"""
Module providing functions that support primary reduction steps.

.. include:: ../include/links.rst
"""
import os
from pathlib import Path
import numpy as np
import copy

from astropy.table import Table

from pypeit import log
from pypeit import PypeItError
from pypeit.calibframe import CalibFrame
from pypeit.images import buildimage
from pypeit import specobjs
from pypeit import find_objects
from pypeit import extraction
from pypeit.manual_extract import ManualExtractionObj
from pypeit import spec2dobj
from pypeit.core import wave

# local_skyregions
from pypeit.core import skysub
from pypeit import io

from pypeit import slittrace
from pypeit import calibrations

from linetools import utils as ltu

from IPython import embed

def get_sci_metadata(spectrograph, fitstbl, frame:int, det):
    """
    Grab the meta data for a given science frame and specific detector

    Args:
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph instance
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            The class holding the metadata for all the frames in this PypeIt run.
        frame (:obj:`int`):
            The index of the frame in the fitstbl
        det (:obj:`int`):
            The detector number (1-indexed)

    Returns:
        5 objects are returned::
            - str: Object type;  science or standard
            - str: Setup/configuration string
            - `astropy.time.Time`_: Time of observation
            - str: Basename of the frame
            - str: Binning of the detector

    """

    # Set binning, obstime, basename, and objtype
    binning = fitstbl['binning'][frame]
    obstime  = fitstbl.construct_obstime(frame)
    basename = fitstbl.construct_basename(frame, obstime=obstime)
    types  = fitstbl['frametype'][frame].split(',')
    if 'science' in types:
        objtype_out = 'science'
    elif 'standard' in types:
        objtype_out = 'standard'
    else:
        raise PypeItError(
            'get_sci_metadata() should only be run on standard or science frames.  Types of this '
            f'frame are: {types}'
        )
    calib_key = CalibFrame.construct_calib_key(fitstbl['setup'][frame],
                                                fitstbl['calib'][frame],
                                                spectrograph.get_det_name(det))
    return objtype_out, calib_key, obstime, basename, binning

def set_bkg_negative(fitstbl, par, bg_frames:list):
    """
    Determine background subtraction settings for a reduction.

    Args:
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            The class holding the metadata for all the frames in this PypeIt run.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the reduction process, 
            including slitmask and object finding parameters.
        bg_frames : list
            A list of indices corresponding to the background frames in the
            `fitstbl`. If empty or None, no background subtraction is performed.

    Returns:
        tuple: Returns three objects:
            has_bg (:obj:`bool`):
                True if background frames are provided and non-empty, False otherwise.
            bkg_redux (:obj:`bool`):
                True if the reduction involves background subtraction, False otherwise.
            find_negative (:obj:`bool`):
                Indicates whether to find negative objects during the reduction. The
                default behavior depends on the type of the background frames
                ('science' or 'sky') unless explicitly overridden by the parameter
                `par['reduce']['findobj']['find_negative']`.
    """
    has_bg = True if bg_frames is not None and len(bg_frames) > 0 else False
    # Is this an b/g subtraction reduction?
    if has_bg:
        bkg_redux = True
        # The default is to find_negative objects if the bg_frames are
        # classified as "science", and to not find_negative objects if the
        # bg_frames are classified as "sky". This can be explicitly
        # overridden if par['reduce']['findobj']['find_negative'] is set to
        # something other than the default of None.
        find_negative = (('science' in fitstbl['frametype'][bg_frames[0]]) |
                                ('standard' in fitstbl['frametype'][bg_frames[0]])) \
                        if par['reduce']['findobj']['find_negative'] is None else \
                            par['reduce']['findobj']['find_negative']
    else:
        bkg_redux = False
        find_negative= False

    # Return
    return has_bg, bkg_redux, find_negative

def calib_one(spectrograph, fitstbl, par, det, calib_ID, calibrations_path:str, 
              reuse_calibs:bool=True,
              qa_path:str=None, show:bool=False, run_state:dict=None,
              stop_at_step:str=None):
    """
    Run Calibration for a single detector, calib_ID pair

    Args:
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph instance
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            The class holding the metadata for all the frames in this PypeIt run.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the reduction, including slitmask and
            object finding parameters.
        det (:obj:`int`):
            Detector number (1-indexed)
        calib_ID (:obj:`str`):
            The calibration group ID
        calibrations_path (:obj:`str`):
            Path to the calibration files
        reuse_calibs (:obj:`bool`, optional):
            If True, reuse existing calibration files if they exist.
        qa_path (:obj:`str`, optional):
            Path to the QA files; if None, use the default path
            defined by the parameters.
        show (:obj:`bool`, optional):
            Show the QA during processing
        run_state (:obj:`dict`, optional):
            A dictionary containing the current state of the reduction.
            If None, a new empty dictionary is created.
        stop_at_step (:obj:`str`, optional):
            Run only up to this calibration step.

    Returns:
        caliBrate (:class:`~pypeit.calibrations.Calibrations`)

    """
    if qa_path is None:
        qa_path = os.path.join(par['rdx']['redux_path'], par['rdx']['qadir'])

    # Handle frames
    in_grp = fitstbl.find_calib_group(calib_ID)
    frame_indx = np.arange(len(fitstbl))
    grp_frames = frame_indx[in_grp]

    # Instantiate Calibrations class
    user_slits = slittrace.merge_user_slit(par['rdx']['slitspatnum'],
                                            par['rdx']['maskIDs'])
    log.info(f'Building/loading calibrations for detector {det}')
    caliBrate = calibrations.Calibrations.get_instance(
        fitstbl, par['calibrations'], spectrograph, calibrations_path, 
        calib_ID, grp_frames[0], det,
        qadir=qa_path,
        reuse_calibs=reuse_calibs, show=show, user_slits=user_slits,
        chk_version=par['rdx']['chk_version'])
        #, state=run_state)

    # Check
    if stop_at_step is not None and stop_at_step not in caliBrate.steps:
        raise PypeItError(
            f"Requested stop_at_step={stop_at_step} is not a valid calibration step.\n Allowed "
            f"steps are: {caliBrate.steps}"
        )
        
    # Run
    caliBrate.run_the_steps(stop_at_step=stop_at_step)

    # Success?
    if not caliBrate.success:
        log.warning(
            f'Calibrations for detector {det} were unsuccessful!  The step that failed was '
            f'{caliBrate.failed_step}.  Continuing to next detector.'
        )

    return caliBrate


def process_one_det(spectrograph, fitstbl, par, frames:list, 
                    det, calib_ID:str, calibrations_path:str, bg_frames:list=None,
                    sci_outfile:str=None, bkg_outfile:str=None): 
    """
    Process a single detector for a given set of frames.

    This function handles the image processing for a specific detector, including
    loading calibrations, building the science image, and optionally subtracting
    a background image.

    Args:
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph instance
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            The class holding the metadata for all the frames in this PypeIt run.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the reduction, including slitmask and
            object finding parameters.
        frames (:obj:`list`):
            List of indices corresponding to the science frames in the
            `fitstbl` to be processed together.
        det (:obj:`int`):
            Detector number (1-indexed)
        calib_ID (:obj:`str`):
            The calibration group ID
        calibrations_path (:obj:`str`):
            Path to the calibration files
        bg_frames (:obj:`list`, optional):
            List of indices corresponding to the background frames in the
            `fitstbl`. If None or empty, no A-B background subtraction is performed.
            Default is None.
        sci_outfile (:obj:`str`, optional):
            The science output file, if this is a science reduction.
            Default is None.
        bkg_outfile (:obj:`str`, optional):
            The background output file, if this is a background
            subtraction reduction. Default is None. 

    Returns:
        tuple:
            - sciImg (:class:`~pypeit.images.pypeitimage.PypeItImage`): 
              The processed science image, with background
              subtraction applied if `bg_frames` is provided.
            - bkg_redux_sciimg (:class:`~pypeit.images.pypeitimage.PypeItImage` or None):
              The science image without
              background subtraction, used to generate a global sky model. This is
              a dictionary with `image` and `ivar` keys if `bg_frames` is provided,
              otherwise it is None.
    """

    # Grab some meta-data needed for the reduction from the fitstbl
    objtype, setup, obstime, basename, binning \
            = get_sci_metadata(spectrograph, fitstbl, frames[0], det)

    # Grab the calibrations
    caliBrate = load_calibrations_for_frame(
        spectrograph, fitstbl, par, frames[0], det, calib_ID, calibrations_path)

    log.info(f"Image processing begins for {basename} on det={det}")

    # Is this a standard star?
    std_redux = objtype == 'standard'
    frame_par = par['calibrations']['standardframe'] \
                    if std_redux else par['scienceframe']

    # Build Science image
    sci_files = fitstbl.frame_paths(frames)
    sciImg = buildimage.buildimage_fromlist(
        spectrograph, det, frame_par,
        sci_files, bias=caliBrate.msbias, bpm=caliBrate.msbpm,
        dark=caliBrate.msdark,
        scattlight=caliBrate.msscattlight,
        flatimages=caliBrate.flatimages,
        slits=caliBrate.slits,  # For flexure correction
        ignore_saturation=False)

    # get no bkg subtracted sciImg to generate a global sky model without bkg subtraction.
    # it's a dictionary with only `image` and `ivar` keys if bkg_redux=False, otherwise it's None
    bkg_redux_sciimg = None

    # Background Image?
    if bg_frames is not None and len(bg_frames) > 0:
        # get no bkg subtracted sciImg
        bkg_redux_sciimg = sciImg
        # Build the background image
        bg_file_list = fitstbl.frame_paths(bg_frames)
        # TODO I think we should create a separate self.par['bkgframe'] parameter set to hold the image
        # processing parameters for the background frames.  This would allow the user to specify different
        # parameters for the background frames than for the science frames.  
        bgimg = buildimage.buildimage_fromlist(spectrograph, det, frame_par, bg_file_list,
                                                bpm=caliBrate.msbpm,
                                                bias=caliBrate.msbias,
                                                dark=caliBrate.msdark,
                                                scattlight=caliBrate.msscattlight,
                                                flatimages=caliBrate.flatimages,
                                                slits=caliBrate.slits,
                                                ignore_saturation=False)

        # NOTE: If the spatial flexure exists for sciImg, the subtraction
        # function propagates that to the subtracted image, ignoring any
        # spatial flexure determined for the background image.
        sciImg = bkg_redux_sciimg.sub(bgimg)

    # Write out the science image?
    if sci_outfile is not None:
        # Generate the folder?
        if not sci_outfile.parent.is_dir():
            sci_outfile.parent.mkdir()
        sciImg.to_file(sci_outfile, overwrite=True)
        log.info(f'Wrote intermediate science image to {sci_outfile}')

    # Write out the background image?
    if bkg_outfile is not None and bkg_redux_sciimg is not None:
        bkg_redux_sciimg.to_file(bkg_outfile, overwrite=True)
        log.info(f'Wrote intermediate background image to {bkg_outfile}')

    # Return
    return sciImg, bkg_redux_sciimg

def findobj_on_det(sciImg, spectrograph, fitstbl, par, frames:list, calib_ID:str, 
                   det, calibrations_path:str, std_outfile:str=None, 
                   bkg_redux=False, find_negative=False, show:bool=False):
    """
    Perform object finding on a specific detector.

    This function is responsible for identifying objects on a given detector
    during the data reduction process. It utilizes calibration data, metadata,
    and reduction parameters to locate and characterize objects in the science image.

    Args:
        sciImg (:class:`~pypeit.images.pypeitimage.PypeItImage`):
            Data container that holds a single image from a
            single detector and its related images (e.g. ivar, mask)
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph instance
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            The class holding the metadata for all the frames in this PypeIt run.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the reduction, including slitmask and
            object finding parameters.
        frames (:obj:`list`):
            List of indices corresponding to the science frames in the
            `fitstbl` to be processed together.
        calib_ID (:obj:`str`):
            The calibration group ID
        det (:obj:`int`):
            Detector number (1-indexed)
        calibrations_path (:obj:`str`):
            Path to the calibration files
        std_outfile (:obj:`str`, optional):
            The standard star output file, if this is a standard star
            reduction. Default is None.
        bkg_redux (:obj:`bool`, optional):
            Indicates whether the reduction involves background subtraction.
            Default is False.
        find_negative (:obj:`bool`, optional):
            Indicates whether to find negative objects during the reduction.
            Default is False.
        show (:obj:`bool`, optional):
            Show the QA during processing. Default is False.

    Returns:
        tuple: A tuple containing:
            - initial_sky (object): The initial sky model.
            - sobjs_obj (object): The detected objects.
            - objFind (object): The object finding instance used for the reduction.
    """

    # Grab some meta-data needed for the reduction from the fitstbl
    objtype, setup, obstime, basename, binning \
            = get_sci_metadata(spectrograph, fitstbl, frames[0], det)

    # Is this a standard star?
    std_redux = objtype == 'standard'
    if std_redux is False and std_outfile is not None:
        std_trace = specobjs.get_std_trace(spectrograph.get_det_name(det), 
                                        std_outfile)
    else:
        std_trace = None
    log.info("Object finding begins for {} on det={}".format(basename, det))

    # Grab the calibrations
    caliBrate = load_calibrations_for_frame(
        spectrograph, fitstbl, par, frames[0], det, calib_ID, calibrations_path)

    log.info(f'Reducing detector {det}')

    # Instantiate Reduce object
    # Required for pypeline specific object
    # At instantiaton, the fullmask in self.sciImg is modified
    objFind = instantiate_objfind(sciImg, spectrograph, fitstbl,
                                  par, frames, det, caliBrate,
                                  bkg_redux,
                                  find_negative, show=show)
    # Do it
    initial_sky, sobjs_obj = objFind.run(std_trace=std_trace, 
                                         show_peaks=show)

    return initial_sky, sobjs_obj, objFind

def finalize_sky_det(spectrograph, fitstbl, par, frame,
                     det, objFind, initial_sky, all_specobjs_objfind,
                     bkg_redux_sciimg=None, bkg_redux=False, show=False):
    """
    Finalize sky subtraction for a specific detector.

    This function performs the final global sky subtraction for a given detector.
    It also updates the slit mask based on bad sky subtraction flags.

    Args:
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph instance.
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            The class holding the metadata for all the frames in this PypeIt run.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the reduction process.
        frame (:obj:`int`):
            The index of the frame in the fitstbl.
        det (:obj:`int`):
            Detector number (1-indexed).
        objFind (:class:`~pypeit.find_objects.FindObjects`):
            The object finding instance for this detector.
        initial_sky (`numpy.ndarray`_):
            The initial sky model for this detector.
        all_specobjs_objfind (:class:`~pypeit.specobjs.SpecObjs`):
            All spectral objects found during object finding.
        bkg_redux_sciimg (:class:`~pypeit.images.pypeitimage.PypeItImage`, optional):
            Background-reduced science image, or None if background reduction is
            not being performed.
        bkg_redux (:obj:`bool`, optional):
            Indicates whether background reduction is being performed. Default is False.
        show (:obj:`bool`, optional):
            Show the QA during processing. Default is False.

    Returns:
        tuple: A tuple containing:
            - final_global_sky (`numpy.ndarray`_): The final global sky model for this detector.
            - bkg_redux_global_sky (`numpy.ndarray`_ or None): The background-reduced global sky
              model, or None if bkg_redux is False.
            - objFind (:class:`~pypeit.find_objects.FindObjects`): The updated object finding
                class for this detector.
    """

    objtype, setup, obstime, basename, binning \
            = get_sci_metadata(spectrograph, fitstbl, frame, det)
    detname = spectrograph.get_det_name(det)

    # Get objects on this detector
    if all_specobjs_objfind.nobj > 0:
        all_specobjs_on_det = all_specobjs_objfind[all_specobjs_objfind.DET == detname]
    else:
        all_specobjs_on_det = all_specobjs_objfind

    # Determine if final global sky subtraction should be performed
    skymask = None
    if 'standard' in objtype or \
            par['reduce']['findobj']['skip_skysub'] or \
            par['reduce']['findobj']['skip_final_global'] or \
            par['reduce']['skysub']['user_regions'] is not None:
        final_global_sky = initial_sky
    else:
        # Update the skymask
        skymask = objFind.create_skymask(all_specobjs_on_det)
        final_global_sky = objFind.global_skysub(previous_sky=initial_sky,
                                                 skymask=skymask, show=show,
                                                 reinit_bpm=False)

    # Get the bkg_redux_global_sky
    bkg_redux_global_sky = None
    if bkg_redux and bkg_redux_sciimg is not None:
        skymask = objFind.create_skymask(all_specobjs_on_det) if skymask is None else skymask
        # DO NOT reinit_bpm, nor update_crmask
        bkg_redux_global_sky = objFind.global_skysub(skymask=skymask,
                                                     bkg_redux_sciimg=bkg_redux_sciimg,
                                                     reinit_bpm=False, update_crmask=False, show=show)

    return final_global_sky, bkg_redux_global_sky, objFind



def load_calibrations_for_frame(spectrograph, fitstbl, par, frame, det, 
                                calib_ID, calibrations_path:str):
    """
    Load calibrations for a specific frame and detector.

    This function initializes the Calibrations class, runs the calibration steps,
    and ensures that the calibrations are successfully loaded for the given frame
    and detector.

    Args:
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph instance
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            The class holding the metadata for all the frames in this PypeIt run.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the reduction, including slitmask and
            object finding parameters.
        frame (:obj:`int`):
            The index of the frame in the fitstbl
        det (:obj:`int`):
            Detector number (1-indexed)
        calib_ID (:obj:`str`):
            The calibration group ID
        calibrations_path (:obj:`str`):
            Path to the calibration files

    Returns:
        :class:`~pypeit.calibrations.Calibrations`:
            The loaded calibrations for the specified frame and detector.

    Raises:
        PypeItError: If the calibrations for the specified detector are unsuccessful.
    """

    # Instantiate Calibrations class
    user_slits = slittrace.merge_user_slit(par['rdx']['slitspatnum'],
                                            par['rdx']['maskIDs'])
    caliBrate = calibrations.Calibrations.get_instance(
        fitstbl, par['calibrations'], spectrograph,
        calibrations_path, 
        calib_ID, frame, det,
        reuse_calibs=True, user_slits=user_slits,
        chk_version=par['rdx']['chk_version'])
    #caliBrate.set_config(frame, det, par['calibrations'])
    caliBrate.run_the_steps(reload_only=True)

    if not caliBrate.success:
        raise PypeItError(
            f'Calibrations for detector {det} were unsuccessful!  The step that failed was '
            f'{caliBrate.failed_step}.'
        )

    return caliBrate


def load_skyregions(spectrograph, fitstbl, par, frame, det, caliBrate,
                    calibrations_path:str, scifile:str=None, initial_slits=False):
    """
    Generate or load sky regions, if defined by the user.

    Sky regions are defined by the internal provided parameters; see
    ``user_regions`` in :class:`~pypeit.par.pypeitpar.SkySubPar`.  If
    included in the pypeit file like so,

    .. code-block:: ini

        [reduce]
            [[skysub]]
                user_regions = :25,75:

    The first and last 25% of all slits are used as sky.  If the user has
    used the ``pypeit_skysub_regions`` GUI to generate a sky mask for a
    given frame, this can be searched for and loaded by setting the
    parameter to ``user``:

    .. code-block:: ini

        [reduce]
            [[skysub]]
                user_regions = user

    Parameters
    ----------
    initial_slits : :obj:`bool`, optional
        Flag to use the initial slits before any tweaking based on the
        slit-illumination profile; see
        :func:`~pypeit.slittrace.SlitTraceSet.select_edges`.

    Returns
    -------
    skymask : `numpy.ndarray`_
        A boolean array used to select sky pixels; i.e., True is a pixel
        that corresponds to a sky region.  If the ``user_regions`` parameter
        is not set (or an empty string), the returned value is None.
    """
    if par['reduce']['skysub']['user_regions'] in [None, '']:
        return None

    # Flexure
    spat_flexure = None

    # First priority given to user_regions first
    if par['reduce']['skysub']['user_regions'] == 'user':
        # Build the file name
        calib_key = CalibFrame.construct_calib_key(
                            fitstbl['setup'][frame],
                            CalibFrame.ingest_calib_id(fitstbl['calib'][frame]),
                            spectrograph.get_det_name(det))
        regfile = buildimage.SkyRegions.construct_file_name(
            calib_key, calib_dir=calibrations_path, 
            basename=io.remove_suffix(scifile))
        regfile = Path(regfile).absolute()
        if not regfile.exists():
            raise PypeItError(
                f'Unable to find SkyRegions file: {regfile} . Create a SkyRegions frame using '
                'pypeit_skysub_regions, or change the user_regions to the percentage format.  '
                'See documentation.'
            )
        log.info(f'Loading SkyRegions file: {regfile}')
        return buildimage.SkyRegions.from_file(regfile).image.astype(bool)

    skyregtxt = par['reduce']['skysub']['user_regions']
    if isinstance(skyregtxt, list):
        skyregtxt = ",".join(skyregtxt)
    log.info(f'Generating skysub mask based on the user defined regions: {skyregtxt}')
    # NOTE : Do not include spatial flexure here!
    #        It is included when generating the mask in the return statement below
    slits_left, slits_right, _ \
        = caliBrate.slits.select_edges(initial=initial_slits, flexure=None)

    maxslitlength = np.max(slits_right-slits_left)
    # Get the regions
    status, regions = skysub.read_userregions(skyregtxt, caliBrate.slits.nslits, maxslitlength)
    if status == 1:
        raise PypeItError(
            "Unknown error in sky regions definition. Please check the value:\n" + skyregtxt
        )
    elif status == 2:
        raise PypeItError(
            "Sky regions definition must contain a percentage range, and therefore must "
            "contain a ':'"
        )
    # Generate and return image
    return skysub.generate_mask(spectrograph.pypeline, regions, caliBrate.slits,
                                slits_left, slits_right, spat_flexure=spat_flexure)


def extract_det(spectrograph, fitstbl, par,
                frames, det, calib_ID:str, calibrations_path:str,
                sciImg, final_sky, sobjs_obj, calib_slits,
                bkg_redux_final_sky=None, bkg_redux:bool=False,
                find_negative:bool=False,
                show:bool=False):
    """
    Extract Objects in a single exposure/detector pair

    sci_ID and det need to have been set internally prior to calling this method

    Args:
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph instance
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            The class holding the metadata for all the frames in this PypeIt run.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the reduction, including slitmask and
            object finding parameters.
        frames (:obj:`list`):
            List of frames to extract; stacked if more than one
            is provided
        det (:obj:`int`):
            Detector number (1-indexed)
        calib_ID (:obj:`str`):
            The calibration group ID
        calibrations_path (:obj:`str`):
            Path to the calibration files
        sciImg (:class:`~pypeit.images.pypeitimage.PypeItImage`):
            Data container that holds a single image from a
            single detector and its related images (e.g. ivar, mask)
        final_sky (`numpy.ndarray`_):
            Final global sky model
        sobjs_obj (:class:`~pypeit.specobjs.SpecObjs`):
            List of objects found during `run_objfind`
        calib_slits (:class:`~pypeit.slittrace.SlitTraceSet`):
            If provided, use these slits instead of those from the
            calibrations.
        bkg_redux_final_sky (`numpy.ndarray`_, optional):
            Final global sky model for the background-reduced science image.
            Default is None.
        bkg_redux (:obj:`bool`, optional):
            Indicates whether the reduction involves background subtraction.
            Default is False.
        find_negative (:obj:`bool`, optional):
            Indicates whether to find negative objects during the reduction.
            Default is False.
        show (:obj:`bool`, optional):
            Show the QA during processing. Default is False.0

    Returns:
        tuple: A tuple containing:
            - spec2DObj (:class:`~pypeit.spec2dobj.Spec2DObj`):
              The 2D spectrum object containing the 2D spectral images
            - sobj (:class:`~pypeit.specobjs.SpecObjs`):
                The extracted 1D spectra
            
    """
    # Grab some meta-data needed for the reduction from the fitstbl
    #self.objtype, self.setup, self.obstime, self.basename, self.binning \
    #        = self.get_sci_metadata(frames[0], det)
    objtype, setup, obstime, basename, binning \
            = get_sci_metadata(spectrograph, fitstbl, frames[0], det)

    # Grab the calibrations
    caliBrate = load_calibrations_for_frame(
        spectrograph, fitstbl, par, frames[0], det, calib_ID, calibrations_path)
    # update slits
    caliBrate.slits = calib_slits

    # Is this a standard star?
    std_redux = 'standard' in objtype

    # Each spec2d file includes the slits object with unique flagging
    #  for extraction failures.  So we make a copy here before those flags
    #  are modified.
    maskdef_designtab = caliBrate.slits.maskdef_designtab
    slits = copy.deepcopy(caliBrate.slits)
    slits.maskdef_designtab = None
    # this is only used for IFU reductions currently
    scaleImg = sciImg.rel_scaleImg

    if not par['reduce']['extraction']['skip_extraction']:
        log.info(f"Extraction begins for {basename} on det={det}")
        # Instantiate Reduce object
        # Required for pipeline specific object
        # At instantiation, the fullmask in self.sciImg is modified
        # TODO Are we repeating steps in the init for FindObjects and Extract??
        exTract = extraction.Extract.get_instance(
            sciImg, slits, sobjs_obj, spectrograph,
            par, objtype, global_sky=final_sky,
            bkg_redux_global_sky=bkg_redux_final_sky,
            waveTilts=caliBrate.wavetilts, wv_calib=caliBrate.wv_calib, 
            flatimages=caliBrate.flatimages,
            bkg_redux=bkg_redux, 
            return_negative=par['reduce']['extraction']['return_negative'],
            std_redux=std_redux, basename=basename, show=show)
        # Perform the extraction
        skymodel, bkg_redux_skymodel, objmodel, ivarmodel, outmask, sobjs, waveImg,\
            tilts, slits = exTract.run()
        slitgpm = np.logical_not(exTract.extract_bpm)
        slitshift = exTract.slitshift
    else:
        log.info(f"Extraction skipped for {basename} on det={det}")
        # Since the extraction was not performed, fill the arrays with the best available information
        skymodel, bkg_redux_skymodel, objmodel, ivarmodel, outmask, sobjs = \
            final_sky, \
            bkg_redux_final_sky, \
            np.zeros_like(sciImg.image), \
            np.copy(sciImg.ivar), \
            sciImg.fullmask, \
            sobjs_obj
        slitgpm = (slits.mask == 0)
        slitshift = sciImg.flex_shift
        # get slitmask (same as in extract)
        slitmask = slits.slit_img(flexure=sciImg.spat_flexure, exclude_flag=slits.bitmask.exclude_for_reducing)
        # get spat_flexure and tilts (same as in extract)
        _spat_flexure = 0. if sciImg.spat_flexure is None else sciImg.spat_flexure
        _tilts_spat_flexure = 0. if caliBrate.wavetilts.spat_flexure is None else caliBrate.wavetilts.spat_flexure
        tilts = caliBrate.wavetilts.fit2tiltimg(slitmask, flexure=_tilts_spat_flexure)
        waveImg = caliBrate.wv_calib.build_waveimg(tilts, slits, spat_flexure=sciImg.spat_flexure, spec_flexure=slitshift)

    # Apply a reference frame correction to each object and the waveimg
    vel_corr, waveImg = refframe_correct(spectrograph, par, slits, 
                                              fitstbl["ra"][frames[0]], 
                                              fitstbl["dec"][frames[0]],
                                              obstime, slitgpm=slitgpm, 
                                              waveimg=waveImg, sobjs=sobjs)

    # TODO -- Do this upstream
    # Tack on wavelength RMS
    for sobj in sobjs:
        iwv = np.where(caliBrate.wv_calib.spat_ids == sobj.SLITID)[0][0]
        sobj.WAVE_RMS =caliBrate.wv_calib.wv_fits[iwv].rms

    # Construct table of spectral flexure
    spec_flex_table = Table()
    spec_flex_table['spat_id'] = slits.spat_id
    spec_flex_table['sci_spec_flexure'] = slitshift

    # Construct the Spec2DObj
    spec2DObj = spec2dobj.Spec2DObj(sciimg=sciImg.image,
                                    ivarraw=sciImg.ivar,
                                    skymodel=skymodel,
                                    bkg_redux_skymodel=bkg_redux_skymodel,
                                    objmodel=objmodel,
                                    ivarmodel=ivarmodel,
                                    scaleimg=scaleImg,
                                    waveimg=waveImg,
                                    bpmmask=outmask,
                                    detector=sciImg.detector,
                                    sci_spat_flexure=sciImg.spat_flexure,
                                    sci_spec_flexure=spec_flex_table,
                                    vel_corr=vel_corr,
                                    vel_type=par['calibrations']['wavelengths']['refframe'],
                                    tilts=tilts,
                                    slits=slits,
                                    wavesol=caliBrate.wv_calib.wave_diagnostics(print_diag=False),
                                    maskdef_designtab=maskdef_designtab)
    spec2DObj.process_steps = sciImg.process_steps

    spec2DObj.calibs = calibrations.Calibrations.get_association(
                                fitstbl, spectrograph, 
                                caliBrate.calib_dir, #calibrations_path,
                                fitstbl[frames[0]]['setup'],
                                fitstbl.find_frame_calib_groups(frames[0])[0], det,
                                must_exist=True, proc_only=True)

    # QA
    spec2DObj.gen_qa()

    # Return
    return spec2DObj, sobjs

def instantiate_objfind(sciImg, spectrograph, fitstbl, par, frames, det,
                        caliBrate, bkg_redux, find_negative, 
                        show:bool=False):
    """
    Instantiate the FindObjects class for object finding in spectroscopic data.

    Args:
        sciImg (:class:`~pypeit.images.pypeitimage.PypeItImage`):
            Data container that holds a single image from a
            single detector and its related images (e.g. ivar, mask)
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph instance
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            The class holding the metadata for all the frames in this PypeIt run.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the reduction, including slitmask and
            object finding parameters.
        frames (:obj:`list`):
            List of indices corresponding to the science frames in the
            `fitstbl` to be processed together.
        det (:obj:`int`):
            Detector number (1-indexed)
        caliBrate (:class:`~pypeit.calibrations.Calibrations`):
            The calibration data for the current frame and detector.
        bkg_redux (:obj:`bool`):
            Indicates whether the reduction involves background subtraction.
        find_negative (:obj:`bool`):
            Indicates whether to find negative objects during the reduction.
        show (:obj:`bool`, optional):
            Show the QA during processing. Default is False.

    Returns:
        FindObjects: An instance of the FindObjects class configured for object finding.

    Notes:
        - This function handles manual extraction if specified in the FITS table.
        - It also builds an initial sky mask if user-defined sky regions are provided.
        - The FindObjects instance is initialized with the relevant parameters for object finding.
    """
    objtype, setup, obstime, basename, binning \
            = get_sci_metadata(spectrograph, fitstbl, frames[0], det)
    std_redux = objtype == 'standard'

    # Deal with manual extraction
    row = fitstbl[frames[0]]
    manual_obj = ManualExtractionObj.by_fitstbl_input(
        row['filename'], row['manual'], spectrograph) if len(row['manual'].strip()) > 0 else None

    # TODO - Move this into FindObjects class
    if par['reduce']['skysub']['user_regions'] in [None, '']:
        initial_skymask = None
    else:
        # Build the initial sky mask
        initial_skymask = load_skyregions(
            spectrograph, fitstbl, par, frames[0], det,
            caliBrate, str(caliBrate.calib_dir), initial_slits=spectrograph.pypeline != 'SlicerIFU',
            scifile=fitstbl.frame_paths(frames[0]))
            

    objFind = find_objects.FindObjects.get_instance(
        sciImg, caliBrate.slits,
        spectrograph, par, objtype,
        wv_calib=caliBrate.wv_calib,
        waveTilts=caliBrate.wavetilts,
        initial_skymask=initial_skymask,
        bkg_redux=bkg_redux,
        manual=manual_obj,
        find_negative=find_negative,
        std_redux=std_redux,
        basename=basename,
        show=show)

    return objFind


def refframe_correct(spectrograph, par, slits, ra, dec, obstime, slitgpm=None, 
                     waveimg=None, sobjs=None):
    """ Correct the calibrated wavelength to the user-supplied reference frame

    Args:
        slits (:class:`~pypeit.slittrace.SlitTraceSet`):
            Slit trace set object
        ra (float, str):
            Right Ascension
        dec (float, str):
            Declination
        obstime (`astropy.time.Time`_):
            Observation time
        slitgpm (`numpy.ndarray`_, None, optional):
            1D boolean array indicating the good slits (True). If None, the gpm will be taken from slits
        waveimg (`numpy.ndarray`_, optional)
            Two-dimensional image specifying the wavelength of each pixel
        sobjs (:class:`~pypeit.specobjs.SpecObjs`, None, optional):
            Spectrally extracted objects

    """
    if slitgpm is None:
        slitgpm = (slits.mask == 0)
    # Correct Telescope's motion
    refframe = par['calibrations']['wavelengths']['refframe']
    vel_corr = 0.0
    if refframe in ['heliocentric', 'barycentric'] \
            and par['calibrations']['wavelengths']['reference'] != 'pixel':
        log.info(f"Performing a {par['calibrations']['wavelengths']['refframe']} correction")
        # Calculate correction
        radec = ltu.radec_to_coord((ra, dec))
        vel, vel_corr = wave.geomotion_correct(radec, obstime,
                                                spectrograph.telescope['longitude'],
                                                spectrograph.telescope['latitude'],
                                                spectrograph.telescope['elevation'],
                                                refframe)
        # Apply correction to objects
        log.info(f'Applying {refframe} correction = {vel:0.5f} km/s')
        if (sobjs is not None) and (sobjs.nobj != 0):
            # Loop on slits to apply
            gd_slitord = slits.slitord_id[slitgpm]
            for slitord in gd_slitord:
                indx = sobjs.slitorder_indices(slitord)
                this_specobjs = sobjs[indx]
                # Loop on objects
                for specobj in this_specobjs:
                    if specobj is None:
                        continue
                    specobj.apply_helio(vel_corr, refframe)

        # Apply correction to wavelength image
        if waveimg is not None:
            waveimg *= vel_corr
    else:
        log.info('A wavelength reference frame correction will not be performed.')

    # Return the value of the correction and the corrected wavelength image
    return vel_corr, waveimg


