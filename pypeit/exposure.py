import numpy as np

from astropy.io import fits

from pypeit import log
from pypeit import PypeItError
from pypeit import outputfiles
from pypeit.core import parse
from pypeit.display import display
from pypeit.history import History
from pypeit import slittrace
from pypeit import specobjs
from pypeit import spec2dobj
from pypeit import calibrations

from pypeit import pypeit_steps

from IPython import embed


def adjust_for_slitmask(sciImg_dict:dict, spectrograph, fitstbl, par, 
                        frame0:int, all_specobjs_objfind,
                        calib_slits:list):
    """
    Adjust slitmask information for a given set of science images.

    This function processes slitmask design information, applies spatial 
    flexure corrections, computes dither offsets, and updates the slitmask 
    calibration data. It also matches slitmask design to detected objects 
    and adds undetected objects based on the design.

    Args:
        sciImg_dict (:obj:`dict`): A dict of science image objects, one for each 
            detector, containing information such as spatial flexure and 
            detector properties.
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph object
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            The class holding the metadata for all the frames in this PypeIt run.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the reduction process, 
            including slitmask and object finding parameters.
        frame0 (:obj:`int`): The index of the current frame in the FITS table.
        binning (:obj:`str`): The binning string (e.g., '2,2') specifying the 
            spectral and spatial binning.
        all_specobjs_objfind (:obj:`list`): A list of SpecObj objects representing 
            detected objects for all detectors.
        calib_slits (:obj:`list`): A list of SlitTraceSet objects containing 
            slitmask calibration data for all detectors.

    Returns:
        tuple: 
            list: Updated list of SlitTraceSet objects with adjusted slitmask 
            information, including object positions and offsets.
            SpecObjs: Updated SpecObjs object with matched and added objects
    """
    # get object positions from slitmask design and slitmask offsets for all the detectors
    spat_flexure = np.array([sciImg_dict[ss].spat_flexure for ss in sciImg_dict])
    # Grab platescale with binning
    bin_spec, bin_spat = parse.parse_binning(fitstbl['binning'][frame0])
    platescale = np.array([sciImg_dict[ss].detector.platescale*bin_spat for ss in sciImg_dict])
    # get the dither offset if available and if desired
    dither_off = None
    if par['reduce']['slitmask']['use_dither_offset']:
        if 'dithoff' in fitstbl.keys():
            dither_off = fitstbl['dithoff'][frame0]

    calib_slits = slittrace.get_maskdef_objpos_offset_alldets(
        all_specobjs_objfind, calib_slits, spat_flexure, platescale,
        par['calibrations']['slitedges']['det_buffer'],
        par['reduce']['slitmask'], dither_off=dither_off)
    # determine if slitmask offsets exist and compute an average offsets over all the detectors
    calib_slits = slittrace.average_maskdef_offset(
        calib_slits, platescale[0], spectrograph.list_detectors(mosaic='MSC' in calib_slits[0].detname))
    # slitmask design matching and add undetected objects
    all_specobjs_objfind = slittrace.assign_addobjs_alldets(
        all_specobjs_objfind, calib_slits, spat_flexure, platescale,
        par['reduce']['slitmask'], par['reduce']['findobj']['find_fwhm'])

    return calib_slits, all_specobjs_objfind
    
                            

def process_exposure(spectrograph, fitstbl, par, frames:list, 
                     calib_ID:str, detectors:list, calibrations_path:str, 
                     bg_frames:list=None): 
    """
    Process all detectors for a given exposure.

    Calls :func:`~pypeit.pypeit_steps.process_one_det` for each detector.
    

    This function processes exposure data for a list of detectors by performing
    the necessary reduction steps and generating science images and background
    reduced science images.

    Args:
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph object
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            The class holding the metadata for all the frames in this PypeIt run.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the reduction process, 
            including slitmask and object finding parameters.
        frames (:obj:`list`): A list of frame indices to process.
        calib_ID (:obj:`str`): The calibration group ID.
        detectors (:obj:`list`): A list of detector indices to process.
        calibrations_path (:obj:`str`): The path to the calibration files.
        bg_frames (:obj:`list`, optional): A list of background frame indices. Defaults to None.

    Returns:
        tuple: A tuple containing:
            - sciImg_dict (:obj:`dict`): A dictionary where the keys are detector indices
              and the values are the corresponding science images.
            - bkg_redux_sciimg_dict (:obj:`dict`): A dictionary where the keys are detector
              indices and the values are the corresponding background reduced science images.
    """

    # dict of sciImg
    sciImg_dict = {}
    # list of bkg_redux_sciimg
    bkg_redux_sciimg_dict = {}

    # Loop on the detectors
    for det in detectors:
        log.info(f'Reducing detector {det}')

        # Process
        sciImg, bkg_redux_sciimg = pypeit_steps.process_one_det(
            spectrograph, fitstbl, par, frames, det, calib_ID, calibrations_path, 
            bg_frames=bg_frames)

        # List em up
        sciImg_dict[det] = sciImg
        bkg_redux_sciimg_dict[det] = bkg_redux_sciimg

    # Return
    return sciImg_dict, bkg_redux_sciimg_dict

def findobj_on_exposure(sciImg_dict:dict, bkg_redux_sciimg_dict:dict,
                        spectrograph, fitstbl, par,
                        frames:list, detectors:list, calib_ID:str, 
                        calibrations_path:str, 
                        std_outfile:str=None, bkg_redux=False, 
                        find_negative=False, show=False):
    """
    Identifies objects on a set of exposures for the specified detectors.
    This function loops over the provided detectors, 
    science images, and identifies objects using the `findobj_on_det` method.
    It returns the initial sky model for each detector and a collection of
    identified spectral objects.

    Calls :func:`~pypeit.pypeit_steps.process_one_det` for each detector.
    
    Args:
        sciImg_dict (:obj:`dict`): A dict of science image objects, one for each 
            detector, containing information such as spatial flexure and 
            detector properties.
        bkg_redux_sciimg_dict (dict): A dict of background image objects, one for each
            detector, containing information such as spatial flexure and
            detector properties.
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph object
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            The class holding the metadata for all the frames in this PypeIt run.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the reduction process, 
        frames (:obj:`list`): A list of frame indices to process.
        calib_ID (:obj:`str`): The calibration group ID.
        detectors (list): List of detectors to process.
        calibrations_path (:obj:`str`): The path to the calibration files.
        std_outfile (str, optional): Path to the standard star output file. Defaults to None.
        bkg_redux (bool, optional): If True, perform A-B background subtraction. Defaults to False.
        find_negative (bool, optional): If True, search for negative objects. Defaults to False.
        show (bool, optional): If True, display intermediate results. Defaults to False.

    Returns:
        tuple:
            - final_sky_dict (dict): Dictionary containing the final sky model; keys are each detector.
            - bkg_redux_final_sky_dict (dict): Dictionary containing the final bkg_redux sky model;
              keys are each detector.
            - all_specobjs_objfind (SpecObjs): Collection of all identified spectral objects.
            - all_silts (list): List of Slits objects, detector by detector
            - sciImg_dict (dict): Dictionary containing updated sciImg objects with global spectral
              flexure and scaleimg information.
    """
    
    # Output
    initial_sky_dict = {}
    final_sky_dict = {}
    bkg_redux_final_sky_dict = {}

    # container for specobjs during first loop (objfind)
    all_specobjs_objfind = specobjs.SpecObjs()
    all_slits = []
    # container for the ObjFind for each detector
    all_objfinds = []

    # #####################################
    # find objects + initial sky subtraction
    # Loop on the detectors
    for det in detectors:
        # Grab the science image
        sciImg = sciImg_dict[det]

        # Run
        initial_sky, sobjs_obj, objFind = \
            pypeit_steps.findobj_on_det(
                sciImg, spectrograph, fitstbl, par, frames, calib_ID, det, 
                calibrations_path, 
                bkg_redux=bkg_redux, find_negative=find_negative, show=show, 
                std_outfile=std_outfile)

        # Slits
        all_slits.append(objFind.slits)
        # objFind
        all_objfinds.append(objFind)

        # Store em
        initial_sky_dict[det] = initial_sky
        if len(sobjs_obj)>0:
            all_specobjs_objfind.add_sobj(sobjs_obj)

    # #####################################
    # slitmask stuff
    if par['reduce']['slitmask']['assign_obj']:
        frame0 = frames[0]
        all_slits, all_specobjs_objfind = adjust_for_slitmask(
            sciImg_dict, spectrograph, fitstbl,
            par, frame0,
            all_specobjs_objfind, all_slits)

    # #####################################
    # final sky subtraction
    for i,det in enumerate(detectors):
        # Load some useful objects
        this_objfind = all_objfinds[i]
        bkg_redux_sciimg = bkg_redux_sciimg_dict[det]
        initial_sky = initial_sky_dict[det]
        # update the slits in objfind
        this_objfind.slits = all_slits[i]

        final_global_sky, bkg_redux_global_sky, this_objfind = \
            pypeit_steps.finalize_sky_det(spectrograph, fitstbl, par, frames[0],
                     det, this_objfind, initial_sky, all_specobjs_objfind,
                     bkg_redux_sciimg=bkg_redux_sciimg, bkg_redux=bkg_redux, show=show)

        # store the final skies
        final_sky_dict[det] = final_global_sky
        bkg_redux_final_sky_dict[det] = bkg_redux_global_sky
        # Update the sciImg with the scaleImg information
        sciImg_dict[det].rel_scaleImg = this_objfind.scaleimg
        # and the global spectral flexure shift
        sciImg_dict[det].flex_shift = this_objfind.slitshift
        # TODO: RJC please check if slitshift here should be assigned or added to sciImg_dict[det].flex_shift

        # update here slits.mask since global_skysub modify reduce_bpm, and we need to propagate it into extraction
        flagged_slits = np.where(this_objfind.reduce_bpm)[0]
        if len(flagged_slits) > 0:
            all_slits[i].mask[flagged_slits] = \
                all_slits[i].bitmask.turn_on(all_slits[i].mask[flagged_slits], 'BADSKYSUB')

    # Return
    return final_sky_dict, bkg_redux_final_sky_dict, all_specobjs_objfind, all_slits, sciImg_dict

def extract_exposure(sciImg_dict:dict, spectrograph, fitstbl, par, frames:list, detectors,
                     calib_ID:str, calibrations_path:str, all_specobjs_objfind, 
                     final_sky_dict:dict, bkg_redux_final_sky_dict:dict,
                     calib_slits, bkg_redux:bool=False, find_negative:bool=False):

    """
    Extracts spectral data from a set of science images and performs background subtraction, 
    sky subtraction, and object extraction.

    Calls :func:`~pypeit.pypeit_steps.extract_det` for each detector.

    Parameters:
        sciImg_dict (:obj:`dict`): A dict of science image objects, one for each 
            detector, containing information such as spatial flexure and 
            detector properties.
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph object
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            The class holding the metadata for all the frames in this PypeIt run.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the reduction process, 
        frames (:obj:`list`): A list of frame indices to process.
        detectors: List of detectors to process.
        calib_ID (:obj:`str`): The calibration group ID.
        calibrations_path (:obj:`str`): The path to the calibration files.
        all_specobjs_objfind: SpecObjs object containing objects found during object finding.
        final_sky_dict (:obj:`dict`): Dictionary containing final sky models for each detector.
        bkg_redux_final_sky_dict (:obj:`dict`): Dictionary containing final bkg_redux sky models for each detector.
        calib_slits (:obj:`list`): A list of SlitTraceSet objects for all detectors updated
        after slitmask adjustment and findobj+sky subtraction.
        bkg_redux (bool, optional): If True, perform background reduction. Default is False.
        find_negative (bool, optional): If True, search for negative objects. Default is False.

    Returns:
        tuple: A tuple containing:
            - all_spec2d (AllSpec2DObj): Container for all extracted 2D spectra.
            - all_specobjs_extract (SpecObjs): Container for all extracted objects.
    """

    # Container for all the Spec2DObj
    all_spec2d = spec2dobj.AllSpec2DObj()
    all_spec2d['meta']['bkg_redux'] = bkg_redux
    all_spec2d['meta']['find_negative'] = find_negative
    # container for specobjs during second loop (extraction)
    all_specobjs_extract = specobjs.SpecObjs()

    # Extract
    for i,det in enumerate(detectors):
        detname = sciImg_dict[det].detector.name

        # TODO: pass back the background frame, pass in background
        # files as an argument. extract one takes a file list as an
        # argument and instantiates science within
        if all_specobjs_objfind.nobj > 0:
            all_specobjs_on_det = all_specobjs_objfind[all_specobjs_objfind.DET == detname]
        else:
            all_specobjs_on_det = all_specobjs_objfind

        # Extract
        all_spec2d[detname], tmp_sobjs = pypeit_steps.extract_det(
            spectrograph, fitstbl, par, frames, det,
            calib_ID, calibrations_path,
            sciImg_dict[det],
            final_sky_dict[det],
            all_specobjs_on_det,
            calib_slits[i],
            bkg_redux_final_sky=bkg_redux_final_sky_dict[det],
            bkg_redux=bkg_redux,
            find_negative=find_negative)

        # Hold em
        if tmp_sobjs.nobj > 0:
            all_specobjs_extract.add_sobj(tmp_sobjs)

        # Add calibration associations to the SpecObjs object
        all_specobjs_extract.calibs = calibrations.Calibrations.get_association(
                                fitstbl, spectrograph, calibrations_path,
                                fitstbl[frames[0]]['setup'],
                                fitstbl.find_frame_calib_groups(frames[0])[0], det,
                                must_exist=True, proc_only=True)

    # Return
    return all_spec2d, all_specobjs_extract

def reduce_exposure(spectrograph, fitstbl, par, frames, calib_ID, 
                    calibrations_path: str, bg_frames=None, 
                    reuse_calibs: bool = True,
                    run_state: dict = None, std_outfile=None,
                    show: bool = False):
    """
    Reduce a set of exposures for a given spectrograph and calibration setup.

    This function performs the full reduction process for a set of science frames,
    including background subtraction, object finding, sky subtraction, and extraction.

    Args:
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph object
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            The class holding the metadata for all the frames in this PypeIt run.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the reduction process, 
            including slitmask and object finding parameters.
        frames (:obj:`list`): A list of frame indices to process.
        calib_ID (:obj:`str`): The calibration group ID.
        calibrations_path (:obj:`str`): The path to the calibration files.
        bg_frames (:obj:`list`, optional): A list of background frame indices. Defaults to None.
        reuse_calibs (bool, optional): Whether to reuse existing calibrations. Defaults to True.
        run_state (dict, optional): Dictionary to track the state of the reduction process. Defaults to None.
        std_outfile (str, optional): Path to the standard star output file. Defaults to None.
        show (bool, optional): Whether to display intermediate results (e.g., using Ginga). Defaults to False.

    Returns:
        tuple:
            - all_spec2d (dict): Dictionary containing the 2D spectral data for all detectors.
            - all_specobjs_extract (list): List of extracted spectral objects.

    Notes:
        - The function handles background subtraction and finding negative traces if applicable.
        - Calibrations are performed for each detector, and unsuccessful calibrations are skipped.
        - Object finding, sky subtraction, and extraction are performed for the specified frames.
        - Slitmask adjustments are applied if enabled in the parameters.
    """

    # if show is set, clear the ginga channels at the start of each new sci_ID
    if show:
        # TODO: Put this in a try/except block?
        display.clear_all(allow_new=True)

    # Prep for background subtraction and finding negative traces
    has_bg, bkg_redux, find_negative = pypeit_steps.set_bkg_negative(
        fitstbl, par, bg_frames)

    # Print status message
    lstr = f'Reducing target {fitstbl["target"][frames[0]]}\n'
    # TODO: Print these when the frames are actually combined,
    # backgrounds are used, etc?
    lstr += 'Combining frames:\n'
    for iframe in frames:
        lstr += f'{fitstbl["filename"][iframe]}\n'
    log.info(lstr)
    if has_bg:
        bg_lstr = ''
        for iframe in bg_frames:
            bg_lstr += f'{fitstbl["filename"][iframe]}\n'
        bg_lstr = '\nUsing background from frames:\n' + bg_lstr
        log.info(bg_lstr)

    # Find the detectors to reduce
    detectors = spectrograph.select_detectors(subset=par['rdx']['detnum'] if par['rdx']['slitspatnum'] is None 
                                              else par['rdx']['slitspatnum'])
    log.info(f'Detectors to work on: {detectors}')

    # #####################################
    # Calibrations
    for det in detectors:
        log.info(f'Calibrating detector {det}')
        # run/load calibration
        caliBrate =  pypeit_steps.calib_one(spectrograph, fitstbl, par, det, calib_ID, calibrations_path,
              show=show, run_state=run_state, reuse_calibs=reuse_calibs)
        if not caliBrate.success:
            log.warning(
                f'Calibrations for detector {det} were unsuccessful!  The step that failed was '
                f'{caliBrate.failed_step}.  Continuing by skipping this detector.'
            )
            # Remove from list of detectors
            detectors.remove(det)
            continue

    # #####################################
    # Process or load processed frames
    sciImg_dict, bkg_redux_sciimg_dict = process_exposure(
            spectrograph, fitstbl, par, frames, calib_ID,
                detectors, calibrations_path, 
                bg_frames=bg_frames) 

    # #####################################
    # Find objects +  sky
    final_sky_dict, bkg_redux_final_sky_dict, all_specobjs_find, calib_slits, sciImg_dict = \
        findobj_on_exposure(sciImg_dict, bkg_redux_sciimg_dict,
                            spectrograph,
                            fitstbl,
                            par, frames, detectors,
                            calib_ID, calibrations_path,
                            std_outfile=std_outfile,
                            bkg_redux=bkg_redux,
                            find_negative=find_negative,
                            show=show)

    # #####################################
    # Extract
    all_spec2d, all_specobjs_extract = extract_exposure(
        sciImg_dict, spectrograph, fitstbl,
        par, frames, detectors, calib_ID,
        calibrations_path, all_specobjs_find,
        final_sky_dict, bkg_redux_final_sky_dict,
        calib_slits, bkg_redux=bkg_redux,
        find_negative=find_negative)

    # Return
    return all_spec2d, all_specobjs_extract

def save_exposure(spectrograph, fitstbl, par, 
                  frame:int, all_spec2d:spec2dobj.AllSpec2DObj,
                  all_specobjs:specobjs.SpecObjs, 
                  calibrations_path:str,
                  history:History=None,
                  in_update_det:int=None,
                  skip_write_2d:bool=False):
    """
    Save the outputs from extraction for a given exposure

    Args:
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph object
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            The class holding the metadata for all the frames in this PypeIt run.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the reduction process, 
            including slitmask and object finding parameters.
        frame (:obj:`int`):
            0-indexed row in the metadata table with the frame
            that has been reduced.
        all_spec2d(:class:`~pypeit.spec2dobj.AllSpec2DObj`):
            The 2D reduced spectrum objects.
        all_specobjs (:class:`~pypeit.specobjs.SpecObjs`):
            The 1D spectral extraction objects.
        calibrations_path (:obj:`str`): The path to the calibration files.
        history (:class:`~pypeit.history.History`, optional):
            History entries to be added to fits header
        in_update_det (:obj:`int`, optional):
            Detector number to use when writing the output files.
            If not None, overwrite the value of par['rdx']['detnum']
        skip_write_2d (:obj:`bool`, optional):
            Skip writing the 2D spectrum to disk
    """
    # Check for the Science/ directory
    science_path = outputfiles.science_path(par)
    if not science_path.is_dir():
        science_path.mkdir()

    # Determine the headers
    row_fitstbl = fitstbl[frame]
    # Need raw file header information
    rawfile = fitstbl.frame_paths(frame)
    head2d = fits.getheader(rawfile, ext=spectrograph.primary_hdrext)


    # NOTE: There are some gymnastics here to keep from altering
    # self.par['rdx']['detnum'].  I.e., I can't just set update_det =
    # self.par['rdx']['detnum'] because that can alter the latter if I don't
    # deepcopy it...
    if in_update_det is not None:
        update_det = in_update_det
    elif par['rdx']['detnum'] is None:
        update_det = None
    elif isinstance(par['rdx']['detnum'], list):
        update_det = [spectrograph.allowed_mosaics.index(d)+1 
                        if isinstance(d, tuple) else d for d in par['rdx']['detnum']]
    else:
        update_det = par['rdx']['detnum']

    subheader = spectrograph.subheader_for_spec(row_fitstbl, head2d)
    # 1D spectra
    if all_specobjs.nobj > 0 and not par['reduce']['extraction']['skip_extraction']:
        # Spectra
        outfile1d = outputfiles.spec_output_file(fitstbl, par, frame)
        # TODO
        #embed(header='deal with the following for maskIDs;  713 of pypeit')
        all_specobjs.write_to_fits(subheader, outfile1d,
                                    update_det=update_det,
                                    slitspatnum=par['rdx']['slitspatnum'],
                                    history=history)
        # Info
        outfiletxt = outputfiles.spec_output_file(fitstbl, par,
                                            frame, ext='.txt')
        # TODO: Note we re-read in the specobjs from disk to deal with situations where
        # only a single detector is run in a second pass but in the same reduction directory.
        # This was to address Issue #1116 in PR #1154. Slightly inefficient, but only other
        # option is to re-work write_info to also "append"
        sobjs = specobjs.SpecObjs.from_fitsfile(outfile1d, chk_version=False)
        sobjs.write_info(outfiletxt, spectrograph.pypeline)

    if skip_write_2d:
        return

    # 2D spectra
    outfile2d = outputfiles.spec_output_file(fitstbl, par, frame, twod=True)

    # Build header
    pri_hdr = all_spec2d.build_primary_hdr(head2d, spectrograph,
                                            redux_path=par['rdx']['redux_path'],
                                            calib_dir=calibrations_path,
                                            subheader=subheader,
                                            history=history)

    # Write
    all_spec2d.write_to_fits(outfile2d, pri_hdr=pri_hdr,
                                update_det=update_det,
                                slitspatnum=par['rdx']['slitspatnum'])
