""" Primary module for guiding the reduction of long/multi-slit data
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import os

from collections import OrderedDict

from pypeit import msgs
from pypeit.core import load
from pypeit import utils
from pypeit.core import save
from pypeit.core import wave
from pypeit.core import pypsetup
from pypeit.core import fsort
from pypeit import calibrations
from pypeit import fluxspec
from pypeit import scienceimage
from pypeit.par import pypeitpar

#from pypeit import arparse as settings

from pypeit.spectrographs.spectrograph import Spectrograph
from pypeit.spectrographs.util import load_spectrograph

from pypeit import debugger


def ARMS(fitstbl, setup_dict, par=None, spectrograph=None):
    """
    Automatic Reduction of Multislit Data

    .. todo::
        - improve docstring

    Parameters
    ----------
    spectrograph : str
    fitstbl : Table
      Contains relevant information from fits header files
    setup_dict : dict

    par (`pypeit.par.pypeitpar.PypitPar`): Uber top-level parameter set

    Returns
    -------
    status : int
      Status of the reduction procedure
      0 = Successful full execution
      1 = Successful processing of setup or calcheck
    """
    # TODO: Provide meaningful status values upon return
    status = 0

    # Generate sciexp list, if need be (it will be soon)
    #sv_std_idx = []
    std_dict = {}
    #sciexp = []
    all_sci_ID = fitstbl['sci_ID'].data[fitstbl['science']]  # Binary system: 1,2,4,8, etc.
    numsci = len(all_sci_ID)
    basenames = [None]*numsci  # For fluxing at the very end

    # Spectrometer class
    if spectrograph is None:
        # Set spectrograph from FITS table instrument header
        # keyword.
        if par is not None and par['rdx']['spectrograph'] != fitstbl['instrume'][0]:
            msgs.error('Specified spectrograph does not match instrument in the fits table!')
        _spectrograph = load_spectrograph(spectrograph=fitstbl['instrume'][0])
    elif isinstance(spectrograph, str):
        _spectrograph = load_spectrograph(spectrograph=spectrograph)
    elif isinstance(spectrograph, Spectrograph):
        _spectrograph = spectrograph
    else:
        raise TypeError('Could not instantiate Spectrograph!')

    # Instantiate the parameters
    _par = _spectrograph.default_pypeit_par() if par is None else par
    if not isinstance(_par, pypeitpar.PypeItPar):
        raise TypeError('Input parameters must be a PypitPar instance.')
    required = [ 'rdx', 'calibrations', 'scienceframe', 'objects', 'extract', 'skysubtract',
                 'flexure', 'fluxcalib' ]
    can_be_None = [ 'skysubtract', 'flexure', 'fluxcalib' ]
    _par.validate_keys(required=required, can_be_None=can_be_None)

    # Init calib dict
    caliBrate = calibrations.MultiSlitCalibrations(fitstbl, spectrograph=_spectrograph,
                                                   par=_par['calibrations'],
                                                   save_masters=True, write_qa=True)

    # Loop on science exposure first
    #  calib frames, e.g. arcs)
    for sc in range(numsci):
        # Init
        sci_dict = OrderedDict()  # This needs to be ordered
        sci_dict['meta'] = {}
        sci_dict['meta']['vel_corr'] = 0.
        sci_ID = all_sci_ID[sc]
        scidx = np.where((fitstbl['sci_ID'] == sci_ID) & fitstbl['science'])[0][0]
        msgs.info("Reducing file {0:s}, target {1:s}".format(fitstbl['filename'][scidx],
                                                             fitstbl['target'][scidx]))
                                                             #slf._target_name))
        # Loop on Detectors
        for kk in range(_spectrograph.ndet):
            det = kk + 1  # Detectors indexed from 1
            if _par['rdx']['detnum'] is not None:
                if det not in map(int, _par['rdx']['detnum']):
                    msgs.warn("Skipping detector {:d}".format(det))
                    continue
                else:
                    msgs.warn("Restricting the reduction to detector {:d}".format(det))
            # Setup
            msgs.info("Working on detector {0}".format(det))
            sci_dict[det] = {}

            setup = pypsetup.instr_setup(sci_ID, det, fitstbl, setup_dict,
                                        _spectrograph.detector[det-1]['numamplifiers'],
                                        must_exist=True)

            #-----------------------------------------------------------
            # Calibrations
            #-----------------------------------------------------------
            caliBrate.reset(setup, det, sci_ID, _par['calibrations'])
            # This also instantiates the datasec_img in
            # caliBrate.spectrograph
            datasec_img = caliBrate.get_datasec_img()
            # Bias frame or command
            msbias = caliBrate.get_bias()
            # Arc Image
            msarc = caliBrate.get_arc()
            # Bad pixel mask
            msbpm = caliBrate.get_bpm()
            # Physical pixel locations on the detector
            pixlocn = caliBrate.get_pixlocn()
            # Slit Tracing
            tslits_dict, maskslits = caliBrate.get_slits(arms=True)
            if tslits_dict is None: # No slits
                msgs.warn('No slits found!')
                continue
            # Generate the 1D wavelength solution
            wv_calib, maskslits = caliBrate.get_wv_calib()
            # Derive the spectral tilt
            mstilts, maskslits = caliBrate.get_tilts()
            # Prepare the pixel flat field frame
            mspixflatnrm, slitprof = caliBrate.get_pixflatnrm()
            # Generate/load a master wave frame
            mswave = caliBrate.get_wave()

            # Above could also be run with::
            #   caliBrate.run_the_steps()
            # And we could pull the items we need out of it or just use the attributes
            #-----------------------------------------------------------

            #-----------------------------------------------------------
            # Science frames
            #-----------------------------------------------------------
            msgs.info("Working on the science frame")
            sci_image_files = fsort.list_of_files(fitstbl, 'science', sci_ID)

            # Instantiate
            sciI = scienceimage.ScienceImage(_spectrograph, file_list=sci_image_files,
                                             frame_par=_par['scienceframe'],
                                             trace_objects_par=_par['objects'],
                                             extract_objects_par=_par['extract'],
                                             tslits_dict=tslits_dict, tilts=mstilts, det=det,
                                             setup=setup, datasec_img=datasec_img, bpm=msbpm,
                                             maskslits=maskslits, pixlocn=pixlocn,
                                             fitstbl=fitstbl, scidx=scidx)

            msgs.sciexp = sciI  # For QA on crash

            # Names and time
            obstime, basename = sciI.init_time_names(_spectrograph.camera,
                                                     timeunit=_spectrograph.timeunit)
            if basenames[sc] is None:
                basenames[sc] = basename

            # Process (includes Variance image and CRs)
            sciframe, rawvarframe, crmask = sciI.process(msbias, mspixflatnrm, apply_gain=True,
                                                         trim=caliBrate.par['trim'])

            # Global skysub
            if _par['skysubtract'] is None:
                # These are set as attributes of sciI
                sciI.global_sky = np.zeros_like(sciframe)
                sciI.modelvarframe = np.zeros_like(sciframe)
            else:
                # The call to global_skysub initializes the attributes
                # of sciI directly
                global_sky, modelvarframe = sciI.global_skysub(
                                        bspline_spacing=_par['skysubtract']['bspline_spacing'])

            # Find objects
            nobj = sciI.find_objects()[1]
            if nobj == 0:
                msgs.warn('No objects to extract for science frame' + msgs.newline()
                          + fitstbl['filename'][scidx])
                specobjs, flg_objs = [], None
            else:
                flg_objs = True  # Objects were found

            # Another round of sky sub
            if _par['skysubtract'] is not None and flg_objs:
                global_sky, modelvarframe = sciI.global_skysub(
                                        bspline_spacing=_par['skysubtract']['bspline_spacing'],
                                                               use_tracemask=True)

            # Another round of finding objects
            if flg_objs:
                nobj = sciI.find_objects()[1]
                if nobj == 0:
                    msgs.warn('No objects to extract for science frame' + msgs.newline()
                              + fitstbl['filename'][scidx])
                    specobjs, flg_objs = [], None

            # Extraction
            if flg_objs:
                specobjs, finalvar, finalsky = sciI.extraction(mswave)

            # Flexure correction?
            if _par['flexure'] is not None and flg_objs and _par['flexure']['method'] is not None:
                sky_file, sky_spectrum = _spectrograph.archive_sky_spectrum()
                flex_list = wave.flexure_obj(specobjs, maskslits, _par['flexure']['method'],
                                               sky_spectrum, sky_file=sky_file,
                                               mxshft=_par['flexure']['maxshift'])
                # QA
                wave.flexure_qa(specobjs, maskslits, basename, det, flex_list)

            # Helio
            # Correct Earth's motion
            vel_corr = -999999.9
            if (caliBrate.par['wavelengths']['frame'] in ['heliocentric', 'barycentric']) and \
                        (caliBrate.par['wavelengths']['reference'] != 'pixel') and flg_objs:
                if _par['extract']['reuse']:
                    msgs.warn('{0} correction'.format(caliBrate.par['wavelengths']['frame'])
                              + 'will not be applied if an extracted science frame exists, '
                              + 'and is used')
                if specobjs is not None:
                    msgs.info("Performing a {0} correction".format(
                                            caliBrate.par['wavelengths']['frame']))

                    vel, vel_corr = wave.geomotion_correct(specobjs, maskslits, fitstbl, scidx,
                                                             obstime,
                                                             _spectrograph.telescope['longitude'],
                                                             _spectrograph.telescope['latitude'],
                                                             _spectrograph.telescope['elevation'],
                                                             caliBrate.par['wavelengths']['frame'])
                else:
                    msgs.info('There are no objects on detector {0} to perform a '.format(det)
                              + '{1} correction'.format(caliBrate.par['wavelengths']['frame']))
            else:
                msgs.info('A wavelength reference-frame correction will not be performed.')

            # Save for outputing (after all detectors are done)
            sci_dict[det]['sciframe'] = sciframe
            if vel_corr > -999999.9:
                sci_dict['meta']['vel_corr'] = vel_corr
            if flg_objs:
                sci_dict[det]['specobjs'] = utils.unravel_specobjs([specobjs])
                sci_dict[det]['finalvar'] = finalvar
                sci_dict[det]['finalsky'] = finalsky
            else:  # Nothing extracted
                sci_dict[det]['specobjs'] = []
                sci_dict[det]['finalvar'] = sciI.modelvarframe
                sci_dict[det]['finalsky'] = sciI.global_sky
            #-----------------------------------------------------------

            #-----------------------------------------------------------
            # Standard star frames
            #-----------------------------------------------------------
            # Can only reduce these frames if the mask is the same
            std_idx = fsort.ftype_indices(fitstbl, 'standard', sci_ID)
            if len(std_idx) > 0:
                std_idx = std_idx[0]
            else:
                msgs.info("No standard star associated with this science frame")
                continue
            #
            msgs.info("Processing standard star")
            std_image_files = fsort.list_of_files(fitstbl, 'standard', sci_ID)
            if std_idx in std_dict.keys():
                if det in std_dict[std_idx].keys():
                    continue
            else:
                std_dict[std_idx] = {}

            if _par['calibrations']['standardframe'] is None:
                msgs.warn('No standard frame parameters provided.  Using default parameters.')

            # Instantiate for the Standard
            # TODO: Uses the same trace and extraction parameter sets used for the science
            # frames.  Should these be different for the standards?
            stdI = scienceimage.ScienceImage(_spectrograph, file_list=std_image_files,
                                             frame_par=_par['calibrations']['standardframe'],
                                             trace_objects_par=_par['objects'],
                                             extract_objects_par=_par['extract'],
                                             tslits_dict=tslits_dict, tilts=mstilts, det=det,
                                             setup=setup, datasec_img=datasec_img, bpm=msbpm,
                                             maskslits=maskslits, pixlocn=pixlocn,
                                             fitstbl=fitstbl, scidx=std_idx, objtype='standard')

            # Names and time
            std_basename = stdI.init_time_names(_spectrograph.camera,
                                                timeunit=_spectrograph.timeunit)[1]
            # Process (includes Variance image and CRs)
            stdframe = stdI.process(msbias, mspixflatnrm, apply_gain=True,
                                    trim=caliBrate.par['trim'])[0]
            # Sky
            stdI.global_skysub(bspline_spacing=_par['skysubtract']['bspline_spacing'])
            # Find objects
            nobj = stdI.find_objects()[1]
            if nobj == 0:
                msgs.warn('No objects to extract for standard frame' + msgs.newline()
                          + fitstbl['filename'][scidx])
                continue
            stdI.global_skysub(bspline_spacing=_par['skysubtract']['bspline_spacing'],
                               use_tracemask=True)
            # Extract
            stdobjs = stdI.extraction(mswave)[0]
            # Save for fluxing and output later
            std_dict[std_idx][det] = {}
            std_dict[std_idx][det]['basename'] = std_basename
            std_dict[std_idx][det]['specobjs'] = utils.unravel_specobjs([stdobjs])
            #-----------------------------------------------------------

        #---------------------------------------------------------------
        # Write the output for this exposure
        #---------------------------------------------------------------
        # Build the final list of specobjs and vel_corr
        all_specobjs = []
        for key in sci_dict:
            if key in ['meta']:
                continue
            #
            try:
                all_specobjs += sci_dict[key]['specobjs']
            except KeyError:  # No object extracted
                continue

        if len(all_specobjs) == 0:
            msgs.warn('No objects!')
            continue

        # Write 1D spectra
        save_format = 'fits'
        if save_format == 'fits':
            outfile = os.path.join(_par['rdx']['scidir'], 'spec1d_{:s}.fits'.format(basename))
            helio_dict = dict(refframe='pixel'
                              if caliBrate.par['wavelengths']['reference'] == 'pixel' 
                              else caliBrate.par['wavelengths']['frame'],
                              vel_correction=sci_dict['meta']['vel_corr'])
            save.save_1d_spectra_fits(all_specobjs, fitstbl[scidx], outfile,
                                        helio_dict=helio_dict, telescope=_spectrograph.telescope)
#        elif save_format == 'hdf5':
#            debugger.set_trace()  # NEEDS REFACTORING
#            arsave.save_1d_spectra_hdf5(None)
        else:
            msgs.error(save_format + ' is not a recognized output format!')
        # Obj info
        save.save_obj_info(all_specobjs, fitstbl, _spectrograph, basename, _par['rdx']['scidir'])
        # Write 2D images for the Science Frame
        save.save_2d_images(sci_dict, fitstbl, scidx, _spectrograph.primary_hdrext,
                              setup, caliBrate.master_dir, _par['rdx']['scidir'], basename)
        #---------------------------------------------------------------

    #-------------------------------------------------------------------
    # Write standard stars at the very end
    #-------------------------------------------------------------------
    for std_idx in std_dict.keys():
        all_std_objs = []
        for det in std_dict[std_idx].keys():
            all_std_objs += std_dict[std_idx][det]['specobjs']
        outfile = os.path.join(_par['rdx']['scidir'],
                               'spec1d_{:s}.fits'.format(std_dict[std_idx][det]['basename']))
        save.save_1d_spectra_fits(all_std_objs, fitstbl[std_idx], outfile,
                                    telescope=_spectrograph.telescope)
    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    # Flux calibrate
    #-------------------------------------------------------------------
    if _par['fluxcalib'] is None or len(std_dict) == 0:
        msgs.info('Flux calibration is not performed.')
        return status

    if _par['fluxcalib'] is None and len(std_dict) > 0:
        msgs.info('Flux calibration parameters not provided.  Standards not used.')
        return status

    # Standard star (is this a calibration, e.g. goes above?)
    msgs.info("Taking one star per detector mosaic")
    msgs.info("Waited until very end to work on it")
    msgs.warn("You should probably consider using the pypeit_flux_spec script anyhow...")

    # Get the sensitivity function
    if _par['fluxcalib']['sensfunc'] is None:
        # Take the first standard
        std_idx = list(std_dict.keys())[0]
        # Build the list of stdobjs
        all_std_objs = []
        for det in std_dict[std_idx].keys():
            all_std_objs += std_dict[std_idx][det]['specobjs']
        FxSpec = fluxspec.FluxSpec(std_specobjs=all_std_objs, spectrograph=_spectrograph,
                                   setup=setup, root_path=caliBrate.master_root,
                                   mode=_par['calibrations']['masters'])
        sensfunc = FxSpec.master(fitstbl[std_idx])
    else:
        # User provided it
        FxSpec = fluxspec.FluxSpec(sens_file=_par['fluxcalib']['sensfunc'],
                                   spectrograph=_spectrograph, root_path=caliBrate.master_root,
                                   mode=_par['calibrations']['masters'])
        sensfunc = FxSpec.sensfunc

    # Apply the flux calibration
    msgs.info("Fluxing with {:s}".format(sensfunc['std']['name']))
    for kk, sci_ID in enumerate(all_sci_ID):
        # Load from disk (we zero'd out the object to free memory)
        if save_format == 'fits':
            sci_spec1d_file = os.path.join(_par['rdx']['scidir'],
                                           'spec1d_{:s}.fits'.format(basenames[kk]))

        # Load
        sci_specobjs, sci_header = load.load_specobj(sci_spec1d_file)
        # TODO: (KBW) I'm wary of this kind of approach.  We want
        # FluxSpec to check that its internals make sense and this
        # bypasses any of that checking.
        FxSpec.sci_specobjs = sci_specobjs
        FxSpec.sci_header = sci_header

        # Flux
        FxSpec.flux_science()
        # Over-write
        FxSpec.write_science(sci_spec1d_file)
    #-------------------------------------------------------------------

    return status

