""" Primary module for guiding the reduction of long/multi-slit data
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import os

from collections import OrderedDict

from pypit import msgs
from pypit import arparse as settings
from pypit import arcalib
from pypit import arload
from pypit import arutils
from pypit.core import arprocimg
from pypit import arsave
from pypit.core import arwave
from pypit import arsciexp
from pypit.core import arsetup
from pypit import arpixels
from pypit.core import arsort
from pypit import fluxspec
from pypit import scienceimage

from pypit import ardebug as debugger


def ARMS(spectrograph, fitstbl, setup_dict):
    """
    Automatic Reduction of Multislit Data

    Parameters
    ----------
    fitsdict : dict
      Contains relevant information from fits header files
    reuseMaster : bool
      If True, a master frame that will be used for another science frame
      will not be regenerated after it is first made.
      This setting comes with a price, and if a large number of science frames are
      being generated, it may be more efficient to simply regenerate the master
      calibrations on the fly.

    Returns
    -------
    status : int
      Status of the reduction procedure
      0 = Successful full execution
      1 = Successful processing of setup or calcheck
    """
    status = 0

    # Generate sciexp list, if need be (it will be soon)
    sv_std_idx = []
    std_dict = {}
    basenames = []  # For fluxing at the very end
    sciexp = []
    all_sci_ID = fitstbl['sci_ID'].data[fitstbl['science']]  # Binary system: 1,2,4,8, etc.
    for sci_ID in all_sci_ID:
        sciexp.append(arsciexp.ScienceExposure(sci_ID, fitstbl, settings.argflag,
                                               settings.spect, do_qa=True))
        basenames.append(sciexp[-1]._basename)
        std_idx = arsort.ftype_indices(fitstbl, 'standard', sci_ID)
        if (len(std_idx) > 0):
            if len(std_idx) > 1:
                msgs.warn("Will only reduce the first, unique standard for each standard frame!")
            if std_idx[0] not in sv_std_idx:  # Only take the first one
                sv_std_idx.append(std_idx[0])
                # Standard stars
                #std_dict[std_idx[0]] = arsciexp.ScienceExposure(sci_ID, fitstbl, settings.argflag,
                #                                                settings.spect, do_qa=False, idx_sci=std_idx[0])
    numsci = len(sciexp)

    # Init calib dict
    calib_dict = {}


    # Loop on science exposure first
    #  calib frames, e.g. arcs)
    for sc in range(numsci):
        # Init
        sci_output = OrderedDict()
        sci_ID = all_sci_ID[sc]
        scidx = np.where((fitstbl['sci_ID'] == sci_ID) & fitstbl['science'])[0][0]
        msgs.info("Reducing file {0:s}, target {1:s}".format(fitstbl['filename'][scidx],
                                                             fitstbl['target'][scidx])) #slf._target_name))
        # Loop on Detectors
        for kk in range(settings.spect['mosaic']['ndet']):
            det = kk + 1  # Detectors indexed from 1
            if settings.argflag['reduce']['detnum'] is not None:
                if det not in map(int, settings.argflag['reduce']['detnum']):
                    msgs.warn("Skipping detector {:d}".format(det))
                    continue
                else:
                    msgs.warn("Restricting the reduction to detector {:d}".format(det))


            dnum = settings.get_dnum(det)
            msgs.info("Working on detector {:s}".format(dnum))
            sci_output[det] = {}

            # Setup
            namp = settings.spect[dnum]["numamplifiers"]
            setup = arsetup.instr_setup(sci_ID, det, fitstbl, setup_dict, namp, must_exist=True)
            settings.argflag['reduce']['masters']['setup'] = setup

            ###############
            # Get data sections (Could avoid doing this for every sciexp, but it is quick)
            # TODO -_ Clean this up!
            scifile = os.path.join(fitstbl['directory'][scidx],fitstbl['filename'][scidx])
            settings_det = settings.spect[dnum].copy()  # Should include naxis0, naxis1 in this
            datasec_img, naxis0, naxis1 = arprocimg.get_datasec_trimmed(
                settings.argflag['run']['spectrograph'], scifile, det, settings_det,
                naxis0=fitstbl['naxis0'][scidx],
                naxis1=fitstbl['naxis1'][scidx])
            # Binning
            settings_det['binning'] = fitstbl['binning'][0]
            # Yes, this looks goofy.  Is needed for LRIS and DEIMOS for now
            settings.spect[dnum] = settings_det.copy()  # Used internally..
            fitstbl['naxis0'][scidx] = naxis0
            fitstbl['naxis1'][scidx] = naxis1

            # Calib dict
            if setup not in calib_dict.keys():
                calib_dict[setup] = {}

            # TODO -- Update/avoid the following with new settings
            tsettings = settings.argflag.copy()
            tsettings['detector'] = settings.spect[settings.get_dnum(det)]
            try:
                tsettings['detector']['dataext'] = tsettings['detector']['dataext01']  # Kludge; goofy named key
            except KeyError: # LRIS, DEIMOS
                tsettings['detector']['dataext'] = None
            tsettings['detector']['dispaxis'] = settings.argflag['trace']['dispersion']['direction']

            ###############################################################################
            # Prepare for Bias subtraction
            if 'bias' in calib_dict[setup].keys():
                msbias = calib_dict[setup]['bias']
            else:
                # Grab it
                #   Bias will either be an image (ndarray) or a command (str, e.g. 'overscan') or none
                msbias, _ = arcalib.get_msbias(det, setup, sci_ID, fitstbl, tsettings)
                # Save
                calib_dict[setup]['bias'] = msbias

            ###############################################################################
            # Generate a master arc frame
            if 'arc' in calib_dict[setup].keys():
                msarc = calib_dict[setup]['arc']
            else:
                msarc, _ = arcalib.get_msarc(det, setup, sci_ID, spectrograph, fitstbl, tsettings, msbias)
                # Save
                calib_dict[setup]['arc'] = msarc

            ###############################################################################
            # Generate a bad pixel mask (should not repeat)
            if 'bpm' in calib_dict[setup].keys():
                msbpm = calib_dict[setup]['bpm']
            else:
                # Grab it
                msbpm, _ = arcalib.get_mspbm(det, spectrograph, tsettings, msarc.shape,
                                      binning=fitstbl['binning'][scidx],
                                      reduce_badpix=settings.argflag['reduce']['badpix'],
                                      msbias=msbias)
                # Save
                calib_dict[setup]['bpm'] = msbpm

            ###############################################################################
            # Generate an array that provides the physical pixel locations on the detector
            pixlocn = arpixels.gen_pixloc(msarc.shape, det, settings.argflag)

            ###############################################################################
            # Slit Tracing
            if 'trace' in calib_dict[setup].keys():  # Internal
                tslits_dict = calib_dict[setup]['trace']
            else:
                # Setup up the settings (will be Refactored with settings)
                ts_settings = dict(trace=settings.argflag['trace'], masters=settings.argflag['reduce']['masters'])
                ts_settings['masters']['directory'] = settings.argflag['run']['directory']['master']+'_'+ settings.argflag['run']['spectrograph']
                # Get it
                tslits_dict, _ = arcalib.get_tslits_dict( det, setup, spectrograph, sci_ID, ts_settings, tsettings, fitstbl, pixlocn, msbias, msbpm, trim=settings.argflag['reduce']['trim'])
                # Save in calib
                calib_dict[setup]['trace'] = tslits_dict

            ###############################################################################
            # Initialize maskslits
            nslits = tslits_dict['lcen'].shape[1]
            maskslits = np.zeros(nslits, dtype=bool)

            ###############################################################################
            # Generate the 1D wavelength solution
            if 'wavecalib' in calib_dict[setup].keys():
                wv_calib = calib_dict[setup]['wavecalib']
                wv_maskslits = calib_dict[setup]['wvmask']
            elif settings.argflag["reduce"]["calibrate"]["wavelength"] == "pixel":
                msgs.info("A wavelength calibration will not be performed")
                wv_calib = None
                wv_maskslits = np.zeros_like(maskslits, dtype=bool)
            else:
                # Setup up the settings (will be Refactored with settings)
                wvc_settings = dict(calibrate=settings.argflag['arc']['calibrate'], masters=settings.argflag['reduce']['masters'])
                wvc_settings['masters']['directory'] = settings.argflag['run']['directory']['master']+'_'+ settings.argflag['run']['spectrograph']
                nonlinear = settings.spect[settings.get_dnum(det)]['saturation'] * settings.spect[settings.get_dnum(det)]['nonlinear']
                # Get it
                wv_calib, wv_maskslits, _ = arcalib.get_wv_calib(
                    det, setup, spectrograph, sci_ID, wvc_settings, fitstbl, tslits_dict, pixlocn,
                    msarc, nonlinear=nonlinear)
                # Save in calib
                calib_dict[setup]['wavecalib'] = wv_calib
                calib_dict[setup]['wvmask'] = wv_maskslits
            # Mask me
            maskslits += wv_maskslits

            ###############################################################################
            # Derive the spectral tilt
            if 'tilts' in calib_dict[setup].keys():
                mstilts = calib_dict[setup]['tilts']
                wt_maskslits = calib_dict[setup]['wtmask']
            else:
                # Settings kludges
                tilt_settings = dict(tilts=settings.argflag['trace']['slits']['tilts'].copy(),
                                     masters=settings.argflag['reduce']['masters'])
                tilt_settings['tilts']['function'] = settings.argflag['trace']['slits']['function']
                tilt_settings['masters']['directory'] = settings.argflag['run']['directory']['master']+'_'+ settings.argflag['run']['spectrograph']
                # Get it
                mstilts, wt_maskslits, _ = arcalib.get_wv_tilts(
                    det, setup, tilt_settings, settings_det, tslits_dict, pixlocn,
                    msarc, wv_calib, maskslits)
                # Save
                calib_dict[setup]['tilts'] = mstilts
                calib_dict[setup]['wtmask'] = wt_maskslits

            # Mask me
            maskslits += wt_maskslits

            ###############################################################################
            # Prepare the pixel flat field frame
            if settings.argflag['reduce']['flatfield']['perform']:  # Only do it if the user wants to flat field
                if 'normpixelflat' in calib_dict[setup].keys():
                    mspixflatnrm = calib_dict[setup]['normpixelflat']
                    slitprof = calib_dict[setup]['slitprof']
                else:
                    # Settings
                    flat_settings = dict(flatfield=settings.argflag['reduce']['flatfield'].copy(),
                                         slitprofile=settings.argflag['reduce']['slitprofile'].copy(),
                                         combine=settings.argflag['pixelflat']['combine'].copy(),
                                         masters=settings.argflag['reduce']['masters'].copy(),
                                         detector=settings.spect[dnum])
                    flat_settings['masters']['directory'] = settings.argflag['run']['directory']['master']+'_'+ settings.argflag['run']['spectrograph']
                    # Get it
                    mspixflatnrm, slitprof, _ = arcalib.get_msflat(
                        det, setup, sci_ID, fitstbl, tslits_dict, datasec_img,
                        flat_settings, msbias, mstilts)
                    # Save internallly
                    calib_dict[setup]['normpixelflat'] = mspixflatnrm
                    calib_dict[setup]['slitprof'] = slitprof
            else:
                mspixflatnrm = None
                slitprof = None


            ###############################################################################
            # Generate/load a master wave frame
            if 'wave' in calib_dict[setup].keys():
                mswave = calib_dict[setup]['wave']
            else:
                if settings.argflag["reduce"]["calibrate"]["wavelength"] == "pixel":
                    mswave = mstilts * (mstilts.shape[0]-1.0)
                else:
                    # Settings
                    wvimg_settings = dict(masters=settings.argflag['reduce']['masters'].copy())
                    wvimg_settings['masters']['directory'] = settings.argflag['run']['directory']['master']+'_'+ settings.argflag['run']['spectrograph']
                    # Get it
                    mswave, _ = arcalib.get_mswave(
                        setup, tslits_dict, wvimg_settings, mstilts, wv_calib, maskslits)
                # Save internally
                calib_dict[setup]['wave'] = mswave

            # CALIBS END HERE
            ###############################################################################


            ###############
            #  Process and extract the science frame
            msgs.info("Loading science frame")
            sci_image_files = arsort.list_of_files(fitstbl, 'science', sci_ID)
            sci_settings = tsettings.copy()
            # Instantiate
            sciI = scienceimage.ScienceImage(file_list=sci_image_files, datasec_img=datasec_img,
                                             bpm=msbpm, det=det, setup=setup, settings=sci_settings,
                                             maskslits=maskslits, pixlocn=pixlocn, tslits_dict=tslits_dict,
                                             tilts=mstilts, fitstbl=fitstbl, scidx=scidx)
            msgs.sciexp = sciI  # For QA on crash
            # Names and time
            obstime, basename = sciI.init_time_names(settings.spect['mosaic']['camera'],
                            timeunit=settings.spect["fits"]["timeunit"])

            # Process (includes Variance image and CRs)
            dnoise = (settings_det['darkcurr'] * float(fitstbl["exptime"][scidx])/3600.0)
            sciframe, rawvarframe, crmask = sciI._process(
                msbias, mspixflatnrm, apply_gain=True, dnoise=dnoise)

            # Global skysub
            setting_skysub = {}
            setting_skysub['skysub'] = settings.argflag['reduce']['skysub'].copy()
            if settings.argflag['reduce']['skysub']['perform']:
                _ = sciI.global_skysub(setting_skysub)
            else:
                sciI.global_sky = np.zeros_like(sciframe)
                sciI.modelvarframe = np.zeros_like(sciframe)

            # Find objects
            flg_objs = True
            _, nobj = sciI.find_objects()
            if nobj == 0:
                msgs.warn("No objects to extract for science frame" + msgs.newline() + fitstbl['filename'][scidx])
                specobjs = []
                flg_objs = False

            # Another round of sky sub
            if settings.argflag['reduce']['skysub']['perform'] and flg_objs:
                _ = sciI.global_skysub(setting_skysub, use_tracemask=True)

            # Another round of finding objects
            if flg_objs:
                _, nobj = sciI.find_objects()
                if nobj == 0:
                    msgs.warn("No objects to extract for science frame" + msgs.newline() + fitstbl['filename'][scidx])
                    specobjs = []
                    flg_objs = False

            # Extraction
            if flg_objs:
                specobjs, finalvar, finalsky = sciI.extraction(mswave)

            # Flexure correction?
            if settings.argflag['reduce']['flexure']['perform'] and flg_objs:
                if settings.argflag['reduce']['flexure']['method'] is not None:
                    flex_list = arwave.flexure_obj(
                        specobjs, maskslits, settings.argflag['reduce']['flexure']['method'],
                        spectrograph,
                        skyspec_fil = settings.argflag['reduce']['flexure']['spectrum'],
                        mxshft = settings.argflag['reduce']['flexure']['maxshift'])
                    #if not msgs._debug['no_qa']:
                    arwave.flexure_qa(specobjs, maskslits, basename,
                                      det, flex_list)

            # Helio
            # Correct Earth's motion
            if (settings.argflag['reduce']['calibrate']['refframe'] in ['heliocentric', 'barycentric']) and \
                    (settings.argflag['reduce']['calibrate']['wavelength'] != "pixel") and flg_objs:
                if settings.argflag['science']['extraction']['reuse']:
                    msgs.warn("{0:s} correction will not be applied if an extracted science frame exists, and is used".format(settings.argflag['reduce']['calibrate']['refframe']))
                if specobjs is not None:
                    msgs.info("Performing a {0:s} correction".format(settings.argflag['reduce']['calibrate']['refframe']))
                    settings_mosaic = {}
                    settings_mosaic['mosaic'] = settings.spect['mosaic'].copy()
                    vel, vel_corr = arwave.geomotion_correct(specobjs, maskslits, fitstbl, scidx,
                                                             obstime, settings_mosaic,
                                                             settings.argflag['reduce']['calibrate']['refframe'])
                else:
                    msgs.info("There are no objects on detector {0:d} to perform a {1:s} correction".format(
                        det, settings.argflag['reduce']['calibrate']['refframe']))
            else:
                msgs.info("A heliocentric correction will not be performed")

            # Save for output after all detectors are done
            sci_output[det]['sciframe'] = sciframe
            if flg_objs:
                sci_output[det]['specobjs'] = arutils.unravel_specobjs([specobjs])
                sci_output[det]['finalvar'] = finalvar
                sci_output[det]['finalsky'] = finalsky
            else:  # Nothing extracted
                sci_output[det]['specobjs'] = []
                sci_output[det]['finalvar'] = sciI.modelvarframe
                sci_output[det]['finalsky'] = sciI.global_sky

            ######################################################
            # Reduce standard here; only legit if the mask is the same
            std_idx = arsort.ftype_indices(fitstbl, 'standard', sci_ID)
            if len(std_idx) > 0:
                std_idx = std_idx[0]
            else:
                msgs.info("No standard star associated with this science frame")
                continue
            #
            std_image_files = arsort.list_of_files(fitstbl, 'standard', sci_ID)
            if std_idx in std_dict.keys():
                if det in std_dict[std_idx].keys():
                    continue
            else:
                std_dict[std_idx] = {}

            '''
            if stdslf.extracted[det-1] is False:
                # Fill up the necessary pieces
                for iattr in ['pixlocn', 'lordloc', 'rordloc', 'pixcen', 'pixwid', 'lordpix', 'rordpix',
                              'slitpix', 'satmask', 'maskslits', 'mswave']:
                    setattr(stdslf, '_'+iattr, getattr(slf, '_'+iattr))  # Brings along all the detectors, but that is ok
                # Load
                stdframe = arload.load_frames(fitstbl, [std_idx], det, frametype='standard', msbias=msbias)
                stdframe = stdframe[:, :, 0]
                # Reduce
                msgs.info("Processing standard frame")
                arproc.reduce_multislit(stdslf, mstilts, stdframe, msbpm, datasec_img, std_idx, fitstbl, det,
                                        mswave, mspixelflatnrm=mspixflatnrm, standard=True, slitprof=slitprof)
                # Finish
                stdslf.extracted[det-1] = True
            '''
            # Instantiate
            stdI = scienceimage.ScienceImage(file_list=std_image_files, datasec_img=datasec_img,
                                             bpm=msbpm, det=det, setup=setup, settings=sci_settings,
                                             maskslits=maskslits, pixlocn=pixlocn, tslits_dict=tslits_dict,
                                             tilts=mstilts, fitstbl=fitstbl, scidx=std_idx,
                                             objtype='standard')
            # Names and time
            _, std_basename = stdI.init_time_names(settings.spect['mosaic']['camera'],
                                                     timeunit=settings.spect["fits"]["timeunit"])
            # Process (includes Variance image and CRs)
            stdframe, _, _ = stdI._process(msbias, mspixflatnrm, apply_gain=True, dnoise=dnoise)
            # Sky
            _ = stdI.global_skysub(setting_skysub)
            # Find objects
            _, nobj = stdI.find_objects()
            if nobj == 0:
                msgs.warn("No objects to extract for standard frame" + msgs.newline() + fitstbl['filename'][scidx])
                continue
            _ = stdI.global_skysub(setting_skysub, use_tracemask=True)
            # Extract
            stdobjs, _, _ = sciI.extraction(mswave)
            # Save
            std_dict[std_idx][det] = {}
            std_dict[std_idx][det]['basename'] = std_basename
            std_dict[std_idx][det]['specobjs'] = arutils.unravel_specobjs([stdobjs])

        ###########################
        # Write
        # Build the final list of specobjs
        all_specobjs = []
        for key in sci_output:
            all_specobjs += sci_output[key]['specobjs']

        # Write 1D spectra
        save_format = 'fits'
        if save_format == 'fits':
            outfile = settings.argflag['run']['directory']['science']+'/spec1d_{:s}.fits'.format(basename)
            helio_dict = dict(refframe=settings.argflag['reduce']['calibrate']['refframe'],
                              vel_correction=vel_corr)
            arsave.save_1d_spectra_fits(all_specobjs, fitstbl[scidx], outfile,
                                            helio_dict=helio_dict, obs_dict=settings.spect['mosaic'])
            #arsave.save_1d_spectra_fits(slf, fitstbl)
        elif save_format == 'hdf5':
            debugger.set_trace()  # NEEDS REFACTORINGj
            arsave.save_1d_spectra_hdf5(None)
        else:
            msgs.error(save_format + ' is not a recognized output format!')
        # Obj info
        arsave.save_obj_info(all_specobjs, fitstbl, settings.spect, basename,
            settings.argflag['run']['directory']['science'])
        # Write 2D images for the Science Frame
        arsave.save_2d_images(
            sci_output, fitstbl, scidx,
            settings.spect['fits']['headext{0:02d}'.format(1)], setup,
            settings.argflag['run']['directory']['master']+'_'+spectrograph, # MFDIR
            settings.argflag['run']['directory']['science'], basename)

    # Write standard stars
    for std_idx in std_dict.keys():
        all_std_objs = []
        for det in std_dict[std_idx].keys():
            all_std_objs += std_dict[std_idx][det]['specobjs']
        outfile = settings.argflag['run']['directory']['science']+'/spec1d_{:s}.fits'.format(std_dict[std_idx][det]['basename'])
        arsave.save_1d_spectra_fits(all_std_objs, fitstbl[std_idx], outfile,
                                        obs_dict=settings.spect['mosaic'])

    #########################
    # Flux towards the very end..
    #########################
    if settings.argflag['reduce']['calibrate']['flux'] and (len(std_dict) > 0):
        # Standard star (is this a calibration, e.g. goes above?)
        msgs.info("Processing standard star")
        msgs.info("Taking one star per detector mosaic")
        msgs.info("Waited until very end to work on it")
        msgs.warn("You should probably consider using the pypit_flux_spec script anyhow...")

        # Kludge settings
        fsettings = settings.spect.copy()
        fsettings['run'] = settings.argflag['run']
        fsettings['reduce'] = settings.argflag['reduce']
        # Generate?
        if settings.argflag['reduce']['calibrate']['sensfunc']['archival'] == 'None':
            # Take the first standard
            std_idx = list(std_dict.keys())[0]
            # Build the list of stdobjs
            all_std_objs = []
            for det in std_dict[std_idx].keys():
                all_std_objs += std_dict[std_idx][det]['specobjs']
            FxSpec = fluxspec.FluxSpec(settings=fsettings, std_specobjs=all_std_objs,
                                       setup=setup)  # This takes the last setup run, which is as sensible as any..
            sensfunc = FxSpec.master(fitstbl[std_idx])
        else:  # Input by user
            FxSpec = fluxspec.FluxSpec(settings=fsettings,
                                       sens_file=settings.argflag['reduce']['calibrate']['sensfunc']['archival'])
            sensfunc = FxSpec.sensfunc
        # Flux
        msgs.info("Fluxing with {:s}".format(sensfunc['std']['name']))
        for kk, sci_ID in enumerate(all_sci_ID):
            # Load from disk (we zero'd out the object to free memory)
            if save_format == 'fits':
                sci_spec1d_file = settings.argflag['run']['directory']['science']+'/spec1d_{:s}.fits'.format(
                    basenames[kk])
            # Load
            sci_specobjs, sci_header = arload.load_specobj(sci_spec1d_file)
            FxSpec.sci_specobjs = sci_specobjs
            FxSpec.sci_header = sci_header
            # Flux
            FxSpec.flux_science()
            # Over-write
            FxSpec.write_science(sci_spec1d_file)

    return status
