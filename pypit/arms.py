""" Primary module for guiding the reduction of long/multi-slit data
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import os

from collections import OrderedDict

from pypit import msgs
from pypit import arparse as settings
from pypit import arload
from pypit import arutils
from pypit.core import arprocimg
from pypit.core import arsave
from pypit.core import arwave
from pypit.core import arsetup
from pypit import arpixels
from pypit.core import arsort
from pypit import calibrations
from pypit.spectrographs import bpmimage
from pypit import flatfield
from pypit import fluxspec
from pypit import traceslits
from pypit import wavecalib
from pypit import waveimage
from pypit import wavetilts
from pypit import scienceimage
from pypit.spectrographs import io

from pypit import ardebug as debugger


def ARMS(spectrograph, fitstbl, setup_dict):
    """
    Automatic Reduction of Multislit Data

    Parameters
    ----------
    spectrograph : str
    fitstbl : Table
      Contains relevant information from fits header files
    setup_dict : dict

    Returns
    -------
    status : int
      Status of the reduction procedure
      0 = Successful full execution
      1 = Successful processing of setup or calcheck
    """
    status = 0

    # Generate sciexp list, if need be (it will be soon)
    #sv_std_idx = []
    std_dict = {}
    #sciexp = []
    all_sci_ID = fitstbl['sci_ID'].data[fitstbl['science']]  # Binary system: 1,2,4,8, etc.
    numsci = len(all_sci_ID)
    basenames = [None]*numsci  # For fluxing at the very end

    # Init calib dict
    calib_dict = {}
    caliBrate = calibrations.MultiSlitCalibrations(fitstbl)


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
            # Setup
            dnum = settings.get_dnum(det)
            msgs.info("Working on detector {:s}".format(dnum))
            sci_dict[det] = {}

            namp = settings.spect[dnum]["numamplifiers"]
            setup = arsetup.instr_setup(sci_ID, det, fitstbl, setup_dict, namp, must_exist=True)


            settings.argflag['reduce']['masters']['setup'] = setup

            ###############
            # Get data sections (Could avoid doing this for every sciexp, but it is quick)
            # TODO -_ Clean this up!
            scifile = os.path.join(fitstbl['directory'][scidx],fitstbl['filename'][scidx])
            settings_det = settings.spect[dnum].copy()  # Should include naxis0, naxis1 in this
            # Binning
            settings_det['binning'] = fitstbl['binning'][0]
            settings_det['dispaxis'] = settings.argflag['trace']['dispersion']['direction']
            # Yes, this looks goofy.  Is needed for LRIS and DEIMOS for now
            datasec, _, naxis0, naxis1 = io.get_datasec(spectrograph, scifile, det, settings_det)
            settings.spect[dnum]['naxis0'] = naxis0
            settings.spect[dnum]['naxis1'] = naxis1
            # Build the datasec_img
            datasec_img = arpixels.pix_to_amp(naxis0, naxis1, datasec, settings_det['numamplifiers'])
            settings.spect[dnum] = settings_det.copy()  # Used internally..

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

            # New ones
            ts_settings = dict(trace=settings.argflag['trace'], masters=settings.argflag['reduce']['masters'])
            ts_settings['masters']['directory'] = settings.argflag['run']['directory']['master']+'_'+ settings.argflag['run']['spectrograph']
            tsettings['trace'] = ts_settings['trace'].copy()
            tsettings['masters'] = ts_settings['masters'].copy()

            ###############################################################################
            # Begin calibrations
            caliBrate.reset(setup, det, sci_ID, tsettings, datasec_img)


            msbias = caliBrate.get_bias() # Bias frame or command
            msarc = caliBrate.get_arc() # Arc Image
            msbpm = caliBrate.get_bpm() # Bad pixel mask
            pixlocn = caliBrate.get_pixlocn()  # Physical pixel locations on the detector
            tslits_dict, maskslits = caliBrate.get_slits() # Slit Tracing
            wv_calib, maskslits = caliBrate.get_wv_calib() # Generate the 1D wavelength solution
            mstilts, maskslits = caliBrate.get_tilts() # Derive the spectral tilt
            mspixflatnrm, slitprof = caliBrate.get_pixflatnrm() # Prepare the pixel flat field frame
            mswave = caliBrate.get_wave() # Generate/load a master wave frame

            ''' Could also be run with 
            caliBrate.run_the_steps()
            '''

            # CALIBS END HERE
            ###############################################################################


            ###############
            #  Process and extract the science frame
            msgs.info("Working on the science frame")
            sci_image_files = arsort.list_of_files(fitstbl, 'science', sci_ID)
            # Settings
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
            if basenames[sc] is None:
                basenames[sc] = basename

            # Process (includes Variance image and CRs)
            dnoise = (settings_det['darkcurr'] * float(fitstbl["exptime"][scidx])/3600.0)
            sciframe, rawvarframe, crmask = sciI._process(
                msbias, mspixflatnrm, apply_gain=True, dnoise=dnoise)

            # Global skysub
            settings_skysub = {}
            settings_skysub['skysub'] = settings.argflag['reduce']['skysub'].copy()
            if settings.argflag['reduce']['skysub']['perform']:
                global_sky, modelvarframe = sciI.global_skysub(settings_skysub)
            else:
                sciI.global_sky = np.zeros_like(sciframe)
                sciI.modelvarframe = np.zeros_like(sciframe)

            # Find objects
            _, nobj = sciI.find_objects()
            if nobj == 0:
                msgs.warn("No objects to extract for science frame" + msgs.newline() + fitstbl['filename'][scidx])
                specobjs, flg_objs = [], None
            else:
                flg_objs = True  # Objects were found

            # Another round of sky sub
            if settings.argflag['reduce']['skysub']['perform'] and flg_objs:
                global_sky, modelvarframe = sciI.global_skysub(settings_skysub,
                                                               use_tracemask=True)
            # Another round of finding objects
            if flg_objs:
                _, nobj = sciI.find_objects()
                if nobj == 0:
                    msgs.warn("No objects to extract for science frame" + msgs.newline() + fitstbl['filename'][scidx])
                    specobjs, flg_objs = [], None

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
                    # QA
                    arwave.flexure_qa(specobjs, maskslits, basename, det, flex_list)

            # Helio
            # Correct Earth's motion
            vel_corr = -999999.9
            if (settings.argflag['reduce']['calibrate']['refframe'] in ['heliocentric', 'barycentric']) and \
                    (settings.argflag['reduce']['calibrate']['wavelength'] != "pixel") and flg_objs:
                if settings.argflag['science']['extraction']['reuse']:
                    msgs.warn("{0:s} correction will not be applied if an extracted science frame exists, and is used".format(settings.argflag['reduce']['calibrate']['refframe']))
                if specobjs is not None:
                    msgs.info("Performing a {0:s} correction".format(settings.argflag['reduce']['calibrate']['refframe']))
                    settings_mosaic = {}  # For long, lat of observatory
                    settings_mosaic['mosaic'] = settings.spect['mosaic'].copy()
                    vel, vel_corr = arwave.geomotion_correct(specobjs, maskslits, fitstbl, scidx,
                                                             obstime, settings_mosaic,
                                                             settings.argflag['reduce']['calibrate']['refframe'])
                else:
                    msgs.info("There are no objects on detector {0:d} to perform a {1:s} correction".format(
                        det, settings.argflag['reduce']['calibrate']['refframe']))
            else:
                msgs.info("A heliocentric correction will not be performed")

            # Save for outputing (after all detectors are done)
            sci_dict[det]['sciframe'] = sciframe
            if vel_corr > -999999.9:
                sci_dict['meta']['vel_corr'] = vel_corr
            if flg_objs:
                sci_dict[det]['specobjs'] = arutils.unravel_specobjs([specobjs])
                sci_dict[det]['finalvar'] = finalvar
                sci_dict[det]['finalsky'] = finalsky
            else:  # Nothing extracted
                sci_dict[det]['specobjs'] = []
                sci_dict[det]['finalvar'] = sciI.modelvarframe
                sci_dict[det]['finalsky'] = sciI.global_sky


            ######################################################
            # Reduce standard here; only legit if the mask is the same
            std_idx = arsort.ftype_indices(fitstbl, 'standard', sci_ID)
            if len(std_idx) > 0:
                std_idx = std_idx[0]
            else:
                msgs.info("No standard star associated with this science frame")
                continue
            #
            msgs.info("Processing standard star")
            std_image_files = arsort.list_of_files(fitstbl, 'standard', sci_ID)
            if std_idx in std_dict.keys():
                if det in std_dict[std_idx].keys():
                    continue
            else:
                std_dict[std_idx] = {}

            # Instantiate for the Standard
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
            _ = stdI.global_skysub(settings_skysub)
            # Find objects
            _, nobj = stdI.find_objects()
            if nobj == 0:
                msgs.warn("No objects to extract for standard frame" + msgs.newline() + fitstbl['filename'][scidx])
                continue
            _ = stdI.global_skysub(settings_skysub, use_tracemask=True)
            # Extract
            stdobjs, _, _ = stdI.extraction(mswave)
            # Save for fluxing and output later
            std_dict[std_idx][det] = {}
            std_dict[std_idx][det]['basename'] = std_basename
            std_dict[std_idx][det]['specobjs'] = arutils.unravel_specobjs([stdobjs])

        ###########################
        # Write
        # Build the final list of specobjs and vel_corr
        all_specobjs = []
        for key in sci_dict:
            if key in ['meta']:
                continue
            #
            all_specobjs += sci_dict[key]['specobjs']

        # Write 1D spectra
        save_format = 'fits'
        if save_format == 'fits':
            outfile = settings.argflag['run']['directory']['science']+'/spec1d_{:s}.fits'.format(basename)
            helio_dict = dict(refframe=settings.argflag['reduce']['calibrate']['refframe'],
                              vel_correction=sci_dict['meta']['vel_corr'])
            arsave.save_1d_spectra_fits(all_specobjs, fitstbl[scidx], outfile,
                                            helio_dict=helio_dict, obs_dict=settings.spect['mosaic'])
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
            sci_dict, fitstbl, scidx,
            settings.spect['fits']['headext{0:02d}'.format(1)], setup,
            settings.argflag['run']['directory']['master']+'_'+spectrograph, # MFDIR
            settings.argflag['run']['directory']['science'], basename)

    # Write standard stars at the very end
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
