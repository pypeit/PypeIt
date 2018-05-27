""" Primary module for guiding the reduction of long/multi-slit data
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import yaml
import numpy as np

from astropy import units

from pypit import msgs
from pypit import arparse as settings
from pypit import arload
from pypit import armasters
from pypit import armbase
from pypit import arproc
from pypit.core import arprocimg
from pypit import arsave
from pypit import arsciexp
from pypit.core import arsetup
from pypit import arpixels
from pypit.core import arsort
from pypit import wavetilts
from pypit import artrace
from pypit import arcimage
from pypit import bpmimage
from pypit import biasframe
from pypit import fluxspec
from pypit import traceslits
from pypit import traceimage
from pypit import wavecalib

from pypit import ardebug as debugger


def ARMS(fitstbl, setup_dict, reuseMaster=False, reloadMaster=True, sciexp=None):
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
    if sciexp is None:
        sciexp = []
        all_sci_ID = fitstbl['sci_ID'].data[fitstbl['science']]  # Binary system: 1,2,4,8, etc.
        for sci_ID in all_sci_ID:
            sciexp.append(arsciexp.ScienceExposure(sci_ID, fitstbl, settings.argflag,
                                                   settings.spect, do_qa=True))
            std_idx = arsort.ftype_indices(fitstbl, 'standard', sci_ID)
            if (len(std_idx) > 0):
                if len(std_idx) > 1:
                    msgs.warn("Will only reduce the first, unique standard for each standard frame!")
                if std_idx[0] not in sv_std_idx:  # Only take the first one
                    sv_std_idx.append(std_idx[0])
                    # Standard stars
                    std_dict[std_idx[0]] = arsciexp.ScienceExposure(sci_ID, fitstbl, settings.argflag,
                                                                    settings.spect, do_qa=False, idx_sci=std_idx[0])
    numsci = len(sciexp)

    # Init calib dict
    calib_dict = {}

    # Loop on science exposure
    for sc in range(numsci):

        slf = sciexp[sc]
        sci_ID = slf.sci_ID
        scidx = slf._idx_sci[0]
        msgs.info("Reducing file {0:s}, target {1:s}".format(fitstbl['filename'][scidx], slf._target_name))
        msgs.sciexp = slf  # For QA writing on exit, if nothing else.  Could write Masters too

        #if reloadMaster and (sc > 0):
        #    settings.argflag['reduce']['masters']['reuse'] = True

        # Loop on Detectors
        for kk in range(settings.spect['mosaic']['ndet']):
            det = kk + 1  # Detectors indexed from 1
            if settings.argflag['reduce']['detnum'] is not None:
                if det not in map(int, settings.argflag['reduce']['detnum']):
                    msgs.warn("Skipping detector {:d}".format(det))
                    continue
                else:
                    msgs.warn("Restricting the reduction to detector {:d}".format(det))

            slf.det = det
            dnum = settings.get_dnum(det)
            msgs.info("Working on detector {:s}".format(dnum))

            # Setup
            namp = settings.spect[dnum]["numamplifiers"]
            setup = arsetup.instr_setup(sci_ID, det, fitstbl, setup_dict, namp, must_exist=True)
            settings.argflag['reduce']['masters']['setup'] = setup
            slf.setup = setup

            ###############
            # Get data sections (Could avoid doing this for every sciexp, but it is quick)
            datasec_img = arprocimg.get_datasec_trimmed(fitstbl, det, scidx, settings.argflag, settings.spect)
            #slf._datasec[det-1] = pix_to_amp(naxis0, naxis1, datasec, numamplifiers)

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

            ###############
            # Prepare for Bias subtraction
            #   bias will either be an image (ndarray) or a command (str, e.g. 'overscan') or none
            if 'bias' in calib_dict[setup].keys():
                msbias = calib_dict[setup]['bias']
            else:
                # Init
                bias = biasframe.BiasFrame(settings=tsettings, setup=setup, det=det, fitstbl=fitstbl, sci_ID=sci_ID)
                # Load the MasterFrame (if it exists and is desired) or the command (e.g. 'overscan')
                msbias = bias.master()
                if msbias is None:  # Build it and save it
                    msbias = bias.build_image()
                    bias.save_master(msbias, raw_files=bias.file_list, steps=bias.steps)
                # Save
                calib_dict[setup]['bias'] = msbias

            ###############
            # Generate a master arc frame
            if 'arc' in calib_dict[setup].keys():
                msarc = calib_dict[setup]['arc']
            else:
                # Instantiate with everything needed to generate the image (in case we do)
                AImage = arcimage.ArcImage([], spectrograph=settings.argflag['run']['spectrograph'],
                                           settings=tsettings, det=det, setup=setup, sci_ID=sci_ID,
                                           msbias=msbias, fitstbl=fitstbl)
                # Load the MasterFrame (if it exists and is desired)?
                msarc = AImage.master()
                if msarc is None:  # Otherwise build it
                    msgs.info("Preparing a master {0:s} frame".format(AImage.frametype))
                    msarc = AImage.build_image()
                    # Save to Masters
                    AImage.save_master(msarc, raw_files=AImage.file_list, steps=AImage.steps)
                # Save
                calib_dict[setup]['arc'] = msarc

            ###############
            # Generate a bad pixel mask (should not repeat)
            if 'bpm' in calib_dict[setup].keys():
                msbpm = calib_dict[setup]['bpm']
            else:
                bpmImage = bpmimage.BPMImage(spectrograph=settings.argflag['run']['spectrograph'],
                                             settings=tsettings, det=det,
                                             shape=msarc.shape,
                                             binning=fitstbl['binning'][scidx],
                                             reduce_badpix=settings.argflag['reduce']['badpix'],
                                             msbias=msbias)
                msbpm = bpmImage.build()
                # Save
                calib_dict[setup]['bpm'] = msbpm

            ###############
            # Generate an array that provides the physical pixel locations on the detector
            pixlocn = arpixels.gen_pixloc(msarc.shape, det, settings.argflag)
            # TODO -- Deprecate using slf for this
            slf.SetFrame(slf._pixlocn, pixlocn, det)

            ###############
            # Slit Tracing
            if 'trace' in calib_dict[setup].keys():  # Internal
                Tslits = calib_dict[setup]['trace']
            else:
                # Setup up the settings (will be Refactored with settings)
                tmp = dict(trace=settings.argflag['trace'], masters=settings.argflag['reduce']['masters'])
                tmp['masters']['directory'] = settings.argflag['run']['directory']['master']+'_'+ settings.argflag['run']['spectrograph']

                # Instantiate (without mstrace)
                Tslits = traceslits.TraceSlits(None, slf._pixlocn[det-1], settings=tmp, det=det, setup=setup, binbpx=msbpm)

                # Load via masters, as desired
                if not Tslits.master():
                    # Build the trace image first
                    trace_image_files = arsort.list_of_files(fitstbl, 'trace', sci_ID)
                    Timage = traceimage.TraceImage(trace_image_files,
                                                   spectrograph=settings.argflag['run']['spectrograph'],
                                                   settings=tsettings, det=det)
                    mstrace = Timage.process(bias_subtract=msbias, trim=settings.argflag['reduce']['trim'])

                    # Load up and get ready
                    Tslits.mstrace = mstrace
                    _ = Tslits.make_binarr()
                    # Now we go forth
                    Tslits.run(armlsd=True)#, ignore_orders=ignore_orders, add_user_slits=add_user_slits)
                    # QA
                    Tslits._qa()
                    # Save to disk
                    Tslits.save_master()

                # Save in calib
                calib_dict[setup]['trace'] = Tslits

            # Save in slf
            # TODO -- Deprecate this means of holding the info (e.g. just pass around Tslits)
            slf.SetFrame(slf._lordloc, Tslits.lcen, det)
            slf.SetFrame(slf._rordloc, Tslits.rcen, det)
            slf.SetFrame(slf._pixcen, Tslits.pixcen, det)
            slf.SetFrame(slf._pixwid, Tslits.pixwid, det)
            slf.SetFrame(slf._lordpix, Tslits.lordpix, det)
            slf.SetFrame(slf._rordpix, Tslits.rordpix, det)
            slf.SetFrame(slf._slitpix, Tslits.slitpix, det)

            # Initialize maskslit
            slf._maskslits[det-1] = np.zeros(slf._lordloc[det-1].shape[1], dtype=bool)

            ###############
            # Generate the 1D wavelength solution
            if 'wavecalib' in calib_dict[setup].keys():
                wv_calib = calib_dict[setup]['wavecalib']
                wv_maskslits = calib_dict[setup]['wvmask']
            elif settings.argflag["reduce"]["calibrate"]["wavelength"] == "pixel":
                msgs.info("A wavelength calibration will not be performed")
                pass
            else:
                # Setup up the settings (will be Refactored with settings)
                tmp = dict(calibrate=settings.argflag['arc']['calibrate'], masters=settings.argflag['reduce']['masters'])
                tmp['masters']['directory'] = settings.argflag['run']['directory']['master']+'_'+ settings.argflag['run']['spectrograph']

                # Instantiate
                Wavecalib = wavecalib.WaveCalib(msarc, spectrograph=settings.argflag['run']['spectrograph'],
                                                settings=tmp, det=det, setup=setup, fitstbl=fitstbl, sci_ID=sci_ID)
                # Load from disk (MasterFrame)?
                wv_calib = Wavecalib.master()
                # Build?
                if wv_calib is None:
                    nonlinear = settings.spect[settings.get_dnum(det)]['saturation'] * settings.spect[settings.get_dnum(det)]['nonlinear']
                    wv_calib, _ = Wavecalib.run(Tslits.lcen, Tslits.rcen, pixlocn, nonlinear=nonlinear)
                    # Save to Masters
                    Wavecalib.save_master(Wavecalib.wv_calib)
                else:
                    Wavecalib.wv_calib = wv_calib
                # Mask
                wv_maskslits = Wavecalib._make_maskslits(Tslits.lcen.shape[1])

                # Save in calib
                calib_dict[setup]['wavecalib'] = wv_calib
                calib_dict[setup]['wvmask'] = wv_maskslits

            # Mask me
            slf._maskslits[det-1] += wv_maskslits

            ###############
            # Derive the spectral tilt
            if 'tilts' in calib_dict[setup].keys():
                mstilts = calib_dict[setup]['tilts']
                wt_maskslits = calib_dict[setup]['wtmask']
            else:
                # Settings kludges
                tilt_settings = dict(tilts=settings.argflag['trace']['slits']['tilts'].copy())
                tilt_settings['tilts']['function'] = settings.argflag['trace']['slits']['function']
                tilt_settings['masters'] = settings.argflag['reduce']['masters']
                tilt_settings['masters']['directory'] = settings.argflag['run']['directory']['master']+'_'+ settings.argflag['run']['spectrograph']
                settings_det = {}
                settings_det[dnum] = settings.spect[dnum].copy()
                # Instantiate
                wTilt = wavetilts.WaveTilts(msarc, settings=tilt_settings, det=det, setup=setup,
                                            lordloc=Tslits.lcen, rordloc=Tslits.rcen,
                                            pixlocn=Tslits.pixlocn, pixcen=Tslits.pixcen,
                                            slitpix=Tslits.slitpix, settings_det=settings_det)
                # Master
                mstilts = wTilt.master()
                if mstilts is None:
                    mstilts, wt_maskslits = wTilt.run(maskslits=slf._maskslits[det-1])
                    wTilt.save_master()
                else:
                    wt_maskslits = np.zeros(len(slf._maskslits[det-1]), dtype=bool)
                # Save
                calib_dict[setup]['tilts'] = mstilts
                calib_dict[setup]['wtmask'] = wt_maskslits
            slf._maskslits[det-1] += wt_maskslits


            ###############
            # Prepare the pixel flat field frame
            update = slf.MasterFlatField(fitstbl, det, msbias, datasec_img, mstilts)
            if update and reuseMaster: armbase.UpdateMasters(sciexp, sc, det, ftype="flat", chktype="pixelflat")

            ###############
            # Generate/load a master wave frame
            update = slf.MasterWave(det, wv_calib, mstilts)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="arc", chktype="wave")

            ###############
            # Load the science frame and from this generate a Poisson error frame
            msgs.info("Loading science frame")
            sciframe = arload.load_frames(fitstbl, [scidx], det,
                                          frametype='science',
                                          msbias=msbias) # slf._msbias[det-1])
            sciframe = sciframe[:, :, 0]
            # Extract
            msgs.info("Processing science frame")
            arproc.reduce_multislit(slf, mstilts, sciframe, msbpm, datasec_img, scidx, fitstbl, det)

            ######################################################
            # Reduce standard here; only legit todo if the mask is the same
            std_idx = arsort.ftype_indices(fitstbl, 'standard', sci_ID)
            if len(std_idx) > 0:
                std_idx = std_idx[0]
            else:
                continue
            stdslf = std_dict[std_idx]
            if stdslf.extracted is False:
                # Fill up the necessary pieces
                for iattr in ['pixlocn', 'lordloc', 'rordloc', 'pixcen', 'pixwid', 'lordpix', 'rordpix',
                              'slitpix', 'satmask', 'maskslits', 'slitprof',
                              'mspixelflatnrm', 'mswave']:
                    setattr(stdslf, '_'+iattr, getattr(slf, '_'+iattr))  # Brings along all the detectors, but that is ok
                # Load
                stdframe = arload.load_frames(fitstbl, [std_idx], det, frametype='standard', msbias=msbias)
                stdframe = stdframe[:, :, 0]
                # Reduce
                msgs.info("Processing standard frame")
                arproc.reduce_multislit(stdslf, mstilts, stdframe, msbpm, datasec_img, std_idx, fitstbl, det, standard=True)
                # Finish
                stdslf.extracted = True

    # Write standard stars
    for key in std_dict.keys():
        outfile = settings.argflag['run']['directory']['science']+'/spec1d_{:s}.fits'.format(std_dict[key]._basename)
        arsave.new_save_1d_spectra_fits(std_dict[key]._specobjs, fitstbl[std_idx], outfile,
                                        obs_dict=settings.spect['mosaic'])

    #########################
    # Flux towards the very end..
    #########################
    if (settings.argflag['reduce']['calibrate']['flux'] == True) and (len(std_dict) > 0):
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
        if (settings.argflag['reduce']['calibrate']['sensfunc']['archival'] == 'None'):
            std_keys = list(std_dict.keys())
            std_key = std_keys[0] # Take the first extraction
            FxSpec = fluxspec.FluxSpec(settings=fsettings, std_specobjs=std_dict[std_key]._specobjs,
                                       setup=setup)  # This takes the last setup run, which is as sensible as any..
            sensfunc = FxSpec.master(fitstbl[std_key])
        else:  # Input by user
            FxSpec = fluxspec.FluxSpec(settings=fsettings,
                                       sens_file=settings.argflag['reduce']['calibrate']['sensfunc']['archival'])
            sensfunc = FxSpec.sensfunc
        # Flux
        msgs.info("Fluxing with {:s}".format(sensfunc['std']['name']))
        for slf in sciexp:
            scidx = slf._idx_sci[0]
            FxSpec._flux_specobjs(slf._specobjs, fitstbl['airmass'][scidx], fitstbl['exptime'][scidx])


    # Write science
    for sc in range(numsci):
        slf = sciexp[sc]

        # TODO -- When we refactor ScienceExposure, something will need to carry all the individual images until we write them out..
        # Write 1D spectra
        save_format = 'fits'
        if save_format == 'fits':
            outfile = settings.argflag['run']['directory']['science']+'/spec1d_{:s}.fits'.format(slf._basename)
            helio_dict = dict(refframe=settings.argflag['reduce']['calibrate']['refframe'],
                              vel_correction=slf.vel_correction)
            arsave.new_save_1d_spectra_fits(slf._specobjs, fitstbl[slf._idx_sci[0]], outfile,
                                            helio_dict=helio_dict, obs_dict=settings.spect['mosaic'])
            #arsave.save_1d_spectra_fits(slf, fitstbl)
        elif save_format == 'hdf5':
            arsave.save_1d_spectra_hdf5(slf)
        else:
            msgs.error(save_format + ' is not a recognized output format!')
        arsave.save_obj_info(slf, fitstbl)
        # Write 2D images for the Science Frame
        arsave.save_2d_images(slf, fitstbl)
        # Free up some memory by replacing the reduced ScienceExposure class
        #sciexp[sc] = None

    return status
