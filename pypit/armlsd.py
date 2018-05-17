""" Primary module for guiding the reduction of long slit data
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import yaml
import numpy as np

from astropy import units

from pypit import msgs
from pypit import arparse as settings
from pypit import arflux
from pypit import arload
from pypit import armasters
from pypit import armbase
from pypit import arproc
from pypit import arsave
from pypit import arsciexp
from pypit import arsetup
from pypit import ardeimos
from pypit import arpixels
from pypit import arsort
from pypit import artrace
from pypit import biasframe
from pypit.core import artraceslits
from pypit import traceslits
from pypit import traceimage

from pypit import ardebug as debugger


def ARMLSD(fitstbl, setup_dict, reuseMaster=False, reloadMaster=True, sciexp=None, original=True):
    """
    Automatic Reduction and Modeling of Long Slit Data

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

    '''  THIS IS NOW RUN IN pypit.py
    # Create a list of science exposure classes
    original=True
    if original:
        mode, sciexp, setup_dict = armbase.setup_science(fitsdict)
        if mode == 'setup':
            status = 1
            return status
        elif mode == 'calcheck':
            status = 2
            return status
        else:
            numsci = len(sciexp)
    else:
        # This should move inside of PYPIT
        setupc = setupclass.SetupClass(settings, fitstbl=fitstbl)
        mode, fitstbl, setup_dict = setupc.run()
    '''
    # Generate sciexp list, if need be (it will be soon)
    if sciexp is None:
        sciexp = []
        all_sci_ID = fitstbl['sci_ID'].data[fitstbl['science']]  # Binary system: 1,2,4,8, etc.
        for ii in all_sci_ID:
            sciexp.append(arsciexp.ScienceExposure(ii, fitstbl, settings.argflag,
                                                   settings.spect, do_qa=True, original=original))
    numsci = len(sciexp)

    # Start reducing the data

    # Loop on Detectors
    for kk in range(settings.spect['mosaic']['ndet']):
        det = kk + 1  # Detectors indexed from 1
        if settings.argflag['reduce']['detnum'] is not None:
            if det not in map(int, settings.argflag['reduce']['detnum']):
                msgs.warn("Skipping detector {:d}".format(det))
                continue
            else:
                msgs.warn("Restricting the reduction to detector {:d}".format(det))

        # Reset Calib dict
        calib_dict = {}

        for sc in range(numsci):

            slf = sciexp[sc]
            sci_ID = slf.sci_ID
            scidx = slf._idx_sci[0]
            msgs.info("Reducing file {0:s}, target {1:s}".format(fitstbl['filename'][scidx], slf._target_name))
            msgs.sciexp = slf  # For QA writing on exit, if nothing else.  Could write Masters too

            #if reloadMaster and (sc > 0):
            #    settings.argflag['reduce']['masters']['reuse'] = True

            slf.det = det
            dnum = settings.get_dnum(det)
            ###############
            # Get data sections
            arproc.get_datasec_trimmed(slf, fitstbl, det, scidx)
            # Setup
            if original:
                setup = arsetup.instr_setup(slf, det, fitstbl, setup_dict, must_exist=True)
            else:
                namp = settings.spect[dnum]["numamplifiers"]
                setup = arsetup.new_instr_setup(sci_ID, det, fitstbl, setup_dict, namp, must_exist=True)
            settings.argflag['reduce']['masters']['setup'] = setup
            slf.setup = setup

            # Allow for variants in the setup for each Science frame (e.g. Arc pairings)
            calib_dict[setup] = {}

            # TODO -- Update/avoid the following with new settings
            tsettings = settings.argflag.copy()
            tsettings['detector'] = settings.spect[settings.get_dnum(det)]
            try:
                tsettings['detector']['dataext'] = tsettings['detector']['dataext01']  # Kludge; goofy named key
            except KeyError: # LRIS, DEIMOS
                pass
            tsettings['detector']['dispaxis'] = settings.argflag['trace']['dispersion']['direction']

            ###############
            # Prepare for Bias subtraction
            #   bias will either be an image (ndarray) or a command (str, e.g. 'overscan') or none
            if 'bias' in calib_dict.keys():
                msbias = calib_dict[setup]['bias']
            else:
                # Init
                bias = biasframe.BiasFrame(settings=tsettings, setup=setup, ind=slf._idx_bias, det=det, fitstbl=fitstbl)
                # If an image is generated, it will be saved to disk a a MasterFrame
                msbias = bias.build_master()
                # Save in Calib dict -- Will replace Master class
                calib_dict[setup]['bias'] = msbias

            # Set in slf -- Will be Deprecated
            if isinstance(bias, np.ndarray):
                slf.SetMasterFrame(bias, "bias", det)

            #update = slf.MasterBias(fitsdict, det)
            #if update and reuseMaster:
            #    armbase.UpdateMasters(sciexp, sc, det, ftype="bias")

            ###############
            # Generate a bad pixel mask (should not repeat)
            update = slf.BadPixelMask(fitstbl, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="arc")
            ###############
            # Generate a master arc frame
            update = slf.MasterArc(fitstbl, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="arc")

            ###############
            # Set the number of spectral and spatial pixels, and the bad pixel mask is it does not exist
            slf._nspec[det-1], slf._nspat[det-1] = slf._msarc[det-1].shape
            if slf._bpix[det-1] is None:
                bpix = np.zeros((slf._nspec[det-1], slf._nspat[det-1]))
                if settings.argflag['run']['spectrograph'] in ['keck_deimos']:
                    bpix = ardeimos.bpm(det)
                slf.SetFrame(slf._bpix, bpix, det)

            ###############
            # Generate an array that provides the physical pixel locations on the detector
            pixlocn = arpixels.gen_pixloc(slf._msarc[det-1].shape, det, settings.argflag)
            slf.SetFrame(slf._pixlocn, pixlocn, det)

            '''
            ###############
            # Estimate gain and readout noise for the amplifiers
            msgs.work("Estimate Gain and Readout noise from the raw frames...")
            update = slf.MasterRN(fitsdict, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="readnoise")
            '''
            ###############
            # Slit Tracing
            # Load master frame?
            if (settings.argflag['reduce']['masters']['reuse']) or (settings.argflag['reduce']['masters']['force']):
                mstrace_name = armasters.master_name('trace', setup)
                Tslits = traceslits.TraceSlits.from_master_files(mstrace_name, load_pix_obj=True)  # Returns None if none exists
            else:
                Tslits = None
            # Build it?  Beginning with the trace image
            if Tslits is None:
                # Build the trace image first
                trace_image_files = arsort.list_of_files(fitstbl, 'trace', sci_ID)
                Timage = traceimage.TraceImage(trace_image_files,
                                                    spectrograph=settings.argflag['run']['spectrograph'],
                                                    settings=tsettings, det=det)
                mstrace = Timage.process(bias_subtract=msbias, trim=settings.argflag['reduce']['trim'])
                # Setup up the settings (will be Refactored with settings)
                tmp = dict(trace=settings.argflag['trace'], masters=settings.argflag['reduce']['masters'])
                tmp['masters']['directory'] = settings.argflag['run']['directory']['master']+'_'+ settings.argflag['run']['spectrograph']
                # Now trace me some slits
                Tslits = traceslits.TraceSlits(mstrace, slf._pixlocn[det-1], det=det, settings=tmp,
                                               binbpx=slf._bpix[det-1], setup=setup)
                _ = Tslits.run(armlsd=True)
                # QA
                Tslits._qa()
                # Save to disk
                Tslits.save_master()

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

            # Save mstrace (should not be necessary)
            #armbase.UpdateMasters(sciexp, sc, det, ftype="flat", chktype="trace")

            ###############
            # Generate the 1D wavelength solution
            update = slf.MasterWaveCalib(fitstbl, sc, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="arc", chktype="trace")

            ###############
            # Derive the spectral tilt
            if slf._tilts[det-1] is None:
                tilts = armasters.load_master_frame(slf, "tilts")
                if tilts is None:
                    # First time tilts are derived for this arc frame --> derive the order tilts
                    tilts, satmask, outpar = artrace.multislit_tilt(slf, slf._msarc[det-1], det)
                    slf.SetFrame(slf._tilts, tilts, det)
                    slf.SetFrame(slf._satmask, satmask, det)
                    msgs.bug("This outpar is only the last slit!!  JXP doesn't think it matters for now")
                    slf.SetFrame(slf._tiltpar, outpar, det)
                    armasters.save_masters(slf, det, mftype='tilts')
                else:
                    slf.SetFrame(slf._tilts, tilts, det)

            ###############
            # Prepare the pixel flat field frame
            update = slf.MasterFlatField(fitstbl, det)
            if update and reuseMaster: armbase.UpdateMasters(sciexp, sc, det, ftype="flat", chktype="pixelflat")

            ###############
            # Generate/load a master wave frame
            update = slf.MasterWave(fitstbl, sc, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="arc", chktype="wave")

            ###############
            # Check if the user only wants to prepare the calibrations only
            msgs.info("All calibration frames have been prepared")
            if settings.argflag['run']['preponly']:
                msgs.info("If you would like to continue with the reduction, disable the command:" + msgs.newline() +
                          "run preponly False")
                continue

            ###############
            # Write setup
            #setup = arsort.calib_setup(sc, det, fitsdict, setup_dict, write=True)
            # Write MasterFrames (currently per detector)
            #armasters.save_masters(slf, det, setup)

            ###############
            # Load the science frame and from this generate a Poisson error frame
            msgs.info("Loading science frame")
            sciframe = arload.load_frames(fitstbl, [scidx], det,
                                          frametype='science',
                                          msbias=slf._msbias[det-1])
            sciframe = sciframe[:, :, 0]
            # Extract
            msgs.info("Processing science frame")
            arproc.reduce_multislit(slf, sciframe, scidx, fitstbl, det)

            ###############
            # Using model sky, calculate a flexure correction


        # Write 1D spectra
        save_format = 'fits'
        if save_format == 'fits':
            arsave.save_1d_spectra_fits(slf, fitstbl)
        elif save_format == 'hdf5':
            arsave.save_1d_spectra_hdf5(slf)
        else:
            msgs.error(save_format + ' is not a recognized output format!')
        arsave.save_obj_info(slf, fitstbl)
        # Write 2D images for the Science Frame
        arsave.save_2d_images(slf, fitstbl)
        # Free up some memory by replacing the reduced ScienceExposure class
        sciexp[sc] = None

    debugger.set_trace()
    #########################
    # Flux at the very end..
    #########################
    if (settings.argflag['reduce']['calibrate']['flux'] == True):
        # Standard star (is this a calibration, e.g. goes above?)
        msgs.info("Processing standard star")
        msgs.info("Assuming one star per detector mosaic")
        msgs.info("Waited until last detector to process")

        if (settings.argflag['reduce']['calibrate']['sensfunc']['archival'] == 'None'):
            update = slf.MasterStandard(fitstbl)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, 0, ftype="standard")
        else:
            sensfunc = yaml.load(open(settings.argflag['reduce']['calibrate']['sensfunc']['archival']))
            # Yaml does not do quantities, so make the sensfunc min/max wave quantities
            sensfunc['wave_max'] *= units.angstrom
            sensfunc['wave_min'] *= units.angstrom
            slf.SetMasterFrame(sensfunc, "sensfunc", None, mkcopy=False)
            msgs.info(
                "Using archival sensfunc {:s}".format(settings.argflag['reduce']['calibrate']['sensfunc']['archival']))

        msgs.info("Fluxing with {:s}".format(slf._sensfunc['std']['name']))
        for kk in range(settings.spect['mosaic']['ndet']):
            det = kk + 1  # Detectors indexed from 1
            if slf._specobjs[det - 1] is not None:
                arflux.apply_sensfunc(slf, det, scidx, fitstbl)
            else:
                msgs.info("There are no objects on detector {0:d} to apply a flux calibration".format(det))

    return status
