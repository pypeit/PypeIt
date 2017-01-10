""" Primary module for guiding the reduction of echelle data
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
from pypit import arparse as settings
from pypit import arload
from pypit import armasters
from pypit import armbase
from pypit import armsgs
from pypit import arproc
from pypit import arsave
from pypit import arsort
from pypit import artrace
from pypit import arqa

try:
    from xastropy.xutils import xdebug as debugger
except ImportError:
    import pdb as debugger

# Logging
msgs = armsgs.get_logger()


def ARMED(fitsdict, reuseMaster=False, reloadMaster=True):
    """
    Automatic Reduction and Modeling of Echelle Data

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
      0 = Successful execution
      1 = ...
    """
    status = 0

    # Create a list of science exposure classes
    sciexp = armbase.SetupScience(fitsdict)
    numsci = len(sciexp)

    # Create a list of master calibration frames
    masters = armasters.MasterFrames(settings.spect['mosaic']['ndet'])

    # Use Masters?  Requires setup file
    setup_file, nexist = arsort.get_setup_file()
    setup_dict = arsort.load_setup()
    settings.argflag['reduce']['masters']['file'] = setup_file

    # Start reducing the data
    for sc in range(numsci):
        slf = sciexp[sc]
        scidx = slf._idx_sci[0]
        msgs.info("Reducing file {0:s}, target {1:s}".format(fitsdict['filename'][scidx], slf._target_name))
        msgs.sciexp = slf  # For QA writing on exit, if nothing else.  Could write Masters too
        if reloadMaster and (sc > 0):
            settings.argflag['reduce']['masters']['reuse'] = True
        # Loop on Detectors
        for kk in range(settings.spect['mosaic']['ndet']):
            det = kk + 1  # Detectors indexed from 1
            slf.det = det
            ###############
            # Get datasec sections
            arproc.get_datasec_trimmed(slf, fitsdict, det, scidx)
            # Setup
            setup = arsort.instr_setup(sc, det, fitsdict, setup_dict, must_exist=True)
            settings.argflag['reduce']['masters']['setup'] = setup
            ###############
            # Generate master bias frame
            update = slf.MasterBias(fitsdict, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="bias")
            ###############
            # Generate a bad pixel mask (should not repeat)
            update = slf.BadPixelMask(fitsdict, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="arc")
            ###############
            # Estimate gain and readout noise for the amplifiers
            msgs.work("Estimate Gain and Readout noise from the raw frames...")
            ###############
            # Generate a master arc frame
            update = slf.MasterArc(fitsdict, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="arc")
            ###############
            # Set the number of spectral and spatial pixels, and the bad pixel mask is it does not exist
            slf._nspec[det-1], slf._nspat[det-1] = slf._msarc[det-1].shape
            if slf._bpix[det-1] is None:
                slf.SetFrame(slf._bpix, np.zeros((slf._nspec[det-1], slf._nspat[det-1])), det)
            ###############
            # Generate a master trace frame
            update = slf.MasterTrace(fitsdict, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="flat", chktype="trace")
            ###############
            # Generate a master pinhole frame
            update = slf.MasterPinhole(fitsdict, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="flat", chktype="pinhole")
            ###############
            # Generate an array that provides the physical pixel locations on the detector
            slf.GetPixelLocations(det)
            ###############
            # Determine the edges of the spectrum (spatial)
            if ('trace'+settings.argflag['reduce']['masters']['setup'] not in settings.argflag['reduce']['masters']['loaded']):
                ###############
                # Determine the centroid of the spectrum (spatial)
                lordloc, rordloc, extord = artrace.trace_slits(slf, slf._mspinhole[det-1], det,
                                                               pcadesc="PCA trace of the slit edges")

                # Using the order centroid, expand the order edges until the edge of the science slit is found
                if settings.argflag['trace']['slits']['expand']:
                    lordloc, rordloc = artrace.expand_slits(slf, slf._mstrace[det-1], det,
                                                            0.5*(lordloc+rordloc), extord)

                # Save the locations of the order edges
                slf.SetFrame(slf._lordloc, lordloc, det)
                slf.SetFrame(slf._rordloc, rordloc, det)

                # Convert physical trace into a pixel trace
                msgs.info("Converting physical trace locations to nearest pixel")
                pixcen = artrace.phys_to_pix(0.5 * (slf._lordloc[det - 1] + slf._rordloc[det - 1]), slf._pixlocn[det - 1], 1)
                pixwid = (slf._rordloc[det - 1] - slf._lordloc[det - 1]).mean(0).astype(np.int)
                lordpix = artrace.phys_to_pix(slf._lordloc[det - 1], slf._pixlocn[det - 1], 1)
                rordpix = artrace.phys_to_pix(slf._rordloc[det - 1], slf._pixlocn[det - 1], 1)
                slf.SetFrame(slf._pixcen, pixcen, det)
                slf.SetFrame(slf._pixwid, pixwid, det)
                slf.SetFrame(slf._lordpix, lordpix, det)
                slf.SetFrame(slf._rordpix, rordpix, det)
                # Save QA for slit traces
                arqa.slit_trace_qa(slf, slf._mstrace[det - 1], slf._lordpix[det - 1], slf._rordpix[det - 1], extord,
                                   desc="Trace of the slit edges", normalize=False)
                armbase.UpdateMasters(sciexp, sc, det, ftype="flat", chktype="trace")

            ###############
            # Generate the 1D wavelength solution
            update = slf.MasterWaveCalib(fitsdict, sc, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="arc", chktype="trace")
            ###############
            # Derive the spectral tilt
            if slf._tilts[det-1] is None:
                if settings.argflag['reduce']['masters']['reuse']:
                    mstilt_name = armasters.master_name(settings.argflag['run']['directory']['master'],
                                                        'tilts', settings.argflag['reduce']['masters']['setup'])
                    try:
                        tilts, head = arload.load_master(mstilt_name, frametype="tilts")
                    except IOError:
                        pass
                    else:
                        slf.SetFrame(slf._tilts, tilts, det)
                        settings.argflag['reduce']['masters']['loaded'].append('tilts'+settings.argflag['reduce']['masters']['setup'])
                if 'tilts'+settings.argflag['reduce']['masters']['setup'] not in settings.argflag['reduce']['masters']['loaded']:
                    # First time tilts are derived for this arc frame --> derive the order tilts
                    tilts, satmask, outpar = artrace.echelle_tilt(slf, slf._msarc[det-1], det)
                    slf.SetFrame(slf._tilts, tilts, det)
                    slf.SetFrame(slf._satmask, satmask, det)
                    slf.SetFrame(slf._tiltpar, outpar, det)

            ###############
            # Prepare the pixel flat field frame
            update = slf.MasterFlatField(fitsdict, det)
            if update and reuseMaster: armbase.UpdateMasters(sciexp, sc, det, ftype="flat", chktype="pixelflat")

            msgs.error("UP TO HERE!!!")
            ###############
            # Generate/load a master wave frame
            update = slf.MasterWave(fitsdict, sc, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="arc", chktype="wave")

            # Check if the user only wants to prepare the calibrations only
            msgs.info("All calibration frames have been prepared")
            if settings.argflag['run']['preponly']:
                msgs.info("If you would like to continue with the reduction,"
                          +msgs.newline()+"disable the run+preponly command")
                continue

            # Write setup
            #setup = arsort.calib_setup(sc, det, fitsdict, setup_dict, write=True)
            # Write MasterFrames (currently per detector)
            armasters.save_masters(slf, det, setup)

            ###############
            # Load the science frame and from this generate a Poisson error frame
            msgs.info("Loading science frame")
            sciframe = arload.load_frames(fitsdict, [scidx], det,
                                          frametype='science',
                                          msbias=slf._msbias[det - 1])
            sciframe = sciframe[:, :, 0]
            # Extract
            msgs.info("Processing science frame")
            arproc.reduce_frame(slf, sciframe, scidx, fitsdict, det)

        # Close the QA for this object
        slf._qa.close()

        # Write 1D spectra
        save_format = 'fits'
        if save_format == 'fits':
            arsave.save_1d_spectra_fits(slf)
        elif save_format == 'hdf5':
            arsave.save_1d_spectra_hdf5(slf)
        else:
            msgs.error(save_format + ' is not a recognized output format!')
        # Write 2D images for the Science Frame
        arsave.save_2d_images(slf)
        # Free up some memory by replacing the reduced ScienceExposure class
        sciexp[sc] = None
    return status
