import pdb
import numpy as np
import arextract
import arflux
import ario
import arload
import armasters
import armbase
import armsgs
import arproc
import ararc
import arspecobj
import artrace
import arqa

try:
    from xastropy.xutils import xdebug as xdb
except:
    pass

# Logging
msgs = armsgs.get_logger()

def ARMLSD(argflag, spect, fitsdict, reuseMaster=False):
    """
    Automatic Reduction and Modeling of Long Slit Data

    Parameters
    ----------
    argflag : dict
      Arguments and flags used for reduction
    spect : dict
      Properties of the spectrograph.
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
    sciexp = armbase.SetupScience(argflag, spect, fitsdict)
    numsci = len(sciexp)

    # Create a list of master calibration frames
    masters = armasters.MasterFrames(spect['mosaic']['ndet'])

    # Start reducing the data
    for sc in range(numsci):
        slf = sciexp[sc]
        scidx = slf._idx_sci[0]
        msgs.info("Reducing file {0:s}, target {1:s}".format(fitsdict['filename'][scidx], slf._target_name))
        # Loop on Detectors
        for kk in xrange(slf._spect['mosaic']['ndet']):
            det = kk + 1  # Detectors indexed from 1
            ###############
            # Get amplifier sections
            fitsdict = arproc.get_ampsec_trimmed(slf, fitsdict, det, scidx)
            ###############
            # Generate master bias frame
            update = slf.MasterBias(fitsdict, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="bias")
            ###############
            # Generate a bad pixel mask (should not repeat)
            update = slf.BadPixelMask(det)
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
            # Determine the dispersion direction (and transpose if necessary)
            slf.GetDispersionDirection(fitsdict, det, scidx)
            if slf._bpix[det-1] is None:
                slf.SetFrame(slf._bpix, np.zeros((slf._nspec[det-1], slf._nspat[det-1])), det)
            ###############
            # Generate a master trace frame
            update = slf.MasterTrace(fitsdict, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="flat", chktype="trace")
            ###############
            # Generate an array that provides the physical pixel locations on the detector
            slf.GetPixelLocations(det)
            ###############
            # Determine the edges of the spectrum (spatial)
            lordloc, rordloc, extord = artrace.trace_orders(slf, slf._mstrace[det-1], det, singleSlit=True, pcadesc="PCA trace of the slit edges")
            slf.SetFrame(slf._lordloc, lordloc, det)
            slf.SetFrame(slf._rordloc, rordloc, det)
            # Convert physical trace into a pixel trace
            msgs.info("Converting physical trace locations to nearest pixel")
            pixcen = artrace.phys_to_pix(0.5*(slf._lordloc[det-1]+slf._rordloc[det-1]), slf._pixlocn[det-1], 1)
            pixwid = (slf._rordloc[det-1]-slf._lordloc[det-1]).mean(0).astype(np.int)
            lordpix = artrace.phys_to_pix(slf._lordloc[det-1], slf._pixlocn[det-1], 1)
            rordpix = artrace.phys_to_pix(slf._rordloc[det-1], slf._pixlocn[det-1], 1)
            slf.SetFrame(slf._pixcen, pixcen, det)
            slf.SetFrame(slf._pixwid, pixwid, det)
            slf.SetFrame(slf._lordpix, lordpix, det)
            slf.SetFrame(slf._rordpix, rordpix, det)
            # Save QA for slit traces
            arqa.slit_trace_qa(slf, slf._mstrace[det-1], slf._lordpix[det-1], slf._rordpix[det-1], extord, desc="Trace of the slit edges")
            ###############
            # Prepare the pixel flat field frame
            update = slf.MasterFlatField(fitsdict, det)
            if update and reuseMaster: armbase.UpdateMasters(sciexp, sc, det, ftype="flat", chktype="pixflat")
            ###############
            # Derive the spectral tilt
            if slf._tilts[det-1] is None:
                # First time tilts are derived for this arc frame --> derive the order tilts
                tilts = None
                nitertilts = 2
                doQA = False
                for tt in range(nitertilts):
                    msgs.info("Iterating on spectral tilts -- Iteration {0:d}/{1:d}".format(tt+1, nitertilts))
                    if tt == nitertilts-1:
                        doQA = True
                    tilts, satmask = artrace.model_tilt(slf, det, slf._msarc[det-1], guesstilts=tilts, plotQA=doQA)
                slf.SetFrame(slf._tilts, tilts, det)
                slf.SetFrame(slf._satmask, tilts, det)

                # Setup arc parameters (e.g. linelist)
                arcparam = ararc.setup_param(slf, sc, det, fitsdict)
                slf.SetFrame(slf._arcparam, arcparam, det)
                ###############
                # Extract arc and identify lines
                wv_calib = ararc.simple_calib(slf, det)
                slf.SetFrame(slf._wvcalib, wv_calib, det)
                ###############
                # Generate a master wave frame
                update = slf.MasterWave(fitsdict, det)
                if update and reuseMaster:
                    armbase.UpdateMasters(sciexp, sc, det, ftype="arc", chktype="wave")


            ###############
            # Check if the user only wants to prepare the calibrations
            msgs.info("All calibration frames have been prepared")
            if slf._argflag['run']['preponly']:
                msgs.info("If you would like to continue with the reduction,"
                          +msgs.newline()+"disable the run+preponly command")
                continue
            ###############
            # Standard star (is this a calibration, e.g. goes above?)
            msgs.info("Processing standard star")
            msgs.warn("Assuming one star per detector mosaic")
            update = slf.MasterStandard(scidx, fitsdict)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="standard")

            ###############
            # Load the science frame and from this generate a Poisson error frame
            msgs.info("Loading science frame")
            sciframe = arload.load_frames(slf, fitsdict, [scidx], det,
                                          frametype='science',
                                          msbias=slf._msbias[det-1],
                                          transpose=slf._transpose)
            sciframe = sciframe[:, :, 0]
            # Extract
            msgs.info("Extracting science frame")
            reduce_frame(slf, sciframe, scidx, fitsdict, det)


            #continue
            #msgs.error("UP TO HERE")
            ###############
            # Perform a velocity correction
            if (slf._argflag['reduce']['heliocorr'] == True) & False:
                if slf._argflag['science']['load']['extracted'] == True:
                    msgs.warn("Heliocentric correction will not be applied if an extracted science frame exists, and is used")
                msgs.work("Perform a full barycentric correction")
                msgs.work("Include the facility to correct for gravitational redshifts and time delays (see Pulsar timing work)")
                msgs.info("Performing a heliocentric correction")
                # Load the header for the science frame
                slf._waveids = arvcorr.helio_corr(slf, scidx[0])
            else:
                msgs.info("A heliocentric correction will not be performed")

            ###############
            # Using model sky, calculate a flexure correction
            msgs.warn("Implement flexure correction!!")

            ###############
            # Flux
            msgs.work("Need to check for existing sensfunc")
            msgs.work("Consider using archived sensitivity if not found")
            msgs.info("Fluxing with {:s}".format(slf._sensfunc['std']['name']))
            arflux.apply_sensfunc(slf, scidx, fitsdict)

        # Write
        ario.write_1d_spectra(slf)
        # Close the QA for this object
        slf._qa.close()
        # Free up some memory by replacing the reduced ScienceExposure class
        sciexp[sc] = None
    return status


def reduce_frame(slf, sciframe, scidx, fitsdict, det, standard=False):
    """ Run standard extraction steps on a frame
    Parameters
    ----------
    sciframe : image
      Bias subtracted image (using arload.load_frame)
    scidx : int
      Index of the frame
    fitsdict : dict
      Contains relevant information from fits header files
    det : int
      Detector index
    standard : bool, optional
      Standard star frame?
    """
    # Convert ADUs to electrons
    sciframe *= slf._spect['det'][det-1]['gain']
    varframe = arproc.variance_frame(slf, det, sciframe, scidx, fitsdict)
    if not standard:
        slf._varframe[det-1] = varframe
    ###############
    # Subtract off the scattered light from the image
    msgs.work("Scattered light subtraction is not yet implemented...")
    ###############
    # Flat field the science frame
    if slf._argflag['reduce']['flatfield']:
        msgs.info("Flat fielding the science frame")
        sciframe = arproc.flatfield(slf, sciframe, slf._mspixflatnrm[det-1], det)
    else:
        msgs.info("Not performing a flat field calibration")
    if not standard:
        slf._sciframe[det-1] = sciframe
    ###############
    # Identify cosmic rays
    msgs.work("Include L.A.Cosmic arguments in the settings files")
    if True: crmask = arproc.lacosmic(slf, fitsdict, det, sciframe, scidx, grow=1.5)
    else: crmask = np.zeros(sciframe.shape)
    msgs.work("For now, perform extraction -- really should do this after the flexure+heliocentric correction")
    ###############
    # Estimate Sky Background
    if slf._argflag['reduce']['bgsubtraction']:
        # Perform an iterative background/science extraction
        msgs.info("Estimating the sky background")
        bgframe = arproc.bg_subtraction(slf, det, sciframe, varframe, crmask)
        varframe = arproc.variance_frame(slf, det, sciframe, scidx, fitsdict, skyframe=bgframe)
        if not standard: # Need to save
            slf._varframe[det-1] = varframe
            slf._bgframe[det-1] = bgframe
    ###############
    # Estimate trace of science objects
    scitrace = artrace.trace_object(slf, det, sciframe-bgframe, varframe, crmask)
    if scitrace is None:
        msgs.info("Not performing extraction for science frame"+msgs.newline()+slf._fitsdict['filename'][scidx[0]])
        pdb.set_trace()
        #continue
    ###############
    # Finalize the Sky Background image
    if slf._argflag['reduce']['bgsubtraction']:
        # Perform an iterative background/science extraction
        msgs.info("Finalizing the sky background image")
        trcmask = scitrace['object'].sum(axis=2)
        trcmask[np.where(trcmask>0.0)] = 1.0
        bgframe = arproc.bg_subtraction(slf, det, sciframe, varframe, crmask, tracemask=trcmask)
        # Redetermine the variance frame based on the new sky model
        varframe = arproc.variance_frame(slf, det, sciframe, scidx, fitsdict, skyframe=bgframe)
        # Save
        if not standard:
            slf._varframe[det-1] = varframe
            slf._bgframe[det-1] = bgframe
    ###############
    # Determine the final trace of the science objects
    msgs.info("Final trace")
    scitrace = artrace.trace_object(slf, det, sciframe-bgframe, varframe, crmask)
    if standard:
        slf._msstd[det-1]['trace'] = scitrace
        specobjs = arspecobj.init_exp(slf, scidx, det, fitsdict,
                                                         trc_img=scitrace,
                                                         objtype='standard')
        slf._msstd[det-1]['spobjs'] = specobjs
    else:
        slf._scitrace[det-1] = scitrace
        # Generate SpecObjExp list
        specobjs = arspecobj.init_exp(slf, scidx, det, fitsdict,
                                      trc_img=scitrace, objtype='science')
        slf._specobjs[det-1] = specobjs
    ###############
    # Extract
    if scitrace is None:
        msgs.info("Not performing extraction for science frame"+msgs.newline()+slf._fitsdict['filename'][scidx[0]])
        pdb.set_trace()
        #continue
    # Boxcar
    arextract.boxcar(slf, det, specobjs, sciframe-bgframe, varframe, bgframe, crmask, scitrace)
    # Return
    return True
