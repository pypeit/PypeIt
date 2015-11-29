import pdb
import numpy as np
import armasters
import artrace
import arutils
import armbase
import arproc
import ararc
import arqa


def ARMLSD(argflag, spect, fitsdict, msgs, reuseMaster=False):
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
    msgs : class
      Messages class used to log data reduction process
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
    sciexp = armbase.SetupScience(argflag, spect, fitsdict, msgs)
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
            fitsdict = arproc.get_ampsec_trimmed(slf, fitsdict, det, scidx, msgs)
            ###############
            # Generate master bias frame
            update = slf.MasterBias(fitsdict, det, msgs)
            if update and reuseMaster: armbase.UpdateMasters(sciexp, sc, det, msgs, ftype="bias")
            ###############
            # Generate a bad pixel mask (should not repeat)
            update = slf.BadPixelMask(det, msgs)
            if update and reuseMaster: armbase.UpdateMasters(sciexp, sc, det, msgs, ftype="arc")
            ###############
            # Estimate gain and readout noise for the amplifiers
            msgs.work("Estimate Gain and Readout noise from the raw frames...")
            ###############
            # Generate a master arc frame
            update = slf.MasterArc(fitsdict, det, msgs)
            if update and reuseMaster: armbase.UpdateMasters(sciexp, sc, det, msgs, ftype="arc")
            ###############
            # Determine the dispersion direction (and transpose if necessary)
            slf.GetDispersionDirection(fitsdict, det, msgs)
            if slf._bpix[det-1] is None:
                slf.SetFrame(slf._bpix, np.zeros((slf._nspec[det-1], slf._nspat[det-1])), det)
            ###############
            # Generate a master trace frame
            update = slf.MasterTrace(fitsdict, det, msgs)
            if update and reuseMaster: armbase.UpdateMasters(sciexp, sc, det, msgs, ftype="flat", chktype="trace")
            ###############
            # Generate an array that provides the physical pixel locations on the detector
            slf.GetPixelLocations(det, msgs)
            ###############
            # Determine the edges of the spectrum (spatial)
            lordloc, rordloc, extord = artrace.trace_orders(slf, slf._mstrace[det-1], det, msgs, singleSlit=True,
                                                            pcadesc="PCA trace of the slit edges")
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
            arqa.slit_trace_qa(slf, slf._mstrace[det-1], slf._lordpix[det-1], slf._rordpix[det-1],
                               extord, desc="Trace of the slit edges")
            ###############
            # Prepare the pixel flat field frame
            update = slf.MasterFlatField(fitsdict, det, msgs)
            if update and reuseMaster: armbase.UpdateMasters(sciexp, sc, det, msgs, ftype="flat", chktype="pixflat")
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
                    tilts, satmask = artrace.model_tilt(slf, det, slf._msarc[det-1], msgs,
                                                        guesstilts=tilts, plotQA=doQA)
                slf.SetFrame(slf._tilts, tilts, det)
                slf.SetFrame(slf._satmask, tilts, det)


                ################################################################
                # Temporary break until core structure is fixed
                slf._qa.close()
                msgs.error("UP TO HERE -- DELETE THE QA PLOT CLOSE HERE!")


                # Setup arc parameters (e.g. linelist)
                slf._arcparam = ararc.setup_param(slf, sc)
                ###############
                # Extract arc and identify lines
                slf.wv_calib = ararc.simple_calib(slf, det)
                # Generate Wavelength Image
                slf._mswvimg = arutils.func_val(slf.wv_calib['fitc'], slf._tilts, slf.wv_calib['function'], minv=slf.wv_calib['fmin'], maxv=slf.wv_calib['fmax'])
                ind = slf._spect['arc']['index'][sc]
                slf._mswvimg_name = "{0:s}/{1:s}/mswvimg{2:s}_{3:03d}.fits".format(os.getcwd(),slf._argflag['run']['masterdir'],slf._spect["det"][det-1]["suffix"],len(slf._done_arcs)-1)
                ind = slf._spect['arc']['index'][sc]

            ###############
            # Check if the user only wants to prepare the calibrations
            msgs.info("All calibration frames have been prepared")
            if slf._argflag['run']['preponly']:
                msgs.info("If you would like to continue with the reduction,"+msgs.newline()+"disable the run+preponly command")
                continue
            ###############
            # Load the science frame and from this generate a Poisson error frame
            sciframe = arload.load_frames(slf, scidx, det, frametype='science', msbias=slf._msbias, transpose=slf._transpose)
            sciframe = sciframe[:,:,0]
            # Convert ADUs to electrons
            sciframe *= slf._spect['det'][det-1]['gain']
            varframe = arproc.variance_frame(slf, det, sciframe, scidx[0])
            # Write
            msvar_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), slf._argflag['run']['masterdir'], slf.target, 0, "var")
            arsave.save_master(slf, varframe, filename=msvar_name, frametype='variance')
            ###############
            # Subtract off the scattered light from the image
            msgs.work("Scattered light subtraction is not yet implemented...")
            ###############
            # Flat field the science frame
            if slf._argflag['reduce']['flatfield']:
                msgs.info("Flat fielding the science frame")
                sciframe = arproc.flatfield(slf, sciframe, mspixflatnrm)
            else:
                msgs.info("Not performing a flat field calibration")
            ###############
            # Identify cosmic rays
            msgs.work("Include L.A.Cosmic arguments in the settings files")
            #if slf._argflag['reduce']['crmask']:
            if True: crmask = arproc.lacosmic(slf, det, sciframe, grow=1.5)
            else: crmask = np.zeros(sciframe.shape)
            #arutils.ds9plot(crmask)
            #arutils.ds9plot(sciframe*(1.0-crmask))
            msgs.work("For now, perform extraction -- really should do this after the flexure+heliocentric correction")
            ###############
            # Estimate Sky Background
            if slf._argflag['reduce']['bgsubtraction']:
                # Perform an iterative background/science extraction
                msgs.info("Estimating the sky background")
                bgframe = arproc.bg_subtraction(slf, det, sciframe, varframe, crmask)
                #scibgsub, bgframe = arproc.background_subtraction(slf, sciframe, varframe)
                # Derive a suitable name for the master sky background frame
                msgs.work("Include an index suffix for each object frame")# e.g. if you have 3 frames of the same object, include a common integer suffix on the filenames
                msbg_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), slf._argflag['run']['masterdir'], slf.target, 0, "sky")
                # Send the data away to be saved
                arsave.save_master(slf, bgframe, filename=msbg_name, frametype='sky background')
                # Derive a suitable name for the sky-subtracted science frame
                msscibg_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), slf._argflag['run']['masterdir'], slf.target, 0, "skysub")
                # Send the data away to be saved
                arsave.save_master(slf, sciframe-bgframe, filename=msscibg_name, frametype='sky subtracted science')
                # Redetermine the variance frame based on the new sky model
                varframe = arproc.variance_frame(slf, det, sciframe, scidx[0], skyframe=bgframe)
            ###############
            # Estimate trace of science objects
            scitrace = artrace.trace_object(slf, sciframe-bgframe, varframe, crmask)
            if scitrace is None:
                msgs.info("Not performing extraction for science frame"+msgs.newline()+slf._fitsdict['filename'][scidx[0]])
                continue
            ###############
            # Finalize the Sky Background image
            if slf._argflag['reduce']['bgsubtraction']:
                # Perform an iterative background/science extraction
                msgs.info("Finalizing the sky background image")
                trcmask = scitrace['object'].sum(axis=2)
                trcmask[np.where(trcmask>0.0)] = 1.0
                bgframe = arproc.bg_subtraction(slf, det, sciframe, varframe, crmask, tracemask=trcmask)
                # Derive a suitable name for the master sky background frame
                msgs.work("Include an index suffix for each object frame")# e.g. if you have 3 frames of the same object, include a common integer suffix on the filenames
                msbg_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), slf._argflag['run']['masterdir'], slf.target, 0, "sky")
                # Send the data away to be saved
                arsave.save_master(slf, bgframe, filename=msbg_name, frametype='sky background')
                # Derive a suitable name for the sky-subtracted science frame
                msscibg_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), slf._argflag['run']['masterdir'], slf.target, 0, "skysub")
                # Send the data away to be saved
                arsave.save_master(slf, sciframe-bgframe, filename=msscibg_name, frametype='sky subtracted science')
                # Redetermine the variance frame based on the new sky model
                varframe = arproc.variance_frame(slf, det, sciframe, scidx[0], skyframe=bgframe)
            ###############
            # Determine the final trace of the science objects
            scitrace = artrace.trace_object(slf, sciframe-bgframe, varframe, crmask)
            # Write
            mstrc_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), slf._argflag['run']['masterdir'], slf.target, 0, "objtrc")
            hdutrc = fits.PrimaryHDU(scitrace['traces'])
            hduobj = fits.ImageHDU(scitrace['object'])
            hdulist = fits.HDUList([hdutrc, hduobj])
            hdulist.writeto(mstrc_name,clobber=True)
            msgs.info("Wrote object trace file: {:s}".format(mstrc_name))
            # Generate SpecObjExp list
            slf._specobjs += arspecobj.init_exp(slf,sc,det,trc_img=scitrace, objtype=sctype)
            ###############
            # Extract
            if scitrace is None:
                msgs.info("Not performing extraction for science frame"+msgs.newline()+slf._fitsdict['filename'][scidx[0]])
                continue
            # Boxcar Extraction
            arextract.boxcar(slf, sciframe-bgframe, varframe, bgframe, crmask, scitrace)
            #Generate and Write spectra
            if False:
                sig = np.sqrt(slf._specobjs[0].boxcar['var'])
                wave = slf._specobjs[0].boxcar['wave']
                flux = slf._specobjs[0].boxcar['counts']
                sky = slf._specobjs[0].boxcar['sky']
                #
                xspec = XSpectrum1D.from_tuple( (wave,flux,sig) )
                spec_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), slf._argflag['run']['masterdir'], slf.target, 0, "boxcar")
                msgs.info("Writing boxcar spectrum: {:s}".format(spec_name))
                xspec.write_to_fits(spec_name, clobber=True)

                skyspec = XSpectrum1D.from_tuple( (wave,sky) )
                skyspec_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), slf._argflag['run']['masterdir'], slf.target, 0, "skybox")
                msgs.info("Writing sky spectrum: {:s}".format(skyspec_name))
                skyspec.write_to_fits(skyspec_name, clobber=True)

            ###############
            # If standard, generate a sensitivity function
            if (sctype == 'standard') & (det == slf._spect['mosaic']['ndet']):
                if sc > 0:
                    msgs.error("What to do with multiple standard exposures??")
                else:
                    msgs.warn("Need to check for existing sensfunc as with Arc, Trace")
                    slf._sensfunc = arflux.generate_sensfunc(slf,sc)
                    # Write
                    msgs.warn("Need to write sensfunc to hard drive")
                    #sensfunc_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.yaml".format(os.getcwd(), slf._argflag['run']['masterdir'], slf._fitsdict['target'][scidx[0]], 0, "sensfunc")
                    #msgs.info("Writing sensfunc: {:s}".format(sensfunc_name))
                    #with open(sensfunc_name, 'w') as yamlf:
                    #    yamlf.write( yaml.dump(slf._sensfunc))
                    #with io.open(sensfunc_name, 'w', encoding='utf-8') as f:
                    #    f.write(unicode(json.dumps(slf._sensfunc, sort_keys=True, indent=4, separators=(',', ': '))))

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
            if sctype == 'science':
                msgs.warn("Implement flexure correction!!")

            ###############
            # Determine the wavelength scale (i.e. the wavelength of each pixel) to be used during the extraction
            '''
            if sctype == 'science':
                msgs.info("Generating the array of extraction wavelengths")
                slf._wavelength = arproc.get_wscale(slf)
            '''

            ###############
            # Flux
            if sctype == 'science':
                msgs.work("Need to check for existing sensfunc")
                msgs.work("Consider using archived sensitivity if not found")
                msgs.info("Fluxing with {:s}".format(slf._sensfunc['std']['name']))
                arflux.apply_sensfunc(slf,sc)

            ###############
            # Append for later stages (e.g. coadding)
            slf._allspecobjs += slf._specobjs
            if sctype == 'science':
                msgs.error("STOP FOR NOW")



        # Free up some memory by replacing the reduced ScienceExposure class
        # Close the QA for this object
        slf._qa.close()
        sciexp[sc] = None
    return status
