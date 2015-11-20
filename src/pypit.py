#import matplotlib
#matplotlib.use('Agg')  # For Travis

import os
import sys
import getopt
#import json, io, yaml
from signal import SIGINT, signal as sigsignal
from warnings import resetwarnings, simplefilter
from time import time
import traceback
import numpy as np
from astropy.io import fits
from astropy import units as u
# Import PYPIT routines
import armsgs as msgs
import ararc
import arextract
import arflux
import arload
import arsave
import arcomb
import arproc
import arqa
import arspecobj
import artrace
import arutils
import arvcorr

try:
    from linetools.spectra.xspectrum1d import XSpectrum1D
except:
    pass
#else:
#    from astropy import units as u

import pdb

try:
    from xastropy.xutils import xdebug as xdb
except:
    pass

last_updated = "Last updated 09 November 2015"
version = 'v0.2'

def usage(prognm):
    print "\n#####################################################################"
    print msgs.pypitheader()
    print "##  -----------------------------------------------------------------"
    print "##  Options: (default values in brackets)"
    print "##   -c or --cpus      : (all) Number of cpu cores to use"
    print "##   -h or --help      : Print this message"
    print "##   -v or --verbose   : (2) Level of verbosity (0-2)"
    print "##  -----------------------------------------------------------------"
    print "##  %s" % last_updated
    print "#####################################################################\n"
    sys.exit()

def PYPIT(argflag, quick=False):
    """
    Main driver of the PYPIT code. Default settings and
    user-specified changes are made, and passed to the
    appropriate code for data reduction.

    Parameters
    ----------
    argflag : dict
      Arguments and flags used for reduction
    quick : bool
      If True, a quick reduction (but possibly less
      accurate) will be performed. This flag is most
      useful for observing at a telescope, but not
      for publication quality results.
    ---------------------------------------------------
    """

    # First send all signals to messages to be dealt with (i.e. someone hits ctrl+c)
    sigsignal(SIGINT, msgs.signal_handler)

    # Ignore all warnings given by python
    resetwarnings()
    simplefilter("ignore")

    # Record the starting time
    tstart=time()

    # Load the Input file
    argflag, parlines, datlines, spclines = arload.load_input(argflag)

    # If a quick reduction has been requested, make sure the requested pipeline
    # is the quick implementation (if it exists), otherwise run the standard pipeline.
    if quick:
        # Change to a "quick" settings file
        msgs.work("TO BE DONE")

    # Load the Spectrograph settings
    spect = arload.load_spect(argflag)

    # Load any changes to the spectrograph settings
    spect = arload.load_spect(argflag, spect=spect, lines=spclines)

    # Load the important information from the fits headers
    fitsdict = arload.load_headers(argflag, spect, datlines)

    # Reduce the data!
    status = 0
    msgs.work("Make appropriate changes to quick reduction")
    if quick:
        msgs.work("define what is needed here")
    # Send the data away to be reduced
    if spect['mosaic']['reduction']=='ARMLSD':
        msgs.info("Data reduction will be performed using PYPIT-ARMLSD")
        import pypArmlsd
        status = pypArmlsd.ARMLSD(argflag, spect, fitsdict)
    elif spect['mosaic']['reduction']=='ARMED':
        msgs.info("Data reduction will be performed using PYPIT-ARMED")
        import pypArmed
        status = pypArmed.ARMED(argflag, spect, fitsdict)
    # Check for successful reduction
    if status == 0:
        msgs.info("Data reduction complete")
    else:
        msgs.error("Data reduction failed with status ID {0:d}".format(status))
    # Capture the end time and print it to user
    tend = time()
    codetime = tend-tstart
    if codetime < 60.0:
        msgs.info("Data reduction execution time: {0:.2f}s".format(codetime))
    elif codetime/60.0 < 60.0:
        mns = int(codetime/60.0)
        scs = codetime - 60.0*mns
        msgs.info("Data reduction execution time: {0:d}m {1:.2f}s".format(mns, scs))
    else:
        hrs = int(codetime/3600.0)
        mns = int(60.0*(codetime/3600.0 - hrs))
        scs = codetime - 60.0*mns - 3600.0*hrs
        msgs.info("Data reduction execution time: {0:d}h {1:d}m {2:.2f}s".format(hrs, mns, scs))
    return


###################################
# Reduction pipelines
###################################

def ARMLSD(self):
    """
    Automatic Reduction & Modelling of Long Slit Data
    """
    success = False
    # INIT
    self._allspecobjs = []
    # Sort data and match calibrations to science+standard frames
    self.Setup()
    sci = self._filesort['science']
    std = self._filesort['standard']
    numsci = np.size(sci)
    numstd = np.size(std)
    msgs.bug("I suspect the above doesn't organize the standard correctly")
    msgs.bug("We want the ones matched to the science frames.")
    if numstd == 0:
        msgs.bug("We will proceed without a standard star, but don't ask to flux calibrate!")
    if numsci == 0:
        msgs.bug("What to do if no science frames are input? The calibrations should still be processed.")
        msgs.work("Maybe assume that each non-identical arc is a science frame, but force no science extraction")
        msgs.error("No science frames are input!!")
        numsci=1
        self._argflag['run']['preponly'] = True # Prepare the calibration files, but don't do the science frame reduction
    # Reduce each standard+science frame entirely before moving to the next one
    stdsci_idx = range(numstd) + range(numsci)
    stdsci_types = ['standard']*numstd + ['science']*numsci
    if stdsci_idx[0] != 0:
        msgs.error("First index should be 0")
    for qq,sc in enumerate(stdsci_idx):
        # Some book-keeping
        sctype = stdsci_types[qq]
        scidx = self._spect[sctype]['index'][sc]
        self._scidx = scidx[0]
        sciext_name_p, sciext_name_e = os.path.splitext(self._fitsdict['filename'][scidx[0]])
        self._specobjs = []
        self.target = ''.join(self._fitsdict['target'][scidx[0]].split())
        msgs.info("Working on file {:s}, target {:s}".format(self._fitsdict['filename'][scidx[0]],self.target))
        ###############
        # First set the index for the science frame
        # Now loop on Detectors
        for kk in xrange(self._spect['mosaic']['ndet']):
            det = kk + 1 # Detectors indexed from 1
            prefix = "{0:s}{1:s}".format(sciext_name_p, self._spect["det"][kk]["suffix"])
            ###############
            # Get amplifier sections
            self._ampsec = arproc.get_ampsec_trimmed(self, det, scidx[0])
            ###############
            # Generate master bias frame
            self._msbias, self._msbias_name = self.MasterBias(sc, det)
            ###############
            # Generate a bad pixel mask (should not repeat)
            if not hasattr(self,'_bpix'):
                self.BadPixelMask(sc, det)
            else:
                msgs.info("Using previously generated BadPixelMask")
                msgs.work("Make sure a new bad pixel mask is created for a new setup (e.g. different binning)")
            ###############
            # Estimate gain and readout noise for the amplifiers
            msgs.work("Estimate Gain and Readout noise from the raw frames...")
            ###############
            # Generate a master arc frame
            self._msarc, self._msarc_name = self.MasterArc(sc, det)
            if (not hasattr(self,'_bpix')) or (self._bpix is None):
                self._bpix = np.zeros_like(self._msarc)
            ###############
            # Determine the dispersion direction (and transpose if necessary) only on the first pass through
            self.GetDispersionDirection(self._spect['arc']['index'][sc], det)
            ###############
            # Generate a master trace frame
            self._mstrace, self._mstrace_name = self.MasterTrace(sc, det)
            ###############
            # Generate an array that provides the physical pixel locations on the detector
            self.GetPixelLocations(det)
            ###############
            # Determine the edges of the spectrum (spatial)
            self._lordloc, self._rordloc = artrace.trace_orders(self, self._mstrace, ARMLSD=True, prefix=prefix, trcprefix=self._trcprefix)
            arsave.save_ordloc(self, self._mstrace_name)
            arqa.slit_trace_qa(self._mstrace, self._lordloc, self._rordloc, self._mstrace_name)
            # Convert physical trace into a pixel trace
            msgs.info("Converting physical trace locations to nearest pixel")
            self._pixcen  = artrace.phys_to_pix(0.5*(self._lordloc+self._rordloc), self._pixlocn, self._dispaxis, 1-self._dispaxis)
            self._pixwid  = (self._rordloc-self._lordloc).mean(0).astype(np.int)
            self._lordpix = artrace.phys_to_pix(self._lordloc, self._pixlocn, self._dispaxis, 1-self._dispaxis)
            self._rordpix = artrace.phys_to_pix(self._rordloc, self._pixlocn, self._dispaxis, 1-self._dispaxis)
            ###############
            # Prepare the pixel flat field frame
            mspixflatnrm, msblaze = arproc.flatnorm(self, det, self._mstrace.copy(), overpix=0, fname=os.path.basename(os.path.splitext(self._mstrace_name)[0]))
            self._mspixflat, self._mspixflatnrm, self._msblaze, self._mspixflat_name = self.MasterFlatField(sc,det)
            ###############
            # Derive the spectral tilt
            if self._foundarc:
                try:
                    # Load the order locations from file
                    self._tilts, self._satmask = arload.load_tilts(self._msarc_name)
                except:
                    # If this fails, rederive the order tilts
                    self._tilts, self._satmask = artrace.model_tilt(self, self._msarc, prefix=prefix, trcprefix=self._trcprefix, tltprefix=self._tltprefix)
                    # Save the traces
                    arsave.save_tilts(self, self._msarc_name)
            else:
                # First time encountered this set of arc frames --> derive the order tilts
                tilts = None
                nitertilts=1
                doQA = False
                doRefine = False
                for tt in range(nitertilts):
                    msgs.info("Iterating on spectral tilts -- Iteration {0:d}/{1:d}".format(tt+1,nitertilts))
                    if tt==nitertilts-1:
                        doQA=True
                        doRefine=True
                    tilts, satmask = artrace.model_tilt(self, det, self._msarc, prefix=prefix, trcprefix=self._trcprefix, tltprefix=self._tltprefix, guesstilts=tilts, plotQA=doQA, refine_tilts=doRefine)
                self._tilts, self._satmask = tilts, satmask
                # Save the tilts
                arsave.save_tilts(self, self._msarc_name)
                msgs.bug("Need to include the definitions below in the above exception as a failsafe")
                # Setup arc parameters (e.g. linelist)
                self._arcparam = ararc.setup_param(self, sc)
                ###############
                # Extract arc and identify lines
                self.wv_calib = ararc.simple_calib(self, det)
                # Generate Wavelength Image
                self._mswvimg = arutils.func_val(self.wv_calib['fitc'], self._tilts, self.wv_calib['function'], min=self.wv_calib['fmin'], max=self.wv_calib['fmax'])
                ind = self._spect['arc']['index'][sc]
                self._mswvimg_name = "{0:s}/{1:s}/mswvimg{2:s}_{3:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],self._spect["det"][det-1]["suffix"],len(self._done_arcs)-1)
                ind = self._spect['arc']['index'][sc]
                arsave.save_master(self, self._mswvimg, filename=self._mswvimg_name, frametype='wvimg', ind=ind)

            ###############
            # Check if the user only wants to prepare the calibrations
            msgs.info("All calibration frames have been prepared")
            if self._argflag['run']['preponly']:
                msgs.info("If you would like to continue with the reduction,"+msgs.newline()+"disable the run+preponly command")
                continue
            ###############
            # Load the science frame and from this generate a Poisson error frame
            sciframe = arload.load_frames(self, scidx, det, frametype='science', msbias=self._msbias, transpose=self._transpose)
            sciframe = sciframe[:,:,0]
            # Convert ADUs to electrons
            sciframe *= self._spect['det'][det-1]['gain']
            varframe = arproc.variance_frame(self, det, sciframe, scidx[0])
            # Write
            msvar_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), self._argflag['run']['masterdir'], self.target, 0, "var")
            arsave.save_master(self, varframe, filename=msvar_name, frametype='variance')
            ###############
            # Subtract off the scattered light from the image
            msgs.work("Scattered light subtraction is not yet implemented...")
            ###############
            # Flat field the science frame
            if self._argflag['reduce']['flatfield']:
                msgs.info("Flat fielding the science frame")
                sciframe = arproc.flatfield(self, sciframe, mspixflatnrm)
            else:
                msgs.info("Not performing a flat field calibration")
            ###############
            # Identify cosmic rays
            msgs.work("Include L.A.Cosmic arguments in the settings files")
            #if self._argflag['reduce']['crmask']:
            if True: crmask = arproc.lacosmic(self, det, sciframe, grow=1.5)
            else: crmask = np.zeros(sciframe.shape)
            #arutils.ds9plot(crmask)
            #arutils.ds9plot(sciframe*(1.0-crmask))
            msgs.work("For now, perform extraction -- really should do this after the flexure+heliocentric correction")
            ###############
            # Estimate Sky Background
            if self._argflag['reduce']['bgsubtraction']:
                # Perform an iterative background/science extraction
                msgs.info("Estimating the sky background")
                bgframe = arproc.bg_subtraction(self, det, sciframe, varframe, crmask)
                #scibgsub, bgframe = arproc.background_subtraction(self, sciframe, varframe)
                # Derive a suitable name for the master sky background frame
                msgs.work("Include an index suffix for each object frame")# e.g. if you have 3 frames of the same object, include a common integer suffix on the filenames
                msbg_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), self._argflag['run']['masterdir'], self.target, 0, "sky")
                # Send the data away to be saved
                arsave.save_master(self, bgframe, filename=msbg_name, frametype='sky background')
                # Derive a suitable name for the sky-subtracted science frame
                msscibg_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), self._argflag['run']['masterdir'], self.target, 0, "skysub")
                # Send the data away to be saved
                arsave.save_master(self, sciframe-bgframe, filename=msscibg_name, frametype='sky subtracted science')
                # Redetermine the variance frame based on the new sky model
                varframe = arproc.variance_frame(self, det, sciframe, scidx[0], skyframe=bgframe)
            ###############
            # Estimate trace of science objects
            scitrace = artrace.trace_object(self, sciframe-bgframe, varframe, crmask)
            if scitrace is None:
                msgs.info("Not performing extraction for science frame"+msgs.newline()+self._fitsdict['filename'][scidx[0]])
                continue
            ###############
            # Finalize the Sky Background image
            if self._argflag['reduce']['bgsubtraction']:
                # Perform an iterative background/science extraction
                msgs.info("Finalizing the sky background image")
                trcmask = scitrace['object'].sum(axis=2)
                trcmask[np.where(trcmask>0.0)] = 1.0
                bgframe = arproc.bg_subtraction(self, det, sciframe, varframe, crmask, tracemask=trcmask)
                # Derive a suitable name for the master sky background frame
                msgs.work("Include an index suffix for each object frame")# e.g. if you have 3 frames of the same object, include a common integer suffix on the filenames
                msbg_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), self._argflag['run']['masterdir'], self.target, 0, "sky")
                # Send the data away to be saved
                arsave.save_master(self, bgframe, filename=msbg_name, frametype='sky background')
                # Derive a suitable name for the sky-subtracted science frame
                msscibg_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), self._argflag['run']['masterdir'], self.target, 0, "skysub")
                # Send the data away to be saved
                arsave.save_master(self, sciframe-bgframe, filename=msscibg_name, frametype='sky subtracted science')
                # Redetermine the variance frame based on the new sky model
                varframe = arproc.variance_frame(self, det, sciframe, scidx[0], skyframe=bgframe)
            ###############
            # Determine the final trace of the science objects
            scitrace = artrace.trace_object(self, sciframe-bgframe, varframe, crmask)
            # Write
            mstrc_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), self._argflag['run']['masterdir'], self.target, 0, "objtrc")
            hdutrc = fits.PrimaryHDU(scitrace['traces'])
            hduobj = fits.ImageHDU(scitrace['object'])
            hdulist = fits.HDUList([hdutrc, hduobj])
            hdulist.writeto(mstrc_name,clobber=True)
            msgs.info("Wrote object trace file: {:s}".format(mstrc_name))
            # Generate SpecObjExp list
            self._specobjs += arspecobj.init_exp(self,sc,det,trc_img=scitrace, objtype=sctype)
            ###############
            # Extract
            if scitrace is None:
                msgs.info("Not performing extraction for science frame"+msgs.newline()+self._fitsdict['filename'][scidx[0]])
                continue
            # Boxcar Extraction
            arextract.boxcar(self, sciframe-bgframe, varframe, bgframe, crmask, scitrace)
            #Generate and Write spectra
            if False:
                sig = np.sqrt(self._specobjs[0].boxcar['var'])
                wave = self._specobjs[0].boxcar['wave']
                flux = self._specobjs[0].boxcar['counts']
                sky = self._specobjs[0].boxcar['sky']
                #
                xspec = XSpectrum1D.from_tuple( (wave,flux,sig) )
                spec_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), self._argflag['run']['masterdir'], self.target, 0, "boxcar")
                msgs.info("Writing boxcar spectrum: {:s}".format(spec_name))
                xspec.write_to_fits(spec_name, clobber=True)

                skyspec = XSpectrum1D.from_tuple( (wave,sky) )
                skyspec_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), self._argflag['run']['masterdir'], self.target, 0, "skybox")
                msgs.info("Writing sky spectrum: {:s}".format(skyspec_name))
                skyspec.write_to_fits(skyspec_name, clobber=True)

            ###############
            # If standard, generate a sensitivity function
            if (sctype == 'standard') & (det == self._spect['mosaic']['ndet']):
                if sc > 0:
                    msgs.error("What to do with multiple standard exposures??")
                else:
                    msgs.warn("Need to check for existing sensfunc as with Arc, Trace")
                    self._sensfunc = arflux.generate_sensfunc(self,sc)
                    # Write
                    msgs.warn("Need to write sensfunc to hard drive")
                    #sensfunc_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.yaml".format(os.getcwd(), self._argflag['run']['masterdir'], self._fitsdict['target'][scidx[0]], 0, "sensfunc")
                    #msgs.info("Writing sensfunc: {:s}".format(sensfunc_name))
                    #with open(sensfunc_name, 'w') as yamlf:
                    #    yamlf.write( yaml.dump(self._sensfunc))
                    #with io.open(sensfunc_name, 'w', encoding='utf-8') as f:
                    #    f.write(unicode(json.dumps(self._sensfunc, sort_keys=True, indent=4, separators=(',', ': '))))

            #continue
            #msgs.error("UP TO HERE")
            ###############
            # Perform a velocity correction
            if (self._argflag['reduce']['heliocorr'] == True) & False:
                if self._argflag['science']['load']['extracted'] == True:
                    msgs.warn("Heliocentric correction will not be applied if an extracted science frame exists, and is used")
                msgs.work("Perform a full barycentric correction")
                msgs.work("Include the facility to correct for gravitational redshifts and time delays (see Pulsar timing work)")
                msgs.info("Performing a heliocentric correction")
                # Load the header for the science frame
                self._waveids = arvcorr.helio_corr(self, scidx[0])
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
                self._wavelength = arproc.get_wscale(self)
            '''

            ###############
            # Flux
            if sctype == 'science':
                msgs.work("Need to check for existing sensfunc")
                msgs.work("Consider using archived sensitivity if not found")
                msgs.info("Fluxing with {:s}".format(self._sensfunc['std']['name']))
                arflux.apply_sensfunc(self,sc)

            ###############
            # Append for later stages (e.g. coadding)
            self._allspecobjs += self._specobjs
            if sctype == 'science':
                msgs.error("STOP FOR NOW")


    # Insert remaining reduction steps here
    return success

if __name__ == "__main__":
    prognm = sys.argv[0]
    debug = True
    quick = False

    # Load options from command line
    try:
        opt,arg=getopt.getopt(sys.argv[1:],'hqc:v:', ['help',
                                                      'quick'])
        argflag = arload.optarg(sys.argv, last_updated)
    except getopt.GetoptError, err:
        msgs.error(err.msg)
        usage(prognm)
    for o,a in opt:
        if   o in ('-h', '--help')      : usage(argflag)
        elif o in ('-q', '--quick')     : quick = True
#		elif o in ('-c', '--cpus')      : argflag['run']['ncpus']     = a
#		elif o in ('-v', '--verbose')   : argflag['out']['verbose']   = int(a)

    if debug:
        PYPIT(argflag, quick=quick)
    else:
        try:
            PYPIT(argflag, quick=quick)
        except Exception:
            # There is a bug in the code, print the file and line number of the error.
            et, ev, tb = sys.exc_info()
            while tb:
                co = tb.tb_frame.f_code
                filename = str(co.co_filename)
                line_no =  str(traceback.tb_lineno(tb))
                tb = tb.tb_next
            filename=filename.split('/')[-1]
            msgs.bug("There appears to be a bug on Line "+line_no+" of "+filename+" with error:"+msgs.newline()+str(ev)+msgs.newline()+"---> please contact the author")
