#import matplotlib
#matplotlib.use('Agg')  # For Travis

import os
import sys
import getopt
from signal import SIGINT, signal as sigsignal
from warnings import resetwarnings, simplefilter
from time import time
import traceback
import numpy as np
# Import PYPIT routines
import armsgs as msgs
import ararc
import arload
import arsave
import arcomb
import arproc
import arsort
import artrace
import arutils
import arvcorr
import arextract
import arqa

import pdb

try:
    from xastropy.xutils import xdebug as xdb
except:
    pass

last_updated = "Last updated 20th July 2015"

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

class ClassMain:

    def __init__(self, argflag, quick=False):
        """
        argflag :: A list of arguments and flags for the reduction
        quick   :: Results in a quick reduction, if a quicker (but less accurate) reduction exists
        ---------------------------------------------------

        """
        #pdb.set_trace()
        #############################
        # Set some universal parameters
        self._argflag = argflag   # Arguments and Flags
        self._transpose = False   # Determine if the frames need to be transposed
        # Arrays to store the name for the frames that have already been combined
        self._done_bias, self._name_bias = [], []
        self._done_flat, self._name_flat = [], []
        self._done_arcs, self._name_arcs = [], []
        #############################

        # First send all signals to messages to be dealt
        # with (i.e. someone hits ctrl+c)
        sigsignal(SIGINT, msgs.signal_handler)

        # Ignore all warnings given by python
        resetwarnings()
        simplefilter("ignore")

        # Record the starting time
        self._tstart=time()

        # Load the Input file
        self._parlines, self._datlines, self._spclines = arload.load_input(self)

        # Determine the type of data that is being reduced
        msgs.work("TO BE DONE")

        # If a quick reduction has been requested, make sure the requested pipeline
        # is the quick implementation (if it exists), otherwise run the standard pipeline.
        if quick:
            # Change to a "quick" settings file
            msgs.work("TO BE DONE")

        # Load the Spectrograph settings
        self._spect = arload.load_spect(self)

        # Load any changes to the spectrograph settings
        self._spect = arload.load_spect(self, lines=self._spclines)

        # Load the important information from the fits headers
        self._fitsdict = arload.load_headers(self)

        # Load the list of standard stars
        self._standardStars = arload.load_standards(self)

        # Reduce the data!
        msgs.work("Send the data away to a definition of the type of reduction needed")
        status = 0
        if quick:
            msgs.work("define what is needed here")
        else:
            success = self.ARMLSD()
        if status==0:
            msgs.info("Reduction complete")

    ###################################
    # Reduction procedures
    ###################################

    def BadPixelMask(self, sc, det):
        if self._argflag['reduce']['badpix']:
            msgs.info("Preparing a bad pixel mask")
            # Get all of the bias frames for this science frame
            if len(self._spect['bias']['index']) == 0:
                msgs.warn("No bias frames available to determine bad pixel mask")
                msgs.info("Not preparing a bad pixel mask")
                self._bpix = None
                return
            ind = self._spect['bias']['index'][sc]
            # Load the Bias frames
            frames = arload.load_frames(self, ind, det, frametype='bias', trim=False)
            tbpix = arcomb.comb_frames(frames, det, spect=self._spect, frametype='bias', **self._argflag['bias']['comb'])
            self._bpix = arproc.badpix(self, det, tbpix)
            del tbpix
        else:
            self._bpix = None
            msgs.info("Not preparing a bad pixel mask")
        return


    def GetDispersionDirection(self, ind, det):
        '''Sets dispersion axis and transposes image as needed
        '''
        if self._argflag['trace']['disp']['direction'] is None:
            self._dispaxis = artrace.dispdir(self._msarc, dispwin=self._argflag['trace']['disp']['window'], mode=0)
        elif self._argflag['trace']['disp']['direction'] in [0,1]:
            self._dispaxis = int(self._argflag['trace']['disp']['direction'])
        else:
            msgs.error("The argument for the dispersion direction (trace+disp+direction)"+msgs.newline()+
                        "must be either:"+msgs.newline()+"  0 if the dispersion axis is predominantly along a row"+msgs.newline()+
                        "  1 if the dispersion axis is predominantly along a column")
        # Perform a check to warn the user if the longest axis is not equal to the dispersion direction
        if (self._msarc.shape[0] > self._msarc.shape[1]):
            if (self._dispaxis==1): msgs.warn("The dispersion axis is set to the shorter axis, is this correct?")
        else:
            if (self._dispaxis==0): msgs.warn("The dispersion axis is set to the shorter axis, is this correct?")

        ###############
        # Change dispersion direction and files if necessary
        # The code is programmed assuming dispaxis=0
        if (self._dispaxis==1):
            msgs.info("Transposing frames and keywords")
            # Flip the transpose switch
            self._transpose = True
            # Transpose the master bias frame
            if self._msbias_name is not None:
                self._msbias = self._msbias.T
                arsave.save_master(self, self._msbias, filename=self._msbias_name, frametype=self._argflag['reduce']['usebias'], ind=ind)
            # Transpose the master arc, and save it
            self._msarc = self._msarc.T
            arsave.save_master(self, self._msarc, filename=self._msarc_name, frametype='arc', ind=ind)
            # Transpose the bad pixel mask
            self._bpix = self._bpix.T
            # Transpose the amplifier sections frame
            self._ampsec = self._ampsec.T
            # Update the keywords of the fits files
            for i in xrange(len(self._fitsdict['naxis0'])):
                temp = self._fitsdict['naxis0'][i]
                self._fitsdict['naxis0'][i] = self._fitsdict['naxis1'][i]
                self._fitsdict['naxis1'][i] = temp
            # Change the user-specified (x,y) pixel sizes
            tmp = self._spect['det'][det-1]['xgap']
            self._spect['det'][det-1]['xgap'] = self._spect['det'][det-1]['ygap']
            self._spect['det'][det-1]['ygap'] = tmp
            self._spect['det'][det-1]['ysize'] = 1.0/self._spect['det'][det-1]['ysize']
            # Update the amplifier/data/overscan sections
            for i in xrange(self._spect['det'][det-1]['numamplifiers']):
                # Flip the order of the sections
                self._spect['det'][det-1]['ampsec{0:02d}'.format(i+1)] = self._spect['det'][det-1]['ampsec{0:02d}'.format(i+1)][::-1]
                self._spect['det'][det-1]['datasec{0:02d}'.format(i+1)] = self._spect['det'][det-1]['datasec{0:02d}'.format(i+1)][::-1]
                self._spect['det'][det-1]['oscansec{0:02d}'.format(i+1)] = self._spect['det'][det-1]['oscansec{0:02d}'.format(i+1)][::-1]
            # Change the user-specified (x,y) pixel sizes
            msgs.work("Transpose gain and readnoise frames")
            # Set the new dispersion axis
            self._dispaxis=0
        return


    def GetPixelLocations(self, det):
        if self._argflag['reduce']['locations'] is None:
            self._pixlocn = artrace.gen_pixloc(self, self._mstrace, det, gen=True)
        elif self._argflag['reduce']['locations'] in ["mstrace"]:
            self._pixlocn = arutils.gen_pixloc(self._spect,self._mstrace, det, gen=False)
        else:
            self._pixlocn = arload.load_master(self._argflag['run']['masterdir']+'/'+self._argflag['reduce']['locations'], frametype=None)
        return


    def MasterArc(self, sc, det):
        '''Generate MasterBias frame for a given detector
        Parameters:
        -----------
        sc: int 
          Index of the science frame
        det: int 
          Index of the detector
        '''
        if self._argflag['reduce']['usearc'] in ['arc']:
            msgs.info("Preparing a master arc frame")
            ind = self._spect['arc']['index'][sc]
            # Check if an *identical* master arc frame has already been produced
            self.SetFoundArc(False)
            for gb in xrange(len(self._done_arcs)):
                if np.array_equal(ind, self._done_arcs[gb]):
                    msgs.info("An identical master arc frame already exists")
                    msarc = arload.load_master(self._name_arcs[gb], frametype='arc')
                    self._tltprefix = os.path.splitext(os.path.basename(msarc_name))[0]
                    msarc_name = self._name_arcs[gb]
                    self.SetFoundArc(True)
            if not self._foundarc:
                # Load the arc frames
                frames = arload.load_frames(self, ind, det, frametype='arc', msbias=self._msbias)
                if self._argflag['reduce']['arcmatch'] > 0.0:
                    sframes = arsort.match_frames(self, frames, self._argflag['reduce']['arcmatch'], frametype='arc', satlevel=self._spect['det']['saturation']*self._spect['det']['nonlinear'])
                    subframes = np.zeros((frames.shape[0], frames.shape[1], len(sframes)))
                    numarr = np.array([])
                    for i in xrange(len(sframes)):
                        numarr = np.append(numarr,sframes[i].shape[2])
                        msarc = arcomb.comb_frames(sframes[i], det, spect=self._spect, frametype='arc', **self._argflag['arc']['comb'])
                        msarc_name = "{0:s}/{1:s}/sub-msarc{2:s}_{3:03d}_{4:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],self._spect["det"]["suffix"],len(self._done_arcs),i)
                        # Send the data away to be saved
                        arsave.save_master(self, msarc, filename=msarc_name, frametype='arc', ind=ind)
                        subframes[:,:,i]=msarc.copy()
                    del sframes
                    # Combine all sub-frames
                    msarc = arcomb.comb_frames(subframes, det, spect=self._spect, frametype='arc', weights=numarr, **self._argflag['arc']['comb'])
                    del subframes
                else:
                    msarc = arcomb.comb_frames(frames, det, spect=self._spect, frametype='arc', **self._argflag['arc']['comb'])
                del frames
                # Derive a suitable name for the master arc frame
                msarc_name = "{0:s}/{1:s}/msarc{2:s}_{3:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],self._spect["det"][det-1]["suffix"],len(self._done_arcs))
                self._tltprefix = os.path.splitext(os.path.basename(msarc_name))[0]
                # Send the data away to be saved
                arsave.save_master(self, msarc, filename=msarc_name, frametype='arc', ind=ind)
                # Store the files used and the master bias name in case it can be used during the later reduction processes
                self._done_arcs.append(ind)
                self._name_arcs.append(msarc_name)
        else:
            msarc_name = self._argflag['run']['masterdir']+'/'+self._argflag['reduce']['usearc']
            self._tltprefix = os.path.splitext(os.path.basename(msarc_name))[0]
            msarc=arload.load_master(msarc_name, frametype=None)
        return msarc, msarc_name


    def MasterBias(self, sc, det):
        '''Generate MasterBias frame for a given detector
        Parameters:
        -----------
        sc: int 
          Index of the science frame
        det: int 
          Index of the detector
        '''
        if self._argflag['reduce']['usebias'] in ['bias', 'dark']:
            msgs.info("Preparing a master {0:s} frame".format(self._argflag['reduce']['usebias']))
            if self._argflag['reduce']['usebias'] == 'bias':
                # Get all of the bias frames for this science frame
                ind = self._spect['bias']['index'][sc]
            elif self._argflag['reduce']['usebias'] == 'dark':
                # Get all of the dark frames for this science frame
                ind = self._spect['dark']['index'][sc]
            # Check if an *identical* master bias frame has already been produced
            found = False
            for gb in xrange(len(self._done_bias)):
                if np.array_equal(ind, self._done_bias[gb]):
                    msgs.info("An identical master {0:s} frame already exists.".format(self._argflag['reduce']['usebias']))
                    msbias = arload.load_master(self._name_bias[gb], frametype=self._argflag['reduce']['usebias'])
                    found = True
                    break
            if not found:
                # Load the Bias/Dark frames
                frames = arload.load_frames(self, ind, det, frametype=self._argflag['reduce']['usebias'], transpose=self._transpose)
                msbias = arcomb.comb_frames(frames, det, spect=self._spect, frametype=self._argflag['reduce']['usebias'], **self._argflag['bias']['comb'])
                # Derive a suitable name for the master bias frame
                msbias_name = "{0:s}/{1:s}/msbias{2:s}_{3:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],self._spect["det"][det-1]["suffix"],len(self._done_bias))
                # Send the data away to be saved
                arsave.save_master(self, msbias, filename=msbias_name, frametype=self._argflag['reduce']['usebias'], ind=ind)
                # Store the files used and the master bias name in case it can be used during the later reduction processes
                self._done_bias.append(ind)
                self._name_bias.append(msbias_name)
        elif self._argflag['reduce']['usebias'] == 'overscan':
            msbias = 'overscan'
            msbias_name = None
        elif self._argflag['reduce']['usebias'] == 'none':
            msgs.info("Not performing a bias/dark subtraction")
            msbias=None
            msbias_name = None
        else: # It must be the name of a file the user wishes to load
            msbias_name = os.getcwd()+self._argflag['run']['masterdir']+'/'+self._argflag['reduce']['usebias']
            msbias = arload.load_master(self._argflag['run']['masterdir']+'/'+self._argflag['reduce']['usebias'], frametype=None)
        return msbias, msbias_name


    def MasterFlatField(self, sc, det):
        if self._argflag['reduce']['flatfield']: # Only do it if the user wants to flat field
            ###############
            # Generate a master pixel flat frame
            if self._argflag['reduce']['useflat'] in ['pixflat', 'blzflat']:
                msgs.info("Preparing a master pixel flat frame with {0:s}".format(self._argflag['reduce']['useflat']))
                if self._argflag['reduce']['useflat'] == 'pixflat':
                    # Get all of the pixel flat frames for this science frame
                    ind = self._spect['pixflat']['index'][sc]
                elif self._argflag['reduce']['useflat'] == 'blzflat':
                    # Get all of the blzflat frames for this science frame
                    ind = self._spect['blzflat']['index'][sc]
                # Check if an *identical* master pixel flat frame has already been produced
                found = False
                for gb in xrange(len(self._done_flat)):
                    if np.array_equal(ind, self._done_flat[gb]):
                        msgs.info("An identical master flat frame already exists")
                        mspixflat_name = self._name_flat[gb]
                        mspixflat = arload.load_master(mspixflat_name, frametype='pixel flat')
                        found = True
                        break
                if not found:
                    # Load the frames for tracing
                    frames = arload.load_frames(self, ind, frametype='pixel flat', msbias=self._msbias, transpose=self._transpose)
                    if self._argflag['reduce']['flatmatch'] > 0.0:
                        sframes = arsort.match_frames(self, frames, self._argflag['reduce']['flatmatch'], frametype='pixel flat', satlevel=self._spect['det'][det-1]['saturation']*self._spect['det'][det-1]['nonlinear'])
                        subframes = np.zeros((frames.shape[0], frames.shape[1], len(sframes)))
                        numarr = np.array([])
                        for i in xrange(len(sframes)):
                            numarr = np.append(numarr,sframes[i].shape[2])
                            mspixflat = arcomb.comb_frames(sframes[i], spect=self._spect, frametype='pixel flat', **self._argflag['pixflat']['comb'])
                            mspixflat_name = "{0:s}/{1:s}/sub-msflat{2:s}_{3:03d}_{4:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],self._spect["det"][det-1]["suffix"],len(self._done_flat),i)
                            # Send the data away to be saved
                            arsave.save_master(self, mspixflat, filename=mspixflat_name, frametype='pixel flat', ind=ind)
                            subframes[:,:,i]=mspixflat.copy()
                        del sframes
                        # Combine all sub-frames
                        mspixflat = arcomb.comb_frames(subframes, spect=self._spect, frametype='pixel flat', weights=numarr, **self._argflag['pixflat']['comb'])
                        del subframes
                    else:
                        mspixflat = arcomb.comb_frames(frames, spect=self._spect, frametype='pixel flat', **self._argflag['pixflat']['comb'])
                    del frames
                    # Derive a suitable name for the master pixel flat frame
                    mspixflat_name = "{0:s}/{1:s}/msflat{2:s}_{3:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],self._spect["det"][det-1]["suffix"],len(self._done_flat))
                    # Send the data away to be saved
                    arsave.save_master(self, mspixflat, filename=mspixflat_name, frametype='pixel flat', ind=ind)
                    # Store the files used and the master bias name in case it can be used during the later reduction processes
                    self._done_flat.append(ind)
                    self._name_flat.append(mspixflat_name)
            else: # It must be the name of a file the user wishes to load
                mspixflat_name = self._argflag['run']['masterdir']+'/'+self._argflag['reduce']['usepixflat']
                mspixflat=arload.load_master(mspixflat_name, frametype=None)
            # Now that the combined, master flat field frame is loaded...
            # Normalize the flat field
            mspixflatnrm, msblaze = arproc.flatnorm(self, det, mspixflat, overpix=0, fname=os.path.basename(os.path.splitext(mspixflat_name)[0]))
        else:
            msgs.work("Pixel Flat arrays need to be generated when not flat fielding")
            msgs.bug("Blaze is currently undefined")
            mspixflat = np.ones_like(self._msarc)
            mspixflatnrm = np.ones_like(self._msarc)
            msblaze = None
            mspixflat_name = None
        return mspixflat, mspixflatnrm, msblaze, mspixflat_name


    def MasterTrace(self, sc, det):
        if self._argflag['reduce']['usetrace'] in ['trace', 'blzflat']:
            msgs.info("Preparing a master trace frame with {0:s}".format(self._argflag['reduce']['usetrace']))
            if self._argflag['reduce']['usetrace'] == 'trace':
                # Get all of the trace frames for this science frame
                ind = self._spect['trace']['index'][sc]
            elif self._argflag['reduce']['usetrace'] == 'blzflat':
                # Get all of the blzflat frames for this science frame
                ind = self._spect['blzflat']['index'][sc]
            # Check if an *identical* master trace frame has already been produced
            foundtrc = False
            for gb in xrange(len(self._done_flat)):
                if np.array_equal(ind, self._done_flat[gb]):
                    msgs.info("An identical master trace frame already exists")
                    mstrace_name = self._name_flat[gb]
                    mstrace = arload.load_master(self._name_flat[gb], frametype='trace')
                    self._trcprefix = os.path.splitext(os.path.basename(self._name_flat[gb]))[0]
                    foundtrc = True
                    break
            if not foundtrc:
                # Load the frames for tracing
                frames = arload.load_frames(self, ind, det, frametype='trace', msbias=self._msbias, trim=self._argflag['reduce']['trim'], transpose=self._transpose)
                if self._argflag['reduce']['flatmatch'] > 0.0:
                    sframes = arsort.match_frames(self, frames, self._argflag['reduce']['flatmatch'], frametype='trace', satlevel=self._spect['det'][det-1]['saturation']*self._spect['det'][det-1]['nonlinear'])
                    subframes = np.zeros((frames.shape[0], frames.shape[1], len(sframes)))
                    numarr = np.array([])
                    for i in xrange(len(sframes)):
                        numarr = np.append(numarr, sframes[i].shape[2])
                        mstrace = arcomb.comb_frames(sframes[i], spect=self._spect, frametype='trace', **self._argflag['trace']['comb'])
                        mstrace_name = "{0:s}/{1:s}/sub-msflat{2:s}_{3:03d}_{4:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],self._spect["det"][det-1]["suffix"],len(self._done_flat),i)
                        null = os.path.splitext(os.path.basename(mstrace_name))[0]
                        self._trcprefix = "{0:s}{1:s}".format(null,self._spect["det"][det-1]["suffix"])
                        # Send the data away to be saved
                        arsave.save_master(self, mstrace, filename=mstrace_name, frametype='trace', ind=ind)
                        subframes[:,:,i]=mstrace.copy()
                    del sframes
                    # Combine all sub-frames
                    mstrace = arcomb.comb_frames(subframes, det, spect=self._spect, frametype='trace', weights=numarr, **self._argflag['trace']['comb'])
                    del subframes
                else:
                    mstrace = arcomb.comb_frames(frames, det, spect=self._spect, frametype='trace', **self._argflag['trace']['comb'])
                del frames
                # Derive a suitable name for the master trace frame
                mstrace_name = "{0:s}/{1:s}/msflat{2:s}_{3:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],self._spect["det"][det-1]["suffix"],len(self._done_flat))
                self._trcprefix = os.path.splitext(os.path.basename(mstrace_name))[0]
                # Send the data away to be saved
                arsave.save_master(self, mstrace, filename=mstrace_name, frametype='trace', ind=ind)
                # Store the files used and the master bias name in case it can be used during the later reduction processes
                self._done_flat.append(ind)
                self._name_flat.append(mstrace_name)
        elif self._argflag['reduce']['usetrace'] == 'science':
            msgs.work("Tracing with a science frame is not implemented")
        else: # It must be the name of a file the user wishes to load
            mstrace_name = self._argflag['run']['masterdir']+'/'+self._argflag['reduce']['usetrace']
            mstrace=arload.load_master(mstrace_name, frametype=None)
            self._trcprefix = os.path.splitext(os.path.basename(mstrace_name))[0]
        return mstrace, mstrace_name


    def Setup(self):
        # Sort the data
        msgs.bug("Files and folders should not be deleted -- there should be an option to overwrite files automatically if they already exist, or choose to rename them if necessary")
        self._filesort = arsort.sort_data(self)
        # Write out the details of the sorted files
        if self._argflag['out']['sorted'] is not None: arsort.sort_write(self)
        # Match Science frames to calibration frames
        arsort.match_science(self)
        # If the user is only debugging, then exit now
        if self._argflag['run']['calcheck']:
            msgs.info("Calibration check complete. Change the 'calcheck' flag to continue with data reduction")
            sys.exit()
        # Make directory structure for different objects
        self._sci_targs = arsort.make_dirs(self)
        return


    # Bunch of Setters
    def SetFoundArc(self, bool): self._foundarc = bool


    ###################################
    # Reduction pipelines
    ###################################

    def ARMED(self):
        """
        Automatic Reduction & Modelling of Echelle Data
        """
        success = False
        # Insert series of reduction steps here
        return success


    def ARMLSD(self):
        """
        Automatic Reduction & Modelling of Long Slit Data
        """
        success = False
        # Sort data and match calibrations to science frames
        self.Setup()
        sci = self._filesort['science']
        numsci = np.size(sci)
        if numsci == 0:
            msgs.bug("What to do if no science frames are input? The calibrations should still be processed.")
            msgs.work("Maybe assume that each non-identical arc is a science frame, but force no science extraction")
            msgs.error("No science frames are input!!")
            numsci=1
            self._argflag['run']['preponly'] = True # Prepare the calibration files, but don't do the science frame reduction
        # Reduce each science frame entirely before moving to the next one
        for sc in xrange(numsci):
            ###############
            # First set the index for the science frame
            scidx = self._spect['science']['index'][sc]
            self._scidx = scidx[0]
            sciext_name_p, sciext_name_e = os.path.splitext(self._fitsdict['filename'][scidx[0]])
            # Now loop on Detectors
            for kk in xrange(self._spect['mosaic']['ndet']):
                det = kk + 1 # Detectors indexed from 1
                prefix = "{0:s}{1:s}".format(sciext_name_p,self._spect["det"][kk]["suffix"])
                ###############
                # Get amplifier sections
                self._ampsec = arproc.get_ampsec_trimmed(self, det, scidx[0])
                ###############
                # Generate master bias frame
                self._msbias, self._msbias_name = self.MasterBias(sc, det)
                ###############
                # Generate a bad pixel mask
                self.BadPixelMask(sc, det)
                ###############
                # Estimate gain and readout noise for the amplifiers
                msgs.work("Estimate Gain and Readout noise from the raw frames...")
                ###############
                # Generate a master arc frame
                self._msarc, self._msarc_name = self.MasterArc(sc, det)
                if self._bpix is None: self._bpix = np.zeros_like(self._msarc)
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
                # Determine the edges of the spectrum
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
                    msbg_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), self._argflag['run']['masterdir'], self._fitsdict['target'][scidx[0]], 0, "sky")
                    # Send the data away to be saved
                    arsave.save_master(self, bgframe, filename=msbg_name, frametype='sky background')
                    # Derive a suitable name for the sky-subtracted science frame
                    msscibg_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), self._argflag['run']['masterdir'], self._fitsdict['target'][scidx[0]], 0, "skysub")
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
                    msbg_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), self._argflag['run']['masterdir'], self._fitsdict['target'][scidx[0]], 0, "sky")
                    # Send the data away to be saved
                    arsave.save_master(self, bgframe, filename=msbg_name, frametype='sky background')
                    # Derive a suitable name for the sky-subtracted science frame
                    msscibg_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.fits".format(os.getcwd(), self._argflag['run']['masterdir'], self._fitsdict['target'][scidx[0]], 0, "skysub")
                    # Send the data away to be saved
                    arsave.save_master(self, sciframe-bgframe, filename=msscibg_name, frametype='sky subtracted science')
                    # Redetermine the variance frame based on the new sky model
                    varframe = arproc.variance_frame(self, det, sciframe, scidx[0], skyframe=bgframe)
                ###############
                # Determine the final trace of the science objects
                scitrace = artrace.trace_object(self, sciframe-bgframe, varframe, crmask)
                if scitrace is None:
                    msgs.info("Not performing extraction for science frame"+msgs.newline()+self._fitsdict['filename'][scidx[0]])
                    continue

                ###############
                # Extraction
                boxcar = arextract.boxcar(sciframe-bgframe, varframe, crmask, scitrace)
                msgs.error("UP TO HERE")
                ###############
                # Perform a velocity correction
                if self._argflag['reduce']['heliocorr'] == True:
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
                msgs.error("Implement flexure correction")

                ###############
                # Determine the wavelength scale (i.e. the wavelength of each pixel) to be used during the extraction
                msgs.info("Generating the array of extraction wavelengths")
                self._wavelength = arproc.get_wscale(self)

                msgs.error("INSERT EXTRACTION HERE")


        # Insert remaining reduction steps here
        return success

if __name__ == "__main__":
    prognm = sys.argv[0]
    debug = True
    quick = False

    # Load options from command line
    try:
        opt,arg=getopt.getopt(sys.argv[1:],'hqc:v:', ['help',
                                                 'quick'
                                                ])
    except getopt.GetoptError, err:
        msgs.error(err.msg)
        usage(prognm)
    for o,a in opt:
        if   o in ('-h', '--help')      : usage(argflag)
        elif o in ('-q', '--quick')     : quick = True
#		elif o in ('-c', '--cpus')      : argflag['run']['ncpus']     = a
#		elif o in ('-v', '--verbose')   : argflag['out']['verbose']   = int(a)

    if debug:
        argflag = arload.optarg(sys.argv, last_updated)
        ClassMain(argflag, quick=quick)
    else:
        try:
            argflag = arload.optarg(sys.argv, last_updated)
            ClassMain(argflag, quick=quick)
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
