import numpy as np
# Import PYPIT routines
import armsgs as msgs
from matplotlib.backends.backend_pdf import PdfPages
import arsort
import arload

class ScienceExposure:

    def __init__(self, scidx, spect):
        """
        A Science Exposure class that carries all information for a given science exposure
        """
        #############################
        # Set some universal parameters
        self._argflag = argflag   # Arguments and Flags
        self._transpose = False   # Determine if the frames need to be transposed
        # Arrays to store the name for the frames that have already been combined
        self._done_bias, self._name_bias = [], []
        self._done_flat, self._name_flat = [], []
        self._done_arcs, self._name_arcs = [], []

        # Initialize the QA for this science exposure
        qafn = "testname.pdf"
        self._qa = pp = PdfPages(qafn)
        # Set indices used for frame combination
        self._idx_sci = scidx
        self._idx_arcs = None
        self._idx_bias = None
        self._idx_flat = None
        # Initialize Variables
        self._dispaxis = None  # Which direction is the predominant spectral (dispersion) axis
        self._pixlocn = None   # Physical locations of each pixel on the detector
        self._lordloc = None   # Array of slit traces (left side) in physical pixel coordinates
        self._rordloc = None   # Array of slit traces (left side) in physical pixel coordinates
        self._pixcen  = None   # Central slit traces in apparent pixel coordinates
        self._pixwid  = None   # Width of slit (at each row) in apparent pixel coordinates
        self._lordpix = None   # Array of slit traces (left side) in apparent pixel coordinates
        self._rordpix = None   # Array of slit traces (right side) in apparent pixel coordinates
        self._tilts   = None   # Array of spectral tilts at each position on the detector
        self._satmask = None   # Array of Arc saturation streaks
        self._arcparam = None  #
        self._wvcalib = None   #
        self._resnarr = None   # Resolution array
        # Initialize the Master Calibration frames
        self._bpix = None          # Bad Pixel Mask
        self._msarc = None         # Master Arc
        self._msbias = None        # Master Bias
        self._mstrace = None       # Master Trace
        self._mspixflat = None     # Master pixel flat
        self._mspixflatnrm = None  # Normalized Master pixel flat
        self._msblaze = None       # Blaze function
        # Initialize the Master Calibration frame names
        self._msarc_name = None      # Master Arc Name
        self._msbias_name = None     # Master Bias Name
        self._mstrace_name = None    # Master Trace Name
        self._mspixflat_name = None  # Master Pixel Flat Name
        # Initialize the science, variance, and background frames
        self._sciframe = None
        self._varframe = None
        self._bgframe  = None
        # Initialize some extraction products
        self._ext_boxcar = None
        self._ext_optimal = None


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
        '''Generate MasterArc frame for a given detector
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
                    msarc_name = self._name_arcs[gb]
                    msarc = arload.load_master(msarc_name, frametype='arc')
                    self._tltprefix = os.path.splitext(os.path.basename(msarc_name))[0]
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

    # Setters
    def SetMasterFrame(self, frame, ftype):
        if ftype == "arc": self._msarc = frame.copy()
        elif ftype == "badpix": self._bpix = frame.copy()
        elif ftype == "bias": self._msbias = frame.copy()
        elif ftype == "normpixflat": self._mspixflatnrm = frame.copy()
        elif ftype == "pixflat": self._mspixflat = frame.copy()
        elif ftype == "trace": self._mstrace = frame.copy()
        else:
            msgs.bug("I could not set master frame of type: {0:s}".format(ftype))
            msgs.error("Please contact the authors")
        return

    def SetFoundArc(self, bool): self._foundarc = bool

    # Getters
    def GetMasterFrame(self, ftype):
        if ftype == "arc": return self._msarc.copy()
        elif ftype == "badpix": return self._bpix.copy()
        elif ftype == "bias": return self._msbias.copy()
        elif ftype == "normpixflat": return self._mspixflatnrm.copy()
        elif ftype == "pixflat": return self._mspixflat.copy()
        elif ftype == "trace": return self._mstrace.copy()
        else:
            msgs.bug("I could not set master frame of type: {0:s}".format(ftype))
        return None

