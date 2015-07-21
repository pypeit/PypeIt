#!/usr/bin/env python
"""
Automatic Reduction and Modeling of Echelle Data
"""
last_updated='Last updated 16th November 2013'

# Import standard libraries
import os
import sys
import time
import copy
import signal
import warnings
import traceback
import numpy as np
# Import a Chi-Squared minimisation package
#from arcsmin import arfit
# Import useful ARMED definitions:
import arload
import arsave
import arsort
import armsgs as msgs
#import arfunc_base
# Import the ARMED reduction files
import ararc
#import arbias
#import arcalib
import arcomb
import arextract
import arextract_temp
#import arflat
import artrace
import arproc
import arvcorr

class ClassMain:

    def __init__(self, argflag, getinst=False):
        if getinst: return # Just get an instance

        # Set parameters
        self._argflag = argflag
        self._transpose = False   # Determine if the frames need to be transposed

        # First send all signals to messages to be dealt
        # with (i.e. someone hits ctrl+c)
        signal.signal(signal.SIGINT, msgs.signal_handler)

        # Ignore all warnings given by python
        warnings.resetwarnings()
        warnings.simplefilter("ignore")

        # Record the starting time
        self._tstart=time.time()

        # Load the Input file
        self._parlines, self._datlines, self._spclines = arload.load_input(self)

        # Load the Spectrograph settings
        self._spect = arload.load_spect(self)

        # Load any changes to the spectrograph settings
        self._spect = arload.load_spect(self, lines=self._spclines)

        # Load the important information from the fits headers
        self._fitsdict = arload.load_headers(self)

        # Load the list of standard stars
        self._standardStars = arload.load_standards(self)

        # Load the Data
        #arload.load_data(self)

        # Reduce the data!
        self.main()

####################################################################################
####################################################################################
####################################################################################

    # Reduce the data
    def main(self):
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
        # Store the name for the frames that have already been combined
        done_bias, name_bias = [], []
        done_flat, name_flat = [], []
        done_arcs, name_arcs = [], []
        # Reduce each science frame entirely before moving to the next one
        sci = self._filesort['science']
        sfx = self._spect["det"]["suffix"]
        numsci = np.size(sci)
        if numsci == 0:
            msgs.bug("What to do if no science frames are input? The calibrations should still be processed.")
            msgs.work("Maybe assume that each non-identical arc is a science frame, but force no science extraction")
            msgs.error("No science frames are input!!")
            numsci=1
            self._argflag['run']['preponly'] = True # Prepare the calibration files, but don't do the science frame reduction
        for sc in range(numsci):
            ###############
            # First set the index for the science frame
            scidx = self._spect['science']['index'][sc]
            sciext_name_p, sciext_name_e = os.path.splitext(self._fitsdict['filename'][scidx[0]])
            prefix = "{0:s}{1:s}".format(sciext_name_p,sfx)

            ###############
            # Generate master bias frame
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
                for gb in range(len(done_bias)):
                    if np.array_equal(ind, done_bias[gb]):
                        msgs.info("An identical master {0:s} frame already exists.".format(self._argflag['reduce']['usebias']))
                        msbias = arload.load_master(name_bias[gb], frametype=self._argflag['reduce']['usebias'])
                        found = True
                        break
                if not found:
                    # Load the Bias/Dark frames
                    frames = arload.load_frames(self, ind, frametype=self._argflag['reduce']['usebias'], transpose=self._transpose)
                    msbias = arcomb.comb_frames(frames, spect=self._spect, frametype=self._argflag['reduce']['usebias'], **self._argflag['bias']['comb'])
                    # Derive a suitable name for the master bias frame
                    msbias_name = "{0:s}/{1:s}/msbias{2:s}_{3:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],sfx,len(done_bias))
                    # Send the data away to be saved
                    arsave.save_master(self, msbias, filename=msbias_name, frametype=self._argflag['reduce']['usebias'], ind=ind)
                    # Store the files used and the master bias name in case it can be used during the later reduction processes
                    done_bias.append(ind)
                    name_bias.append(msbias_name)
            elif self._argflag['reduce']['usebias'] == 'overscan':
                msbias = 'overscan'
            elif self._argflag['reduce']['usebias'] == 'none':
                msgs.info("Not performing a bias/dark subtraction")
                msbias=None
            else: # It must be the name of a file the user wishes to load
                msbias_name = os.getcwd()+self._argflag['run']['masterdir']+'/'+self._argflag['reduce']['usebias']
                msbias = arload.load_master(self._argflag['run']['masterdir']+'/'+self._argflag['reduce']['usebias'], frametype=None)

            ###############
            # Generate a bad pixel mask
            if self._argflag['reduce']['badpix']:
                msgs.info("Preparing a bad pixel mask")
                # Get all of the bias frames for this science frame
                ind = self._spect['bias']['index'][sc]
                # Load the Bias frames
                frames = arload.load_frames(self, ind, frametype='bias', trim=False)
                tbpix = arcomb.comb_frames(frames, spect=self._spect, frametype='bias', **self._argflag['bias']['comb'])
                self._bpix = arproc.badpix(self,tbpix)
                del tbpix
            else:
                msgs.info("Not preparing a bad pixel mask")


            ###############
            # Estimate gain and readout noise for the amplifiers
            msgs.work("Estimate Gain and Readout noise...")

            ###############
            # Generate a master arc frame
            if self._argflag['reduce']['usearc'] in ['arc']:
                msgs.info("Preparing a master arc frame")
                ind = self._spect['arc']['index'][sc]
                # Check if an *identical* master arc frame has already been produced
                foundarc = False
                for gb in range(len(done_arcs)):
                    if np.array_equal(ind, done_arcs[gb]):
                        msgs.info("An identical master arc frame already exists")
                        msarc = arload.load_master(name_arcs[gb], frametype='arc')
                        tltprefix = os.path.splitext(os.path.basename(msarc_name))[0]
                        msarc_name = name_arcs[gb]
                        foundarc = True
                if not foundarc:
                    # Load the arc frames
                    frames = arload.load_frames(self, ind, frametype='arc', msbias=msbias)
                    if self._argflag['reduce']['arcmatch'] > 0.0:
                        sframes = arsort.match_frames(self, frames, self._argflag['reduce']['arcmatch'], frametype='arc', satlevel=self._spect['det']['saturation']*self._spect['det']['nonlinear'])
                        subframes = np.zeros((frames.shape[0], frames.shape[1], len(sframes)))
                        numarr = np.array([])
                        for i in range(len(sframes)):
                            numarr = np.append(numarr,sframes[i].shape[2])
                            msarc = arcomb.comb_frames(sframes[i], spect=self._spect, frametype='arc', **self._argflag['arc']['comb'])
                            msarc_name = "{0:s}/{1:s}/sub-msarc{2:s}_{3:03d}_{4:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],sfx,len(done_arcs),i)
                            # Send the data away to be saved
                            arsave.save_master(self, msarc, filename=msarc_name, frametype='arc', ind=ind)
                            subframes[:,:,i]=msarc.copy()
                        del sframes
                        # Combine all sub-frames
                        msarc = arcomb.comb_frames(subframes, spect=self._spect, frametype='arc', weights=numarr, **self._argflag['arc']['comb'])
                        del subframes
                    else:
                        msarc = arcomb.comb_frames(frames, spect=self._spect, frametype='arc', **self._argflag['arc']['comb'])
                    del frames
                    # Derive a suitable name for the master arc frame
                    msarc_name = "{0:s}/{1:s}/msarc{2:s}_{3:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],sfx,len(done_arcs))
                    tltprefix = os.path.splitext(os.path.basename(msarc_name))[0]
                    # Send the data away to be saved
                    arsave.save_master(self, msarc, filename=msarc_name, frametype='arc', ind=ind)
                    # Store the files used and the master bias name in case it can be used during the later reduction processes
                    done_arcs.append(ind)
                    name_arcs.append(msarc_name)
            else:
                msarc_name = self._argflag['run']['masterdir']+'/'+self._argflag['reduce']['usearc']
                tltprefix = os.path.splitext(os.path.basename(msarc_name))[0]
                msarc=arload.load_master(msarc_name, frametype=None)

            ###############
            # Determine the dispersion direction only on the first pass through
            if (sc==0):
                ###############
                # Guess dispersion direction
                if self._argflag['trace']['disp']['direction'] is None:
                    self._dispaxis = artrace.dispdir(msarc, dispwin=self._argflag['trace']['disp']['window'], mode=0)
                elif self._argflag['trace']['disp']['direction'] in [0,1]:
                    self._dispaxis = int(self._argflag['trace']['disp']['direction'])
                else:
                    msgs.error("The argument for the dispersion direction (trace+disp+direction)"+msgs.newline()+
                                "must be either:"+msgs.newline()+"  0 if the dispersion axis is predominantly along a row"+msgs.newline()+
                                "  1 if the dispersion axis is predominantly along a column")
                # Perform a check to warn the user if the longest axis is not equal to the dispersion direction
                if (msarc.shape[0] > msarc.shape[1]):
                    if (self._dispaxis==1): msgs.warn("The dispersion axis is set to the shorter axis, is this correct?")
                else:
                    if (self._dispaxis==0): msgs.warn("The dispersion axis is set to the shorter axis, is this correct?")

                ###############
                # Change dispersion direction and files if necessary
                # The code is programmed assuming dispaxis=0
                if (self._dispaxis==1):
                    msgs.info("Transposing frames")
                    # Flip the transpose switch
                    self._transpose = True
                    # Transpose the master bias frame
                    msbias = msbias.T
                    arsave.save_master(self, msbias, filename=msbias_name, frametype=self._argflag['reduce']['usebias'], ind=ind)
                    # Transpose the master arc, and save it
                    msarc = msarc.T
                    arsave.save_master(self, msarc, filename=msarc_name, frametype='arc', ind=ind)
                    # Transpose the bad pixel mask
                    self._bpix = self._bpix.T
                    # Change the user-specified (x,y) pixel sizes
                    tmp = self._spect['det']['xgap']
                    self._spect['det']['xgap'] = self._spect['det']['ygap']
                    self._spect['det']['ygap'] = tmp
                    self._spect['det']['ysize'] = 1.0/self._spect['det']['ysize']
                    # Update the amplifier/data/overscan sections
                    for i in range(self._spect['det']['numamplifiers']):
                        # Flip the order of the sections
                        self._spect['det']['ampsec{0:02d}'.format(i+1)] = self._spect['det']['ampsec{0:02d}'.format(i+1)][::-1]
                        self._spect['det']['datasec{0:02d}'.format(i+1)] = self._spect['det']['datasec{0:02d}'.format(i+1)][::-1]
                        self._spect['det']['oscansec{0:02d}'.format(i+1)] = self._spect['det']['oscansec{0:02d}'.format(i+1)][::-1]
                    # Change the user-specified (x,y) pixel sizes
                    msgs.work("Transpose gain and readnoise frames")
                    # Set the new dispersion axis
                    self._dispaxis=0

            ###############
            # Generate a master trace frame
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
                for gb in range(len(done_flat)):
                    if np.array_equal(ind, done_flat[gb]):
                        msgs.info("An identical master trace frame already exists")
                        mstrace_name = name_flat[gb]
                        mstrace = arload.load_master(name_flat[gb], frametype='trace')
                        trcprefix = os.path.splitext(os.path.basename(name_flat[gb]))[0]
                        foundtrc = True
                        break
                if not foundtrc:
                    # Load the frames for tracing
                    frames = arload.load_frames(self, ind, frametype='trace', msbias=msbias, trim=self._argflag['reduce']['trim'], transpose=self._transpose)
                    if self._argflag['reduce']['flatmatch'] > 0.0:
                        sframes = arsort.match_frames(self, frames, self._argflag['reduce']['flatmatch'], frametype='trace', satlevel=self._spect['det']['saturation']*self._spect['det']['nonlinear'])
                        subframes = np.zeros((frames.shape[0], frames.shape[1], len(sframes)))
                        numarr = np.array([])
                        for i in range(len(sframes)):
                            numarr = np.append(numarr,sframes[i].shape[2])
                            mstrace = arcomb.comb_frames(sframes[i], spect=self._spect, frametype='trace', **self._argflag['trace']['comb'])
                            mstrace_name = "{0:s}/{1:s}/sub-msflat{2:s}_{3:03d}_{4:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],sfx,len(done_flat),i)
                            null = os.path.splitext(os.path.basename(mstrace_name))[0]
                            trcprefix = "{0:s}{1:s}".format(null,sfx)
                            # Send the data away to be saved
                            arsave.save_master(self, mstrace, filename=mstrace_name, frametype='trace', ind=ind)
                            subframes[:,:,i]=mstrace.copy()
                        del sframes
                        # Combine all sub-frames
                        mstrace = arcomb.comb_frames(subframes, spect=self._spect, frametype='trace', weights=numarr, **self._argflag['trace']['comb'])
                        del subframes
                    else:
                        mstrace = arcomb.comb_frames(frames, spect=self._spect, frametype='trace', **self._argflag['trace']['comb'])
                    del frames
                    # Derive a suitable name for the master trace frame
                    mstrace_name = "{0:s}/{1:s}/msflat{2:s}_{3:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],sfx,len(done_flat))
                    trcprefix = os.path.splitext(os.path.basename(mstrace_name))[0]
                    # Send the data away to be saved
                    arsave.save_master(self, mstrace, filename=mstrace_name, frametype='trace', ind=ind)
                    # Store the files used and the master bias name in case it can be used during the later reduction processes
                    done_flat.append(ind)
                    name_flat.append(mstrace_name)
            elif self._argflag['reduce']['usetrace'] == 'science':
                msgs.work("Tracing with a science frame is not implemented")
            else: # It must be the name of a file the user wishes to load
                mstrace_name = self._argflag['run']['masterdir']+'/'+self._argflag['reduce']['usetrace']
                mstrace=arload.load_master(mstrace_name, frametype=None)
                trcprefix = os.path.splitext(os.path.basename(mstrace_name))[0]


            ###############
            # Generate an array that provides the physical pixel locations on the detector
            if self._argflag['reduce']['locations'] is None:
                self._pixlocn = artrace.gen_pixloc(self,mstrace,gen=True)
            elif self._argflag['reduce']['locations'] in ["mstrace"]:
                self._pixlocn = arutils.gen_pixloc(self._spect,mstrace,gen=False)
            else:
                self._pixlocn = arload.load_master(self._argflag['run']['masterdir']+'/'+self._argflag['reduce']['locations'], frametype=None)


            ###############
            # Determine the location of and trace the orders for the trace frame
#			if foundtrc:
#				try:
#					# Load the order locations from file
#					self._lordloc, self._rordloc = arload.load_ordloc(mstrace_name)
#				except:
#					# If this fails, rederive the order locations
#					self._lordloc, self._rordloc = artrace.trace_longslit(self, mstrace, prefix=prefix, trcprefix=trcprefix)
#					# Save the traces
#					arsave.save_ordloc(self, mstrace_name)
#			else:
#				# First time encountered this set of trace frames --> derive the order locations
#				self._lordloc, self._rordloc = artrace.trace_longslit(self, mstrace, prefix=prefix, trcprefix=trcprefix)
#				# Save the traces
#				arsave.save_ordloc(self, mstrace_name)
            # Take the edges of the (trimmed) data to be the order edges
            self._lordloc = self._pixlocn[:,0,1].reshape((self._pixlocn.shape[0],1))
            self._rordloc = self._pixlocn[:,-1,1].reshape((self._pixlocn.shape[0],1))
            # Convert physical trace into a pixel trace
            msgs.info("Converting physical trace locations to nearest pixel")
            self._pixcen  = artrace.phys_to_pix(0.5*(self._lordloc+self._rordloc), self._pixlocn, self._dispaxis, 1-self._dispaxis)
            self._pixwid  = (self._rordloc-self._lordloc).mean(0).astype(np.int)
            self._lordpix = artrace.phys_to_pix(self._lordloc, self._pixlocn, self._dispaxis, 1-self._dispaxis)
            self._rordpix = artrace.phys_to_pix(self._rordloc, self._pixlocn, self._dispaxis, 1-self._dispaxis)


            ###############
            # Prepare the Flat-field frame
            if self._argflag['reduce']['flatfield']:
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
                    for gb in range(len(done_flat)):
                        if np.array_equal(ind, done_flat[gb]):
                            msgs.info("An identical master flat frame already exists")
                            mspixflat_name = name_flat[gb]
                            mspixflat = arload.load_master(mspixflat_name, frametype='pixel flat')
                            found = True
                            break
                    if not found:
                        # Load the frames for tracing
                        frames = arload.load_frames(self, ind, frametype='pixel flat', msbias=msbias)
                        if self._argflag['reduce']['flatmatch'] > 0.0:
                            sframes = arsort.match_frames(self, frames, self._argflag['reduce']['flatmatch'], frametype='pixel flat', satlevel=self._spect['det']['saturation']*self._spect['det']['nonlinear'])
                            subframes = np.zeros((frames.shape[0], frames.shape[1], len(sframes)))
                            numarr = np.array([])
                            for i in range(len(sframes)):
                                numarr = np.append(numarr,sframes[i].shape[2])
                                mspixflat = arcomb.comb_frames(sframes[i], spect=self._spect, frametype='pixel flat', **self._argflag['pixflat']['comb'])
                                mspixflat_name = "{0:s}/{1:s}/sub-msflat{2:s}_{3:03d}_{4:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],sfx,len(done_flat),i)
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
                        mspixflat_name = "{0:s}/{1:s}/msflat{2:s}_{3:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],sfx,len(done_flat))
                        # Send the data away to be saved
                        arsave.save_master(self, mspixflat, filename=mspixflat_name, frametype='pixel flat', ind=ind)
                        # Store the files used and the master bias name in case it can be used during the later reduction processes
                        done_flat.append(ind)
                        name_flat.append(mspixflat_name)
                else: # It must be the name of a file the user wishes to load
                    mspixflat_name = self._argflag['run']['masterdir']+'/'+self._argflag['reduce']['usepixflat']
                    mspixflat=arload.load_master(mspixflat_name, frametype=None)
                # Now that the combined, master flat field frame is loaded...
                # Normalize the flat field
                mspixflatnrm, msblaze = arproc.flatnorm(self, mspixflat, overpix=0, fname=os.path.basename(os.path.splitext(mspixflat_name)[0]))

            ###############
            # Derive the order tilts
            if foundarc:
                try:
                    # Load the order locations from file
                    self._tilts, self._satmask = arload.load_tilts(msarc_name)
                except:
                    # If this fails, rederive the order tilts
                    self._tilts, self._satmask = artrace.model_tilt(self, msarc, prefix=prefix, trcprefix=trcprefix, tltprefix=tltprefix) # NOTE: self._tilts = tan(tilt angle), where "tilt angle" is the angle between (1) the line representing constant wavelength and (2) the column of pixels that is most closely parallel with the spatial direction of the slit.
                    # Save the traces
                    arsave.save_tilts(self, msarc_name)
            else:
                # First time encountered this set of arc frames --> derive the order tilts
                self._tilts, self._satmask = artrace.model_tilt(self, msarc, prefix=prefix, trcprefix=trcprefix, tltprefix=tltprefix) # NOTE: self._tilts = tan(tilt angle), where "tilt angle" is the angle between (1) the line representing constant wavelength and (2) the column of pixels that is most closely parallel with the spatial direction of the slit.
                # Save the tilts
                msgs.error("OK?")
                arsave.save_tilts(self, msarc_name)

            ###############
            # Extract wavelength arrays and save the extracted lines
            msarcext_name_p, msarcext_name_e = os.path.splitext(msarc_name)
            msarcext_name = msarcext_name_p + "_ext1d" + msarcext_name_e
            if self._argflag['arc']['load']['extracted'] == True:
                if os.path.exists(msarcext_name):
                    self._arcext = arload.load_master(msarcext_name, frametype=None)
                else:
                    msgs.warn("Master 1D arc extraction does not exist for file:"+msgs.newline()+msarc_name)
                    msgs.info("Reduction will continue by performing 1D arc extraction")
                    self._arcext, extprops = arextract.extract(self, msarc, mode='mean', frametype='arc')
                    arsave.save_master(self, self._arcext, filename=msarcext_name, frametype='1D arc extraction')
            else:
                self._arcext, extprops = arextract.extract(self, msarc, mode='mean', frametype='arc')
                # Send the data away to be saved
                arsave.save_master(self, self._arcext, filename=msarcext_name, frametype='1D arc extraction')


            ###############
            # Perform wavelength calibration
            if self._argflag['arc']['load']['calibrated'] == True:
                if os.path.exists(somename):
                    self._waveids = arload.load_orderdetail(somename, frametype=None)
                else:
                    msgs.warn("Master arc calibration does not exist for file:"+msgs.newline()+msarc_name)
                    msgs.info("Reduction will continue by calibrating 1D arc extraction")
                    self._waveids = ararc.calibrate(self)
                    arsave.save_orderdetail(self, filename=somename)
            else:
                msarc_id_p, msarc_id_e = os.path.splitext(msarc_name)
                fname = msarc_id_p + "_id" + msarc_id_e
                #pixels=arload.waveids(fname)
                pixels=None
                self._waveids = ararc.calibrate(self,msarc_name,pixtmp=pixels, prefix=prefix)
                # Send the data away to be saved
                msgs.work("Think about the best way to save the order location/tilt/calibration details")
                msgs.bug("Order wavelength IDs not saved")
                #arsave.save_orderdetail(self, filename=somename)

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
            # Determine the wavelength scale (i.e. the wavelength of each pixel) to be used during the extraction
            msgs.info("Generating the array of extraction wavelengths")
            self._wavelength = arproc.get_wscale(self)


            ###############
            # Check if the user only wants to prepare the calibrations
            msgs.info("All calibration frames have been prepared")
            if self._argflag['run']['preponly']:
                msgs.info("If you would like to continue with the reduction,"+msgs.newline()+"disable the run+preponly command")
                continue

            ###############
            # Load the science frame and from this generate a Poisson error frame
            sciframe = arload.load_frames(self, scidx, frametype='science', msbias=msbias)
            sciframe = sciframe[:,:,0]
            # Convert ADUs to electrons
            sciframe *= self._spect['det']['gain']
            varframe = arproc.variance_frame(slf, sciframe, scidx[0])
            snframe = arproc.sn_frame(self, sciframe, scidx[0])

            ###############
            # Subtract off the scattered light from the image
            msgs.work("Scattered light subtraction is not yet implemented...")
            """
            Do we want to do a scattered light subtraction? Probably not, since
            these photons contribute to the noise (error) frame. If a scattered
            light subtraction is implemented, make sure that the scattered light
            frame is passed to the "arproc.error_frame" procedure below when
            calculating the error frame.
            """


            ###############
            # Flat field the science frame
            if self._argflag['reduce']['flatfield']:
                msgs.info("Flat fielding the science frame")
                sciframe, errframe = arproc.flatfield(self, sciframe, mspixflatnrm, snframe=snframe)
            else:
                msgs.info("Not performing a flat field calibration")
                errframe = np.zeros_like(sciframe)
                wnz = np.where(snframe>0.0)
                errframe[wnz] = sciframe[wnz]/snframe[wnz]

            ###############
            # Flat field the science frame
            if self._argflag['reduce']['bgsubtraction']:
                msgs.info("Preparing a sky background frame from the science frame")
                scibgsub, bgframe = arproc.background_subtraction(self, sciframe, errframe)
                # Derive a suitable name for the master sky background frame
                msbg_name = "{0:s}/{1:s}/{2:s}_{3:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],self._fitsdict['target'][scidx[0]],"sky")
                # Send the data away to be saved
                arsave.save_master(self, bgframe, filename=msbg_name, frametype='sky background', ind=ind)
                # Derive a suitable name for the sky-subtracted science frame
                msscibg_name = "{0:s}/{1:s}/{2:s}_{3:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],self._fitsdict['target'][scidx[0]],"skysub")
                # Send the data away to be saved
                arsave.save_master(self, scibgsub, filename=msscibg_name, frametype='sky subtracted science', ind=ind)

#			msgs.error("There is a problem with the 2D spectrum fit  -  the mask flag (-999999.9) is not being correctly loaded in after being extracted."+
#						"You either need to find a way to load this value correctly, or save the order width in the header. This affects line 150 or arextract")

            #msgs.work("Temporary -- until extraction algorithm is improved")
            #sciext1D, scierr1D, skyext1D, skyerr1D = arextract_temp.optimalextract(self, sciframe, prefix=prefix)

            extdim = "ext1D"
            sciext_name = "{0:s}/{1:s}/{2:s}/{3:s}_{4:s}.fits".format(os.getcwd(),self._argflag['run']['scidir'],self._fitsdict['target'][scidx[0]],prefix,extdim)
            # Send the data away to be saved
            arsave.save_extraction(self, sciext1D, scidx, scierr=scierr1D, filename=sciext_name, frametype='science extraction', wave=self._sciwav, sky=skyext1D, skyerr=skyerr1D)
            msgs.work("Really need to work out the blaze/normalized flat when there are more than one amplifier")
            msgs.error("Up to here")
            continue













            ###############
            # Extract and save the science frame
            # DEPRECATED -> self._sciext, self._sciwav = arextract.extract(self, sciframe, mode=self._argflag["science"]["extraction"]["method"], frametype='science', wave=self._wavelength)
#			msgs.error("Havn't completed this yet -- also need to calculate sciwav in arextract before returning")
            # Derive a suitable name for the extracted science frame
            if self._argflag["science"]["extraction"]["method"].lower() in ["2d"]:
                extdim = "ext2D"
            else:
                extdim = "ext1D"
            sciext_name = "{0:s}/{1:s}/{2:s}/{3:s}_{4:s}.fits".format(os.getcwd(),self._argflag['run']['scidir'],self._fitsdict['target'][scidx[0]],prefix,extdim)
            if self._argflag['science']['load']['extracted'] == True:
                if os.path.exists(sciext_name):
                    msgs.work("Loaded extraction arrays are sometimes off by the machine precision (1 part in 10^16)")
                    """
                    This can be fixed by generating a wavelength scale using the pixsize from the header and finding the closest match to lambda_0.
                    """
                    self._sciext, self._sciwav, extprops = arload.load_extraction(sciext_name, frametype="Science", wave=True)
                else:
                    msgs.warn("A file containing the science extraction does not exist for file:"+msgs.newline()+sciext_name)
                    msgs.info("Reduction will continue by performing science extraction")
                    self._sciext, extprops, self._sciwav = arextract.extract(self, sciframe, error=errframe, mode=self._argflag["science"]["extraction"]["method"], frametype='Science', wave=self._wavelength)
                    # Send the data away to be saved
                    arsave.save_extraction(self, self._sciext, scidx, filename=sciext_name, frametype='science extraction', wave=self._sciwav, extprops=extprops)
            else:
                self._sciext, extprops, self._sciwav = arextract.extract(self, sciframe, error=errframe, mode=self._argflag["science"]["extraction"]["method"], frametype='Science', wave=self._wavelength)
                # Send the data away to be saved
                arsave.save_extraction(self, self._sciext, scidx, filename=sciext_name, frametype='science extraction', wave=self._sciwav, extprops=extprops)

            ###############
            # Continue with the next science target if the user only wanted a 1D extraction
            if self._argflag["science"]["extraction"]["method"].lower() not in ["2d"]:
                continue

            ###############
            # Otherwise, fit the 2D extracted science frame, and save the results
            msgs.work("The error frame is currently not self-consistently extracted from the science frame"+msgs.newline()+
                        "Instead, an error frame is generated from the extracted science frame")
            #scierr = arproc.error_frame_postext(self, self._sciext, scidx[0])
            sciext1D, scierr1D, skyext1D, skyerr1D = arextract.fit_spec2D(self, self._sciext, extprops, prefix=prefix)

            extdim = "ext1D"
            sciext_name = "{0:s}/{1:s}/{2:s}/{3:s}_{4:s}.fits".format(os.getcwd(),self._argflag['run']['scidir'],self._fitsdict['target'][scidx[0]],prefix,extdim)
            # Send the data away to be saved
            arsave.save_extraction(self, sciext1D, scidx, scierr=scierr1D, filename=sciext_name, frametype='science extraction', wave=self._sciwav, sky=skyext1D, skyerr=skyerr1D)

            #msgs.error("EXITING NORMALLY...")


            #msgs.error("A list of things to do:")
            # Take a median of the best-fitting sigma for the arc lines and refine the arc line centroids
            # Divide both the science frame and blaze flat frame by the pixel flat frame.
            # Do a PCA on the slit profile of (blaze divided by pixel) flat frame.
            # Apply this result to the science frame when extracting science object. (This will ensure that the gradient or slit profile doesn't affect the object profile).



# WHAT YOU NEED
# 1. Several bias frames (or a region for overscan)
# 2. Several flat field frames (ideally both pixel and blaze flats)
# 3. Several arc frames (with a slit length the same as you flat field)
# 4. Several trace frames (with a slit length less than the flat, and ideally the same length as the seeing)
# 5. A science frame!
#
# THINGS TO THINK ABOUT
# 1. What if a chip has two (or more!) amplifiers  ---  Make the user specify the number of amplifiers and the CCD regions the amplifiers apply to?
#
# THE DATA REDUCTION PROCESS
#
# Generate master bias frame (or identify the o/scan region)
# Estimate gain and readout noise for the amplifiers
# Subtract master bias from all individual trace frames (you might as well subtract the bias from all frames except the arcs)
# Generate several sub-master trace frames (combine all frames with the same instrument setup, weighted by counts --- perhaps if you weight by counts you should get the best result here --- you can go straight to a master trace frame.)
# Estimate dispersion direction (make sure all setups agree, if not, take the most common, or take the largest dispersion one [where the dispersion is measured with std(pixel[1:]-pixel[:-1]) ])
# Determine the location of (and trace) the orders for each submaster trace frame
# Find overlap in results from submaster traces
# PCA to derive the order locations and extrapolate to all possible orders on the chip
# Hence derive the inter-order positions (Not sure if this is necessary now that I'm weighting by counts.)
# Generate master trace frame (where you take the highest S/N sub-master)
# Generate Master flat frame (do the exact same as what you do with the trace)
# Generate Master arc frame (again, you should weight by counts, but leave saturated pixels in)
# Using the location of the orders from the trace frame, derive the slit profile from the master flat field for all orders.
# Extract the arc spectrum for each order (taking into account the slit profile and pixel separation on the chip)
# Automatically fit the Arc spectrum with a polynomial and a series of gaussians (keep adding Gaussians until some chi-squared threshold is reached)
# Transform the wavelength scale onto the 2D master flat-field frame with the user-specified pixel velocity size.
# Extract orders from master flat-field.
# Normalise the orders in the flat-field (i.e. go back to the 2D flat-field frame and for each wavelength bin, divide all pixels in this bin by the value derived in the previous step (times the fraction of this pixel in the appropriate wavelength bin).
# Divide the science frame by the normalised flat-field frame.
# Transform the wavelength scale onto the 2D science frame, and extract the data.
#
##############################################

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

if __name__ == "__main__":
    debug = True
    if debug:
        argflag = arload.optarg(sys.argv, last_updated)
        ClassMain(argflag)
    else:
        try:
            argflag = arload.optarg(sys.argv, last_updated)
            ClassMain(argflag)
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

