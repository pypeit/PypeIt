""" Class for book-keeping the reduction process
"""
from __future__ import absolute_import, division, print_function

import sys
import numpy as np
from astropy.time import Time
import datetime
from matplotlib.backends.backend_pdf import PdfPages
# Import PYPIT routines
from pypit import artrace
from pypit import arload
from pypit import arcomb
from pypit import arflux
from pypit import arlris
from pypit import armasters
from pypit import armsgs
from pypit import arparse
from pypit import arproc
from pypit import arsort
from pypit import arutils

try:
    from xastropy.xutils import xdebug as debugger
except ImportError:
    import pdb as debugger

# Logging and settings
msgs = armsgs.get_logger()
argflag = arparse.get_argflag().__dict__['_argflag']
spect = arparse.get_spect().__dict__['_spect']


class ScienceExposure:
    """
    A Science Exposure class that carries all information for a given science exposure
    """

    def __init__(self, snum, fitsdict, do_qa=True):

        #############################
        # Set some universal parameters
        self._transpose = False   # Determine if the frames need to be transposed

        # Set indices used for frame combination
        self._idx_sci = spect['science']['index'][snum]
        self._idx_arcs = spect['arc']['index'][snum]
        self._idx_trace = spect['trace']['index'][snum]
        self._idx_std = spect['standard']['index'][snum]
        if argflag['bias']['useframe'] == 'bias': self._idx_bias = spect['bias']['index'][snum]
        elif argflag['bias']['useframe'] == 'dark':  self._idx_bias = spect['dark']['index'][snum]
        else: self._idx_bias = []
        if argflag['reduce']['trace']['useframe'] == 'trace': self._idx_trace = spect['trace']['index'][snum]
        elif argflag['reduce']['trace']['useframe'] == 'slitflat': self._idx_trace = spect['slitflat']['index'][snum]
        else: self._idx_trace = []
        if argflag['reduce']['flatfield']['useframe'] == 'pixelflat': self._idx_flat = spect['pixelflat']['index'][snum]
        elif argflag['reduce']['flatfield']['useframe'] == 'slitflat': self._idx_flat = spect['slitflat']['index'][snum]
        else: self._idx_flat = []

        # Set the base name and extract other names that will be used for output files
        self.SetBaseName(fitsdict)

        # Initialize the QA for this science exposure
        qafn = "{0:s}/QA_{1:s}.pdf".format(argflag['run']['directory']['qa'], self._basename)
        if do_qa:
            self._qa = PdfPages(qafn)

        # Initialize Variables
        ndet = spect['mosaic']['ndet']
        self._nonlinear = [spect['det{0:02d}'.format(det+1)]['saturation'] *
                           spect['det{0:02d}'.format(det+1)]['nonlinear']
                           for det in range(ndet)]
        self._dispaxis = None  # Which direction is the predominant spectral (dispersion) axis
        self._nspec    = [None for all in range(ndet)]   # Number of spectral pixels
        self._nspat    = [None for all in range(ndet)]   # Number of spatial pixels
        self._ampsec   = [None for all in range(ndet)]   # Locations of the amplifiers on each detector
        self._pixlocn  = [None for all in range(ndet)]   # Physical locations of each pixel on the detector
        self._lordloc  = [None for all in range(ndet)]   # Array of slit traces (left side) in physical pixel coordinates
        self._rordloc  = [None for all in range(ndet)]   # Array of slit traces (left side) in physical pixel coordinates
        self._pixcen   = [None for all in range(ndet)]   # Central slit traces in apparent pixel coordinates
        self._pixwid   = [None for all in range(ndet)]   # Width of slit (at each row) in apparent pixel coordinates
        self._lordpix  = [None for all in range(ndet)]   # Array of slit traces (left side) in apparent pixel coordinates
        self._rordpix  = [None for all in range(ndet)]   # Array of slit traces (right side) in apparent pixel coordinates
        self._tilts    = [None for all in range(ndet)]   # Array of spectral tilts at each position on the detector
        self._tiltpar  = [None for all in range(ndet)]   # Dict parameters for tilt fitting
        self._satmask  = [None for all in range(ndet)]   # Array of Arc saturation streaks
        self._arcparam = [None for all in range(ndet)]   # Dict guiding wavelength calibration
        self._wvcalib  = [None for all in range(ndet)]   #
        self._resnarr  = [None for all in range(ndet)]   # Resolution array
        # Initialize the Master Calibration frames
        self._bpix = [None for all in range(ndet)]          # Bad Pixel Mask
        self._msarc = [None for all in range(ndet)]         # Master Arc
        self._mswave = [None for all in range(ndet)]         # Master Wavelength image
        self._msbias = [None for all in range(ndet)]        # Master Bias
        self._msrn = [None for all in range(ndet)]          # Master ReadNoise image
        self._mstrace = [None for all in range(ndet)]       # Master Trace
        self._mspixelflat = [None for all in range(ndet)]     # Master pixel flat
        self._mspixelflatnrm = [None for all in range(ndet)]  # Normalized Master pixel flat
        self._msblaze = [None for all in range(ndet)]       # Blaze function
        self._msstd = [{} for all in range(ndet)]           # Master Standard dict
        # Initialize the Master Calibration frame names
        self._msarc_name = [None for all in range(ndet)]      # Master Arc Name
        self._msbias_name = [None for all in range(ndet)]     # Master Bias Name
        self._mstrace_name = [None for all in range(ndet)]    # Master Trace Name
        self._mspixelflat_name = [None for all in range(ndet)]  # Master Pixel Flat Name
        # Initialize the science, variance, and background frames
        self._sciframe = [None for all in range(ndet)]
        self._rawvarframe = [None for all in range(ndet)]    # Variance based on detected counts + RN
        self._modelvarframe = [None for all in range(ndet)]  # Variance from sky and object models
        self._bgframe = [None for all in range(ndet)]
        self._scimask = [None for all in range(ndet)]        # Mask (1=Bad pix; 2=CR)
        self._scitrace = [None for all in range(ndet)]
        self._specobjs = [None for all in range(ndet)]
        # Initialize some extraction products
        self._ext_boxcar = [None for all in range(ndet)]
        self._ext_optimal = [None for all in range(ndet)]
        return

    def SetBaseName(self, fitsdict):
        """
        Set the base name that is used for all outputs

        Parameters
        ----------
        fitsdict : dict
          Contains relevant information from fits header files
        """
        scidx = self._idx_sci[0]
        try:
            "T" in fitsdict['date'][scidx]
        except IndexError:
            debugger.set_trace()
        if "T" in fitsdict['date'][scidx]:
            tbname = fitsdict['date'][scidx]
        else:
            # Not ideal, but convert MJD into a date+time
            timval = Time(fitsdict['time'][scidx]/24.0, scale='tt', format='mjd')
            tbname = timval.isot
        try:
            tval = datetime.datetime.strptime(tbname, '%Y-%m-%dT%H:%M:%S.%f')
        except ValueError:
            tval = datetime.datetime.strptime(tbname, '%Y-%m-%dT%H:%M:%S')
        self._inst_name = spect['mosaic']['camera']
        self._target_name = fitsdict['target'][self._idx_sci[0]].replace(" ", "")
        self._basename = self._target_name+'_'+self._inst_name+'_'+ \
                         datetime.datetime.strftime(tval, '%Y%b%dT') + \
                         tbname.split("T")[1].replace(':','')
        return

    ###################################
    # Reduction procedures
    ###################################

    def BadPixelMask(self, fitsdict, det):
        """
        Generate Bad Pixel Mask for a given detector

        Parameters
        ----------
        fitsdict : dict
          Contains relevant information from fits header files
        det : int
          Index of the detector

        Returns
        -------
        boolean : bool
          Should other ScienceExposure classes be updated?
        """
        bpix = None
        if argflag['reduce']['badpix'] == 'bias':
            if argflag['reduce']['masters']['reuse']:
                # Attempt to load the Master Frame
                bpix_name = armasters.master_name(argflag['run']['directory']['master'],
                                                'badpix', argflag['reduce']['masters']['setup'])
                try:
                    bpix, head = arload.load_master(bpix_name, frametype="badpix")
                except IOError:
                    msgs.warn("No MasterBadPix frame found {:s}".format(bpix_name))
                else:
                    argflag['reduce']['masters']['loaded'].append('badpix'+argflag['reduce']['masters']['setup'])
            if 'badpix'+argflag['reduce']['masters']['setup'] not in argflag['reduce']['masters']['loaded']:
                msgs.info("Preparing a bad pixel mask")
                # Get all of the bias frames for this science frame
                if len(self._idx_bias) == 0:
                    msgs.warn("No bias frames available to determine bad pixel mask")
                    msgs.info("Not preparing a bad pixel mask")
                    #self._bpix = None
                    return False
                # Load the Bias frames
                bpix = arproc.badpix(self, det, self.GetMasterFrame('bias', det))
        else:
            # Instrument dependent
            if argflag['run']['spectrograph'] in ['lris_red']:
                bpix = arlris.bpm(self, 'red', fitsdict, det)
            else:
                msgs.info("Not preparing a bad pixel mask")
                return False
        self.SetFrame(self._bpix, bpix, det)
        del bpix
        return True

    def GetDispersionDirection(self, fitsdict, det, scidx):
        """ Set the dispersion axis.
        If necessary, transpose frames and adjust information as needed

        Parameters
        ----------
        fitsdict : dict
          Contains relevant information from fits header files
        det : int
          Index of the detector
        scidx : int
          Index of frame

        Returns
        -------
        fitsdict : dict
          Updates to the input fitsdict
        """
        if argflag['trace']['dispersion']['direction'] is None:
            self._dispaxis = artrace.dispdir(self._msarc[det-1], dispwin=argflag['trace']['dispersion']['window'], mode=0)
        elif argflag['trace']['dispersion']['direction'] in [0, 1]:
            self._dispaxis = int(argflag['trace']['dispersion']['direction'])
        else:
            msgs.error("The argument for the dispersion direction (trace+disp+direction)"+msgs.newline() +
                       "must be either:"+msgs.newline()+"  0 if the dispersion axis is predominantly along a row" +
                       msgs.newline() + "  1 if the dispersion axis is predominantly along a column")
        # Perform a check to warn the user if the longest axis is not equal to the dispersion direction
        if self._msarc[det-1].shape[0] > self._msarc[det-1].shape[1]:
            if self._dispaxis == 1: msgs.warn("The dispersion axis is set to the shorter axis, is this correct?")
        else:
            if self._dispaxis == 0: msgs.warn("The dispersion axis is set to the shorter axis, is this correct?")

        ###############
        # Change dispersion direction and files if necessary
        # The code is programmed assuming dispaxis=0
        if self._dispaxis == 1:
            msgs.info("Transposing frames and keywords")
            # Flip the transpose switch
            self._transpose = True
            # Transpose the master bias frame
            if self._msbias[det-1] is not None:
                if type(self._msbias[det-1]) is str: pass  # Overscan sub - change the oscansec parameters below
                elif type(self._msbias[det-1]) is np.ndarray:
                    if 'bias'+argflag['reduce']['masters']['setup'] not in argflag['reduce']['masters']['loaded']:
                        self.SetMasterFrame(self._msbias[det-1].T, 'bias', det)
            # Transpose the master arc, and save it
            if 'arc'+argflag['reduce']['masters']['setup'] not in argflag['reduce']['masters']['loaded']:
                self.SetMasterFrame(self._msarc[det-1].T, 'arc', det)
            # Transpose the bad pixel mask
            if self._bpix[det-1] is not None:
                if 'badpix'+argflag['reduce']['masters']['setup'] not in argflag['reduce']['masters']['loaded']:
                    self.SetFrame(self._bpix, self._bpix[det-1].T, det)
            # Transpose the amplifier sections frame
            self.SetFrame(self._ampsec, self._ampsec[det-1].T, det)
            # Update the keywords of the fits files
            temp = fitsdict['naxis0'][scidx]
            fitsdict['naxis0'][scidx] = fitsdict['naxis1'][scidx]
            fitsdict['naxis1'][scidx] = temp
            # Change the user-specified (x,y) pixel sizes
            tmp = spect['det'][det-1]['xgap']
            spect['det'][det-1]['xgap'] = spect['det'][det-1]['ygap']
            spect['det'][det-1]['ygap'] = tmp
            spect['det'][det-1]['ysize'] = 1.0/spect['det'][det-1]['ysize']
            # Update the amplifier/data/overscan sections
            for i in range(spect['det'][det-1]['numamplifiers']):
                # Flip the order of the sections
                spect['det'][det-1]['ampsec{0:02d}'.format(i+1)] = spect['det'][det-1]['ampsec{0:02d}'.format(i+1)][::-1]
                spect['det'][det-1]['datasec{0:02d}'.format(i+1)] = spect['det'][det-1]['datasec{0:02d}'.format(i+1)][::-1]
                spect['det'][det-1]['oscansec{0:02d}'.format(i+1)] = spect['det'][det-1]['oscansec{0:02d}'.format(i+1)][::-1]
            # Change the user-specified (x,y) pixel sizes
            msgs.work("Transpose gain and readnoise frames")
            # Set the new dispersion axis
            self._dispaxis = 0
        else:
            msgs.info("Not transposing")
        # Set the number of spectral and spatial pixels
        self._nspec[det-1], self._nspat[det-1] = self._msarc[det-1].shape

    def GetPixelLocations(self, det):
        """
        Generate or load the physical location of each pixel

        Parameters
        ----------
        det : int
          Index of the detector
        """
        if argflag['reduce']['pixel']['locations'] is None:
            self.SetFrame(self._pixlocn, artrace.gen_pixloc(self._dispaxis, self._mstrace[det-1], det, gen=True), det)
        elif argflag['reduce']['pixel']['locations'] in ["mstrace"]:
            self.SetFrame(self._pixlocn, artrace.gen_pixloc(self._dispaxis, self._mstrace[det-1], det, gen=False), det)
        else:
            mname = argflag['run']['directory']['master']+'/'+argflag['reduce']['pixel']['locations']
            self.SetFrame(self._pixlocn, arload.load_master(mname, frametype=None), det)
        return

    def MasterArc(self, fitsdict, det):
        """
        Generate Master Arc frame for a given detector

        Parameters
        ----------
        fitsdict : dict
          Contains relevant information from fits header files
        det : int
          Index of the detector

        Returns
        -------
        boolean : bool
          Should other ScienceExposure classes be updated?
        """

        if self._msarc[det-1] is not None:
            msgs.info("An identical master arc frame already exists")
            return False
        if argflag['arc']['useframe'] in ['arc']:
            # Attempt to load the Master Frame
            if argflag['reduce']['masters']['reuse']:
                msarc_name = armasters.master_name(argflag['run']['directory']['master'],
                                                    'arc', argflag['reduce']['masters']['setup'])
                try:
                    msarc, head = arload.load_master(msarc_name, frametype="arc")
                except IOError:
                    msgs.warn("No MasterArc frame found {:s}".format(msarc_name))
                else:
                    self._transpose = head['transp']
                    if self._transpose:  # Need to setup for flipping
                        argflag['trace']['dispersion']['direction'] = 1
                    else:
                        argflag['trace']['dispersion']['direction'] = 0
                    # Append as loaded
                    argflag['reduce']['masters']['loaded'].append('arc'+argflag['reduce']['masters']['setup'])
            if 'arc'+argflag['reduce']['masters']['setup'] not in argflag['reduce']['masters']['loaded']:
                msgs.info("Preparing a master arc frame")
                ind = self._idx_arcs
                # Load the arc frames
                frames = arload.load_frames(fitsdict, ind, det, frametype='arc', msbias=self._msbias[det-1])
                if argflag['arc']['combine']['match'] > 0.0:
                    sframes = arsort.match_frames(frames, argflag['arc']['combine']['match'], msgs, frametype='arc',
                                                  satlevel=spect['det']['saturation']*spect['det']['nonlinear'])
                    subframes = np.zeros((frames.shape[0], frames.shape[1], len(sframes)))
                    numarr = np.array([])
                    for i in range(len(sframes)):
                        numarr = np.append(numarr, sframes[i].shape[2])
                        msarc = arcomb.comb_frames(sframes[i], det,
                                                   frametype='arc', **argflag['arc']['combine'])
                        # Send the data away to be saved
                        subframes[:,:,i] = msarc.copy()
                    del sframes
                    # Combine all sub-frames
                    msarc = arcomb.comb_frames(subframes, det, frametype='arc',
                                               weights=numarr, **argflag['arc']['combine'])
                    del subframes
                else:
                    msarc = arcomb.comb_frames(frames, det, frametype='arc', **argflag['arc']['combine'])
                del frames
            # # Derive a suitable name for the master arc frame
            # msarc_name = "{0:s}/{1:s}/msarc{2:s}_{3:03d}.fits".format(os.getcwd(),argflag['run']['directory']['master'],spect["det"][det-1]["suffix"],len(self._done_arcs))
            # self._tltprefix = os.path.splitext(os.path.basename(msarc_name))[0]
            # # Send the data away to be saved
            # arsave.save_master(self, msarc, filename=msarc_name, frametype='arc', ind=ind)
            # # Store the files used and the master bias name in case it can be used during the later reduction processes
            # self._done_arcs.append(ind)
            # self._name_arcs.append(msarc_name)
        else:
            msarc_name = argflag['run']['directory']['master']+'/'+argflag['arc']['useframe']
            msarc, head = arload.load_master(msarc_name, frametype=None)
        # Set and then delete the Master Arc frame
        self.SetMasterFrame(msarc, "arc", det)
        del msarc
        return True

    def MasterBias(self, fitsdict, det):
        """
        Generate Master Bias frame for a given detector

        Parameters
        ----------
        fitsdict : dict
          Contains relevant information from fits header files
        det : int
          Index of the detector

        Returns
        -------
        boolean : bool
          Should other ScienceExposure classes be updated?
        """

        # If the master bias is already made, use it
        if self._msbias[det-1] is not None:
            msgs.info("An identical master {0:s} frame already exists".format(argflag['bias']['useframe']))
            return False
        elif argflag['bias']['useframe'] in ['bias', 'dark']:
            # Load from hard-drive?
            if argflag['reduce']['masters']['reuse']:
                # Attempt to load the Master Frame
                msbias_name = armasters.master_name(argflag['run']['directory']['master'],
                                                    'bias', argflag['reduce']['masters']['setup'])
                try:
                    msbias, head = arload.load_master(msbias_name, frametype="bias")
                except IOError:
                    msgs.warn("No MasterBias frame found {:s}".format(msbias_name))
                else:
                    argflag['reduce']['masters']['loaded'].append('bias'+argflag['reduce']['masters']['setup'])
            if 'bias'+argflag['reduce']['masters']['setup'] not in argflag['reduce']['masters']['loaded']:
                msgs.info("Preparing a master {0:s} frame".format(argflag['bias']['useframe']))
                # Get all of the bias frames for this science frame
                ind = self._idx_bias
                # Load the Bias/Dark frames
                frames = arload.load_frames(fitsdict, ind, det, frametype=argflag['bias']['useframe'], transpose=self._transpose)
                msbias = arcomb.comb_frames(frames, det, frametype=argflag['bias']['useframe'], **argflag['bias']['combine'])
                del frames
        elif argflag['bias']['useframe'] == 'overscan':
            self.SetMasterFrame('overscan', "bias", det, mkcopy=False)
            return False
        elif argflag['bias']['useframe'] == 'none':
            msgs.info("Not performing a bias/dark subtraction")
            self.SetMasterFrame(None, "bias", det, mkcopy=False)
            return False
        else: # It must be the name of a file the user wishes to load
            msbias_name = argflag['run']['directory']['master']+'/'+argflag['bias']['useframe']
            msbias, head = arload.load_master(msbias_name, frametype="bias")
            argflag['reduce']['masters']['loaded'].append('bias')
        # Set and then delete the Master Bias frame
        self.SetMasterFrame(msbias, "bias", det)

        del msbias
        return True

    '''
    def MasterRN(self, fitsdict, det):
        """
        Generate Master ReadNoise frame for a given detector
        Parameters
        ----------
        fitsdict : dict
          Contains relevant information from fits header files
        det : int
          Index of the detector
        Returns
        -------
        boolean : bool
          Should other ScienceExposure classes be updated?
        """

        # If the master bias is already made, use it
        if self._msrn[det-1] is not None:
            msgs.info("An identical master ReadNoise frame already exists")
            return False
        msrn = np.zeros((self._nspec[det-1], self._nspat[det-1]))
        # Systems with multiple amps will need help here
        if spect['det'][det-1]['numamplifiers'] > 1:
            msgs.work("Readnoise needs to be updated for multiple amps")
        # Set
        rnoise = spect['det'][det-1]['ronoise'] #+ (0.5*spect['det'][det-1]['gain'])**2
        msrn[:] = rnoise
        # Save
        self.SetMasterFrame(msrn, "readnoise", det)
        del msrn
        return True
    '''

    def MasterFlatField(self, fitsdict, det):
        """
        Generate Master Flat-field frame for a given detector

        Parameters
        ----------
        fitsdict : dict
          Contains relevant information from fits header files
        det : int
          Index of the detector

        Returns
        -------
        boolean : bool
          Should other ScienceExposure classes be updated?
        """

        if argflag['reduce']['flatfield']['perform']:  # Only do it if the user wants to flat field
        # If the master pixelflat is already made, use it
            if self._mspixelflat[det-1] is not None:
                msgs.info("An identical master pixelflat frame already exists")
                if self._mspixelflatnrm[det-1] is None:
                    # Normalize the flat field
                    msgs.info("Normalizing the pixel flat")
                    mspixelflatnrm, msblaze = arproc.flatnorm(self, det, self.GetMasterFrame("pixelflat", det),
                                                            overpix=0, plotdesc="Blaze function")
                    self.SetFrame(self._msblaze, msblaze, det)
                    self.SetMasterFrame(mspixelflatnrm, "normpixelflat", det)
                return False
            ###############
            # Generate a master pixel flat frame
            if argflag['reduce']['flatfield']['useframe'] in ['pixelflat', 'slitflat']:
                if argflag['reduce']['masters']['reuse']:
                    # Attempt to load the Master Frame
                    msflat_name = armasters.master_name(argflag['run']['directory']['master'],
                                                    'normpixelflat', argflag['reduce']['masters']['setup'])
                    try:
                        mspixelflatnrm, head = arload.load_master(msflat_name, frametype="normpixelflat")
                    except IOError:
                        msgs.warn("No MasterFlatField frame found {:s}".format(msflat_name))
                    else:
                        argflag['reduce']['masters']['loaded'].append('normpixelflat'+argflag['reduce']['masters']['setup'])
                        mspixelflat = mspixelflatnrm
                if 'normpixelflat'+argflag['reduce']['masters']['setup'] not in argflag['reduce']['masters']['loaded']:
                    msgs.info("Preparing a master pixel flat frame with {0:s}".format(argflag['reduce']['flatifeld']['useframe']))
                    # Get all of the pixel flat frames for this science frame
                    ind = self._idx_flat
                    # Load the frames for tracing
                    frames = arload.load_frames(fitsdict, ind, det, frametype='pixel flat',
                                                msbias=self._msbias[det-1], transpose=self._transpose)
                    if argflag['pixelflat']['combine']['match'] > 0.0:
                        sframes = arsort.match_frames(frames, argflag['pixelflat']['combine']['match'],
                                                      frametype='pixel flat', satlevel=self._nonlinear)
                        subframes = np.zeros((frames.shape[0], frames.shape[1], len(sframes)))
                        numarr = np.array([])
                        for i in range(len(sframes)):
                            numarr = np.append(numarr, sframes[i].shape[2])
                            mspixelflat = arcomb.comb_frames(sframes[i], det, frametype='pixel flat',
                                                           **argflag['pixelflat']['combine'])
                            subframes[:,:,i] = mspixelflat.copy()
                        del sframes
                        # Combine all sub-frames
                        mspixelflat = arcomb.comb_frames(subframes, det, frametype='pixel flat',
                                                       weights=numarr, **argflag['pixelflat']['combine'])
                        del subframes
                    else:
                        mspixelflat = arcomb.comb_frames(frames, det, frametype='pixel flat',
                                                       **argflag['pixelflat']['combine'])
                    del frames
                    # Apply gain (instead of ampsec scale)
                    mspixelflat *= arproc.gain_frame(self, det)
                    # Normalize the flat field
                    mspixelflatnrm, msblaze = arproc.flatnorm(self, det, mspixelflat, overpix=0, plotdesc="Blaze function")
                    self.SetFrame(self._msblaze, msblaze, det)
            else:  # It must be the name of a file the user wishes to load
                mspixelflat_name = armasters.user_master_name(argflag['run']['directory']['master'],
                                                            argflag['reduce']['flatfield']['useframe'])
                mspixelflatnrm, head = arload.load_master(mspixelflat_name, exten=det, frametype=None)
                mspixelflat = mspixelflatnrm
            # Now that the combined, master flat field frame is loaded...
        else:
            msgs.work("Pixel Flat arrays need to be generated when not flat fielding")
            msgs.bug("Blaze is currently undefined")
            mspixelflat = np.ones_like(self._msarc)
            mspixelflatnrm = np.ones_like(self._msarc)
        # Set Master Frames
        self.SetMasterFrame(mspixelflat, "pixelflat", det)
        self.SetMasterFrame(mspixelflatnrm, "normpixelflat", det)
        return True

    def MasterTrace(self, fitsdict, det):
        """
        Generate Master Trace frame for a given detector

        Parameters
        ----------
        fitsdict : dict
          Contains relevant information from fits header files
        det : int
          Index of the detector

        Returns
        -------
        boolean : bool
          Should other ScienceExposure classes be updated?
        """

        # If the master trace is already made, use it
        if self._mstrace[det-1] is not None:
            msgs.info("An identical master trace frame already exists")
            return False
        if argflag['reduce']['trace']['useframe'] in ['trace', 'slitflat']:
            if argflag['reduce']['masters']['reuse']:
                # Attempt to load the Master Frame
                mstrace_name = armasters.master_name(argflag['run']['directory']['master'],
                                                   'trace', argflag['reduce']['masters']['setup'])
                try:
                    mstrace, head = arload.load_master(mstrace_name, frametype="trace")
                except IOError:
                    msgs.warn("No MasterTrace frame found {:s}".format(mstrace_name))
                else:
                    # Extras
                    lordloc, _ = arload.load_master(mstrace_name, frametype="trace", exten=1)
                    rordloc, _ = arload.load_master(mstrace_name, frametype="trace", exten=2)
                    pixcen, _ = arload.load_master(mstrace_name, frametype="trace", exten=3)
                    pixwid, _ = arload.load_master(mstrace_name, frametype="trace", exten=4)
                    lordpix, _ = arload.load_master(mstrace_name, frametype="trace", exten=5)
                    rordpix, _ = arload.load_master(mstrace_name, frametype="trace", exten=6)
                    self.SetFrame(self._lordloc, lordloc, det)
                    self.SetFrame(self._rordloc, rordloc, det)
                    self.SetFrame(self._pixcen, pixcen.astype(np.int), det)
                    self.SetFrame(self._pixwid, pixwid.astype(np.int), det)
                    self.SetFrame(self._lordpix, lordpix.astype(np.int), det)
                    self.SetFrame(self._rordpix, rordpix.astype(np.int), det)
                    #
                    argflag['reduce']['masters']['loaded'].append('trace'+argflag['reduce']['masters']['setup'])
            if 'trace'+argflag['reduce']['masters']['setup'] not in argflag['reduce']['masters']['loaded']:
                msgs.info("Preparing a master trace frame with {0:s}".format(argflag['reduce']['trace']['useframe']))
                ind = self._idx_trace
                # Load the frames for tracing
                frames = arload.load_frames(fitsdict, ind, det, frametype='trace', msbias=self._msbias[det-1],
                                            trim=argflag['reduce']['trim'], transpose=self._transpose)
                if argflag['trace']['combine']['match'] > 0.0:
                    sframes = arsort.match_frames(frames, argflag['trace']['combine']['match'], msgs, frametype='trace', satlevel=spect['det'][det-1]['saturation']*spect['det'][det-1]['nonlinear'])
                    subframes = np.zeros((frames.shape[0], frames.shape[1], len(sframes)))
                    numarr = np.array([])
                    for i in range(len(sframes)):
                        numarr = np.append(numarr, sframes[i].shape[2])
                        mstrace = arcomb.comb_frames(sframes[i], det, frametype='trace', **argflag['trace']['combine'])
                        subframes[:,:,i] = mstrace.copy()
                    del sframes
                    # Combine all sub-frames
                    mstrace = arcomb.comb_frames(subframes, det, frametype='trace', weights=numarr, **argflag['trace']['combine'])
                    del subframes
                else:
                    mstrace = arcomb.comb_frames(frames, det, frametype='trace', **argflag['trace']['combine'])
                del frames
        elif argflag['reduce']['trace']['useframe'] == 'science':
            msgs.error("Tracing with a science frame is not yet implemented")
        else: # It must be the name of a file the user wishes to load
            mstrace_name = argflag['run']['directory']['master']+'/'+argflag['reduce']['trace']['useframe']
            mstrace, head = arload.load_master(mstrace_name, frametype=None)
            debugger.set_trace()  # NEED TO LOAD EXTRAS AS ABOVE
        # Set and then delete the Master Trace frame
        self.SetMasterFrame(mstrace, "trace", det)
        del mstrace
        return True

    def MasterWave(self, fitsdict, sc, det):
        """
        Generate Master Wave frame for a given detector

        Parameters
        ----------
        fitsdict : dict
          Contains relevant information from fits header files
        det : int
          Index of the detector

        Returns
        -------
        boolean : bool
          Should other ScienceExposure classes be updated?
        """
        if self._mswave[det-1] is not None:
            msgs.info("An identical master arc frame already exists")
            return False
        # Attempt to load the Master Frame
        if argflag['reduce']['masters']['reuse']:
            mswave_name = armasters.master_name(argflag['run']['directory']['master'],
                                                'wave', argflag['reduce']['masters']['setup'])
            try:
                mswave, head = arload.load_master(mswave_name, frametype="arc")
            except IOError:
                msgs.warn("No MasterWave frame found {:s}".format(mswave_name))
            else:
                argflag['reduce']['masters']['loaded'].append('wave'+argflag['reduce']['masters']['setup'])
        if 'wave'+argflag['reduce']['masters']['setup'] not in argflag['reduce']['masters']['loaded']:
            msgs.info("Preparing a master wave frame")
            wv_calib = self._wvcalib[det-1]
            mswave = arutils.func_val(wv_calib['fitc'], self._tilts[det-1], wv_calib['function'], minv=wv_calib['fmin'], maxv=wv_calib['fmax'])
        # Set and then delete the Master Arc frame
        self.SetMasterFrame(mswave, "wave", det)
        del mswave
        return True

    def MasterWaveCalib(self, fitsdict, sc, det):
        """
        Generate Master 1D Wave Solution (down slit center)

        Parameters
        ----------
        fitsdict : dict
          Contains relevant information from fits header files
        det : int
          Index of the detector

        Returns
        -------
        boolean : bool
          Should other ScienceExposure classes be updated?
        """
        from pypit import ararc

        if self._wvcalib[det-1] is not None:
            msgs.info("An identical master wave calib frame already exists")
            return False
        # Attempt to load the Master Frame
        if argflag['reduce']['masters']['reuse']:
            mswv_calib_name = armasters.master_name(argflag['run']['directory']['master'],
                                                    'wave_calib', argflag['reduce']['masters']['setup'])
            try:
                wv_calib = arload.load_master(mswv_calib_name, frametype="wv_calib")
            except IOError:
                msgs.warn("No MasterWave1D data found {:s}".format(mswv_calib_name))
            else:
                argflag['reduce']['masters']['loaded'].append('wave_calib'+argflag['reduce']['masters']['setup'])
        if 'wave_calib'+argflag['reduce']['masters']['setup'] not in argflag['reduce']['masters']['loaded']:
            # Setup arc parameters (e.g. linelist)
            arcparam = ararc.setup_param(self, sc, det, fitsdict)
            self.SetFrame(self._arcparam, arcparam, det)
            ###############
            # Extract arc and identify lines
            wv_calib = ararc.simple_calib(self, det)
        # Set
        if wv_calib is not None:
            self.SetFrame(self._wvcalib, wv_calib, det)
        del wv_calib
        return True

    def MasterStandard(self, scidx, fitsdict):
        """
        Generate Master Standard frame for a given detector
        and generates a sensitivity function
        Currently only uses first standard star exposure
        Currently takes brightest source on the mosaic

        Parameters
        ----------
        fitsdict : dict
          Contains relevant information from fits header files

        Returns
        -------
        boolean : bool
        """

        if len(self._msstd[0]) != 0:
            msgs.info("Using existing standard frame")
            return False
        #
        msgs.info("Preparing the standard")
        # Get all of the pixel flat frames for this science frame
        ind = self._idx_std
        msgs.warn("Taking only the first standard frame for now")
        ind = [ind[0]]
        # Extract
        all_specobj = []
        for kk in range(spect['mosaic']['ndet']):
            det = kk+1
            # Load the frame(s)
#            set_trace()
            frame = arload.load_frames(fitsdict, ind, det, frametype='standard',
                                       msbias=self._msbias[det-1],
                                       transpose=self._transpose)
#            msgs.warn("Taking only the first standard frame for now")
#            ind = ind[0]
            sciframe = frame[:, :, 0]
            # Save RA/DEC
            if kk == 0:
                self._msstd[det-1]['RA'] = fitsdict['ra'][ind[0]]
                self._msstd[det-1]['DEC'] = fitsdict['dec'][ind[0]]
            #debugger.set_trace()
            arproc.reduce_frame(self, sciframe, ind[0], fitsdict, det, standard=True)

            #
            all_specobj += self._msstd[det-1]['spobjs']
#        debugger.set_trace()
        # If standard, generate a sensitivity function
        sensfunc = arflux.generate_sensfunc(self, scidx, all_specobj, fitsdict)
        # Set the sensitivity function
        self.SetMasterFrame(sensfunc, "sensfunc", None, mkcopy=False)
        return True

    def Setup(self):

        # Sort the data
        msgs.bug("Files and folders should not be deleted -- there should be an option to overwrite files automatically if they already exist, or choose to rename them if necessary")
        self._filesort = arsort.sort_data(self)
        # Write out the details of the sorted files
        if argflag['output']['sorted'] is not None: arsort.sort_write(self)
        # Match Science frames to calibration frames
        arsort.match_science(self)
        # If the user is only debugging, then exit now
        if argflag['run']['calcheck']:
            msgs.info("Calibration check complete. Change the 'calcheck' flag to continue with data reduction")
            sys.exit()
        # Make directory structure for different objects
        self._sci_targs = arsort.make_dirs(self)
        return

    # Setters
    @staticmethod
    def SetFrame(toarray, value, det, mkcopy=True):
        if mkcopy:
            toarray[det-1] = value.copy()
        else:
            toarray[det-1] = value
        return

    def SetMasterFrame(self, frame, ftype, det, mkcopy=True):
        """ Set the Master Frame
        Parameters
        ----------
        frame : object
        ftype : str
          frame type
        det : int
          Detector index
        mkcopy

        Returns
        -------

        """
        if det is not None:
            det -= 1
        if mkcopy:
            cpf = frame.copy()
        else:
            cpf = frame
        # Set the frame
        if ftype == "arc": self._msarc[det] = cpf
        elif ftype == "wave": self._mswave[det] = cpf
        elif ftype == "bias": self._msbias[det] = cpf
        elif ftype == "readnoise": self._msrn[det] = cpf
        elif ftype == "normpixelflat": self._mspixelflatnrm[det] = cpf
        elif ftype == "pixelflat": self._mspixelflat[det] = cpf
        elif ftype == "trace": self._mstrace[det] = cpf
        elif ftype == "standard": self._msstd[det] = cpf
        elif ftype == "sensfunc": self._sensfunc = cpf
        else:
            msgs.bug("I could not set master frame of type: {0:s}".format(ftype))
            msgs.error("Please contact the authors")
        return

    # Getters
    @staticmethod
    def GetFrame(getarray, det, mkcopy=True):
        if mkcopy:
            return getarray[det-1].copy()
        else:
            return getarray[det-1]

    def GetMasterFrame(self, ftype, det, mkcopy=True):

        det -= 1
        # Get the frame
        if mkcopy:
            if ftype == "arc": return self._msarc[det].copy()
            elif ftype == "wave": return self._mswave[det].copy()
            elif ftype == "bias": return self._msbias[det].copy()
            elif ftype == "normpixelflat": return self._mspixelflatnrm[det].copy()
            elif ftype == "pixelflat": return self._mspixelflat[det].copy()
            elif ftype == "trace": return self._mstrace[det].copy()
            elif ftype == "standard": return mkcopy.copy(self._msstd[det])
            elif ftype == "sensfunc": return mkcopy.copy(self._sensfunc)
            else:
                msgs.bug("I could not get master frame of type: {0:s}".format(ftype))
                msgs.error("Please contact the authors")
        else:
            if ftype == "arc": return self._msarc[det]
            elif ftype == "wave": return self._mswave[det]
            elif ftype == "bias": return self._msbias[det]
            elif ftype == "normpixelflat": return self._mspixelflatnrm[det]
            elif ftype == "pixelflat": return self._mspixelflat[det]
            elif ftype == "trace": return self._mstrace[det]
            elif ftype == "standard": return self._msstd[det]
            elif ftype == "sensfunc": return self._sensfunc
            else:
                msgs.bug("I could not get master frame of type: {0:s}".format(ftype))
                msgs.error("Please contact the authors")
        return None

    def update_sci_pixmask(self, det, mask_pix, mask_type):
        """ Update the binary pixel mask for a given science frame

        Parameters
        ----------
        det : int
        mask_pix : ndarray
          Image of pixels to mask (anything >0)
        mask_type : str
          Type of masked pixel
            BadPix = 1
            CR = 2
        """
        mask_dict = dict(BadPix=1, CR=2)
        if mask_type not in mask_dict.keys():
            msgs.error("Bad pixel mask type")
        # Find pixels to mask
        mask = np.where(mask_pix > 0)
        if len(mask[0]) == 0:
            return
        # Update those that need it
        prev_val = self._scimask[det-1][mask]
        upd = np.where((prev_val % 2**(mask_dict[mask_type]+1)) < 2**(mask_dict[mask_type]))[0]
        if len(upd) > 0:
            self._scimask[det-1][mask[0][upd], mask[1][upd]] += 2**mask_dict[mask_type]
        # Return
        return


    def __repr__(self):
        return ('<{:s}: frame={:s} target={:s} index={:d}>'.format(
                self.__class__.__name__, self._basename, self._target_name,
                self._idx_sci[0]))
