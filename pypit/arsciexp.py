""" Class for book-keeping the reduction process
"""
from __future__ import (absolute_import, division, print_function, unicode_literals)

import sys
import numpy as np
from astropy.time import Time
from matplotlib.backends.backend_pdf import PdfPages
# Import PYPIT routines
from pypit import arparse as settings
from pypit import artrace
from pypit import arload
from pypit import arcomb
from pypit import arflux
from pypit import arlris
from pypit import armasters
from pypit import armsgs
from pypit import arproc
from pypit import arsort
from pypit import arutils
from pypit import arsave

from pypit import ardebug as debugger

# Logging
msgs = armsgs.get_logger()


class ScienceExposure:
    """
    A Science Exposure class that carries all information for a given science exposure
    """

    def __init__(self, snum, fitsdict, do_qa=True):

        # Set indices used for frame combination
        self._idx_sci = settings.spect['science']['index'][snum]
        self._idx_arcs = settings.spect['arc']['index'][snum]
        self._idx_std = settings.spect['standard']['index'][snum]
        if settings.argflag['bias']['useframe'] == 'bias': self._idx_bias = settings.spect['bias']['index'][snum]
        elif settings.argflag['bias']['useframe'] == 'dark':  self._idx_bias = settings.spect['dark']['index'][snum]
        else: self._idx_bias = []
        if settings.argflag['reduce']['trace']['useframe'] == 'trace': self._idx_trace = settings.spect['trace']['index'][snum]
        else: self._idx_trace = []
        if settings.argflag['reduce']['flatfield']['useframe'] == 'pixelflat': self._idx_flat = settings.spect['pixelflat']['index'][snum]
        elif settings.argflag['reduce']['flatfield']['useframe'] == 'trace': self._idx_flat = settings.spect['trace']['index'][snum]
        else: self._idx_flat = []
        if settings.argflag['reduce']['slitcen']['useframe'] == 'trace': self._idx_cent = settings.spect['trace']['index'][snum]
        elif settings.argflag['reduce']['slitcen']['useframe'] == 'pinhole': self._idx_cent = settings.spect['pinhole']['index'][snum]
        else: self._idx_cent = []
        self.sc = snum

        # Set the base name and extract other names that will be used for output files
        #  Also parses the time input
        self.SetBaseName(fitsdict)

        # Initialize the QA for this science exposure
        qafn = "{0:s}/QA_{1:s}.pdf".format(settings.argflag['run']['directory']['qa'], self._basename)
        self.qaroot = "{0:s}/PNGs/QA_{1:s}".format(settings.argflag['run']['directory']['qa'], self._basename)
        #if do_qa and not msgs._debug['no_qa']:
        #    self._qa = PdfPages(qafn)

        # Initialize Variables
        ndet = settings.spect['mosaic']['ndet']
        self._nonlinear = [settings.spect[settings.get_dnum(det+1)]['saturation'] *
                           settings.spect[settings.get_dnum(det+1)]['nonlinear']
                           for det in range(ndet)]
        self._nspec    = [None for all in range(ndet)]   # Number of spectral pixels
        self._nspat    = [None for all in range(ndet)]   # Number of spatial pixels
        self._datasec  = [None for all in range(ndet)]   # Locations of the data on each detector
        self._pixlocn  = [None for all in range(ndet)]   # Physical locations of each pixel on the detector
        self._lordloc  = [None for all in range(ndet)]   # Array of slit traces (left side) in physical pixel coordinates
        self._rordloc  = [None for all in range(ndet)]   # Array of slit traces (left side) in physical pixel coordinates
        self._pixcen   = [None for all in range(ndet)]   # Central slit traces in apparent pixel coordinates
        self._pixwid   = [None for all in range(ndet)]   # Width of slit (at each row) in apparent pixel coordinates
        self._lordpix  = [None for all in range(ndet)]   # Array of slit traces (left side) in apparent pixel coordinates
        self._rordpix  = [None for all in range(ndet)]   # Array of slit traces (right side) in apparent pixel coordinates
        self._slitpix  = [None for all in range(ndet)]   # Array identifying if a given pixel belongs to a given slit
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
        self._mspinhole = [None for all in range(ndet)]       # Master Pinhole
        self._mspixelflat = [None for all in range(ndet)]     # Master Pixel Flat
        self._mspixelflatnrm = [None for all in range(ndet)]  # Normalized Master pixel flat
        self._msblaze = [None for all in range(ndet)]       # Blaze function
        self._msstd = [{} for all in range(ndet)]           # Master Standard dict
        self._sensfunc = None                               # Sensitivity function
        # Initialize the Master Calibration frame names
        self._msarc_name = [None for all in range(ndet)]      # Master Arc Name
        self._msbias_name = [None for all in range(ndet)]     # Master Bias Name
        self._mstrace_name = [None for all in range(ndet)]    # Master Trace Name
        self._mspinhole_name = [None for all in range(ndet)]    # Master Pinhole Name
        self._mspixelflat_name = [None for all in range(ndet)]  # Master Pixel Flat Name
        # Initialize the science, variance, and background frames
        self._sciframe = [None for all in range(ndet)]
        self._rawvarframe = [None for all in range(ndet)]    # Variance based on detected counts + RN
        self._modelvarframe = [None for all in range(ndet)]  # Variance from sky and object models
        self._bgframe = [None for all in range(ndet)]
        self._scimask = [None for all in range(ndet)]        # Mask (1=Bad pix; 2=CR)
        self._scitrace = [None for all in range(ndet)]
        self._slitprof = [None for all in range(ndet)]   # Slit profiles at each position on the detector
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
        import datetime
        #
        scidx = self._idx_sci[0]
        tbname = None
        try:
            if "T" in fitsdict['date'][scidx]:
                tbname = fitsdict['date'][scidx]
        except IndexError:
            debugger.set_trace()
        else:
            if tbname is None:
                if settings.spect["fits"]["timeunit"] == "mjd":
                    # Not ideal, but convert MJD into a date+time
                    timval = Time(fitsdict['time'][scidx] / 24.0, scale='tt', format='mjd')
                    tbname = timval.isot
                else:
                    # Really not ideal... just append date and time
                    tbname = fitsdict['date'][scidx] + "T" + str(fitsdict['time'][scidx])
        '''
        if "T" in fitsdict['date'][scidx]:
            tbname = fitsdict['date'][scidx]
        else:
            # Not ideal, but convert MJD into a date+time
            debugger.set_trace() # CANNOT GET HERE
            timval = Time(fitsdict['time'][scidx]/24.0, scale='tt', format='mjd')
            tbname = timval.isot
        '''
        tval = Time(tbname, format='isot')#'%Y-%m-%dT%H:%M:%S.%f')
        dtime = datetime.datetime.strptime(tval.value, '%Y-%m-%dT%H:%M:%S.%f')
        #except ValueError:
            #tval = datetime.datetime.strptime(tbname, '%Y-%m-%dT%H:%M:%S')
        self._inst_name = settings.spect['mosaic']['camera']
        self._target_name = fitsdict['target'][self._idx_sci[0]].replace(" ", "")
        self._basename = self._target_name+'_'+self._inst_name+'_'+ \
                         datetime.datetime.strftime(dtime, '%Y%b%dT') + \
                         tbname.split("T")[1].replace(':','')
        # Save Time object
        self._time = tval
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
        if settings.argflag['reduce']['badpix'] == 'bias':
            if settings.argflag['reduce']['masters']['reuse']:
                # Attempt to load the Master Frame
                bpix_name = armasters.master_name('badpix', settings.argflag['reduce']['masters']['setup'])
                try:
                    bpix, head = arload.load_master(bpix_name, frametype="badpix")
                except IOError:
                    msgs.warn("No MasterBadPix frame found {:s}".format(bpix_name))
                else:
                    settings.argflag['reduce']['masters']['loaded'].append('badpix'+settings.argflag['reduce']['masters']['setup'])
            if 'badpix'+settings.argflag['reduce']['masters']['setup'] not in settings.argflag['reduce']['masters']['loaded']:
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
            if settings.argflag['run']['spectrograph'] in ['lris_red']:
                bpix = arlris.bpm(self, 'red', fitsdict, det)
            else:
                msgs.info("Not preparing a bad pixel mask")
                return False
        self.SetFrame(self._bpix, bpix, det)
        del bpix
        return True

    def GetPixelLocations(self, det):
        """
        Generate or load the physical location of each pixel

        Parameters
        ----------
        det : int
          Index of the detector
        """
        if settings.argflag['reduce']['pixel']['locations'] is None:
            self.SetFrame(self._pixlocn, artrace.gen_pixloc(self._mstrace[det-1], det, gen=True), det)
        elif settings.argflag['reduce']['pixel']['locations'] in ["mstrace"]:
            self.SetFrame(self._pixlocn, artrace.gen_pixloc(self._mstrace[det-1], det, gen=False), det)
        else:
            mname = settings.argflag['run']['directory']['master']+'/'+settings.argflag['reduce']['pixel']['locations']
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
        dnum = settings.get_dnum(det)

        if self._msarc[det-1] is not None:
            msgs.info("An identical master arc frame already exists")
            return False
        if settings.argflag['arc']['useframe'] in ['arc']:
            # Attempt to load the Master Frame
            if settings.argflag['reduce']['masters']['reuse']:
                msarc_name = armasters.master_name('arc', settings.argflag['reduce']['masters']['setup'])
                try:
                    msarc, head = arload.load_master(msarc_name, frametype="arc")
                except IOError:
                    msgs.warn("No MasterArc frame found {:s}".format(msarc_name))
                else:
                    self._transpose = head['transp']
                    if self._transpose:  # Need to setup for flipping
                        settings.argflag['trace']['dispersion']['direction'] = 1
                    else:
                        settings.argflag['trace']['dispersion']['direction'] = 0
                    # Append as loaded
                    settings.argflag['reduce']['masters']['loaded'].append('arc'+settings.argflag['reduce']['masters']['setup'])
            if 'arc'+settings.argflag['reduce']['masters']['setup'] not in settings.argflag['reduce']['masters']['loaded']:
                msgs.info("Preparing a master arc frame")
                ind = self._idx_arcs
                # Load the arc frames
                frames = arload.load_frames(fitsdict, ind, det, frametype='arc', msbias=self._msbias[det-1])
                if settings.argflag['arc']['combine']['match'] > 0.0:
                    sframes = arsort.match_frames(frames, settings.argflag['arc']['combine']['match'], frametype='arc',
                                                  satlevel=settings.spect[dnum]['saturation']*settings.spect[dnum]['nonlinear'])
                    subframes = np.zeros((frames.shape[0], frames.shape[1], len(sframes)))
                    numarr = np.array([])
                    for i in range(len(sframes)):
                        numarr = np.append(numarr, sframes[i].shape[2])
                        msarc = arcomb.comb_frames(sframes[i], det, 'arc')
                        # Send the data away to be saved
                        subframes[:,:,i] = msarc.copy()
                    del sframes
                    # Combine all sub-frames
                    msarc = arcomb.comb_frames(subframes, det, 'arc', weights=numarr)
                    del subframes
                else:
                    msarc = arcomb.comb_frames(frames, det, 'arc')
                del frames
            # # Derive a suitable name for the master arc frame
            # msarc_name = "{0:s}/{1:s}/msarc{2:s}_{3:03d}.fits".format(os.getcwd(),settings.argflag['run']['directory']['master'],settings.spect["det"][det-1]["suffix"],len(self._done_arcs))
            # self._tltprefix = os.path.splitext(os.path.basename(msarc_name))[0]
            # # Send the data away to be saved
            # arsave.save_master(self, msarc, filename=msarc_name, frametype='arc', ind=ind)
            # # Store the files used and the master bias name in case it can be used during the later reduction processes
            # self._done_arcs.append(ind)
            # self._name_arcs.append(msarc_name)
        else:
            msarc_name = settings.argflag['run']['directory']['master']+'/'+settings.argflag['arc']['useframe']
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
            msgs.info("An identical master {0:s} frame already exists".format(settings.argflag['bias']['useframe']))
            return False
        elif settings.argflag['bias']['useframe'] in ['bias', 'dark']:
            # Load from hard-drive?
            if settings.argflag['reduce']['masters']['reuse']:
                # Attempt to load the Master Frame
                msbias_name = armasters.master_name('bias', settings.argflag['reduce']['masters']['setup'])
                try:
                    msbias, head = arload.load_master(msbias_name, frametype="bias")
                except IOError:
                    msgs.warn("No MasterBias frame found {:s}".format(msbias_name))
                else:
                    settings.argflag['reduce']['masters']['loaded'].append('bias'+settings.argflag['reduce']['masters']['setup'])
            if 'bias'+settings.argflag['reduce']['masters']['setup'] not in settings.argflag['reduce']['masters']['loaded']:
                msgs.info("Preparing a master {0:s} frame".format(settings.argflag['bias']['useframe']))
                # Get all of the bias frames for this science frame
                ind = self._idx_bias
                # Load the Bias/Dark frames
                frames = arload.load_frames(fitsdict, ind, det, frametype=settings.argflag['bias']['useframe'])
                msbias = arcomb.comb_frames(frames, det, 'bias', printtype=settings.argflag['bias']['useframe'])
                del frames
        elif settings.argflag['bias']['useframe'] == 'overscan':
            self.SetMasterFrame('overscan', "bias", det, mkcopy=False)
            return False
        elif settings.argflag['bias']['useframe'] == 'none':
            msgs.info("Not performing a bias/dark subtraction")
            self.SetMasterFrame(None, "bias", det, mkcopy=False)
            return False
        else:  # It must be the name of a file the user wishes to load
            msbias_name = settings.argflag['run']['directory']['master']+u'/'+settings.argflag['bias']['useframe']
            msbias, head = arload.load_master(msbias_name, frametype="bias")
            settings.argflag['reduce']['masters']['loaded'].append('bias')
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
        if settings.spect['det'][det-1]['numamplifiers'] > 1:
            msgs.work("Readnoise needs to be updated for multiple amps")
        # Set
        rnoise = settings.spect['det'][det-1]['ronoise'] #+ (0.5*settings.spect['det'][det-1]['gain'])**2
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
        from pypit import arqa

        if settings.argflag['reduce']['flatfield']['perform']:  # Only do it if the user wants to flat field
            # If the master pixelflat is already made, use it
            if self._mspixelflat[det-1] is not None:
                msgs.info("An identical master pixelflat frame already exists")
                if self._mspixelflatnrm[det-1] is None:
                    # Normalize the flat field
                    msgs.info("Normalizing the pixel flat")
                    slit_profiles, mstracenrm, msblaze, flat_ext1d = \
                        arproc.slit_profile(self, self.GetMasterFrame("pixelflat", det),
                                            det, ntcky=settings.argflag['reduce']['flatfield']['params'][0])
                    # mspixelflatnrm, msblaze = arproc.flatnorm(self, det, self.GetMasterFrame("pixelflat", det),
                    #                                           overpix=0, plotdesc="Blaze function")
                    mspixelflatnrm = mstracenrm.copy()
                    winpp = np.where(slit_profiles != 0.0)
                    mspixelflatnrm[winpp] /= slit_profiles[winpp]
                    self.SetMasterFrame(mspixelflatnrm, "normpixelflat", det)
                    if np.array_equal(self._idx_flat, self._idx_trace):
                        # The flat field frame is also being used to trace the slit edges and determine the slit
                        # profile. Avoid recalculating the slit profile and blaze function and save them here.
                        self.SetFrame(self._msblaze, msblaze, det)
                        self.SetFrame(self._slitprof, slit_profiles, det)
                        if settings.argflag["reduce"]["slitprofile"]["perform"]:
                            msgs.info("Preparing QA of each slit profile")
                            arqa.slit_profile(self, mstracenrm, slit_profiles, self._lordloc[det - 1], self._rordloc[det - 1],
                                              self._slitpix[det - 1], desc="Slit profile")
                        msgs.info("Saving blaze function QA")
                        arqa.plot_orderfits(self, msblaze, flat_ext1d, desc="Blaze function")
                return False
            ###############
            # Generate a master pixel flat frame
            if settings.argflag['reduce']['flatfield']['useframe'] in ['pixelflat', 'trace']:
                if settings.argflag['reduce']['masters']['reuse']:
                    # Attempt to load the Master Frame
                    msflat_name = armasters.master_name('normpixelflat', settings.argflag['reduce']['masters']['setup'])
                    try:
                        mspixelflatnrm, head = arload.load_master(msflat_name, frametype="normpixelflat")
                    except IOError:
                        msgs.warn("No MasterFlatField frame found {:s}".format(msflat_name))
                    else:
                        settings.argflag['reduce']['masters']['loaded'].append('normpixelflat'+settings.argflag['reduce']['masters']['setup'])
                        mspixelflat = mspixelflatnrm
                if 'normpixelflat'+settings.argflag['reduce']['masters']['setup'] not in settings.argflag['reduce']['masters']['loaded']:
                    msgs.info("Preparing a master pixel flat frame with {0:s}".format(settings.argflag['reduce']['flatfield']['useframe']))
                    # Get all of the pixel flat frames for this science frame
                    ind = self._idx_flat
                    # Load the frames for tracing
                    frames = arload.load_frames(fitsdict, ind, det, frametype='pixel flat',
                                                msbias=self._msbias[det-1])
                    if settings.argflag['pixelflat']['combine']['match'] > 0.0:
                        sframes = arsort.match_frames(frames, settings.argflag['pixelflat']['combine']['match'],
                                                      frametype='pixel flat', satlevel=self._nonlinear)
                        subframes = np.zeros((frames.shape[0], frames.shape[1], len(sframes)))
                        numarr = np.array([])
                        for i in range(len(sframes)):
                            numarr = np.append(numarr, sframes[i].shape[2])
                            mspixelflat = arcomb.comb_frames(sframes[i], det, 'pixelflat', printtype='pixel flat')
                            subframes[:, :, i] = mspixelflat.copy()
                        del sframes
                        # Combine all sub-frames
                        mspixelflat = arcomb.comb_frames(subframes, det, 'pixelflat', weights=numarr,
                                                         printtype='pixel flat')
                        del subframes
                    else:
                        mspixelflat = arcomb.comb_frames(frames, det, 'pixelflat', printtype='pixel flat')
                    del frames
                    # Apply gain (instead of ampsec scale)
                    mspixelflat *= arproc.gain_frame(self, det)
                    # Normalize the flat field
                    msgs.info("Normalizing the pixel flat")
                    slit_profiles, mstracenrm, msblaze, flat_ext1d = \
                        arproc.slit_profile(self, mspixelflat, det, ntcky=settings.argflag['reduce']['flatfield']['params'][0])
                    # mspixelflatnrm, msblaze = arproc.flatnorm(self, det, self.GetMasterFrame("pixelflat", det),
                    #                                         overpix=0, plotdesc="Blaze function")
                    mspixelflatnrm = mstracenrm.copy()
                    winpp = np.where(slit_profiles != 0.0)
                    mspixelflatnrm[winpp] /= slit_profiles[winpp]
                    if np.array_equal(self._idx_flat, self._idx_trace):
                        # The flat field frame is also being used to trace the slit edges and determine the slit
                        # profile. Avoid recalculating the slit profile and blaze function and save them here.
                        self.SetFrame(self._msblaze, msblaze, det)
                        self.SetFrame(self._slitprof, slit_profiles, det)
                        if settings.argflag["reduce"]["slitprofile"]["perform"]:
                            msgs.info("Preparing QA of each slit profile")
                            arqa.slit_profile(self, mstracenrm, slit_profiles, self._lordloc[det - 1], self._rordloc[det - 1],
                                              self._slitpix[det - 1], desc="Slit profile")
                        msgs.info("Saving blaze function QA")
                        arqa.plot_orderfits(self, msblaze, flat_ext1d, desc="Blaze function")
            else:  # It must be the name of a file the user wishes to load
                mspixelflat_name = armasters.user_master_name(settings.argflag['run']['directory']['master'],
                                                              settings.argflag['reduce']['flatfield']['useframe'])
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

    def MasterPinhole(self, fitsdict, det):
        """
        Generate Master pinhole frame for a given detector

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
        dnum = settings.get_dnum(det)
        # If the master pinhole is already made, use it
        if self._mspinhole[det - 1] is not None:
            msgs.info("An identical master pinhole frame already exists")
            return False
        if settings.argflag['reduce']['slitcen']['useframe'] in ['trace', 'pinhole']:
            if settings.argflag['reduce']['masters']['reuse']:
                # Attempt to load the Master Frame
                mspinhole_name = armasters.master_name('pinhole', settings.argflag['reduce']['masters']['setup'])
                try:
                    mspinhole, head = arload.load_master(mspinhole_name, frametype="pinhole")
                except IOError:
                    msgs.warn("No MasterPinhole frame found {:s}".format(mspinhole_name))
                else:
                    settings.argflag['reduce']['masters']['loaded'].append(
                        'pinhole' + settings.argflag['reduce']['masters']['setup'])
            if 'pinhole' + settings.argflag['reduce']['masters']['setup'] not in \
                    settings.argflag['reduce']['masters']['loaded']:
                msgs.info("Preparing a master pinhole frame with {0:s}".format(
                    settings.argflag['reduce']['slitcen']['useframe']))
                ind = self._idx_cent
                # Load the pinhole frames
                frames = arload.load_frames(fitsdict, ind, det, frametype='pinhole', msbias=self._msbias[det - 1],
                                            trim=settings.argflag['reduce']['trim'])
                if settings.argflag['pinhole']['combine']['match'] > 0.0:
                    sframes = arsort.match_frames(frames, settings.argflag['pinhole']['combine']['match'],
                                                  frametype='pinhole', satlevel=settings.spect[dnum]['saturation'] *
                                                  settings.spect['det'][det - 1]['nonlinear'])
                    subframes = np.zeros((frames.shape[0], frames.shape[1], len(sframes)))
                    numarr = np.array([])
                    for i in range(len(sframes)):
                        numarr = np.append(numarr, sframes[i].shape[2])
                        mspinhole = arcomb.comb_frames(sframes[i], det, 'pinhole')
                        subframes[:, :, i] = mspinhole.copy()
                    del sframes
                    # Combine all sub-frames
                    mspinhole = arcomb.comb_frames(subframes, det, 'pinhole', weights=numarr)
                    del subframes
                else:
                    mspinhole = arcomb.comb_frames(frames, det, 'pinhole')
                del frames
        else:  # It must be the name of a file the user wishes to load
            mspinhole_name = settings.argflag['run']['directory']['master'] + '/' + \
                              settings.argflag['reduce']['slitcen']['useframe']
            mspinhole, head = arload.load_master(mspinhole_name, frametype=None)
            debugger.set_trace()  # NEED TO LOAD EXTRAS AS ABOVE
        # Set and then delete the Master Trace frame
        self.SetMasterFrame(mspinhole, "pinhole", det)
        del mspinhole
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
        dnum = settings.get_dnum(det)
        # If the master trace is already made, use it
        if self._mstrace[det-1] is not None:
            msgs.info("An identical master trace frame already exists")
            return False
        if settings.argflag['reduce']['trace']['useframe'] in ['trace']:
            if settings.argflag['reduce']['masters']['reuse']:
                # Attempt to load the Master Frame
                mstrace_name = armasters.master_name('trace', settings.argflag['reduce']['masters']['setup'])
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
                    slitpix, _ = arload.load_master(mstrace_name, frametype="trace", exten=7)
                    self.SetFrame(self._lordloc, lordloc, det)
                    self.SetFrame(self._rordloc, rordloc, det)
                    self.SetFrame(self._pixcen, pixcen.astype(np.int), det)
                    self.SetFrame(self._pixwid, pixwid.astype(np.int), det)
                    self.SetFrame(self._lordpix, lordpix.astype(np.int), det)
                    self.SetFrame(self._rordpix, rordpix.astype(np.int), det)
                    self.SetFrame(self._slitpix, slitpix.astype(np.int), det)
                    #
                    settings.argflag['reduce']['masters']['loaded'].append('trace'+settings.argflag['reduce']['masters']['setup'])
            if 'trace'+settings.argflag['reduce']['masters']['setup'] not in settings.argflag['reduce']['masters']['loaded']:
                msgs.info("Preparing a master trace frame with {0:s}".format(settings.argflag['reduce']['trace']['useframe']))
                ind = self._idx_trace
                # Load the frames for tracing
                frames = arload.load_frames(fitsdict, ind, det, frametype='trace', msbias=self._msbias[det-1],
                                            trim=settings.argflag['reduce']['trim'])
                if settings.argflag['trace']['combine']['match'] > 0.0:
                    sframes = arsort.match_frames(frames, settings.argflag['trace']['combine']['match'], frametype='trace', satlevel=settings.spect[dnum]['saturation']*settings.spect['det'][det-1]['nonlinear'])
                    subframes = np.zeros((frames.shape[0], frames.shape[1], len(sframes)))
                    numarr = np.array([])
                    for i in range(len(sframes)):
                        numarr = np.append(numarr, sframes[i].shape[2])
                        mstrace = arcomb.comb_frames(sframes[i], det, 'trace')
                        subframes[:,:,i] = mstrace.copy()
                    del sframes
                    # Combine all sub-frames
                    mstrace = arcomb.comb_frames(subframes, det, 'trace', weights=numarr)
                    del subframes
                else:
                    mstrace = arcomb.comb_frames(frames, det, 'trace')
                del frames
        else: # It must be the name of a file the user wishes to load
            mstrace_name = settings.argflag['run']['directory']['master']+'/'+settings.argflag['reduce']['trace']['useframe']
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
        if settings.argflag['reduce']['masters']['reuse']:
            mswave_name = armasters.master_name('wave', settings.argflag['reduce']['masters']['setup'])
            try:
                mswave, head = arload.load_master(mswave_name, frametype="arc")
            except IOError:
                msgs.warn("No MasterWave frame found {:s}".format(mswave_name))
            else:
                settings.argflag['reduce']['masters']['loaded'].append('wave'+settings.argflag['reduce']['masters']['setup'])
        if 'wave'+settings.argflag['reduce']['masters']['setup'] not in settings.argflag['reduce']['masters']['loaded']:
            msgs.info("Preparing a master wave frame")
            if settings.argflag["reduce"]["calibrate"]["wavelength"] == "pixel":
                mswave = self._tilts[det - 1] * (self._tilts[det - 1].shape[0]-1.0)
            else:
                wv_calib = self._wvcalib[det - 1]
                mswave = arutils.func_val(wv_calib['fitc'], self._tilts[det - 1], wv_calib['function'],
                                          minv=wv_calib['fmin'], maxv=wv_calib['fmax'])
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
        else:
            wv_calib = None
        # Attempt to load the Master Frame
        if settings.argflag['reduce']['masters']['reuse']:
            mswv_calib_name = armasters.master_name('wave_calib', settings.argflag['reduce']['masters']['setup'])
            try:
                wv_calib = arload.load_master(mswv_calib_name, frametype="wv_calib")
            except (IOError, ValueError):
                msgs.warn("No MasterWave1D data found {:s}".format(mswv_calib_name))
            else:
                settings.argflag['reduce']['masters']['loaded'].append('wave_calib'+settings.argflag['reduce']['masters']['setup'])
        if settings.argflag["reduce"]["calibrate"]["wavelength"] == "pixel":
            msgs.info("A wavelength calibration will not be performed")
        else:
            if 'wave_calib' + settings.argflag['reduce']['masters']['setup'] not in settings.argflag['reduce']['masters']['loaded']:
                # Setup arc parameters (e.g. linelist)
                arcparam = ararc.setup_param(self, sc, det, fitsdict)
                self.SetFrame(self._arcparam, arcparam, det)
                ###############
                # Extract arc and identify lines
                if settings.argflag['arc']['calibrate']['method'] == 'simple':
                    wv_calib = ararc.simple_calib(self, det)
                elif settings.argflag['arc']['calibrate']['method'] == 'arclines':
                    wv_calib = ararc.calib_with_arclines(self, det)
        # Set
        if wv_calib is not None:
            self.SetFrame(self._wvcalib, wv_calib, det)
        del wv_calib
        return True

    def MasterStandard(self, fitsdict):
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
        if self._sensfunc is not None:
            msgs.info("Using existing sensitivity function.")
            return False
        # Attempt to load the Master Frame
        if settings.argflag['reduce']['masters']['reuse']:
            sfunc_name = armasters.master_name('sensfunc',
                settings.argflag['reduce']['masters']['setup'])
            try:
                sensfunc = arload.load_master(sfunc_name, frametype="sensfunc")
            except (IOError, ValueError):
                msgs.warn("No MasterSensFunc data found {:s}".format(sfunc_name))
            else:
                msgs.info("Loaded sensitivity function from {:s}".format(sfunc_name))
                settings.argflag['reduce']['masters']['loaded'].append('sensfunc'+settings.argflag['reduce']['masters']['setup'][0])
                self._sensfunc = sensfunc.copy()
                return True
        # Grab the standard star frames
        msgs.info("Preparing the standard")
        ind = self._idx_std
        msgs.warn("Taking only the first standard frame for now")
        ind = [ind[0]]
        # Extract
        all_specobj = []
        for kk in range(settings.spect['mosaic']['ndet']):
            det = kk+1
            # Load the frame(s)
            frame = arload.load_frames(fitsdict, ind, det, frametype='standard',
                                       msbias=self._msbias[det-1])
            sciframe = frame[:, :, 0] # First exposure
            # Save RA/DEC
            self._msstd[det-1]['RA'] = fitsdict['ra'][ind[0]]
            self._msstd[det-1]['DEC'] = fitsdict['dec'][ind[0]]
            self._msstd[det-1]['spobjs'] = None
            # Use this detector? Need to check this after setting RA/DEC above
            if settings.argflag['reduce']['detnum'] is not None:
                msgs.warn("If your standard wasnt on this detector, you will have trouble..")
                if det != settings.argflag['reduce']['detnum']:
                    continue
            if settings.spect["mosaic"]["reduction"] == "ARMLSD":
                arproc.reduce_multislit(self, sciframe, ind[0], fitsdict, det, standard=True)
            elif settings.spect["mosaic"]["reduction"] == "ARMED":
                arproc.reduce_echelle(self, sciframe, ind[0], fitsdict, det, standard=True)
            else:
                msgs.error("Not ready for reduction type {0:s}".format(settings.spect["mosaic"]["reduction"]))

            if self._msstd[det-1]['spobjs'] is not None:
                all_specobj += self._msstd[det-1]['spobjs']
        # If standard, generate a sensitivity function
        sensfunc = arflux.generate_sensfunc(self, ind[0], all_specobj, fitsdict)
        # Set the sensitivity function
        self.SetMasterFrame(sensfunc, "sensfunc", None, mkcopy=False)
        # Apply to Standard
        for kk in range(settings.spect['mosaic']['ndet']):
            det = kk + 1  # Detectors indexed from 1
            arflux.apply_sensfunc(self, det, ind[0], fitsdict, standard=True)
        # Save
        armasters.save_sensfunc(self, settings.argflag['reduce']['masters']['setup'])
        # Save standard star spectrum to disk
        outfile = settings.argflag['run']['directory']['science']+'/spec1d_{:s}.fits'.format(
            fitsdict['filename'][ind[0]].split('.')[0])
        arsave.save_1d_spectra_fits(self, fitsdict, standard=True, outfile=outfile)
        return True

    '''
    def Setup(self):

        # Sort the data
        msgs.bug("Files and folders should not be deleted -- there should be an option to overwrite files automatically if they already exist, or choose to rename them if necessary")
        self._filesort = arsort.sort_data(self)
        # Write out the details of the sorted files
        if settings.argflag['output']['sorted'] is not None: arsort.sort_write(self)
        # Match Science frames to calibration frames
        arsort.match_science(self)
        # If the user is only debugging, then exit now
        if settings.argflag['run']['calcheck']:
            msgs.info("Calibration check complete. Change the 'calcheck' flag to continue with data reduction")
            sys.exit()
        # Make directory structure for different objects
        self._sci_targs = arsort.make_dirs(self)
        return
    '''

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
        elif ftype == "pinhole": self._mspinhole[det] = cpf
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
            elif ftype == "pinhole": return self._mspinhole[det].copy()
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
            elif ftype == "pinhole": return self._mspinhole[det]
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
