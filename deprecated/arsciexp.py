""" Class for book-keeping the reduction process
"""
from __future__ import (absolute_import, division, print_function, unicode_literals)

import sys
import datetime

import numpy as np

from astropy.time import Time

# Import PYPIT routines
from pypit import msgs
from pypit import arparse as settings
from pypit import arload
from pypit import arcomb
from pypit import armasters
from pypit.core import arsort
from pypit import arutils
from pypit import ardebug as debugger

class ScienceExposure:
    """
    A Science Exposure class that carries all information for a given science exposure
    """

    def __init__(self, sci_ID, fitstbl, settings_argflag, settings_spect, do_qa=True,
                 idx_sci=None):

        # Set indices used for frame combination
        msgs.error("DEPRECATED")
        self.sci_ID = sci_ID  # Binary 1,2,4,8,..
        self._idx_sci = np.where((fitstbl['sci_ID'] == sci_ID) & fitstbl['science'])[0]
        if idx_sci is not None:
            self._idx_sci = np.array([idx_sci])
        if settings_argflag['reduce']['masters']['force']:
            #self._idx_bias = []
            #self._idx_flat = []
            self._idx_cent = []
            #self._idx_trace = []
            #self._idx_arcs = []
            #self._idx_std = []
        else:
            #self._idx_arcs = arsort.ftype_indices(fitstbl, 'arc', self.sci_ID)
            #self._idx_std = arsort.ftype_indices(fitstbl, 'standard', self.sci_ID)
            # Bias
            #if settings_argflag['bias']['useframe'] == 'bias':
            #    self._idx_bias = arsort.ftype_indices(fitstbl, 'bias', self.sci_ID)
            #elif settings_argflag['bias']['useframe'] == 'dark':
            #    self._idx_bias = arsort.ftype_indices(fitstbl, 'dark', self.sci_ID)
            #else: self._idx_bias = []
            # Trace
            #self._idx_trace = arsort.ftype_indices(fitstbl, 'trace', self.sci_ID)
            # Flat
            #if settings_argflag['reduce']['flatfield']['useframe'] == 'pixelflat':
            #    self._idx_flat = arsort.ftype_indices(fitstbl, 'pixelflat', self.sci_ID)
            #elif settings_argflag['reduce']['flatfield']['useframe'] == 'trace':
            #    self._idx_flat = arsort.ftype_indices(fitstbl, 'trace', self.sci_ID)
            #else: self._idx_flat = []
            # Cent
            if settings_argflag['reduce']['slitcen']['useframe'] == 'trace':
                self._idx_cent = arsort.ftype_indices(fitstbl, 'trace', self.sci_ID)
            elif settings_argflag['reduce']['slitcen']['useframe'] == 'pinhole':  # Not sure this will work
                self._idx_cent = arsort.ftype_indices(fitstbl, 'pinhole', self.sci_ID)
            else: self._idx_cent = []

        # Set the base name and extract other names that will be used for output files
        #  Also parses the time input
        self.SetBaseName(fitstbl)

        # Setup
        self.setup = ''

        # Velocity correction (e.g. heliocentric)
        self.vel_correction = 0.

        # Initialize the QA for this science exposure
        qafn = "{0:s}/QA_{1:s}.pdf".format(settings_argflag['run']['directory']['qa'], self._basename)
        self.qaroot = "{0:s}/PNGs/QA_{1:s}".format(settings_argflag['run']['directory']['qa'], self._basename)

        # Initialize Variables
        ndet = settings_spect['mosaic']['ndet']
        self.extracted = [False for all in range(ndet)]   # Mainly for standard stars
        self._nonlinear = [settings_spect[settings.get_dnum(det+1)]['saturation'] *
                           settings_spect[settings.get_dnum(det+1)]['nonlinear']
                           for det in range(ndet)]
        #self._nspec    = [None for all in range(ndet)]   # Number of spectral pixels
        #self._nspat    = [None for all in range(ndet)]   # Number of spatial pixels
        #self._datasec  = [None for all in range(ndet)]   # Locations of the data on each detector
        self._pixlocn  = [None for all in range(ndet)]   # Physical locations of each pixel on the detector
        self._lordloc  = [None for all in range(ndet)]   # Array of slit traces (left side) in physical pixel coordinates
        self._rordloc  = [None for all in range(ndet)]   # Array of slit traces (left side) in physical pixel coordinates
        self._pixcen   = [None for all in range(ndet)]   # Central slit traces in apparent pixel coordinates
        self._pixwid   = [None for all in range(ndet)]   # Width of slit (at each row) in apparent pixel coordinates
        self._lordpix  = [None for all in range(ndet)]   # Array of slit traces (left side) in apparent pixel coordinates
        self._rordpix  = [None for all in range(ndet)]   # Array of slit traces (right side) in apparent pixel coordinates
        self._slitpix  = [None for all in range(ndet)]   # Array identifying if a given pixel belongs to a given slit
        #self._tilts    = [None for all in range(ndet)]   # Array of spectral tilts at each position on the detector
        #self._tiltpar  = [None for all in range(ndet)]   # Dict parameters for tilt fitting
        self._satmask  = [None for all in range(ndet)]   # Array of Arc saturation streaks
        #self._arcparam = [None for all in range(ndet)]   # Dict guiding wavelength calibration
        #self._wvcalib  = [None for all in range(ndet)]   # List of dict's
        self._resnarr  = [None for all in range(ndet)]   # Resolution array
        self._maskslits = [None for all in range(ndet)]  # Mask for whether to analyze a given slit (True=masked)
        # Initialize the Master Calibration frames
        #self._bpix = [None for all in range(ndet)]          # Bad Pixel Mask
        #self._msarc = [None for all in range(ndet)]         # Master Arc
        self._mswave = [None for all in range(ndet)]         # Master Wavelength image
        #self._msbias = [None for all in range(ndet)]        # Master Bias
        self._msrn = [None for all in range(ndet)]          # Master ReadNoise image
        #self._mstrace = [None for all in range(ndet)]       # Master Trace
        self._mspinhole = [None for all in range(ndet)]       # Master Pinhole
        #self._mspixelflat = [None for all in range(ndet)]     # Master Pixel Flat
        #self._mspixelflatnrm = [None for all in range(ndet)]  # Normalized Master pixel flat
        #self._msblaze = [None for all in range(ndet)]       # Blaze function
        #self._msstd = [{} for all in range(ndet)]           # Master Standard dict
        #self._sensfunc = None                               # Sensitivity function
        # Initialize the Master Calibration frame names
        #self._msarc_name = [None for all in range(ndet)]      # Master Arc Name
        #self._msbias_name = [None for all in range(ndet)]     # Master Bias Name
        #self._mstrace_name = [None for all in range(ndet)]    # Master Trace Name
        self._mspinhole_name = [None for all in range(ndet)]    # Master Pinhole Name
        #self._mspixelflat_name = [None for all in range(ndet)]  # Master Pixel Flat Name
        # Initialize the science, variance, and background frames
        self._sciframe = [None for all in range(ndet)]
        self._rawvarframe = [None for all in range(ndet)]    # Variance based on detected counts + RN
        self._modelvarframe = [None for all in range(ndet)]  # Variance from sky and object models
        self._bgframe = [None for all in range(ndet)]
        self._scimask = [None for all in range(ndet)]        # Mask (1=Bad pix; 2=CR)
        self._scitrace = [None for all in range(ndet)]
        #self._slitprof = [None for all in range(ndet)]   # Slit profiles at each position on the detector
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
    '''
    def BadPixelMask(self, fitsdict, det, msbias):
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
            try:
                bpix = armasters.get_master_frame(self, "badpix")
            except IOError:
                msgs.info("Preparing a bad pixel mask")
                # Get all of the bias frames for this science frame
                if len(self._idx_bias) == 0:
                    msgs.warn("No bias frames available to determine bad pixel mask")
                    msgs.info("Not preparing a bad pixel mask")
                    return False
                # Load the Bias frames
                bpix = arproc.badpix(self, det, msbias) # self.GetMasterFrame('bias', det))
        else:
            # Instrument dependent
            if settings.argflag['run']['spectrograph'] in ['keck_lris_red']:
                bpix = arlris.bpm(self, 'red', fitsdict, det)
            else:
                msgs.info("Not preparing a bad pixel mask")
                return False
        # Save
        self.SetFrame(self._bpix, bpix, det)
        armasters.save_masters(self, det, mftype='badpix')
        del bpix
        return True
    '''

    '''
    def GetPixelLocations(self, det):
        """
        Generate or load the physical location of each pixel

        Parameters
        ----------
        det : int
          Index of the detector
        """
        if settings.argflag['reduce']['pixel']['locations'] is None:
            self.SetFrame(self._pixlocn, arpixels.gen_pixloc(self._mstrace[det-1], det, gen=True), det)
        elif settings.argflag['reduce']['pixel']['locations'] in ["mstrace"]:
            self.SetFrame(self._pixlocn, arpixels.gen_pixloc(self._mstrace[det-1], det, gen=False), det)
        else:
            mname = settings.argflag['run']['directory']['master']+'/'+settings.argflag['reduce']['pixel']['locations']
            self.SetFrame(self._pixlocn, armasters.load_master(mname, frametype=None), det)
        return
    '''

    '''
    def MasterArc(self, fitsdict, det, msbias):
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
            msgs.info("A master arc frame already exists for this frame")
            return False
        if settings.argflag['arc']['useframe'] in ['arc']:
            # Master Frame
            msarc = armasters.load_master_frame(self, "arc")
            if msarc is None:
                msgs.info("Preparing a master arc frame")
                ind = self._idx_arcs
                # Load the arc frames
                frames = arload.load_frames(fitsdict, ind, det, frametype='arc', msbias=msbias) #self._msbias[det-1])
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
        else: # Use input frame name located in MasterFrame directory
            msarc_name = settings.argflag['run']['directory']['master']+'/'+settings.argflag['arc']['useframe']
            msarc, _ = armasters.load_master(msarc_name, frametype=None)

        # Set and then delete the Master Arc frame
        self.SetMasterFrame(msarc, "arc", det)
        armasters.save_masters(self, det, mftype='arc')
        del msarc
        return True
    '''

    '''
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
        msgs.error("DEPRECATED)
        # If the master bias is already made, use it
        if self._msbias[det-1] is not None:
            msgs.info("An identical master {0:s} frame already exists".format(settings.argflag['bias']['useframe']))
            return False
        elif settings.argflag['bias']['useframe'] in ['bias', 'dark']:
            try:
                msbias = armasters.get_master_frame(self, "bias")
            except IOError:
                msgs.info("Preparing a master {0:s} frame".format(settings.argflag['bias']['useframe']))
                # Get all of the bias frames for this science frame
                ind = self._idx_bias
                # Load the Bias/Dark frames
                frames = arload.load_frames(fitsdict, ind, det,
                                            frametype=settings.argflag['bias']['useframe'], trim=False)
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
            msbias, head = armasters.load_master(msbias_name, frametype="bias")
            settings.argflag['reduce']['masters']['loaded'].append('bias')
        # Set and then delete the Master Bias frame
        self.SetMasterFrame(msbias, "bias", det)
        armasters.save_masters(self, det, mftype='bias')

        del msbias
        return True
        '''

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

    def MasterFlatField(self, fitsdict, det, msbias, datasec_img, tilts):
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
        msgs.error("SHOULD NOT GET HERE")
        '''
        dnum = settings.get_dnum(det)
        if settings.argflag['reduce']['flatfield']['perform']:  # Only do it if the user wants to flat field
            # If the master pixelflat is already made, use it
            if self._mspixelflat[det-1] is not None:
                msgs.info("An identical master pixelflat frame already exists")
                if self._mspixelflatnrm[det-1] is None:
                    # Normalize the flat field
                    msgs.info("Normalizing the pixel flat")
                    slit_profiles, mstracenrm, msblaze, flat_ext1d, extrap_slit = \
                        arflat.norm_slits(
                            self.GetMasterFrame("pixelflat", det), datasec_img, self._lordloc[det-1], self._rordloc[det-1],
                            self._pixwid[det-1], self._slitpix[det-1], det, tilts,
                            settings.argflag, settings.spect,
                            ntcky=settings.argflag['reduce']['flatfield']['params'][0])
                    #arflat.norm_slits(self, self.GetMasterFrame("pixelflat", det),
                    #                  det, ntcky=settings.argflag['reduce']['flatfield']['params'][0])
                    # If some slit profiles/blaze functions need to be extrapolated, do that now
                    if settings.spect['mosaic']['reduction'] == 'AMRED':
                        if np.sum(extrap_slit) != 0.0:
                            slit_profiles, mstracenrm, msblaze = arflat.slit_profile_pca(
                                self.GetMasterFrame("pixelflat", det),
                                tilts, msblaze, extrap_slit, slit_profiles,
                                self._lordloc[det-1], self._rordloc[det-1], self._pixwid[det-1],
                                self._slitpix[det-1], self.setup)
                    mspixelflatnrm = mstracenrm.copy()
                    winpp = np.where(slit_profiles != 0.0)
                    mspixelflatnrm[winpp] /= slit_profiles[winpp]
                    self.SetMasterFrame(mspixelflatnrm, "normpixelflat", det)
                    armasters.save_masters(self, det, mftype='normpixelflat')
                    if np.array_equal(self._idx_flat, self._idx_trace):
                        # The flat field frame is also being used to trace the slit edges and determine the slit
                        # profile. Avoid recalculating the slit profile and blaze function and save them here.
                        self.SetFrame(self._msblaze, msblaze, det)
                        self.SetFrame(self._slitprof, slit_profiles, det)
                        armasters.save_masters(self, det, mftype='slitprof')
                        if settings.argflag["reduce"]["slitprofile"]["perform"]:
                            msgs.info("Preparing QA of each slit profile")
#                            arqa.slit_profile(self, mstracenrm, slit_profiles, self._lordloc[det - 1], self._rordloc[det - 1],
#                                              self._slitpix[det - 1], desc="Slit profile")
                            arflat.slit_profile_qa(self, mstracenrm, slit_profiles,
                                                   self._lordloc[det - 1], self._rordloc[det - 1],
                                                   self._slitpix[det - 1], desc="Slit profile")
                        msgs.info("Saving blaze function QA")
#                        arqa.plot_orderfits(self, msblaze, flat_ext1d, desc="Blaze function")
                        artracewave.plot_orderfits(self, msblaze, flat_ext1d, desc="Blaze function")
                return False
            ###############
            # Generate/load a master pixel flat frame
            if settings.argflag['reduce']['flatfield']['useframe'] in ['pixelflat', 'trace']:
                mspixelflatnrm = armasters.load_master_frame(self, "normpixelflat")
                if mspixelflatnrm is None:
                    msgs.info("Preparing a master pixel flat frame with {0:s}".format(settings.argflag['reduce']['flatfield']['useframe']))
                    # Get all of the pixel flat frames for this science frame
                    ind = self._idx_flat
                    # Load the frames for tracing
                    frames = arload.load_frames(fitsdict, ind, det, frametype='pixel flat', msbias=msbias)
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
                    mspixelflat *= arprocimg.gain_frame(datasec_img, settings.spect[dnum]['numamplifiers'],
                                             settings.spect[dnum]['gain'])
                    # Normalize the flat field
                    msgs.info("Normalizing the pixel flat")
                    slit_profiles, mstracenrm, msblaze, flat_ext1d, extrap_slit = \
                        arflat.norm_slits(mspixelflat, datasec_img, self._lordloc[det-1], self._rordloc[det-1],
                                          self._pixwid[det-1], self._slitpix[det-1], det, tilts,
                                          settings.argflag, settings.spect,
                                          ntcky=settings.argflag['reduce']['flatfield']['params'][0])
                    # If some slit profiles/blaze functions need to be extrapolated, do that now
                    if settings.spect['mosaic']['reduction'] == 'AMRED':
                        if np.sum(extrap_slit) != 0.0:
                            slit_profiles, mstracenrm, msblaze = arflat.slit_profile_pca(
                                mspixelflat, tilts, msblaze, extrap_slit, slit_profiles,
                            self._lordloc[det-1], self._rordloc[det-1], self._pixwid[det-1],
                                self._slitpix[det-1], self.setup)
                    mspixelflatnrm = mstracenrm.copy()
                    winpp = np.where(slit_profiles != 0.0)
                    mspixelflatnrm[winpp] /= slit_profiles[winpp]
                    if np.array_equal(self._idx_flat, self._idx_trace):
                        # The flat field frame is also being used to trace the slit edges and determine the slit
                        # profile. Avoid recalculating the slit profile and blaze function and save them here.
                        self.SetFrame(self._msblaze, msblaze, det)
                        self.SetFrame(self._slitprof, slit_profiles, det)
                        armasters.save_masters(self, det, mftype='slitprof')
                        if settings.argflag["reduce"]["slitprofile"]["perform"]:
                            msgs.info("Preparing QA of each slit profile")
#                            arqa.slit_profile(self, mstracenrm, slit_profiles, self._lordloc[det - 1], self._rordloc[det - 1],
#                                              self._slitpix[det - 1], desc="Slit profile")
                            arflat.slit_profile_qa(mstracenrm, slit_profiles,
                                                   self._lordloc[det - 1], self._rordloc[det - 1],
                                                   self._slitpix[det - 1], desc="Slit profile",
                                                   setup=self.setup)
                        msgs.info("Saving blaze function QA")
#                        arqa.plot_orderfits(self, msblaze, flat_ext1d, desc="Blaze function")
                        artracewave.plot_orderfits(self.setup, msblaze, flat_ext1d, desc="Blaze function")
                else:
                    mspixelflat = mspixelflatnrm
            else:  # It must be the name of a file the user wishes to load
                mspixelflat_name = armasters.user_master_name(settings.argflag['run']['directory']['master'],
                                                              settings.argflag['reduce']['flatfield']['useframe'])
                mspixelflatnrm, head, _ = armasters._load(mspixelflat_name, exten=det, frametype=None,
                      force=settings.argflag['reduce']['masters']['force'])
                mspixelflat = mspixelflatnrm
            # Now that the combined, master flat field frame is loaded...
        else:
            msgs.work("Pixel Flat arrays need to be generated when not flat fielding")
            msgs.bug("Blaze is currently undefined")
            debugger.set_trace()
            mspixelflat = np.ones_like(self._msarc)
            mspixelflatnrm = np.ones_like(self._msarc)
        # Set Master Frames
        self.SetMasterFrame(mspixelflat, "pixelflat", det)
        self.SetMasterFrame(mspixelflatnrm, "normpixelflat", det)
        armasters.save_masters(self, det, mftype='normpixelflat')
        return True
        '''

    def MasterPinhole(self, fitsdict, det, msbias):
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
            try:
                mspinhole = armasters.get_master_frame(self, "pinhole")
            except IOError:
                msgs.info("Preparing a master pinhole frame with {0:s}".format(
                    settings.argflag['reduce']['slitcen']['useframe']))
                ind = self._idx_cent
                # Load the pinhole frames
                frames = arload.load_frames(fitsdict, ind, det, frametype='pinhole', msbias=msbias, # self._msbias[det - 1],
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
            mspinhole, head = armasters.load_master(mspinhole_name, frametype=None)
            debugger.set_trace()  # NEED TO LOAD EXTRAS AS ABOVE
        # Set and then delete the Master Trace frame
        self.SetMasterFrame(mspinhole, "pinhole", det)
        #armasters.save_masters(self, det, mftype='pinhole')
        del mspinhole
        return True

    '''
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
        # DEPRECATED
        msgs.error("DEPRECATED!")
        dnum = settings.get_dnum(det)
        # If the master trace is already made, use it
        if self._mstrace[det-1] is not None:
            msgs.info("An identical master trace frame already exists")
            return False
        if settings.argflag['reduce']['trace']['useframe'] in ['trace']:
            mstrace = armasters.load_master_frame(self, "trace", det=det)  # Also loads up the various arrays
            if mstrace is None:
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
            mstrace, _ = armasters.load_master(mstrace_name, frametype=None)
            debugger.set_trace()  # NEED TO LOAD EXTRAS AS ABOVE;  SHOULD MODIFY get_master_frame()
        # Set and then delete the Master Trace frame
        self.SetMasterFrame(mstrace, "trace", det)
        del mstrace
        return True
    '''

    def MasterWave(self, det, all_wvcalib, tilts):
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
        msgs.error("SHOULD NOT BE HERE ANYMORE")
        '''
        if self._mswave[det-1] is not None:
            msgs.info("An identical master arc frame already exists")
            return False
        mswave = armasters.load_master_frame(self, "wave")
        if mswave is None:
            msgs.info("Preparing a master wave frame")
            if settings.argflag["reduce"]["calibrate"]["wavelength"] == "pixel":
                mswave = tilts * (tilts.shape[0]-1.0)
            else:
                ok_slits = np.where(~self._maskslits[det-1])[0]
                mswave = np.zeros_like(tilts)
                for slit in ok_slits:
                    iwv_calib = all_wvcalib[str(slit)]
                    tmpwv = arutils.func_val(iwv_calib['fitc'], tilts, iwv_calib['function'],
                                          minv=iwv_calib['fmin'], maxv=iwv_calib['fmax'])
                    word = np.where(self._slitpix[det - 1] == slit+1)
                    mswave[word] = tmpwv[word]
        # Set and then delete the Master Arc frame
        self.SetMasterFrame(mswave, "wave", det)
        armasters.save_masters(self, det, mftype='wave')
        del mswave
        return True
        '''

    '''
    def MasterWaveCalib(self, fitstbl, det, msarc):
        """
        Generate Master 1D Wave Solution (down slit/order centers)

        Parameters
        ----------
        fitstbl : dict
          Contains relevant information from fits header files
        det : int
          Index of the detector

        Returns
        -------
        boolean : bool
          Should other ScienceExposure classes be updated?
        """

        if self._wvcalib[det-1] is not None:
            msgs.info("An identical master wave calib frame already exists")
            return False
        else:
            wv_calib = None
        # Attempt to load the Master Frame
        wv_calib = armasters.load_master_frame(self, "wv_calib")
        if wv_calib is None:
            if settings.argflag["reduce"]["calibrate"]["wavelength"] == "pixel":
                msgs.info("A wavelength calibration will not be performed")
            else:
                # Setup arc parameters (e.g. linelist)
                arc_idx = arsort.ftype_indices(fitstbl, 'arc', self.sci_ID)
                arcparam = ararc.setup_param(msarc.shape, fitstbl, arc_idx[0])
                self.SetFrame(self._arcparam, arcparam, det)
                ###############
                # Extract an arc down each slit
                arccen, maskslit, _ = artrace.get_censpec(self, msarc, det, gen_satmask=False)
                ok_mask = np.where(maskslit == 0)[0]

                # Fill up the calibrations
                wv_calib = {}
                for kk,slit in enumerate(ok_mask):
                    ###############
                    # Extract arc and identify lines
                    if settings.argflag['arc']['calibrate']['method'] == 'simple':
                        iwv_calib = ararc.simple_calib(self, det, msarc, censpec=arccen[:,kk], slit=slit)
                    elif settings.argflag['arc']['calibrate']['method'] == 'arclines':
                        iwv_calib = ararc.calib_with_arclines(self, det, msarc, slit, arcparam, censpec=arccen[:,kk])
                    wv_calib[str(slit)] = iwv_calib.copy()
        # Set
        if wv_calib is not None:
            self.SetFrame(self._wvcalib, wv_calib, det)
            armasters.save_masters(self, det, mftype='wv_calib')

            # Set mask based on wv_calib existing
            mask = np.array([True]*self._lordloc[det-1].shape[1])
            for key in self._wvcalib[det-1].keys():
                mask[int(key)] = False
            self._maskslits[det-1] = mask
            #
            del wv_calib
        return True
    '''

    def MasterStandard(self, fitsdict, msbias):
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
        msgs.error("THIS IS DEPRECATED")
        '''
        if self._sensfunc is not None:
            msgs.info("Using existing sensitivity function.")
            return False
        # Attempt to load the Master Frame
        try:
            sensfunc = armasters.get_master_frame(self, "sensfunc")
        except IOError:
            # Grab the standard star frames
            msgs.info("Preparing the standard")
            ind = self._idx_std
            msgs.warn("Taking only the first standard frame for now")
            ind = [ind[0]]
            # Extract
            all_specobj = []
            for kk in range(settings.spect['mosaic']['ndet']):
                det = kk+1
                # Use this detector? Need to check this after setting RA/DEC above
                if settings.argflag['reduce']['detnum'] is not None:
                    msgs.warn("If your standard wasnt on this detector, you will have trouble..")
                    if det not in map(int, settings.argflag['reduce']['detnum']):
                        continue
                # Load the frame(s)
                frame = arload.load_frames(fitsdict, ind, det, frametype='standard',
                                           msbias=msbias) # self._msbias[det-1])
                sciframe = frame[:, :, 0] # First exposure
                # Save RA/DEC
                self._msstd[0]['RA'] = fitsdict['ra'][ind[0]]    # Yes, this needs to be index 0
                self._msstd[0]['DEC'] = fitsdict['dec'][ind[0]]  # Yes, this needs to be index 0
                self._msstd[det-1]['spobjs'] = None
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
        else:
            self._sensfunc = sensfunc.copy()
            return True
        '''

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
        #if ftype == "arc": self._msarc[det] = cpf
        if ftype == "wave": self._mswave[det] = cpf
        #elif ftype == "bias": self._msbias[det] = cpf
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
            #if ftype == "arc": return self._msarc[det].copy()
            if ftype == "wave": return self._mswave[det].copy()
            #elif ftype == "bias": return self._msbias[det].copy()
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
            #if ftype == "arc": return self._msarc[det]
            if ftype == "wave": return self._mswave[det]
            #elif ftype == "bias": return self._msbias[det]
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


def dummy_self(inum=1, fitstbl=None, nfile=10):
    """ Generate a dummy self class for testing
    Parameters:
    -----------
    inum : int, optional
      Index in sciexp
    Returns:
    --------
    slf
    """
    # Dummy fitsdict
    if fitstbl is None:
        fitstbl = arsort.dummy_fitstbl(nfile=nfile)
    # Dummy Class
    return ScienceExposure(inum, fitstbl, settings.argflag, settings.spect, do_qa=False)




