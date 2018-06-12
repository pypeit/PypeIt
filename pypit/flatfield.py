# Module for generating the Flat Field
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

#from importlib import reload

from pypit import msgs
from pypit import processimages
from pypit import armasters
from pypit import masterframe
from pypit.core import arprocimg
from pypit.core import arflat
from pypit import ginga

from pypit import ardebug as debugger

# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

# Does not need to be global, but I prefer it
frametype = 'pixelflat'

def default_settings():
    default_settings = dict(flatfield={'method': 'bspline',
                                   "params": [20],
                                   },
                        slitprofile={'perform': True,
                                     },
                        )
    return default_settings

class FlatField(processimages.ProcessImages, masterframe.MasterFrame):
    """
    This class will generat the pixel-level FlatField
      The master() method returns the image

    Parameters
    ----------
    file_list : list
      List of raw files to produce the flat field
    spectrograph : str
    settings : dict-like
    msbias : ndarray or str or None
    tslits_dict : dict
      dict from TraceSlits class (e.g. slitpix)
    tilts : ndarray
      tilts from WaveTilts class
    det : int
    setup : str

    Attributes
    ----------
    frametype : str
      Set to 'pixelflat'
    mspixelflat : ndarray
      Stacked image
    mspixelflatnrm : ndarray
      Normalized flat
    extrap_slit
    msblaze : ndarray
      Blaze function fit to normalize
    blazeext :
    slit_profiles : ndarray
      Slit profile(s)
    self.ntckx : int
      Number of knots in the spatial dimension
    self.ntcky : int
      Number of knots in the spectral dimension

    """
    # Keep order same as processimages (or else!)
    def __init__(self, file_list=[], spectrograph=None, settings=None, msbias=None,
                 tslits_dict=None, tilts=None, det=None, setup=None, datasec_img=None):

        # Parameters unique to this Object
        self.msbias = msbias
        self.det = det
        self.setup = setup
        self.tslits_dict = tslits_dict
        self.tilts = tilts

        # Start us up
        processimages.ProcessImages.__init__(self, file_list, spectrograph=spectrograph,
                                             settings=settings, det=det, datasec_img=datasec_img)

        # Attributes (set after init)
        self.frametype = frametype

        # Key outputs
        self.mspixelflat = None
        self.mspixelflatnrm = None

        # Settings
        # The copy allows up to update settings with user settings without changing the original
        if settings is None:
            # Defaults
            self.settings = processimages.default_settings()
        else:
            self.settings = settings.copy()
            # The following is somewhat kludgy and the current way we do settings may
            #   not touch all the options (not sure .update() would help)
            if 'combine' not in settings.keys():
                if self.frametype in settings.keys():
                    self.settings['combine'] = settings[self.frametype]['combine']
        if 'flatfield' not in settings.keys():
            self.settings.update(default_settings())

        # Child-specific Internals
        #    See ProcessImages and MasterFrame for others
        self.extrap_slit = None
        self.msblaze = None
        self.blazeext = None
        self.slit_profiles = None
        self.ntckx = None
        self.ntcky = None

        # MasterFrames
        masterframe.MasterFrame.__init__(self, self.frametype, self.setup, self.settings)

    @property
    def nslits(self):
        """
        Number of slits

        Returns
        -------
        nslits : int

        """
        if self.tslits_dict is not None:
            return self.tslits_dict['lcen'].shape[1]
        else:
            return 0

    def build_pixflat(self, trim=True):
        """
        # Generate the flat image

        Parameters
        ----------
        trim : bool, optional

        Returns
        -------
        self.mspixelflat

        """
        self.mspixelflat = self.process(bias_subtract=self.msbias, trim=trim, apply_gain=True)
        # Step
        self.steps.append(inspect.stack()[0][3])
        #
        return self.mspixelflat

    def _prep_tck(self):
        """
        Setup for the bspline fitting

        Wrapper to arflat.prep_ntck

        Returns
        -------
        self.ntckx -- set internally
          Number of knots in the spatial dimension
        self.ntcky -- set internally
          Number of knots in the spectral dimension
        """
        # Step
        self.steps.append(inspect.stack()[0][3])
        self.ntckx, self.ntcky = arflat.prep_ntck(self.tslits_dict['pixwid'], self.settings)

    def load_master_slitprofile(self):
        """
        Load the slit profile from a saved Master file

        Returns
        -------
        self.slit_profiles

        """
        return armasters.core_load_master_frame('slitprof', self.setup, self.mdir)

    def slit_profile(self, slit):
        """
        Generate the slit profile for a given slit

        Wrapper to arflat.slit_profile()

        Parameters
        ----------
        slit : int

        Returns
        -------
        modvals : ndarray
          Pixels in the slit
        nrmvals : ndarray
          Pixels in the slit
        msblaze_slit : ndarray (nwave)
        blazeext_slit : ndarray (nwave)
        iextrap_slit : float
          0 = Do not extrapolate
          1 = Do extrapolate

        """
        # Check
        if self.ntckx is None:
            msgs.warn("Need to set self.ntckx with _prep_tck first!")
            return [None]*5
        # Wrap me
        slordloc = self.tslits_dict['lcen'][:,slit]
        srordloc = self.tslits_dict['rcen'][:,slit]
        modvals, nrmvals, msblaze_slit, blazeext_slit, iextrap_slit = arflat.slit_profile(
            slit, self.mspixelflat, self.tilts, slordloc, srordloc,
            self.tslits_dict['slitpix'], self.tslits_dict['pixwid'],
            ntckx=self.ntckx, ntcky=self.ntcky)
        # Step
        step = inspect.stack()[0][3]
        if step not in self.steps:  # Only do it once
            self.steps.append(step)
        # Return
        return modvals, nrmvals, msblaze_slit, blazeext_slit, iextrap_slit

    def run(self, armed=False):
        """
        Main driver to generate normalized flat field

        Code flow:
          1.  Generate the pixelflat image (if necessary)
          2.  Prepare b-spline knot spacing
          3.  Loop on slits/orders
             a. Calculate the slit profile
             b. Normalize
             c. Save

        Parameters
        ----------
        datasec_img
        armed : bool, optional

        Returns
        -------
        self.mspixelflatnrm
        self.slit_profiles

        """

        # Build the pixel flat (as needed)
        if self.mspixelflat is None:
            self.mspixelflat = self.build_pixflat()

        # Prep tck (sets self.ntckx, self.ntcky)
        self._prep_tck()

        # Setup
        self.extrap_slit = np.zeros(self.nslits, dtype=np.int)
        self.mspixelflatnrm = self.mspixelflat.copy()
        self.msblaze = np.ones_like(self.tslits_dict['lcen'])
        self.blazeext = np.ones_like(self.tslits_dict['lcen'])
        self.slit_profiles = np.ones_like(self.mspixelflat)

        # Loop on slits
        for slit in range(self.nslits):
            # Normalize a single slit
            modvals, nrmvals, msblaze_slit, blazeext_slit, iextrap_slit = self.slit_profile(slit)

            word = np.where(self.tslits_dict['slitpix'] == slit+1)
            self.extrap_slit[slit] = iextrap_slit
            if modvals is None:
                continue
            # Fill
            if self.settings["slitprofile"]["perform"]:
                # Leave slit_profiles as ones if the slitprofile is not being determined, otherwise, set the model.
                self.slit_profiles[word] = modvals/nrmvals
            self.mspixelflatnrm[word] /= nrmvals
            # Fill
            self.msblaze[:,slit] = msblaze_slit
            self.blazeext[:,slit] = blazeext_slit

        # If some slit profiles/blaze functions need to be extrapolated, do that now
        if armed:
            if np.sum(self.extrap_slit) != 0.0:
                slit_profiles, mstracenrm, msblaze = arflat.slit_profile_pca(
                    self.mspixelflat, self.tilts, self.msblaze, self.extrap_slit, self.slit_profiles,
                    self.tslits_dict['lcen'], self.tslits_dict['rcen'], self.tslits_dict['pixwid'],
                    self.tslits_dict['slitpix'])

        # Apply slit profile
        winpp = np.where(self.slit_profiles != 0.0)
        self.mspixelflatnrm[winpp] /= self.slit_profiles[winpp]

        # Set pixels not in slits to 1.
        msgs.info("Setting pixels outside of slits to 1. in the flat.")
        inslit = self.tslits_dict['slitpix'] >= 1.
        self.mspixelflatnrm[~inslit] = 1.

        # QA
        '''
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
                arproc.slit_profile_qa(self, mstracenrm, slit_profiles,
                                       self._lordloc[det - 1], self._rordloc[det - 1],
                                       self._slitpix[det - 1], desc="Slit profile")
            msgs.info("Saving blaze function QA")
            #                        arqa.plot_orderfits(self, msblaze, flat_ext1d, desc="Blaze function")
            artracewave.plot_orderfits(self.setup, msblaze, flat_ext1d, desc="Blaze function")
        #
        '''

        # Return
        return self.mspixelflatnrm, self.slit_profiles

    def show(self, attr, display='ginga'):
        """
        Show one of the internal images

        Parameters
        ----------
        attr : str
          mspixelflat -- Show the combined flat image, unnormalized
          norm -- Show the combined normalized flat image
        display : str, optional

        Returns
        -------

        """
        if attr == 'mspixelflat':
            if self.mspixelflat is not None:
                ginga.show_image(self.mspixelflat)
        elif attr == 'norm':
            if self.mspixelflatnrm is not None:
                ginga.show_image(self.mspixelflatnrm)
