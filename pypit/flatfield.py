# Module for generating the Flat Field
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np

#from importlib import reload

from pypit import msgs
from pypit import processimages
from pypit.core import armasters
from pypit import masterframe
from pypit.core import arsort
from pypit.core import arflat
from pypit import ginga
from pypit.par import pypitpar

from pypit import ardebug as debugger

# TODO, JFH: I do not understand why ArcImage is its own class whereas
# there is no class called FlatImage. In principle this is because
# several codes like tilts, and wv_calib use the msarc so it made sense
# to make reading the image a separate class. However, that is also true
# for flats, since the flat fielding code and the tracing code both use
# the flat image
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

    # Frame type is a class attribute
    frametype = 'pixelflat'

    def __init__(self, spectrograph, file_list=[], det=1, par=None, setup=None, root_path=None,
                 mode=None, flatpar=None, msbias=None, tslits_dict=None, tilts=None):

        # Image processing parameters
        self.par = pypitpar.FrameGroupPar(self.frametype) if par is None else par

        # Start us up
        processimages.ProcessImages.__init__(self, spectrograph, file_list=file_list, det=det,
                                             par=self.par['process'])

        # MasterFrames: Specifically pass the ProcessImages-constructed
        # spectrograph even though it really only needs the string name
        directory_path = None if root_path is None \
                                else root_path+'_'+self.spectrograph.spectrograph
        masterframe.MasterFrame.__init__(self, self.frametype, setup,
                                         directory_path=directory_path, mode=mode)

        # Parameters unique to this Object
        self.msbias = msbias
        self.tslits_dict = tslits_dict
        self.tilts = tilts

        # Key outputs
        self.mspixelflat = None
        self.mspixelflatnrm = None

        # Child-specific Internals
        self.extrap_slit = None
        self.msblaze = None
        self.blazeext = None
        self.slit_profiles = None
        self.ntckx = None
        self.ntcky = None

        # FieldFlattening parameters
        self.flatpar = pypitpar.FlatFieldPar() if flatpar is None else flatpar


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
        self.mspixelflat (points at self.stack)

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
        self.ntckx, self.ntcky = arflat.prep_ntck(self.tslits_dict['pixwid'],
                                                  method=self.flatpar['method'],
                                                  params=self.flatpar['params'],
                                                  get_slitprofile=self.flatpar['slitprofile'])

    def load_master_slitprofile(self):
        """
        Load the slit profile from a saved Master file

        Returns
        -------
        self.slit_profiles

        """
        return armasters.load_master_frame('slitprof', self.setup, self.mdir)

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
        modvals, nrmvals, msblaze_slit, blazeext_slit, iextrap_slit \
                = arflat.slit_profile(slit, self.mspixelflat, self.tilts, slordloc, srordloc,
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
            if self.flatpar['slitprofile']:
                # Leave slit_profiles as ones if the slitprofile is not
                # being determined, otherwise, set the model.
                self.slit_profiles[word] = modvals/nrmvals
            self.mspixelflatnrm[word] /= nrmvals
            # Fill
            self.msblaze[:,slit] = msblaze_slit
            self.blazeext[:,slit] = blazeext_slit

        # If some slit profiles/blaze functions need to be extrapolated, do that now
        if armed:
            if np.sum(self.extrap_slit) != 0.0:
                slit_profiles, mstracenrm, msblaze \
                            = arflat.slit_profile_pca(self.mspixelflat, self.tilts, self.msblaze,
                                                      self.extrap_slit, self.slit_profiles,
                                                      self.tslits_dict['lcen'],
                                                      self.tslits_dict['rcen'],
                                                      self.tslits_dict['pixwid'],
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

    def show(self, attr='norm', display='ginga'):
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


