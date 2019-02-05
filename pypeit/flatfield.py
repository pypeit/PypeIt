# Module for generating the Flat Field
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

#from importlib import reload

from pypeit import msgs

from pypeit import processimages
#from pypeit.core import masters
from pypeit import masterframe
from pypeit.core import flat
from pypeit import ginga
from pypeit.par import pypeitpar
from pypeit.core import pixels
from pypeit.core import trace_slits
from pypeit.core import tracewave


from pypeit import debugger

class FlatField(processimages.ProcessImages, masterframe.MasterFrame):
    """
    This class will generate the pixel-level FlatField
      The master() method returns the image

    Parameters
    ----------
    spectrograph : str or Spectrograph
    file_list : list
      List of raw files to produce the flat field
    settings : dict-like
    msbias : ndarray or str or None
    tslits_dict : dict
      dict from TraceSlits class (e.g. slitpix)
    tilts : ndarray
      tilts from WaveTilts class
    det : int
    master_key : str

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

    def __init__(self, spectrograph, files=None, binning=None, det=1, par=None, master_key=None,
                 master_dir=None, reuse_masters=False, flatpar=None, msbias=None, msbpm=None,
                 tslits_dict=None, tilts_dict=None):

        # Image processing parameters
        self.par = pypeitpar.FrameGroupPar(self.frametype) if par is None else par

        # Start us up
        processimages.ProcessImages.__init__(self, spectrograph, files=files, det=det,
                                             par=self.par['process'])

        # MasterFrames: Specifically pass the ProcessImages-constructed
        # spectrograph even though it really only needs the string name
        masterframe.MasterFrame.__init__(self, self.frametype, master_key,
                                         master_dir=master_dir, reuse_masters=reuse_masters)

        # Parameters unique to this Object
        self.msbias = msbias
        self.tslits_dict = tslits_dict
        self.tilts_dict = tilts_dict
        self.msbpm = msbpm
        self.binning = binning
        if master_dir is None:
            self.master_dir = os.getcwd()
        else:
            self.master_dir = master_dir

        # Key outputs
        self.rawflatimg = None
        self.mspixelflat = None
        self.msillumflat = None
        self.flat_model = None

        # Child-specific Internals
        self.extrap_slit = None
        self.msblaze = None
        self.blazeext = None
        self.slit_profiles = None
        self.ntckx = None
        self.ntcky = None

        # FieldFlattening parameters
        self.flatpar = pypeitpar.FlatFieldPar() if flatpar is None else flatpar


    @property
    def nslits(self):
        """
        Number of slits

        Returns
        -------
        nslits : int

        """
        if self.tslits_dict is not None:
            return self.tslits_dict['slit_left'].shape[1]
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
        self.mspixelflat = self.process(bias_subtract=self.msbias, bpm = self.msbpm, trim=trim, apply_gain=True)
        # Step
        self.steps.append(inspect.stack()[0][3])
        #
        return self.mspixelflat

    def _prep_tck(self):
        """
        Setup for the bspline fitting

        Wrapper to flat.prep_ntck

        Returns
        -------
        self.ntckx -- set internally
          Number of knots in the spatial dimension
        self.ntcky -- set internally
          Number of knots in the spectral dimension
        """
        # Step
        self.steps.append(inspect.stack()[0][3])
        self.ntckx, self.ntcky = flat.prep_ntck(self.tslits_dict['pixwid'],
                                                  method=self.flatpar['method'],
                                                  params=self.flatpar['params'],
                                                  get_slitprofile=self.flatpar['slitprofile'])

    # ToDO JFH:
    # This load_master_illumflat code is a bit of a kludge. Usually one reads in masters with the master and load_master method, but there are
    # technically two master files for the flats, i.e. a pixelflat and illumination flat. Perhaps the better way to deal with this
    # would be to package them into one output file and just change the load_master and save_master methods to deal with the
    # possible existence of an illumination flat
    def load_master_illumflat(self, force=False):
        """
        Load the slit illumination profile from a saved Master file

        Returns
        -------
        self.slit_profiles

        """
        ms_name = masterframe.master_name('illumflat', self.master_key, self.mdir)
        msframe = self.load_master(ms_name)
        if msframe is None:
            msgs.warn("No Master frame found of type {:s}: {:s}".format('illumflat', ms_name))

        return msframe

    # TODO Need to add functionality to use a different frame for the ilumination flat, e.g. a sky flat
    def run(self, debug=False, show=False):
        """
        Main driver to generate normalized flat field and illumination flats
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

        Returns
        -------
        self.mspixelflatnrm
        self.slit_profiles

        """

        # Build the pixel flat (as needed)
        if self.rawflatimg is None:
            self.rawflatimg = self.build_pixflat()

        # Prep tck (sets self.ntckx, self.ntcky)
        #self._prep_tck()

        # Setup
        self.mspixelflat = np.ones_like(self.rawflatimg)
        self.msillumflat = np.ones_like(self.rawflatimg)
        self.flat_model = np.zeros_like(self.rawflatimg)
        self.slitmask = pixels.tslits2mask(self.tslits_dict)


        final_tilts = np.zeros_like(self.rawflatimg)

        # Loop on slits
        for slit in range(self.nslits):
            msgs.info('Computing flat field image for slit: {:d}/{:d}'.format(slit,self.nslits-1))
            if self.msbpm is not None:
                inmask = ~self.msbpm
            else:
                inmask = np.ones_like(self.rawflatimg,dtype=bool)

            # Fit flats for a single slit
            this_tilts_dict = {'tilts':self.tilts_dict['tilts'], 'coeffs':self.tilts_dict['coeffs'][:,:,slit].copy(),
                               'slitcen':self.tilts_dict['slitcen'][:,slit].copy(),
                               'func2d':self.tilts_dict['func2d']}
            nonlinear_counts = self.spectrograph.detector[self.det - 1]['nonlinear']*\
                               self.spectrograph.detector[self.det - 1]['saturation']
            pixelflat, illumflat, flat_model, tilts_out, thismask_out, slit_left_out, slit_righ_out = \
                flat.fit_flat(self.rawflatimg, this_tilts_dict, self.tslits_dict, slit,
                              inmask=inmask,
                              nonlinear_counts=nonlinear_counts,
                              spec_samp_fine=self.flatpar['spec_samp_fine'], spec_samp_coarse=self.flatpar['spec_samp_coarse'],
                              spat_samp=self.flatpar['spat_samp'], tweak_slits=self.flatpar['tweak_slits'],
                              tweak_slits_thresh=self.flatpar['tweak_slits_thresh'],
                              tweak_slits_maxfrac=self.flatpar['tweak_slits_maxfrac'],debug=debug)
            self.mspixelflat[thismask_out] = pixelflat[thismask_out]
            self.msillumflat[thismask_out] = illumflat[thismask_out]
            self.flat_model[thismask_out] = flat_model[thismask_out]
            # Did we tweak slit boundaries? If so, update the tslits_dict and the tilts_dict
            if self.flatpar['tweak_slits']:
                self.tslits_dict['slit_left'][:, slit] = slit_left_out
                self.tslits_dict['slit_righ'][:, slit] = slit_righ_out
                final_tilts[thismask_out] = tilts_out[thismask_out]

        # If we tweaked the slits update the tilts_dict
        if self.flatpar['tweak_slits']:
            self.tilts_dict['tilts'] = final_tilts

        if show:
            # Global skysub is the first step in a new extraction so clear the channels here
            self.show(slits=True, wcs_match = True)

        # If illumination flat fielding is turned off, set the illumflat to be None.
        if not self.flatpar['illumflatten']:
            msgs.warn('You have set illumflatten=False. No illumination flat will be applied to your data.')
            self.msillumflat = None


        # Return
        return self.mspixelflat, self.msillumflat



    def show(self, slits = True, wcs_match = True):

        viewer, ch = ginga.show_image(self.mspixelflat, chname='pixeflat', cuts=(0.9, 1.1), wcs_match=wcs_match, clear=True)
        viewer, ch = ginga.show_image(self.msillumflat, chname='illumflat', cuts=(0.9, 1.1), wcs_match=wcs_match)
        viewer, ch = ginga.show_image(self.rawflatimg, chname='flat', wcs_match=wcs_match)
        viewer, ch = ginga.show_image(self.flat_model, chname='flat_model', wcs_match=wcs_match)


        if slits:
            if self.tslits_dict is not None:
                slit_ids = [trace_slits.get_slitid(self.rawflatimg.shape, self.tslits_dict['slit_left'], self.tslits_dict['slit_righ'], ii)[0]
                            for ii in range(self.tslits_dict['slit_left'].shape[1])]
                ginga.show_slits(viewer, ch, self.tslits_dict['slit_left'], self.tslits_dict['slit_righ'], slit_ids)


