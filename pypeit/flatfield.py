"""
Implements the flat-field class.

.. _numpy.ndarray: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
"""
import inspect
import numpy as np
import os

from pypeit import msgs

from pypeit import processimages
from pypeit import masterframe
from pypeit.core import flat
from pypeit import ginga
from pypeit.par import pypeitpar
from pypeit.core import pixels
from pypeit.core import trace_slits

from pypeit import debugger

class FlatField(processimages.ProcessImages, masterframe.MasterFrame):
    """
    Builds pixel-level flat-field and the illumination flat-field.

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The `Spectrograph` instance that sets the instrument used to
            take the observations.  See usage by
            :class:`pypeit.processimages.ProcessImages` base class.
        par (:class:`pypeit.par.pypeitpar.FrameGroupPar`):
            The parameters used to type and process the flat frames.
        files (:obj:`list`, optional):
            The list of files to process.  Can be an empty list.
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.
        master_key (:obj:`str`, optional):
            The string identifier for the instrument configuration.  See
            :class:`pypeit.masterframe.MasterFrame`.
        master_dir (:obj:`str`, optional):
            Path to master frames
        msbias (`numpy.ndarray`_, :obj:`str`, optional):
            Either an image with the bias to be subtracted, or a string
            providing the method to use for bias correction.
        msbpm (`numpy.ndarray`_, optional):
            Bad pixel mask image
        flatpar (:class:`pypeit.par.pypeitpar.FlatFieldPar`, optional):
            User-level parameters for constructing the flat-field
            corrections.  If None, the default parameters are used.
        tslits_dict (:obj:`dict`):
            The current slit traces; see
            :class:`pypeit.traceslits.TraceSlits`.
        tilts_dict (:obj:`dict`, optional):
            The current wavelength tilt traces; see
            :class:`pypeit.wavetilts.WaveTilts`.
        reuse_masters (:obj:`bool`, optional):
            Reuse already created master files from disk.

    Attributes:
        TODO: Fill this in...
        mspixelflat (ndarray): Stacked image
        mspixelflatnrm (ndarray): Normalized flat
        ntckx (int): Number of knots in the spatial dimension
        ntcky (int): Number of knots in the spectral dimension
    """

    # Frame type is a class attribute
    frametype = 'pixelflat'
    master_type = 'Flat'

    def __init__(self, spectrograph, par, files=None, det=1, master_key=None,
                 master_dir=None, reuse_masters=False, flatpar=None, msbias=None, msbpm=None,
                 tslits_dict=None, tilts_dict=None):

        # Image processing parameters
        if not isinstance(par, pypeitpar.FrameGroupPar):
            msgs.error('Must provide a FrameGroupPar instance as the parameters for FlatField.')
        self.par = par

        # Instantiate the base classes
        #   - Basic processing of the raw images
        processimages.ProcessImages.__init__(self, spectrograph, par=self.par['process'],
                                             files=files, det=det)
        #   - Construction and interface as a master frame
        masterframe.MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
                                         master_key=master_key, reuse_masters=reuse_masters)

        # FieldFlattening parameters
        self.flatpar = pypeitpar.FlatFieldPar() if flatpar is None else flatpar

        # Input master frame data
        self.msbias = msbias
        self.msbpm = msbpm
        self.tslits_dict = tslits_dict
        self.tilts_dict = tilts_dict

        # Parameters unique to this Object
        self.rawflatimg = None      # Un-normalized pixel flat
        self.mspixelflat = None     # Normalized pixel flat
        self.msillumflat = None     # Slit illumination flat
        self.flat_model = None      # Model flat

        # Completed steps
        self.steps = []

        # Child-specific Internals
#        self.extrap_slit = None
#        self.msblaze = None
#        self.blazeext = None
#        self.slit_profiles = None
#        self.ntckx = None
#        self.ntcky = None

    @property
    def nslits(self):
        """
        Number of slits from the :attr:`tslits_dict` key=slit_left

        Returns:
            int: Number of slits

        """
        if self.tslits_dict is not None:
            return self.tslits_dict['slit_left'].shape[1]
        else:
            return 0

    def build_pixflat(self, trim=True, force=False):
        """
        Generate the flat image.

        Args:
            trim (:obj:`bool`, optional):
                Trim the image down to just the data section.
            force (:obj:`bool`, optional):
                Force the flat to be reconstructed if it already exists

        Returns:
            `numpy.ndarray`_:  The image with the unnormalized
            pixel-flat data.
        """
        if self.rawflatimg is None or force:
            self.rawflatimg = self.process(bias_subtract=self.msbias, bpm=self.msbpm, trim=trim,
                                           apply_gain=True)
            self.steps.append(inspect.stack()[0][3])
        return self.rawflatimg

#    def _prep_tck(self):
#        """
#        Setup for the bspline fitting
#
#        Wrapper to flat.prep_ntck
#
#        :attr:`ntckx` -- Number of knots in the spatial dimension
#        :attr:`ntcky` -- Number of knots in the spectral dimension
#
#        """
#        # Step
#        self.steps.append(inspect.stack()[0][3])
#        self.ntckx, self.ntcky = flat.prep_ntck(self.tslits_dict['pixwid'],
#                                                  method=self.flatpar['method'],
#                                                  params=self.flatpar['params'],
#                                                  get_slitprofile=self.flatpar['slitprofile'])

    # TODO Need to add functionality to use a different frame for the ilumination flat, e.g. a sky flat
    def run(self, debug=False, show=False):
        """
        Generate normalized pixel and illumination flats

        Code flow::
           1.  Generate the pixelflat image (if necessary)
           2.  Prepare b-spline knot spacing
           3.  Loop on slits/orders
               a. Calculate the slit profile
               b. Normalize
               c. Save

        Args:
            debug (:obj:`bool`, optional):
                Run in debug mode.
            show (:obj:`bool`, optional):
                Show the results in the ginga viewer.

        Returns:
            `numpy.ndarray`_: Two arrays are returned, the normalized
            pixel flat data and the slit illumination correction data.
        """
        # Build the pixel flat
        self.build_pixflat()

        # Prep tck (sets self.ntckx, self.ntcky)
        #self._prep_tck()

        # Setup
        self.mspixelflat = np.ones_like(self.rawflatimg)
        self.msillumflat = np.ones_like(self.rawflatimg)
        self.flat_model = np.zeros_like(self.rawflatimg)
        self.slitmask = pixels.tslits2mask(self.tslits_dict)

        final_tilts = np.zeros_like(self.rawflatimg)

        # If we are tweaking slits allocate the new aray to hold tweaked slit boundaries
        if self.flatpar['tweak_slits']:
            self.tslits_dict['slit_left_tweak'] = np.zeros_like(self.tslits_dict['slit_left'])
            self.tslits_dict['slit_righ_tweak'] = np.zeros_like(self.tslits_dict['slit_righ'])

        # Loop on slits
        for slit in range(self.nslits):
            msgs.info('Computing flat field image for slit: {:d}/{:d}'.format(slit,self.nslits-1))
            if self.msbpm is not None:
                inmask = np.invert(self.msbpm)
            else:
                inmask = np.ones_like(self.rawflatimg,dtype=bool)

            # Fit flats for a single slit
            this_tilts_dict = {'tilts':self.tilts_dict['tilts'],
                               'coeffs':self.tilts_dict['coeffs'][:,:,slit].copy(),
                               'slitcen':self.tilts_dict['slitcen'][:,slit].copy(),
                               'func2d':self.tilts_dict['func2d']}
            nonlinear_counts = self.spectrograph.nonlinear_counts(det=self.det)

            pixelflat, illumflat, flat_model, tilts_out, thismask_out, slit_left_out, \
                    slit_righ_out \
                            = flat.fit_flat(self.rawflatimg, this_tilts_dict, self.tslits_dict,
                                           slit, inmask=inmask, nonlinear_counts=nonlinear_counts,
                                           spec_samp_fine=self.flatpar['spec_samp_fine'],
                                           spec_samp_coarse=self.flatpar['spec_samp_coarse'],
                                           spat_samp=self.flatpar['spat_samp'],
                                           tweak_slits=self.flatpar['tweak_slits'],
                                           tweak_slits_thresh=self.flatpar['tweak_slits_thresh'],
                                           tweak_slits_maxfrac=self.flatpar['tweak_slits_maxfrac'],
                                           debug=debug)

            self.mspixelflat[thismask_out] = pixelflat[thismask_out]
            self.msillumflat[thismask_out] = illumflat[thismask_out]
            self.flat_model[thismask_out] = flat_model[thismask_out]

            # Did we tweak slit boundaries? If so, update the tslits_dict and the tilts_dict
            if self.flatpar['tweak_slits']:
                self.tslits_dict['slit_left'][:,slit] = slit_left_out
                self.tslits_dict['slit_righ'][:,slit] = slit_righ_out
                self.tslits_dict['slit_left_tweak'][:,slit] = slit_left_out
                self.tslits_dict['slit_righ_tweak'][:,slit] = slit_righ_out
                final_tilts[thismask_out] = tilts_out[thismask_out]

        # If we tweaked the slits update the tilts_dict
        if self.flatpar['tweak_slits']:
            self.tilts_dict['tilts'] = final_tilts

        if show:
            # Global skysub is the first step in a new extraction so clear the channels here
            self.show(slits=True, wcs_match = True)

        # If illumination flat fielding is turned off, set the illumflat to be None.
        if not self.flatpar['illumflatten']:
            msgs.warn('No illumination flat will be applied to your data (illumflatten=False).')
            self.msillumflat = None

        # Return
        return self.mspixelflat, self.msillumflat

    def show(self, slits=True, wcs_match=True):
        """
        Show all of the flat field products

        Args:
            slits (bool, optional):
            wcs_match (bool, optional):

        """
        viewer, ch = ginga.show_image(self.mspixelflat, chname='pixeflat', cuts=(0.9, 1.1),
                                      wcs_match=wcs_match, clear=True)
        viewer, ch = ginga.show_image(self.msillumflat, chname='illumflat', cuts=(0.9, 1.1),
                                      wcs_match=wcs_match)
        viewer, ch = ginga.show_image(self.rawflatimg, chname='flat', wcs_match=wcs_match)
        viewer, ch = ginga.show_image(self.flat_model, chname='flat_model', wcs_match=wcs_match)

        if slits:
            if self.tslits_dict is not None:
                slit_ids = [trace_slits.get_slitid(self.rawflatimg.shape,
                                                   self.tslits_dict['slit_left'],
                                                   self.tslits_dict['slit_righ'], ii)[0]
                                for ii in range(self.tslits_dict['slit_left'].shape[1])]
                ginga.show_slits(viewer, ch, self.tslits_dict['slit_left'],
                                 self.tslits_dict['slit_righ'], slit_ids)


#    # TODO JFH: This load_master_illumflat code is a bit of a kludge.
#    # Usually one reads in masters with the master and load_master
#    # method, but there are technically two master files for the flats,
#    # i.e. a pixelflat and illumination flat. Perhaps the better way to
#    # deal with this would be to package them into one output file and
#    # just change the load_master and save_master methods to deal with
#    # the possible existence of an illumination flat
#    def load_master_illumflat(self):
#        """
#        Load the slit illumination profile from a saved Master file
#
#        Returns:
#            ndarray: Image from disk
#
#        """
#        ms_name = masterframe.master_name('illumflat', self.master_key, self.mdir)
#        msframe, head = self.load_master(ms_name)
#        if msframe is None:
#            msgs.warn("No Master frame found of type {:s}: {:s}".format('illumflat', ms_name))
#
#        return msframe, head

    def save(self, outfile=None, overwrite=True):
        """
        Save the flat-field master data.

        Args:
            outfile (:obj:`str`, optional):
                Name for the output file.  Defaults to
                :attr:`file_path`.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file.
        """
        super(FlatField, self).save([self.rawflatimg, self.mspixelflat, self.msillumflat],
                                    ['RAWFLAT', 'PIXELFLAT', 'ILLUMFLAT'], outfile=outfile,
                                    overwrite=overwrite, raw_files=self.files, steps=self.steps)

    # TODO: it would be better to have this instantiate the full class
    # as a classmethod.
    def load(self, ifile=None, return_header=False):
        """
        Load the flat-field data from a save master frame.

        Args:
            ifile (:obj:`str`, optional):
                Name of the master frame file.  Defaults to
                :attr:`file_path`.
            return_header (:obj:`bool`, optional):
                Return the header

        Returns:
            tuple: Returns three `numpy.ndarray`_ objects with the raw
            flat-field image, the normalized pixel flat, and the
            illumination flat.  Also returns the primary header, if
            requested.
        """
        return super(FlatField, self).load(['RAWFLAT', 'PIXELFLAT', 'ILLUMFLAT'], ifile=ifile,
                                           return_header=return_header)

