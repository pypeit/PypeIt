"""
Implements the flat-field class.

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import inspect
import numpy as np

from scipy import interpolate, ndimage

from matplotlib import pyplot as plt

from IPython import embed

from pypeit import msgs
from pypeit import utils
from pypeit import ginga
from pypeit import masterframe
from pypeit import bspline

from pypeit.par import pypeitpar
from pypeit.images import calibrationimage
from pypeit.images import pypeitimage
from pypeit.core import flat
from pypeit.core import save
from pypeit.core import load
from pypeit.core import pixels
from pypeit.core import procimg
from pypeit.core import tracewave
from pypeit.core import basis


class FlatField(calibrationimage.CalibrationImage, masterframe.MasterFrame):
    """
    Builds pixel-level flat-field and the illumination flat-field.

    For the primary methods, see :func:`run`.

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
        slits (:class:`pypeit.edgetrace.SlitTraceSet`):
            The current slit traces.
        tilts_dict (:obj:`dict`, optional):
            The current wavelength tilt traces; see
            :class:`pypeit.wavetilts.WaveTilts`.
        reuse_masters (:obj:`bool`, optional):
            Reuse already created master files from disk.

    Attributes:
        rawflatimg (PypeItImage):
        mspixelflat (ndarray):
            Normalized flat
        msillumflat (ndarray):
            Illumination flat
    """

    # Frame type is a class attribute
    frametype = 'pixelflat'
    master_type = 'Flat'

    @classmethod
    def from_master_file(cls, master_file, par=None):
        """
        Instantiate the class from a master file

        Args:
            master_file (str):
            par (:class:`pypeit.par.pypeitpar.PypeItPar`, optional):
                Full par set

        Returns:
            :class:`pypeit.flatfield.FlatField`:
                With the flat images loaded up

        """
        # Spectrograph
        spectrograph, extras = masterframe.items_from_master_file(master_file)
        head0 = extras[0]
        # Par
        if par is None:
            par = spectrograph.default_pypeit_par()
        # Master info
        master_dir = head0['MSTRDIR']
        master_key = head0['MSTRKEY']
        # Instantiate
        slf = cls(spectrograph, par['calibrations']['pixelflatframe'], master_dir=master_dir,
                  master_key=master_key, reuse_masters=True)
        # Load
        rawflatimg, slf.mspixelflat, slf.msillumflat = slf.load(ifile=master_file)
        # Convert rawflatimg to a PypeItImage
        slf.rawflatimg = pypeitimage.PypeItImage(rawflatimg)
        # Return
        return slf

    def __init__(self, spectrograph, par, files=None, det=1, master_key=None,
                 master_dir=None, reuse_masters=False, flatpar=None, msbias=None, msbpm=None,
                 slits=None, tilts_dict=None):

        # Image processing parameters
        if not isinstance(par, pypeitpar.FrameGroupPar):
            msgs.error('Must provide a FrameGroupPar instance as the parameters for FlatField.')
        self.par = par

        # Instantiate the base classes
        #   - Basic processing of the raw images
        calibrationimage.CalibrationImage.__init__(self, spectrograph, det, self.par['process'],
                                                   files=files)
        #   - Construction and interface as a master frame
        masterframe.MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
                                         master_key=master_key, reuse_masters=reuse_masters)

        # FieldFlattening parameters
        self.flatpar = pypeitpar.FlatFieldPar() if flatpar is None else flatpar

        # Input master frame data
        self.msbias = msbias
        self.msbpm = msbpm
        self.slits = slits
        self.tilts_dict = tilts_dict

        # Attributes unique to this Object
        self.rawflatimg = None      # Un-normalized pixel flat as a PypeItImage
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
        Return the number of slits.  Pulled directly from :attr:`slits`, if it exists.
        """
        return 0 if self.slits is None else self.slits.nslits

    def build_pixflat(self, trim=True, force=False):
        """
        Process the flat flat images.

        Processing steps are the result of
        :func:`pypeit.core.procimg.init_process_steps`, ``trim``
        (based on the input argument), ``apply_gain``, and ``orient``.
        Currently, cosmic-ray rejection (``cr_reject``) is not done.

        Args:
            trim (:obj:`bool`, optional):
                Trim the image down to just the data section.
            force (:obj:`bool`, optional):
                Force the flat to be reconstructed if it already exists

        Returns:
            pypeitimage.PypeItImage:  The image with the unnormalized pixel-flat data.
        """
        if self.rawflatimg is None or force:
            # Process steps
            self.process_steps = procimg.init_process_steps(self.msbias, self.par['process'])
            ## JFH We never need untrimmed images. Why is this even an option?
            if trim:
                self.process_steps += ['trim']
            self.process_steps += ['apply_gain']
            self.process_steps += ['orient']
            # Turning this on leads to substantial edge-tracing problems when last tested
            #     JXP November 22, 2019
            #if self.par['cr_reject']:
            #    self.process_steps += ['crmask']
            self.steps.append(inspect.stack()[0][3])
            # Do it
            self.rawflatimg = super(FlatField, self).build_image(bias=self.msbias, bpm=self.msbpm,
                                                                 ignore_saturation=True)
        return self.rawflatimg

#    # TODO Need to add functionality to use a different frame for the ilumination flat, e.g. a sky flat
#    def run(self, debug=False, show=False):
#        """
#        Generate normalized pixel and illumination flats
#
#        Code flow::
#           1.  Generate the pixelflat image (if necessary)
#           2.  Prepare b-spline knot spacing
#           3.  Loop on slits/orders
#               a. Calculate the slit profile
#               b. Normalize
#               c. Save
#
#        Args:
#            debug (:obj:`bool`, optional):
#                Run in debug mode.
#            show (:obj:`bool`, optional):
#                Show the results in the ginga viewer.
#
#        Returns:
#            `numpy.ndarray`_: Two arrays are returned, the normalized
#            pixel flat data and the slit illumination correction data.
#        """
#        # Build the pixel flat (as needed)
#        self.build_pixflat()
#
#        # Prep tck (sets self.ntckx, self.ntcky)
#        #self._prep_tck()
#
#        # Setup
#        self.mspixelflat = np.ones_like(self.rawflatimg.image)
#        self.msillumflat = np.ones_like(self.rawflatimg.image)
#        self.flat_model = np.zeros_like(self.rawflatimg.image)
#        self.tslits_dict = self.slits.to_tslits_dict()
#        self.slitmask = pixels.tslits2mask(self.tslits_dict)
##        self.slitmask = self.slits.slit_img()
#
#        final_tilts = np.zeros_like(self.rawflatimg.image)
#
#        # If we are tweaking slits allocate the new aray to hold tweaked slit boundaries
#        if self.flatpar['tweak_slits']:
##            self.slits.init_tweaks()
#            self.tslits_dict['slit_left_tweak'] = np.zeros_like(self.tslits_dict['slit_left'])
#            self.tslits_dict['slit_righ_tweak'] = np.zeros_like(self.tslits_dict['slit_righ'])
#
#        inmask = np.ones_like(self.rawflatimg.image, dtype=bool) if self.msbpm is None \
#                    else np.invert(self.msbpm)
#        nonlinear_counts = self.spectrograph.nonlinear_counts(det=self.det)
#
#        # Loop on slits
#        for slit in range(self.nslits):
#            # Is this a good slit??
#            if self.slits.mask[slit]:
#                msgs.info('Skipping bad slit: {}'.format(slit))
#                continue
#
#            msgs.info('Computing flat field image for slit: {:d}/{:d}'.format(slit,self.nslits-1))
#
#            # Fit flats for a single slit
#            this_tilts_dict = {'tilts':self.tilts_dict['tilts'],
#                               'coeffs':self.tilts_dict['coeffs'][:,:,slit].copy(),
#                               'slitcen':self.tilts_dict['slitcen'][:,slit].copy(),
#                               'func2d':self.tilts_dict['func2d']}
#
#            pixelflat, illumflat, flat_model, tilts_out, thismask_out, slit_left_out, \
#                    slit_righ_out \
#                            = flat.fit_flat(self.rawflatimg.image, this_tilts_dict, self.tslits_dict,
#                                           slit, inmask=inmask, nonlinear_counts=nonlinear_counts,
#                                           spec_samp_fine=self.flatpar['spec_samp_fine'],
#                                           spec_samp_coarse=self.flatpar['spec_samp_coarse'],
#                                           spat_samp=self.flatpar['spat_samp'],
#                                           tweak_slits=self.flatpar['tweak_slits'],
#                                           tweak_slits_thresh=self.flatpar['tweak_slits_thresh'],
#                                           tweak_slits_maxfrac=self.flatpar['tweak_slits_maxfrac'],
#                                           debug=debug)
#
#            self.mspixelflat[thismask_out] = pixelflat[thismask_out]
#            self.msillumflat[thismask_out] = illumflat[thismask_out]
#            self.flat_model[thismask_out] = flat_model[thismask_out]
#
#            # Did we tweak slit boundaries? If so, update the tslits_dict and the tilts_dict
#            if self.flatpar['tweak_slits']:
#                # TODO: Why do we need slit_left, slit_left_orig, and
#                # slit_left_tweak? Shouldn't we only need two of these?
#                
##                self.slits.left[:,slit] = slit_left_out
##                self.slits.left_tweak[:,slit] = slit_left_out
##                self.slits.right[:,slit] = slit_righ_out
##                self.slits.right_tweak[:,slit] = slit_righ_out
#                self.tslits_dict['slit_left'][:,slit] = slit_left_out
#                self.tslits_dict['slit_righ'][:,slit] = slit_righ_out
#                self.tslits_dict['slit_left_tweak'][:,slit] = slit_left_out
#                self.tslits_dict['slit_righ_tweak'][:,slit] = slit_righ_out
#                final_tilts[thismask_out] = tilts_out[thismask_out]
#
#        # If we tweaked the slits update the tilts_dict
#        if self.flatpar['tweak_slits']:
#            self.tilts_dict['tilts'] = final_tilts
#
#        if show:
#            # Global skysub is the first step in a new extraction so clear the channels here
#            self.show(slits=True, wcs_match = True)
#
#        # If illumination flat fielding is turned off, set the illumflat to be None.
#        if not self.flatpar['illumflatten']:
#            msgs.warn('No illumination flat will be applied to your data (illumflatten=False).')
#            self.msillumflat = None
#
#        # Return
#        return self.mspixelflat, self.msillumflat

    # TODO: Need to add functionality to use a different frame for the
    # ilumination flat, e.g. a sky flat
    def run(self, debug=False, show=False):
        """
        Generate normalized pixel and illumination flats.

        This is a simple wrapper for the main flat-field methods:
        
            - Flat-field images are processed using
              :func:`build_pixflat`.
            - Full 2D model, illumination flat, and pixel flat images
              are constructed by :func:`fit`.
            - The results can be shown in a ginga window using
              :func:`show`.

        The method is a simple wrapper for :func:`build_pixflat`,
        :func:`fit`, and :func:`show`.

        Args:
            debug (:obj:`bool`, optional):
                Run in debug mode.
            show (:obj:`bool`, optional):
                Show the results in the ginga viewer.

        Returns:
            `numpy.ndarray`_: Two arrays are returned, the normalized
            pixel flat data and the slit illumination correction data.
        """
        # Build the pixel flat (as needed)
        self.build_pixflat()

        # Fit it
        # NOTE: Tilts do not change and self.slits is updated
        # internally.
        self.fit(debug=debug)

        if show:
            # Global skysub is the first step in a new extraction so clear the channels here
            self.show(slits=True, wcs_match = True)

#        # If illumination flat fielding is turned off, set the
#        # illumflat to be None.
#        # TODO: Why are we removing this? The illumflat is constructed
#        # by the fit algorithm regardless. It seems like a waste to get
#        # rid of it. Can we just use this parameter to decide if the
#        # data should be illumination flattened at the relevant place
#        # in the codee instead basing it on whether or not msillumflat
#        # is None?
#        if not self.flatpar['illumflatten']:
#            msgs.warn('No illumination flat will be applied to your data (illumflatten=False).')
#            self.msillumflat = None

        # Return
        return self.mspixelflat, self.msillumflat

    def show(self, slits=True, wcs_match=True):
        """
        Show all of the flat field products

        Args:
            slits (bool, optional):
            wcs_match (bool, optional):

        """
        # TODO: Add an option that shows the relevant stuff in a
        # matplotlib window.
        viewer, ch = ginga.show_image(self.mspixelflat, chname='pixeflat', cuts=(0.9, 1.1),
                                      wcs_match=wcs_match, clear=True)
        viewer, ch = ginga.show_image(self.msillumflat, chname='illumflat', cuts=(0.9, 1.1),
                                      wcs_match=wcs_match)
        viewer, ch = ginga.show_image(self.rawflatimg, chname='flat', wcs_match=wcs_match)
        viewer, ch = ginga.show_image(self.flat_model, chname='flat_model', wcs_match=wcs_match)

        if slits and self.slits is not None:
            ginga.show_slits(viewer, ch, self.slits.left, self.slits.right, self.slits.id)

    def save(self, outfile=None, overwrite=True):
        """
        Save the flat-field master data to a FITS file

        Extensions are:
            RAWFLAT
            PIXELFLAT
            ILLUMFLAT
            MODEL

        Args:
            outfile (:obj:`str`, optional):
                Name for the output file.  Defaults to
                :attr:`file_path`.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file.
        """
        _outfile = self.master_file_path if outfile is None else outfile
        # Check if it exists
        if os.path.exists(_outfile) and not overwrite:
            msgs.warn('Master file exists: {0}'.format(_outfile) + msgs.newline()
                      + 'Set overwrite=True to overwrite it.')
            return

        # Setup the items
        hdr = self.build_master_header(steps=self.steps, raw_files=self.file_list)
        data = [self.rawflatimg.image, self.mspixelflat, self.msillumflat, self.flat_model]
        extnames = ['RAWFLAT', 'PIXELFLAT', 'ILLUMFLAT', 'MODEL']

        # Save to a multi-extension FITS
        save.write_fits(hdr, data, _outfile, extnames=extnames)
        msgs.info('Master frame written to {0}'.format(_outfile))

    def load(self, ifile=None):
        """
        Load the flat-field data from a save master frame.

        Args:
            ifile (:obj:`str`, optional):
                Name of the master frame file.  Defaults to
                :attr:`file_path`.

        Returns:
            tuple: Returns three `numpy.ndarray`_ objects with the raw
            flat-field image, the normalized pixel flat, and the
            illumination flat.
        """
        # Check on whether to reuse and whether the file exists
        master_file = self.chk_load_master(ifile)
        if master_file is None:
            return None, None, None
        # Load
        ext = ['RAWFLAT', 'PIXELFLAT', 'ILLUMFLAT', 'MODEL']
        self.rawflatimg, self.mspixelflat, self.msillumflat, self.flat_model, head0 \
                = load.load_multiext_fits(master_file, ext)
        # Return
        return self.rawflatimg, self.mspixelflat, self.msillumflat

    def fit(self, debug=False):
        """
        Construct a model of the flat-field image.

        For this method to work, :attr:`rawflatimg` must have been
        previously constructed; see :func:`build_pixflat`.

        The method loops through all slits provided by the :attr:`slits`
        object, except those that have been masked (i.e., slits with
        ``self.slits.mask == True`` are skipped).  For each slit:

            - Collapse the flat-field data spatially using the
              wavelength coordinates provided by the fit to the arc-line
              traces (:class:`pypeit.wavetilts.WaveTilts`), and fit the
              result with a bspline.  This provides the
              spatially-averaged spectral response of the instrument.
              The data used in the fit is trimmed toward the clit
              spatial center via the ``slit_trim`` parameter in
              :attr:`flatpar`.
            - Use the bspline fit to construct and normalize out the
              spectral response.
            - Collapse the normalized flat-field data spatially using a
              coordinate system defined by the left slit edge.  The data
              included in the spatial (illumination) profile calculation
              is expanded beyond the nominal slit edges using the
              ``slit_pad`` parameter in :attr:`flatpar`.  The raw,
              collapsed data is then median filtered (see ``spat_samp``
              in :attr:`flatpar`) and Gaussian filtered; see
              :func:`pypeit.core.flat.illum_filter`.  This creates an
              empirical, highly smoothed representation of the
              illumination profile that is fit with a bspline.  The
              construction of the empirical illumination profile (i.e.,
              before the bspline fitting) can be done iteratively, where
              each iteration sigma-clips outliers; see the
              ``illum_iter`` and ``illum_rej`` parameters in
              :attr:`flatpar` and
              :func:`pypeit.core.flat.construct_illum_profile`.
            - Use the bspline fit to construct the 2D illumination image
              (:attr:`msillumflat`) and normalize out the spatial
              response.
            - If requested, the 1D illumination profile is used to
              "tweak" the slit edges by offsetting them to a threshold
              of the illumination peak to either side of the slit center
              (see ``tweak_slits_thresh`` in :attr:`flatpar`), up to a
              maximum allowed shift from the existing slit edge (see
              ``tweak_slits_maxfrac`` in :attr:`flatpar`).  See
              :func:`pypeit.core.tweak_slit_edges`.
            - Fit the residuals of the flat-field data that has been
              independently normalized for its spectral and spatial
              response with a 2D bspline-polynomial fit.  The order of
              the polynomial has been optimized via experimentation; it
              can be changed but you should use extreme caution when
              doing so (see ``twod_fit_npoly``).  The multiplication of
              the 2D spectral response, 2D spatial response, and joint
              2D fit to the high-order residuals define the final flat
              model (:attr:`flat_model`).
            - Finally, the pixel-to-pixel response of the instrument is
              defined as the ratio of the raw flat data to the
              best-fitting flat-field model (:attr:`mspixelflat`)

        This method is the primary method that builds the
        :class:`FlatField` instance, constructing :attr:`mspixelflat`,
        :attr:`msillumflat`, and :attr:`flat_model`.  All of these
        attributes are altered internally.  If the slit edges are to be
        tweaked using the 1D illumination profile (``tweak_slits`` in
        :attr:`flatpar`), the tweaked slit edge arrays in the internal
        :class:`pypeit.edgetrace.SlitTraceSet` object, :attr:`slits`,
        are also altered.

        Used parameters from :attr:`flatpar`
        (:class:`pypeit.par.pypeitpar.FlatFieldPar`) are
        ``spec_samp_fine``, ``spec_samp_coarse``, ``spat_samp``,
        ``tweak_slits``, ``tweak_slits_thresh``,
        ``tweak_slits_maxfrac``, ``rej_sticky``, ``slit_trim``,
        ``slit_pad``, ``illum_iter``, ``illum_rej``, and
        ``twod_fit_npoly``.

        Revision History
        ----------------

            - 11-Mar-2005  First version written by Scott Burles.
            - 2005-2018    Improved by J. F. Hennawi and J. X. Prochaska
            - 3-Sep-2018 Ported to python by J. F. Hennawi and significantly improved

        Args:
            debug (:obj:`bool`, optional):
                Show plots useful for debugging. This will block
                further execution of the code until the plot windows
                are closed.

        """
        # TODO: break up this function!  Can it be partitioned into a
        # series of "core" methods?

        # TODO: The difference between run() and fit() is pretty minimal
        # if we just built rawflatimg here...

        # Flat must have been constructed
        if self.rawflatimg is None:
            raise ValueError('The flat-field image has not been built: run build_pixflat first.')

        # Set parameters (for convenience; get rid of this and just use
        # the parameter values directly?)
        spec_samp_fine = self.flatpar['spec_samp_fine']
        spec_samp_coarse = self.flatpar['spec_samp_coarse']
        spat_samp = self.flatpar['spat_samp']
        tweak_slits = self.flatpar['tweak_slits']
        tweak_slits_thresh = self.flatpar['tweak_slits_thresh']
        tweak_slits_maxfrac = self.flatpar['tweak_slits_maxfrac']
        # If sticky, points rejected at each stage (spec, spat, 2d) are
        # propagated to the next stage
        sticky = self.flatpar['rej_sticky']
        trim = self.flatpar['slit_trim']
        pad = self.flatpar['slit_pad']
        # Iteratively construct the illumination profile by rejecting outliers
        illum_iter = self.flatpar['illum_iter']
        illum_rej = self.flatpar['illum_rej']
        npoly = self.flatpar['twod_fit_npoly']

        # Setup images
        nspec, nspat = self.rawflatimg.image.shape
        # TODO: The above should be the same as self.slits.nspec, self.slits.nspat
        rawflat = self.rawflatimg.image
        gpm = np.ones_like(rawflat, dtype=bool) if self.msbpm is None else np.invert(self.msbpm)

        # Flat-field modeling is done in the log of the counts
        flat_log = np.log(np.fmax(rawflat, 1.0))
        gpm_log = (rawflat > 1.0) & gpm
        # set errors to just be 0.5 in the log
        ivar_log = gpm_log.astype(float)/0.5**2

        # Other setup
        nonlinear_counts = self.spectrograph.nonlinear_counts(det=self.det)

        median_slit_width = np.median(self.slits.right - self.slits.left, axis=0)

        if tweak_slits:
            # NOTE: This copies the input slit edges to a set that can
            # be tweaked. Because these are copied immediately before
            # the calls to slit_img below, all the slit images use the
            # original slit edges, not any pre-existing tweaked ones.
            self.slits.init_tweaked()

        # TODO: This needs to include a padding check

        # Construct three versions of the slit ID image
        #   - an image that uses the padding defined by self.slits
        slitid_img = self.slits.slit_img()
        #   - an image that uses the padding defined by self.flatpar.
        #     This was always 5 pixels in the previous version.
        padded_slitid_img = self.slits.slit_img(pad=pad)
        #   - and an image that trims the width of the slit using the
        #     parameter in self.flatpar. This was always 3 pixels in
        #     the previous version.
        # TODO: Fix this for when trim is a tuple
        trimmed_slitid_img = self.slits.slit_img(pad=-trim)

        # Prep for results
        self.mspixelflat = np.ones_like(rawflat)
        self.msillumflat = np.ones_like(rawflat)
        self.flat_model = np.zeros_like(rawflat)

        # Model each slit independently
        for slit in range(self.slits.nslits):
            # Is this a good slit??
            if self.slits.mask[slit]:
                msgs.info('Skipping bad slit: {}'.format(slit))
                continue

            msgs.info('Modeling the flat-field response for slit: {0}/{1}'.format(
                        slit+1, self.slits.nslits))

            # Find the pixels on the slit
            onslit = slitid_img == slit

            # Check for saturation of the flat. If there are not enough
            # pixels do not attempt a fit, and continue to the next
            # slit.  TODO: set the threshold to a parameter?
            good_frac = np.sum(onslit & (rawflat < nonlinear_counts))/np.sum(onslit)
            if good_frac < 0.5:
                # TODO: Used slit ID in these print statments instead of slit index
                msgs.warn('Only {:4.2f}'.format(100*good_frac)
                          + '% of the pixels on this slit are not saturated.' + msgs.newline()
                          + 'Consider raising nonlinear_counts={:5.3f}'.format(nonlinear_counts) +
                          msgs.newline() + 'Not attempting to flat field slit {:d}'.format(slit))
                continue

            # Demand at least 10 pixels per row (on average) per degree
            # of the polynomial.
            # NOTE: This is not used until the 2D fit. Defined here to
            # be close to the definition of ``onslit``.
            if npoly is None:
                # Approximate number of pixels sampling each spatial pixel
                # for this (original) slit.
                npercol = np.fmax(np.floor(np.sum(onslit)/nspec),1.0)
                #npoly  = max(1, int(np.ceil(npercol/10.)))
                # TODO: Added an upper limit.  Needs testing.
                npoly  = int(np.clip(np.ceil(npercol/10.), 1, 10))

            # TODO: Always calculate the optimized `npoly` and warn the
            # user if npoly is provided but higher than the nominal
            # calculation?

            # Create an image with the spatial coordinates relative to the left edge of this slit
            spat_coo = self.slits.spatial_coordinate_image(slitids=slit, full=True)

            # Find pixels on the padded and trimmed slit coordinates
            onslit_padded = padded_slitid_img == slit
            onslit_trimmed = trimmed_slitid_img == slit

            # ----------------------------------------------------------
            # Collapse the slit spatially and fit the spectral function

            # Create the tilts image for this slit
            # TODO: Is the copy needed?
            tilts = tracewave.fit2tilts(rawflat.shape, self.tilts_dict['coeffs'][:,:,slit].copy(),
                                        self.tilts_dict['func2d'])
            # Convert the tilt image to an image with the spectral pixel index
            spec_coo = tilts * (nspec-1)

            # Only include the trimmed set of pixels in the flat-field
            # fit along the spectral direction.
            spec_gpm = onslit_trimmed & gpm_log # & (rawflat < nonlinear_counts)
            spec_nfit = np.sum(spec_gpm)
            spec_ntot = np.sum(onslit)
            msgs.info('Spectral fit of flatfield for {0}/{1} '.format(spec_nfit, spec_ntot)
                      + ' pixels in the slit.')
            # Set this to a parameter?
            if spec_nfit/spec_ntot < 0.5:
                # TODO: Shouldn't this raise an exception or continue to the next slit instead?
                msgs.warn('Spectral fit includes only {:.1f}'.format(100*spec_nfit/spec_ntot)
                          + '% of the pixels on this slit.' + msgs.newline()
                          + '          Either the slit has many bad pixels or the number of '
                            'trimmed pixels is too large.')

            # Sort the pixels by their spectral coordinate.
            # TODO: Include ivar and sorted gpm in outputs?
            spec_gpm, spec_srt, spec_coo_data, spec_flat_data \
                    = flat.sorted_flat_data(flat_log, spec_coo, gpm=spec_gpm)
            # NOTE: By default np.argsort sorts the data over the last
            # axis. Just to avoid the possibility (however unlikely) of
            # spec_coo[spec_gpm] returning an array, all the arrays are
            # explicitly flattened.
            spec_ivar_data = ivar_log[spec_gpm].ravel()[spec_srt]
            spec_gpm_data = gpm_log[spec_gpm].ravel()[spec_srt]

            # Rejection threshold for spectral fit in log(image)
            # TODO: Make this a parameter?
            logrej = 0.5

            # Fit the spectral direction of the blaze.
            # TODO: Figure out how to deal with the fits going crazy at
            # the edges of the chip in spec direction
            # TODO: Can we add defaults to bspline_profile so that we
            # don't have to instantiate invvar and profile_basis
            spec_bspl, spec_gpm_fit, spec_flat_fit, _, exit_status \
                    = utils.bspline_profile(spec_coo_data, spec_flat_data, spec_ivar_data,
                                            np.ones_like(spec_coo_data), ingpm=spec_gpm_data,
                                            nord=4, upper=logrej, lower=logrej,
                                            kwargs_bspline={'bkspace': spec_samp_fine},
                                            kwargs_reject={'groupbadpix': True, 'maxrej': 5})

# TODO: Used for testing bspline
#            np.savez_compressed('slit{0}_spec.npz'.format(slit+1),
#                                spec_coo_data=spec_coo_data,
#                                spec_flat_data=spec_flat_data,
#                                spec_ivar_data=spec_ivar_data,
#                                spec_gpm_data=spec_gpm_data)

            if exit_status > 1:
                # TODO: Should this fault?
                msgs.warn('Problem in bspline fit!')

            # Debugging/checking spectral fit
            if debug:
                utils.bspline_profile_qa(spec_coo_data, spec_flat_data, spec_bspl, spec_gpm_fit,
                                         spec_flat_fit, xlabel='Spectral Pixel',
                                         ylabel='log(flat counts)',
                                         title='Spectral Fit for slit={:d}'.format(slit))

            if sticky:
                # Add rejected pixels to gpm
                gpm[spec_gpm] = (spec_gpm_fit & spec_gpm_data)[np.argsort(spec_srt)]

            # Construct the model of the flat-field spectral shape
            spec_model = np.ones_like(rawflat)
            spec_model[onslit_padded] = np.exp(spec_bspl.value(spec_coo[onslit_padded])[0])
            # ----------------------------------------------------------

            # ----------------------------------------------------------
            # To fit the spatial response, first normalize out the
            # spectral response, and then collapse the slit spectrally.

            # Normalize out the spectral shape of the flat
            norm_spec = np.ones_like(rawflat)
            norm_spec[onslit_padded] = rawflat[onslit_padded] \
                                            / np.fmax(spec_model[onslit_padded],1.0)

            # Find pixels fot fit in the spatial direction:
            #   - Fit pixels in the padded slit that haven't been masked
            #     by the BPM
            spat_gpm = onslit_padded & gpm #& (rawflat < nonlinear_counts)
            #   - Fit pixels with non-zero flux and less than 70% above
            #     the average spectral profile.
            spat_gpm &= (norm_spec > 0.0) & (norm_spec < 1.7)
            #   - Determine maximum counts in median filtered flat
            #     spectrum model.
            spec_interp = interpolate.interp1d(spec_coo_data, spec_flat_fit, kind='linear',
                                               assume_sorted=True, bounds_error=False,
                                               fill_value=-np.inf)
            spec_sm = utils.fast_running_median(np.exp(spec_interp(np.arange(nspec))),
                                                np.fmax(np.ceil(0.10*nspec).astype(int),10))
            #   - Only fit pixels with at least values > 10% of this maximum and no less than 1.
            spat_gpm &= (spec_model > 0.1*np.amax(spec_sm)) & (spec_model > 1.0)

            # Report
            spat_nfit = np.sum(spat_gpm)
            spat_ntot = np.sum(onslit_padded)
            msgs.info('Spatial fit of flatfield for {0}/{1} '.format(spat_nfit, spat_ntot)
                      + ' pixels in the slit.')
            if spat_nfit/spat_ntot < 0.5:
                # TODO: Shouldn't this raise an exception or continue to the next slit instead?
                msgs.warn('Spatial fit includes only {:.1f}'.format(100*spat_nfit/spat_ntot)
                          + '% of the pixels on this slit.' + msgs.newline()
                          + '          Either the slit has many bad pixels, the model of the '
                          'spectral shape is poor, or the illumination profile is very irregular.')

            # Construct the empirical illumination profile
            _spat_gpm, spat_srt, spat_coo_data, spat_flat_data_raw, spat_flat_data \
                    = flat.construct_illum_profile(norm_spec, spat_coo, median_slit_width[slit],
                                                   spat_gpm=spat_gpm, spat_samp=spat_samp,
                                                   illum_iter=illum_iter, illum_rej=illum_rej,
                                                   debug=debug)

            if sticky:
                # Add rejected pixels to gpm
                gpm[spat_gpm] &= (spat_gpm & _spat_gpm)[spat_gpm]

            # Make sure that the normalized and filtered flat is finite!
            if np.any(np.invert(np.isfinite(spat_flat_data))):
                msgs.error('Inifinities in slit illumination function computation!')

            # Determine the breakpoint spacing from the sampling of the
            # spatial coordinates. Use breakpoints at a spacing of a
            # 1/10th of a pixel, but do not allow a bsp smaller than
            # the typical sampling. Use the bspline class to determine
            # the breakpoints:
#            spat_bspl = pydl.bspline(spat_coo_data, nord=4,
#                                     bkspace=np.fmax(1.0/median_slit_width[slit]/10.0,
#                                                     1.2*np.median(np.diff(spat_coo_data))))
            spat_bspl = bspline.bspline(spat_coo_data, nord=4,
                                        bkspace=np.fmax(1.0/median_slit_width[slit]/10.0,
                                                     1.2*np.median(np.diff(spat_coo_data))))
            # TODO: Can we add defaults to bspline_profile so that we
            # don't have to instantiate invvar and profile_basis
            spat_bspl, spat_gpm_fit, spat_flat_fit, _, exit_status \
                    = utils.bspline_profile(spat_coo_data, spat_flat_data,
                                            np.ones_like(spat_flat_data),
                                            np.ones_like(spat_flat_data), nord=4, upper=5.0,
                                            lower=5.0, fullbkpt=spat_bspl.breakpoints)

# TODO: Used for testing bspline
#            np.savez_compressed('slit{0}_spat.npz'.format(slit+1),
#                                spat_coo_data=spat_coo_data,
#                                spat_flat_data=spat_flat_data,
#                                median_slit_width=median_slit_width[slit])

            # Compute the model illumination profile
            self.msillumflat[onslit_padded] = spat_bspl.value(spat_coo[onslit_padded])[0]
            # ----------------------------------------------------------

            # ----------------------------------------------------------
            # Tweak the slit edges based on the empirical slit
            # illumination profiles, if requested
            if tweak_slits:
                # TODO: Should the tweak be based on the bspline fit?
                left_thresh, left_shift, self.slits.left_tweak[:,slit], right_thresh, \
                    right_shift, self.slits.right_tweak[:,slit] \
                        = flat.tweak_slit_edges(self.slits.left[:,slit], self.slits.right[:,slit],
                                                spat_coo_data, spat_flat_data,
                                                thresh=tweak_slits_thresh,
                                                maxfrac=tweak_slits_maxfrac)
                # TODO: Because the padding doesn't consider adjacent
                # slits, calling slit_img for individual slits can be
                # different from the result when you construct the
                # image for all slits. Fix this...

                # Update the onslit mask
                _slitid_img = self.slits.slit_img(slitids=slit)
                onslit = _slitid_img == slit
                spat_coo = self.slits.spatial_coordinate_image(slitids=slit,
                                                               slitid_img=_slitid_img)
                # Note that nothing changes with the tilts, since these were
                # already extrapolated across the whole image.
            else:
                _slitid_img = slitid_img

            # Add an approximate pixel axis at the top
            if debug:
                # TODO: Move this into a qa plot that gets saved
                ax = utils.bspline_profile_qa(spat_coo_data, spat_flat_data, spat_bspl,
                                              spat_gpm_fit, spat_flat_fit, show=False)
                ax.scatter(spat_coo_data, spat_flat_data_raw, marker='.', s=1, zorder=0, color='k',
                           label='raw data')
                # Force the center of the slit to be at the center of the plot for the hline
                ax.set_xlim(-0.1,1.1)
                ax.axvline(0.0, color='lightgreen', linestyle=':', linewidth=2.0,
                           label='original left edge', zorder=8)
                ax.axvline(1.0, color='red', linestyle=':', linewidth=2.0,
                           label='original right edge', zorder=8)
                if tweak_slits and left_shift > 0:
                    label = 'threshold = {:5.2f}'.format(tweak_slits_thresh) \
                                + ' % of max of left illumprofile'
                    ax.axhline(left_thresh, xmax=0.5, color='lightgreen', linewidth=3.0,
                               label=label, zorder=10)
                    ax.axvline(left_shift, color='lightgreen', linestyle='--', linewidth=3.0,
                               label='tweaked left edge', zorder=11)
                if tweak_slits and right_shift > 0:
                    label = 'threshold = {:5.2f}'.format(tweak_slits_thresh) \
                                + ' % of max of right illumprofile'
                    ax.axhline(right_thresh, xmin=0.5, color='red', linewidth=3.0, label=label,
                               zorder=10)
                    ax.axvline(1-right_shift, color='red', linestyle='--', linewidth=3.0,
                               label='tweaked right edge', zorder=20)
                ax.legend()
                ax.set_xlabel('Normalized Slit Position')
                ax.set_ylabel('Normflat Spatial Profile')
                ax.set_title('Illumination Function Fit for slit={:d}'.format(slit))
                plt.show()

            # ----------------------------------------------------------

            # ----------------------------------------------------------
            # Fit the 2D residuals of the 1D spectral and spatial fits.

            msgs.info('Performing 2D illumination + scattered light flat field fit')

            # Construct the spectrally and spatially normalized flat
            norm_spec_spat = np.ones_like(rawflat)
            norm_spec_spat[onslit_padded] = rawflat[onslit_padded] \
                                                / np.fmax(spec_model[onslit_padded], 1.0) \
                                                / np.fmax(self.msillumflat[onslit_padded], 0.01)

            # Sort the pixels by their spectral coordinate. The mask
            # uses the nominal padding defined by the slits object.
            twod_gpm, twod_srt, twod_spec_coo_data, twod_flat_data \
                    = flat.sorted_flat_data(norm_spec_spat, spec_coo, gpm=onslit)
            # Also apply the sorting to the spatial coordinates
            twod_spat_coo_data = spat_coo[twod_gpm].ravel()[twod_srt]
            # TODO: Reset back to origin gpm if sticky is true?
            twod_gpm_data = gpm[twod_gpm].ravel()[twod_srt]
            # Only fit data with less than 30% variations
            # TODO: Make 30% a parameter?
            twod_gpm_data &= np.absolute(twod_flat_data - 1) < 0.3
            # Here we ignore the formal photon counting errors and
            # simply assume that a typical error per pixel. This guess
            # is somewhat aribtrary. We then set the rejection
            # threshold with sigrej_twod
            # TODO: Make twod_sig and twod_sigrej parameters?
            twod_sig = 0.01
            twod_ivar_data = twod_gpm_data.astype(float)/(twod_sig**2)
            twod_sigrej = 4.0

#            poly_basis = pydl.fpoly(2.0*twod_spat_coo_data - 1.0, npoly).T
            poly_basis = basis.fpoly(2.0*twod_spat_coo_data - 1.0, npoly)

            # Perform the full 2d fit
            twod_bspl, twod_gpm_fit, twod_flat_fit, _ , exit_status \
                    = utils.bspline_profile(twod_spec_coo_data, twod_flat_data, twod_ivar_data,
                                            poly_basis, ingpm=twod_gpm_data, nord=4,
                                            upper=twod_sigrej, lower=twod_sigrej,
                                            kwargs_bspline={'bkspace': spec_samp_coarse},
                                            kwargs_reject={'groupbadpix': True, 'maxrej': 10})

# TODO: Used for testing bspline
#            np.savez_compressed('slit{0}_twod.npz'.format(slit+1),
#                                twod_spat_coo_data=twod_spat_coo_data,
#                                twod_spec_coo_data=twod_spec_coo_data,
#                                twod_flat_data=twod_flat_data,
#                                twod_ivar_data=twod_ivar_data,
#                                poly_basis=poly_basis,
#                                twod_gpm_data=twod_gpm_data)

            if debug:
                # TODO: Make a plot that shows the residuals in the 2D
                # image
                resid = twod_flat_data - twod_flat_fit
                goodpix = twod_gpm_fit & twod_gpm_data
                badpix = np.invert(twod_gpm_fit) & twod_gpm_data

                plt.clf()
                ax = plt.gca()
                ax.plot(twod_spec_coo_data[goodpix], resid[goodpix], color='k', marker='o',
                        markersize=0.2, mfc='k', fillstyle='full', linestyle='None',
                        label='good points')
                ax.plot(twod_spec_coo_data[badpix], resid[badpix], color='red', marker='+',
                        markersize=0.5, mfc='red', fillstyle='full', linestyle='None',
                        label='masked')
                ax.axhline(twod_sigrej*twod_sig, color='lawngreen', linestyle='--',
                           label='rejection thresholds', zorder=10, linewidth=2.0)
                ax.axhline(-twod_sigrej*twod_sig, color='lawngreen', linestyle='--', zorder=10,
                           linewidth=2.0)
#                ax.set_ylim(-0.05, 0.05)
                ax.legend()
                ax.set_xlabel('Spectral Pixel')
                ax.set_ylabel('Residuals from pixelflat 2-d fit')
                ax.set_title('Spectral Residuals for slit={:d}'.format(slit))
                plt.show()

                plt.clf()
                ax = plt.gca()
                ax.plot(twod_spat_coo_data[goodpix], resid[goodpix], color='k', marker='o',
                        markersize=0.2, mfc='k', fillstyle='full', linestyle='None',
                        label='good points')
                ax.plot(twod_spat_coo_data[badpix], resid[badpix], color='red', marker='+',
                        markersize=0.5, mfc='red', fillstyle='full', linestyle='None',
                        label='masked')
                ax.axhline(twod_sigrej*twod_sig, color='lawngreen', linestyle='--',
                           label='rejection thresholds', zorder=10, linewidth=2.0)
                ax.axhline(-twod_sigrej*twod_sig, color='lawngreen', linestyle='--', zorder=10,
                           linewidth=2.0)
#                ax.set_ylim((-0.05, 0.05))
#                ax.set_xlim(-0.02, 1.02)
                ax.legend()
                ax.set_xlabel('Normalized Slit Position')
                ax.set_ylabel('Residuals from pixelflat 2-d fit')
                ax.set_title('Spatial Residuals for slit={:d}'.format(slit))
                plt.show()

            # Save the 2D residual model
            twod_model = np.ones_like(rawflat)
            twod_model[twod_gpm] = twod_flat_fit[np.argsort(twod_srt)]

            # Construct the full flat-field model
            # TODO: Why is the 0.05 here for the illumflat compared to the 0.01 above?
            self.flat_model[onslit] = twod_model[onslit] \
                                        * np.fmax(self.msillumflat[onslit],0.05) \
                                        * np.fmax(spec_model[onslit],1.0)

            # Construct the pixel flat
            self.mspixelflat[onslit] = rawflat[onslit]/self.flat_model[onslit]
            # TODO: Add some code here to treat the edges and places where fits
            # go bad?

        # Set the pixelflat to 1.0 wherever the flat was nonlinear
        self.mspixelflat[rawflat >= nonlinear_counts] = 1.0
        # Do not apply pixelflat field corrections that are greater than
        # 100% to avoid creating edge effects, etc.
        # TODO: Should we do the same for the illumflat?
        self.mspixelflat = np.clip(self.mspixelflat, 0.5, 2.0)



