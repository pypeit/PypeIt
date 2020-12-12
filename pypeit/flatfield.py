"""
Implements the flat-field class.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import copy
import inspect
import numpy as np

from scipy import interpolate

from matplotlib import pyplot as plt

from IPython import embed

from pypeit import msgs
from pypeit import utils
from pypeit import bspline

from pypeit import datamodel
from pypeit import masterframe
from pypeit.display import display
from pypeit.core import flat
from pypeit.core import tracewave
from pypeit.core import basis
from pypeit.core import fitting
from pypeit.core import coadd
from pypeit import slittrace


class FlatImages(datamodel.DataContainer):
    """
    Simple DataContainer for the output from Flatfield

    All of the items in the datamodel are required for instantiation,
      although they can be None (but shouldn't be)

    """
    minimum_version = '1.1.0'
    version = '1.1.0'

    # I/O
    output_to_disk = None  # This writes all items that are not None
    hdu_prefix = None

    # Master fun
    master_type = 'Flat'
    master_file_format = 'fits'

    datamodel = {'pixelflat_raw': dict(otype=np.ndarray, atype=np.floating,
                                       descr='Processed, combined pixel flats'),
                 'pixelflat_norm': dict(otype=np.ndarray, atype=np.floating,
                                        descr='Normalized pixel flat'),
                 'pixelflat_model': dict(otype=np.ndarray, atype=np.floating, descr='Model flat'),
                 'pixelflat_spat_bsplines': dict(otype=np.ndarray, atype=bspline.bspline,
                                                 descr='B-spline models for pixel flat'),
                 'pixelflat_bpm': dict(otype=np.ndarray, atype=np.integer,
                                       descr='Mirrors SlitTraceSet mask for Flat-specific flags'),
                 'pixelflat_spec_illum': dict(otype=np.ndarray, atype=np.floating,
                                              descr='Relative spectral illumination'),
                 'illumflat_raw': dict(otype=np.ndarray, atype=np.floating,
                                       descr='Processed, combined illum flats'),
                 'illumflat_spat_bsplines': dict(otype=np.ndarray, atype=bspline.bspline,
                                                 descr='B-spline models for illum flat'),
                 'illumflat_bpm': dict(otype=np.ndarray, atype=np.integer,
                                       descr='Mirrors SlitTraceSet mask for Flat-specific flags'),
                 'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 'spat_id': dict(otype=np.ndarray, atype=np.integer, descr='Slit spat_id')}

    def __init__(self, pixelflat_raw=None, pixelflat_norm=None, pixelflat_bpm=None,
                 pixelflat_model=None, pixelflat_spat_bsplines=None, pixelflat_spec_illum=None,
                 illumflat_raw=None, illumflat_spat_bsplines=None, illumflat_bpm=None,
                 PYP_SPEC=None, spat_id=None):
        # Parse
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=d)

    def _init_internals(self):
        self.filename = None
        # Master stuff
        self.master_key = None
        self.master_dir = None

    def _validate(self):
        #
        if self.pixelflat_spat_bsplines is not None and len(self.pixelflat_spat_bsplines) > 0:
            if len(self.spat_id) != len(self.pixelflat_spat_bsplines):
                msgs.error("Pixelflat Bsplines are out of sync with the slit IDs")
        if self.illumflat_spat_bsplines is not None and len(self.illumflat_spat_bsplines) > 0:
            if len(self.spat_id) != len(self.illumflat_spat_bsplines):
                msgs.error("Illumflat Bsplines are out of sync with the slit IDs")

    def is_synced(self, slits):
        """
        Confirm the slits in WaveTilts are aligned to that in SlitTraceSet

        Barfs if not

        Args:
            slits (:class:`pypeit.slittrace.SlitTraceSet`):

        """
        if not np.array_equal(self.spat_id, slits.spat_id):
            msgs.error("Your flat solutions are out of sync with your slits.  Remove Masters and start from scratch")

    def _bundle(self):
        """
        Over-write default _bundle() method to write one
        HDU per image.  Any extras are in the HDU header of
        the primary image.

        Returns:
            :obj:`list`: A list of dictionaries, each list element is
            written to its own fits extension. See the description
            above.
        """
        d = []
        for key in self.keys():
            # Skip None
            if self[key] is None:
                continue
            if self.datamodel[key]['otype'] == np.ndarray and 'bsplines' not in key:
                d += [{key: self[key]}]
            elif 'bsplines' in key:
                flattype = 'pixelflat' if 'pixelflat' in key else 'illumflat'
                d += [{'{0}_spat_id-{1}_bspline'.format(flattype, self.spat_id[ss]): self[key][ss]}
                            for ss in range(len(self[key]))]
            else:
                if len(d) > 0:
                    d[0][key] = self[key]
                else:
                    d += [{key: self[key]}]
        # Return
        return d

    @classmethod
    def _parse(cls, hdu, ext=None, transpose_table_arrays=False, hdu_prefix=None):

        # Grab everything but the bsplines. The bsplines are not parsed
        # because the tailored extension names do not match any of the
        # datamodel keys.
        d, version_passed, type_passed, parsed_hdus = super(FlatImages, cls)._parse(hdu)

        # Find bsplines, if they exist
        nspat = len(d['spat_id'])
        hdunames = [h.name for h in hdu]
        for flattype in ['pixelflat', 'illumflat']:
            ext = ['{0}_SPAT_ID-{1}_BSPLINE'.format(flattype.upper(), d['spat_id'][i])
                   for i in range(nspat)]
            indx = np.isin(ext, hdunames)
            if np.any(indx) and not np.all(indx):
                msgs.error('Expected {0} {1} bspline extensions, but only found {1}.'.format(
                           nspat, flattype, np.sum(indx)))
            if np.all(indx):
                key = '{0}_spat_bsplines'.format(flattype)
                try:
                    d[key] = np.array([bspline.bspline.from_hdu(hdu[k]) for k in ext])
                except Exception as e:
                    msgs.warn('Error in bspline extension read:\n {0}: {1}'.format(
                                e.__class__.__name__, str(e)))
                    # Assume this is because the type failed
                    type_passed = False
                else:
                    version_passed &= np.all([d[key][i].version == bspline.bspline.version 
                                              for i in range(nspat)])
                    parsed_hdus += ext

        return d, version_passed, type_passed, parsed_hdus

    def shape(self):
        if self.pixelflat_raw is not None:
            return self.pixelflat_raw.shape
        elif self.illumflat_raw is not None:
            return self.illumflat_raw.shape
        else:
            msgs.error("Shape of FlatImages could not be determined")

    def get_procflat(self, frametype='pixel'):
        if frametype == 'illum':
            return self.illumflat_raw
        else:
            return self.pixelflat_raw

    def get_bpmflats(self, frametype='pixel'):
        # Check if both BPMs are none
        if self.pixelflat_bpm is None and self.illumflat_bpm is None:
            msgs.warn("FlatImages contains no BPM - trying to generate one")
            return np.zeros(self.shape(), dtype=np.int)
        # Now return the requested case, checking for None
        if frametype == 'illum':
            if self.illumflat_bpm is not None:
                return self.illumflat_bpm
            else:
                msgs.warn("illumflat has no BPM - using the pixelflat BPM")
                return self.pixelflat_bpm
        else:
            if self.pixelflat_bpm is not None:
                return self.pixelflat_bpm
            else:
                msgs.warn("pixelflat has no BPM - using the illumflat BPM")
                return self.illumflat_bpm

    def get_spat_bsplines(self, frametype='illum'):
        """
        Grab a list of bspline fits

        Args:
            frametype (str):

        Returns:
            list:

        """
        # Check if both spat bsplines are none
        if self.pixelflat_spat_bsplines is None and self.illumflat_spat_bsplines is None:
            msgs.error("FlatImages contains no spatial bspline fit")
        # Now return the requested case, checking for None
        if frametype == 'illum':
            if self.illumflat_spat_bsplines is not None:
                return self.illumflat_spat_bsplines
            else:
                msgs.warn("illumflat has no spatial bspline fit - using the pixelflat")
                return self.pixelflat_spat_bsplines
        else:
            if self.pixelflat_spat_bsplines is not None:
                return self.pixelflat_spat_bsplines
            else:
                msgs.warn("pixelflat has no spatial bspline fit - using the illumflat")
                return self.illumflat_spat_bsplines

    def get_pixelflat(self):
        return self.pixelflat_norm

    def get_spec_illum(self):
        return self.pixelflat_spec_illum

    def get_flat_model(self):
        return self.pixelflat_model

    def fit2illumflat(self, slits, frametype='illum', initial=False, flexure_shift=None):
        """

        Args:
            slits (:class:`pypeit.slittrace.SlitTraceSet`):
            frametype (str):
                Should the pixel or illum flat spatial profile be generated? (options include: 'illum' or 'pixel').
                The default is to use 'illum' unless frametype='pixel'.
            initial (bool, optional):
            flexure_shift (float, optional):

        Returns:

        """
        illumflat = np.ones(self.shape())
        # Load spatial bsplines
        spat_bsplines = self.get_spat_bsplines(frametype=frametype)

        # Loop
        for slit_idx in range(slits.nslits):
            # Skip masked
            if slits.mask[slit_idx] != 0:
                continue
            # Skip those without a bspline
            # DO it
            _slitid_img = slits.slit_img(slitidx=slit_idx, initial=initial, flexure=flexure_shift)
            onslit = _slitid_img == slits.spat_id[slit_idx]
            spat_coo = slits.spatial_coordinate_image(slitidx=slit_idx,
                                                      initial=initial,
                                                      slitid_img=_slitid_img,
                                                      flexure_shift=flexure_shift)
            illumflat[onslit] = spat_bsplines[slit_idx].value(spat_coo[onslit])[0]
        # TODO -- Update the internal one?  Or remove it altogether??
        return illumflat

    def show(self, frametype='all', slits=None, wcs_match=True):
        """
        Simple wrapper to show_flats()

        Args:
            frametype (str):
                Which flats should be displayed? Must be one of 'illum', 'pixel', or 'all' (default)
            slits:
            wcs_match:

        Returns:

        """
        illumflat_pixel, illumflat_illum = None, None
        # Try to grab the slits
        if slits is None:
            # Warning: This parses the filename, not the Header!
            master_key, master_dir = masterframe.grab_key_mdir(self.filename, from_filename=True)
            try:
                slit_masterframe_name = masterframe.construct_file_name(slittrace.SlitTraceSet, master_key,
                                                                        master_dir=master_dir)
                slits = slittrace.SlitTraceSet.from_file(slit_masterframe_name)
            except:
                msgs.warn('Could not load slits to show with flat-field images. Did you provide the master info??')
        if slits is not None:
            slits.mask_flats(self)
            illumflat_pixel = self.fit2illumflat(slits, frametype='pixel')
            if self.illumflat_spat_bsplines is not None:
                illumflat_illum = self.fit2illumflat(slits, frametype='illum')
        # Decide which frames should be displayed
        if frametype == 'pixel':
            image_list = zip([self.pixelflat_norm, illumflat_pixel, self.pixelflat_raw,
                              self.pixelflat_model, self.pixelflat_spec_illum],
                             ['pixelflat_norm', 'pixelflat_spat_illum', 'pixelflat_raw',
                              'pixelflat_model', 'pixelflat_spec_illum'],
                             [(0.9, 1.1), (0.9, 1.1), None, None,
                              (0.8, 1.2)])
        elif frametype == 'illum':
            image_list = zip([illumflat_illum, self.illumflat_raw],
                             ['illumflat_spat_illum', 'illumflat_raw'],
                             [(0.9, 1.1), None])
        else:
            # Show everything that's available (anything that is None will not be displayed)
            image_list = zip([self.pixelflat_norm, illumflat_pixel, self.pixelflat_raw,
                              self.pixelflat_model, self.pixelflat_spec_illum,
                              illumflat_illum, self.illumflat_raw],
                             ['pixelflat_norm', 'pixelflat_spat_illum', 'pixelflat_raw',
                              'pixelflat_model', 'pixelflat_spec_illum',
                              'illumflat_spat_illum', 'illumflat_raw'],
                             [(0.9, 1.1), (0.9, 1.1), None, None,
                              (0.8, 1.2), (0.9, 1.1), None])
        # Display frames
        show_flats(image_list, wcs_match=wcs_match, slits=slits)



class FlatField(object):
    """
    Builds pixel-level flat-field and the illumination flat-field.

    For the primary methods, see :func:`run`.

    Args:
        rawflatimg (:class:`pypeit.images.pypeitimage.PypeItImage`):
            Processed, combined set of pixelflat images
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The `Spectrograph` instance that sets the instrument used to
            take the observations.  See usage by
            :class:`pypeit.processimages.ProcessImages` base class.
        flatpar (:class:`pypeit.par.pypeitpar.FlatFieldPar`):
            User-level parameters for constructing the flat-field
            corrections.  If None, the default parameters are used.
        slits (:class:`pypeit.slittrace.SlitTraceSet`):
            The current slit traces.
        wavetilts (:class:`pypeit.wavetilts.WaveTilts`):
            The current wavelength tilt traces; see
        spat_illum_only (boolean):
            Only perform the spatial illumination calculation, and ignore
            the 2D bspline fit. This should only be set to true if you
            want the spatial illumination profile only. If you want to
            simultaneously generate a pixel flat and a spatial
            illumination profile from the same input, this should be
            False (which is the default).

    Attributes:
        rawflatimg (:class:`pypeit.images.pypeitimage.PypeItImage`):
        mspixelflat (`numpy.ndarray`_):
            Normalized flat
        msillumflat (`numpy.ndarray`_):
            Illumination flat
        flat_model (`numpy.ndarray`_):
            Model of the flat
        list_of_spat_bsplines (list):
        spec_illum (`numpy.ndarray`_):
            Image of the relative spectral illumination for a multislit spectrograph

    """

    # Frame type is a class attribute
    frametype = 'pixelflat'
    master_type = 'Flat'


    def __init__(self, rawflatimg, spectrograph, flatpar, slits, wavetilts, wv_calib, spat_illum_only=False):

        # Defaults
        self.spectrograph = spectrograph
        # FieldFlattening parameters
        self.flatpar = flatpar

        # Input data
        self.slits = slits
        self.wavetilts = wavetilts
        self.wv_calib = wv_calib

        # Worth a check
        self.wavetilts.is_synced(self.slits)

        # Attributes unique to this Object
        self.rawflatimg = rawflatimg      # Un-normalized pixel flat as a PypeItImage
        self.mspixelflat = None     # Normalized pixel flat
        self.msillumflat = None     # Illumination flat
        self.flat_model = None      # Model flat
        self.list_of_spat_bsplines = None
        self.spat_illum_only = spat_illum_only
        self.spec_illum = None      # Relative spectral illumination image

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
            :class:`FlatImages`:
        """
        # Build the pixel flat (as needed)
        #self.build_pixflat()

        # Fit it
        # NOTE: Tilts do not change and self.slits is updated internally.
        self.fit(spat_illum_only=self.spat_illum_only, debug=debug)

        if show:
            # Global skysub is the first step in a new extraction so clear the channels here
            self.show(wcs_match=True)

        # Build the mask
        bpmflats = self.build_mask()

        # Return
        if self.spat_illum_only:
            # Illumination correction only
            return FlatImages(illumflat_raw=self.rawflatimg.image,
                              illumflat_spat_bsplines=np.asarray(self.list_of_spat_bsplines),
                              illumflat_bpm=bpmflats, PYP_SPEC=self.spectrograph.name,
                              spat_id=self.slits.spat_id)
        else:
            # Pixel and illumination correction only
            return FlatImages(pixelflat_raw=self.rawflatimg.image,
                              pixelflat_norm=self.mspixelflat,
                              pixelflat_model=self.flat_model,
                              pixelflat_spat_bsplines=np.asarray(self.list_of_spat_bsplines),
                              pixelflat_bpm=bpmflats, pixelflat_spec_illum=self.spec_illum,
                              PYP_SPEC=self.spectrograph.name, spat_id=self.slits.spat_id)

    def build_mask(self):
        """
        Generate bad pixel mask.

        Returns:
            :obj:`numpy.ndarray` : bad pixel mask
        """
        bpmflats = np.zeros_like(self.slits.mask, dtype=self.slits.bitmask.minimum_dtype())
        for flag in ['SKIPFLATCALIB', 'BADFLATCALIB']:
            bpm = self.slits.bitmask.flagged(self.slits.mask, flag)
            if np.any(bpm):
                bpmflats[bpm] = self.slits.bitmask.turn_on(bpmflats[bpm], flag)
        return bpmflats

    def show(self, wcs_match=True):
        """
        Show all of the flat field products in ginga.

        Args:
            wcs_match (:obj:`bool`, optional):
                Match the WCS of the flat-field images
        """
        # Prepare the images to show, their names and their cuts
        image_list = zip([self.mspixelflat, self.msillumflat, self.rawflatimg.image, self.flat_model],
                         ['pixelflat', 'spat_illum', 'raw', 'model', 'spec_illum'],
                         [(0.9, 1.1), (0.9, 1.1), None, None, (0.8, 1.2)])
        show_flats(image_list, wcs_match=wcs_match, slits=self.slits)

    def fit(self, spat_illum_only=False, debug=False):
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
              The data used in the fit is trimmed toward the slit
              spatial center via the ``slit_trim`` parameter in
              :attr:`flatpar`.
            - Use the bspline fit to construct and normalize out the
              spectral response.
            - Collapse the normalized flat-field data spatially using a
              coordinate system defined by the left slit edge.  The data
              included in the spatial (illumination) profile calculation
              is expanded beyond the nominal slit edges using the
              ``slit_illum_pad`` parameter in :attr:`flatpar`.  The raw,
              collapsed data is then median filtered (see ``spat_samp``
              in :attr:`flatpar`) and Gaussian filtered; see
              :func:`pypeit.core.flat.illum_filter`.  This creates an
              empirical, highly smoothed representation of the
              illumination profile that is fit with a bspline using
              the :func:`spatial_fit` method.  The
              construction of the empirical illumination profile (i.e.,
              before the bspline fitting) can be done iteratively, where
              each iteration sigma-clips outliers; see the
              ``illum_iter`` and ``illum_rej`` parameters in
              :attr:`flatpar` and
              :func:`pypeit.core.flat.construct_illum_profile`.
            - If requested, the 1D illumination profile is used to
              "tweak" the slit edges by offsetting them to a threshold
              of the illumination peak to either side of the slit center
              (see ``tweak_slits_thresh`` in :attr:`flatpar`), up to a
              maximum allowed shift from the existing slit edge (see
              ``tweak_slits_maxfrac`` in :attr:`flatpar`).  See
              :func:`pypeit.core.tweak_slit_edges`.  If tweaked, the
              :func:`spatial_fit` is repeated to place it on the tweaked
              slits reference frame.
            - Use the bspline fit to construct the 2D illumination image
              (:attr:`msillumflat`) and normalize out the spatial
              response.
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
        ``slit_illum_pad``, ``illum_iter``, ``illum_rej``, and
        ``twod_fit_npoly``, ``saturated_slits``.

        **Revision History**:

            - 11-Mar-2005  First version written by Scott Burles.
            - 2005-2018    Improved by J. F. Hennawi and J. X. Prochaska
            - 3-Sep-2018 Ported to python by J. F. Hennawi and significantly improved

        Args:
            spat_illum_only (:obj:`bool`, optional):
                If true, only the spatial illumination profile will be calculated.
                The 2D bspline fit will not be performed. This is primarily used
                to build an illumflat.
            debug (:obj:`bool`, optional):
                Show plots useful for debugging. This will block
                further execution of the code until the plot windows
                are closed.

        """
        # TODO: break up this function!  Can it be partitioned into a series of "core" methods?
        # TODO: JFH I wrote all this code and will have to maintain it and I don't want to see it broken up.
        # TODO: JXP This definitely needs breaking up..

        # Initialise with a series of bad splines (for when slits go wrong)
        if self.list_of_spat_bsplines is None:
            self.list_of_spat_bsplines = [bspline.bspline(None) for all in self.slits.spat_id]

        # Set parameters (for convenience;
        spec_samp_fine = self.flatpar['spec_samp_fine']
        spec_samp_coarse = self.flatpar['spec_samp_coarse']
        tweak_slits = self.flatpar['tweak_slits']
        tweak_slits_thresh = self.flatpar['tweak_slits_thresh']
        tweak_slits_maxfrac = self.flatpar['tweak_slits_maxfrac']
        # If sticky, points rejected at each stage (spec, spat, 2d) are
        # propagated to the next stage
        sticky = self.flatpar['rej_sticky']
        trim = self.flatpar['slit_trim']
        pad = self.flatpar['slit_illum_pad']
        # Iteratively construct the illumination profile by rejecting outliers
        npoly = self.flatpar['twod_fit_npoly']
        saturated_slits = self.flatpar['saturated_slits']

        # Setup images
        nspec, nspat = self.rawflatimg.image.shape
        rawflat = self.rawflatimg.image
        # Good pixel mask
        gpm = np.ones_like(rawflat, dtype=bool) if self.rawflatimg.bpm is None else (
                1-self.rawflatimg.bpm).astype(bool)

        # Flat-field modeling is done in the log of the counts
        flat_log = np.log(np.fmax(rawflat, 1.0))
        gpm_log = (rawflat > 1.0) & gpm
        # set errors to just be 0.5 in the log
        ivar_log = gpm_log.astype(float)/0.5**2

        # Other setup
        nonlinear_counts = self.spectrograph.nonlinear_counts(self.rawflatimg.detector)

        # TODO -- JFH -- CONFIRM THIS SHOULD BE ON INIT
        # It does need to be *all* of the slits
        median_slit_widths = np.median(self.slits.right_init - self.slits.left_init, axis=0)

        if tweak_slits:
            # NOTE: This copies the input slit edges to a set that can be tweaked.
            self.slits.init_tweaked()

        # TODO: This needs to include a padding check
        # Construct three versions of the slit ID image, all of unmasked slits!
        #   - an image that uses the padding defined by self.slits
        slitid_img_init = self.slits.slit_img(initial=True)
        #   - an image that uses the extra padding defined by
        #     self.flatpar. This was always 5 pixels in the previous
        #     version.
        padded_slitid_img = self.slits.slit_img(initial=True, pad=pad)
        #   - and an image that trims the width of the slit using the
        #     parameter in self.flatpar. This was always 3 pixels in
        #     the previous version.
        # TODO: Fix this for when trim is a tuple
        trimmed_slitid_img = self.slits.slit_img(pad=-trim, initial=True)

        # Prep for results
        self.mspixelflat = np.ones_like(rawflat)
        self.msillumflat = np.ones_like(rawflat)
        self.flat_model = np.zeros_like(rawflat)

        # Allocate work arrays only once
        spec_model = np.ones_like(rawflat)
        norm_spec = np.ones_like(rawflat)
        norm_spec_spat = np.ones_like(rawflat)
        twod_model = np.ones_like(rawflat)
        twod_gpm_out = np.ones_like(rawflat, dtype=np.bool)

        # #################################################
        # Model each slit independently
        for slit_idx, slit_spat in enumerate(self.slits.spat_id):
            # Is this a good slit??
            if self.slits.mask[slit_idx] != 0:
                msgs.info('Skipping bad slit: {}'.format(slit_spat))
                continue

            msgs.info('Modeling the flat-field response for slit spat_id={}: {}/{}'.format(
                        slit_spat, slit_idx+1, self.slits.nslits))

            # Find the pixels on the initial slit
            onslit_init = slitid_img_init == slit_spat

            # Check for saturation of the flat. If there are not enough
            # pixels do not attempt a fit, and continue to the next
            # slit.
            # TODO: set the threshold to a parameter?
            good_frac = np.sum(onslit_init & (rawflat < nonlinear_counts))/np.sum(onslit_init)
            if good_frac < 0.5:
                common_message = 'To change the behavior, use the \'saturated_slits\' parameter ' \
                                 'in the \'flatfield\' parameter group; see here:\n\n' \
                                 'https://pypeit.readthedocs.io/en/latest/pypeit_par.html \n\n' \
                                 'You could also choose to use a different flat-field image ' \
                                 'for this calibration group.'
                if saturated_slits == 'crash':
                    msgs.error('Only {:4.2f}'.format(100*good_frac)
                               + '% of the pixels on slit {0} are not saturated.  '.format(slit_spat)
                               + 'Selected behavior was to crash if this occurred.  '
                               + common_message)
                elif saturated_slits == 'mask':
                    self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(self.slits.mask[slit_idx], 'BADFLATCALIB')
                    msgs.warn('Only {:4.2f}'.format(100*good_frac)
                                                + '% of the pixels on slit {0} are not saturated.  '.format(slit_spat)
                              + 'Selected behavior was to mask this slit and continue with the '
                              + 'remainder of the reduction, meaning no science data will be '
                              + 'extracted from this slit.  ' + common_message)
                elif saturated_slits == 'continue':
                    self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(self.slits.mask[slit_idx], 'SKIPFLATCALIB')
                    msgs.warn('Only {:4.2f}'.format(100*good_frac)
                              + '% of the pixels on slit {0} are not saturated.  '.format(slit_spat)
                              + 'Selected behavior was to simply continue, meaning no '
                              + 'field-flatting correction will be applied to this slit but '
                              + 'pypeit will attempt to extract any objects found on this slit.  '
                              + common_message)
                else:
                    # Should never get here
                    raise NotImplementedError('Unknown behavior for saturated slits: {0}'.format(
                                              saturated_slits))
                continue

            # Demand at least 10 pixels per row (on average) per degree
            # of the polynomial.
            # NOTE: This is not used until the 2D fit. Defined here to
            # be close to the definition of ``onslit``.
            if npoly is None:
                # Approximate number of pixels sampling each spatial pixel
                # for this (original) slit.
                npercol = np.fmax(np.floor(np.sum(onslit_init)/nspec),1.0)
                npoly  = np.clip(7, 1, int(np.ceil(npercol/10.)))
            
            # TODO: Always calculate the optimized `npoly` and warn the
            #  user if npoly is provided but higher than the nominal
            #  calculation?

            # Create an image with the spatial coordinates relative to the left edge of this slit
            spat_coo_init = self.slits.spatial_coordinate_image(slitidx=slit_idx, full=True, initial=True)

            # Find pixels on the padded and trimmed slit coordinates
            onslit_padded = padded_slitid_img == slit_spat
            onslit_trimmed = trimmed_slitid_img == slit_spat

            # ----------------------------------------------------------
            # Collapse the slit spatially and fit the spectral function
            # TODO: Put this stuff in a self.spectral_fit method?

            # Create the tilts image for this slit
            # TODO -- JFH Confirm the sign of this shift is correct!
            _flexure = 0. if self.wavetilts.spat_flexure is None else self.wavetilts.spat_flexure
            tilts = tracewave.fit2tilts(rawflat.shape, self.wavetilts['coeffs'][:,:,slit_idx],
                                        self.wavetilts['func2d'], spat_shift=-1*_flexure)
            # Convert the tilt image to an image with the spectral pixel index
            spec_coo = tilts * (nspec-1)

            # Only include the trimmed set of pixels in the flat-field
            # fit along the spectral direction.
            spec_gpm = onslit_trimmed & gpm_log  # & (rawflat < nonlinear_counts)
            spec_nfit = np.sum(spec_gpm)
            spec_ntot = np.sum(onslit_init)
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
            #  the edges of the chip in spec direction
            # TODO: Can we add defaults to bspline_profile so that we
            #  don't have to instantiate invvar and profile_basis
            try:
                spec_bspl, spec_gpm_fit, spec_flat_fit, _, exit_status \
                    = fitting.bspline_profile(spec_coo_data, spec_flat_data, spec_ivar_data,
                                            np.ones_like(spec_coo_data), ingpm=spec_gpm_data,
                                            nord=4, upper=logrej, lower=logrej,
                                            kwargs_bspline={'bkspace': spec_samp_fine},
                                            kwargs_reject={'groupbadpix': True, 'maxrej': 5})
            except:
                embed(header='808 of flatfield')

            if exit_status > 1:
                # TODO -- MAKE A FUNCTION
                msgs.warn('Flat-field spectral response bspline fit failed!  Not flat-fielding '
                          'slit {0} and continuing!'.format(slit_spat))
                self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(self.slits.mask[slit_idx], 'BADFLATCALIB')
                continue

            # Debugging/checking spectral fit
            if debug:
                fitting.bspline_qa(spec_coo_data, spec_flat_data, spec_bspl, spec_gpm_fit,
                                 spec_flat_fit, xlabel='Spectral Pixel', ylabel='log(flat counts)',
                                 title='Spectral Fit for slit={:d}'.format(slit_spat))

            if sticky:
                # Add rejected pixels to gpm
                gpm[spec_gpm] = (spec_gpm_fit & spec_gpm_data)[np.argsort(spec_srt)]

            # Construct the model of the flat-field spectral shape
            # including padding on either side of the slit.
            spec_model[...] = 1.
            spec_model[onslit_padded] = np.exp(spec_bspl.value(spec_coo[onslit_padded])[0])
            # ----------------------------------------------------------

            # ----------------------------------------------------------
            # To fit the spatial response, first normalize out the
            # spectral response, and then collapse the slit spectrally.

            # Normalize out the spectral shape of the flat
            norm_spec[...] = 1.
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

            # First fit -- With initial slits
            exit_status, spat_coo_data,  spat_flat_data, spat_bspl, spat_gpm_fit, \
                spat_flat_fit, spat_flat_data_raw \
                        = self.spatial_fit(norm_spec, spat_coo_init, median_slit_widths[slit_idx],
                                           spat_gpm, gpm, debug=debug)

            if tweak_slits:
                # TODO: Should the tweak be based on the bspline fit?
                # TODO: Will this break if
                left_thresh, left_shift, self.slits.left_tweak[:,slit_idx], right_thresh, \
                    right_shift, self.slits.right_tweak[:,slit_idx] \
                        = flat.tweak_slit_edges(self.slits.left_init[:,slit_idx],
                                                self.slits.right_init[:,slit_idx],
                                                spat_coo_data, spat_flat_data,
                                                thresh=tweak_slits_thresh,
                                                maxfrac=tweak_slits_maxfrac, debug=debug)
                # TODO: Because the padding doesn't consider adjacent
                #  slits, calling slit_img for individual slits can be
                #  different from the result when you construct the
                #  image for all slits. Fix this...

                # Update the onslit mask
                _slitid_img = self.slits.slit_img(slitidx=slit_idx, initial=False)
                onslit_tweak = _slitid_img == slit_spat
                spat_coo_tweak = self.slits.spatial_coordinate_image(slitidx=slit_idx,
                                                               slitid_img=_slitid_img)

                # Construct the empirical illumination profile
                # TODO This is extremely inefficient, because we only need to re-fit the illumflat, but
                #  spatial_fit does both the reconstruction of the illumination function and the bspline fitting.
                #  Only the b-spline fitting needs be reddone with the new tweaked spatial coordinates, so that would
                #  save a ton of runtime. It is not a trivial change becauase the coords are sorted, etc.
                exit_status, spat_coo_data, spat_flat_data, spat_bspl, spat_gpm_fit, \
                    spat_flat_fit, spat_flat_data_raw = self.spatial_fit(
                    norm_spec, spat_coo_tweak, median_slit_widths[slit_idx], spat_gpm, gpm, debug=False)

                spat_coo_final = spat_coo_tweak
            else:
                _slitid_img = slitid_img_init
                spat_coo_final = spat_coo_init
                onslit_tweak = onslit_init

            # Add an approximate pixel axis at the top
            if debug:
                # TODO: Move this into a qa plot that gets saved
                ax = fitting.bspline_qa(spat_coo_data, spat_flat_data, spat_bspl, spat_gpm_fit,
                                      spat_flat_fit, show=False)
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
                ax.set_title('Illumination Function Fit for slit={:d}'.format(slit_spat))
                plt.show()

            # ----------------------------------------------------------
            # Construct the illumination profile with the tweaked edges
            # of the slit
            if exit_status <= 1:
                # TODO -- JFH -- Check this is ok for flexure!!
                self.msillumflat[onslit_tweak] = spat_bspl.value(spat_coo_final[onslit_tweak])[0]
                self.list_of_spat_bsplines[slit_idx] = spat_bspl
                # No need to proceed further if we just need the illumination profile
                if spat_illum_only:
                    continue
            else:
                # Save the nada
                msgs.warn('Slit illumination profile bspline fit failed!  Spatial profile not '
                          'included in flat-field model for slit {0}!'.format(slit_spat))
                self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(self.slits.mask[slit_idx], 'BADFLATCALIB')
                continue

            # ----------------------------------------------------------
            # Fit the 2D residuals of the 1D spectral and spatial fits.
            msgs.info('Performing 2D illumination + scattered light flat field fit')

            # Construct the spectrally and spatially normalized flat
            norm_spec_spat[...] = 1.
            norm_spec_spat[onslit_tweak] = rawflat[onslit_tweak] / np.fmax(spec_model[onslit_tweak], 1.0) \
                                                    / np.fmax(self.msillumflat[onslit_tweak], 0.01)

            # Sort the pixels by their spectral coordinate. The mask
            # uses the nominal padding defined by the slits object.
            twod_gpm, twod_srt, twod_spec_coo_data, twod_flat_data \
                    = flat.sorted_flat_data(norm_spec_spat, spec_coo, gpm=onslit_tweak)
            # Also apply the sorting to the spatial coordinates
            twod_spat_coo_data = spat_coo_final[twod_gpm].ravel()[twod_srt]
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

            poly_basis = basis.fpoly(2.0*twod_spat_coo_data - 1.0, npoly)

            # Perform the full 2d fit
            twod_bspl, twod_gpm_fit, twod_flat_fit, _, exit_status \
                    = fitting.bspline_profile(twod_spec_coo_data, twod_flat_data, twod_ivar_data,
                                            poly_basis, ingpm=twod_gpm_data, nord=4,
                                            upper=twod_sigrej, lower=twod_sigrej,
                                            kwargs_bspline={'bkspace': spec_samp_coarse},
                                            kwargs_reject={'groupbadpix': True, 'maxrej': 10})
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
                ax.set_title('Spectral Residuals for slit={:d}'.format(slit_spat))
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
                ax.set_title('Spatial Residuals for slit={:d}'.format(slit_spat))
                plt.show()

            # Save the 2D residual model
            twod_model[...] = 1.
            if exit_status > 1:
                msgs.warn('Two-dimensional fit to flat-field data failed!  No higher order '
                          'flat-field corrections included in model of slit {0}!'.format(slit_spat))
            else:
                twod_model[twod_gpm] = twod_flat_fit[np.argsort(twod_srt)]
                twod_gpm_out[twod_gpm] = twod_gpm_fit[np.argsort(twod_srt)]


            # Construct the full flat-field model
            # TODO: Why is the 0.05 here for the illumflat compared to the 0.01 above?
            self.flat_model[onslit_tweak] = twod_model[onslit_tweak] \
                                        * np.fmax(self.msillumflat[onslit_tweak], 0.05) \
                                        * np.fmax(spec_model[onslit_tweak], 1.0)

            # Construct the pixel flat
            #self.mspixelflat[onslit] = rawflat[onslit]/self.flat_model[onslit]
            #self.mspixelflat[onslit_tweak] = 1.
            #trimmed_slitid_img_anew = self.slits.slit_img(pad=-trim, slitidx=slit_idx)
            #onslit_trimmed_anew = trimmed_slitid_img_anew == slit_spat
            self.mspixelflat[onslit_tweak] = rawflat[onslit_tweak]/self.flat_model[onslit_tweak]
            # TODO: Add some code here to treat the edges and places where fits
            #  go bad?

        # No need to continue if we're just doing the spatial illumination
        if spat_illum_only:
            return

        # Set the pixelflat to 1.0 wherever the flat was nonlinear
        self.mspixelflat[rawflat >= nonlinear_counts] = 1.0
        # Set the pixelflat to 1.0 within trim pixels of all the slit edges
        trimmed_slitid_img_new = self.slits.slit_img(pad=-trim, initial=False)
        tweaked_slitid_img = self.slits.slit_img(initial=False)
        self.mspixelflat[(trimmed_slitid_img_new < 0) & (tweaked_slitid_img > 0)] = 1.0

        # Do not apply pixelflat field corrections that are greater than
        # 100% to avoid creating edge effects, etc.
        self.mspixelflat = np.clip(self.mspixelflat, 0.5, 2.0)

        # Finally, using the above products, calculate the relative spectral illumination, if requested
        if self.flatpar['slit_illum_relative']:
            self.spec_illum = self.spectral_illumination(twod_gpm_out, debug=debug)

    def spatial_fit(self, norm_spec, spat_coo, median_slit_width, spat_gpm, gpm, debug=False):
        """
        Perform the spatial fit

        Args:
            norm_spec (`numpy.ndarray`_):
            spat_coo (`numpy.ndarray`_):
                Spatial coordinate array
            median_slit_width (:obj:`float`):
            spat_gpm (`numpy.ndarray`_):
            gpm (`numpy.ndarray`_):
            debug (bool, optional):

        Returns:
            tuple: 7 objects
                 - exit_status (int):
                 - spat_coo_data
                 - spat_flat_data
                 - spat_bspl (:class:`pypeit.bspline.bspline.bspline`): Bspline model of the spatial fit.  Used for illumflat
                 - spat_gpm_fit
                 - spat_flat_fit
                 - spat_flat_data_raw
        """

        # Construct the empirical illumination profile
        _spat_gpm, spat_srt, spat_coo_data, spat_flat_data_raw, spat_flat_data \
            = flat.construct_illum_profile(norm_spec, spat_coo, median_slit_width,
                                           spat_gpm=spat_gpm,
                                           spat_samp=self.flatpar['spat_samp'],
                                           illum_iter=self.flatpar['illum_iter'],
                                           illum_rej=self.flatpar['illum_rej'],
                                           debug=debug)

        if self.flatpar['rej_sticky']:
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
        spat_bspl = bspline.bspline(spat_coo_data, nord=4,
                                    bkspace=np.fmax(1.0 / median_slit_width / 10.0,
                                                    1.2 * np.median(np.diff(spat_coo_data))))
        # TODO: Can we add defaults to bspline_profile so that we
        #  don't have to instantiate invvar and profile_basis
        spat_bspl, spat_gpm_fit, spat_flat_fit, _, exit_status \
            = fitting.bspline_profile(spat_coo_data, spat_flat_data,
                                    np.ones_like(spat_flat_data),
                                    np.ones_like(spat_flat_data), nord=4, upper=5.0,
                                    lower=5.0, fullbkpt=spat_bspl.breakpoints)
        # Return
        return exit_status, spat_coo_data, spat_flat_data, spat_bspl, spat_gpm_fit, \
               spat_flat_fit, spat_flat_data_raw

    def spectral_illumination(self, gpm=None, debug=False):
        """
        Generate a relative scaling image for a slit-based IFU. All
        slits are scaled relative to the zeroth slit. There are three
        stages in this approach:

            1. Get a quick, rough scaling between the orders using a
               low order polynomial

            2. Using this rough scale, perform a joint b-spline fit
               to all slits. This step ensures that a single
               functional form is used in step 3 to fit all slits. It
               also ensures that the model covers the min and max
               wavelength range of all slits.

            3. Calculate the relative scale of each slit, using the
               joint model calculated in step (2).

        Parameters
        ----------
        gpm : `numpy.ndarray`_, None
            Good pixel mask
        debug : bool
            Debug the routine

        Returns
        -------
        scale_model: `numpy.ndarray`_
            An image containing the appropriate scaling
        """
        msgs.info("Deriving spectral illumination profile")
        # Generate a wavelength image
        msgs.info("Generating wavelength image")
        flex = self.wavetilts.spat_flexure
        slitmask = self.slits.slit_img(initial=True, flexure=flex)
        tilts = self.wavetilts.fit2tiltimg(slitmask, flexure=flex)
        #waveimg = wavecalib.build_waveimg(self.spectrograph, tilts, self.slits, self.wv_calib, spat_flexure=flex)
        waveimg = self.wv_calib.build_waveimg(tilts, self.slits, spat_flexure=flex)
        msgs.info('Performing a joint fit to the flat-field response')
        # Grab some parameters
        trim = self.flatpar['slit_trim']
        spec_samp_fine = self.flatpar['spec_samp_coarse']
        rawflat = self.rawflatimg.image.copy() / self.msillumflat.copy()
        # Grab the GPM and the slit images
        if gpm is None:
            gpm = np.ones_like(rawflat, dtype=bool) if self.rawflatimg.bpm is None else (
                    1 - self.rawflatimg.bpm).astype(bool)

        slitid_img_init = self.slits.slit_img(pad=0, initial=True)
        slitid_img_trim = self.slits.slit_img(pad=-trim, initial=True)
        # Find all good slits, and create a mask of pixels to include (True=include)
        wgd = self.slits.spat_id[np.where(self.slits.mask == 0)]
        # Obtain the minimum and maximum wavelength of all slits
        mnmx_wv = np.zeros((self.slits.nslits, 2))
        for slit_idx, slit_spat in enumerate(self.slits.spat_id):
            onslit_init = (slitid_img_init == slit_spat)
            mnmx_wv[slit_idx, 0] = np.min(waveimg[onslit_init])
            mnmx_wv[slit_idx, 1] = np.max(waveimg[onslit_init])
        # Sort by increasing minimum wavelength
        swslt = np.argsort(mnmx_wv[:, 0])

        ### STEP 1
        relscl_model = illum_profile_spectral(rawflat, waveimg, self.slits, model=None, gpmask=gpm, skymask=None,
                                              trim=trim, flexure=flex)

        ### STEP 2
        # Perform a simultaneous fit to all pixels in all slits to get a "global" shape of the flat spectrum.
        # This ensures that the final fit smoothly covers the full wavelength range covered on the detector.
        # Get the pixels containing good slits
        spec_tot = np.isin(slitid_img_init, wgd)  # & (rawflat < nonlinear_counts)
        # Apply the relative scaling
        rawflatscl = rawflat / relscl_model
        # Flat-field modeling is done in the log of the counts
        flat_log = np.log(np.fmax(rawflatscl, 1.0))
        gpm_log = (rawflatscl > 1.0) & gpm
        # set errors to just be 0.5 in the log
        ivar_log = gpm_log.astype(float) / 0.5 ** 2
        # Only include the trimmed set of pixels in the flat-field
        # fit along the spectral direction.
        spec_gpm = np.isin((slitid_img_trim), wgd) & gpm_log  # & (rawflat < nonlinear_counts)
        spec_nfit = np.sum(spec_gpm)
        spec_ntot = np.sum(spec_tot)
        msgs.info('Spectral fit of flatfield for {0}/{1} '.format(spec_nfit, spec_ntot)
                  + ' pixels on all slits.')
        # Sort the pixels by their spectral coordinate.
        # TODO: Include ivar and sorted gpm in outputs?
        spec_gpm, spec_srt, spec_coo_data, spec_flat_data \
            = flat.sorted_flat_data(flat_log, waveimg, gpm=spec_gpm)
        spec_ivar_data = ivar_log[spec_gpm].ravel()[spec_srt]
        spec_gpm_data = gpm_log[spec_gpm].ravel()[spec_srt]

        # Fit the spectral direction of the blaze.
        logrej = 0.5
        spec_bspl, spec_gpm_fit, spec_flat_fit, _, exit_status \
            = fitting.bspline_profile(spec_coo_data, spec_flat_data, spec_ivar_data,
                                    np.ones_like(spec_coo_data), ingpm=spec_gpm_data,
                                    nord=4, upper=logrej, lower=logrej,
                                    kwargs_bspline={'bkspace': spec_samp_fine},
                                    kwargs_reject={'groupbadpix': True, 'maxrej': 5})

        ### STEP 3
        # Redo the scale model, now using the bspline fit
        scale_model = np.ones_like(self.rawflatimg.image)
        for slit_idx in range(0, self.slits.spat_id.size):
            msgs.info("Generating model relative response image for slit {0:d}".format(slit_idx))
            # Only use the overlapping regions of the slits, where the same wavelength range is covered
            onslit = (slitid_img_trim == self.slits.spat_id[swslt[slit_idx]])
            onslit_init = (slitid_img_init == self.slits.spat_id[swslt[slit_idx]])
            onslit_gpm = onslit & gpm
            # Fit a low order polynomial
            minw, maxw = mnmx_wv[slit_idx, 0], mnmx_wv[slit_idx, 1]
            xfit = (waveimg[onslit_gpm] - minw) / (maxw - minw)
            yfit = rawflat[onslit_gpm] / np.exp(spec_bspl.value(waveimg[onslit_gpm])[0])
            srtd = np.argsort(xfit)
            # Rough outlier rejection
            med = np.median(yfit)
            mad = 1.4826*np.median(np.abs(med-yfit))
            inmsk = (yfit-med > -10*mad) & (yfit-med < 10*mad)
            slit_bspl, _, _, _, exit_status \
                = fitting.bspline_profile(xfit[srtd], yfit[srtd], np.ones_like(xfit)/mad**2, np.ones_like(xfit),
                                        nord=4, upper=3, lower=3, ingpm=inmsk[srtd],
                                        kwargs_bspline={'bkspace': spec_samp_fine},
                                        kwargs_reject={'groupbadpix': True, 'maxrej': 5})
            # TODO : Perhaps mask a slit if it fails...
            if exit_status > 1:
                msgs.warn("b-spline fit of relative scale failed for slit {0:d}".format(slit_idx))
            else:
                scale_model[onslit_init] = 1/slit_bspl.value((waveimg[onslit_init] - minw) / (maxw - minw))[0]

        if debug:
            embed()
            pltflat = self.rawflatimg.image.copy() / self.msillumflat.copy()
            censpec = np.round(0.5 * (self.slits.left_init + self.slits.right_init)).astype(np.int)
            for ss in range(self.slits.nslits):
                #plt.plot(waveimg[(np.arange(censpec.shape[0]), censpec[:,ss].flatten())], scale_model[(np.arange(censpec.shape[0]), censpec[:,ss].flatten())])
                plt.plot(waveimg[(np.arange(censpec.shape[0]), censpec[:, ss].flatten())],
                         pltflat[(np.arange(censpec.shape[0]), censpec[:, ss].flatten())] *
                         scale_model[(np.arange(censpec.shape[0]), censpec[:, ss].flatten())])
            plt.show()
            # This code generates the wavy patterns seen in KCWI
            debug_model = np.ones_like(self.rawflatimg.image)
            blaze_model = np.ones_like(self.rawflatimg.image)
            if exit_status > 1:
                msgs.warn("Joint blaze fit failed")
            else:
                blaze_model[...] = 1.
                blaze_model[spec_tot] = np.exp(spec_bspl.value(waveimg[spec_tot])[0])
                # Now take out the relative scaling
                blaze_model /= scale_model
                # Now, we want to use the raw flat image, corrected for spatial illumination and pixel-to-pixel variations
                corr_model = self.msillumflat
                corr_model *= self.mspixelflat
                debug_model = self.rawflatimg.image.copy() / corr_model
                debug_model /= blaze_model
            import astropy.io.fits as fits
            hdu = fits.PrimaryHDU(debug_model)
            hdu.writeto('debug_model.fits', overwrite=True)

            # Shift to approximately constant wavelength
            shift_image = np.ones_like(self.rawflatimg.image.copy())
            # ratio = ratio of twilight to internal flats
            ratio = np.ones_like(self.rawflatimg.image.copy())  # placeholder... need to load "ratio" image from file
            for slit_idx in range(0, self.slits.spat_id.size):
                # Only use the overlapping regions of the slits, where the same wavelength range is covered
                onslit_init = (slitid_img_init == self.slits.spat_id[swslt[slit_idx]])
                onslit_olap = np.where(onslit_init & (waveimg >= minw) & (waveimg <= maxw))
                shifted = (onslit_olap[0] - onslit_olap[0].min(), onslit_olap[1],)
                shift_image[shifted] = ratio[onslit_olap]

        return scale_model


def show_flats(image_list, wcs_match=True, slits=None):
    """
    Interface to ginga to show a set of flat images

    Args:
        pixelflat (`numpy.ndarray`_):
        illumflat (`numpy.ndarray`_ or None):
        procflat (`numpy.ndarray`_):
        flat_model (`numpy.ndarray`_):
        spec_illum (`numpy.ndarray`_ or None):
        wcs_match (bool, optional):
        slits (:class:`pypeit.slittrace.SlitTraceSet`, optional):

    Returns:

    """
    display.connect_to_ginga(raise_err=True, allow_new=True)
    if slits is not None:
        left, right, mask = slits.select_edges()
        gpm = mask == 0
    # Loop me
    clear = True
    for img, name, cut in image_list:
        if img is None:
            continue
        # TODO: Add an option that shows the relevant stuff in a
        # matplotlib window.
        viewer, ch = display.show_image(img, chname=name, cuts=cut, wcs_match=wcs_match,
                                        clear=clear)
        if slits is not None:
            display.show_slits(viewer, ch, left[:, gpm], right[:, gpm],
                               slit_ids=slits.spat_id[gpm])
        # Turn off clear
        if clear:
            clear = False


def illum_profile_spectral(rawimg, waveimg, slits, model=None, gpmask=None, skymask=None, trim=3, flexure=None):
    """
    Generate a rough estimate of the relative spectral scaling of slits
    using a low order polynomial. This routine is for slit-based IFUs.

    Parameters
    ----------
    rawimg : `numpy.ndarray`_
        Image data that will be used to estimate the spectral relative sensitivity
    waveimg : `numpy.ndarray`_
        Wavelength image
    slits : :class:`pypeit.slittrace.SlitTraceSet`
        Information stored about the slits
    model : `numpy.ndarray`_, None
        A model of the rawimg data. If None, rawimg will be used.
    gpmask : `numpy.ndarray`_, None
        Good pixel mask
    skymask : `numpy.ndarray`_, None
        Sky mask
    trim : int
        Number of pixels to trim from the edges of the slit
        when deriving the spectral illumination
    flexure : float, None
        Spatial flexure

    Returns
    -------
    scale_model: `numpy.ndarray`_
        An image containing the appropriate scaling
    """
    msgs.info("Performing relative spectral sensitivity correction")
    # Setup some helpful parameters
    skymask_now = skymask if (skymask is not None) else np.ones_like(rawimg, dtype=bool)
    gpm = gpmask if (gpmask is not None) else np.ones_like(rawimg, dtype=bool)
    modelimg = model if (model is not None) else rawimg.copy()
    # Setup the slits
    slitid_img_init = slits.slit_img(pad=0, initial=True, flexure=flexure)
    slitid_img_trim = slits.slit_img(pad=-trim, initial=True, flexure=flexure)
    scaleImg = np.ones_like(rawimg)
    rawimg_copy = rawimg.copy()
    # Obtain the minimum and maximum wavelength of all slits
    mnmx_wv = np.zeros((slits.nslits, 2))
    for slit_idx, slit_spat in enumerate(slits.spat_id):
        onslit_init = (slitid_img_init == slit_spat)
        mnmx_wv[slit_idx, 0] = np.min(waveimg[onslit_init])
        mnmx_wv[slit_idx, 1] = np.max(waveimg[onslit_init])

    # Prepare reference spectrum
    specmin = np.argmin(mnmx_wv[:, 0])
    specmax = np.argmax(mnmx_wv[:, 1])
    dwav = np.max((mnmx_wv[:, 1] - mnmx_wv[:, 0])/slits.nspec)
    numsamp = int((np.max(mnmx_wv) - np.min(mnmx_wv)) / dwav)
    bins = np.linspace(np.min(mnmx_wv), np.max(mnmx_wv), numsamp)
    # Ease the minimum and maximum spectra into each other to create a smooth reference
    ww = np.where((bins > mnmx_wv[specmax, 0]) & (bins < mnmx_wv[specmin, 1]))  # Yes, this is correct
    easing = np.ones(numsamp)
    easing[ww] = 1 - np.linspace(0, 1, ww[0].size)
    easing[ww[0].max():] = 0
    onslit_specmin = (slitid_img_trim == slits.spat_id[specmin])
    onslit_specmax = (slitid_img_trim == slits.spat_id[specmax])
    weights = np.zeros(rawimg.shape)
    weights[onslit_specmin] = interpolate.interp1d(bins, easing, kind='linear', bounds_error=False,
                                                   fill_value="extrapolate")(waveimg[onslit_specmin])
    weights[onslit_specmax] = interpolate.interp1d(bins, 1-easing, kind='linear', bounds_error=False,
                                                   fill_value="extrapolate")(waveimg[onslit_specmax])
    # Generate a reference spectrum
    hist, edge = np.histogram(waveimg, bins=bins, weights=modelimg * weights)
    cntr, edge = np.histogram(waveimg, bins=bins, weights=weights)
    wave_ref = 0.5 * (edge[1:] + edge[:-1])
    spec_ref = hist / cntr

    # Go through the slits and calculate the overlapping flux
    maxiter = 10
    lo_prev, hi_prev = 1.0E-32, 1.0E32
    sn_smooth_npix = int(np.round(wave_ref.size / 10))
    for rr in range(maxiter):
        # Reset the relative scaling for this iteration
        relscl_model = np.ones_like(rawimg)

        # Temporary code
        for slit_idx in range(0, slits.spat_id.size):
            # Only use the overlapping regions of the slits, where the same wavelength range is covered
            onslit_b = (slitid_img_trim == slits.spat_id[slit_idx])
            onslit_b_init = (slitid_img_init == slits.spat_id[slit_idx])
            onslit_b_olap = onslit_b & gpm & (waveimg >= mnmx_wv[slit_idx, 0]) & (waveimg <= mnmx_wv[slit_idx, 1]) & skymask_now
            hist, edge = np.histogram(waveimg[onslit_b_olap], bins=bins, weights=rawimg_copy[onslit_b_olap])
            cntr, edge = np.histogram(waveimg[onslit_b_olap], bins=bins)
            cntr = cntr.astype(np.float)
            cntr *= spec_ref
            norm = (cntr != 0)/(cntr + (cntr == 0))
            arr = hist*norm
            gdmask = (arr != 0)
            relscale = coadd.smooth_weights(arr, gdmask, sn_smooth_npix)
            rescale_model = interpolate.interp1d(wave_ref, relscale, kind='linear', bounds_error=False,
                                                 fill_value="extrapolate")(waveimg[onslit_b_init])
            # Store the result
            relscl_model[onslit_b_init] = rescale_model.copy()

        minv, maxv = np.min(relscl_model), np.max(relscl_model)
        if 1/minv + maxv > lo_prev+hi_prev:
            # Adding noise, so break
            # NOTE : THe best precision one might hope for is about:
            # 1.4826 * MAD(arr) / np.sqrt(sn_smooth_npix/ 10)  # /10 comes from the coadd.smooth_weights function
            break
        else:
            lo_prev, hi_prev = 1/minv, maxv
        msgs.info("Iteration {0:d} :: Minimum/Maximum scales = {1:.5f}, {2:.5f}".format(rr + 1, minv, maxv))
        # Store rescaling
        scaleImg *= relscl_model
        rawimg_copy /= relscl_model
        if max(abs(1/minv), abs(maxv)) < 1.001:  # Relative accruacy of 0.1% is sufficient
            break
    return scaleImg


def merge(init_cls, merge_cls):
    """
    Merge merge_cls into init_cls, and return a merged :class:`pypeit.flatfield.FlatImages` class.
    If an element exists in both init_cls and merge_cls, the merge_cls value is taken

    Parameters
    ----------
    init_cls : :class:`pypeit.flatfield.FlatImages`
        Initial class (the elements of this class will be considered the default)
    merge_cls : :class:`pypeit.flatfield.FlatImages`
        The non-zero elements will be merged into init_cls.

    Returns
    -------
    :class:`pypeit.flatfield.FlatImages` : A new instance of the FlatImages class with merged properties.
    """
    # Check the class to be merged in is not None
    if merge_cls is None:
        return init_cls
    # Initialise variables
    # extract all elements that are prefixed with 'pixelflat_' or 'illumflat_'
    keys = [a for a in list(init_cls.__dict__.keys()) if '_' in a and a.split('_')[0] in ['illumflat', 'pixelflat']]
    dd = dict()
    for key in keys:
        dd[key] = None
    # Cherry pick the values from each class
    dd['PYP_SPEC'] = merge_cls.PYP_SPEC if init_cls.PYP_SPEC is None else init_cls.PYP_SPEC
    dd['spat_id'] = merge_cls.spat_id if init_cls.spat_id is None else init_cls.spat_id
    for key in keys:
        mrg = False
        val = None
        namespace = dict({'val': val, 'init_cls':init_cls, 'merge_cls':merge_cls, 'mrg':mrg})
        exec("val = init_cls.{0:s}".format(key), namespace)
        exec("mrg = merge_cls.{0:s} is not None".format(key), namespace)
        if namespace['mrg']:
            exec("val = merge_cls.{0:s}".format(key), namespace)
        dd[key] = namespace['val']
    # Construct the merged class
    return FlatImages(**dd)
