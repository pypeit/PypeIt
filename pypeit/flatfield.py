"""
Implements the flat-field class.

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import inspect
import numpy as np

from scipy import interpolate

from matplotlib import pyplot as plt

from IPython import embed

from pypeit import msgs
from pypeit import utils
from pypeit import ginga
from pypeit import bspline

from pypeit.par import pypeitpar
from pypeit import datamodel
from pypeit import masterframe
from pypeit.core import flat
from pypeit.core import tracewave
from pypeit.core import basis
from pypeit import slittrace


class FlatImages(datamodel.DataContainer):
    """
    Simple DataContainer for the output from BuildWaveTilts

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

    datamodel = {
        'procflat':  dict(otype=np.ndarray, atype=np.floating, desc='Processed, combined flats'),
        'pixelflat': dict(otype=np.ndarray, atype=np.floating, desc='Pixel normalized flat'),
        'flat_model': dict(otype=np.ndarray, atype=np.floating, desc='Model flat'),
        'PYP_SPEC': dict(otype=str, desc='PypeIt spectrograph name'),
        'bpmflats': dict(otype=np.ndarray, atype=np.integer,
                         desc='Mirrors SlitTraceSet mask for the Flat-specific flags'),
        'spat_bsplines': dict(otype=np.ndarray, atype=bspline.bspline,
                              desc='B-spline models for Illumination flat'),
        'spat_id': dict(otype=np.ndarray, atype=np.integer, desc='Slit spat_id '),
    }

    def __init__(self, procflat=None, pixelflat=None, bpmflats=None,
                 flat_model=None, spat_bsplines=None, PYP_SPEC=None, spat_id=None):
        # Parse
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=d)

    def _init_internals(self):
        # Master stuff
        self.master_key = None
        self.master_dir = None

    def _validate(self):
        #
        if self.spat_bsplines is not None and len(self.spat_bsplines) > 0:
            if len(self.spat_id) != len(self.spat_bsplines):
                msgs.error("Bsplines are out of sync with the slit IDs")

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

        # Rest of the datamodel
        for key in self.keys():
            # Skip None
            if self[key] is None:
                continue
            # Array?
            if self.datamodel[key]['otype'] == np.ndarray and key != 'spat_bsplines':
                d.append({key: self[key]})
            elif key == 'spat_bsplines':
                for ss, spat_bspl in enumerate(self[key]):
                    # Naming
                    spat_bspl.hdu_prefix = 'SPAT_ID-{}_'.format(self.spat_id[ss])
                    dkey = 'bspline-{}'.format(self.spat_id[ss])
                    # Save
                    d.append({dkey: spat_bspl})
            else: # Add to header of the primary image
                d[0][key] = self[key]
        # Return
        return d

    def fit2illumflat(self, slits, flexure_shift=None):
        """

        Args:
            slits (:class:`pypeit.slittrace.SlitTraceSet`):
            flexure_shift (float, optional):

        Returns:

        """
        illumflat = np.ones_like(self.procflat)
        # Loop
        for slit_idx in range(slits.nslits):
            # Skip masked
            if slits.mask[slit_idx] != 0:
                continue
            # Skip those without a bspline
            # DO it
            _slitid_img = slits.slit_img(slitidx=slit_idx, flexure=flexure_shift)
            onslit = _slitid_img == slits.spat_id[slit_idx]
            spat_coo = slits.spatial_coordinate_image(slitidx=slit_idx,
                                                      slitid_img=_slitid_img,
                                                      flexure_shift=flexure_shift)

            illumflat[onslit] = self.spat_bsplines[slit_idx].value(spat_coo[onslit])[0]
        # TODO -- Update the internal one?  Or remove it altogether??
        return illumflat

    @classmethod
    def _parse(cls, hdu, ext=None, transpose_table_arrays=False, debug=False,
               hdu_prefix=None):
        # Grab everything but the bspline's
        _d, dm_version_passed, dm_type_passed = super(FlatImages, cls)._parse(hdu)
        # Now the bsplines
        list_of_bsplines = []
        spat_ids = []
        for ihdu in hdu:
            if 'BSPLINE' in ihdu.name:
                ibspl = bspline.bspline.from_hdu(ihdu)
                if ibspl.version != bspline.bspline.version:
                    msgs.warn("Your bspline is out of date!!")
                list_of_bsplines.append(ibspl)
                # Grab SPAT_ID for checking
                i0 = ihdu.name.find('ID-')
                i1 = ihdu.name.find('_BSP')
                spat_ids.append(int(ihdu.name[i0+3:i1]))
        # Check
        if spat_ids != _d['spat_id'].tolist():
            msgs.error("Bad parsing of the MasterFlat")
        # Finish
        _d['spat_bsplines'] = np.asarray(list_of_bsplines)
        return _d, dm_version_passed, dm_type_passed

    def show(self, slits=None, wcs_match=True):
        """
        Simple wrapper to show_flats()

        Args:
            slits:
            wcs_match:

        Returns:

        """
        illumflat = None
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
            illumflat = self.fit2illumflat(slits)
        # Show
        show_flats(self.pixelflat, illumflat, self.procflat, self.flat_model,
                   wcs_match=wcs_match, slits=slits)

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

    Attributes:
        rawflatimg (:class:`pypeit.images.pypeitimage.PypeItImage`):
        mspixelflat (`numpy.ndarray`_):
            Normalized flat
        msillumflat (`numpy.ndarray`_):
            Illumination flat
        flat_model (`numpy.ndarray`_):
            Model of the flat
        list_of_spat_bsplines (list):

    """

    # Frame type is a class attribute
    frametype = 'pixelflat'
    master_type = 'Flat'


    def __init__(self, rawflatimg, spectrograph, flatpar, slits, wavetilts):

        # Defatuls
        self.spectrograph = spectrograph
        # FieldFlattening parameters
        self.flatpar = flatpar

        # Input data
        self.slits = slits
        self.wavetilts = wavetilts

        # Worth a check
        self.wavetilts.is_synced(self.slits)

        # Attributes unique to this Object
        self.rawflatimg = rawflatimg      # Un-normalized pixel flat as a PypeItImage
        self.mspixelflat = None     # Normalized pixel flat
        self.msillumflat = None     # Illumination flat
        self.flat_model = None      # Model flat
        self.list_of_spat_bsplines = None

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

#    def build_pixflat(self, trim=True, force=False):
#        """
#        Process the flat flat images.
#
#        Processing steps are the result of
#        :func:`pypeit.core.procimg.init_process_steps`, ``trim``
#        (based on the input argument), ``apply_gain``, and ``orient``.
#        Currently, cosmic-ray rejection (``cr_reject``) is not done.
#
#        Args:
#            trim (:obj:`bool`, optional):
#                Trim the image down to just the data section.
#            force (:obj:`bool`, optional):
#                Force the flat to be reconstructed if it already exists
#
#        Returns:
#            pypeitimage.PypeItImage:  The image with the unnormalized pixel-flat data.
#        """
#        if self.rawflatimg is None or force:
#            # Process steps
#            self.process_steps = procimg.init_process_steps(self.msbias, self.par['process'])
#            if trim:
#                self.process_steps += ['trim']
#            self.process_steps += ['apply_gain']
#            self.process_steps += ['orient']
#            # Turning this on leads to substantial edge-tracing problems when last tested
#            #     JXP November 22, 2019
#            #if self.par['cr_reject']:
#            #    self.process_steps += ['crmask']
#            self.steps.append(inspect.stack()[0][3])
#            # Do it
#            self.rawflatimg = super(FlatField, self).build_image(bias=self.msbias, bpm=self.msbpm,
#                                                                 ignore_saturation=True)
#        return self.rawflatimg

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
        # NOTE: Tilts do not change and self.slits is updated
        # internally.
        # TODO -- Move the infinitely long fit() code back here?
        self.fit(debug=debug)

        if show:
            # Global skysub is the first step in a new extraction so clear the channels here
            self.show(wcs_match=True)

        # Build the mask
        bpmflats = np.zeros_like(self.slits.mask, dtype=self.slits.bitmask.minimum_dtype())
        for flag in ['SKIPFLATCALIB', 'BADFLATCALIB']:
            bpm = self.slits.bitmask.flagged(self.slits.mask, flag)
            if np.any(bpm):
                bpmflats[bpm] = self.slits.bitmask.turn_on(bpmflats[bpm], flag)

        # Return
        return FlatImages(procflat=self.rawflatimg.image,
                          pixelflat=self.mspixelflat,
                          flat_model=self.flat_model,
                          spat_bsplines=np.asarray(self.list_of_spat_bsplines),
                          bpmflats=bpmflats,
                          PYP_SPEC=self.spectrograph.spectrograph,
                          spat_id=self.slits.spat_id)

    def show(self, wcs_match=True):
        """
        Show all of the flat field products in ginga.

        Args:
            wcs_match (:obj:`bool`, optional):
                Match the WCS of the flat-field images
        """
        # Get the slits
        show_flats(self.mspixelflat, self.msillumflat, self.rawflatimg.image, self.flat_model,
                   wcs_match=wcs_match, slits=self.slits)

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
            debug (:obj:`bool`, optional):
                Show plots useful for debugging. This will block
                further execution of the code until the plot windows
                are closed.

        """
        # TODO: break up this function!  Can it be partitioned into a series of "core" methods?
        # TODO: JFH I wrote all this code and will have to maintain it and I don't want to see it broken up.
        # TODO: JXP This definitely needs breaking up..

        # Init
        self.list_of_spat_bsplines = []

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

        # #################################################
        # Model each slit independently
        for slit_idx, slit_spat in enumerate(self.slits.spat_id):
            # Is this a good slit??
            if self.slits.mask[slit_idx] != 0:
                msgs.info('Skipping bad slit: {}'.format(slit_spat))
                self.list_of_spat_bsplines.append(bspline.bspline(None))
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
                self.list_of_spat_bsplines.append(bspline.bspline(None))
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
            spec_bspl, spec_gpm_fit, spec_flat_fit, _, exit_status \
                    = utils.bspline_profile(spec_coo_data, spec_flat_data, spec_ivar_data,
                                            np.ones_like(spec_coo_data), ingpm=spec_gpm_data,
                                            nord=4, upper=logrej, lower=logrej,
                                            kwargs_bspline={'bkspace': spec_samp_fine},
                                            kwargs_reject={'groupbadpix': True, 'maxrej': 5})

            if exit_status > 1:
                # TODO -- MAKE A FUNCTION
                msgs.warn('Flat-field spectral response bspline fit failed!  Not flat-fielding '
                          'slit {0} and continuing!'.format(slit_spat))
                self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(self.slits.mask[slit_idx], 'BADFLATCALIB')
                self.list_of_spat_bsplines.append(bspline.bspline(None))
                continue

            # Debugging/checking spectral fit
            if debug:
                utils.bspline_qa(spec_coo_data, spec_flat_data, spec_bspl, spec_gpm_fit,
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
                ax = utils.bspline_qa(spat_coo_data, spat_flat_data, spat_bspl, spat_gpm_fit,
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
                self.list_of_spat_bsplines.append(spat_bspl)
            else:
                # Save the nada
                msgs.warn('Slit illumination profile bspline fit failed!  Spatial profile not '
                          'included in flat-field model for slit {0}!'.format(slit_spat))
                self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(self.slits.mask[slit_idx], 'BADFLATCALIB')
                self.list_of_spat_bsplines.append(bspline.bspline(None))
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
            twod_bspl, twod_gpm_fit, twod_flat_fit, _ , exit_status \
                    = utils.bspline_profile(twod_spec_coo_data, twod_flat_data, twod_ivar_data,
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

        # Set the pixelflat to 1.0 wherever the flat was nonlinear
        self.mspixelflat[rawflat >= nonlinear_counts] = 1.0
        # Set the pixelflat to 1.0 within trim pixels of all the slit edges
        trimmed_slitid_img_new = self.slits.slit_img(pad=-trim, initial=False)
        tweaked_slitid_img = self.slits.slit_img(initial=False)
        self.mspixelflat[(trimmed_slitid_img_new < 0) & (tweaked_slitid_img > 0)] = 1.0


        # Do not apply pixelflat field corrections that are greater than
        # 100% to avoid creating edge effects, etc.
        self.mspixelflat = np.clip(self.mspixelflat, 0.5, 2.0)

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
        # don't have to instantiate invvar and profile_basis
        spat_bspl, spat_gpm_fit, spat_flat_fit, _, exit_status \
            = utils.bspline_profile(spat_coo_data, spat_flat_data,
                                    np.ones_like(spat_flat_data),
                                    np.ones_like(spat_flat_data), nord=4, upper=5.0,
                                    lower=5.0, fullbkpt=spat_bspl.breakpoints)
        # Return
        return exit_status, spat_coo_data, spat_flat_data, spat_bspl, spat_gpm_fit, \
               spat_flat_fit, spat_flat_data_raw


def show_flats(pixelflat, illumflat, procflat, flat_model, wcs_match=True, slits=None):
    """
    Interface to ginga to show a set of flat images

    Args:
        pixelflat (`numpy.ndarray`_):
        illumflat (`numpy.ndarray`_ or None):
        procflat (`numpy.ndarray`_):
        flat_model (`numpy.ndarray`_):
        wcs_match (bool, optional):
        slits (:class:`pypeit.slittrace.SlitTraceSet`, optional):

    Returns:

    """
    ginga.connect_to_ginga(raise_err=True, allow_new=True)
    if slits is not None:
        left, right, mask = slits.select_edges()
        gpm = mask == 0
    # Loop me
    clear = True
    for img, name, cut in zip([pixelflat, illumflat, procflat, flat_model],
                         ['pixelflat', 'illumflat', 'flat', 'flat_model'],
                         [(0.9, 1.1), (0.9, 1.1), None, None]):
        if img is None:
            continue
        # TODO: Add an option that shows the relevant stuff in a
        # matplotlib window.
        viewer, ch = ginga.show_image(img, chname=name, cuts=cut,
                                  wcs_match=wcs_match, clear=clear)
        if slits is not None:
            ginga.show_slits(viewer, ch, left[:, gpm], right[:, gpm],
                             slit_ids=slits.spat_id[gpm])
        # Turn off clear
        if clear:
            clear = False

