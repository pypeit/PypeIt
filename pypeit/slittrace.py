"""
Implements the objects used to hold slit edge data.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import inspect

from IPython import embed

import numpy as np

from pypeit import msgs
from pypeit import datamodel
from pypeit.bitmask import BitMask


class SlitTraceBitMask(BitMask):
    """
    Mask bits used during slit tracing.
    """
    version = '1.0.0'

    def __init__(self):
        mask = dict([
            ('SHORTSLIT', 'Slit formed by left and right edge is too short. Not ignored for flexure'),
            ('BOXSLIT', 'Slit formed by left and right edge is valid (large enough to be a valid '
                        'slit), but too short to be a science slit'),
            ('USERIGNORE', 'User has specified to ignore this slit. Not ignored for flexure.'),
            ('BADWVCALIB', 'Wavelength calibration failed for this slit'),
            ('BADTILTCALIB', 'Tilts analysis failed for this slit'),
            ('SKIPFLATCALIB', 'Flat field generation failed for this slit. Skip flat fielding'),
            ('BADFLATCALIB', 'Flat field generation failed for this slit. Ignore it fully.'),
            ('BADREDUCE', 'Skysub/extraction failed for this slit'),
        ])
        super(SlitTraceBitMask, self).__init__(list(mask.keys()), descr=list(mask.values()))

    @property
    def exclude_for_reducing(self):
        # Ignore these flags when reducing or considering reduced slits
        return ['SKIPFLATCALIB']

    @property
    def exclude_for_flexure(self):
        # Ignore these flags when performing a flexure calculation
        #  Currently they are *all* of the flags..
        return ['SHORTSLIT', 'USERIGNORE', 'BADWVCALIB', 'BADTILTCALIB',
                'SKIPFLATCALIB', 'BADFLATCALIB', 'BADREDUCE']



class SlitTraceSet(datamodel.DataContainer):
    """
    Defines a generic class for holding and manipulating image traces
    organized into left-right slit pairs.

    Instantiation arguments map directly to the object
    :attr:`datamodel`. The only additional argument is ``load``,
    described below.

    :class:`SlitTraceSet` objects can be instantiated with only the
    master-frame arguments. If this is done, it's expected that
    you'll attempt to load an existing master frame, either using
    ``load=True`` on instantiation or by a call to :func:`load`.
    Otherwise, all the elements of the data model will be empty.

    Args:
        load (:obj:`bool`, optional):
            Attempt to load an existing master frame with the slit
            trace data. WARNING: If ``load`` is True, all of the slit
            information provided to the initialization is ignored!
            Only the arguments relevant to the
            :class:`pypeit.masterframe.MasterFrame` components of
            :class:`SlitTraceSet` are used.

    Attributes:
        left_flexure (`numpy.ndarray`_):  Convenient spot to hold flexure corrected left
        right_flexure (`numpy.ndarray`_):  Convenient spot to hold flexure corrected right
        master_key (:obj:`str`):
        master_dir (:obj:`str`):

    """
    frametype = 'slits'
    master_type = 'Slits'
    """Name for type of master frame."""
    master_file_format = 'fits.gz'
    """File format for the master frame file."""
    minimum_version = '1.1.0'
    version = '1.1.0'
    """SlitTraceSet data model version."""

    hdu_prefix = None
    output_to_disk = None

    bitmask = SlitTraceBitMask()

    # Define the data model
    datamodel = {'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 'pypeline': dict(otype=str, descr='PypeIt pypeline name'),
                 'nspec': dict(otype=int,
                               descr='Number of pixels in the image spectral direction.'),
                 'nspat': dict(otype=int,
                               descr='Number of pixels in the image spatial direction.'),
                 'binspec': dict(otype=int,
                                 descr='Number of pixels binned in the spectral direction.'),
                 'binspat': dict(otype=int,
                                 descr='Number of pixels binned in the spatial direction.'),
                 'pad': dict(otype=int,
                             descr='Integer number of pixels to consider beyond the slit edges.'),
                 'spat_id': dict(otype=np.ndarray, atype=(int,np.integer),
                                 descr='Slit ID number from SPAT measured at half way point.'),
                 'maskdef_id': dict(otype=np.ndarray, atype=(int,np.integer),
                                    descr='Slit ID number slitmask'),
                 'ech_order': dict(otype=np.ndarray, atype=(int,np.integer),
                                   descr='Slit ID number echelle order'),
                 'nslits': dict(otype=int,
                                descr='Total number of slits, derived from shape of left_init.'),
                 'left_init': dict(otype=np.ndarray, atype=np.floating,
                              descr='Spatial coordinates (pixel indices) of all left edges, one '
                                    'per slit.  Derived from the TraceImage. Shape is Nspec by '
                                    'Nslits.'),
                 'right_init': dict(otype=np.ndarray, atype=np.floating,
                              descr='Spatial coordinates (pixel indices) of all right edges, one '
                                    'per slit.  Derived from the TraceImage. Shape is Nspec by '
                                    'Nslits.'),
                 'left_tweak': dict(otype=np.ndarray, atype=np.floating,
                                    descr='Spatial coordinates (pixel indices) of all left '
                                          'edges, one per slit.  These traces have been adjusted '
                                          'by the flat-field.  Shape is Nspec by Nslits.'),
                 'right_tweak': dict(otype=np.ndarray, atype=np.floating,
                                     descr='Spatial coordinates (pixel indices) of all right '
                                           'edges, one per slit.  These traces have been adjusted '
                                           'by the flat-field.  Shape is Nspec by Nslits.'),
                 'center': dict(otype=np.ndarray, atype=np.floating,
                               descr='Spatial coordinates of the slit centers from left_init, '
                                     'right_init.  Shape is Nspec by Nslits.'),
                 'mask_init': dict(otype=np.ndarray, atype=np.integer,
                                   descr='Bit mask for slits at instantiation.  Used to reset'),
                 'mask': dict(otype=np.ndarray, atype=np.integer,
                              descr='Bit mask for slits (fully good slits have 0 value).  Shape '
                                    'is Nslits.'),
                'slitbitm': dict(otype=str, desc='List of BITMASK keys from SlitTraceBitMask'),
                'specmin': dict(otype=np.ndarray, atype=np.floating,
                                descr='Minimum spectral position allowed for each slit/order.  '
                                      'Shape is Nslits.'),
                'specmax': dict(otype=np.ndarray, atype=np.floating,
                                descr='Maximum spectral position allowed for each slit/order.  '
                                      'Shape is Nslits.')}
    """Provides the class data model."""

    # TODO: Allow tweaked edges to be arguments?
    # TODO: May want nspat to be a required argument.
    # The INIT must contain every datamodel item or risk fail on I/O when it is a nested container
    def __init__(self, left_init, right_init, pypeline, nspec=None, nspat=None, PYP_SPEC=None,
                 mask_init=None, specmin=None, specmax=None, binspec=1, binspat=1, pad=0,
                 spat_id=None, maskdef_id=None,
                 ech_order=None, nslits=None, left_tweak=None,
                 right_tweak=None, center=None, mask=None, slitbitm=None):

        # Instantiate the DataContainer
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = dict([(k,values[k]) for k in args[1:]])
        # The dictionary passed to DataContainer.__init__ does not
        # contain self or the MasterFrame arguments.
        # TODO: Does it matter if the calling function passes the
        # keyword arguments in a different order?
        datamodel.DataContainer.__init__(self, d=_d)

    def _validate(self):
        """
        Validate the slit traces.
        """
        # Allow the object to be empty
        if self.left_init is None or self.right_init is None:
            return
        if self.left_init.shape != self.right_init.shape:
            raise ValueError('Input left_init and right_init traces should have the same shape.')
        if self.left_init.ndim == 1:
            # Object contains a single slit.  Expand the dimensions.
            self.left_init = np.expand_dims(self.left_init, 1)
            self.right_init = np.expand_dims(self.right_init, 1)

        # Do the same for the tweaked traces. NOTE: Tweaked traces will
        # only exist if they were read from an output file; i.e., the
        # current init does not allow for the tweaked traces to be
        # defined directly.
        if self.left_tweak is not None:
            if self.left_tweak.shape != self.right_tweak.shape:
                # TODO: Shouldn't get here, but just in case the init changes...
                raise ValueError('Tweaked left and right traces should have the same shape.')
            if self.left_tweak.ndim == 1:
                # Object contains a single slit.  Expand the dimensions.
                self.left_tweak = np.expand_dims(self.left_tweak, 1)
                self.right_tweak = np.expand_dims(self.right_tweak, 1)

        self.nspec, self.nslits = self.left_init.shape

        # Center is always defined by the original traces, not the
        # tweaked ones. This means that the slit IDs are always tied
        # to the original traces, not the tweaked ones.
        self.center = (self.left_init+self.right_init)/2

        if self.nspat is None:
            # TODO: May want nspat to be a required argument given the
            # only other option is this kludge, which should basically
            # never be useful.
            self.nspat = int(np.amax(np.append(self.left_init, self.right_init)))
        if self.spat_id is None:
            self.spat_id = np.round(self.center[int(np.round(0.5 * self.nspec)), :]).astype(int)
        if self.PYP_SPEC is None:
            self.PYP_SPEC = 'unknown'
        if self.mask_init is None:
            self.mask_init = np.zeros(self.nslits, dtype=self.bitmask.minimum_dtype())
        if self.specmin is None:
            self.specmin = np.full(self.nslits, -1, dtype=float)
        if self.specmax is None:
            self.specmax = np.full(self.nslits, self.nspec, dtype=float)

        # If the echelle order is provided, check that the number of
        # orders matches the number of provided "slits"
        if self.ech_order is not None and len(self.ech_order) != self.nslits:
            msgs.error('Number of provided echelle orders does not match the number of '
                       'order traces.')

        # Make sure mask, specmin, and specmax are at least 1D arrays.
        # TODO: Is there a way around this?
        self.mask_init = np.atleast_1d(self.mask_init)
        self.specmin = np.atleast_1d(self.specmin)
        self.specmax = np.atleast_1d(self.specmax)
        if self.slitbitm is None:
            self.slitbitm = ','.join(list(self.bitmask.keys()))
        else:
            # Validate
            if self.slitbitm != ','.join(list(self.bitmask.keys())):
                msgs.error("Input BITMASK keys differ from current data model!")
        # Mask
        if self.mask is None:
            self.mask = self.mask_init.copy()

    def _init_internals(self):
        self.left_flexure = None
        self.right_flexure = None
        # Master stuff
        self.master_key = None
        self.master_dir = None

    def _bundle(self):
        """
        Bundle the data in preparation for writing to a fits file.

        See :func:`pypeit.datamodel.DataContainer._bundle`. Data is
        always written to a 'SLITS' extension.
        """
        return super(SlitTraceSet, self)._bundle(ext='SLITS', transpose_arrays=True)

    @classmethod
    def _parse(cls, hdu, hdu_prefix=None):
        """
        Parse the data that was previously written to a fits file.

        See :func:`pypeit.datamodel.DataContainer._parse`. Data is
        always read from the 'SLITS' extension.
        """
        return super(SlitTraceSet, cls)._parse(hdu, ext='SLITS', transpose_table_arrays=True)

    def init_tweaked(self):
        """
        Initialize the tweaked slits.
        """
        self.left_tweak = self.left_init.copy()
        self.right_tweak = self.right_init.copy()

    def rm_tweaked(self):
        """
        Delete (set to None) the tweaked traces.
        """
        self.left_tweak = None
        self.right_tweak = None

    @property
    def slit_info(self):
        """

        Returns:
            `numpy.ndarray`_:

        """
        #
        info = np.vstack([self.spat_id, self.mask])
        info = np.vstack([info, np.zeros_like(self.spat_id)]) \
                    if self.maskdef_id is None else np.vstack([info, self.maskdef_id])
        return info.astype(int).T

    # TODO: Do we need both of these? I.e., can the 'spat_id' for
    # echelle spectrographs just be the echelle order?
    @property
    def slitord_id(self):
        """
        Return array of slit_spatId (MultiSlit, IFU) or ech_order (Echelle) values

        Returns:
            `numpy.ndarray`_:

        """
        if self.pypeline in ['MultiSlit', 'IFU']:
            return self.spat_id
        elif self.pypeline == 'Echelle':
            return self.ech_order
        else:
            msgs.error('Unrecognized Pypeline {:}'.format(self.pypeline))

    def spatid_to_zero(self, spat_id):
        """
        Convert input spat_id into a zero-based index

        Args:
            spat_id (int):

        Returns:
            int: zero-based index of the input spat_id

        """
        return np.where(self.spat_id == spat_id)[0][0]

    def slitord_to_zero(self, slitord):
        """
        Convert input slitord into a zero-based index

        Args:
            slitord (int):

        Returns:
            int: zero-based index of the input spat_id

        """
        if self.pypeline in ['MultiSlit', 'IFU']:
            return np.where(self.spat_id == slitord)[0][0]
        elif self.pypeline in ['Echelle']:
            return np.where(self.ech_order == slitord)[0][0]
        else:
            msgs.error('Unrecognized Pypeline {:}'.format(self.pypeline))

    def get_slitlengths(self, initial=False, median=False):
        """
        Get the length of each slit in pixels.

        By default, the method will return the tweaked slit lengths
        if they have been defined. If they haven't been defined the
        nominal edges (:attr:`left` and :attr:`right`) are returned.
        Use ``initial=True`` to return the nominal edges regardless
        of the presence of the tweaked edges.

        Args:
            initial (:obj:`bool`, optional):
                To use the initial edges regardless of the presence of
                the tweaked edges, set this to True.
            median (:obj:`bool`, optional):
                The default is to return the slit length as a function
                of the spectral coordinate. If median is set to true,
                the median slit length of each slit is returned.

        Returns:
            `numpy.ndarray`_: Slit lengths.
        """
        left, right, _ = self.select_edges(initial=initial)
        slitlen = right - left
        if median is True:
            slitlen = np.median(slitlen, axis=1)
        return slitlen

    def get_radec_image(self, wcs, initial=True, flexure=None, trace_cen=None):
        """Generate an RA and DEC image for every pixel in the frame

        Parameters
        ----------
        wcs : astropy.wcs
            The World Coordinate system of a science frame
        maxslitlen : int
            This is the slit length in pixels, and it should be the same
            value that was passed to get_wcs() to generate the WCS that
            is passed into this function as an argument.
        initial : bool
            Select the initial slit edges?
        flexure : float, optional
            If provided, offset each slit by this amount.
        trace_cen : `numpy.ndarray`_, optional
            Central traces of each slit. Shape should be (slits.nspec, slits.nslits).
            If None, the average of the left and right slit edges will be used

        Returns
        -------
        tuple : There are three elements in the tuple. The first two are 2D numpy arrays
                of shape (nspec, nspat), where the first ndarray is the RA image, and the
                second ndarray is the DEC image. RA and DEC are in units degrees. The third
                element of the tuple stores the minimum and maximum difference (in pixels)
                between the WCS reference (usually the centre of the slit) and the edge of
                the slits. The third array has a shape of (nslits, 2).
        """
        # Grab the central trace, if none was provided
        if trace_cen is None:
            left, right, _ = self.select_edges(initial=initial, flexure=flexure)
            trace_cen = 0.5 * (left + right)

        # Initialise the output
        raimg = np.zeros((self.nspec, self.nspat))
        decimg = np.zeros((self.nspec, self.nspat))
        minmax = np.zeros((self.nslits, 2))

        # Get the slit information
        slitid_img_init = self.slit_img(pad=0, initial=initial, flexure=flexure)
        for slit_idx, spatid in enumerate(self.spat_id):
            onslit = (slitid_img_init == spatid)
            onslit_init = np.where(onslit)
            evalpos = onslit_init[1] - trace_cen[onslit_init[0], slit_idx]
            minmax[:, 0] = np.min(evalpos)
            minmax[:, 1] = np.max(evalpos)
            slitID = np.ones(evalpos.size) * slit_idx - wcs.wcs.crpix[0]
            world_ra, world_dec, _ = wcs.wcs_pix2world(slitID, evalpos, onslit_init[0], 0)
            # Set the RA first and DEC next
            raimg[onslit] = world_ra.copy()
            decimg[onslit] = world_dec.copy()
        return raimg, decimg, minmax

    def select_edges(self, initial=False, flexure=None):
        """
        Select between the initial or tweaked slit edges and allow for
        flexure correction.

        By default, the method will return the tweaked slits if they
        have been defined. If they haven't been defined the nominal
        edges (:attr:`left` and :attr:`right`) are returned. Use
        ``original=True`` to return the nominal edges regardless of
        the presence of the tweaked edges.

        Args:
            initial (:obj:`bool`, optional):
                To use the initial edges regardless of the presence of
                the tweaked edges, set this to True.
            flexure (:obj:`float`, optional):
                If provided, offset each slit by this amount

        Returns:
            tuple: Returns the full arrays containing the left and right
            edge coordinates and the mask, respectively.
            These are returned as copies.
        """
        # TODO: Add a copy argument?
        if self.left_tweak is not None and self.right_tweak is not None and not initial:
            left, right = self.left_tweak, self.right_tweak
        else:
            left, right = self.left_init, self.right_init

        # Add in spatial flexure?
        if flexure:
            self.left_flexure = left + flexure
            self.right_flexure = right + flexure
            left, right = self.left_flexure, self.right_flexure

        # Return
        return left.copy(), right.copy(), self.mask.copy()

    def slit_img(self, pad=None, slitidx=None, initial=False, flexure=None,
                 exclude_flag=None, use_spatial=True):
        r"""
        Construct an image identifying each pixel with its associated
        slit.

        The output image has the same shape as the original trace
        image. Each pixel in the image is set to the index
        of its associated slit (i.e, the pixel value is
        :math:`0..N_{\rm slit}-1`). Pixels not associated with any
        slit are given values of -1.

        The width of the slit is extended at either edge by a fixed
        number of pixels using the `pad` parameter in :attr:`par`.
        This value can be overridden using the method keyword
        argument.

        .. warning::

            - The function does not check that pixels end up in
              multiple pixels or that the padding is sensible given
              the separation between slit edges!
            - All slits are identified in the image, even if they are
              masked with :attr:`mask`.

        Args:
            pad (:obj:`float`, :obj:`int`, :obj:`tuple`, optional):
                The number of pixels used to pad (extend) the edge of
                each slit. This can be a single scale to pad both
                left and right edges equally or a 2-tuple that
                provides separate padding for the left (first
                element) and right (2nd element) edges separately. If
                not None, this overrides the value in :attr:`par`.
                The value can be negative, which means that the
                widths are **trimmed** instead of padded.
            slitidx (:obj:`int`, array_like, optional):
                List of indexes (zero-based) to include in the image.
                If None, all slits not flagged are included.
            initial (:obj:`bool`, optional):
                By default, the method will use the tweaked slit
                edges if they have been defined. If they haven't
                been, the initial edges (:attr:`left_init` and
                :attr:`right_init`) are used. To use the nominal edges
                regardless of the presence of the tweaked edges, set
                this to True. See :func:`select_edges`.
            exclude_flag (:obj:`str`, optional):
                Bitmask flag to ignore when masking
                Warning -- This could conflict with input slitids, i.e. avoid using both
            use_spatial (bool, optional):
                If True, use self.spat_id value instead of 0-based indices


        Returns:
            `numpy.ndarray`_: The image with the slit index
            identified for each pixel.
        """
        #
        if slitidx is not None and exclude_flag is not None:
            msgs.error("Cannot pass in both slitidx and exclude_flag!")
        # Check the input
        if pad is None:
            pad = self.pad
        _pad = pad if isinstance(pad, tuple) else (pad,pad)
        if len(_pad) != 2:
            msgs.error('Padding for both left and right edges should be provided as a 2-tuple!')

        # Pixel coordinates
        spat = np.arange(self.nspat)
        spec = np.arange(self.nspec)

        left, right, _ = self.select_edges(initial=initial, flexure=flexure)

        # Choose the slits to use
        if slitidx is not None:
            slitidx = np.atleast_1d(slitidx).ravel()
        else:
            bpm = self.mask.astype(bool)
            if exclude_flag:
                bpm &= np.invert(self.bitmask.flagged(self.mask, flag=exclude_flag))
            slitidx = np.where(np.invert(bpm))[0]

        # TODO: When specific slits are chosen, need to check that the
        # padding doesn't lead to slit overlap.

        # Find the pixels in each slit, limited by the minimum and
        # maximum spectral position.
        slitid_img = np.full((self.nspec,self.nspat), -1, dtype=int)
        for i in slitidx:
            slit_id = self.spat_id[i] if use_spatial else i
            indx = (spat[None,:] > left[:,i,None] - _pad[0]) \
                        & (spat[None,:] < right[:,i,None] + _pad[1]) \
                        & (spec > self.specmin[i])[:,None] & (spec < self.specmax[i])[:,None]
            slitid_img[indx] = slit_id
        # Return
        return slitid_img

    def spatial_coordinate_image(self, slitidx=None, full=False, slitid_img=None,
                                 pad=None, initial=False, flexure_shift=None):
        r"""
        Generate an image with the normalized spatial coordinate
        within each slit.

        Args:
            slitidx (:obj:`int`, array_like, optional):
                List of indices to include in the image. If None,
                all slits are included.
            full (:obj:`bool`, optional):
                If True, return a full image with the coordinates
                based on the single provided slit. In this case,
                slitids must be provided and it can be only one slit!
                If False, only those coordinates falling in the slit
                (as determined by ``slitid_image`` or by a call to
                :func:`slit_img`) will be returned.
            slitid_img (`numpy.ndarray`_, optional):
                Image identifying the slit associated with each
                pixel. Should have the same shape as expected by the
                object (:attr:`nspec` by :attr:`nspat`). Pixels not
                associated with any pixel have a value of -1. This is
                a convenience parameter that allows previous
                executions of :func:`slit_img` to be passed in
                directly instead of having to repeat the operation. A
                slit ID image is only required if ``full`` is False.
            pad (:obj:`float`, :obj:`int`, optional):
                When constructing the slit ID image, use this value
                to override the value of `pad` in :attr:`par`. Only
                used if ``slitid_img`` is not provided directly and
                ``full`` is False.
            initial (:obj:`bool`, optional):
                By default, the method will use the tweaked slit
                edges if they have been defined. If they haven't
                been, the nominal edges (:attr:`left` and
                :attr:`right`) are used. To use the nominal edges
                regardless of the presence of the tweaked edges, set
                this to True. See :func:`select_edges`.

        Returns:
            `numpy.ndarray`_: Array specifying the spatial coordinate
            of pixel in its own slit, scaled to go from 0 to 1. If
            ``full`` is True, the image provides the coordinates
            relative to the left edge for the provided slit over the
            full image.
        """
        # Slit indices to include
        _slitidx = np.arange(self.nslits) if slitidx is None else np.atleast_1d(slitidx).ravel()
        if full and len(_slitidx) > 1:
            msgs.error('For a full image with the slit coordinates, must select a single slit.')

        # Generate the slit ID image if it wasn't provided
        if not full:
            if slitid_img is None:
                slitid_img = self.slit_img(pad=pad, slitidx=_slitidx, initial=initial)
            if slitid_img.shape != (self.nspec,self.nspat):
                msgs.error('Provided slit ID image does not have the correct shape!')

        # Choose the slit edges to use
        left, right, _ = self.select_edges(initial=initial, flexure=flexure_shift)

        # Slit width
        slitwidth = right - left

        # TODO: This check should go in it's own function and/or
        # checked when instantiating the object.
        # Slit widths have to be larger than 0
        indx = slitwidth <= 0.
        if np.any(indx[:,_slitidx]):
            bad_slits = np.where(np.any(indx, axis=0))[0]
            # TODO: Shouldn't this fault?
            msgs.warn('Slits {0} have negative (or 0) slit width!'.format(bad_slits))

        # Output image
        coo_img = np.zeros((self.nspec,self.nspat), dtype=float)
        spat = np.arange(self.nspat)
        for i in _slitidx:
            coo = (spat[None,:] - left[:,i,None])/slitwidth[:,i,None]
            if not full:
                indx = slitid_img == self.spat_id[i]
                coo_img[indx] = coo[indx]
            else:
                coo_img = coo
        return coo_img

    def spatial_coordinates(self, initial=False, flexure=None):
        """
        Return a fiducial coordinate for each slit.

        This is a simple wrapper for :func:`select_edges` and
        :func:`slit_spat_pos`.

        Args:
            original (:obj:`bool`, optional):
                By default, the method will use the tweaked slit
                edges if they have been defined. If they haven't
                been, the nominal edges (:attr:`left` and
                :attr:`right`) are used. To use the nominal edges
                regardless of the presence of the tweaked edges, set
                this to True. See :func:`select_edges`.

        Returns:
            `numpy.ndarray`_: Vector with the list of floating point
            spatial coordinates.
        """
        # TODO -- Confirm it makes sense to pass in flexure
        left, right, _ = self.select_edges(initial=initial, flexure=flexure)
        return SlitTraceSet.slit_spat_pos(left, right, self.nspat)

    @staticmethod
    def slit_spat_pos(left, right, nspat):
        r"""
        Return a fidicial, normalized spatial coordinate for each slit.

        The fiducial coordinates are given by::
   
            nspec = left.shape[0]
            (left[nspec//2,:] + right[nspec//2,:])/2/nspat

        Args:
            left (`numpy.ndarray`_):
                Array with left slit edges. Shape is :math:`(N_{\rm
                spec},N_{\rm slits})`.
            right (`numpy.ndarray`_):
                Array with right slit edges. Shape is :math:`(N_{\rm
                spec},N_{\rm slits})`.
            nspat (:obj:`int`):
                Number of pixels in the spatial direction in the image
                used to trace the slit edges.

        Returns:
            `numpy.ndarray`_: Vector with the list of floating point
            spatial coordinates.
        """
        if left.shape != right.shape:
            msgs.error('Left and right traces must have the same shape.')
        nspec = left.shape[0]
        return (left[nspec//2,:] + right[nspec//2,:])/2/nspat

    def user_mask(self, det, slitspatnum):
        """
        Mask all but the input slit

        Args:
            det (:obj:`int`): Detector number
            slitspatnum (:obj:`str` or :obj:`list`):
        """
        # Parse
        dets, spat_ids = parse_slitspatnum(slitspatnum)
        if det not in dets:
            return
        # Cut down for convenience
        indet = dets == det
        spat_ids = spat_ids[indet]
        #
        msk = np.ones(self.nslits, dtype=bool)
        for slit_spat in spat_ids:
            #TODO -- Consider putting in a tolerance which if not met causes a crash
            idx = np.argmin(np.abs(self.spat_id - slit_spat))
            msk[idx] = False
        self.mask[msk] = self.bitmask.turn_on(self.mask[msk], 'USERIGNORE')

    def mask_flats(self, flatImages):
        """
        Mask from a :class:`pypeit.flatfield.Flats` object

        Args:
            flatImages (:class:`pypeit.flatfield.FlatImages`):

        """
        # Loop on all the FLATFIELD BPM keys
        for flag in ['SKIPFLATCALIB', 'BADFLATCALIB']:
            bad_flats = self.bitmask.flagged(flatImages.get_bpmflats(), flag)
            if np.any(bad_flats):
                self.mask[bad_flats] = self.bitmask.turn_on(self.mask[bad_flats], flag)

    def mask_wvcalib(self, wv_calib):
        """
        Mask from a WaveCalib object

        Args:
            wv_calib (:obj:`dict`):

        """
        for islit in range(self.nslits):
            if wv_calib.wv_fits[islit] is None or wv_calib.wv_fits[islit].pypeitfit is None:
                self.mask[islit] = self.bitmask.turn_on(self.mask[islit], 'BADWVCALIB')

    def mask_wavetilts(self, waveTilts):
        """
        Mask from a :class:`pypeit.wavetilts.WaveTilts` object

        Args:
            waveTilts (:class:`pypeit.wavetilts.WaveTilts`):

        """
        # There is only one BPM for Tilts (so far)
        bad_tilts = waveTilts.bpmtilts > 0
        if np.any(bad_tilts):
            self.mask[bad_tilts] = self.bitmask.turn_on(self.mask[bad_tilts], 'BADTILTCALIB')


def parse_slitspatnum(slitspatnum):
    """
    Parse the slitspatnum into a list of detectors and SPAT_IDs

    Args:
        slitspatnum (:obj:`str` or :obj:`list`:

    Returns:
        tuple:  dets, spat_ids  (each is an `numpy.ndarray`_ of int's)

    """
    dets = []
    spat_ids = []
    for item in slitspatnum.split(','):
        spt = item.split(':')
        dets.append(int(spt[0]))
        spat_ids.append(int(spt[1]))
    # Return
    return np.array(dets).astype(int), np.array(spat_ids).astype(int)
