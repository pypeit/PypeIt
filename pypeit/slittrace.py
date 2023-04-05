"""
Implements the objects used to hold slit edge data.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import inspect

from IPython import embed

import numpy as np

from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
from astropy import units
from astropy.stats import sigma_clipped_stats
from scipy.interpolate import RegularGridInterpolator, interp1d

from pypeit import msgs
from pypeit import datamodel
from pypeit import calibframe
from pypeit import specobj
from pypeit.bitmask import BitMask
from pypeit.core import parse


class SlitTraceBitMask(BitMask):
    """
    Mask bits used during slit tracing.
    """
    version = '1.0.1'

    def __init__(self):
        # Only ever append new bits (and don't remove old ones)
        mask = dict([
            ('SHORTSLIT', 'Slit formed by left and right edge is too short. Not ignored for flexure'),
            ('BOXSLIT', 'Slit formed by left and right edge is valid (large enough to be a valid '
                        'slit), but too short to be a science slit'),
            ('USERIGNORE', 'User has specified to ignore this slit. Not ignored for flexure.'),
            ('BADWVCALIB', 'Wavelength calibration failed for this slit'),
            ('BADTILTCALIB', 'Tilts analysis failed for this slit'),
            ('SKIPFLATCALIB', 'Flat field generation failed for this slit. Skip flat fielding'),
            ('BADFLATCALIB', 'Flat field generation failed for this slit. Ignore it fully.'),
            ('BADREDUCE', 'Reduction failed for this slit'), # THIS IS DEPRECATED (we may remove in v1.13) BUT STAYS HERE TO ALLOW FOR BACKWARDS COMPATIBILITY
            ('BADSKYSUB', 'Skysub failed for this slit'),
            ('BADEXTRACT', 'Extraction failed for this slit'),
        ])
        super().__init__(list(mask.keys()), descr=list(mask.values()))

    @property
    def exclude_for_reducing(self):
        # Ignore these flags when reducing or considering reduced slits
        return ['SKIPFLATCALIB']

    @property
    def exclude_for_flexure(self):
        # Ignore these flags when performing a flexure calculation
        #  Currently they are *all* of the flags..
        return ['SHORTSLIT', 'USERIGNORE', 'BADWVCALIB', 'BADTILTCALIB',
                'SKIPFLATCALIB', 'BADFLATCALIB', 'BADSKYSUB', 'BADEXTRACT']



class SlitTraceSet(calibframe.CalibFrame):
    """
    Defines a generic class for holding and manipulating image traces
    organized into left-right slit pairs.

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_slittraceset.rst

    Attributes:
        left_flexure (`numpy.ndarray`_):  Convenient spot to hold flexure corrected left
        right_flexure (`numpy.ndarray`_):  Convenient spot to hold flexure corrected right
    """
    calib_type = 'Slits'
    """Name for type of calibration frame."""

    calib_file_format = 'fits.gz'
    """File format for the calibration frame file."""

    version = '1.1.4'
    """SlitTraceSet data model version."""

    bitmask = SlitTraceBitMask()
    """
    Bit interpreter for slit masks.
    """

    internals = calibframe.CalibFrame.internals + ['left_flexure', 'right_flexure']
    """
    Attributes kept separate from the datamodel.
    """

    # Define the data model
    datamodel = {'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 'pypeline': dict(otype=str, descr='PypeIt pypeline name'),
                 'detname': dict(otype=str, descr='Identifier for detector or mosaic'),
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
                 'maskdef_designtab': dict(otype=Table, descr='Table with slitmask design and object info'),
                 'maskfile': dict(otype=str, descr='Data file that yielded the slitmask info'),
                 'maskdef_posx_pa': dict(otype=float, descr='PA that aligns with spatial dimension of the detector'),
                 'maskdef_offset': dict(otype=float, descr='Slitmask offset (pixels) from position expected '
                                                           'by the slitmask design'),
                 'maskdef_objpos': dict(otype=np.ndarray, atype=np.floating,
                                         descr='Object positions expected by the slitmask design [relative pixels]'),
                 'maskdef_slitcen': dict(otype=np.ndarray, atype=np.floating,
                                         descr='Slit centers expected by the slitmask design'),
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
                'slitbitm': dict(otype=str, descr='List of BITMASK keys from SlitTraceBitMask'),
                'specmin': dict(otype=np.ndarray, atype=np.floating,
                                descr='Minimum spectral position (pixel units) allowed for each slit/order.  '
                                      'Shape is Nslits.'),
                'specmax': dict(otype=np.ndarray, atype=np.floating,
                                descr='Maximum spectral position (pixel units) allowed for each slit/order.  '
                                      'Shape is Nslits.')}
    """Provides the class data model."""

    # TODO: Allow tweaked edges to be arguments?
    # TODO: May want nspat to be a required argument.
    # The INIT must contain every datamodel item or risk fail on I/O when it is a nested container
    def __init__(self, left_init, right_init, pypeline, detname=None, nspec=None, nspat=None,
                 PYP_SPEC=None, mask_init=None, specmin=None, specmax=None, binspec=1, binspat=1,
                 pad=0, spat_id=None, maskdef_id=None, maskdef_designtab=None, maskfile=None,
                 maskdef_posx_pa=None, maskdef_offset=None, maskdef_objpos=None,
                 maskdef_slitcen=None, ech_order=None, nslits=None, left_tweak=None,
                 right_tweak=None, center=None, mask=None, slitbitm=None):

        # Instantiate the DataContainer
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = dict([(k,values[k]) for k in args[1:]])
        # The dictionary passed to DataContainer.__init__ does not
        # contain self.
        # TODO: Does it matter if the calling function passes the
        # keyword arguments in a different order? No.
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
            # Validate -- All of the keys must be present and in current order, but new ones can exist
            bitms = self.slitbitm.split(',')
            curbitm = list(self.bitmask.keys())
            for kk, bit in enumerate(bitms):
                if curbitm[kk] != bit:
                    msgs.error("Input BITMASK keys differ from current data model!")
            # Update to current, no matter what
            self.slitbitm = ','.join(list(self.bitmask.keys()))
        # Mask
        if self.mask is None:
            self.mask = self.mask_init.copy()

    def _bundle(self):
        """
        Bundle the data in preparation for writing to a fits file.

        See :func:`pypeit.datamodel.DataContainer._bundle`. Data is
        always written to a 'SLITS' extension.
        """
        bndl = super()._bundle(ext='SLITS', transpose_arrays=True)
        if self.maskdef_designtab is not None:
            # save the table
            tab_detached = bndl[0]['SLITS']['maskdef_designtab']
            # remove `tab_detached` from the dict
            bndl[0]['SLITS'].pop('maskdef_designtab')
            # create a dict for the `tab_detached`
            tab_dict = {'maskdef_designtab': tab_detached}
            return [bndl[0], tab_dict]
        return bndl

    # TODO: Although I don't like doing it, kwargs is here to catch the
    # extraneous keywords that can be passed to _parse from the base class but
    # won't be used.
    @classmethod
    def _parse(cls, hdu, hdu_prefix=None, **kwargs):
        """
        Parse the data that was previously written to a fits file.

        See :func:`pypeit.datamodel.DataContainer._parse`. Data is
        always read from the 'SLITS' extension.
        """
        if not hasattr(hdu, '__len__'):
            return super()._parse(hdu, transpose_table_arrays=True)

        # TODO: My edit to the code causes the code below to fault in some cases
        # because of a consistency limitation that I put on the values that ext
        # could have.  The if statement above fixes the issue.  But I think the
        # code in the except block will always fault now, and I don't remember
        # why we needed this try/except block in the first place.
        try:
            return super()._parse(hdu, ext=['SLITS', 'MASKDEF_DESIGNTAB'],
                                  transpose_table_arrays=True)
        except KeyError:
            return super()._parse(hdu, ext='SLITS', transpose_table_arrays=True)

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
        THIS NEEDS A DESCRIPTION

        Returns:
            `numpy.ndarray`_:

        """
        #
        info = np.vstack([self.spat_id, self.mask])
        info = np.vstack([info, np.zeros_like(self.spat_id)]) \
                    if self.maskdef_id is None else np.vstack([info, self.maskdef_id])
        return info.astype(int).T

    @property
    def slitord_id(self):
        """
        Return array of slit_spatId (MultiSlit, IFU) or ech_order (Echelle) values

        Returns:
            `numpy.ndarray`_:

        """
        if self.pypeline in ['MultiSlit', 'IFU']:
            return self.spat_id
        if self.pypeline == 'Echelle':
            return self.ech_order
        msgs.error(f'Unrecognized Pypeline {self.pypeline}')

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

    def get_radec_image(self, wcs, alignSplines, tilts, initial=True, flexure=None):
        """Generate an RA and DEC image for every pixel in the frame
        NOTE: This function is currently only used for IFU reductions.

        Parameters
        ----------
        wcs : astropy.wcs
            The World Coordinate system of a science frame
        alignSplines : :class:`pypeit.alignframe.AlignmentSplines`
            An instance of the AlignmentSplines class that allows one to build and
            transform between detector pixel coordinates and WCS pixel coordinates.
        tilts : `numpy.ndarray`
            Spectral tilts.
        initial : bool
            Select the initial slit edges?
        flexure : float, optional
            If provided, offset each slit by this amount.

        Returns
        -------
        tuple : There are three elements in the tuple. The first two are 2D numpy arrays
                of shape (nspec, nspat), where the first ndarray is the RA image, and the
                second ndarray is the DEC image. RA and DEC are in units degrees. The third
                element of the tuple stores the minimum and maximum difference (in pixels)
                between the WCS reference (usually the centre of the slit) and the edges of
                the slits. The third array has a shape of (nslits, 2).
        """
        # Initialise the output
        raimg = np.zeros((self.nspec, self.nspat))
        decimg = np.zeros((self.nspec, self.nspat))
        minmax = np.zeros((self.nslits, 2))
        # Get the slit information
        slitid_img_init = self.slit_img(pad=0, initial=initial, flexure=flexure)
        for slit_idx, spatid in enumerate(self.spat_id):
            onslit = (slitid_img_init == spatid)
            onslit_init = np.where(onslit)
            if self.mask[slit_idx] != 0:
                msgs.error(f"Slit {spatid} ({slit_idx+1}/{self.spat_id.size}) is masked. Cannot generate RA/DEC image.")
            # Retrieve the pixel offset from the central trace
            evalpos = alignSplines.transform(slit_idx, onslit_init[1], onslit_init[0])
            minmax[slit_idx, 0] = np.min(evalpos)
            minmax[slit_idx, 1] = np.max(evalpos)
            # Calculate the WCS from the pixel positions
            slitID = np.ones(evalpos.size) * slit_idx - wcs.wcs.crpix[0]
            world_ra, world_dec, _ = wcs.wcs_pix2world(slitID, evalpos, tilts[onslit_init]*(self.nspec-1), 0)
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

    def slit_img(self, pad=None, slitidx=None, initial=False, 
                 flexure=None,
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
            flexure (:obj:`float`, optional):
                If provided, offset each slit by this amount
                Done in select_edges()

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

    def mask_add_missing_obj(self, sobjs, spat_flexure, fwhm, boxcar_rad):
        """
        Generate new SpecObj and add them into the SpecObjs object for any slits missing the targeted source.

        Args:
            sobjs (:class:`pypeit.specobjs.SpecObjs`):
                List of SpecObj that have been found and traced
            spat_flexure (:obj:`float`):
                Shifts, in spatial pixels, between this image and SlitTrace
            fwhm (:obj:`float`):
                FWHM in pixels to be used in the optimal extraction
            boxcar_rad (:obj:`float`):
                BOX_RADIUS in pixels to be used in the boxcar extraction

        Returns:
            :class:`pypeit.specobjs.SpecObjs`: Updated list of SpecObj that have been found and traced

        """
        msgs.info('Add undetected objects at the expected location from slitmask design.')

        if fwhm is None:
            msgs.error('A FWHM for the optimal extraction must be provided. See `find_fwhm` in `FindObjPar`.')

        if self.maskdef_objpos is None:
            msgs.error('An array with the object positions expected from slitmask design is missing.')

        if self.maskdef_offset is None:
            msgs.error('A value for the slitmask offset must be provided.')

        # Restrict to objects on this detector
        if sobjs.nobj > 0:
            on_det = (sobjs.DET == self.detname) & (sobjs.OBJID > 0) # use only positive detections
            cut_sobjs = sobjs[on_det]
        else:
            cut_sobjs = sobjs

        # get slits edges init
        left_init, _, _ = self.select_edges(initial=True, flexure=spat_flexure)  # includes flexure
        # get slits edges tweaked
        left_tweak, right_tweak, _ = self.select_edges(initial=False, flexure=spat_flexure)  # includes flexure

        # midpoint in the spectral direction
        specmid = left_init[:,0].size//2

        # Loop on all the good slits
        gd_slit = self.mask == 0
        for islit in np.where(gd_slit)[0]:
            # Check for assigned obj
            if len(cut_sobjs) > 0:
                target_obj_in_slit = (cut_sobjs.MASKDEF_ID == self.maskdef_id[islit]) & \
                                     (cut_sobjs.MASKDEF_OBJNAME != 'SERENDIP')
                if np.any(target_obj_in_slit):
                    continue

            # Object index
            oidx = np.where(self.maskdef_designtab['MASKDEF_ID'].data == self.maskdef_id[islit])[0][0]
            #
            # Do it
            SPAT_PIXPOS = float((self.maskdef_objpos[oidx]+self.maskdef_offset) + left_init[specmid, islit])

            # TODO -- DO THIS RIGHT FOR SLITS (LIKELY) OFF THE DETECTOR/SLIT
            #  If we keep what follows, probably should add some tolerance to be off the edge
            #  Otherwise things break in skysub
            if (SPAT_PIXPOS > right_tweak[specmid, islit]) or (SPAT_PIXPOS < left_tweak[specmid, islit]):
                msgs.warn("Targeted object is off the detector")
                continue

            # Generate a new specobj
            specobj_dict = {'SLITID': self.spat_id[islit], # Confirm this
                            'DET': self.detname,
                            'OBJTYPE': 'science',  # Confirm this is ok
                            'PYPELINE': self.pypeline}
            thisobj = specobj.SpecObj(**specobj_dict)

            # Fill me up
            if cut_sobjs.nobj == 0:
                # This uses slit edge trace
                xoff = SPAT_PIXPOS - self.center[specmid, islit]
                thisobj.TRACE_SPAT = self.center[:, islit] + xoff
                thisobj.SPAT_PIXPOS = SPAT_PIXPOS
                thisobj.SPAT_FRACPOS = (SPAT_PIXPOS - left_tweak[specmid, islit]) / \
                                       (right_tweak[specmid, islit]-left_tweak[specmid, islit])
                thisobj.trace_spec = np.arange(left_tweak.shape[0])
                # OBJID
                thisobj.OBJID = 1
            else:
                # This uses an object
                idx_nearest = np.argmin(np.abs(cut_sobjs.SPAT_PIXPOS-SPAT_PIXPOS))
                specmid = cut_sobjs[idx_nearest].TRACE_SPAT.size//2
                xoff = SPAT_PIXPOS - cut_sobjs[idx_nearest].TRACE_SPAT[specmid]
                thisobj.TRACE_SPAT = cut_sobjs[idx_nearest].TRACE_SPAT + xoff
                thisobj.SPAT_PIXPOS = SPAT_PIXPOS
                thisobj.trace_spec = cut_sobjs[idx_nearest].trace_spec
                thisobj.SPAT_FRACPOS = (SPAT_PIXPOS - left_tweak[specmid, islit]) / \
                                       (right_tweak[specmid, islit]-left_tweak[specmid, islit])
                # OBJID
                any_obj_in_slit = cut_sobjs.MASKDEF_ID == self.maskdef_id[islit]
                if np.any(any_obj_in_slit):
                    thisobj.OBJID = np.max(cut_sobjs[any_obj_in_slit].OBJID) + 1
                else:
                    thisobj.OBJID = 1

            # FWHM
            thisobj.FWHM = fwhm  # pixels
            thisobj.BOX_RADIUS = boxcar_rad  # pixels
            thisobj.maskwidth = 4. * fwhm  # matches objfind() in extract.py
            thisobj.smash_snr = 0.
            thisobj.smash_peakflux = 0.
            #thisobj.THRESHOLD = 0.
            # Finishing up
            thisobj.set_name()
            # Mask info
            thisobj.RA = self.maskdef_designtab['OBJRA'][oidx]
            thisobj.DEC = self.maskdef_designtab['OBJDEC'][oidx]
            thisobj.MASKDEF_OBJNAME = self.maskdef_designtab['OBJNAME'][oidx]
            thisobj.MASKDEF_OBJMAG = self.maskdef_designtab['OBJMAG'][oidx]
            thisobj.MASKDEF_OBJMAG_BAND = self.maskdef_designtab['OBJMAG_BAND'][oidx]
            thisobj.MASKDEF_ID = self.maskdef_designtab['MASKDEF_ID'][oidx]
            thisobj.MASKDEF_EXTRACT = True
            thisobj.hand_extract_flag = False
            # Add to SpecObjs
            sobjs.add_sobj(thisobj)

        if sobjs.nobj > 0:
            # Sort objects according to their spatial location
            spat_pixpos = sobjs.SPAT_PIXPOS
            sobjs = sobjs[spat_pixpos.argsort()]

            # Vette
            for sobj in sobjs:
                if not sobj.ready_for_extraction():
                    msgs.error("Bad SpecObj.  Can't proceed")

        # Return
        return sobjs

    def assign_maskinfo(self, sobjs, plate_scale, spat_flexure, TOLER=1.):
        """
        Assign RA, DEC, Name to objects.
        Modified in place.

        Args:
            sobjs (:class:`pypeit.specobjs.SpecObjs`): List of SpecObj that have been found and traced.
            plate_scale (:obj:`float`): platescale for the current detector.
            spat_flexure (:obj:`float`): Shifts, in spatial pixels, between this image and SlitTrace.
            det_buffer (:obj:`int`): Minimum separation between detector edges and a slit edge.
            TOLER (:obj:`float`, optional): Matching tolerance in arcsec.

        Returns:
            :class:`pypeit.specobjs.SpecObjs`: Updated list of SpecObj that have been found and traced.

        Returns:
            :class:`pypeit.specobjs.SpecObjs`: Updated list of SpecObj that have been found and traced

        """

        if self.maskdef_objpos is None:
            msgs.error('An array of object positions predicted by the slitmask design must be provided.')
        if self.maskdef_slitcen is None:
            msgs.error('An array of slit centers predicted by the slitmask design must be provided.')
        if self.maskdef_offset is None:
            msgs.error('A value for the slitmask offset must be provided.')

        # Unpack -- Remove this once we have a DataModel
        obj_maskdef_id = self.maskdef_designtab['MASKDEF_ID'].data
        obj_slit_coords = SkyCoord(ra=self.maskdef_designtab['SLITRA'],
                                   dec=self.maskdef_designtab['SLITDEC'], frame='fk5', unit='deg')
        obj_slit_pa = self.maskdef_designtab['SLITPA']
        # RMS (in arcsec) of the x-correlation between predicted and traced left slits edges
        cc_rms = self.maskdef_designtab.meta['MASKRMSL'] * plate_scale

        # Restrict to objects on this detector
        if sobjs.nobj > 0:
            on_det = (sobjs.DET == self.detname) & (sobjs.OBJID > 0) # use only positive detections
            cut_sobjs = sobjs[on_det]
            if cut_sobjs.nobj == 0:
                msgs.warn('NO detected objects.')
                return sobjs
        else:
            msgs.warn('NO detected objects.')
            return sobjs

        msgs.info('Assign slitmask design info to detected objects. '
                  'Matching tolerance includes user-provided tolerance, slit tracing uncertainties and object size.')

        # get slits edges init
        left_init, right_init, _ = self.select_edges(initial=True, flexure=spat_flexure)  # includes flexure

        # midpoint in the spectral direction
        specmid = left_init[:, 0].size // 2

        # First pass
        measured, expected = [], []
        for sobj in cut_sobjs:
            # Set MASKDEF_ID
            sobj.MASKDEF_ID = self.maskdef_id[self.spat_id == sobj.SLITID][0]
            # object ID
            # TODO -- Add to SpecObj DataModel?
            # There is small chance that self.maskdef_id=-99. This would definitely happen if the user
            # add a custom slit. If maskdef_id=-99 for a certain object, we cannot assign OBJECT, RA, DEC
            oidx = np.where(obj_maskdef_id == sobj.MASKDEF_ID)[0]
            if oidx.size == 0:
                # In this case I have to give a fake value for index reasons later in the code.
                measured.append(-9999.9)
                expected.append(9999.9)
            else:
                oidx = oidx[0]
                # Expected object position (distance from left edge) in pixels, corrected for edge loss
                expected_objpos = self.maskdef_objpos[oidx]
                # Measured object position (distance from left edge) in pixels
                measured_objpos = sobj.SPAT_PIXPOS - left_init[specmid, self.maskdef_id == sobj.MASKDEF_ID][0]
                # Finish
                measured.append(measured_objpos)
                expected.append(expected_objpos)
        measured = np.array(measured)
        expected = np.array(expected)

        # Assign
        # Loop on slits to deal with multiple sources within TOLER
        # Exclude the objects that have maskdef_id=-99
        uni_maskid = np.unique(cut_sobjs.MASKDEF_ID[cut_sobjs.MASKDEF_ID!=-99])
        for maskid in uni_maskid:
            # Index for SpecObjs on this slit
            idx = np.where(cut_sobjs.MASKDEF_ID == maskid)[0]
            # Index for slitmask
            sidx = np.where(self.maskdef_designtab['MASKDEF_ID'] == maskid)[0][0]
            # Within TOLER?
            # separation in pixels
            separ = measured[idx] - (expected[idx] + self.maskdef_offset)
            msgs.info('MASKDEF_ID:{}'.format(maskid))
            msgs.info('Difference between expected and detected object '
                      'positions: {} arcsec'.format(np.round(separ*plate_scale, 2)))
            # we include in the tolerance the rms of the slit edges matching and the size of
            # the detected object with the highest peak flux
            if np.any((cut_sobjs[idx].smash_peakflux != None) & (cut_sobjs[idx].smash_peakflux != 0.)):
                ipeak = np.argmax(cut_sobjs[idx].smash_peakflux)
                obj_fwhm = cut_sobjs[ipeak].FWHM*plate_scale
            else:
                obj_fwhm = 0.
            in_toler = np.abs(separ*plate_scale) < (TOLER + cc_rms + obj_fwhm/2)
            if np.any(in_toler):
                # Find positive peakflux
                peak_flux = cut_sobjs[idx].smash_peakflux[in_toler]
                pos_peak_flux = np.where(peak_flux>0)[0]
                if np.any(pos_peak_flux):
                    # find the object with the smallest separation
                    closest_idx = np.argmin(np.abs(separ[in_toler][pos_peak_flux]))
                else:
                    closest_idx = np.argmin(np.abs(separ[in_toler]))
                imx_idx = idx[in_toler][closest_idx]
                # Object in Mask Definition
                oidx = np.where(obj_maskdef_id == maskid)[0][0]
                # Assign
                sobj = cut_sobjs[imx_idx]
                sobj.RA = self.maskdef_designtab['OBJRA'][oidx]
                sobj.DEC = self.maskdef_designtab['OBJDEC'][oidx]
                sobj.MASKDEF_OBJNAME = self.maskdef_designtab['OBJNAME'][oidx]
                sobj.MASKDEF_OBJMAG = self.maskdef_designtab['OBJMAG'][oidx]
                sobj.MASKDEF_OBJMAG_BAND = self.maskdef_designtab['OBJMAG_BAND'][oidx]
                sobj.MASKDEF_EXTRACT = False
                # Remove that idx value
                idx = idx.tolist()
                idx.remove(imx_idx)
                idx = np.array(idx)
            # Fill in the rest
            for ss in idx:
                sobj = cut_sobjs[ss]
                # Measured coordinates
                offset = measured[ss] - (self.maskdef_slitcen + self.maskdef_offset - left_init)[specmid, self.spat_id == sobj.SLITID][0]
                new_obj_coord = obj_slit_coords[sidx].directional_offset_by(
                    np.radians(obj_slit_pa[sidx]), (offset*plate_scale)*units.arcsec)
                # Assign
                sobj.RA = new_obj_coord.ra.value
                sobj.DEC = new_obj_coord.dec.value
                sobj.MASKDEF_OBJNAME = 'SERENDIP'
                sobj.MASKDEF_EXTRACT = False
        # Give fake values of RA, DEC, and MASKDEF_OBJNAME for object with maskdef_id=-99.
        noidx = np.where(cut_sobjs.MASKDEF_ID == -99)[0]
        if noidx.size > 0:
            for sobj in cut_sobjs[noidx]:
                # Assign
                sobj.RA = 0.0
                sobj.DEC = 0.0
                sobj.MASKDEF_OBJNAME = 'NONE'
                sobj.MASKDEF_EXTRACT = False

        # Return
        return sobjs

    def det_of_slit(self, spat_id:int, det_img:np.ndarray,
                    slit_img:np.ndarray=None):
        """ Identify the 'best' detector for this slit/order
        The metric is the detector where the slit appears
        the most frequently.

        Only sensibly used for mosaic images

        Args:
            spat_id (`int`): 
                spat_id value for the slit of interest
            det_img (`numpy.ndarray`_): 
                :obj:`int` image specifying the detector number
                (1-based) for each pixel in the mosaic
            slit_img (`numpy.ndarray`_, optional): 
                image identifying each pixel with its associated slit.

        Returns:
            :obj:`int`: Detector number for the slit (1-based)
        """
        # Grab slit image?
        if slit_img is None:
            slit_img = self.slit_img()
        # Find the most common detector value
        det, cnt = np.unique(
            det_img[(slit_img == spat_id) & (det_img > 0)], 
            return_counts=True)
        return det[np.argmax(cnt)]

    def get_maskdef_objpos(self, plate_scale, det_buffer):
        """
        Determine the object positions expected by the slitmask design

        Args:
            plate_scale (:obj:`float`): platescale for the current detector
            det_buffer (:obj:`int`): Minimum separation between detector edges and a slit edge

        """

        # midpoint in the spectral direction
        specmid = self.left_init[:,0].size//2

        # Unpack -- Remove this once we have a DataModel
        obj_maskdef_id = self.maskdef_designtab['MASKDEF_ID'].data
        # Distance (arcsec) of the object from the left edge
        obj_topdist = self.maskdef_designtab['OBJ_TOPDIST'].data
        obj_botdist = self.maskdef_designtab['OBJ_BOTDIST'].data

        # slit lengths
        expected_slitlen = (obj_topdist + obj_botdist) / plate_scale  # binned pixels
        measured_slitlen = self.right_init[specmid, :] - self.left_init[specmid, :]  # binned pixels
        # difference between measured and expected slit length (but only for the left side).
        left_edgeloss = np.zeros(self.maskdef_id.size)
        # define a new slit center to take into account the slits that are smaller than what
        # should be because they are cut by the detector edges.
        new_slitcen = self.center.copy()
        for i, maskid in enumerate(self.maskdef_id):
            if maskid == -99:
                left_edgeloss[i] = -9999.9
            else:
                if (self.left_init[specmid, i] != det_buffer) & (self.right_init[specmid, i] != self.nspat - det_buffer):
                    # for all but the slits that fall outside the detector, we assume the edge loss happens
                    # equally in both sides of the slit
                    left_edgeloss[i] = (expected_slitlen[obj_maskdef_id == maskid][0] - measured_slitlen[i]) / 2.  # pixels
                    # because edge loss happens equally in both sides, we assume no shift in slit center

        # Stats (typical edge loss across the detector)
        if left_edgeloss[(left_edgeloss != 0)&(left_edgeloss != -9999.9)].size > 3:
            _, median_edgeloss, _ = sigma_clipped_stats(left_edgeloss[(left_edgeloss != 0)&(left_edgeloss != -9999.9)],
                                                        sigma=2.)
        else:
            median_edgeloss = 0.

        expected_objpos_all = np.zeros(obj_maskdef_id.size)
        for i, maskid in enumerate(obj_maskdef_id):
            # fix the edge loss for slits that fall off the detector
            if self.left_init[specmid, i] == det_buffer:
                # for the leftmost slit that is cut by the detector edge, we assume that the edge loss
                # happens mostly in the left side and we assume a value (median_edgeloss) for the right side
                left_edgeloss[i] = expected_slitlen[obj_maskdef_id == maskid][0] - measured_slitlen[i] \
                                   - median_edgeloss  # pixels
                # shift in slit center due to slit partially falling outside the detector
                censhift = left_edgeloss / 2.
                new_slitcen[i] -= censhift
            if self.right_init[specmid, i] == self.nspat - det_buffer:
                # for the rightmost slit that is cut by the detector edge, we assume that the edge loss
                # happens mostly in the right side, and `left_edgeloss` is assumed to be equal to `median_edgeloss`.
                left_edgeloss[i] = median_edgeloss
                # shift in slit center due to slit partially falling outside the detector
                censhift = (expected_slitlen[obj_maskdef_id == maskid][0] - measured_slitlen[i] - median_edgeloss) / 2.
                new_slitcen[i] += censhift

            # Expected objects position (distance from left edge) in pixels, corrected for edge loss
            expected_objpos_all[i] = obj_topdist[i] / plate_scale - left_edgeloss[self.maskdef_id == maskid][0]

        self.maskdef_objpos = expected_objpos_all
        self.maskdef_slitcen = new_slitcen

        return

    def get_maskdef_offset(self, sobjs, platescale, spat_flexure, slitmask_off, bright_maskdefid,
                           snr_thrshd, use_alignbox, dither_off=None):
        """
        Determine the Slitmask offset (pixels) from position expected by the slitmask design

        Args:
            sobjs (:class:`pypeit.specobjs.SpecObjs`): List of SpecObj that have been found and traced
            platescale (:obj:`float`): Platescale
            spat_flexure (:obj:`float`): Shifts, in spatial pixels, between this image and SlitTrace
            slitmask_off (:obj:`float`): User provided slitmask offset in pixels
            bright_maskdefid (:obj:`str`): User provided maskdef_id of a bright object to be used to measure offset
            snr_thrshd (:obj:`float`): Objects detected above this S/N ratio threshold will be use to
                                        compute the slitmask offset
            use_alignbox (:obj:`bool`): Flag that determines if the alignment boxes are used to measure the offset
            dither_off (:obj:`float`, optional): dither offset recorded in the header of the observations


        """
        if self.maskdef_objpos is None:
            msgs.error('An array of object positions predicted by the slitmask design must be provided.')
        if self.maskdef_slitcen is None:
            msgs.error('An array of slit centers predicted by the slitmask design must be provided.')

        # If slitmask offset provided by the user, just save it and return
        if slitmask_off is not None:
            self.maskdef_offset = slitmask_off
            msgs.info('User-provided slitmask offset: {} pixels ({} arcsec)'.format(round(self.maskdef_offset, 2),
                                                                            round(self.maskdef_offset*platescale, 2)))
            return
        # If using the dither offeset recorde in the header, just save it and return
        if dither_off is not None:
            self.maskdef_offset = -dither_off/platescale
            msgs.info('Slitmask offset from the dither pattern: {} pixels ({} arcsec)'.
                      format(round(self.maskdef_offset, 2), round(self.maskdef_offset*platescale, 2)))
            return

        # Restrict to objects on this detector
        if sobjs.nobj > 0:
            on_det = (sobjs.DET == self.detname) & (sobjs.OBJID > 0) # use only positive detections
            cut_sobjs = sobjs[on_det]
            if cut_sobjs.nobj == 0:
                msgs.warn('NO detected objects. Slitmask offset cannot be estimated in '
                          f'{self.detname}.')
                self.maskdef_offset = 0.0
                return
        else:
            msgs.warn('NO detected objects. Slitmask offset cannot be estimated in '
                      f'{self.detname}.')
            self.maskdef_offset = 0.0
            return

        # Maskdef ID
        obj_maskdef_id = self.maskdef_designtab['MASKDEF_ID'].data
        # Flag for slits used for alignment (1-yes; 0-no)
        flag_align = self.maskdef_designtab['ALIGN'].data
        align_maskdef_ids = obj_maskdef_id[flag_align == 1]

        # get slits edges init
        left_init, _, _ = self.select_edges(initial=True, flexure=spat_flexure)  # includes flexure
        # midpoint in the spectral direction
        specmid = left_init[:, 0].size // 2

        # First pass
        measured, expected = [], []
        for sobj in cut_sobjs:
            # Set MASKDEF_ID
            sobj.MASKDEF_ID = self.maskdef_id[self.spat_id == sobj.SLITID][0]
            # object ID
            # TODO -- Add to SpecObj DataModel?
            # There is small chance that self.maskdef_id=-99. This would definitely happen if the user
            # add a custom slit. If maskdef_id=-99 for a certain object, we cannot assign OBJECT, RA, DEC
            oidx = np.where(obj_maskdef_id == sobj.MASKDEF_ID)[0]
            if oidx.size == 0:
                # In this case I have to give a fake value for index reasons later in the code.
                measured.append(-9999.9)
                expected.append(9999.9)
            else:
                oidx = oidx[0]
                # Expected object position (distance from left edge) in pixels, corrected for edge loss
                expected_objpos = self.maskdef_objpos[oidx]
                # Measured object position (distance from left edge) in pixels
                measured_objpos = sobj.SPAT_PIXPOS - left_init[specmid, self.maskdef_id == sobj.MASKDEF_ID][0]
                # Finish
                measured.append(measured_objpos)
                expected.append(expected_objpos)
        measured = np.array(measured)
        expected = np.array(expected)

        # If the users want to use the alignment boxes to trace the offset they will set the flag `use_alignbox` to True
        if use_alignbox:
            align_offs = []
            for align_id in align_maskdef_ids:
                sidx = np.where(cut_sobjs.MASKDEF_ID == align_id)[0]
                if sidx.size > 0:
                    # Take the brightest source as the star
                    peak_flux = cut_sobjs[sidx].smash_peakflux
                    imx_peak = np.argmax(peak_flux)
                    imx_sidx = sidx[imx_peak]
                    bright_measured = measured[imx_sidx]
                    bright_expected = expected[imx_sidx]
                    align_offs.append(bright_measured - bright_expected)
            align_offs = np.array(align_offs)
            if align_offs.size > 0:
                mean, median_off, std = sigma_clipped_stats(align_offs, sigma=2.)
                self.maskdef_offset = median_off
                msgs.info(f'Slitmask offset estimated using ALIGN BOXES in {self.detname}: '
                          f'{round(self.maskdef_offset, 2)} pixels ('
                          f'{round(self.maskdef_offset*platescale, 2)} arcsec).')
            else:
                self.maskdef_offset = 0.0
                msgs.info('NO objects detected in ALIGN BOXES. Slitmask offset '
                          f'cannot be estimated in {self.detname}.')
            return

        # if the maskdef_id of a bright object is provided by the user, check if it is in
        # this detector and use it to compute the offset
        if bright_maskdefid is not None:
            if bright_maskdefid in obj_maskdef_id:
                sidx = np.where(cut_sobjs.MASKDEF_ID == bright_maskdefid)[0]
                if sidx.size == 0:
                    self.maskdef_offset = 0.0
                    msgs.info(f'Object in slit {bright_maskdefid} not detected. Slitmask offset '
                              f'cannot be estimated in {self.detname}.')
                else:
                    # Parse the peak fluxes
                    peak_flux = cut_sobjs[sidx].smash_peakflux
                    imx_peak = np.argmax(peak_flux)
                    imx_sidx = sidx[imx_peak]
                    bright_measured = measured[imx_sidx]
                    bright_expected = expected[imx_sidx]
                    self.maskdef_offset = bright_measured - bright_expected
                    msgs.info('Slitmask offset computed using bright object in slit '
                              f'{bright_maskdefid} ({self.detname}): '
                              f'{round(self.maskdef_offset, 2)} pixels ('
                              f'{round(self.maskdef_offset*platescale, 2)} arcsec)')
            else:
                self.maskdef_offset = 0.0
            return

        # Determine offsets using only detections with the highest signal
        # objects added in manual extraction have smash_snr = None
        nonone = cut_sobjs.smash_snr != None
        if len(cut_sobjs[nonone]) > 0:
            highsnr_measured = measured[nonone][cut_sobjs[nonone].smash_snr > snr_thrshd]
            highsnr_expected = expected[nonone][cut_sobjs[nonone].smash_snr > snr_thrshd]
            if len(highsnr_measured) >= 3:
                off = highsnr_measured - highsnr_expected
                mean, median_off, std = sigma_clipped_stats(off, sigma=2.)
                self.maskdef_offset = median_off
                msgs.info(f'Slitmask offset estimated in {self.detname}: '
                          f'{round(self.maskdef_offset, 2)} pixels ('
                          f'{round(self.maskdef_offset*platescale, 2)} arcsec)')
            else:
                msgs.warn(f'Less than 3 objects detected above {snr_thrshd} sigma threshold. '
                          f'Slitmask offset cannot be estimated in {self.detname}.')
                self.maskdef_offset = 0.0
        else:
            msgs.warn(f'Less than 3 objects detected above {snr_thrshd} sigma threshold. '
                      f'Slitmask offset cannot be estimated in {self.detname}.')
            self.maskdef_offset = 0.0

        return

    def get_maskdef_extract_fwhm(self, sobjs, platescale, fwhm_parset, find_fwhm):
        """
        This method determines the fwhm to use for the optimal extraction
        of maskdef_extract (i.e., undetected) objects.
        If the user provides a fwhm, it would be used. Otherwise fwhm
        will be computed using the average fwhm of the detected objects.

        Args:
            sobjs (:class:`pypeit.specobjs.SpecObjs`):
                List of SpecObj that have been found and traced.
            platescale (:obj:`float`):
                Platescale.
            fwhm_parset (:obj:`float`, optional):
                Parset that guides the determination of the fwhm of the maskdef_extract objects.
                If None (default) the fwhm are computed as the averaged from the detected objects,
                if it is a number it will be adopted as the fwhm.
            find_fwhm (:obj:`float`):
            Initial guess of the objects fwhm in pixels (used in object finding)

        Returns:
            :obj:`float`: FWHM in pixels to be used in the optimal extraction

        """
        msgs.info('Determining the FWHM to be used for the optimal extraction of `maskdef_extract` objects')
        fwhm = None
        if fwhm_parset is not None:
            msgs.info(f'Using user-provided FWHM = {fwhm_parset}"')
            fwhm = fwhm_parset/platescale
        elif sobjs.nobj > 0:
            # Use average FWHM of detected objects, but remove the objects in the alignment boxes
            # Find align boxes
            maskdef_id = self.maskdef_designtab['MASKDEF_ID'].data
            # Flag to identify alignment boxes (1-yes; 0-no)
            flag_align = self.maskdef_designtab['ALIGN'].data
            align_maskdef_ids = maskdef_id[flag_align == 1]
            all_fwhm = np.array([])
            for ss in sobjs:
                # append only the FWHM of objects that were detected not in the alignment boxes
                if ss.MASKDEF_ID not in align_maskdef_ids:
                    all_fwhm = np.append(all_fwhm, ss.FWHM)
            if all_fwhm.size > 0:
                # compute median
                _, fwhm, _ = sigma_clipped_stats(all_fwhm, sigma=2.)
                msgs.info('Using median FWHM = {:.3f}" from detected objects.'.format(fwhm*platescale))
        if fwhm is None:
            fwhm = find_fwhm
            msgs.warn('The median FWHM cannot be determined because no objects were detected. '
                      'Using `find_fwhm` = {:.3f}". if the user wants to provide a value '
                      'set parameter `missing_objs_fwhm` in `SlitMaskPar`'.format(fwhm*platescale))

        return fwhm

    def user_mask(self, det, user_slits):
        """
        Mask all but the input slit

        Args:
            det (:obj:`int`): Detector number
            user_slits (:obj:`dict`):
        """
        if user_slits['method'] == 'slitspat':
            # Parse
            dets, spat_ids = parse.parse_slitspatnum(
                user_slits['slit_info'])
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
            self.mask[msk] = self.bitmask.turn_on(self.mask[msk],
                                                  'USERIGNORE')
        elif user_slits['method'] == 'maskIDs':
            # Mask only the good one
            msk = np.logical_not(np.isin(self.maskdef_id, user_slits['slit_info']))
            # Set
            self.mask[msk] = self.bitmask.turn_on(self.mask[msk],
                                                  'USERIGNORE')
        else:
            msgs.error('Not ready for this method: {:s}'.format(
                user_slits['method']))

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
                # This condition is added to avoid to give a 'BADWVCALIB' flag to alignment boxes,
                # which are purposely not wavelength calibrated, but are used during find object
                if np.logical_not(self.bitmask.flagged(self.mask[islit], flag='BOXSLIT')):
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




def merge_user_slit(slitspatnum, maskIDs):
    # Not set?
    if slitspatnum is None and maskIDs is None:
        return None
    #
    if slitspatnum is not None and maskIDs is not None:
        msgs.error("These should not both have been set")
    # MaskIDs
    user_slit_dict = {}
    if maskIDs is not None:
        user_slit_dict['method'] = 'maskIDs'
        user_slit_dict['slit_info'] = maskIDs
    else:
        user_slit_dict['method'] = 'slitspat'
        user_slit_dict['slit_info'] = slitspatnum
    # Return
    return user_slit_dict


def get_maskdef_objpos_offset_alldets(sobjs, calib_slits, spat_flexure, platescale, det_buffer, slitmask_par,
                                      dither_off=None):
    """
    Loop around all the calibrated detectors to extract information on the object positions
    expected by the slitmask design and the offsets between the expected and measure slitmask position.
    This info is recorded in the `SlitTraceSet` datamodel.

    Args:
        sobjs (:class:`pypeit.specobjs.SpecObjs`): List of SpecObj that have been found and traced
        calib_slits (:obj:`list`): List of `SlitTraceSet` with information on the traced slit edges
        spat_flexure (:obj:`list`): List of shifts, in spatial pixels, between this image and SlitTrace
        platescale (:obj:`list`): List of platescale for every detector
        det_buffer (:obj:`int`): Minimum separation between detector edges and a slit edge
        slitmask_par (:class:`pypeit.par.pypeitpar.PypeItPar`): slitmask PypeIt parameters
        dither_off (:obj:`float`, optional): dither offset recorded in the header of the observations

    Returns:
        List of `SlitTraceSet` with updated information on the traced slit edges

    """

    # grab corresponding detectors
    calib_dets = np.array([ss.detname for ss in calib_slits])
    for i in range(calib_dets.size):
        if calib_slits[i].maskdef_designtab is not None:
            # get object positions expected by slitmask design
            calib_slits[i].get_maskdef_objpos(platescale[i], det_buffer)

            # get slitmask offset in each single detector
            calib_slits[i].get_maskdef_offset(sobjs, platescale[i], spat_flexure[i],
                                              slitmask_par['slitmask_offset'],
                                              slitmask_par['bright_maskdef_id'],
                                              slitmask_par['snr_thrshd'],
                                              slitmask_par['use_alignbox'],
                                              dither_off=dither_off)

    return calib_slits


def average_maskdef_offset(calib_slits, platescale, list_detectors):
    """
    Loop around all the calibrated detectors to compute the median offset between
    the expected and measure slitmask position. This info is recorded in the `SlitTraceSet` datamodel.

    Args:
        calib_slits (:obj:`list`):
            List of :class:`~pypeit.slittrace.SlitTraceSet` objects with
            information on the traced slit edges.
        platescale (:obj:`float`):
            Platescale, must be the same for every detector.
        list_detectors (`numpy.ndarray`_):
            An array that lists the detector numbers of the current
            spectrograph; see
            :func:`~pypeit.spectrographs.spectrograph.Spectrograph.list_detectors`.
            If there are multiple detectors along the dispersion direction,
            there are ordered along the first axis.  For example, all the
            "bluest" detectors would be in ``list_detectors[0]``.

    Returns:
        `numpy.ndarray`_: Array of :class:`~pypeit.slittrace.SlitTraceSet`
        objects with updated information on the traced slit edges.
    """

    calib_slits = np.array(calib_slits)
    if list_detectors is None:
        msgs.warn('No average slitmask offset computed')
        return calib_slits

    # unpack list_detectors
    blue_and_red = list_detectors.ndim > 1
    spectrograph_dets = list_detectors if blue_and_red else np.expand_dims(list_detectors, 0)

    # determine if a slitmask offset exist and use the average offset over all the detectors
    # grab slitmask offsets from slits calibrations
    slitmask_offsets = np.array([ss.maskdef_offset for ss in calib_slits])
    # grab corresponding detectors
    calib_dets = np.array([ss.detname for ss in calib_slits])

    # remove eventual None and zeros (zero is assigned when no offset could be measured.)
    calib_dets = calib_dets[(slitmask_offsets != None) & (slitmask_offsets != 0)]
    slitmask_offsets = slitmask_offsets[(slitmask_offsets != None) & (slitmask_offsets != 0)].astype('float')

    if slitmask_offsets.size == 0:
        # If all detectors have maskdef_offset=0 give a warning
        msgs.warn('No slitmask offset could be measured. Assumed to be zero. ')
        msgs.warn('RA, DEC, OBJNAME assignment and forced extraction of undetected objects MAY BE WRONG! '
                  'Especially for dithered observations!')
        msgs.warn('To provide a value set `slitmask_offset` in `SlitMaskPar`')

        return calib_slits

    # are there dets from calib_slits that are blue?
    indx_b = np.where(np.in1d(calib_dets, spectrograph_dets[0]))[0]
    # if this spectrograph is not split into blue and red detectors
    # or if it is but there are no available offsets in the blue
    if not blue_and_red or indx_b.size == 0:
        # use all the available offsets to compute the median
        _, median_off, _ = sigma_clipped_stats(slitmask_offsets, sigma=2.)
        for cs in calib_slits:
            # assign median to each det
            cs.maskdef_offset = median_off
        msgs.info('Average Slitmask offset: {:.2f} pixels ({:.2f} arcsec).'.format(median_off, median_off * platescale))

        return calib_slits

    if indx_b.size > 0:
        # compute median if these blue dets have values of slitmask_offsets
        _, median_off, _ = sigma_clipped_stats(slitmask_offsets[indx_b], sigma=2.)
        for cs in calib_slits:
            if cs.detname in spectrograph_dets[0]:
                # assign median to each blue det
                cs.maskdef_offset = median_off
        msgs.info('Average Slitmask offset for the blue detectors: '
                  '{:.2f} pixels ({:.2f} arcsec).'.format(median_off, median_off * platescale))

        # which dets from calib_slits are red?
        indx_r = np.where(np.in1d(calib_dets, spectrograph_dets[1]))[0]
        if indx_r.size > 0:
            # compute median if these red dets have values of slitmask_offsets
            _, median_off, _ = sigma_clipped_stats(slitmask_offsets[indx_r], sigma=2.)

        # assign median to each red det (median would be the one computed for red dets if exists
        # or the median computed for blue dets)
        for cs in calib_slits:
            if cs.detname in spectrograph_dets[1]:
                cs.maskdef_offset = median_off
        msgs.info('Average Slitmask offset for the red detectors: '
                  '{:.2f} pixels ({:.2f} arcsec).'.format(median_off, median_off * platescale))

    return calib_slits


def assign_addobjs_alldets(sobjs, calib_slits, spat_flexure, platescale, slitmask_par, find_fwhm):
    """
    Loop around all the calibrated detectors to assign RA, DEC and OBJNAME to
    extracted object and to force extraction of undetected objects.

    Args:
        sobjs (:class:`~pypeit.specobjs.SpecObjs`):
            List of SpecObj that have been found and traced.
        calib_slits (`numpy.ndarray`_):
            Array of `SlitTraceSet` with information on the traced slit edges.
        spat_flexure (:obj:`list`):
            List of shifts, in spatial pixels, between this image and SlitTrace.
        platescale (:obj:`list`):
            List of platescale for every detector.
        slitmask_par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            Slitmask PypeIt parameters.
        find_fwhm (:obj:`float`):
            Initial guess of the objects fwhm in pixels (used in object finding)

    Returns:
        :class:`~pypeit.specobjs.SpecObjs`:
            Updated list of spectra that have been found and traced.
    """

    # grab corresponding detectors
    calib_dets = np.array([ss.detname for ss in calib_slits])
    for i in range(calib_dets.size):
        msgs.info('DET: {}'.format(calib_dets[i]))
        # Assign RA,DEC, OBJNAME to detected objects and add undetected objects
        if calib_slits[i].maskdef_designtab is not None:
            # Assign slitmask design information to detected objects
            sobjs = calib_slits[i].assign_maskinfo(sobjs, platescale[i], spat_flexure[i],
                                                   TOLER=slitmask_par['obj_toler'])

            if slitmask_par['extract_missing_objs']:
                # Set the FWHM for the extraction of missing objects
                fwhm = calib_slits[i].get_maskdef_extract_fwhm(sobjs, platescale[i],
                                                               slitmask_par['missing_objs_fwhm'], find_fwhm)
                # Assign undetected objects
                sobjs = calib_slits[i].mask_add_missing_obj(sobjs, spat_flexure[i], fwhm,
                                                            slitmask_par['missing_objs_boxcar_rad']/platescale[i])

    return sobjs



