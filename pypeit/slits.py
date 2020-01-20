"""
Implements the objects used to hold slit edge data.

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import inspect

from IPython import embed

import numpy as np

from pypeit import msgs
from pypeit import masterframe
from pypeit import datamodel

class SlitTraceSet(datamodel.DataContainer):
    """
    Defines a generic class for holding and manipulating image traces
    organized into left-right slit pairs.

    Instantiation arguments map directly to the object
    :attr:`datamodel`.
    """
    # Set the version of this class
    version = '1.0.0'
    # Define the data model
    datamodel = {'spectrograph': dict(otype=str, descr='Spectrograph used to take the data.'),
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
                 'nslits': dict(otype=int, descr='Number of slits.'),
                 'id': dict(otype=np.ndarray, atype=int, descr='Slit ID number'),
                 'left': dict(otype=np.ndarray, atype=float,
                              descr='Spatial coordinates (pixel indices) of all left edges, one '
                                    'per slit.  Shape is Nspec by Nslits.'),
                 'right': dict(otype=np.ndarray, atype=float,
                              descr='Spatial coordinates (pixel indices) of all right edges, one '
                                    'per slit.  Shape is Nspec by Nslits.'),
                 'left_tweak': dict(otype=np.ndarray, atype=float,
                                    descr='Spatial coordinates (pixel indices) of all left '
                                          'edges, one per slit.  These traces have been adjusted '
                                          'by the flat-field.  Shape is Nspec by Nslits.'),
                 'right_tweak': dict(otype=np.ndarray, atype=float,
                                     descr='Spatial coordinates (pixel indices) of all right '
                                           'edges, one per slit.  These traces have been adjusted '
                                           'by the flat-field.  Shape is Nspec by Nslits.'),
                 'center': dict(otype=np.ndarray, atype=float,
                               descr='Spatial coordinates of the slit centers.  Shape is Nspec '
                                     'by Nslits.'),
                 'mask': dict(otype=np.ndarray, atype=bool,
                              descr='Bad-slit mask (good slits are False).  Shape is Nslits.'),
                 'specmin': dict(otype=np.ndarray, atype=float,
                                 descr='Minimum spectral position allowed for each slit/order.  '
                                       'Shape is Nslits.'),
                 'specmax': dict(otype=np.ndarray, atype=float,
                                 descr='Maximum spectral position allowed for each slit/order.  '
                                       'Shape is Nslits.')}
    """Provides the class data model."""
    # NOTE: The docstring above is for the ``datamodel`` attribute.

    # TODO: Allow tweaked edges to be arguments?
    # TODO: May want nspat to be a required argument.
    def __init__(self, left, right, nspat=None, spectrograph=None, mask=None, specmin=None,
                 specmax=None, binspec=1, binspat=1, pad=0):
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        super(SlitTraceSet, self).__init__({k: values[k] for k in args[1:]})

    def _validate(self):
        """
        Validate the slit traces.
        """
        if self.left.shape != self.right.shape:
            raise ValueError('Input left and right traces should have the same shape.')
        if self.left.ndim == 1:
            # Object contains a single slit.  Expand the dimensions.
            self.left = np.expand_dims(self.left, 1)
            self.right = np.expand_dims(self.right, 1)

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

        self.nspec, self.nslits = self.left.shape

        # Center is always defined by the original traces, not the
        # tweaked ones. TODO: Is that the behavior we want? An argument
        # in favor is that this means that the slit IDs are always tied
        # to the original traces, not the tweaked ones.
        self.center = (self.left+self.right)/2

        if self.nspat is None:
            # TODO: May want nspat to be a required argument given the
            # only other option is this kludge, which should basically
            # never be useful.
            self.nspat = np.amax(np.append(self.left, self.right))
        if self.id is None:
            self._set_slitids()
        if self.spectrograph is None:
            self.spectrograph = 'unknown'
        if self.mask is None:
            self.mask = np.zeros(self.nslits, dtype=bool)
        if self.specmin is None:
            self.specmin = np.full(self.nslits, -1, dtype=float)
        if self.specmax is None:
            self.specmax = np.full(self.nslits, self.nspec, dtype=float)

        # Make sure mask, specmin, and specmax are at least 1D arrays.
        # TODO: Is there a way around this?
        self.mask = np.atleast_1d(self.mask)
        self.specmin = np.atleast_1d(self.specmin)
        self.specmax = np.atleast_1d(self.specmax)

    def _bundle(self):
        """
        Bundle the data in preparation for writing to a fits file.

        See :func:`pypeit.datamodel.DataContainer._bundle`. Data is
        always written to a 'SLITS' extension.
        """
        return super(SlitTraceSet, self)._bundle(ext='SLITS', transpose_arrays=True)

    @classmethod
    def _parse(cls, hdu):
        """
        Parse the data that was previously written to a fits file.

        See :func:`pypeit.datamodel.DataContainer._parse`. Data is
        always read from the 'SLITS' extension.
        """
        return super(SlitTraceSet, cls)._parse(hdu, ext='SLITS', transpose_table_arrays=True)

    def _set_slitids(self, specfrac=0.5):
        """
        Assign the slit ID numbers.

        Slit IDs are set based on the fractional location of the slit
        center at the provided fractional location along the spectral
        direction, and it is computed as follows::

            self.id = np.round(self.center[int(np.round(specfrac*self.nspec)),:]
                                / self.nspat * 1e4).astype(int)

        I.e., if the slit center is 1024 for an image with 2048
        spatial pixels, the slit ID is 5000. The slit ID numbers are
        assigned to :attr:`id`.

        Args:
            specfrac (:obj:`float`, optional):
                The fractional spectral position of the slit center
                to use for the ID calculation.

        Raises:
            ValueError:
                Raised if :attr:`center`, :attr:`nspat`, or
                :attr:`nspec` is not defined (None).
        """
        # TODO: This doesn't need to be a separate private function. It
        # could just be a line in _validate().
        if self.center is None or self.nspec is None or self.nspat is None:
            raise ValueError('Object does not have center, nspec, or nspat; '
                             'cannot assign slit IDs.')
        self.id = np.round(self.center[int(np.round(specfrac*self.nspec)),:]
                                / self.nspat * 1e4).astype(int)

    def init_tweaked(self):
        """
        Initialize the tweaked slits.
        """
        self.left_tweak = self.left.copy()
        self.right_tweak = self.right.copy()

    def rm_tweaked(self):
        """
        Delete (set to None) the tweaked traces.
        """
        self.left_tweak = None
        self.right_tweak = None

    def select_edges(self, original=False):
        """
        Select between the original or tweaked slit edges.

        By default, the method will return the tweaked slits if they
        have been defined. If they haven't been defined the nominal
        edges (:attr:`left` and :attr:`right`) are returned. Use
        ``original=True`` to return the nominal edges regardless of
        the presence of the tweaked edges.

        Args:
            original (:obj:`bool`, optional):
                To use the nominal edges regardles of the presence of
                the tweaked edges, set this to True.

        Returns:
            tuple: Returns the arrays containing the left and right
            edge coordinates, respectively. These are returned as
            pointers to the internal attributes, **not** copies.
        """
        # TODO: Add a copy argument?
        if self.left_tweak is not None and self.right_tweak is not None and not original:
            return self.left_tweak, self.right_tweak
        return self.left, self.right

    def slit_img(self, pad=None, slitids=None, original=False):
        r"""
        Construct an image identifying each pixel with its associated
        slit.

        The output image has the same shape as the original trace
        image. Each pixel Each pixel in the image is set to the index
        of its associated slit (i.e, the pixel value is
        :math:`0..N_{\rm slit}-1`). Pixels not associated with any
        slit are given values of -1.

        The width of the slit is extended at either edge by a fixed
        number of pixels using the `pad` parameter in :attr:`par`.
        This value can be overridden using the method keyword
        argument.

        .. warning::

            The function does not check that pixels end up in
            multiple pixels or that the padding is sensible given the
            separation between slit edges!

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
            slitids (:obj:`int`, array_like, optional):
                List of slit IDs to include in the image. If None,
                all slits are included.
            original (:obj:`bool`, optional):
                By default, the method will use the tweaked slit
                edges if they have been defined. If they haven't
                been, the nominal edges (:attr:`left` and
                :attr:`right`) are used. To use the nominal edges
                regardless of the presence of the tweaked edges, set
                this to True. See :func:`select_edges`.

        Returns:
            `numpy.ndarray`_: The image with the slit index
            identified for each pixel.
        """
        # Check the input
        if pad is None:
            pad = self.pad
        _pad = pad if isinstance(pad, tuple) else (pad,pad)
        if len(_pad) != 2:
            msgs.error('Padding for both left and right edges should be provided as a 2-tuple!')

        # Slit IDs to include in the image
        slitids = np.arange(self.nslits) if slitids is None else np.atleast_1d(slitids).ravel()

        # Pixel coordinates
        spat = np.arange(self.nspat)
        spec = np.arange(self.nspec)

        # Choose the slit edges to use
        left, right = self.select_edges(original=original)

        # TODO: When specific slits are chosen, need to check that the
        # padding doesn't lead to slit overlap.

        # Find the pixels in each slit, limited by the minimum and
        # maximum spectral position.
        slitid_img = np.full((self.nspec,self.nspat), -1, dtype=int)
        for i in slitids:
            indx = (spat[None,:] > left[:,i,None] - _pad[0]) \
                        & (spat[None,:] < right[:,i,None] + _pad[1]) \
                        & (spec > self.specmin[i])[:,None] & (spec < self.specmax[i])[:,None]
            slitid_img[indx] = i
        return slitid_img

    def spatial_coordinate_image(self, slitids=None, full=False, slitid_img=None, pad=None,
                                 original=False):
        r"""
        Generate an image with the normalized spatial coordinate
        within each slit.

        Args:
            slitids (:obj:`int`, array_like, optional):
                List of slit IDs to include in the image. If None,
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
            original (:obj:`bool`, optional):
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
        # Slit IDs to include
        _slitids = np.arange(self.nslits) if slitids is None else np.atleast_1d(slitids).ravel()
        if full and len(_slitids) > 1:
            msgs.error('For a full image with the slit coordinates, must select a single slit.')

        # Generate the slit ID if it wasn't provided
        if not full:
            if slitid_img is None:
                slitid_img = self.slit_img(pad=pad, slitids=_slitids, original=original)
            if slitid_img.shape != (self.nspec,self.nspat):
                msgs.error('Provided slit ID image does not have the correct shape!')

        # Choose the slit edges to use
        left, right = self.select_edges(original=original)

        # Slit width
        slitwidth = right - left

        # TODO: This check should go in it's own function and/or
        # checked when instantiating the object.
        # Slit widths have to be larger than 0
        indx = slitwidth <= 0.
        if np.any(indx[:,_slitids]):
            bad_slits = np.where(np.any(indx, axis=0))[0]
            # TODO: Shouldn't this fault?
            msgs.warn('Slits {0} have negative (or 0) slit width!'.format(bad_slits))

        # Output image
        coo_img = np.zeros((self.nspec,self.nspat), dtype=float)
        spat = np.arange(self.nspat)
        for i in _slitids:
            coo = (spat[None,:] - left[:,i,None])/slitwidth[:,i,None]
            if not full:
                indx = slitid_img == i
                coo_img[indx] = coo[indx]
            else:
                coo_img = coo
        return coo_img

    def spatial_coordinates(self, original=False):
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
        left, right = self.select_edges(original=original)
        return slit_spat_pos(left, right, self.nspat)

    # TODO: Temporary until it's fully removed
    def to_tslits_dict(self):
        tslits_dict = {}
        tslits_dict['slit_left'] = self.left.copy()
        tslits_dict['slit_righ'] = self.right.copy()
        tslits_dict['slitcen'] = self.center.copy()
        tslits_dict['maskslits'] = self.mask.copy()
        tslits_dict['nspec'] = self.nspec
        tslits_dict['nspat'] = self.nspat
        tslits_dict['nslits'] = self.nslits
        tslits_dict['binspectral'] = self.binspec
        tslits_dict['binspatial'] = self.binspat
        tslits_dict['spectrograph'] = self.spectrograph
        tslits_dict['spec_min'] = self.specmin.copy()
        tslits_dict['spec_max'] = self.specmax.copy()
        tslits_dict['pad'] = self.pad
        return tslits_dict


