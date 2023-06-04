"""
Implements a general purpose object used to decompose and predict
traces using principle-component analysis.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import warnings
from IPython import embed

import numpy as np

from astropy.io import fits

from pypeit import msgs
from pypeit import utils
from pypeit.io import hdu_iter_by_ext
from pypeit.core import trace
from pypeit.core import pca
from pypeit.core.fitting import PypeItFit
from pypeit.datamodel import DataContainer

# TODO: This is even more general than a "trace" PCA. Rename to
# "VectorPCA"?

class TracePCA(DataContainer):
    r"""
    Class to build and interact with PCA model of traces.

    This is primarily a container class for the results of
    :func:`pypeit.core.pca.pca_decomposition`,
    :func:`pypeit.core.pca.fit_pca_coefficients`, and
    :func:`pypeit.core.pca.pca_predict`.

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_tracepca.rst

    Args:
        trace_cen (`numpy.ndarray`_, optional):
            A floating-point array with the spatial location of each
            each trace. Shape is :math:`(N_{\rm spec}, N_{\rm
            trace})`. If None, the object is "empty" and all of the
            other keyword arguments are ignored.
        npca (:obj:`bool`, optional):
            The number of PCA components to keep. See
            :func:`pypeit.core.pca.pca_decomposition`.
        pca_explained_var (:obj:`float`, optional):
            The percentage (i.e., not the fraction) of the variance
            in the data accounted for by the PCA used to truncate the
            number of PCA coefficients to keep (see `npca`). Ignored
            if `npca` is provided directly. See
            :func:`pypeit.core.pca.pca_decomposition`.
        reference_row (:obj:`int`, optional):
            The row (spectral position) in `trace_cen` to use as the
            reference coordinate system for the PCA. If None, set to
            the :math:`N_{\rm spec}/2`.
        coo (`numpy.ndarray`_, optional):
            Floating-point array with the reference coordinates for
            each trace. If provided, the shape must be :math:`(N_{\rm
            trace},)`. If None, the reference coordinate system is
            defined by the value of `trace_cen` at the spectral
            position defined by `reference_row`. See the `mean`
            argument of :func:`pypeit.core.pca.pca_decomposition`.
    """

    version = '1.1.0'
    """Datamodel version."""

    datamodel = {'reference_row': dict(otype=(int, np.integer),
                                       descr='The row (spectral position) used as the reference ' \
                                             'coordinate system for the PCA.'),
                 'trace_coo': dict(otype=np.ndarray, atype=(float,np.floating),
                                   descr='Trace coordinates.  Shape must be '
                                         r':math:`(N_{\rm spec},N_{\rm trace})`.'),
                 'nspec': dict(otype=int,
                               descr='Number of pixels in the image spectral direction.'),
                 'ntrace': dict(otype=int, descr='Number of traces used to construct the PCA.'),
                 'input_npca': dict(otype=int,
                                    descr='Requested number of PCA components if provided.'),
                 'npca': dict(otype=int, descr='Number of PCA components used.'),
                 'input_pcav': dict(otype=(float,np.floating),
                                    descr='Requested variance accounted for by PCA decomposition.'),
                 'pca_coeffs': dict(otype=np.ndarray, atype=(float,np.floating),
                                    descr=r'PCA component coefficients. If the PCA decomposition '
                                          r'used :math:`N_{\rm comp}` components for '
                                          r':math:`N_{\rm vec}` vectors, the shape of this array '
                                          r'must be :math:`(N_{\rm vec}, N_{\rm comp})`. The '
                                          r'array can be 1D with shape :math:`(N_{\rm vec},)` if '
                                          r'there was only one PCA component.'),
                 'pca_components': dict(otype=np.ndarray, atype=(float,np.floating),
                                        descr=r'Vectors with the PCA components.  Shape must be '
                                              r':math:`(N_{\rm comp}, N_{\rm spec})`.'),
                 'pca_mean': dict(otype=np.ndarray, atype=(float,np.floating),
                                  descr=r'The mean offset of the PCA decomposotion for each '
                                        r' spectral pixel. Shape is :math:`(N_{\rm spec},)`.'),
                 'pca_coeffs_model': dict(otype=np.ndarray, atype=PypeItFit,
                                          descr='An array of PypeItFit objects, one per PCA '
                                                'component, that models the trend of the PCA '
                                                'component coefficients with the reference '
                                                'coordinate of each vector.  These models are '
                                                r'used by :func:`predict` to model the expected '
                                                'coefficients at a new reference coordinate.')}
    """Object datamodel."""

    internals = ['is_empty']

    # TODO: Add a show method that plots the pca coefficients and the
    # current fit, if there is one
    def __init__(self, trace_cen=None, npca=None, pca_explained_var=99.0, reference_row=None,
                 coo=None):

        # Instantiate as an empty DataContainer
        super(TracePCA, self).__init__()
        self.is_empty = True

        # Only do the decomposition if the trace coordinates are provided.
        if trace_cen is not None:
            self.decompose(trace_cen, npca=npca, pca_explained_var=pca_explained_var,
                           reference_row=reference_row, coo=coo)

    def decompose(self, trace_cen, npca=None, pca_explained_var=99.0, reference_row=None,
                  coo=None):
        r"""
        Construct the PCA from scratch.
    
        Args:
            trace_cen (`numpy.ndarray`_):
                A floating-point array with the spatial location of
                each each trace. Shape is :math:`(N_{\rm spec},
                N_{\rm trace})`.  Cannot be None.
            npca (:obj:`bool`, optional):
                The number of PCA components to keep. See
                :func:`pypeit.core.pca.pca_decomposition`.
            pca_explained_var (:obj:`float`, optional):
                The percentage (i.e., not the fraction) of the
                variance in the data accounted for by the PCA used to
                truncate the number of PCA coefficients to keep (see
                `npca`). Ignored if `npca` is provided directly. See
                :func:`pypeit.core.pca.pca_decomposition`.
            reference_row (:obj:`int`, optional):
                The row (spectral position) in `trace_cen` to use as
                the reference coordinate system for the PCA. If None,
                set to the :math:`N_{\rm spec}/2`.
            coo (`numpy.ndarray`_, optional):
                Floating-point array with the reference coordinates
                for each trace. If provided, the shape must be
                :math:`(N_{\rm trace},)`. If None, the reference
                coordinate system is defined by the value of
                `trace_cen` at the spectral position defined by
                `reference_row`. See the `mean` argument of
                :func:`pypeit.core.pca.pca_decomposition`.
        """
        if trace_cen is None:
            raise ValueError('Must provide traces to construct the PCA.')
        self._reinit()
        self.is_empty = False

        # Set the reference row to use for the coordinates of the trace
        self.reference_row = trace_cen.shape[0]//2 if reference_row is None else reference_row
        self.ntrace = trace_cen.shape[1]

        self.trace_coo = trace_cen[self.reference_row,:] if coo is None else np.atleast_1d(coo)
        if self.trace_coo.size != self.ntrace:
            raise ValueError('Provided reference coordinates have incorrect shape.')

        # Save the input
        self.input_npca = npca
        self.input_pcav = pca_explained_var

        # Perform the PCA decomposition of the traces
        self.pca_coeffs, self.pca_components, self.pca_mean, _ \
                = pca.pca_decomposition(trace_cen.T, npca=npca,
                                        pca_explained_var=pca_explained_var, mean=self.trace_coo)
        self.npca = self.pca_coeffs.shape[1]
        self.nspec = self.pca_components.shape[1]

    def _reinit(self):
        """
        Erase and/or define all the attributes of the class.
        """
        # The following are *all* objects assigned to self.
        self.is_empty = True
        self.reference_row = None
        self.ntrace = None
        self.trace_coo = None
        self.input_npca = None
        self.input_pcav = None
        self.pca_coeffs = None
        self.pca_components = None
        self.pca_mean = None
        self.npca = None
        self.nspec = None
        self.pca_coeffs_model = None

    def build_interpolator(self, order, ivar=None, weights=None, function='polynomial', lower=3.0,
                           upper=3.0, maxrej=1, maxiter=25, minx=None, maxx=None, debug=False):
        """
        Wrapper for :func:`fit_pca_coefficients` that uses class
        attributes and saves the input parameters.
        """
        if self.is_empty:
            raise ValueError('TracePCA object is empty; re-instantiate or run decompose().')
        self.pca_coeffs_model \
                = pca.fit_pca_coefficients(self.pca_coeffs, order, ivar=ivar, weights=weights,
                                           function=function, lower=lower, upper=upper,
                                           minx=minx, maxx=maxx, maxrej=maxrej, maxiter=maxiter,
                                           coo=self.trace_coo, debug=debug)

    def predict(self, x):
        r"""
        Predict one or more traces given the functional forms for the
        PCA coefficients.

        .. warning::
            The PCA coefficients must have first been modeled by a
            function before using this method. An error will be
            raised if :attr:`fit_coeff` is not defined.

        Args:
            x (:obj:`float`, `numpy.ndarray`_):
                One or more spatial coordinates (at the PCA reference
                row) at which to sample the PCA coefficients and
                produce the PCA model for the trace spatial position
                as a function of spectral pixel.

        Returns:
            `numpy.ndarray`_: The array with the predicted spatial
            locations of the trace. If the provided coordinate is a
            single value, the returned shape is :math:`(N_{\rm
            pix},)`; otherwise it is :math:`(N_{\rm pix}, N_{\rm
            x})`.
        """
        if self.is_empty:
            raise ValueError('TracePCA object is empty; re-instantiate or run decompose().')
        if self.pca_coeffs_model is None:
            msgs.error('PCA coefficients have not been modeled; run build_interpolator first.')
        return pca.pca_predict(x, self.pca_coeffs_model, self.pca_components, self.pca_mean, x).T

    def _bundle(self, ext='PCA'):
        """Bundle the data for writing."""
        d = super(TracePCA, self)._bundle(ext=ext)
        if self.pca_coeffs_model is None:
            return d

        # Fix the default bundling to handle the fact that
        # pca_coeffs_model is an array of PypeItFit instances.
        _d = d[0] if ext is None else d[0][ext]
        del _d['pca_coeffs_model']
        d += [{'{0}_MODEL_{1}'.format(ext,i+1): self.pca_coeffs_model[i]} 
                for i in range(self.npca)]
        return d

    # TODO: Although I don't like doing it, kwargs is here to catch the
    # extraneous keywords that can be passed to _parse from the base class but
    # won't be used.
    @classmethod
    def _parse(cls, hdu, hdu_prefix=None, **kwargs):
        """
        Parse the data from the provided HDU.

        See :func:`pypeit.datamodel.DataContainer._parse` for the
        argument descriptions.
        """
        # Run the default parser to get most of the data
        d, version_passed, type_passed, parsed_hdus \
                = super(TracePCA, cls)._parse(hdu, hdu_prefix=hdu_prefix)

        # This should only ever read one hdu!
        if len(parsed_hdus) > 1:
            msgs.error('CODING ERROR: Parsing saved TracePCA instances should only parse 1 HDU, '
                       'independently of the PCA PypeItFit models.')

        # Check if any models exist
        if hasattr(hdu, '__len__') \
                and np.any(['{0}_MODEL'.format(parsed_hdus[0]) in h.name for h in hdu]):
            # Parse the models
            model_ext = ['{0}_MODEL_{1}'.format(parsed_hdus[0],i+1) for i in range(d['npca'])]
            d['pca_coeffs_model'] = np.array([PypeItFit.from_hdu(hdu[e]) for e in model_ext])
            parsed_hdus += model_ext

        return d, version_passed, type_passed, parsed_hdus

    @classmethod
    def from_dict(cls, d=None):
        """
        Instantiate from a dictionary.

        This is a basic wrapper for
        :class:`pypeit.datamodel.DataContainer.from_dict` that
        appropriately toggles :attr:`is_empty`.
        """
        self = super(TracePCA, cls).from_dict(d=d)
        self.is_empty = False
        return self


# TODO: Like with the use of TracePCA in EdgeTraceSet, we should
# integrate the elements of the function below into classes that trace
# objects and tilts so that the PCA can be called and used later
def pca_trace_object(trace_cen, order=None, trace_bpm=None, min_length=0.6, npca=None,
                     pca_explained_var=99.0, reference_row=None, coo=None, minx=None, maxx=None,
                     trace_wgt=None, function='polynomial', lower=3.0, upper=3.0, maxrej=1,
                     maxiter=25, debug=False):
    r"""
    Decompose and reconstruct the provided traces using
    principle-component analysis.

    Args:
        trace_cen (`numpy.ndarray`_):
            A floating-point array with the spatial location of each
            each trace. Shape is :math:`(N_{\rm spec}, N_{\rm
            trace})`.
        order (:obj:`int`, :obj:`list`, optional):
            The order of the polynomial to use fit each PCA
            coefficient as a function of trace position. If None,
            `order` is set to :math:`3.3 N_{\rm use}/N_{\rm trace}`,
            where :math:`N_{\rm use}` is the number of traces used to
            construct the PCA and :math:`N_{\rm trace}` is the number
            of total traces provided. If an integer (determined
            automatically if the argument is `None`), the order per
            PCA component (see `npca`) is set to cascade from
            high-to-low order as follows::

                _order = np.clip(order - np.arange(npca), 1, None).astype(int)

        trace_bpm (`numpy.ndarray`_, optional):
            Bad-pixel mask for the trace data (True is bad; False is
            good). Must match the shape of `trace_cen`.
        min_length (:obj:`float`, optional):
            The good position of the trace must cover at least this
            fraction of the spectral dimension for use in the PCA
            decomposition.
        npca (:obj:`bool`, optional):
            The number of PCA components to keep. See
            :func:`pypeit.core.pca.pca_decomposition`.
        pca_explained_var (:obj:`float`, optional):
            The percentage (i.e., not the fraction) of the variance
            in the data accounted for by the PCA used to truncate the
            number of PCA coefficients to keep (see `npca`). Ignored
            if `npca` is provided directly. See
            :func:`pypeit.core.pca.pca_decomposition`.
        reference_row (:obj:`int`, optional):
            The row (spectral position) in `trace_cen` to use as the
            reference coordinate system for the PCA. If None, set to
            the :math:`N_{\rm spec}/2` or based on the spectral
            position that crosses the most number of valid trace
            positions.
        coo (`numpy.ndarray`_, optional):
            Floating-point array with the reference coordinates to
            use for each trace. If None, coordinates are defined at
            the reference row of `trace_cen`. Shape must be
            :math:`(N_{\rm trace},)`.
        minx, maxx (:obj:`float`, optional):
            Minimum and maximum values used to rescale the
            independent axis data. If None, the minimum and maximum
            values of `coo` are used. See
            :func:`utils.robust_polyfit_djs`.
        trace_wgt (`numpy.ndarray`_, optional):
            Weights to apply to the PCA coefficient of each trace
            during the fit. Weights are independent of the PCA
            component. See `weights` parameter of
            :func:`pypeit.core.pca.fit_pca_coefficients`. Shape must
            be :math:`(N_{\rm trace},)`.
        function (:obj:`str`, optional):
            Type of function used to fit the data.
        lower (:obj:`float`, optional):
            Number of standard deviations used for rejecting data
            **below** the mean residual. If None, no rejection is
            performed. See :func:`utils.robust_polyfit_djs`.
        upper (:obj:`float`, optional):
            Number of standard deviations used for rejecting data
            **above** the mean residual. If None, no rejection is
            performed. See :func:`utils.robust_polyfit_djs`.
        maxrej (:obj:`int`, optional):
            Maximum number of points to reject during fit iterations.
            See :func:`utils.robust_polyfit_djs`.
        maxiter (:obj:`int`, optional):
            Maximum number of rejection iterations allows. To force
            no rejection iterations, set to 0.
        debug (:obj:`bool`, optional):
            Show plots useful for debugging.
    """
    # Check the input
    if trace_bpm is None:
        use_trace = np.ones(trace_cen.shape[1], dtype=bool)
        _reference_row = trace_cen.shape[0]//2 if reference_row is None else reference_row
    else:
        use_trace = np.sum(np.invert(trace_bpm), axis=0)/trace_cen.shape[0] > min_length
        _reference_row = trace.most_common_trace_row(trace_bpm[:,use_trace]) \
                                if reference_row is None else reference_row
    _coo = None if coo is None else coo[use_trace]

    # Instantiate the PCA
    cenpca = TracePCA(trace_cen[:,use_trace], npca=npca, pca_explained_var=pca_explained_var,
                      reference_row=_reference_row, coo=_coo)

    # Set the order of the function fit to the PCA coefficients:
    # Order is set to cascade down to lower order for components
    # that account for a smaller percentage of the variance.
    if order is None:
        # TODO: Where does this come from?
        order = int(np.clip(np.floor(3.3*np.sum(use_trace)/trace_cen.shape[1]),1.0,3.0))
    _order = np.atleast_1d(order)
    if _order.size == 1:
        _order = np.clip(order - np.arange(cenpca.npca), 1, None).astype(int)
    if _order.size != cenpca.npca:
        msgs.error('Number of polynomial orders does not match the number of PCA components.')
    msgs.info('Order of function fit to each component: {0}'.format(_order))

    # Apply a 10% relative error to each coefficient. This performs
    # better than use_mad, since larger coefficients will always be
    # considered inliers, if the coefficients vary rapidly with
    # order as they sometimes do.

    # TODO: This inverse variance usage has performance issues and
    # tends to lead to rejection of coefficients that are near 0.
    # Instead of setting the floor to an absolute value 0.1, why not a
    # relative value like the mean or median of the coefficients? I.e.
#    ivar = utils.inverse(numpy.square(np.fmax(0.1*np.absolute(cenpca.pca_coeffs),
#                                              0.1*np.median(cenpca.pca_coeffs))))
    ivar = utils.inverse(np.square(np.fmax(0.1*np.absolute(cenpca.pca_coeffs), 0.1)))

    # Set any additional weights for each trace
    weights = np.ones(np.sum(use_trace), dtype=float) \
                if trace_wgt is None else trace_wgt[use_trace]

    # TODO: This combination of ivar and weights is as it has been
    # previously. However, we recently changed the slit-edge tracing
    # code to fit the PCA coefficients with unity weights (the default
    # when passing weights=None to build_interpolator) and ivar=None.
    # The means the PCA coefficients are fit with uniform weighting and
    # the rejection is based on the median absolute deviation of the
    # data with respect to the fitted model.
    # 
    # This current call below to build_interpolator will instead weight
    # the fit according to the weights set above, and it will reject
    # points based on the inverse variance set above. We need to check
    # that this makes sense!

    # Build the interpolator that allows prediction of new traces
    cenpca.build_interpolator(_order, ivar=ivar, weights=weights, function=function,
                              lower=lower, upper=upper, maxrej=maxrej, maxiter=maxiter,
                              minx=minx, maxx=maxx, debug=debug)

    # Return the traces predicted for all input traces
    return cenpca.predict(trace_cen[_reference_row,:] if coo is None else coo)

