"""
Implements a general purpose object used to decompose and predict
traces using principle-component analysis.

.. _numpy.ndarray: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
"""

import numpy as np

from pypeit import utils
from pypeit.core import trace
from pypeit.core import pca

# TODO: This is even more general than a "trace" PCA. Could think of a
# more general name.

class TracePCA:
    r"""
    Class to build and interact with PCA model of traces.

    This is primarily a container class for the results of
    :func:`pypeit.core.pca.pca_decomposition`,
    :func:`pypeit.core.pca.fit_pca_coefficients`, and
    :func:`pypeit.core.pca.pca_predict`.

    Args:
        trace_cen (`numpy.ndarray`_):
            A floating-point array with the spatial location of each
            each trace. Shape is :math:`(N_{\rm spec}, N_{\rm
            trace})`.
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
    """
    # TODO: Add a show method that plots the pca coefficients and the
    # current fit, if there is one
    def __init__(self, trace_cen, npca=None, pca_explained_var=99.0, reference_row=None):

        # Set the reference row to use for the coordinates of the trace
        self.reference_row = trace_cen.shape[0]//2 if reference_row is None else reference_row
        self.trace_coo = trace_cen[self.reference_row,:]
        # Save the input
        self.input_npca = npca
        self.input_pcav = pca_explained_var

        # Perform the PCA decomposition of the traces
        self.pca_coeffs, self.pca_components, self.pca_mean, self.trace_mean \
                = pca.pca_decomposition(trace_cen.T, npca=npca,
                                        pca_explained_var=pca_explained_var, mean=self.trace_coo)
        self.npca = self.pca_coeffs.shape[1]
        self.nspec = self.pca_components.shape[1]

        # Instantiate the remaining attributes
        self.pca_mask = np.zeros(self.pca_coeffs.shape, dtype=bool)
        self.function = None
        self.fit_coeff = None
        self.minx = None
        self.maxx = None
        self.lower = None
        self.upper = None
        self.maxrej = None
        self.maxiter = None

    def build_interpolator(self, order, ivar=None, weights=None, function='polynomial', lower=3.0,
                           upper=3.0, minx=None, maxx=None, maxrej=1, maxiter=25, debug=False):
        """
        Wrapper for :func:`fit_pca_coefficients` that uses class
        attributes and saves the input parameters.
        """
        # Save the input
        # TODO: Keep the ivar and weights?
        self.function = function
        self.lower = lower
        self.upper = upper
        self.maxrej = maxrej
        self.maxiter = maxiter
        self.pca_mask, self.fit_coeff, self.minx, self.maxx \
                = pca.fit_pca_coefficients(self.pca_coeffs, order, ivar=ivar, weights=weights,
                                           function=self.function, lower=lower, upper=upper,
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
        if self.fit_coeff is None:
            msgs.error('PCA coefficients have not been modeled; run model_coeffs first.')
        return pca.pca_predict(x, self.fit_coeff, self.pca_components, self.pca_mean, x,
                               function=self.function).T

# TODO: This is a place holder; it's never called and will fault if it is!
def pca_trace_object(trace_cen, order=None, trace_bpm=None, min_length=0.6, npca=None,
                     pca_explained_var=99.0, coeff_weights=None, function='polynomial', lower=3.0,
                     upper=3.0, minx=None, maxx=None, maxrej=1, maxiter=25, debug=False):

        if trace_bpm is None:
            use_trace = np.ones(trace_cen.shape[1], dtype=bool)
            reference_row = trace_cen.shape[0]//2
        else:
            use_trace = np.sum(np.invert(trace_bpm), axis=0)/trace_cen.shape[0] > min_length
            reference_row = trace.most_common_trace_row(trace_bpm[:,use_trace])

        # Instantiate the PCA
        cenpca = TracePCA(trace_cen[:,use_trace], npca=npca, pca_explained_var=pca_explained_var,
                          reference_row=reference_row)

        # Set the order of the function fit to the PCA coefficiencts:
        # Order is set to cascade down to lower order for components
        # that account for a smaller percentage of the variance.
        if order is None:
            order = int(np.clip(np.floor(3.3*np.sum(use_trace)/trace_cen.shape[1]),1.0,3.0))
        _order = np.clip(order - np.arange(self.pca[i].npca), 1, None).astype(int)
        msgs.info('Order of function fit to each component: {0}'.format(_order))

        # Apply a 10% relative error to each coefficient. This performs
        # better than use_mad, since larger coefficients will always be
        # considered inliers, if the coefficients vary rapidly with
        # order as they sometimes do.
        ivar = utils.inverse(numpy.square(np.fmax(0.1*np.absolute(cenpca.pca_coeffs), 0.1)))
        weights = np.fmax(sobjs_final[indx_obj_id].ech_snr, 1.0)

        # Run the fit
        cenpca.build_interpolator(_order, ivar=ivar, weights=weights, function=function,
                                  lower=lower, upper=upper, maxrej=maxrej, maxiter=maxiter,
                                  debug=debug)

        return cenpca.predict(trace_cen[:,use_trace])
