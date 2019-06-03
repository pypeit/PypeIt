"""
Implements edge tracing.

.. _numpy.ndarray: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
"""
import os
import time
import inspect
from collections import Counter, OrderedDict

import numpy as np

from scipy import ndimage, signal, interpolate

import matplotlib
from matplotlib import pyplot as plt
from matplotlib import ticker, rc

from sklearn.decomposition import PCA

from astropy.io import fits
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy import table

from pypeit.bitmask import BitMask
from pypeit.core import parse, pydl, procimg, arc
from pypeit.par.parset import ParSet

from pypeit import utils
from pypeit import ginga
from pypeit import masterframe
from pypeit import msgs
from pypeit.traceimage import TraceImage
from pypeit.spectrographs.util import load_spectrograph
from pypeit.par.pypeitpar import TraceSlitsPar
from pypeit import io

from pypeit import sampling
from pypeit import moment
from pypeit.spectrographs import slitmask

class TraceBitMask(BitMask):
    """
    Mask bits used during slit tracing.
    """
    # TODO: Create a script that will dynamically write the used bits
    # to a doc for readthedocs.
    def __init__(self):
        # TODO: This needs to be an OrderedDict for now to ensure that
        # the bit assigned to each key is always the same. As of python
        # 3.7, normal dict types are guaranteed to preserve insertion
        # order as part of its data model. When/if we require python
        # 3.7, we can remove this (and other) OrderedDict usage in
        # favor of just a normal dict.
        mask = OrderedDict([
                       ('NOEDGE', 'No edge found/input for this trace in this column.'),
                    ('MATHERROR', 'A math error occurred during the calculation (e.g., div by 0)'),
                  ('MOMENTERROR', 'Recentering moment calculation had a large error'),
                   ('LARGESHIFT', 'Recentering resulted in a large shift'),
                ('DISCONTINUOUS', 'Pixel included in a trace but part of a discontinuous segment'),
                    ('DUPLICATE', 'Trace is a duplicate based on trace matching tolerance'),
                   ('SHORTRANGE', 'Trace does not meet the minimum spectral range criterion'),
                       ('HITMIN', 'Trace crosses the minimum allowed column'),
                       ('HITMAX', 'Trace crosses the maximum allowed column'),
                  ('OFFDETECTOR', 'Trace lands off the detector'),
                   ('USERINSERT', 'Trace was inserted as requested by user'),
                   ('SYNCINSERT', 'Trace was inserted during left and right edge sync'),
                   ('MASKINSERT', 'Trace was inserted based on drilled slit-mask locations'),
                    ('SHORTSLIT', 'Slit formed by left and right edge is too short'),
                 ('ABNORMALSLIT', 'Slit formed by left and right edge has abnormal length')
                           ])
        super(TraceBitMask, self).__init__(list(mask.keys()), descr=list(mask.values()))

    def bad_flags(self):
        """
        List the flags that mean the trace is bad.
        """
        return list(set(self.bits.keys()) - set(['SYNCINSERT', 'MASKINSERT']))


class EdgeTracePCA:
    r"""
    Class to build and interact with PCA model of edge traces.

    This is primarily a container class for the results of
    :func:`pca_decomposition`, :func:`fit_pca_coefficients`, and
    :func:`pca_predict`.

    Args:
        trace (`numpy.ndarray`_):
            A floating-point array with the spatial location of each
            each trace edge. Shape is :math:`(N_{\rm spec}, N_{\rm
            trace})`.
        npca (:obj:`bool`, optional):
            The number of PCA components to keep. See
            :func:`pca_decomposition`.
        pca_explained_var (:obj:`float`, optional):
            The percentage (i.e., not the fraction) of the variance
            in the data accounted for by the PCA used to truncate the
            number of PCA coefficients to keep (see `npca`). Ignored
            if `npca` is provided directly. See
            :func:`pca_decomposition`.
        reference_row (:obj:`int`, optional):
            The row (spectral position) in `trace` to use as the
            reference coordinate system for the PCA. If None, set to
            the :math:`N_{\rm spec}/2`.
    """
    # TODO: Add a show method that plots the pca coefficients and the
    # current fit, if there is one
    def __init__(self, trace, npca=None, pca_explained_var=99.0, reference_row=None):

        # Set the reference row to use for the coordinates of the trace
        self.reference_row = trace.shape[0]//2 if reference_row is None else reference_row
        self.trace_coo = trace[self.reference_row,:]
        # Save the input
        self.input_npca = npca
        self.input_pcav = pca_explained_var

        # Perform the PCA decomposition of the traces
        self.pca_coeffs, self.pca_components, self.pca_mean, self.trace_mean \
                = pca_decomposition(trace.T, npca=npca, pca_explained_var=pca_explained_var,
                                    mean=self.trace_coo)
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

    def build_interpolator(self, order, function='polynomial', lower=3.0, upper=3.0, minx=None,
                           maxx=None, maxrej=1, maxiter=25, debug=False):
        """
        Wrapper for :func:`fit_pca_coefficients` that uses class
        attributes and saves the input parameters.
        """
        # Save the input
        self.function = function
        self.lower = lower
        self.upper = upper
        self.maxrej = maxrej
        self.maxiter = maxiter
        self.pca_mask, self.fit_coeff, self.minx, self.maxx \
                = fit_pca_coefficients(self.pca_coeffs, order, function=self.function, lower=lower,
                                       upper=upper, minx=minx, maxx=maxx, maxrej=maxrej,
                                       maxiter=maxiter, coo=self.trace_coo, debug=debug)

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
            raise ValueError('PCA coefficients have not been modeled; run model_coeffs first.')
        return pca_predict(x, self.fit_coeff, self.pca_components, self.pca_mean, x,
                           function=self.function).T

#    def trace_map(self, nspat, trace=None):
#        """
#        Construct a map of the spatial position of a full map of
#        trace coordinates.
#
#        Args:
#            nspat (:obj:`int`):
#                Number of spatial pixels (second axis) in the image.
#            trace (`numpy.ndarray`_, optional):
#                A subset of traces to replace the PCA predictions at
#                specific spatial positions. If None, this function is
#                equivalent to::
#
#                    self.predict(np.arange(nspat))
#
#        Returns:
#            `numpy.ndarray`_: Spatial position of each trace starting
#            from the reference row.
#        """
#        # Use the PCA to predict traces at all spatial pixels
#        trace_full = self.predict(np.arange(nspat))
#
#        if trace is None:
#            # No anchor traces to insert
#            return trace_full
#
#        # Force a 2D trace
#        _trace = trace.reshape(-1,1) if trace.ndim == 1 else trace
#        # Anchor the PCA predictions to the provided traces at the
#        # relevant spatial positions
#        trace_full[:,np.round(_trace[self.reference_row,:]).astype(int)] = _trace
#        return trace_full


class EdgeTracePar(ParSet):
    """
    Parameters used for slit edge tracing.
    
    For a table with the current keywords, defaults, and descriptions,
    see :ref:`pypeitpar`.
    """
    prefix = 'ETP'  # Prefix for writing parameters to a header is a class attribute
    def __init__(self, filt_iter=None, sobel_mode=None, edge_thresh=None, follow_span=None,
                 minimum_spec_length=None, valid_flux_thresh=None, max_spat_shift=None,
                 max_spat_error=None, match_tol=None, fit_function=None, fit_order=None,
                 fit_maxdev=None, fit_maxiter=None, fit_niter=None, pca_n=None,
                 pca_var_percent=None, pca_function=None, pca_order=None, pca_sigrej=None,
                 pca_maxrej=None, pca_maxiter=None, smash_range=None, edge_detect_clip=None,
                 trace_median_frac=None, trace_thresh=None, fwhm_uniform=None, niter_uniform=None,
                 fwhm_gaussian=None, niter_gaussian=None, det_buffer=None, max_nudge=None,
                 sync_trace=None, sync_center=None, sync_to_edge=None, min_slit_gap=None,
                 minimum_slit_length=None, length_range=None, clip=None, pad=None, add_slits=None,
                 rm_slits=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])      # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['filt_iter'] = 0
        dtypes['filt_iter'] = int
        descr['filt_iter'] = 'Number of median-filtering iterations to perform on sqrt(trace) ' \
                             'image before applying to Sobel filter to detect slit/order edges.'

        defaults['sobel_mode'] = 'nearest'
        options['sobel_mode'] = EdgeTracePar.valid_sobel_modes()
        dtypes['sobel_mode'] = str
        descr['sobel_mode'] = 'Mode for Sobel filtering.  Default is \'nearest\'; note we find' \
                              '\'constant\' works best for DEIMOS.'

        defaults['edge_thresh'] = 20.
        dtypes['edge_thresh'] = [int, float]
        descr['edge_thresh'] = 'Threshold for finding edges in the Sobel-filtered significance' \
                               ' image.'

        defaults['follow_span'] = 20
        dtypes['follow_span'] = int
        descr['follow_span'] = 'In the initial connection of spectrally adjacent edge ' \
                               'detections, this sets the number of previous spectral rows ' \
                               'to consider when following slits forward.'

        defaults['minimum_spec_length'] = 1/3.
        dtypes['minimum_spec_length'] = [int, float]
        descr['minimum_spec_length'] = 'The minimum spectral length of a trace allowed as a ' \
                                       'fraction of the full detector spectral span.'
        
        defaults['valid_flux_thresh'] = 500.
        dtypes['valid_flux_thresh'] = [int, float]
        descr['valid_flux_thresh'] = 'The flux in the image used to construct the edge traces ' \
                                     'is valid if its median value is above this threshold.  ' \
                                     'Any edge tracing issues are then assumed not to be an ' \
                                     'issue with the trace image itself.'

        defaults['max_spat_shift'] = 0.15
        dtypes['max_spat_shift'] = [int, float]
        descr['max_spat_shift'] = 'Maximum spatial shift in pixels between the edges in ' \
                                  'adjacent spectral positions.'

        defaults['max_spat_error'] = 0.2
        dtypes['max_spat_error'] = [int, float]
        descr['max_spat_error'] = 'Maximum error in the spatial position of edges in pixels.'

        defaults['match_tol'] = 3.
        dtypes['match_tol'] = [int, float]
        descr['match_tol'] = 'Same-side slit edges below this separation in pixels are ' \
                             'considered part of the same edge.'

        defaults['fit_function'] = 'legendre'
        options['fit_function'] = EdgeTracePar.valid_functions()
        dtypes['fit_function'] = str
        descr['fit_function'] = 'Function fit to edge measurements.  ' \
                                'Options are: {0}'.format(', '.join(options['fit_function']))

        defaults['fit_order'] = 5
        dtypes['fit_order'] = int
        descr['fit_order'] = 'Order of the function fit to edge measurements.'

        defaults['fit_maxdev'] = 5.0
        dtypes['fit_maxdev'] = [int, float]
        descr['fit_maxdev'] = 'Maximum deviation between the fitted and measured edge position ' \
                              'for rejection in spatial pixels.'

        defaults['fit_maxiter'] = 25
        dtypes['fit_maxiter'] = int
        descr['fit_maxiter'] = 'Maximum number of rejection iterations during edge fitting.'

        defaults['fit_niter'] = 9
        dtypes['fit_niter'] = int
        descr['fit_niter'] = 'Number of iterations of re-measuring and re-fitting the edge ' \
                             'data; see :func:`pypeit.new_trace.fit_trace`.'

        dtypes['pca_n'] = int
        descr['pca_n'] = 'The number of PCA components to keep, which must be less than the ' \
                         'number of detected traces.  If not provided, determined by ' \
                         'calculating the minimum number of components required to explain a ' \
                         'given percentage of variance in the edge data; see `pca_var_percent`.'
            
        defaults['pca_var_percent'] = 99.8
        dtypes['pca_var_percent'] = [int, float]
        descr['pca_var_percent'] = 'The percentage (i.e., not the fraction) of the variance in ' \
                                   'the edge data accounted for by the PCA used to truncate ' \
                                   'the number of PCA coefficients to keep (see `pca_n`).  ' \
                                   'Ignored if `pca_n` is provided directly.'
        
        defaults['pca_function'] = 'polynomial'
        dtypes['pca_function'] = str
        options['pca_function'] = EdgeTracePar.valid_functions()
        descr['pca_function'] = 'Type of function fit to the PCA coefficients for each ' \
                                'component.  Options are: {0}'.format(
                                    ', '.join(options['pca_function']))
        
        defaults['pca_order'] = 3
        dtypes['pca_order'] = int
        descr['pca_order'] = 'Order of the function fit to the PCA coefficients.'
        
        defaults['pca_sigrej'] = [2., 2.]
        dtypes['pca_sigrej'] = [int, float, list]
        descr['pca_sigrej'] = 'Sigma rejection threshold for fitting PCA components. Individual ' \
                              'numbers are used for both lower and upper rejection. A list of ' \
                              'two numbers sets these explicitly (e.g., [2., 3.]).'

        defaults['pca_maxrej'] = 1
        dtypes['pca_maxrej'] = int
        descr['pca_maxrej'] = 'Maximum number of PCA coefficients rejected during a given fit ' \
                              'iteration.'

        defaults['pca_maxiter'] = 25
        dtypes['pca_maxiter'] = int
        descr['pca_maxiter'] = 'Maximum number of rejection iterations when fitting the PCA ' \
                               'coefficients.'

        defaults['smash_range'] = [0., 1.]
        dtypes['smash_range'] = list
        descr['smash_range'] = 'Range of the slit in the spectral direction (in fractional ' \
                               'units) to smash when searching for slit edges.  If the ' \
                               'spectrum covers only a portion of the image, use that range.'

        dtypes['edge_detect_clip'] = [int, float]
        descr['edge_detect_clip'] = 'Sigma clipping level for peaks detected in the collapsed, ' \
                                    'Sobel-filtered significance image.'

        defaults['trace_median_frac'] = 0.01
        dtypes['trace_median_frac'] = [int, float]
        descr['trace_median_frac'] = 'After detection of peaks in the rectified Sobel-filtered ' \
                                     'image and before refitting the edge traces, the rectified ' \
                                     'image is median filtered with a kernel width of ' \
                                     '`trace_median_frac*nspec` along the spectral dimension.'
        
        defaults['trace_thresh'] = 10.
        dtypes['trace_thresh'] = [int, float]
        descr['trace_thresh'] = 'After rectification and median filtering of the Sobel-filtered ' \
                                'image (see `trace_median_frac`), values in the median-filtered ' \
                                'image *below* this threshold are masked in the refitting of ' \
                                'the edge trace data.'

        defaults['fwhm_uniform'] = 3.0
        dtypes['fwhm_uniform'] = [int, float]
        descr['fwhm_uniform'] = 'The `fwhm` parameter to use when using uniform weighting in ' \
                                ':func:`fit_trace` when refining the PCA predictions of edges.' \
                                'See description :func:`peak_trace`.'

        defaults['niter_uniform'] = 9
        dtypes['niter_uniform'] = int
        descr['niter_uniform'] = 'The number of iterations of :func:`fit_trace` to use when ' \
                                 'using uniform weighting.'

        defaults['fwhm_gaussian'] = 3.0
        dtypes['fwhm_gaussian'] = [int, float]
        descr['fwhm_gaussian'] = 'The `fwhm` parameter to use when using Gaussian weighting in ' \
                                 ':func:`fit_trace` when refining the PCA predictions of edges.' \
                                 'See description :func:`peak_trace`.'

        defaults['niter_gaussian'] = 6
        dtypes['niter_gaussian'] = int
        descr['niter_gaussian'] = 'The number of iterations of :func:`fit_trace` to use when ' \
                                  'using Gaussian weighting.'

        defaults['det_buffer'] = 5
        dtypes['det_buffer'] = int
        descr['det_buffer'] = 'The minimum separation between the detector edges and a slit ' \
                              'edge for any added edge traces.  Must be positive.'

        defaults['max_nudge'] = 100
        dtypes['max_nudge'] = int
        descr['max_nudge'] = 'If parts of the (predicted) trace fall off the detector edge, ' \
                             'allow them to be nudged away from the detector edge up to and ' \
                             'including this maximum number of pixels. Can be 0 or larger.'

        defaults['sync_trace'] = 'pca'
        options['sync_trace'] = EdgeTracePar.valid_trace_modes()
        dtypes['sync_trace'] = str
        descr['sync_trace'] = 'Mode to use when predicting the form of the trace to insert.  ' \
                              'Use `pca` to use the PCA decomposition or `nearest` to reproduce ' \
                              'the shape of the nearest trace.'
                      
        defaults['sync_center'] = 'median'
        options['sync_center'] = EdgeTracePar.valid_center_modes()
        dtypes['sync_center'] = str
        descr['sync_center'] = 'Mode to use for determining the location of traces to insert.  ' \
                               'Use `median` to use the median of the matched left and right ' \
                               'edge pairs or `nearest` to ue the length of the nearest slit`.'

        defaults['sync_to_edge'] = True
        dtypes['sync_to_edge'] = bool
        descr['sync_to_edge'] = 'If adding a first left edge or a last right edge, ignore ' \
                                '`center_mode` for these edges and place them at the edge of ' \
                                'the detector (with the relevant shape).'

        defaults['min_slit_gap'] = 1.
        dtypes['min_slit_gap'] = [int, float]
        descr['min_slit_gap'] = 'Minimum allowed gap in pixels between the mean spatial ' \
                                'location of left and right edges.'

#        defaults['minimum_slit_length'] = 6.
        dtypes['minimum_slit_length'] = [int, float]
        descr['minimum_slit_length'] = 'Minimum slit length in arcsec.  Shorter slits are ' \
                                       'masked or clipped.'

#        defaults['length_range'] = 0.3
        dtypes['length_range'] = [int, float]
        descr['length_range'] = 'Range in relative slit length.  For example, a value of 0.3 ' \
                                'means that slit lengths should not vary more than 30%.  ' \
                                'Relatively shorter or longer slits are masked or clipped.'

        defaults['clip'] = True
        dtypes['clip'] = bool
        descr['clip'] = 'Instead of just masking bad slit trace edges, remove them.'

#        # Force trim to be a tuple
#        if pars['trim'] is not None and not isinstance(pars['trim'], tuple):
#            try:
#                pars['trim'] = tuple(pars['trim'])
#            except:
#                raise TypeError('Could not convert provided trim to a tuple.')
#        defaults['trim'] = (0,0)
#        dtypes['trim'] = tuple
#        descr['trim'] = 'How much to trim off each edge of each slit.  Each number should be 0 ' \
#                        'or positive'

        defaults['pad'] = 0
        dtypes['pad'] = int
        descr['pad'] = 'Integer number of pixels to consider beyond the slit edges.'

        dtypes['add_slits'] = [str, list]
        descr['add_slits'] = 'Add one or more user-defined slits.  The syntax to define a ' \
                             'slit to add is: \'det:spec:spat_left:spat_right\' where ' \
                             'det=detector, spec=spectral pixel, spat_left=spatial pixel of ' \
                             'left slit boundary, and spat_righ=spatial pixel of right slit ' \
                             'boundary.  For example, \'2:2000:2121:2322,3:2000:1201:1500\' ' \
                             'will add a slit to detector 2 passing through spec=2000 ' \
                             'extending spatially from 2121 to 2322 and another on detector 3 ' \
                             'at spec=2000 extending from 1201 to 1500.'

        dtypes['rm_slits'] = [str, list]
        descr['rm_slits'] = 'Remove one or more user-specified slits.  The syntax used to ' \
                            'define a slit to remove is: \'det:spec:spat\' where det=detector, ' \
                            'spec=spectral pixel, spat=spatial pixel.  For example, ' \
                            '\'2:2000:2121,3:2000:1500\' will remove the slit on detector 2 ' \
                            'that contains pixel (spat,spec)=(2000,2121) and on detector 3 ' \
                            'that contains pixel (2000,2121).'

        # Instantiate the parameter set
        super(EdgeTracePar, self).__init__(list(pars.keys()), values=list(pars.values()),
                                           defaults=list(defaults.values()),
                                           options=list(options.values()),
                                           dtypes=list(dtypes.values()),
                                           descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = ['filt_iter', 'sobel_mode', 'edge_thresh', 'follow_span', 'minimum_spec_length',
                   'valid_flux_thresh', 'max_spat_shift', 'max_spat_error', 'match_tol',
                   'fit_function', 'fit_order', 'fit_maxdev', 'fit_maxiter', 'fit_niter', 'pca_n',
                   'pca_var_percent', 'pca_function', 'pca_order', 'pca_sigrej', 'pca_maxrej',
                   'pca_maxiter', 'smash_range', 'edge_detect_clip', 'trace_median_frac',
                   'trace_thresh', 'fwhm_uniform', 'niter_uniform', 'fwhm_gaussian',
                   'niter_gaussian', 'det_buffer', 'max_nudge', 'sync_trace', 'sync_center',
                   'sync_to_edge', 'minimum_slit_length', 'length_range', 'clip', 'pad',
                   'add_slits', 'rm_slits']
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_functions():
        """
        Return the list of valid functions to use for slit tracing.
        """
        return [ 'polynomial', 'legendre', 'chebyshev' ]

    @staticmethod
    def valid_sobel_modes():
        """Return the valid sobel modes."""
        return [ 'nearest', 'constant' ]

    @staticmethod
    def valid_trace_modes():
        """Return the valid trace prediction modes."""
        return ['pca', 'nearest']

    @staticmethod
    def valid_center_modes():
        """Return the valid center prediction modes."""
        return ['median', 'nearest']

    def validate(self):
        pass


class EdgeTraceSet(masterframe.MasterFrame):
    r"""
    Core class that identifies, traces, and pairs edges in an image
    to define the slit apertures.

    

    TODO: MORE

    trace_img, mask, and det should be considered mutually exclusive compare to load

    load takes precedence.  I.e., if both trace_img and load are provided, trace_img is ignored!

    if trace_img is provided, the initialization will also run
    :func:`initial_trace` or :func:`auto_trace`, depending on the
    value of `auto`, *and* save the output.

    Nominal run:
        - initial_trace
        - moment_refine
        - fit_refine (calls fit_trace)
        - build_pca
        - peak_refine (calls peak_trace, which uses both pca_trace and fit_trace)

    Final trace is based on a run of fit_refine that pass through the detected peaks

    design and object data are only available if the spectrograph
    class has a get_slitmask function, and that the slit mask data
    includes the object information.

    TODO: Write a script that uses the empty function to document the
    data model for these two tables. And for EdgeTraceSet as a whole?

    TODO: Talk about loading; PCA data not currently saved; rebuilt if loaded

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The object that sets the instrument used to take the
            observations. Used to set :attr:`spectrograph`.
        par (:class:`pypeit.par.pypeitpar.EdgeTracePar`):
            The parameters used to guide slit tracing
        master_key (:obj:`str`, optional):
            The string identifier for the instrument configuration.  See
            :class:`pypeit.masterframe.MasterFrame`.
        master_dir (:obj:`str`, optional):
            Path to master frames.
        qa_path (:obj:`str`, optional):
            Directory for QA output.
        trace_img (`numpy.ndarray`_,
            :class:`pypeit.traceimage.TraceImage`, optional):
            Two-dimensional image used to trace slit edges. If a
            :class:`pypeit.traceimage.TraceImage` is provided, the
            raw files used to construct the image are saved.
        mask (`numpy.ndarray`_, optional):
            Mask for the trace image. Must have the same shape as
            `trace_img`. If None, all pixels are assumed to be valid.
        det (:obj:`int`, optional):
            The 1-indexed detector number that provided the trace
            image. This is *only* used to determine whether or not
            bad columns in the image are actually along columns or
            along rows, as determined by :attr:`spectrograph` and the
            result of a call to
            :func:`pypeit.spectrograph.Spectrograph.raw_is_transposed`.
        binning (`str`, optional):
            Comma-separated binning along the spectral and spatial
            directions following the PypeIt convention (e.g., '2,1').
            This is used to set the pixel scale of the image in
            arcsec per pixel, as needed for some assessments of the
            edge traces.
        auto (:obj:`bool`, optional):
            If a trace image is provided (`trace_img`), run
            :func:`auto_trace` instead of :func:`initial_trace`.
        load (:obj:`bool`, optional):
            Attempt to load existing output. If True and the file
            does not exist, an error is raised.

    Attributes:
        spectrograph
            (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            See argument list.
        par (:class:`pypeit.par.pypeitpar.EdgeTracePar`):
            See argument list.
        files (:obj:`list`):
            The list of raw files used to constuct the image used to
            detect and trace slit edges. Only defined if argument
            `trace_img` in :func:`initial_trace` or
            :func:`auto_trace` is a
            :class:`pypeit.traceimage.TraceImage` object.
        trace_img (`numpy.ndarray`_):
            Image data used to trace slit edges.
        trace_msk (`numpy.ndarray`_):
            A boolean array with the mask for pixels to ignore in the
            trace image (True is bad).
        det (:obj:`int`):
            1-indexed detector number.
        binning (:obj:`str`):
            Comma-separated binning along the spectral and spatial
            directions following the PypeIt convention (e.g., '2,1').
        sobel_sig (`numpy.ndarray`_)):
            Sobel-filtered image used to detect left and right edges
            of slits.
        sobel_sig_left (`numpy.ndarray`_):
            Lazy-loaded version of `sobel_sig` that clips the
            features related to right edges. Only kept for
            convenience.
        sobel_sig_right (`numpy.ndarray`_):
            Lazy-loaded version of `sobel_sig` that clips the
            features related to left edges. Only kept for
            convenience.
        nspec (:obj:`int`):
            Number of spectral pixels (rows) in the trace image
            (`axis=0`).
        nspat (:obj:`int`):
            Number of spatial pixels (columns) in the trace image
            (`axis=1`).
        traceid (`numpy.ndarray`_):
            The list of unique trace IDs.
        spat_img (`numpy.ndarray`_):
            An integer array with the spatial pixel nearest to each
            trace edge. This is identically::

                self.spat_img = np.round(self.spat_cen
                                         if self.spat_fit is None
                                         else self.spat_fit).astype(int)

        spat_cen (`numpy.ndarray`_):
            A floating-point array with the location of the slit edge
            for each spectral pixel *as measured* from the trace
            image. Shape is :math:`(N_{\rm spec},N_{\rm trace})`.
        spat_err (`numpy.ndarray`_):
            Error in slit edge locations; measurments without errors
            have their errors set to -1.
        spat_msk (`numpy.ndarray`_):
            An integer array with the mask bits assigned to each
            trace centroid; see :class:`TraceBitMask`.
        spat_fit (`numpy.ndarray`_):
            A model fit to the `spat_cen` data.
        spat_fit_type (:obj:`str`):
            An informational string identifier for the type of model
            used to fit the trace data.
        pca (:class:`EdgeTracePCA`):
            Result of a PCA decomposition of the edge traces and used
            to predict new traces.
        pca_type (:obj:`str`)
            An informational string indicating which data were used
            in the PCA decomposition, 'center' for `spat_cen` or
            'fit' for `spat_fit`.
        design (`astropy.table.Table`_):
            Collated slit-mask design data matched to the edge
            traces.
        objects (`numpy.recarray`_):
            Collated object ID and coordinate information matched to
            the design table.
        qa_path (:obj:`str`):
            Root path for QA output files.
        log (:obj:`list`):
            A list of strings indicating the main methods applied
            when tracing.
    """
    master_type = 'Trace'   # For MasterFrame base
    bitmask = TraceBitMask()            # Object used to define and toggle tracing mask bits
    def __init__(self, spectrograph, par, master_key=None, master_dir=None, qa_path=None,
                 trace_img=None, mask=None, det=1, binning=None, auto=False, load=False):

        # TODO: It's possible for the master key and the detector
        # number to be inconsistent...
        masterframe.MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
                                         master_key=master_key)

        self.spectrograph = spectrograph    # Spectrograph used to take the data
        self.par = par                      # Parameters used for slit edge tracing

        self.files = None               # Files used to construct the trace image
        self.trace_img = None           # The image used to find the slit edges
        self.trace_msk = None           # Mask for the trace image
        # TODO: Need a separate mask for the sobel image?
        self.det = None                 # Detector used for the trace image
        self.binning = None             # Detector ordered spectral then spatial
        self.sobel_sig = None           # Sobel filtered image used to detect edges
        self.sobel_sig_left = None      # Sobel filtered image used to trace left edges
        self.sobel_sig_right = None     # Sobel filtered image used to trace right edges
        self.nspec = None               # The shape of the trace image is (nspec,nspat)
        self.nspat = None

        self.traceid = None             # The ID numbers for each trace
        self.spat_img = None            # (Integer) Pixel nearest the slit edge for each trace
        self.spat_cen = None            # (Floating-point) Spatial coordinate of the slit edges
                                        # for each spectral pixel
        self.spat_err = None            # Error in the slit edge spatial coordinate
        self.spat_msk = None            # Mask for the slit edge position for each spectral pixel

        self.spat_fit = None            # The result of modeling the slit edge positions
        self.spat_fit_type = None       # The type of fitting performed
        
        self.pca = None                 # EdgeTracePCA with the PCA decomposition
        self.pca_type = None            # Measurements used to construct the PCA (center or fit)

        self.design = None              # Table that collates slit-mask design data matched to
                                        # the edge traces
        self.objects = None             # Table that collates object information, if available
                                        # in the slit-mask design, matched to the `design` table.

        self.qa_path = qa_path          # Directory for QA output

        self.log = None                 # Log of methods applied

        if trace_img is not None and load:
            raise ValueError('Arguments trace_img and load are mutually exclusive.  Choose to '
                             'either trace a new image or load a previous trace.')

        if load:
            # Attempt to load an existing master frame
            self.load()
        elif trace_img is not None:
            # Provided a trace image so instantiate the object.
            if auto:
                self.auto_trace(trace_img, mask=mask, det=det, binning=binning)
            else:
                self.initial_trace(trace_img, mask=mask, det=det, binning=binning)

    @property
    def file_path(self):
        """
        Overwrite MasterFrame default to force the file to be gzipped.
        """
        # TODO: Change the MasterFrame default to a compressed file?
        return '{0}.gz'.format(os.path.join(self.master_dir, self.file_name))

    @property
    def ntrace(self):
        """
        The number of edges (left and right) traced.
        """
        return self.traceid.size

    @staticmethod
    def empty_design_table(rows=None):
        """
        Construct an empty `design` table.

        Args:
            rows (:obj:`int`, optional):
                Number of rows for each column. If None, the table
                has empty columns.

        Returns:
            `astropy.table.Table`_: Instance of the empty design
            table.
        """
        length = 0 if rows is None else rows
        return table.Table([
                    table.Column(name='TRACEID', dtype=int, length=length,
                                 description='Trace ID Number'),
                    table.Column(name='TRACESROW', dtype=int, length=length,
                                 description='Spectral row for provided left and right edges.'),
                    table.Column(name='TRACELPIX', dtype=int, length=length,
                                 description='Spatial pixel coordinate for left edge'),
                    table.Column(name='TRACERPIX', dtype=int, length=length,
                                 description='Spatial pixel coordinate for right edge'),
                    table.Column(name='SLITID', dtype=int, length=length,
                                 description='Slit ID Number'),
                    table.Column(name='SLITLFOC', dtype=int, length=length,
                                 description='Left edge of the slit in mm at the focal plane'),
                    table.Column(name='SLITRFOC', dtype=int, length=length,
                                 description='Right edge of the slit in mm at the focal plane'),
                    table.Column(name='SLITRA', dtype=int, length=length,
                                 description='Right ascension of the slit center (deg)'),
                    table.Column(name='SLITDEC', dtype=int, length=length,
                                 description='Declination of the slit center (deg)'),
                    table.Column(name='SLITLEN', dtype=int, length=length,
                                 description='Slit length (arcsec)'),
                    table.Column(name='SLITWID', dtype=int, length=length,
                                 description='Slit width (arcsec)'),
                    table.Column(name='SLITPA', dtype=int, length=length,
                                 description='Slit position angle onsky (deg from N through E)')
                           ])

    @staticmethod
    def empty_objects_table(rows=None):
        """
        Construct an empty `objects` table.

        Args:
            rows (:obj:`int`, optional):
                Number of rows for each column. If None, the table
                has empty columns.

        Returns:
            `astropy.table.Table`_: Instance of the empty object
            table.
        """
        length = 0 if rows is None else rows
        return table.Table([
                    table.Column(name='OBJID', dtype=int, length=length,
                                 description='Object ID Number'),
                    table.Column(name='OBJRA', dtype=int, length=length,
                                 description='Right ascension of the object (deg)'),
                    table.Column(name='OBJDEC', dtype=int, length=length,
                                 description='Declination of the object (deg)'),
                    table.Column(name='SLITID', dtype=int, length=length,
                                 description='Slit ID Number'),
                    table.Column(name='SLITINDX', dtype=int, length=length,
                                 description='Row index of relevant slit in the design table')
                           ])

    def rectify(self, flux, mask=None, extract_width=None, mask_threshold=0.5):
        r""""
        Rectify the provided image based on the current edge trace
        PCA model.

        The is primarily a wrapper for :func:`rectify_image`; see its
        documentation for more detail.

        Args:
            flux (`numpy.ndarray`_):
                The 2D image to rectify. Its shape should match the
                image used to construct the edge traces:
                :math:`(N_{\rm spec}, N_{\rm spat})`.
            mask (`numpy.ndarray`_, optional):
                Boolean mask for pixels to ignore in input image. If
                None, no pixels are masked in the rectification. If
                provided, shape must match `flux`.
            extract_width (:obj:`float`, optional):
                The width of the extraction aperture to use for the
                image rectification. When using extraction to rectify
                the image, flux conservation is not as accurate. If
                None, the image recification is performed using
                :class:`pypeit.sampling.Resample` along each row.
            mask_threshold (:obj:`float`, optional):
                Either due to `mask` or the bounds of the provided
                `flux`, pixels in the rectified image may not be fully
                covered by valid pixels in `flux`. Pixels in the
                output image with less than this fractional coverage
                by input pixels are flagged in the output.

        Returns:
             Two `numpy.ndarray`_ objects are returned both with
             shape :math:`(N_{\rm spec}, N_{\rm spat})`, the rectified
             image and its boolean mask.
        """
        if self.pca is None:
            raise ValueError('Must first run the PCA analysis for the traces; run build_pca.')

        # Get the traces that cross the reference row at the first and last pixels of the image
        first_last_trace = self.pca.predict(np.array([0,self.nspat-1]))
        # Use these two traces to define the spatial pixel coordinates to sample
        start = np.ceil(np.amax(np.amin(first_last_trace, axis=1))).astype(int)
        buffer = self.nspat - np.floor(np.amin(np.amax(first_last_trace, axis=1))).astype(int) \
                    + start
        # Rectify the image
        ocol = np.arange(self.nspat+buffer)-start
        return rectify_image(flux, self.pca.predict(ocol), mask=mask, ocol=ocol,
                             max_ocol=self.nspat-1, extract_width=extract_width,
                             mask_threshold=mask_threshold)

    def auto_trace(self, trace_img, mask=None, det=1, binning=None, save=True, debug=False):
        r"""
        Execute a fixed series of methods to automatically identify
        and trace slit edges.

        Args:
            trace_img (`numpy.ndarray`_, :class:`pypeit.traceimage.TraceImage`):
                2D image used to trace slit edges. If a
                :class:`pypeit.traceimage.TraceImage` is provided,
                the raw files used to construct the image and on-chip
                binning are saved; the latter overrides any directly
                provided `binning`. The array should have shape
                :math:`(N_{\rm spec},N_{\rm spat})`; i.e., spectra
                are ordered along columns.
            mask (`numpy.ndarray`_, optional):
                Mask for the trace image. Must have the same shape as
                `trace_img`. If None, all pixels are assumed to be
                valid.
            det (:obj:`int`, optional):
                The 1-indexed detector number that provided the trace
                image. This is *only* used to determine whether or
                not bad columns in the image are actually along
                columns or along rows, as determined by
                :attr:`spectrograph` and the result of a call to
                :func:`pypeit.spectrograph.Spectrograph.raw_is_transposed`.
            binning (:obj:`str`, optional):
                On-detector binning of the data ordered spectral then
                spatial with format, e.g., `2,1`. Ignored if
                `trace_img` is an instance of
                :class:`pypeit.traceimage.TraceImage`.
            save (:obj:`bool`, optional):
                Save the result to the master frame.
            debug (:obj:`bool`, optional):
                Run in debug mode.
        """
        self.initial_trace(trace_img, mask=mask, det=det, binning=binning, save=False)
        self.moment_refine()
        self.fit_refine(debug=debug)
        self.pca_refine(debug=debug)
        self.peak_refine(rebuild_pca=True, debug=debug)
        self.sync()
        self.log += [inspect.stack()[0][3]]
        if save:
            self.save()

    def initial_trace(self, trace_img, mask=None, det=1, binning=None, save=True):
        r"""
        Initialize the object for tracing a new image.

        This effectively reinstantiates the object and must be the
        first method called for tracing an image.  The algorithm:
            - Lightly boxcar smooths the trace image spectrally.
            - Replaces bad pixel columns, if a mask is provided.
            - Applies a Sobel filter to the trace image along columns
              to detect slit edges using steep positive gradients
              (left edges) and steep negative gradients (right
              edges). See :func:`detect_slit_edges`.
            - Follows the detected left and right edges along
              spectrally adjacent pixels to identify coherent traces.
              See :func:`identify_traces`.
            - Performs basic handling of orphaned left or right
              edges. See :func:`handle_orphan_edges`.
            - Initializes the attributes that provide the trace
              position for each spectral pixel based on these
              results.

        Used parameters from :attr:`par` (:class:`EdgeTracePar`) are
        `filt_iter`, `sobel_mode`, `edge_thresh`,
        `minimum_spec_length`, `follow_span`, and
        `valid_flux_thresh`.

        The results of this are, by default, saved to the master
        frame; see `save` argument and :func:`save`.

        Args:
            trace_img (`numpy.ndarray`_, :class:`pypeit.traceimage.TraceImage`):
                2D image used to trace slit edges. If a
                :class:`pypeit.traceimage.TraceImage` is provided,
                the raw files used to construct the image and on-chip
                binning are saved; the latter overrides any directly
                provided `binning`. The array should have shape
                :math:`(N_{\rm spec},N_{\rm spat})`; i.e., spectra
                are ordered along columns.
            mask (`numpy.ndarray`_, optional):
                Mask for the trace image. Must have the same shape as
                `trace_img`. If None, all pixels are assumed to be
                valid.
            det (:obj:`int`, optional):
                The 1-indexed detector number that provided the trace
                image. This is *only* used to determine whether or
                not bad columns in the image are actually along
                columns or along rows, as determined by
                :attr:`spectrograph` and the result of a call to
                :func:`pypeit.spectrograph.Spectrograph.raw_is_transposed`.
            binning (:obj:`str`, optional):
                On-detector binning of the data ordered spectral then
                spatial with format, e.g., `2,1`. Ignored if
                `trace_img` is an instance of
                :class:`pypeit.traceimage.TraceImage`.
            save (:obj:`bool`, optional):
                Save the result to the master frame.
        """
        # Parse the input based on its type
        if isinstance(trace_img, TraceImage):
            self.files = trace_img.files
            _trace_img = trace_img.stack

            # TODO: Waiting for changes in X's upcoming PR to know what
            # to do. CombinedImage class should have a single binning
            # value.
            _binning = trace_img.binning[0]
            if binning is not None:
                msgs.warn('Using binning of the trace image provided by the object directly, '
                          'not the value provided to this function.')
            # TODO: does TraceImage have a mask?
            # TODO: instead keep the TraceImage object instead of
            # deconstructing it...
        else:
            _trace_img = trace_img
            _binning = binning

        # Check the input
        if _trace_img.ndim != 2:
            raise ValueError('Trace image must be 2D.')
        self.trace_img = _trace_img
        self.nspec, self.nspat = self.trace_img.shape
        self.trace_msk = np.zeros((self.nspec, self.nspat), dtype=bool) if mask is None else mask
        if self.trace_msk.shape != self.trace_img.shape:
            raise ValueError('Mask is not the same shape as the trace image.')
        self.det = det
        # TODO: The default PypeIt format for setting `binning` for
        # unbinned data should be function. It shouldn't be hardcoded
        # in case it changes.
        self.binning = '1,1' if _binning is None else _binning

        # Lightly smooth the image before using it to trace edges
        _trace_img = ndimage.uniform_filter(self.trace_img, size=(3, 1), mode='mirror')

        # Replace bad-pixel columns if they exist
        # TODO: Do this before passing the image to this function?
        # Instead of just replacing columns, replace all bad pixels...
        if np.any(self.trace_msk):
            # Do we need to replace bad *rows* instead of bad columns?
            flip = self.spectrograph.raw_is_transposed(det=self.det)
            axis = 1 if flip else 0

            # Replace bad columns that cover more than half the image
            bad_cols = np.sum(self.trace_msk, axis=axis) > (self.trace_msk.shape[axis]//2)
            if flip:
                # Deal with the transposes
                _trace_img = procimg.replace_columns(_trace_img.T, bad_cols, copy=True,
                                                     replace_with='linear').T
            else:
                _trace_img = procimg.replace_columns(_trace_img, bad_cols, copy=True,
                                                     replace_with='linear')

        # Filter the trace image and use the filtered image to detect
        # slit edges
        # TODO: Decide if mask should be passed to this or not,
        # currently not because of issues when masked pixels happen to
        # land in slit gaps.
        self.sobel_sig, edge_img = detect_slit_edges(_trace_img,
                                                     median_iterations=self.par['filt_iter'],
                                                     sobel_mode=self.par['sobel_mode'],
                                                     sigdetect=self.par['edge_thresh'])

        # Empty out the images prepared for left and right tracing
        # until they're needed.
        self.sobel_sig_left = None
        self.sobel_sig_right = None

        # Identify traces by following the detected edges in adjacent
        # spectral pixels.
        minimum_spec_pix = self.nspec * self.par['minimum_spec_length']
        trace_id_img = identify_traces(edge_img, follow_span=self.par['follow_span'],
                                       minimum_spec_length=minimum_spec_pix)

        # Update the traces by handling single orphan edges and/or
        # traces without any left or right edges.
        flux_valid = np.median(_trace_img) > self.par['valid_flux_thresh']
        trace_id_img = handle_orphan_edge(trace_id_img, self.sobel_sig, mask=self.trace_msk,
                                          flux_valid=flux_valid, copy=True)

        # Set the ID image to a MaskedArray to ease some subsequent
        # calculations; pixels without a detected edge are masked.
        trace_id_img = np.ma.MaskedArray(trace_id_img, mask=trace_id_img == 0)

        # Find the set of trace IDs; left traces are negative, right
        # traces are positive
        self.traceid = np.unique(trace_id_img.compressed())

        # Initialize the mask bits for the trace coordinates and
        # initialize them all as having no edge
        self.spat_msk = np.zeros((self.nspec, self.ntrace), dtype=self.bitmask.minimum_dtype())
        self.spat_msk = self.bitmask.turn_on(self.spat_msk, 'NOEDGE')

        # Save the input trace edges and remove the mask for valid edge
        # coordinates
        self.spat_img = np.zeros((self.nspec, self.ntrace), dtype=int)
        for i in range(self.ntrace):
            row, col = np.where(np.invert(trace_id_img.mask)
                                    & (trace_id_img.data == self.traceid[i]))
            self.spat_img[row,i] = col
            self.spat_msk[row,i] = 0            # Turn-off the mask

        # Instantiate objects to store the floating-point trace
        # centroids and errors.
        self.spat_cen = self.spat_img.astype(float)   # This makes a copy
        self.spat_err = np.zeros((self.nspec, self.ntrace), dtype=float)

        # No fitting has been done yet
        self.spat_fit_type = None
        self.spat_fit = None

        # No PCA has been constructed yet
        self.pca = None
        self.pca_type = None

        # No design or object data
        self.design = None
        self.objects = None

        # Restart the log
        self.log = [inspect.stack()[0][3]]

        # Save if requested
        if save:
            self.save()

    def save(self, outfile=None, overwrite=True, checksum=True):
        """
        Save the trace object to a file for full recall.

        Args:
            outfile (:obj:`str`, optional):
                Name for the output file.  Defaults to
                :attr:`file_path`.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file.
            checksum (:obj:`bool`, optional):
                Passed to `astropy.io.fits.HDUList.writeto` to add
                the DATASUM and CHECKSUM keywords fits header(s).
        """
        _outfile = self.file_path if outfile is None else outfile
        # Check if it exists
        if os.path.exists(_outfile) and not overwrite:
            # TODO: Raise an error instead?
            msgs.warn('Master file exists: {0}'.format(_outfile) + msgs.newline()
                      + 'Set overwrite=True to overwrite it.')
            return

        # Build the primary header
        #   - Initialize with basic metadata
        prihdr = self.initialize_header()
        #   - Add the qa path
        prihdr['QADIR'] = (self.qa_path, 'PypeIt: QA directory')
        #   - Add metadata specific to this class
        prihdr['SPECT'] = (self.spectrograph.spectrograph, 'PypeIt: Spectrograph Name')
        #   - List the processed raw files, if available
        if self.files is not None:
            nfiles = len(self.files)
            ndig = int(np.log10(nfiles))+1
            for i in range(nfiles):
                prihdr['RAW{0}'.format(str(i+1).zfill(ndig))] \
                            = (self.files[i], 'PypeIt: Processed raw file')
        #   - Add the binning
        prihdr['BINNING'] = (self.binning, 'PypeIt: Binning')
        #   - Add the detector number
        prihdr['DET'] = (self.det, 'PypeIt: Detector')
        #   - Add the tracing parameters
        self.par.to_header(prihdr)
        #   - List the completed methods, if there are any
        if self.log is not None:
            ndig = int(np.log10(len(self.log)))+1
            for i,m in enumerate(self.log):
                prihdr['LOG{0}'.format(str(i+1).zfill(ndig))] \
                        = (m, '{0}: Completed method'.format(self.__class__.__name__))
        #   - PCA type, used for rebuilding the PCA when loading
        prihdr['PCATYPE'] = ('None' if self.pca is None else self.pca_type,
                             'PypeIt: Edge trace PCA type')
        #   - Indicate the type if fit (TODO: Keep the fit parameters?)
        fithdr = fits.Header()
        fithdr['FITTYP'] = 'None' if self.spat_fit_type is None else self.spat_fit_type

        # Only put the definition of the bits in the trace mask in the
        # header of the appropriate extension.
        mskhdr = fits.Header()
        self.bitmask.to_header(mskhdr)

        # Determine if the file should be compressed
        compress = False
        if _outfile.split('.')[-1] == 'gz':
            _outfile = _outfile[:_outfile.rfind('.')] 
            compress = True
    
        # Write the fits file; note not everything is written. Some
        # arrays are reconstruced by the load function.
        fits.HDUList([fits.PrimaryHDU(header=prihdr),
                      fits.ImageHDU(data=self.trace_img, name='TRACEIMG'),
                      fits.ImageHDU(data=self.trace_msk.astype(np.int16), name='TRACEMSK'),
                      fits.ImageHDU(data=self.sobel_sig, name='SOBELSIG'),
                      fits.ImageHDU(data=self.traceid, name='TRACEID'),
                      fits.ImageHDU(data=self.spat_cen, name='CENTER'),
                      fits.ImageHDU(data=self.spat_err, name='CENTER_ERR'),
                      fits.ImageHDU(header=mskhdr, data=self.spat_msk, name='CENTER_MASK'),
                      fits.ImageHDU(header=fithdr, data=self.spat_fit, name='CENTER_FIT'),
                    ]).writeto(_outfile, overwrite=True, checksum=checksum)
        msgs.info('Master frame written to {0}'.format(_outfile))

        # Compress the file if the output filename has a '.gz'
        # extension; this is slow but still faster than if you have
        # astropy.io.fits do this directly
        if compress:
            msgs.info('Compressing file to: {0}.gz'.format(_outfile))
            io.compress_file(_outfile, overwrite=overwrite)

    def load(self, validate=True, rebuild_pca=True):
        """
        Load and reinitialize the trace data.

        Data is read from :attr:`file_path` and used to overwrite any
        internal data. Specific comparisons of the saved data are
        performed to ensure the file is consistent with having been
        written by an identical instantiation; see :func:`_reinit`.

        To load a full :class:`EdgeTraceSet` directly from a file,
        use :func:`from_file`.

        Args:
            validate (:obj:`bool`, optional):
                Validate that the spectrograph, parameter set, and
                bitmask have not changed between the current internal
                values and the values read from the fits file. The
                method raises an error if the spectrograph or
                parameter set are different. If the bitmask is
                different, a warning is issued and the bitmask
                defined by the header is used instead of the existing
                :attr:`bitmask`.
            rebuild_pca (:obj:`bool`, optional):
                If the primary header indicates that the PCA
                decompostion had been performed on the save object,
                use the saved parameters to rebuild that PCA.
        Raises:
            FileNotFoundError:
                Raised if no data has been written for this master
                frame.
            ValueError:
                Raised if validation of the data fails (actually
                raised by :func:`_reinit`).
        """
        filename = self.file_path
        # Check the file exists
        if not os.path.isfile(filename):
            raise FileNotFoundError('File does not exit: {0}'.format(filename))
        with fits.open(filename) as hdu:
            # Re-initialize and validate
            self._reinit(hdu, validate=validate, rebuild_pca=rebuild_pca)

    @classmethod
    def from_file(cls, filename, rebuild_pca=True):
        """
        Instantiate using data from a file.

        To reload data that has been saved for an existing
        instantiation, use :func:`load`.

        Args:
            filename (:obj:`str`):
                Fits file produced by :func:`save`.
            rebuild_pca (:obj:`bool`, optional):
                If the primary header indicates that the PCA
                decompostion had been performed on the save object,
                use the saved parameters to rebuild that PCA.
        """
        # Check the file exists
        if not os.path.isfile(filename):
            raise FileNotFoundError('File does not exit: {0}'.format(filename))
        msgs.info('Loading EdgeTraceSet data from: {0}'.format(filename))
        with fits.open(filename) as hdu:
            this = cls(load_spectrograph(hdu[0].header['SPECT']),
                       EdgeTracePar.from_header(hdu[0].header),
                       master_key=hdu[0].header['MSTRKEY'], master_dir=hdu[0].header['MSTRDIR'],
                       qa_path=hdu[0].header['QADIR'])

            # Re-initialize and validate
            this._reinit(hdu, rebuild_pca=rebuild_pca)
        return this

    def _reinit(self, hdu, validate=True, rebuild_pca=True):
        """
        Reinitialize the internals based on the provided fits HDU.

        Args:
            hdu (`astropy.io.fits.Header`):
                The fits data used to reinitialize the object written
                by :func:`save`.
            validate (:obj:`bool`, optional):
                Validate that the spectrograph, parameter set, and
                bitmask have not changed between the current internal
                values and the values read from the fits file. The
                method raises an error if the spectrograph or
                parameter set are different. If the bitmask is
                different, a warning is issued and the bitmask
                defined by the header is used instead of the existing
                :attr:`bitmask`.
            rebuild_pca (:obj:`bool`, optional):
                If the primary header indicates that the PCA
                decompostion had been performed on the save object,
                use the saved parameters to rebuild that PCA.
        """
        # Read and assign data from the fits file
        self.files = io.parse_hdr_key_group(hdu[0].header, prefix='RAW')
        if len(self.files) == 0:
            self.files = None
        self.trace_img = hdu['TRACEIMG'].data
        self.nspec, self.nspat = self.trace_img.shape
        self.trace_msk = hdu['TRACEMSK'].data.astype(bool)
        self.det = hdu[0].header['DET']
        self.binning = hdu[0].header['BINNING']
        self.sobel_sig = hdu['SOBELSIG'].data
        self.traceid = hdu['TRACEID'].data
        self.spat_cen = hdu['CENTER'].data
        self.spat_err = hdu['CENTER_ERR'].data
        self.spat_msk = hdu['CENTER_MASK'].data
        self.spat_fit = hdu['CENTER_FIT'].data
        self.spat_fit_type = None if hdu['CENTER_FIT'].header['FITTYP'] == 'None' \
                                else hdu['CENTER_FIT'].header['FITTYP']

        self.spat_img = np.round(self.spat_cen if self.spat_fit is None
                                 else self.spat_fit).astype(int)

        # Rebuild the PCA if it existed previously and requested
        self.pca_type = None if hdu[0].header['PCATYPE'] == 'None' else hdu[0].header['PCATYPE']
        self._reset_pca(rebuild_pca and self.pca_type is not None)

        self.log = io.parse_hdr_key_group(hdu[0].header, prefix='LOG')

        # TODO: Recalculate Sobel left and right images instead of
        # setting them to None?
        self.sobel_sig_left = None
        self.sobel_sig_right = None

        # Finished, if not validating
        if not validate:
            return

        # Test the bitmask has the same keys and key values
        hdr_bitmask = BitMask.from_header(hdu['CENTER_MASK'].header)
        if hdr_bitmask.bits != self.bitmask.bits:
            msgs.warn('The bitmask in this fits file appear to be out of date!  Will continue '
                      'by using old bitmask but errors may occur.  You should recreate this '
                      'master frame.')
            self.bitmask = hdr_bitmask

        # Test the spectrograph is the same
        if self.spectrograph.spectrograph != hdu[0].header['SPECT']:
            raise ValueError('Data used for this master frame was from a different spectrograph!')

        # Test the parameters used are the same
        par = EdgeTracePar.from_header(hdu[0].header)
        if self.par.data != par.data:
            # TODO: The above inequality works for non-nested ParSets,
            # but will need to be more careful for nested ones, or just
            # avoid writing nested ParSets to headers...
            raise ValueError('Parameters used to construct this master used different parameters!')

    def show(self, include_error=False, thin=1, in_ginga=False, include_img=False):
        """
        Show a scatter plot of the current trace data and fit, if
        it's available.

        Args:
            include_error (:obj:`bool`, optional):
                Show the errors on the measurements
            thin (:obj:`int`, optional):
                Thin the data plotted by plotting every `thin`
                measurement in the spectral direction. Default is to
                plot all data; to show every other datum, set
                `thin=2`.
            in_ginga (:obj:`bool`, optional):
                Show the trace against the trace image in a ginga
                viewer instead of a line and scatter plot.
            include_img (:obj:`bool`, optional):
                Overlay the trace data on the trace image.
        """
        if in_ginga:
            raise NotImplementedError('Someone want to do this?')

        # Show the traced image
        if include_img:
            plt.imshow(self.trace_img.T, origin='lower', interpolation='nearest', aspect='auto')
        # Spectral position
        spec = np.tile(np.arange(self.nspec), (self.ntrace,1)).T
        if include_error:
            plt.errorbar(spec[::thin,:], self.spat_cen[::thin,:], yerr=self.spat_err, fmt='none',
                         ecolor='k', elinewidth=0.5, alpha=0.3, capthick=0, zorder=3)
        left = self.traceid < 0
        plt.scatter(spec[::thin,left], self.spat_cen[::thin,left], marker='.', color='k', s=30,
                    lw=0, zorder=4, label='left edge measurements', alpha=0.8)
        right = self.traceid > 0
        plt.scatter(spec[::thin,right], self.spat_cen[::thin,right], marker='.', color='0.7',
                    s=30, lw=0, zorder=4, label='right edge measurements', alpha=0.8)
        if self.spat_fit is None:
            plt.legend()
            plt.show()
            return

        for i in range(self.ntrace):
            # If statement structure primarily for the labels. Only
            # difference between left and right is the color.
            if left[i]:
                left_line = plt.plot(spec[::thin,i], self.spat_fit[::thin,i], color='C3', lw=0.5,
                                     zorder=5)
            else:
                right_line = plt.plot(spec[::thin,i], self.spat_fit[::thin,i], color='C1', lw=0.5,
                                      zorder=5)
        left_line[0].set_label('left edge fit')
        right_line[0].set_label('right edge fit')
        plt.legend()
        plt.show()

    def qa_plot(self, fileroot=None, min_spat=20):
        """
        Build a series of QA plots showing the edge traces.

        Args:
            fileroot (:obj:`str`, optional):
                Root name for the output files. The number of output
                files depends on the layout and the number of traces
                found. If None, plots are displayed interactively.
            min_spat (:obj:`int`, optional):
                Minimum number of spectral pixels to plot for each
                trace. If None, set to twice the difference between
                the minimum and maximum centroid of the plotted
                trace.
        """
        # Restore matplotlib defaults
        # TODO: Is this going to screw up later plots?
        matplotlib.rcParams.update(matplotlib.rcParamsDefault)

        # Set font size
        rc('font', size=8)

        # Spectral pixel coordinate vector and global plot limits
        spec = np.arange(self.nspec)
        xlim = [-1,self.nspec]
        img_zlim = utils.growth_lim(self.trace_img, 0.95, fac=1.05)
        sob_zlim = utils.growth_lim(self.sobel_sig, 0.95, fac=1.05)

        # Set figure
        w,h = plt.figaspect(1)
        fig = plt.figure(figsize=(1.5*w,1.5*h))

        # Grid for plots
        n = np.array([2,3])
        buff = np.array([0.05, 0.03])
        strt = np.array([0.07, 0.04])
        end = np.array([0.99, 0.99])
        delt = (end-(n-1)*buff-strt)/n

        # Determine the number of plot pages
        npages = self.ntrace//int(np.prod(n))
        if npages * np.prod(n) < self.ntrace:
            npages += 1
        ndig = int(np.log10(npages))+1

        # Make plots
        j = 0
        page = 0
        msgs.info('Constructing Trace QA plots')
        for i in range(self.ntrace):

            # Plot index
            jj = j//n[0]
            ii = j - jj*n[0]

            # Plot coordinates
            ax_x = strt[0]+ii*(buff[0]+delt[0])
            ax_y0 = strt[1]+(n[1]-jj-1)*(buff[1]+delt[1])

            # Spatial pixel plot limits for this trace
            indx = np.invert(self.bitmask.flagged(self.spat_msk[:,i],
                                                  flag=self.bitmask.bad_flags()))
            ylim = utils.growth_lim(self.spat_cen[indx,i], 1.0, fac=2.0)
            if min_spat is not None and np.diff(ylim) < min_spat:
                ylim = np.sum(ylim)/2 + np.array([-1,1])*min_spat/2

            # Plot the trace image and the fit (if it exists)
            ax = fig.add_axes([ax_x, ax_y0 + 2*delt[1]/3, delt[0], delt[1]/3.])
            ax.minorticks_on()
            ax.tick_params(which='major', length=10, direction='in', top=True, right=True)
            ax.tick_params(which='minor', length=5, direction='in', top=True, right=True)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.xaxis.set_major_formatter(ticker.NullFormatter())
            ax.imshow(self.trace_img.T, origin='lower', interpolation='nearest', vmin=img_zlim[0],
                      vmax=img_zlim[1], aspect='auto')
            if self.spat_fit is not None:
                ax.plot(spec, self.spat_fit[:,i], color='C3' if self.traceid[i] < 0 else 'C1')
            ax.text(0.95, 0.8, 'Trace {0}'.format(self.traceid[i]), ha='right', va='center',
                    transform=ax.transAxes, fontsize=12)

            # Plot the filtered image and the fit (if it exists)
            ax = fig.add_axes([ax_x, ax_y0 + delt[1]/3, delt[0], delt[1]/3.])
            ax.minorticks_on()
            ax.tick_params(which='major', length=10, direction='in', top=True, right=True)
            ax.tick_params(which='minor', length=5, direction='in', top=True, right=True)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.xaxis.set_major_formatter(ticker.NullFormatter())
            ax.imshow(self.sobel_sig.T, origin='lower', interpolation='nearest', vmin=sob_zlim[0],
                      vmax=sob_zlim[1], aspect='auto')
            if self.spat_fit is not None:
                ax.plot(spec, self.spat_fit[:,i], color='C3' if self.traceid[i] < 0 else 'C1')
            if ii == 0:
                ax.text(-0.13, 0.5, 'Spatial Coordinate (pix)', ha='center', va='center',
                        transform=ax.transAxes, rotation='vertical')

            # Plot the trace centroids and the fit (if it exists)
            ax = fig.add_axes([ax_x, ax_y0, delt[0], delt[1]/3.])
            ax.minorticks_on()
            ax.tick_params(which='major', length=10, direction='in', top=True, right=True)
            ax.tick_params(which='minor', length=5, direction='in', top=True, right=True)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.scatter(spec[indx], self.spat_cen[indx,i], marker='.', s=50, color='k', lw=0)
            nindx = np.invert(indx)
            if np.any(nindx):
                ax.scatter(spec[nindx], self.spat_cen[nindx,i], marker='x', s=30, color='0.5',
                           lw=0.5)
            if self.spat_fit is not None:
                ax.plot(spec, self.spat_fit[:,i], color='C3' if self.traceid[i] < 0 else 'C1')
            if jj == n[1]-1:
                ax.text(0.5, -0.3, 'Spectral Coordinate (pix)', ha='center', va='center',
                        transform=ax.transAxes)

            # Prepare for the next trace plot
            j += 1
            if j == np.prod(n) or i == self.ntrace-1:
                j = 0
                if fileroot is None:
                    plt.show()
                else:
                    page += 1
                    ofile = os.path.join(self.qa_path,
                                         '{0}_{1}.png'.format(fileroot, str(page).zfill(ndig)))
                    fig.canvas.print_figure(ofile, bbox_inches='tight')
                    msgs.info('Finished page {0}/{1}'.format(page, npages))
                fig.clear()
                plt.close(fig)
                fig = plt.figure(figsize=(1.5*w,1.5*h))

    def _side_dependent_sobel(self, side):
        """
        Return the Sobel-filtered image relevant to tracing the given
        side.

        The calculation of the side-dependent Sobel image should only
        need to be done once per side per instantiation. Unless they
        are reset to None (such as when the object is
        reinstantiated), multiple calls to this function allow for
        the data to be "lazy loaded" by performing the calculation
        once and then keeping the result in memory.

        Args:
            side (:obj:`str`):
                The side to return; must be 'left' or 'right'
                (case-sensitive).
    
        Returns:
            `numpy.ndarray`_: The manipulated Sobel image relevant to
            tracing the specified edge side.
        """
        # TODO: Add boxcar to TraceSlitsPar
        boxcar = 5
        if side == 'left':
            if self.sobel_sig_left is None:
                self.sobel_sig_left = prepare_sobel_for_trace(self.sobel_sig, boxcar=boxcar,
                                                              side='left')
            return self.sobel_sig_left
        if side == 'right':
            if self.sobel_sig_right is None:
                self.sobel_sig_right = prepare_sobel_for_trace(self.sobel_sig, boxcar=boxcar,
                                                               side='right')
            return self.sobel_sig_right

    def moment_refine(self, maxshift_start=0.5, continuous=False):
        """
        Refine the edge positions using a moment analysis and assess
        the results.

        For each set of edges (left and right), this method uses
        :func:`follow_trace_moment` to refine the centroids of the
        currently identified traces. The resulting traces are then
        checked that they cover at least a minimum fraction of the
        detector, whether or not they hit the detector edge, and
        whether or not they cross one another; see
        :func:`check_traces`. These two operations are done
        iteratively until all input traces are either refined or
        flagged for deletion.

        Nominally, this method should be run directly after
        :func:`initial_trace`.

        Used parameters from :attr:`par` (:class:`EdgeTracePar`) are
        `fwhm_uniform`, `max_spat_shift`, `max_spat_error`, and
        `clip`.

        .. warning::
            - This function modifies the internal trace arrays **in
              place**.
            - Because this changes :attr:`spat_cen` and
              :attr:`spat_err`, any model fitting of these data are
              erased by this function! I.e., :attr:`spat_fit` and
              :attr:`spat_fit_type` are set to None.
            - This *always* removes the PCA if it exists.

        Args:
            maxshift_start (:obj:`float`, optional):
                Maximum shift in pixels allowed for the adjustment of
                the first row analyzed, which is the row that has the
                most slit edges that cross through it.
            continuous (:obj:`bool`, optional):
                Keep only the continuous part of the traces from the
                starting row.
        """
        # Parse parameters and report
        width = 2 * self.par['fwhm_uniform']
        maxshift_follow = self.par['max_spat_shift']
        maxerror = self.par['max_spat_error']

        msgs.info('Width of window for centroiding the edges: {0:.1f}'.format(width))
        msgs.info('Max shift between spectrally adjacent pixels: {0:.2f}'.format(maxshift_follow))
        msgs.info('Max centroid error: {0:.2f}'.format(maxerror))
    
        # To improve performance, generate bogus ivar and mask once
        # here so that they don't have to be generated multiple times.
        # TODO: Keep these as work space as class attributes?
        ivar = np.ones_like(self.sobel_sig, dtype=float)
        _mask = np.zeros_like(self.sobel_sig, dtype=bool) \
                    if self.trace_msk is None else self.trace_msk
        fwgt = np.ones_like(self.sobel_sig, dtype=float)

        # Book-keeping objects to keep track of which traces have been
        # analyzed and which ones should be removed
        untraced = np.ones(self.ntrace, dtype=bool)
        rmtrace = np.zeros(self.ntrace, dtype=bool)

        # To hold the refined traces and mask
        cen = np.zeros_like(self.spat_cen)
        err = np.zeros_like(self.spat_err)
        msk = np.zeros_like(self.spat_msk)

        # Refine left then right
        for side in ['left', 'right']:
            
            # Get the image relevant to tracing this side
            _sobel_sig = self._side_dependent_sobel(side)

            # Identify the traces on the correct side: Traces on the
            # left side are negative.
            this_side = self.traceid < 0 if side == 'left' else self.traceid > 0

            # Loop continues until all traces are refined
            # TODO: Not sure why this while loop is necessary...
            i = 0
            while np.any(this_side & untraced):
                msgs.info('Iteration {0} for {1} side'.format(i+1, side))

                # TODO: Deal with single untraced edge

                # Get the traces to refine
                indx = this_side & untraced
                msgs.info('Number to retrace: {0}'.format(np.sum(indx)))

                # Find the most common row index
                trace_mask = self.bitmask.flagged(self.spat_msk, flag=self.bitmask.bad_flags())
                start_row = most_common_trace_row(trace_mask)
                msgs.info('Starting row is: {0}'.format(start_row))

                # Trace starting from this row
                msgs.info('Sequentially tracing first moment of Sobel-filtered image to higher '
                          'and lower spectral positions.')
                cen[:,indx], err[:,indx], msk[:,indx] \
                        = follow_trace_moment(_sobel_sig, start_row, self.spat_img[start_row,indx],
                                              ivar=ivar, mask=_mask, fwgt=fwgt, width=width,
                                              maxshift_start=maxshift_start,
                                              maxshift_follow=maxshift_follow, maxerror=maxerror,
                                              continuous=continuous, bitmask=self.bitmask)

                # Check the traces
                mincol = None if side == 'left' else 0
                maxcol = _sobel_sig.shape[1]-1 if side == 'left' else None
                good, bad = self.check_traces(cen, err, msk, subset=indx, mincol=mincol,
                                              maxcol=maxcol)

                # Save the results and update the book-keeping
                self.spat_cen[:,good] = cen[:,good]
                self.spat_err[:,good] = err[:,good]
                self.spat_msk[:,good | bad] |= msk[:,good | bad]
                untraced[good | bad] = False
                rmtrace[bad] = True

                # Increment the iteration counter
                i += 1

        # Update the image coordinates
        self.spat_img = np.round(self.spat_cen).astype(int)

        # Erase any previous fitting, PCA, and slit-mask-match results
        self.spat_fit_type = None
        self.spat_fit = None
        self.pca_type = None
        self.pca = None
        self.design = None
        self.objects = None

        # Remove bad traces and re-order the trace IDs
        if self.par['clip']:
            self.remove_traces(rmtrace)

        # Add to the log
        self.log += [inspect.stack()[0][3]]

    def check_traces(self, cen, err, msk, subset=None, mincol=None, maxcol=None):
        r"""
        Validate new trace data to be added.

        Steps are:
            - Remove duplicates based on the provided matching
              tolerance.
            - Remove traces that do not cover at least some fraction of
              the detector.
            - Remove traces that are at a minimum or maximum column
              (typically the edge of the detector).

        Used parameters from :attr:`par` (:class:`EdgeTracePar`) are
        `match_tol`, and `minimum_spec_length`.

        .. warning::
            `msk` is edited in-place

        Args:
            cen (`numpy.ndarray`_):
                The adjusted center of the refined traces. Shape is
                :math:`(N_{\rm spec}, N_{\rm refine},)`.
            err (`numpy.ndarray`_):
                The errors in the adjusted center of the refined
                traces. Shape is :math:`(N_{\rm spec}, N_{\rm
                refine},)`.
            msk (`numpy.ndarray`_):
                The mask bits for the adjusted center of the refined
                traces. Shape is :math:`(N_{\rm spec}, N_{\rm
                refine},)`.  This is edited in-place!
            subset (`numpy.ndarray`_, optional):
                Boolean array selecting the traces to compare. Shape
                is :math:`(N_{\rm trace},)`, with :math:`N_{\rm
                refine}` True values. It is expected that all the
                traces selected by a subset must be from the same
                slit side (left or right). If None, no repeat traces
                can be identified.
            mincol (:obj:`int`, optional):
                Clip traces that hit this minimum column value at the
                center row (`self.nspec//2`). If None, no traces
                clipped.
            maxcol (:obj:`int`, optional):
                Clip traces that hit this maximum column value at the
                center row (`self.nspec//2`). If None, no traces
                clipped.

        Returns:
            Returns two boolean arrays selecting the good and bad
            traces.  Shapes are :math:`(N_{\rm trace},)`.
        """
        # The nearest image column
        col = np.round(cen).astype(int)

        indx = np.ones(self.ntrace, dtype=bool) if subset is None else subset

        # Find repeat traces; comparison of traces must include
        # unmasked trace data, be traces of the same edge (left or
        # right), and be within the provided matching tolerance
        repeat = np.zeros_like(indx, dtype=bool)
        if subset is not None:
            msgs.info('Tolerance for finding repeat traces: {0:.1f}'.format(self.par['match_tol']))
            s = -1 if np.all(self.traceid[indx] < 0) else 1
            compare = (s*self.traceid > 0) & np.invert(indx)
            if np.any(compare):
                # Use masked arrays to ease exclusion of masked data
                _col = np.ma.MaskedArray(np.round(cen).astype(int), mask=msk > 0)
                trace_mask = self.bitmask.flagged(self.spat_msk, flag=self.bitmask.bad_flags())
                spat_img = np.ma.MaskedArray(self.spat_img, mask=trace_mask)
                mindiff = np.ma.amin(np.absolute(_col[:,indx,None]-spat_img[:,None,compare]),
                                    axis=(0,2))
                # TODO: This tolerance uses the integer image
                # coordinates, not the floating-point centroid
                # coordinates....
                repeat[indx] =  (mindiff.data < self.par['match_tol']) & np.invert(mindiff.mask)
                if np.any(repeat):
                    msk[:,repeat] = self.bitmask.turn_on(msk[:,repeat], 'DUPLICATE')
                    msgs.info('Found {0} repeat trace(s).'.format(np.sum(repeat)))

        # Find short traces
        short = np.zeros_like(indx, dtype=bool)
        if self.par['minimum_spec_length'] is not None:
            msgs.info('Minimum trace length (detector fraction): {0:.2f}'.format(
                        self.par['minimum_spec_length']))
            short[indx] = (np.sum(msk[:,indx] == 0, axis=0)
                                < self.par['minimum_spec_length']*self.nspec)
            if np.any(short):
                msk[:,short] = self.bitmask.turn_on(msk[:,short], 'SHORTRANGE')
                msgs.info('Found {0} short trace(s).'.format(np.sum(short)))

        # Find traces that are at the minimum column at the center row
        # TODO: Why only the center row?
        hit_min = np.zeros_like(indx, dtype=bool)
        if mincol is not None:
            hit_min[indx] = (col[self.nspec//2,indx] <= mincol) & (msk[self.nspec//2,indx] == 0)
            if np.any(hit_min):
                msk[:,hit_min] = self.bitmask.turn_on(msk[:,hit_min], 'HITMIN')
                msgs.info('{0} trace(s) hit the minimum centroid value.'.format(np.sum(hitmin)))
            
        # Find traces that are at the maximum column at the center row
        # TODO: Why only the center row?
        hit_max = np.zeros_like(indx, dtype=bool)
        if maxcol is not None:
            hit_max[indx] = (col[self.nspec//2,indx] >= maxcol) & (msk[self.nspec//2,indx] == 0)
            if np.any(hit_max):
                msk[:,hit_max] = self.bitmask.turn_on(msk[:,hit_max], 'HITMAX')
                msgs.info('{0} trace(s) hit the maximum centroid value.'.format(np.sum(hitmax)))

        # Good traces
        bad = indx & (repeat | short | hit_min | hit_max)
        msgs.info('Identified {0} bad trace(s) in all.'.format(np.sum(bad)))
        good = indx & np.invert(bad)
        return good, bad

    def is_synced(self):
        """
        Check if the slit edges are synced.
        """
        side = np.clip(self.traceid, -1, 1)
        return np.all(side[1:] + side[:-1] == 0)

    def check_synced(self, rebuild_pca=False):
        """
        Quality check and masking of the synchronized edges.

        Before executing this method, the slit edges must be
        synchronized (see :func:`sync`) and ordered spatially
        in left-right pairs (see :func:`spatial_sort`). The former is
        checked explicitly.

        Used parameters from :attr:`par` (:class:`EdgeTracePar`) are
        `minimum_slit_length` and `length_range`.

        Checks are:
            - Any trace falling off the edge of the detector is
              masked (see :class:`TraceBitMask`). This is the only
              check performed by default (i.e., no keywords are
              provided).
            - Traces that form a slit with a length (the difference
              between the left and right edges) below an absolute
              tolerance (i.e., `right-left < atol`) are masked. The
              absolute tolerance is set using the platescale provided
              by the spectrograph class, the spatial binning (from
              :attr:`binning`), and the minimum slit length in arcsec
              (`minimum_slit_length` in :attr:`par`).
            - Traces that form a slit with a length that is abnormal
              relative to the median width are masked; i.e.,
              `abs(log((right[0]-left[0])/median(right-left))) >
              log(1+rtol)`, where `rtol` is identically
              `length_range` in :attr:`par`.

        Args:
            rebuild_pca (:obj:`bool`, optional):
                If the pca exists, rebuild the PCA using the new
                traces and trace masks and the previous parameter
                set.

        """
        # Parse parameters and report
        atol = None
        rtol = self.par['length_range']
        if self.par['minimum_slit_length'] is not None:
            platescale = parse.parse_binning(self.binning)[1] \
                            * self.spectrograph.detector[self.det-1]['platescale']
            atol = self.par['minimum_slit_length']/platescale
            msgs.info('Binning: {0}'.format(self.binning))
            msgs.info('Platescale per binned pixel: {0}'.format(platescale))
            msgs.info('Minimum slit length (binned pixels): {0}'.format(atol))

        # Use the fit data if available
        trace = self.spat_cen if self.spat_fit is None else self.spat_fit
        # Keep track of whether or not any new masks were applied
        new_masks = False

        # Flag trace locations falling off the detector
        # TODO: This is a general check that is independent of whether
        # or not the traces are synced. Make check_traces more general
        # so that it can effectively be called at any time.
        indx = (trace < 0) | (trace >= self.nspat)
        if np.any(indx):
            new_masks = True
            self.spat_msk[indx] = self.bitmask.turn_on(self.spat_msk[indx], 'OFFDETECTOR')

        # Check the slits are synced
        if not self.is_synced():
            raise ValueError('Edge traces are not yet (or improperly) synced; run sync().')

        # Calculate the slit length
        slit_length = np.mean(np.squeeze(np.diff(trace.reshape(self.nspec,-1,2), axis=-1)), axis=0)

        if atol is not None:
            # Find any short slits (flag both edges of the slit)
            indx = np.repeat(slit_length < atol, 2)
            if np.sum(indx) == self.ntrace:
                msgs.error('All slits are too short!')
            if np.any(indx):
                new_masks = True
                msgs.info('Rejecting {0} slits that are too short.'.format(np.sum(indx)))
                self.spat_msk[:,indx] = self.bitmask.turn_on(self.spat_msk[:,indx], 'SHORTSLIT')

        if rtol is not None:
            msgs.info('Relative range in slit length limited to +/-{0:.1f}%'.format(rtol*100))
            # Find slits that are not within the provided fraction of
            # the median length
            indx = np.repeat(np.absolute(np.log(slit_length/np.median(slit_length)))
                             > np.log(1+rtol), 2)
            if np.any(indx):
                new_masks = True
                msgs.info('Rejecting {0} abnormally long or short slits.'.format(np.sum(indx)))
                self.spat_msk[:,indx] = self.bitmask.turn_on(self.spat_msk[:,indx], 'ABNORMALSLIT')

        # TODO: Check that slit edges meet minimum slit gap?

        if self.par['clip']:
            # Remove traces that have been fully flagged as bad
            rmtrace = np.all(self.bitmask.flagged(self.spat_msk, flag=self.bitmask.bad_flags()),
                             axis=0)
            if np.sum(rmtrace) == self.ntrace:
                msgs.error('All slit edges are fully masked!')
            self.remove_traces(rmtrace, rebuild_pca=rebuild_pca)
        elif new_masks:
            # Reset the PCA if new masks are applied
            self._reset_pca(rebuild_pca and self.pca is not None)

    def remove_traces(self, indx, resort=True, rebuild_pca=False):
        r"""
        Remove a set of traces.

        Args:
            indx (array-like):
                The boolean array with the traces to remove. Length
                must be :math:`(N_{\rm trace},)`.
            resort (:obj:`bool`, optional):
                Re-sort the traces and trace IDs to be sequential in
                the spatial direction. See :func:`spatial_sort`.
            rebuild_pca (:obj:`bool`, optional):
                If the pca exists, rebuild it using the new traces
                and the previous parameter set.
        """
        # Make sure there are traces to remove
        if not np.any(indx):
            msgs.warn('No traces removed.')
            return
        msgs.info('Removing {0} edge traces.'.format(np.sum(indx)))

        # Reset the trace data
        keep = np.invert(indx)
        self.spat_img = self.spat_img[:,keep]
        self.spat_cen = self.spat_cen[:,keep]
        self.spat_err = self.spat_err[:,keep]
        self.spat_msk = self.spat_msk[:,keep]
        if self.spat_fit is not None:
            self.spat_fit = self.spat_fit[:,keep]
        self.traceid = self.traceid[keep]

        if resort:
            # Resort by the spatial dimension
            self.spatial_sort()

        # Reset the PCA
        self._reset_pca(rebuild_pca and self.pca is not None)

    def spatial_sort(self):
        """
        Sort the traces spatially.

        The sorting is based on the mean measured center
        (:attr:`spat_cen`) of each trace. The trace IDs are also
        reassigned to be sorted spatially; i.e., the trace IDs for
        three synced slits would be `[-1, 1, -2, 2, -3, 3]`.
        
        All attributes are edited in-place.
        """
        msgs.info('Re-sorting edge traces by their mean spatial position.')

        # Sort the traces by their spatial position (always use
        # measured positions even if fit positions are available)
        # TODO: Instead do this based on a reference row?
        srt = np.argsort(np.mean(self.spat_cen, axis=0))

        self.traceid = self.traceid[srt]
        self.spat_img = self.spat_img[:,srt]
        self.spat_cen = self.spat_cen[:,srt]
        self.spat_err = self.spat_err[:,srt]
        self.spat_msk = self.spat_msk[:,srt]
        if self.spat_fit is not None:
            self.spat_fit = self.spat_fit[:,srt]

        # Reorder the trace numbers
        indx = self.traceid < 0
        self.traceid[indx] = -1-np.arange(np.sum(indx))
        indx = np.invert(indx)
        self.traceid[indx] = 1+np.arange(np.sum(indx))

    def _reset_pca(self, rebuild):
        """"
        Reset the PCA decomposition.

        The PCA is reset by either rebuilding it with the previous
        set of parameters (`rebuild` is True) or removing it (setting
        the relevant attributes to `None`; `rebuild` is False).
        """
        if rebuild:
            # Rebuild the PCA using the previous parameters
            return self.build_pca(use_center=self.pca_type == 'center') 
        # Remove the existing PCA
        self.pca_type = None
        self.pca = None

    def current_trace_img(self):
        """
        Return an image with the trace IDs at the locations of each
        edge in the original image.
        """
        edge_img = np.zeros((self.nspec, self.nspat), dtype=int)
        i = np.tile(np.arange(self.nspec), (self.ntrace,1)).T.ravel()
        edge_img[i, self.spat_img.ravel()] = np.tile(self.traceid, (self.nspec,1)).ravel()
        return edge_img

    def fit_refine(self, weighting='uniform', debug=False, idx=None):
        """
        Fit the edge location data with a function.

        Primarily a wrapper for :func:`fit_trace`, run once per edge
        side (left and right). See documentation of
        :func:`fit_trace`.

        Used parameters from :attr:`par` (:class:`EdgeTracePar`) are
        `fit_function`, `fit_order`, `fwhm_uniform`, `fwhm_gaussian`,
        `fit_maxdev`, `fit_maxiter`, and `fit_niter`.

        Args:
            weighting (:obj:`str`, optional):
                The weighting to apply to the position within each
                integration window (see :func:`fit_trace` and
                :func:`pypeit.moment.moment1d`).
            debug (:obj:`bool`, optional):
                Run in debug mode.
            idx (`numpy.ndarray`_, optional):
                Array of strings with the IDs for each object. Used
                only if show_fits is true for the plotting. Default
                is just a running number.
        """
        # Parse parameters and report
        function = self.par['fit_function']
        order = self.par['fit_order']
        fwhm = self.par['fwhm_uniform'] if weighting == 'uniform' else self.par['fwhm_gaussian']
        maxdev = self.par['fit_maxdev']
        maxiter = self.par['fit_maxiter']
        niter = self.par['fit_niter']
        xmin = 0.
        xmax = self.nspec-1.

        msgs.info('Trace fitting function: {0}'.format(function))
        msgs.info('Trace fitting order: {0}'.format(order))
        msgs.info('Weighting for remeasuring edge centroids: {0}'.format(weighting))
        msgs.info('FWHM parameter for remeasuring edge centroids: {0:.1f}'.format(fwhm))
        msgs.info('Maximum deviation for fitted data: {0:.1f}'.format(maxdev))
        msgs.info('Maximum number of rejection iterations: {0}'.format(maxiter))
        msgs.info('Number of remeasuring and refitting iterations: {0}'.format(niter))

        # Generate bogus ivar and mask once here so that they don't
        # have to be generated multiple times.
        # TODO: Keep these as work space as class attributes?
        ivar = np.ones_like(self.sobel_sig, dtype=float)
        mask = np.zeros_like(self.sobel_sig, dtype=bool) \
                    if self.trace_msk is None else self.trace_msk

        # Parameters
        _order = self.par['trace_npoly'] if order is None else order
        _function = self.par['function'] if function is None else function

        trace_fit = np.zeros_like(self.spat_cen, dtype=float)
        trace_cen = np.zeros_like(self.spat_cen, dtype=float)
        trace_err = np.zeros_like(self.spat_cen, dtype=float)
        bad_trace = np.zeros_like(self.spat_cen, dtype=bool)

        # Flag bad traces; excludes inserted traces
        trace_mask = self.bitmask.flagged(self.spat_msk, flag=self.bitmask.bad_flags())

        # Fit both sides
        for side in ['left', 'right']:
            # Get the image relevant to tracing this side
            _sobel_sig = self._side_dependent_sobel(side)
            # Select traces on this side
            this_side = self.traceid < 0 if side == 'left' else self.traceid > 0
            # Perform the fit
            trace_fit[:,this_side], trace_cen[:,this_side], trace_err[:,this_side], \
                bad_trace[:,this_side], _ \
                    = fit_trace(_sobel_sig, self.spat_cen[:,this_side], _order, ivar=ivar,
                                mask=mask, trace_mask=trace_mask[:,this_side],
                                weighting=weighting, fwhm=fwhm, function=_function,
                                maxdev=maxdev, maxiter=maxiter, niter=niter, show_fits=debug,
                                idx=idx, xmin=xmin, xmax=xmax)

        # Save the results of the edge measurements ...
        self.spat_cen = trace_cen
        self.spat_err = trace_err
        # ... and the model fits
        self.spat_fit = trace_fit
        self.spat_fit_type = '{0} : order={1}'.format(_function, _order)
        self.spat_img = np.round(self.spat_fit).astype(int)
        # TODO: Flag pixels with a bad_trace (bad moment measurement)?
        # moment1d replaces any bad measurement with the input center
        # value, so may not want to flag them. Flagging would likely
        # mean that the traces would be ignored in subsequent
        # analysis...
        self.log += [inspect.stack()[0][3]]

    def build_pca(self, use_center=False, debug=False):
        """
        Build a PCA model of the current trace data using both the
        left and right traces.

        Primarily a wrapper for instantiation of :attr:`pca`, which
        has type :class:`EdgeTracePCA`. After executing this, traces
        can be predicted using the pca by calling
        `self.pca.predict(spat)`; see :func:`EdgeTracePCA.predict`.

        If no parametrized function has been fit to the trace data or
        if specifically requested (see `use_center`), the PCA is
        based on the measured trace centroids; othwerwise, the PCA
        uses the parametrized trace fits.

        Used parameters from :attr:`par` (:class:`EdgeTracePar`) are
        `pca_n`, `pca_var_percent`, `pca_function`, `pca_order`,
        `pca_sigrej`, `pca_maxrej`, and `pca_maxiter`.

        Args:
            use_center (:obj:`bool`, optional):
                Use the center measurements for the PCA decomposition
                instead of the functional fit to those data. This is
                only relevant if both are available. If not fits have
                been performed, the function will automatically use
                the center measurements.
            debug (:obj:`bool`, optional):
                Run in debug mode.
        """
        # Parse parameters and report
        npca = self.par['pca_n']
        pca_explained_var = self.par['pca_var_percent']
        function = self.par['pca_function']
        order = self.par['pca_order']
        lower, upper = self.par['pca_sigrej'] if hasattr(self.par['pca_sigrej'], '__len__') \
                        else (self.par['pca_sigrej'],)*2
        maxrej = self.par['pca_maxrej']
        maxiter = self.par['pca_maxiter']

        if npca is not None:
            msgs.info('Restricted number of PCA components: {0}'.format(npca))
        if pca_explained_var is not None:
            msgs.info('Requested pecentage of variance explained by PCA: {0:.1f}'.format(
                        pca_explained_var))
        msgs.info('Function fit to PCA coefficients: {0}'.format(function))
        msgs.info('Lower sigma rejection: {0:.1f}'.format(lower))
        msgs.info('Upper sigma rejection: {0:.1f}'.format(upper))
        msgs.info('Maximum number of rejections per iteration: {0}'.format(maxrej))
        msgs.info('Maximum number of rejection iterations: {0}'.format(maxiter))

        # Check the state of the current object
        if self.pca is not None:
            msgs.warn('PCA model already exists and will be overwritten.')
        if self.spat_fit is None and not use_center:
            msgs.warn('No trace fits exits.  PCA based on trace centroid measurements.')
        self.pca_type = 'center' if self.spat_fit is None or use_center else 'fit'

        # When constructing the PCA, ignore bad trace component *and*
        # any traces inserted by hand.
        trace_mask = self.bitmask.flagged(self.spat_msk)

        # Only use traces that are unmasked for at least some fraction of the detector
        # TODO: Make this minimum length check a paramter
        minimum_spec_length=0.9
        # TODO: Is there a way to propagate the mask to the PCA?
        # TODO: Keep a separate mask specifically for the fit data?
        use_trace = np.sum(np.invert(trace_mask), axis=0)/self.nspec > minimum_spec_length
        if np.sum(use_trace) < 2:
            msgs.error('Insufficient traces for PCA decomposition.')
        trace_inp = self.spat_cen[:,use_trace] if self.spat_fit is None or use_center \
                        else self.spat_fit[:,use_trace]
        msgs.info('Using {0}/{1} traces in the PCA analysis (omitting short traces).'.format(
                  np.sum(use_trace), self.ntrace))

        # Instantiate the PCA
        self.pca = EdgeTracePCA(trace_inp, npca=npca, pca_explained_var=pca_explained_var,
                                reference_row=most_common_trace_row(trace_mask[:,use_trace]))

        # Set the order of the function fit to the PCA coefficiencts:
        # Order is set to cascade down to lower order for components
        # that account for a smaller percentage of the variance.
        _order = np.clip(order - np.arange(self.pca.npca), 1, None).astype(int)
        msgs.info('Order of function fit to each component: {0}'.format(_order))
        # Run the fit
        self.pca.build_interpolator(_order, function=function, lower=lower, upper=upper, minx=0.,
                                    maxx=self.nspat-1., maxrej=maxrej, maxiter=maxiter,
                                    debug=debug)

    def pca_refine(self, use_center=False, debug=False, force=False):
        """
        Use a PCA decomposition to refine the traces.

        If no parametrized function has been fit to the trace data or
        if specifically requested (see `use_center`), the PCA is
        based on the measured trace centroids; othwerwise, the PCA
        uses the parametrized trace fits.

        If needed or forced to, this first executes :func:`build_pca`
        and then uses :func:`EdgeTracePCA.predict` to use the PCA to
        reset the trace data.

        Args:
            use_center (:obj:`bool`, optional):
                Use the center measurements for the PCA decomposition
                instead of the functional fit to those data. This is
                only relevant if both are available. If not fits have
                been performed, the function will automatically use
                the center measurements.
            debug (:obj:`bool`, optional):
                Run in debug mode.
            force (:obj:`bool`, optional):
                Force the recalculation of the PCA even if it has
                already been done.
        """
        # NOTE: All parameters parsed by build_pca
        # Perform the PCA decomposition if necessary
        _pca_type = 'center' if use_center or self.spat_fit is None else 'fit'
        if force or self.pca is None or self.pca_type != _pca_type:
            self.build_pca(use_center=use_center, debug=debug)

        # Get the spatial positions of each trace at the reference row
        trace_ref = self.spat_cen[self.pca.reference_row,:] if self.pca_type == 'center' \
                        else self.spat_fit[self.pca.reference_row,:]

        # Predict the traces
        self.spat_fit = self.pca.predict(trace_ref)
        self.spat_fit_type = 'pca'
        self.log += [inspect.stack()[0][3]]

    def peak_refine(self, rebuild_pca=False, debug=False):
        """
        Refine the trace by isolating peaks and troughs in the
        Sobel-filtered image.

        This function *requires* that the PCA model exists; see
        :func:`build_pca` or :func:`pca_refine`. It is also primarily
        a wrapper for :func:`peak_trace`.

        Used parameters from :attr:`par` (:class:`EdgeTracePar`) are
        `edge_thresh`, `smash_range`, `edge_detect_clip`, trace_median_frac`,
        `trace_thresh`, `fit_function`, `fit_order`, `fwhm_uniform`,
        `fwhm_uniform`, `niter_gaussian`, `niter_gaussian`,
        `fit_maxdev`, and `fit_maxiter`.

        Args:
            rebuild_pca (:obj:`bool`, optional):
                This method fundamentally resets the trace data,
                meaning that the PCA is no longer valid. Use this
                boolean to have the method rebuild the PCA based on
                the refined traces. Note that the PCA is *not* then
                used to reset the fitted trace data; i.e.,
                :attr:`spat_fit` remains based on the output of
                :func:`peak_trace`.
            debug (:obj:`bool`, optional):
                Run in debug mode.
        """
        if self.pca is None:
            raise ValueError('Must first run the PCA analysis fo the traces; run build_pca.')

        # Parse parameters and report
        peak_thresh = self.par['edge_thresh']
        smash_range = self.par['smash_range']
        peak_clip = self.par['edge_detect_clip']
        trace_median_frac = self.par['trace_median_frac']
        trace_thresh = self.par['trace_thresh']
        function = self.par['fit_function']
        order = self.par['fit_order']
        fwhm_uniform = self.par['fwhm_uniform']
        niter_uniform = self.par['niter_uniform']
        fwhm_gaussian = self.par['fwhm_gaussian']
        niter_gaussian = self.par['niter_gaussian']
        maxdev = self.par['fit_maxdev']
        maxiter = self.par['fit_maxiter']

        msgs.info('Threshold for peak detection: {0:.1f}'.format(peak_thresh))
        msgs.info('Detector range (spectral axis) collapsed: {0}'.format(smash_range))
        msgs.info('Image fraction filtered for trace masking: {0:.2f}'.format(trace_median_frac))
        msgs.info('Threshold for trace masking: {0:.1f}'.format(trace_thresh))
        msgs.info('Trace fitting function: {0}'.format(function))
        msgs.info('Trace fitting order: {0}'.format(order))
        msgs.info('FWHM parameter for uniform-weighted centroids: {0:.1f}'.format(fwhm_uniform))
        msgs.info('Number of uniform-weighted iterations: {0:.1f}'.format(niter_uniform))
        msgs.info('FWHM parameter for Gaussian-weighted centroids: {0:.1f}'.format(fwhm_gaussian))
        msgs.info('Number of Gaussian-weighted iterations: {0:.1f}'.format(niter_gaussian))
        msgs.info('Maximum deviation for fitted data: {0:.1f}'.format(maxdev))
        msgs.info('Maximum number of rejection iterations: {0}'.format(maxiter))

        # TODO: Much of this is identical to fit_refine; abstract to a
        # single function that selects the type of refinement to make?
        # Also check traces after fitting or PCA?

        # Generate bogus ivar and mask once here so that they don't
        # have to be generated multiple times.
        # TODO: Keep these as work space as class attributes? so that
        # they don't need to be reinstantiated.
        ivar = np.ones_like(self.sobel_sig, dtype=float)
        mask = np.zeros_like(self.sobel_sig, dtype=bool) \
                    if self.trace_msk is None else self.trace_msk

        # Get the image relevant to tracing
        _sobel_sig = prepare_sobel_for_trace(self.sobel_sig, boxcar=5, side=None)

        # Find and trace both peaks and troughs in the image. The input
        # trace data (`trace` argument) is the PCA prediction of the
        # trace that passes through each spatial position at the
        # reference spectral pixel.
        trace_fit, trace_cen, trace_err, bad_trace, nleft \
                = peak_trace(_sobel_sig, ivar=ivar, mask=mask,
                             trace=self.pca.predict(np.arange(self.nspat)),
                             smash_range=smash_range, peak_thresh=peak_thresh, peak_clip=peak_clip,
                             trough=True, trace_median_frac=trace_median_frac,
                             trace_thresh=trace_thresh, fwhm_uniform=fwhm_uniform,
                             fwhm_gaussian=fwhm_gaussian, function=function, order=order,
                             maxdev=maxdev, maxiter=maxiter, niter_uniform=niter_uniform,
                             niter_gaussian=niter_gaussian, debug=debug)

        # Assess the output
        ntrace = trace_fit.shape[1]
        if ntrace < self.ntrace:
            msgs.warn('Found fewer traces using peak finding than originally available.  '
                      'May want to reset peak threshold.')

        # Reset the trace data
        self.spat_msk = np.zeros_like(bad_trace, dtype=self.bitmask.minimum_dtype())
        if np.any(bad_trace):
            self.spat_msk[bad_trace] = self.bitmask.turn_on(self.spat_msk[bad_trace], 'MATHERROR')
        self.traceid = np.zeros(ntrace, dtype=int)
        self.traceid[:nleft] = -1-np.arange(nleft)
        self.traceid[nleft:] = 1+np.arange(ntrace-nleft)
        self.spat_fit = trace_fit
        self.spat_fit_type = '{0} : order={1}'.format(function, order)
        self.spat_cen = trace_cen
        self.spat_err = trace_err
        self.spat_img = np.round(self.spat_fit).astype(int)

        # Spatially sort the traces
        self.spatial_sort()
        # Reset the PCA
        self._reset_pca(rebuild_pca)
        self.log += [inspect.stack()[0][3]]

    # TODO: Make this a core function?
    def _get_insert_locations(self):
        """
        Find where edges need to be inserted.

        See :func:`sync`.

        Returns:
            Three `numpy.ndarray`_ objects are returned:
                - An integer vector identifying the type of side for
                  the fully synchronized edge set. Elements of the
                  vector should alternate left (-1) and right (1).
                - A boolean array selecting the edges in the returned
                  list of sides that should be added to the existing
                  trace set.
                - An array with the indices in the existing trace
                  arrays where the new traces should be inserted.
        """
        side = np.clip(self.traceid, -1, 1)
        add_edge = np.zeros(self.ntrace, dtype=bool)
        add_indx = np.zeros(self.ntrace, dtype=int)

        if side[0] > 0:
            # First side is a right, so add a left
            side = np.insert(side, 0, -1)
            add_edge = np.insert(add_edge, 0, True)
            add_indx = np.insert(add_indx, 0, 0)
        if side[-1] < 0:
            # Last side is a left, so add a right
            side = np.append(side, 1)
            add_edge = np.append(add_edge, True)
            add_indx = np.append(add_indx, self.ntrace)

        # Find missing lefts and rights
        diff = side[1:] + side[:-1]
        if np.all(diff == 0):
            # All edges are paired left-right
            return side, add_edge, add_indx

        # Missing lefts have diff == 2, rights have diff == -2
        missing = np.where(diff != 0)[0] + 1
        # Set the full side vector
        side = np.insert(side, missing, -np.sign(diff[missing-1]))
        # Keep track of which edges will have been added
        add_edge = np.insert(add_edge, missing, np.ones(missing.size, dtype=bool))
        # Keep track of where the new edges should be inserted in the
        # existing set
        add_indx = np.insert(add_indx, missing, missing)
        # Return the edges to add, their side, and where to insert them 
        return side, add_edge, add_indx

    def _get_reference_locations(self, trace, add_edge):
        """
        Insert the reference locations for the traces to add.

        Used parameters from :attr:`par` (:class:`EdgeTracePar`) are
        `sync_center`, `sync_to_edge`, and `min_slit_gap`. See
        :func:`sync`.

        Args:
            trace (`numpy.ndarray`_):
                Trace data to use for determining new edge locations.
            add_edge (`numpy.ndarray`_):
                Boolean array indicating that a trace in the new
                array is an added trace. The number of *False*
                entries in `add_edge` should match the length of the
                2nd axis of `trace`.
        """
        # Parse parameters and report
        center_mode = self.par['sync_center']
        to_edge = self.par['sync_to_edge']
        min_slit_gap = self.par['min_slit_gap']

        msgs.info('Mode used to set spatial position of new traces: {0}'.format(center_mode))
        msgs.info('For first left and last right, set trace to the edge: {0}'.format(to_edge))

        if center_mode not in ['median', 'nearest']:
            raise ValueError('Unknown centering mode: {0}'.format(center_mode))

        # Get the reference row for the placement calculation; allow
        # the use of inserted traces.
        trace_mask = self.bitmask.flagged(self.spat_msk, flag=self.bitmask.bad_flags())
        reference_row = most_common_trace_row(trace_mask) if self.pca is None \
                            else self.pca.reference_row

        # Check that the trace data are at this row
        if not np.array_equal(np.arange(trace.shape[1]), np.argsort(trace[reference_row,:])):
            raise ValueError('Trace data must be spatially sorted.')

        # Build a masked array with the trace positions at that row,
        # masked where new traces are supposed to go.
        trace_ref = np.ma.masked_all(add_edge.size)
        trace_ref[np.invert(add_edge)] = trace[reference_row,:]
        trace_ref = trace_ref.reshape(-1,2)

        # Get the length and center of each slit in pixels
        nslits = trace_ref.shape[0]
        slit_length = np.ma.diff(trace_ref, axis=1).ravel()
        slit_center = np.ma.mean(trace_ref, axis=1)

        # Mask any bad calculations
        missing_a_side = np.ma.any(add_edge.reshape(-1,2), axis=1)
        slit_length[missing_a_side] = np.ma.masked
        slit_center[missing_a_side] = np.ma.masked

        # Determine how to offset and get the reference locations of the new edges to add
        if center_mode == 'median':
            offset = np.full(nslits, np.ma.median(slit_length), dtype=float)
        if center_mode == 'nearest':
            # Find the index of the nearest slit with both existing
            # edges (i.e. has an unmasked slit length)
            nearest = utils.nearest_unmasked(slit_center, use_indices=True)
            # The offset is the slit length of the nearest valid slit
            offset = slit_length.data[nearest]

        # Set the new edge trace reference locations
        for slit in range(nslits):
            if not slit_length.mask[slit]:
                # Both slit edges already defined
                continue
            if trace_ref.mask[slit,0]:
                # Add the left edge
                trace_ref[slit,0] = 0 if slit == 0 and to_edge \
                                        else trace_ref[slit,1] - offset[slit]
                continue
            # Need to add the right edge
            trace_ref[slit,1] = self.nspat - 1 if slit == nslits-1 and to_edge \
                                    else trace_ref[slit,0] + offset[slit]

        # Check for slit overlaps
        overlap = trace_ref[1:,0] - trace_ref[:-1,1] < min_slit_gap
        if np.any(overlap):
            msgs.warn('Found {0} overlapping slit(s) using edge offsets.  '.format(np.sum(overlap))
                      + 'Moving left edges to minimum gap of {0} pixel(s).'.format(min_slit_gap))
            indx = np.where(overlap)[0]+1
            trace_ref[indx,0] = trace_ref[indx-1,1] + min_slit_gap

        # TODO: Nothing should now be masked. Get rid of this once
        # satisfied that the coding is correct.
        if np.any(trace_ref.mask):
            raise ValueError('Coding error: this should not happen')
        return trace_ref.data.ravel()

    def _nudge_traces(self, trace):
        r"""
        Nudge traces away from the detector edge.

        Traces are shifted spatially, up to a maximum value set by
        `max_nudge`, to be no closer than a minimum of `det_buffer`
        pixels from the detector edges. Both parameters are pulled
        from :attr:`par` (:class:`EdgeTracePar`).

        Args:
            trace (`numpy.ndarray`_):
                Array with trace locations to adjust. Must be 2D with
                shape :math:`(N_{\rm spec}, N_{\rm trace})`.

        Returns:
            `numpy.ndarray`_: The nudged traces.
        """
        # Check input
        if self.par['max_nudge'] <= 0:
            # Nothing to do
            return trace
        if trace.shape[0] != self.nspec:
            raise ValueError('Traces have incorrect length.')
        _buffer = self.par['det_buffer']
        if _buffer < 0:
            msgs.warn('Buffer must be greater than 0; ignoring.')
            _buffer = 0

        msgs.info('Nudging traces, by at most {0} pixel(s)'.format(self.par['max_nudge'])
                  + ', to be no closer than {0} pixel(s) from the detector edge.'.format(_buffer))

        # Should never happen, but this makes a compromise if a trace
        # crosses both the left and right spatial edge of the
        # detector...
        offset = np.clip(_buffer - np.amin(trace, axis=0), 0, self.par['max_nudge']) \
                    + np.clip(self.nspat - 1 - _buffer - np.amax(trace, axis=0),
                              -self.par['max_nudge'], 0)
        # Offset and return the traces
        return trace + offset[None,:]

    def sync(self, rebuild_pca=True):
        """
        Match left and right edge traces to construct slit edge pairs.

        First, the method ensures that the edge traces are sorted
        spatiall; see :func:`spatial_sort`. Then it determines where
        traces need to be inserted to create a left-right pair (using
        :func:`_get_insert_locations`).

        The next steps are determine the reference positions where
        the traces should be inserted and the shape the trace should
        take. The former is determined by the `sync_center` parameter
        in :attr:`par`, and the latter is determined by `sync_trace`.

        Current the position of the trace is determine either by the
        median of or the nearest slit length for the already
        identified left-right edge pairs. The form of the trace with
        spectral pixel can either be predicted by the PCA
        decomposition or take exactly the same shape as the nearest
        left or right edge.

        After the new traces are generated, they are added to the
        edge traces using :func:`insert_traces` and flagged as having
        been inserted by the sync operation. The full set of
        synchronized edge traces are then checked using
        :func:`check_synced`.

        Used parameters from :attr:`par` (:class:`EdgeTracePar`) are
        `det_buffer` and `sync_trace`.

        Args:
            rebuild_pca (:obj:`bool`, optional):
                If the pca exists and traces are removed (see
                :func:`check_synced`), rebuild the PCA using the new
                traces and the previous parameter set. Note that
                inserted traces are *not* included in the PCA
                decomposition.
        """
        # Check input
        if self.par['sync_trace'] not in ['pca', 'nearest']:
            raise ValueError('Unknown trace mode: {0}'.format(self.par['sync_trace']))
        if self.par['sync_trace'] == 'pca' and self.pca is None:
            raise ValueError('The PCA decomposition does not exist.  Either run self.build_pca '
                             'or use a different trace_mode.')

        # Make sure that the traces are sorted spatially
        # TODO: This should be the convention of the class and should
        # *always* be true; instead check for this and raise an error
        # if it's not?
        self.spatial_sort()

        # Find the edges to add, what side they're on, and where to
        # insert them into the existing trace array
        side, add_edge, add_indx = self._get_insert_locations()
        if not np.any(add_edge):
            # No edges to add
            return

        # TODO: Report to the user what's going on
        msgs.info('Found {0} left and {1} right trace(s) to add.'.format(
                    np.sum((side == -1) & add_edge), np.sum((side == 1) & add_edge)))

        # Allow the edges to be synced, even if a fit hasn't been done yet
        trace = self.spat_cen if self.spat_fit is None else self.spat_fit

        # Instantiate the traces to add
        trace_add = np.zeros((self.nspec, np.sum(add_edge)), dtype=float)

        # If there was only one edge, just add the other one
        if side.size == 2:
            msgs.warn('Only one edge traced.  Ignoring center_mode and adding edge at the '
                      'opposite edge of the detector.')
            msgs.info('Detector edge buffer: {0}'.format(self.par['det_buffer']))
            # TODO: PCA would have failed because there needs to be at least
            # two traces.  Get rid of this test eventually.
            if self.par['sync_trace'] == 'pca':
                raise ValueError('Coding error: this should not happen.')
            # Set the offset to add to the existing trace
            offset = self.par['det_buffer'] - np.amin(trace[:,0]) if add_edge[0] \
                        else self.nspat - np.amax(trace[:,0]) - self.par['det_buffer']
            # Construct the trace to add and insert it
            trace_add[:,0] = trace[:,0] + offset
            self.insert_traces(side[add_edge], trace_add, loc=add_indx, mode='sync')
            return

        # Get the reference locations for the new edges
        trace_ref = self._get_reference_locations(trace, add_edge)

        # Predict the traces either using the PCA or using the nearest slit edge
        if self.par['sync_trace'] == 'pca':
            trace_add = self.pca.predict(trace_ref[add_edge])
        if self.par['sync_trace'] == 'nearest':
            # Index of trace nearest the ones to add
            # TODO: Force it to use the nearest edge of the same side;
            # i.e., when inserting a new right, force it to use the
            # nearest right instead of the nearest left?
            nearest = utils.nearest_unmasked(np.ma.MaskedArray(trace_ref, mask=add_edge))
            # Indices of the original traces
            indx = np.zeros(len(add_edge), dtype=int)
            indx[np.invert(add_edge)] = np.arange(self.ntrace)
            # Offset the original traces by a constant based on the
            # reference row to construct the new traces.
            trace_add = trace[:,indx[nearest[add_edge]]] + trace_ref[add_edge] \
                            - trace_ref[nearest[add_edge]]

        # Insert the new traces, check the full synchronized list and
        # log completion of the method
        self.insert_traces(side[add_edge], trace_add, loc=add_indx[add_edge], mode='sync')
        self.check_synced(rebuild_pca=rebuild_pca)
        self.log += [inspect.stack()[0][3]]

    def insert_traces(self, side, trace, loc=None, mode='user', resort=True):
        r"""
        Insert/append a set of edge traces.

        New traces to add are first nudged away from the detector
        edge (see :func:`_nudge_traces`) according to parameters
        `max_nudge` and `det_buffer` from :attr:`par`
        (:class:`EdgeTracePar`). They are then inserted or appended
        to the existing traces and masked according to the provided
        `mode`. The traces are added to *both* the measured centroid
        list and the fitted model data. Then the full list of traces
        can be resorted spatially, according to the provided
        `resort`.

        Typically, the inserted traces will be masked, which means
        that any existing PCA decomposition will be unchanged.
        However, if `mode` is None, these inserted traces would be
        used in the construction of the PCA.

        Args:
            side (:obj:`int`, `numpy.ndarray`_):
                Side for each trace to be added: -1 for left, 1 for
                right. Shape is :math:`(N_{\rm new},)`.
            trace (`numpy.ndarray`_):
                Array with one or more vectors of trace locations.
                Can be 1D or 2D with shape :math:`(N_{\rm spec},)` or
                :math:`(N_{\rm spec}, N_{\rm new})`, respectively.
            loc (:obj:`int`, `numpy.ndarray`_, optional):
                Indices in the current trace arrays at which to
                insert the new traces; see `numpy.insert`. If None,
                traces are appended.
            mode (:obj:`str`, optional):
                Mode used for generating the traces to insert used to
                flag the traces. Options are:
                    - None: Traces are simply inserted without
                    flagging.
                    - 'user': Traces are the result of a user
                    request.
                    - 'sync': Traces were generated by synchronizing
                    left and right traces.
                    - 'mask': Traces were generated based on the
                    expected slit positions from mask design data.
            resort (:obj:`bool`, optional):
                Resort the traces in the spatial dimension; see
                :func:`spatial_sort`.
        """
        # Check input
        _side = np.atleast_1d(side)
        ntrace = _side.size
        _trace = trace.reshape(-1,1) if trace.ndim == 1 else trace
        if _trace.shape[1] != ntrace:
            raise ValueError('Number of sides does not match the number of traces to insert.')
        if loc is None:
            # Insertion locations not provided so append
            loc = np.full(ntrace, self.ntrace, dtype=int)
        if loc.size != ntrace:
            raise ValueError('Number of sides does not match the number of insertion locations.')

        # Nudge the traces
        _trace = self._nudge_traces(_trace)

        # Set the mask
        mask = np.zeros(_trace.shape, dtype=self.bitmask.minimum_dtype())
        # Flag the traces pixels that fall off the detector
        indx = (_trace < 0) | (_trace >= self.nspat)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'OFFDETECTOR')
        # Flag the mode used to insert these traces, if provided
        if mode == 'user':
            mask = self.bitmask.turn_on(mask, 'USERINSERT')
        elif mode == 'sync':
            mask = self.bitmask.turn_on(mask, 'SYNCINSERT')
        elif mode == 'mask':
            mask = self.bitmask.turn_on(mask, 'MASKINSERT')

        # Set the ID numbers for the new traces
        _traceid = np.empty(ntrace, dtype=int)
        indx = _side < 0
        if np.any(indx):
            _traceid[indx] = np.amin(self.traceid) - 1 - np.arange(np.sum(indx))
        indx = _side > 0
        if np.any(indx):
            _traceid[indx] = np.amax(self.traceid) + 1 + np.arange(np.sum(indx))

        # Add the new traces. The new traces are added to both the
        # fitted list and the center list!
        self.traceid = np.insert(self.traceid, loc, _traceid)
        self.spat_img = np.insert(self.spat_img, loc, np.round(_trace).astype(int), axis=1)
        self.spat_cen = np.insert(self.spat_cen, loc, _trace, axis=1)
        self.spat_err = np.insert(self.spat_err, loc, np.zeros(_trace.shape, dtype=float), axis=1)
        self.spat_msk = np.insert(self.spat_msk, loc, mask, axis=1)
        self.spat_fit = np.insert(self.spat_fit, loc, _trace, axis=1)

        if resort:
            self.spatial_sort()

    def mask_refine(self, design_file=None):
        """
        Use the mask design data to refine the edge trace positions.

        Use of this method requires:
            - a PCA decomposition is available,
            - the traces are synchronized into left-right pairs, and
            - :attr:`spectrograph` has a viable `get_slitmask` method
              to read slit mask design data from a file. That file is
              either provided directly or pulled from one of the
              files used to construct the trace image; see
              `design_file`. The result of the `get_slitmask` method
              must provide a
              :class:`pypeit.spectrographs.slitmask.SlitMask` object
              with the slit-mask design data.

        Also useful, but not required, is for :attr:`spectrograph` to
        have a viable `get_detector_map` method that provides a
        :class:`pypeit.spectrograph.opticalmodel.DetectorMap` object,
        which is used to provide a guess offset between the slit-mask
        focal-plane positions and the trace pixel positions. If no
        such `get_detector_method` exists, the guess offset is::

            this

        and the match between expected and traced slit positions may
        be unstable.

        The method uses
        :class:`pypeit.spectrographs.slitmask.SlitRegister` to match
        the expected and traced position and identify both missing
        and erroneous trace locations. The former are used to add new
        traces and the latter are removed. The method also constructs
        the :attr:`design` and :attr:`objects` tables, depending on
        the data accessible via the
        :class:`pypeit.spectrographs.slitmask.SlitMask` instance.

        Args:
            design_file (:obj:`str`, optional):
                A file with the mask design data. If None, the method
                will use the first file in :attr:`files`; if
                :attr:`files` is also None, the method will raise an
                exception.
        """
        # The PCA decomposition must have already been determined
        if self.pca is None:
            msgs.error('Must first run the PCA analysis for the traces; run build_pca.')

        # Get the file to use when parsing the mask design information
        _design_file = (None if self.files is None else self.files[0]) if design_file is None \
                            else design_file
        if _design_file is None or not os.path.isfile(_design_file):
            msgs.error('Design file not found or none provided.')

        # Read the design data
        msgs.info('Reading slit-mask design information from: {0}'.format(_design_file))
        if self.spectrograph.get_slitmask(_design_file) is None:
            msgs.error('Unable to read design file or unable no slit-mask design reader '
                       'defined for {0}.'.format(self.spectrograph.spectrograph))

        # Try to match to both the left and right edges simultaneously
        x_mask = np.array([np.amin(self.spectrograph.slitmask.corners[:,:,0], axis=1),
                           np.amax(self.spectrograph.slitmask.corners[:,:,0], axis=1)]).T.ravel()
        x_trace = self.spat_fit[self.pca.reference_row,:]

        # Estimate the scale in pixels/mm as the telescope platescale
        # in arcsec/mm divided by the detector platescale in
        # arcsec/pixel
        pix_per_mm = self.spectrograph.telescope.platescale() \
                        / self.spectrograph.detector[self.det-1]['platescale']

        # Take a guess at the offset and set the bounds
#        try:
        # Try using the spectrograph detector map
        self.spectrograph.get_detector_map()
        # Set the offset based on the location of this detector
        offset = self.spectrograph.detector_map.npix[0]/2 \
                    - self.spectrograph.detector_map.ccd_center[self.det-1,0]
        # Set the bounds to some nominal fraction of the detector size and pix/mm scale
        bounds = [[offset-0.1*self.spectrograph.detector_map.npix[0],
                    offset+0.1*self.spectrograph.detector_map.npix[0]],
                    [pix_per_mm/1.1, pix_per_mm*1.1]]
        # except:
            # # No detector map
            # msgs.warn('No detector map available for {0}'.format(self.spectrograph.spectrograph)
                      # + '; attempting to match to slit-mask design anyway.')
            # # Set the guess offset such that two sets of coordinates
            # # are offset to their mean
            # offset = np.mean(x_trace) - np.mean(pix_per_mm * x_mask)
            # # Assume that the range of mask locations is larger than
            # # the range of trace locations; set the bounds to force the
            # # trace coordinates to overlap with the mask by ~90%
            # bounds = [[offset-(np.amin(x_trace)-np.amin(pix_per_mm*x_mask))*1.1,
                       # offset+(np.amax(x_trace) - np.amax(pix_per_mm*x_mask))*1.1], 
                      # [pix_per_mm/1.1, pix_per_mm*1.1]]

        import pdb
        pdb.set_trace()

        # Register the data
        register = slitmask.SlitRegister(x_trace, x_mask, guess=[offset, pix_per_mm],
                                         bounds=bounds, fit=True)

        pdb.set_trace()










#-----------------------------------------------------------------------
# Done with EdgeTraceSet class.  Below are "pypeit.core.trace" routines.
#-----------------------------------------------------------------------
def detect_slit_edges(flux, mask=None, median_iterations=0, min_sqm=30., sobel_mode='nearest',
                      sigdetect=30.):
    """
    Find slit edges using the input image.

    The primary algorithm is to run a Sobel filter on the image and
    then trigger on all significant gradients. Positive gradients are
    left edges, negative gradients are right edges.

    Args:
        flux (`numpy.ndarray`_):
            Calibration frame used to identify slit edges.  Likely a
            flat-field image that has been lightly smoothed in the
            spectral direction.  The image should also have its bad
            pixels replaced (see
            :func:`pypeit.core.procimg.replace_columns`).  Its
            orientation *must* have spectra dispersed along rows.
        mask (`numpy.ndarray`_, optional):
            A boolean or integer bad-pixel mask.  If None, all pixels
            are assumed valid.  This is used to ignore features in the
            image that may be due to bad pixels.
        median_iterations (:obj:`int`, optional):
            Number of median smoothing iteration to perform on the trace
            image.  The size of the smoothing is always (7,3).  For
            long-slit data, we recommend `median_iterations=0`.
        min_sqm (:obj:`float`, optional):
            Minimum error used when detecting a slit edge.  TODO: This
            needs a better description.
        sobel_mode (:obj:`str`, optional):
            Mode to use with the Sobel filter.  See
            `scipy.ndimage.sobel`_.
        sigdetect (:obj:`float`, optional):
            Threshold for edge detection.

    Returns:
        Returns two `numpy.ndarray`_ objects: (1) The image of the
        significance of the edge detection in sigma and (2) the array
        isolating the slit edges. In the latter, left edges have a
        value of -1 and right edges have a value of 1.
    """
    # Checks
    if flux.ndim != 2:
        msgs.error('Trace image must be 2D.')
    _mask = np.zeros_like(flux, dtype=int) if mask is None else mask.astype(int)
    if _mask.shape != flux.shape:
        msgs.error('Mismatch in mask and trace image shapes.')

    # Specify how many times to repeat the median filter.  Even better
    # would be to fit the filt/sqrt(abs(binarr)) array with a Gaussian
    # near the maximum in each column
    msgs.info("Detecting slit edges in the trace image")

    # Generate sqrt image
    sqmstrace = np.sqrt(np.abs(flux))

    # Median filter
    # TODO: Add size to parameter list
    for ii in range(median_iterations):
        sqmstrace = ndimage.median_filter(sqmstrace, size=(7, 3))

    # Make sure there are no spuriously low pixels
    sqmstrace[(sqmstrace < 1.0) & (sqmstrace >= 0.0)] = 1.0
    sqmstrace[(sqmstrace > -1.0) & (sqmstrace <= 0.0)] = -1.0

    # Filter with a Sobel
    filt = ndimage.sobel(sqmstrace, axis=1, mode=sobel_mode)
    # Apply the bad-pixel mask
    filt *= (1.0 - _mask)
    # Significance of the edge detection
    sobel_sig = np.sign(filt)*np.power(filt,2)/np.maximum(sqmstrace, min_sqm)

    # First edges assigned according to S/N
    # TODO: why not match the sign of the Sobel image to the edge it
    # traces? I.e., why is the sign flipped?
    tedges = np.zeros(flux.shape, dtype=np.float)
    tedges[np.where(sobel_sig > sigdetect)] = -1.0  # A positive gradient is a left edge
    tedges[np.where(sobel_sig < -sigdetect)] = 1.0  # A negative gradient is a right edge
    
    # Clean the edges
    wcl = np.where((ndimage.maximum_filter1d(sobel_sig, 10, axis=1) == sobel_sig) & (tedges == -1))
    wcr = np.where((ndimage.minimum_filter1d(sobel_sig, 10, axis=1) == sobel_sig) & (tedges == 1))
    edge_img = np.zeros(sobel_sig.shape, dtype=np.int)
    edge_img[wcl] = -1
    edge_img[wcr] = 1

    if mask is not None:
        msgs.info("Applying bad pixel mask")
        edge_img *= (1-_mask)
        sobel_sig *= (1-_mask)

    return sobel_sig, edge_img


def identify_traces(edge_img, max_spatial_separation=4, follow_span=10, minimum_spec_length=50):
    """
    Follow slit edges to identify unique slit traces.

    Args:
        edge_img (`numpy.ndarray`_):
            An array marked with -1 for left slit edges and +1 for right
            slit edges and 0 everywhere else.  The image *must* be
            oriented with the spatial dimension primarily along the
            first axis and spectral dimension primarily along the
            second.  See :func:`detect_slit_edges`.
        max_spatial_separation (:obj:`int`, optional):
            The maximum spatial separation between two edges in proximal
            spectral rows before they become separated into different
            slit traces.
        follow_span (:obj:`int`, optional):
            The number of previous spectral rows to consider when
            following slits forward.
        minimum_spec_length (:obj:`int`, optional):
            The minimum number of spectral rows in an edge trace.
            Traces that do not meet this criterion are ignored.

    Returns:
        `numpy.ndarray`_: An integer array with trace ID numbers at
        pixels locating the edge of that trace. Negative traces are
        for left edges, positive for right edges. The number of left
        and right traces can be determined using
        :func:`count_edge_traces`. Pixels not associated to any edge
        have a value of 0.
    """
    msgs.info('Finding unique traces among detected edges.')
    # Check the input
    if edge_img.ndim > 2:
        raise ValueError('Provided edge image must be 2D.')
    if not np.array_equal(np.unique(edge_img), [-1,0,1]):
        raise ValueError('Edge image must only have -1, 0, or 1 values.')

    # Find the left and right coordinates
    lx, ly = np.where(edge_img == -1)
    rx, ry = np.where(edge_img == 1)
    x = np.concatenate((lx, rx))
    # Put left traces at negative y
    y = np.concatenate((-ly, ry))

    # The trace ID to associate with each coordinate
    trace = np.full_like(x, -1)

    # Loop over spectral channels
    last = 0
    for row in range(np.amin(x), np.amax(x)+1):
        # Find the slit edges in this row
        indx = x == row
        in_row = np.sum(indx)
        if in_row == 0:
            # No slits found in this row
            continue

        # Find the unique edge y positions in the selected set of
        # previous rows and their trace IDs
        prev_indx = np.logical_and(x < row, x > row - follow_span)
        if not np.any(prev_indx):
            # This is likely the first row or the first row with any
            # slit edges
            trace[indx] = np.arange(in_row)+last
            last += in_row
            continue
        uniq_y, uniq_i = np.unique(y[prev_indx], return_index=True)
        uniq_t = trace[prev_indx][uniq_i]

        # Assign trace IDs to this row
        #   - First match to any previous IDs
        row_trace = np.full(in_row, -1)
        for i, _y in enumerate(y[indx]):
            dist = np.absolute(uniq_y-_y)
            mindist = np.argmin(dist)
            if dist[mindist] < max_spatial_separation:
                row_trace[i] = uniq_t[mindist]
        #   - Assign new trace IDs to unmatched edges
        unassigned = row_trace == -1
        n_unassigned = np.sum(unassigned)
        row_trace[unassigned] = np.arange(n_unassigned)+last
        last += n_unassigned
        #   - Assign all edges and continue
        trace[indx] = row_trace

    # Reorder the traces and remove any that do not meet the specified
    # length.
    #   - Left edges.  Given negative IDs starting with -1
    indx = y < 0
    left, reconstruct, counts = np.unique(trace[indx], return_inverse=True,
                                             return_counts=True)
#    if np.any(counts > edge_img.shape[0]):
#        warnings.warn('Some traces have more pixels than allowed by the image.  The maximum '
#                      'spatial separation for the edges in a given trace may be too large.')
    good_trace = counts > minimum_spec_length
    left[:] = 0
    left[good_trace] = -1-np.arange(np.sum(good_trace))
    trace[indx] = left[reconstruct]
    #   - Right edges.  Given positive IDs starting with 1
    indx = np.invert(indx)
    right, reconstruct, counts = np.unique(trace[indx], return_inverse=True,
                                              return_counts=True)
#    if np.any(counts > edge_img.shape[0]):
#        warnings.warn('Some traces have more pixels than allowed by the image.  The maximum '
#                      'spatial separation for the edges in a given trace may be too large.')
    good_trace = counts > minimum_spec_length
    right[:] = 0
    right[good_trace] = 1+np.arange(np.sum(good_trace))
    trace[indx] = right[reconstruct]

    # Construct the image with the trace IDs and return
    trace_id_img = np.zeros_like(edge_img, dtype=int)
    trace_id_img[x,np.absolute(y)] = trace
    return trace_id_img


def count_edge_traces(edge_img):
    """
    Count the number of left and right edges traced.

    Args:
        edge_img (`numpy.ndarray`_):
            Image with edge trace pixels numbered by their associated
            trace.  Pixels with positive numbers follow right slit edges
            and negative numbers follow left slit edges.
    
    Returns:
        Two integers with the number of left and right edges,
        respectively.
    """
    # Avoid returning -0
    nleft = np.amin(edge_img)
    return 0 if nleft == 0 else -nleft, np.amax(edge_img)


# TODO: This needs to be better tested
def atleast_one_edge(edge_img, mask=None, flux_valid=True, copy=False):
    """
    Ensure that there is at least one left and one right slit edge
    identified.

    This is especially useful for long slits that fill the full
    detector, e.g. Shane Kast.

    Args:
        edge_img (`numpy.ndarray`_):
            Image with edge trace pixels numbered by their associated
            trace.  Pixels with positive numbers follow right slit edges
            and negative numbers follow left slit edges.
        mask (`numpy.ndarray`_, optional):
            Integer (0 unmasked; 1 masked) or boolean array indicating
            bad pixels in the image.  If None, all pixels are considered
            good.
        flux_valid (:obj:`bool`, optional):
            The flux in the image used to construct the edge traces is
            valid meaning that any problems should not be an issue with
            the trace image itself.
        copy (:obj:`bool`, optional):
            Copy `edge_img` to a new array before making any
            modifications.  Otherwise, `edge_img` is modified in-place.
   
    Returns:
        `numpy.ndarray`_: The modified trace image, which is either a
        new array or points to the in-place modification of `edge_img`
        according to the value of `copy`.  If no slit edges were found
        and the flux in the trace image is invalid (`flux_valid=False`),
        function returns `None`.
    """
    # Get the number of traces
    nleft, nright = count_edge_traces(edge_img)

    # Determine whether or not to edit the image in place
    _edge_img = edge_img.copy() if copy else edge_img

    if nleft != 0 and nright != 0:
        # Don't need to add anything
        return _edge_img

    if nleft == 0 and nright == 0 and not flux_valid:
        # No traces and fluxes are invalid.  Warn the user and continue.
        msgs.warn('Unable to trace any edges!  Image flux is low; check trace image is correct.')
        return None

    # Use the mask to determine the first and last valid pixel column
    sum_bpm = np.ones(edge_img.shape[1]) if mask is None else np.sum(mask, axis=0) 

    if nleft == 0:
        # Add a left edge trace at the first valid column
        msgs.warn('No left edge found. Adding one at the detector edge.')
        gdi0 = np.min(np.where(sum_bpm == 0)[0])
        _edge_img[:,gdi0] = -1

    if nright == 0:
        # Add a right edge trace at the last valid column
        msgs.warn('No right edge found. Adding one at the detector edge.')
        gdi1 = np.max(np.where(sum_bpm == 0)[0])
        _edge_img[:,gdi1] = 1

    return _edge_img


# TODO: This needs to be better tested
def handle_orphan_edge(edge_img, sobel_sig, mask=None, flux_valid=True, copy=False):
    """
    In the case of single left/right traces and multiple matching
    traces, pick the most significant matching trace and remove the
    others.

    If *no* left and/or right edge is present, this will add one using
    :func:`atleast_one_edge`.

    Args:
        edge_img (`numpy.ndarray`_):
            Image with edge trace pixels numbered by their associated
            trace.  Pixels with positive numbers follow right slit edges
            and negative numbers follow left slit edges.
        sobel_sig (`numpy.ndarray`_):
            Image with the significance of the edge detection.  See
            :func:`detect_slit_edges`.
        mask (`numpy.ndarray`_, optional):
            Integer (0 unmasked; 1 masked) or boolean array indicating
            bad pixels in the image.  If None, all pixels are considered
            good.
        flux_valid (:obj:`bool`, optional):
            The flux in the image used to construct the edge traces is
            valid meaning that any problems should not be an issue with
            the trace image itself.
        copy (:obj:`bool`, optional):
            Copy `edge_img` to a new array before making any
            modifications. Otherwise, `edge_img` is modified
            in-place.

    Returns:
        `numpy.ndarray`_: The modified trace image, which is either a
        new array or points to the in-place modification of
        `edge_img` according to the value of `copy`.
    """
    # Get the number of traces
    nleft, nright = count_edge_traces(edge_img)

    if nleft == 0 or nright == 0:
        # Deal with no left or right edges
        _edge_img = atleast_one_edge(edge_img, mask=mask, flux_valid=flux_valid, copy=copy)
    else:
        # Just do basic setup
        _edge_img = edge_img.copy() if copy else edge_img

    if nleft != 1 and nright != 1 or nleft == 1 and nright == 1:
        # Nothing to do
        return _edge_img
    
    if nright > 1:
        # To get here, nleft must be 1.  This is mainly in here for
        # LRISb, which is a real pain..
        msgs.warn('Only one left edge, and multiple right edges.')
        msgs.info('Restricting right edge detection to the most significantly detected edge.')
        # Find the most significant right trace
        best_trace = np.argmin([-np.median(sobel_sig[_edge_img==t]) for t in range(nright)])+1
        # Remove the other right traces
        indx = _edge_img == best_trace
        _edge_img[(_edge_img > 0) & np.invert(indx)] = 0
        # Reset the number to a single right trace
        _edge_img[indx] = 1
        return _edge_img

    # To get here, nright must be 1.
    msgs.warn('Only one right edge, and multiple left edges.')
    msgs.info('Restricting left edge detection to the most significantly detected edge.')
    # Find the most significant left trace
    best_trace = np.argmax([np.median(sobel_sig[_edge_img == -t]) for t in range(nleft)])+1
    # Remove the other left traces
    indx = _edge_img == best_trace
    _edge_img[(_edge_img > 0) & np.invert(indx)] = 0
    # Reset the number to a single left trace
    _edge_img[indx] = 1

    return _edge_img


def most_common_trace_row(trace_mask):
    """
    Find the spectral position (row) that crosses the most traces.

    If provided the mask for a single trace, this just returns the
    median of the unmasked rows.

    Args:
        trace_mask (`numpy.ndarray`_):
            Mask for the trace data (True is bad; False is good). Can
            be a 1D array for a single trace or a 2D array with shape
            (nspec, ntrace) for multiple traces.

    Returns:
        :obj:`int`: The row that crosses the most valid trace data.
    """
    if trace_mask.ndim == 1:
        rows = np.where(np.invert(trace_mask))[0]
        return rows[rows.size//2]
    return Counter(np.where(np.invert(trace_mask))[0]).most_common(1)[0][0]


def prepare_sobel_for_trace(sobel_sig, boxcar=5, side='left'):
    """
    Prepare the Sobel filtered image for tracing.

    This method:
        - Flips the value of the Sobel image for the right traces;
        the pixels along right traces are negative.
        - Smooths along rows (spatially)

    Args:           
        sobel_sig (`numpy.ndarray`_):
            Image with the significance of the edge detection; see
            :func:`detect_slit_edges`.
        boxcar (:obj:`int`, optional):
            Boxcar smooth the detection image along rows before
            recentering the edge centers; see
            :func:`pypeit.utils.boxcar_smooth_rows`. If `boxcar` is
            less than 1, no smoothing is performed.
        side (:obj:`str`, optional):
            The side that the image will be used to trace. In the
            Sobel image, positive values are for left traces,
            negative for right traces. If 'left', the image is
            clipped at a minimum value of -0.1. If 'right', the image
            sign is flipped and then clipped at a minimum of -0.1. If
            None, the image is not flipped or clipped, only smoothed.

    Returns:
        `numpy.ndarray`_: The smoothed image.
    """
    if side not in ['left', 'right', None]:
        raise ValueError('Side must be left, right, or None.')
    if side is None:
        img = sobel_sig
    else:
        img = np.maximum(sobel_sig, -0.1) if side == 'left' \
                    else np.maximum(-1*sobel_sig, -0.1)
    return utils.boxcar_smooth_rows(img, boxcar) if boxcar > 1 else img


def follow_trace_moment(flux, start_row, start_cen, ivar=None, mask=None, fwgt=None, width=6.0,
                        maxshift_start=0.5, maxshift_follow=0.15, maxerror=0.2, continuous=True,
                        bitmask=None):
    """
    Follow a set of traces using a moment analysis of the provided
    image.

    Starting from a specified row and input centers along each
    column, attempt to follow a set of traces to both lower and
    higher rows in the provided image.

    Importantly, this function does not treat each row independently
    (as would be the case for a direct call to
    :func:`pypeit.moment.moment1d` for trace positions along all
    rows), but treats the calculation of the centroids sequentially
    where the result for each row is dependent and starts from the
    result from the previous row. The only independent measurement is
    the one performed at the input `start_row`. This function is much
    slower than :func:`pypeit.moment.moment1d` because of this
    introduced dependency.

    Args:
        flux (`numpy.ndarray`_):
            Image used to weight the column coordinates when
            recentering. This should typically be the Sobel filtered
            trace image after adjusting for the correct side and
            performing any smoothing; see
            :func:`prepare_sobel_for_trace`.
        start_row (:obj:`int`):
            Row at which to start the recentering. The function
            begins with this row and then traces first to higher
            indices and then to lower indices. This will almost
            certainly need to be different than the default, which is
            to start at the first row.
        start_cen (:obj:`int`, `numpy.ndarray`_, optional):
            One or more trace coordinates to recenter. If an array,
            must be 1D.
        ivar (`numpy.ndarray`_, optional):
            Inverse variance in the weight image. If not provided,
            unity variance is assumed. Used for the calculation of
            the errors in the moment analysis. If this is not
            provided, be careful with the value set for
            `maxerror_center` (see below).
        mask (`numpy.ndarray`_, optional):
            A boolean mask used to ignore pixels in the weight image.
            Pixels to ignore are masked (`mask==True`), pixels to
            analyze are not masked (`mask==False`)
        fwgt (`numpy.ndarray`_, optional):
            An additional weight to apply to each pixel in `flux`.  If
            None, weights are uniform.
        width (:obj:`float`, `numpy.ndarray`_, optional):
            The size of the window about the provided starting center
            for the moment integration window. See
            :func:`pypeit.moment.moment1d`.
        maxshift_start (:obj:`float`, optional):
            Maximum shift in pixels allowed for the adjustment of the
            first row analyzed, which is the row that has the most
            slit edges that cross through it.
        maxshift_follow (:obj:`float`, optional):
            Maximum shift in pixels between traces in adjacent rows
            as the routine follows the trace away from the first row
            analyzed.
        maxerror (:obj:`float`, optional):
            Maximum allowed error in the adjusted center of the trace
            returned by :func:`pypeit.moment.moment1d`.
        continuous (:obj:`bool`, optional):
            Keep only the continuous part of the traces from the
            starting row.
        bitmask (:class:`pypeit.bitmask.BitMask`, optional):
            Object used to flag traces. If None, assessments use
            boolean to flag traces. If not None, errors will be
            raised if the object cannot interpret the correct flag
            names defined. In addition to flags used by
            :func:`_recenter_trace_row`, this function uses the
            DISCONTINUOUS flag.

    Returns:
        Two numpy arrays are returned, the optimized center and an
        estimate of the error; the arrays are masked arrays if
        `start_cen` is provided as a masked array.
    """
    if flux.ndim != 2:
        raise ValueError('Input image must be 2D.')
    # Shape of the image with pixel weights
    nr, nc = flux.shape

    # Instantiate theses supplementary arrays here to speed things up
    # in iterative calling of moment1d. moment1d will check the array
    # sizes.
    _ivar = np.ones_like(flux, dtype=float) if ivar is None else ivar
    _mask = np.zeros_like(flux, dtype=bool) if mask is None else mask
    _fwgt = np.ones_like(flux, dtype=float) if fwgt is None else fwgt

    # Number of starting coordinates
    _cen = np.atleast_1d(start_cen)
    nt = _cen.size
    # Check coordinates are within the image
    if np.any((_cen > nc-1) | (_cen < 0)) or start_row < 0 or start_row > nr-1:
        raise ValueError('Starting coordinates incompatible with input image!')
    # Check the dimensionality
    if _cen.ndim != 1:
        raise ValueError('Input coordinates to be at most 1D.')

    # Instantiate output; just repeat input for all image rows.
    xc = np.tile(_cen, (nr,1)).astype(float)
    xe = np.zeros_like(xc, dtype=float)
    xm = np.zeros_like(xc, dtype=bool) if bitmask is None \
                else np.zeros_like(xc, dtype=bitmask.minimum_dtype())

    # Recenter the starting row
    i = start_row
    xc[i,:], xe[i,:], xm[i,:] = _recenter_trace_row(i, xc[i,:], flux, _ivar, _mask, _fwgt, width,
                                                    maxshift=maxshift_start, maxerror=maxerror,
                                                    bitmask=bitmask)

    # Go to higher indices using the result from the previous row
    for i in range(start_row+1,nr):
        xc[i,:], xe[i,:], xm[i,:] = _recenter_trace_row(i, xc[i-1,:], flux, _ivar, _mask, _fwgt,
                                                        width, maxshift=maxshift_follow,
                                                        maxerror=maxerror, bitmask=bitmask)

    # Go to lower indices using the result from the previous row
    for i in range(start_row-1,-1,-1):
        xc[i,:], xe[i,:], xm[i,:] = _recenter_trace_row(i, xc[i+1,:], flux, _ivar, _mask, _fwgt,
                                                        width, maxshift=maxshift_follow,
                                                        maxerror=maxerror, bitmask=bitmask)

    # NOTE: In edgearr_tcrude, skip_bad (roughly opposite of continuous
    # here) was True by default, meaning continuous would be False by
    # default.
    if not continuous:
        # Not removing discontinuous traces
        return xc, xe, xm

    # Keep only the continuous part of the trace starting from the
    # initial row
    # TODO: Instead keep the longest continuous segment?
    # TODO: Convert continuous from T/F to a tolerance that, e.g.,
    # allows for few-pixel discontinuities, but allows the trace to
    # pick back up again
    bad = xm > 0
    p = np.arange(nr)
    for i in range(nt):
        indx = bad[:,i] & (p > start_row)
        if np.any(indx):
            s = np.amin(p[indx])
            xm[s:,i] = True if bitmask is None else bitmask.turn_on(xm[s:,i], 'DISCONTINUOUS')
        indx = bad[:,i] & (p < start_row)
        if np.any(indx):
            e = np.amax(p[indx])-1
            xm[:e,i] = True if bitmask is None else bitmask.turn_on(xm[:e,i], 'DISCONTINUOUS')

    # Return centers, errors, and mask
    return xc, xe, xm


def _recenter_trace_row(row, cen, flux, ivar, mask, fwgt, width, maxshift=None, maxerror=None,
                        bitmask=None, fill_error=-1):
    """
    Recenter the trace along a single row.

    This method is not meant for general use; it is a support method
    for :func:`follow_trace_moment`. The method executes
    :func:`pypeit.moment.moment1d` (with uniform weighting) for the
    specified row (spectral position) in the image, imposing a
    maximum difference with the input and centroid error, and asseses
    the result. The assessments are provided as either boolean flags
    or mask bits, depending on the value of `bitmask` (see below).

    Measurements flagged by :func:`pypeit.moment.moment1d` or flagged
    as having a centroid error that is larger than `maxerror` (and
    `maxerror` is not None) are replaced by their input value.

    If `maxshift`, `maxerror`, and `bitmask` are all None, this is
    equivalent to::

        return moment.moment1d(flux, cen, width, ivar=ivar, mask=mask, row=row,
                               order=1, fill_error=fill_error)

    Args:
        row (:obj:`int`):
            Row (index along the first axis; spectral position) in
            `flux` at which to recenter the trace position. See
            `row` in :func:`pypeit.moment.moment1d`.
        cen (`numpy.ndarray`_):
            Current estimate of the trace center.
        flux (`numpy.ndarray`_):
            Array used for the centroid calculations.
        ivar (`numpy.ndarray`_):
            Inverse variance in `flux`; passed directly to
            :func:`pypeit.moment.moment1d` and can be None.
        mask (`numpy.ndarray`_):
            Boolean mask for `flux`; passed directly to
            :func:`pypeit.moment.moment1d` and can be None.
        fwgt (`numpy.ndarray`_):
            A weight to apply to each pixel in `flux`.
        width (:obj:`float`):
            Passed directly to :func:`pypeit.moment.moment1d`; see
            the documentation there.
        maxshift (:obj:`float`, optional):
            Maximum shift allowed between the input and recalculated
            centroid.  If None, no limit is applied.
        maxerror (:obj:`flaot`, optional):
            Maximum error allowed in the calculated centroid.
            Measurements with errors larger than this value are
            returned as the input center value. If None, no limit is
            applied.
        bitmask (:class:`pypeit.bitmask.BitMask`, optional):
            Object used to toggle the returned bit masks. If
            provided, must be able to interpret MATHERROR,
            MOMENTERROR, and LARGESHIFT flags. If None, the function
            returns boolean flags set to True if there was an error
            in :func:`pypeit.moment.moment1d` or if the error is
            larger than `maxerror` (and `maxerror` is not None);
            centroids that have been altered by the maximum shift are
            *not* flagged.

    Returns:
        Returns three `numpy.ndarray`_ objects: the new centers, the
        center errors, and the measurement flags with a data type
        depending on `bitmask`.
    """
    xfit, xerr, matherr = moment.moment1d(flux, cen, width, ivar=ivar, mask=mask, fwgt=fwgt,
                                          row=row, order=1, fill_error=fill_error)
#    xfit, xerr, matherr = recenter_moment(flux, cen, ivar=ivar, mask=mask, ycen=row,
#                                          width=width, fill_error=fill_error)
    if maxshift is None and maxerror is None and bitmask is None:
        # Nothing else to do
        return xfit, xerr, matherr

    # Toggle the mask bits
    if bitmask is not None:
        xmsk = np.zeros_like(xfit, dtype=bitmask.minimum_dtype())
        xmsk[matherr] = bitmask.turn_on(xmsk[matherr], 'MATHERROR')
        if maxerror is not None:
            indx = xerr > maxerror
            xmsk[indx] = bitmask.turn_on(xmsk[indx], 'MOMENTERROR')
        if maxshift is not None:
            indx = np.absolute(xfit - cen) > maxshift
            xmsk[indx] = bitmask.turn_on(xmsk[indx], 'LARGESHIFT')

    # Impose the maximum shift
    if maxshift is not None:
        xfit = np.clip(xfit - cen, -maxshift, maxshift)+cen

    # Reset 'bad' values to the input
    indx = matherr
    if maxerror is not None:
        indx |= (xerr > maxerror)
    xfit[indx] = cen[indx]
    xerr[indx] = fill_error

    # Return the new centers, errors, and flags
    return xfit, xerr, indx if bitmask is None else xmsk


def fit_trace(flux, trace, order, ivar=None, mask=None, trace_mask=None, weighting='uniform',
              fwhm=3.0, function='legendre', maxdev=5.0, maxiter=25, niter=9, show_fits=False,
              idx=None, xmin=None, xmax=None):
    """
    Iteratively fit the trace of a feature in the provided image.

    Each iteration performs two steps:
        - Redetermine the trace data using
          :func:`pypeit.moment.moment1d`. The size of the integration
          window (see the definition of the `width` parameter for
          :func:`pypeit.moment.moment1d`)depends on the type of
          weighting: For *uniform weighting*, the code does a third
          of the iterations with window `width = 2*1.3*fwhm`, a third
          with `width = 2*1.1*fhwm`, and a third with `width =
          2*fwhm`. For
          *Gaussian weighting*, all iterations use `width =
          fwhm/2.3548`.
        - Fit the centroid measurements with a 1D function of the
          provided order. See :func:`pypeit.core.pydl.TraceSet`.

    The number of iterations performed is set by the keyword argument
    `niter`. There is no convergence test, meaning that this number
    of iterations is *always* performed.

    History:
        23-June-2018  Written by J. Hennawi

    Args:
        flux (`numpy.ndarray`_):
            Image to use for tracing. Must be 2D with shape (nspec,
            nspat).
        trace (`numpy.ndarray`_):
            Initial guesses for spatial direction trace. This can
            either be an 2-d array with shape (nspec, nTrace) array,
            or a 1-d array with shape (nspec) for the case of a
            single trace.
        order (:obj:`int`):
            Order of function to fit to each trace.  See `function`.
        ivar (`numpy.ndarray`_, optional):
            Inverse variance of the image intensity.  If not provided,
            unity variance is used.  If provided, must have the same
            shape as `flux`.
        mask (`numpy.ndarray`_, optional):
            Boolean array with the input mask for the image. If not
            provided, all values in `flux` are considered valid. If
            provided, must have the same shape as `flux`.
        trace_mask (`numpy.ndarray`_, optional):
            Boolean array with the trace mask; i.e., places where you
            know the trace is going to be bad that you always want to
            mask in the fits. Shape must match `trace`.
        weighting (:obj:`str`, optional):
            The weighting to apply to the position within each
            integration window (see :func:`pypeit.moment.moment1d`).
        fwhm (:obj:`float`, optional):
            The expected width of the feature to trace, which is used
            to define the size of the integration window during the
            centroid calculation; see description above.
        function (:obj:`str`, optional):
            The name of the function to fit. Must be a valid
            selection. See :class`pypeit.core.pydl.TraceSet`.
        maxdev (:obj:`float`, optional):
            If provided, reject points with `abs(data-model) >
            maxdev` during the fitting. If None, no points are
            rejected. See :func:`pypeit.utils.robust_polyfit_djs`.
        maxiter (:obj:`int`, optional):
            Maximum number of rejection iterations allowed during the
            fitting. See :func:`pypeit.utils.robust_polyfit_djs`.
        niter (:obj:`int`, optional):
            The number of iterations for this method; i.e., the
            number of times the two-step fitting algorithm described
            above is performed.
        show_fits (:obj:`bool`, optional):
            Plot the data and the fits.
        idx (`numpy.ndarray`_, optional):
            Array of strings with the IDs for each object. Used only
            if show_fits is true for the plotting. Default is just a
            running number.
        xmin (:obj:`float`, optional):
            Lower reference for robust_polyfit polynomial fitting.
            Default is to use zero
        xmax (:obj:`float`, optional):
            Upper refrence for robust_polyfit polynomial fitting.
            Default is to use the image size in nspec direction

    Returns:
        Returns four `numpy.ndarray`_ objects and the
        :class:`pypeit.core.pydl.TraceSet` object with the
        best-fitting polynomial parameters for the traces. The four
        `numpy.ndarray`_ objects all have the same shape as the input
        positions (`trace`) and provide:

            - The best-fitting positions of each trace determined by
              the polynomial fit.
            - The centroids of the trace determined by either flux-
              or Gaussian-weighting, to which the polynomial is fit.
            - The errors in the centroids.
            - Boolean flags for each centroid measurement (see
              :func:`pypeit.moment.moment1d`).
    """
    # Ensure setup is correct
    if flux.ndim != 2:
        raise ValueError('Input image must be 2D.')
    if ivar is None:
        ivar = np.ones_like(flux, dtype=float)
    if ivar.shape != flux.shape:
        raise ValueError('Inverse variance array shape is incorrect.')
    if mask is None:
        mask = np.zeros_like(flux, dtype=bool)
    if mask.shape != flux.shape:
        raise ValueError('Mask array shape is incorrect.')
    if trace_mask is None:
        trace_mask = np.zeros_like(trace, dtype=bool)

    # Allow for single vectors as input as well:
    _trace = trace.reshape(-1,1) if trace.ndim == 1 else trace
    _trace_mask = trace_mask.reshape(-1, 1) if trace.ndim == 1 else trace_mask
    nspec, ntrace = _trace.shape
    if _trace.shape != _trace_mask.shape:
        raise ValueError('Trace data and its mask do not have the same shape.')

    # Define the fitting limits
    if xmin is None:
        xmin = 0.0
    if xmax is None:
        xmax = float(nspec-1)

    # Abscissa for fitting; needs to be float type when passed to
    # TraceSet
    trace_coo = np.tile(np.arange(nspec), (ntrace,1)).astype(float)

    # Setup the width to use for each iteration depending on the weighting used
    width = np.full(niter, 2*fwhm if weighting == 'uniform' else fwhm/2.3548, dtype=float)
    if weighting == 'uniform':
        width[:niter//3] *= 1.3
        width[niter//3:2*niter//3] *= 1.1

    trace_fit = np.copy(_trace)
    # Uniform weighting during the fit
    trace_fit_ivar = np.ones_like(trace_fit)

    for i in range(niter):
        # First recenter the trace using the previous trace fit/data
        trace_cen, trace_err, bad_trace = moment.moment1d(flux, trace_fit, width[i], ivar=ivar,  
                                                          mask=mask, weighting=weighting, order=1)

        # TODO: Update trace_mask with bad_trace?

        # Do not do any kind of masking based on the trace recentering
        # errors. Trace fitting is much more robust when masked pixels
        # are simply replaced by the input trace values.
        
        # Do not do weighted fits, i.e. uniform weights but set the
        # error to 1.0 pixel
        traceset = pydl.TraceSet(trace_coo, trace_cen.T, inmask=np.invert(_trace_mask.T),
                                 function=function, ncoeff=order, maxdev=maxdev, maxiter=maxiter,
                                 invvar=trace_fit_ivar.T, xmin=xmin, xmax=xmax)

        # TODO: Report iteration number and mean/stddev in difference
        # of coefficients with respect to previous iteration
        trace_fit = traceset.yfit.T

    # Plot the final fit if requested
    if show_fits:
        # Set the title based on the type of weighting used
        title_text = 'Flux Weighted' if weighting == 'uniform' else 'Gaussian Weighted'
        if idx is None:
            idx = np.arange(1,ntrace+1).astype(str)

        # Bad pixels have errors set to 999 and are returned to lie on
        # the input trace. Use this only for plotting below.
        for i in range(ntrace):
            plt.scatter(trace_coo[i,:], trace_cen[:,i], marker='o', color='k', s=30,
                        label=title_text + ' Centroid')
            plt.plot(trace_coo[i,:], _trace[:,i], color='g', zorder=25, linewidth=2.0,
                     linestyle='--', label='initial guess')
            plt.plot(trace_coo[i,:], trace_fit[:,i], c='red', zorder=30, linewidth=2.0,
                     label ='fit to trace')
            if np.any(bad_trace[:,i]):
                plt.scatter(trace_coo[i,bad_trace[:,i]], trace_fit[bad_trace[:,i],i], c='blue',
                            marker='+', s=50, zorder=20, label='masked points, set to init guess')
            if np.any(_trace_mask[:,i]):
                plt.scatter(trace_coo[i,_trace_mask[:,i]], trace_fit[_trace_mask[:,i],i],
                            c='orange', marker='s', s=30, zorder=20,
                            label='input masked points, not fit')

            plt.title(title_text + ' Centroid to object {0}.'.format(idx[i]))
            plt.ylim((0.995*np.amin(trace_fit[:,i]), 1.005*np.amax(trace_fit[:,i])))
            plt.xlabel('Spectral Pixel')
            plt.ylabel('Spatial Pixel')
            plt.legend()
            plt.show()

    # Returns the fit, the actual weighted traces and flag, and the TraceSet object
    return trace_fit, trace_cen, trace_err, bad_trace, traceset


def build_trace_mask(flux, trace, mask=None, boxcar=None, thresh=None, median_kernel=None):
    """
    Construct a trace mask.

    If no keyword arguments are provided, the traces are only masked
    when they land outside the bounds of the image.

    If both `boxcar` and `thresh` are provided, traces are also
    masked by extracting the provided image along the trace (see
    :func:`pypeit.moment.moment1d`) and flagging extracted values
    below the provided threshold.

    Args:
        flux (`numpy.ndarray`_):
            Image to use for tracing. Shape is expected to be (nspec,
            nspat).
        trace (`numpy.ndarray`_):
            Trace locations. Can be a 1D array for a single trace or a
            2D array with shape (nspec, ntrace) for multiple traces.
        mask (`numpy.ndarray`_, optional):
            Boolean array with the input mask for the image. If not
            provided, all values in `flux` are considered valid. If
            provided, must have the same shape as `flux`.
        boxcar (:obj:`float`, optional):
            The width of the extraction window used for all traces
            and spectral rows. If None, the trace mask will not
            consider the extracted flux.
        thresh (:obj:`float`, optional):
            The minimum valid value of the extraced flux used to mask
            the traces. If None, the trace mask will not consider the
            extracted flux.
        median_kernel (:obj:`int`, optional):
            The spectral width of the kernel to use with
            `scipy.signal.medfilt` to filter the *extracted* data
            before setting the trace mask based on the provided
            threshold. If None, the extracted data are not filtered
            before flagging data below the threshold.

    Returns:
        `numpy.ndarray`_: The boolean mask for the traces.
    """
    # Setup and ensure input is correct
    if flux.ndim != 2:
        raise ValueError('Input image must be 2D.')
    nspec, nspat = flux.shape
    if mask is None:
        mask = np.zeros_like(flux, dtype=bool)
    if mask.shape != flux.shape:
        raise ValueError('Mask array shape is incorrect.')
    _trace = trace.reshape(-1,1) if trace.ndim == 1 else trace
    if _trace.shape[0] != nspec:
        raise ValueError('Must provide trace position for each spectral pixel.')

    # Flag based on the trace positions
    trace_mask = (_trace < 0) | (_trace > nspat - 1)

    if boxcar is None or thresh is None:
        # Only flagging based on the trace positions
        return trace_mask

    # Get the extracted flux
    extract_flux = moment.moment1d(flux, _trace, boxcar, mask=mask)[0]
    if median_kernel is not None:
        # Median filter the extracted data
        extract_flux = signal.medfilt(extract_flux, kernel_size=(median_kernel,1))
    return trace_mask | (extract_flux < thresh)


def rectify_image(img, col, mask=None, ocol=None, max_ocol=None, extract_width=None,
                  mask_threshold=0.5):
    r"""
    Rectify the image by shuffling flux along columns using the
    provided column mapping.

    The image recification is one dimensional, treating each image
    row independently. It can be done either by a direct resampling
    of the image columns using the provided mapping of output to
    input column location (see `col` and
    :class:`pypeit.sampling.Resample`) or by an extraction along the
    provided column locations (see `extract_width`). The latter is
    generally faster; however, when resampling each row, the flux is
    explicitly conserved (see the `conserve` argument of
    :class:`pypeit.sampling.Resample`.

    Args:
        img (`numpy.ndarray`_):
            The 2D image to rectify. Shape is :math:`(N_{\rm row},
            N_{\rm col})`.
        col (`numpy.ndarray`_):
            The array mapping each output column to its location in
            the input image. That is, e.g., `col[:,0]` provides the
            column coordinate in `img` that should be rectified to
            column 0 in the output image. Shape is :math:`(N_{\rm
            row}, N_{\rm map})`.
        mask (`numpy.ndarray`_, optional):
            Boolean mask for pixels to ignore in input image. If
            None, no pixels are masked in the rectification. If
            provided, shape must match `img`.
        ocol (`numpy.ndarray`_, optional):
            The column in the output image for each column in `col`.
            If None, assume::

                ocol = np.arange(col.shape[1])

            These coordinates can fall off the output image (i.e.,
            :math:`<0` or :math:`\geq N_{\rm out,col}`), but those
            columns are removed from the output).
        max_ocol (:obj:`int`, optional):
            The last viable column *index* to include in the output
            image; ie., for an image with `ncol` columns, this should
            be `ncol-1`. If None, assume `max(ocol)`.
        extract_width (:obj:`float`, optional):
            The width of the extraction aperture to use for the image
            rectification. If None, the image recification is
            performed using :class:`pypeit.sampling.Resample` along
            each row.
        mask_threshold (:obj:`float`, optional):
            Either due to `mask` or the bounds of the provided `img`,
            pixels in the rectified image may not be fully covered by
            valid pixels in `img`. Pixels in the output image with
            less than this fractional coverage of an input pixel are
            flagged in the output.

    Returns:
        Two `numpy.ndarray`_ objects are returned both with shape
        `(nrow,max_ocol+1)`, the rectified image and its boolean
        mask.
    """
    # Check the input
    if img.ndim != 2:
        raise ValueError('Input image must be 2D.')
    if mask is not None and mask.shape != img.shape:
        raise ValueError('Image mask must match image shape.')
    _img = np.ma.MaskedArray(img, mask=mask)
    nrow, ncol = _img.shape
    if col.ndim != 2:
        raise ValueError('Column mapping array must be 2D.')
    if col.shape[0] != nrow:
        raise ValueError('Number of rows in column mapping array must match image to rectify.')
    _ocol = np.arange(col.shape[1]) if ocol is None else np.atleast_1d(ocol)
    if _ocol.ndim != 1:
        raise ValueError('Output column indices must be provided as a vector.')
    if _ocol.size != col.shape[1]:
        raise ValueError('Output column indices must match columns in column mapping array.')
    _max_ocol = np.amax(_ocol) if max_ocol is None else max_ocol

    # Use an aperture extraction to rectify the image
    if extract_width is not None:
        # Select viable columns
        indx = (_ocol >= 0) & (_ocol <= _max_ocol)
        # Initialize the output image as all masked
        out_img = np.ma.masked_all((nrow,_max_ocol+1), dtype=float)
        # Perform the extraction
        out_img[:,_ocol[indx]] = moment.moment1d(_img.data, col[:,indx], extract_width,
                                                 mask=_img.mask)[0]
        # Determine what fraction of the extraction fell off the image
        coo = col[:,indx,None] + np.arange(np.ceil(extract_width)).astype(int)[None,None,:] \
                    - extract_width/2
        in_image = (coo >= 0) | (coo < ncol)
        out_msk = np.sum(in_image, axis=2)/extract_width < mask_threshold
        # Return the filled numpy.ndarray and boolean mask
        return out_img.filled(0.0)/extract_width, out_img.mask | out_msk

    # Directly resample the image
    # 
    # `col` provides the column in the input image that should
    # resampled to a given position in the output image: the value of
    # the flux at img[:,col[:,0]] should be rectified to `outimg[:,0]`.
    # To run the resampling algorithm we need to invert this. That is,
    # instead of having a function that provides the output column as a
    # function of the input column, we want a function that provies the
    # input column as a function of the output column.

    # Instantiate the output image
    out_img = np.zeros((nrow,_max_ocol+1), dtype=float)
    out_msk = np.zeros((nrow,_max_ocol+1), dtype=bool)
    icol = np.arange(ncol)
    for i in range(nrow):
        # Get the coordinate vector of the output column
        _icol = interpolate.interp1d(col[i,:], _ocol, copy=False, bounds_error=False,
                                     fill_value='extrapolate', assume_sorted=True)(icol)
        # Resample it
        r = sampling.Resample(_img[i,:], x=_icol, newRange=[0,_max_ocol], newpix=_max_ocol+1,
                              newLog=False, conserve=True)
        # Save the resampled data
        out_img[i,:] = r.outy
        # Flag pixels
        out_msk[i,:] = r.outf < mask_threshold
    return out_img, out_msk


# TODO: Add an option where the user specifies the number of slits, and
# so it takes only the highest peaks from detect_lines
def peak_trace(flux, ivar=None, mask=None, trace=None, extract_width=None, smash_range=None,
               peak_thresh=100.0, peak_clip=None, trough=False, trace_median_frac=0.01,
               trace_thresh=10.0, fwhm_uniform=3.0, fwhm_gaussian=3.0, function='legendre',
               order=5, maxdev=5.0, maxiter=25, niter_uniform=9, niter_gaussian=6, debug=False):
    """
    Find and trace features in an image by identifying peaks/troughs
    after collapsing along the spectral axis.

    The image is either compressed directly or after rectification
    using the supplied `trace`. The provided trace *must* have the
    same shape as the input `flux` image and map each spatial
    position as a function of spectral position. This can be the
    output of :func:`pca_predict` where the provided coordinates are
    `np.arange(flux.shape[1])`; see also
    :func:`EdgeTracePCA.predict`. The rectification of the input
    `flux` is done using a boxcar extraction along the provided
    traces with a width of `extract_width`; see
    :func:`pypeit.moment.moment1d`.

    The (rectified) image is collapsed spectrally (see `smash_range`)
    giving the sigma-clipped mean flux as a function of spatial
    position. Peaks are then isolated in this vector (see
    :func:`pypeit.core.arc.detect_lines`).

    Traces that pass through these peak positions are then passed to
    two iterations of :func:`fit_trace`, which both remeasures the
    centroids of the trace and fits a polynomial to those trace data.
    The first iteration determines the centroids with uniform
    weighting, passing `fwhm=fwhm_uniform` to :func:`fit_trace`, and
    the second uses Gaussian weighting for the centroid measurements
    (passing `fwhm=fwhm_gaussian` to :func:`fit_trace`). The results
    of this second iteration of :func:`fit_trace` are the data
    returned.

    Troughs in the image can also be traced, which is done by
    flipping the sign of the image about its median and then
    repeating the "peak" finding and :func:`fit_trace` iterations. If
    troughs are fit, the traces are order with the set of peak traces
    first (the number of which is given by the last returned object
    of the function), followed by the trough traces.

    Args:
        flux (`numpy.ndarray`_):
            Image to use for tracing.
        ivar (`numpy.ndarray`_, optional):
            Inverse variance of the image intensity.  If not provided,
            unity variance is used.  If provided, must have the same
            shape as `flux`.
        mask (`numpy.ndarray`_, optional):
            Boolean array with the input mask for the image. If not
            provided, all values in `flux` are considered valid. If
            provided, must have the same shape as `flux`.
        trace (`numpy.ndarray`_, optional):
            Trace data that maps the spatial position of all spectra
            as a function of spectral row. For example, this can be
            the output of :func:`pca_predict` where the provided
            coordinates are `np.arange(flux.shape[1])`; see also
            :func:`EdgeTracePCA.predict`. This is used to rectify the
            input image so that spectra are identically organized
            along image rows. Shape *must* be identical to `flux`. If
            None, `flux` is assumed to be rectified on input.
        extract_width (:obj:`float`, optional):
            The width of the extract aperture to use when rectifying
            the flux image. If None, set to `fwhm_gaussian`.
        smash_range (:obj:`tuple`, optional):
            Spectral range to over which to collapse the input image
            into a 1D flux as a function of spatial position. This 1D
            vector is used to detect features for tracing. This is
            useful (and recommended) for definining the relevant
            detector range for data with spectra that do not span the
            length of the detector. The tuple gives the minimum and
            maximum in the fraction of the full spectral length
            (nspec). If None, the full image is collapsed.
        peak_thresh (:obj:`float, optional):
            The threshold for detecting peaks in the image. See the
            `input_thresh` parameter for
            :func:`pypeit.core.arc.detect_lines`.
        peak_clip (:obj:`float, optional):
            Sigma-clipping threshold used to clip peaks outside with
            aberrant amplitudes. If None, no clipping is performed.
        trough (:obj:`bool`, optional):
            Trace both peaks **and** troughs in the input image. This
            is done by flipping the value of the smashed image about
            its median value, such that troughs can be identified as
            peaks.
        trace_median_frac (:obj:`float`, optional):
            After rectification of the image and before refitting the
            traces, the rectified image is median filtered with a
            kernel width of trace_median_frac*nspec along the
            spectral dimension.
        trace_thresh (:obj:`float`, optional):
            After rectification and median filtering of the image
            (see `trace_median_frac`), values in the resulting image
            that are *below* this threshold are masked in the
            refitting of the trace using :func:`fit_trace`.
        fwhm_uniform (:obj:`float`, optional):
            The `fwhm` parameter to use when using uniform weighting
            in the calls to :func:`fit_trace`. See description of the
            algorithm above.
        fwhm_gaussian (:obj:`float`, optional):
            The `fwhm` parameter to use when using Gaussian weighting
            in the calls to :func:`fit_trace`. See description of the
            algorithm above.
        function (:obj:`str`, optional):
            The type of polynomial to fit to the trace data. See
            :func:`fit_trace`.
        order (:obj:`int`, optional):
            Order of the polynomial to fit to each trace. See
            :func:`fit_trace`.
        maxdev (:obj:`float`, optional):
            See :func:`fit_trace`. If provided, reject points with
            `abs(data-model) > maxdev` when fitting the trace. If
            None, no points are rejected.
        maxiter (:obj:`int`, optional):
            Maximum number of rejection iterations allowed during the
            fitting. See :func:`fit_trace`.
        niter_uniform (:obj:`int`, optional):
            The number of iterations for :func:`fit_trace` when edge
            measurements are based on uniform weighting. See
            description above.
        niter_gaussian (:obj:`int`, optional):
            The number of iterations for :func:`fit_trace` when edge
            measurements are based on Gaussian weighting. See
            description above.
        debug (:obj:`bool`, optional):
            Show plots useful for debugging.

    Returns:
        Returns four `numpy.ndarray`_ objects and the number of peak
        traces. The number of peak traces should be used to separate
        peak from trough traces; if `trough` is False, this will just
        be the total number of traces. The four
        `numpy.ndarray`_ objects provide:

            - The best-fitting positions of each trace determined by
              the polynomial fit.
            - The centroids of the trace determined by the
              Gaussian-weighting iteration, to which the polynomial
              is fit.
            - The errors in the Gaussian-weighted centroids.
            - Boolean flags for each centroid measurement (see
              :func:`pypeit.moment.moment1d`).
    """
    # Setup and ensure input is correct
    if flux.ndim != 2:
        raise ValueError('Input image must be 2D.')
    nspec, nspat = flux.shape
    if ivar is None:
        ivar = np.ones_like(flux, dtype=float)
    if ivar.shape != flux.shape:
        raise ValueError('Inverse variance array shape is incorrect.')
    if mask is None:
        mask = np.zeros_like(flux, dtype=bool)
    if mask.shape != flux.shape:
        raise ValueError('Mask array shape is incorrect.')

    # Define the region to collapse
    if smash_range is None:
        smash_range = (0,1)

    # Set the image to collapse
    if trace is None:
        # Assume image already rectified
        flux_extract = flux
        # Just set the trace to the follow the spatial columns
        trace = np.tile(np.arange(nspat), (nspec,1))
    else:
        # Check there is a trace for each image pixel
        if trace.shape != flux.shape:
            raise ValueError('Provided trace data must match the image shape.')
        msgs.info('Rectifying image by extracting along trace for each spatial pixel')
        # TODO: JFH What should this aperture size be? I think fwhm=3.0
        # since that is the width of the sobel filter
        flux_extract = rectify_image(flux, trace, mask=mask,
                                     extract_width=fwhm_gaussian if extract_width is None
                                                     else extract_width)[0]
#        if debug:
#            ginga.show_image(flux_extract, chname ='rectified image')

    # Collapse the image along the spectral direction to isolate peaks/troughs
    start, end = np.clip(np.asarray(smash_range)*nspec, 0, nspec).astype(int)
    msgs.info('Collapsing image spectrally between pixels {0}:{1}'.format(start, end))
    flux_smash_mean, flux_smash_median, flux_smash_sig \
            = sigma_clipped_stats(flux_extract[start:end,:], axis=0, sigma=4.0)

    # Offset by the median
    # TODO: If tracing Sobel-filtered image, this should be close to,
    # or identically, 0
    flux_median = np.median(flux_smash_mean)
    flux_smash_mean -= flux_median

    # Trace peak or both peaks and troughs
    label = ['peak', 'trough'] if trough else ['peak']
    sign = [1, -1] if trough else [1]

    # Instantiate output
    npeak = 0
    trace_fit = np.empty((nspec,0), dtype=float)
    trace_cen = np.empty((nspec,0), dtype=float)
    trace_err = np.empty((nspec,0), dtype=float)
    bad_trace = np.empty((nspec,0), dtype=bool)

    # Get the smoothing kernel width and ensure it is odd
    median_kernel = int(np.ceil(nspec*trace_median_frac))//2 * 2 + 1

    # Identify and trace features in the image
    for i,(l,s) in enumerate(zip(label,sign)):

        # Identify the peaks
        peak, _, cen, _, _, best, _, _ \
                = arc.detect_lines(s*flux_smash_mean, cont_subtract=False, fwhm=fwhm_gaussian,
                                   input_thresh=peak_thresh, max_frac_fwhm=4.0,
                                   min_pkdist_frac_fwhm=5.0, debug=debug)

        if len(cen) == 0 or not np.any(best):
            msgs.warn('No good {0}s found!'.format(l))
            continue
        msgs.info('Found {0} good {1}(s) in the rectified, collapsed image'.format(
                    len(cen[best]),l))

        # Set the reference spatial locations to use for tracing the detected peaks
        cen = cen[best]
        loc = np.round(cen).astype(int) 
        if peak_clip is not None:
            # Clip the peaks based on their amplitude as a stop-gap for
            # a detection threshold that may be too low.
            # TODO: Only clip low values?
            clipped_peak = sigma_clip(peak[best])
            peak_mask = np.ma.getmaskarray(clipped_peak)
            if np.any(peak_mask):
                msgs.warn('Clipping {0} detected peak(s) with aberrant amplitude(s).'.format(
                                np.sum(peak_mask)))
                loc = loc[np.invert(peak_mask)]
                cen = cen[np.invert(peak_mask)]

        # As the starting point for the iterative trace fitting, use
        # the input trace data at the positions of the detected peaks.
        # The trace at column `loc` is expected to pass through `loc`
        # at the reference spectral pixel. The offset of the trace is
        # to allow for non-integer measurements of the peak centroid.
        trace_peak = trace[:,loc] + (cen-loc)[None,:]

        # Image to trace: flip when tracing the troughs and clip low
        # values
        # TODO: This -1 is drawn out of the ether
        _flux = np.clip(s*(flux - flux_median), -1, None)

        # Construct the trace mask
        trace_peak_mask = build_trace_mask(_flux, trace_peak, mask=mask, boxcar=fwhm_gaussian,
                                           thresh=trace_thresh, median_kernel=median_kernel)

        # Remeasure and fit the trace using uniform weighting
        trace_peak, cen, err, bad, _ \
                = fit_trace(_flux, trace_peak, order, ivar=ivar, mask=mask,
                            trace_mask=trace_peak_mask, fwhm=fwhm_uniform, function=function,
                            maxdev=maxdev, maxiter=maxiter, niter=niter_uniform, show_fits=debug)

        # Reset the mask
        # TODO: Use or include `bad` resulting from fit_trace()?
        trace_peak_mask = build_trace_mask(_flux, trace_peak, mask=mask, boxcar=fwhm_gaussian,
                                           thresh=trace_thresh, median_kernel=median_kernel)

        # Redo the measurements and trace fitting with Gaussian
        # weighting
        trace_peak, cen, err, bad, _ \
                = fit_trace(_flux, trace_peak, order, ivar=ivar, mask=mask,
                            trace_mask=trace_peak_mask, weighting='gaussian', fwhm=fwhm_gaussian,
                            function=function, maxdev=maxdev, maxiter=maxiter,
                            niter=niter_gaussian, show_fits=debug)

        # Save the results
        trace_fit = np.append(trace_fit, trace_peak, axis=1)
        trace_cen = np.append(trace_cen, cen, axis=1)
        trace_err = np.append(trace_err, err, axis=1)
        bad_trace = np.append(bad_trace, bad, axis=1)

        if i == 0:
            # Save the number of peaks (troughs are appended, if they're located)
            npeak = cen.shape[1]

    return trace_fit, trace_cen, trace_err, bad_trace, npeak

#-----------------------------------------------------------------------
# TODO: These should go in a pypeit.core.pca module
def pca_decomposition(vectors, npca=None, pca_explained_var=99.0, mean=None):
    r"""
    Perform principle-component analysis (PCA) for a set of 1D vectors.

    The vectors are first passed to an unconstrained PCA to determine
    the growth curve of the accounted variance as a function of the
    PCA component. If specifying a number of PCA components to use
    (see `npca`), this yields the percentage of the variance
    accounted for in the analysis. If instead specifying the target
    variance percentage (see `pca_explained_var`), this is used to
    determine the number of PCA components to use in the final
    analysis.

    Args:
        vectors (`numpy.ndarray`_):
            A 2D array with vectors to analyze with shape
            :math:`(N_{\rm vec}, N_{\rm pix})`. All vectors must be
            the same length and cannot be masked.
        npca (:obj:`bool`, optional):
            The number of PCA components to keep, which must be less
            than :math:`N_{\rm vec}`. If `npca==nvec`, no PCA
            compression occurs. If None, `npca` is automatically
            determined by calculating the minimum number of
            components required to explain a given percentage of
            variance in the data. (see `pca_explained_var`).
        pca_explained_var (:obj:`float`, optional):
            The percentage (i.e., not the fraction) of the variance
            in the data accounted for by the PCA used to truncate the
            number of PCA coefficients to keep (see `npca`). Ignored
            if `npca` is provided directly.
        mean (`numpy.ndarray`_, optional):
            The mean value of each vector to subtract from the data
            before performing the PCA. If None, this is determined
            directly from the data. Shape must be :math:`N_{\rm
            vec}`.

    Returns:
        Returns four `numpy.ndarray`_ objects:
            - The coefficients of each PCA component, `coeffs`. Shape
              is :math:`(N_{\rm vec},N_{\rm comp})`.
            - The PCA component vectors, `components`. Shape is
              :math:`(N_{\rm comp},N_{\rm pix})`.
            - The mean offset of each PCA for each pixel, `pca_mean`.
              Shape is :math:`(N_{\rm pix},)`.
            - The mean offset applied to each vector before the PCA,
              `vec_mean`. Shape is :math:`(N_{\rm vec},)`.

        To reconstruct the PCA representation of the input vectors, compute::

            np.dot(coeffs, components) + pca_mean[None,:] + vec_mean[:,None]

    """
    # Check input
    if vectors.ndim != 2:
        raise ValueError('Input trace data must be a 2D array')
    nvec = vectors.shape[0]
    if nvec < 2:
        raise ValueError('There must be at least 2 vectors for the PCA analysis.')

    # Take out the mean value of each vector
    if mean is None:
        mean = np.mean(vectors, axis=1)
    vec_pca = vectors - mean[:,None]

    # Perform unconstrained PCA of the vectors
    pca = PCA()
    pca.fit(vec_pca)

    # Compute the cumulative distribution of the variance explained by the PCA components.
    # TODO: Why round to 6 decimals?  Why work in percentages?
    var_growth = np.cumsum(np.round(pca.explained_variance_ratio_, decimals=6) * 100)
    # Number of components for a full decomposition
    npca_tot = var_growth.size

    msgs.info('The unconstrained PCA yields {0} components.'.format(npca_tot))
    if npca is None:
        # Assign the number of components to use based on the variance
        # percentage
        if pca_explained_var is None:
            raise ValueError('Must provide percentage explained variance.')
        npca = int(np.ceil(np.interp(pca_explained_var, var_growth, np.arange(npca_tot)+1))) \
                    if var_growth[0] < pca_explained_var else 1
    elif npca_tot < npca:
        raise ValueError('Too few vectors for a PCA of the requested dimensionality.  '
                         'The full (uncompressing) PCA has {0} component(s)'.format(npca_tot)
                         + ', which is less than the requested {0} component(s).'.format(npca)
                         + '  Lower the number of requested PCA component(s) or turn off the PCA.')

    msgs.info('PCA will include {0} component(s), '.format(npca)
              + 'containing {0:.3f}% of the total variance.'.format(var_growth[npca-1]))

    # Determine the PCA coefficients with the revised number of
    # components, and return the results
    pca = PCA(n_components=npca)
    pca_coeffs = pca.fit_transform(vec_pca)
    return pca_coeffs, pca.components_, pca.mean_, mean


def fit_pca_coefficients(coeff, order, function='legendre', lower=3.0, upper=3.0, minx=None,
                         maxx=None, maxrej=1, maxiter=25, coo=None, debug=False):
    r"""
    Fit a parameterized function to a set of PCA coefficients,
    primarily for the purpose of predicting coefficients at
    intermediate locations.

    The coefficients of each PCA component are fit by a low-order
    polynomial, where the abscissa is set by the `coo` argument (see
    :func:`pypeit.utils.robust_polyfit_djs`).

    .. note::
        This is a general function, not really specific to the PCA;
        and is really just a wrapper for
        :func:`pypeit.utils.robust_polyfit_djs`.

    Args:
        coeff (`numpy.ndarray`_):
            PCA component coefficients. If the PCA decomposition used
            :math:`N_{\rm comp}` components for :math:`N_{\rm vec}`
            vectors, the shape of this array must be :math:`(N_{\rm
            comp}, N_{\rm vec})`. The array can be 1D with shape
            :math:`(N_{\rm vec},)` if there was only one PCA
            component.
        order (:obj:`int`, `numpy.ndarray`_):
            The order, :math:`o`, of the function used to fit the PCA
            coefficients. Can be a single number for all PCA
            components, or an array with an order specific to each
            component. If the latter, the shape must be
            :math:`(N_{\rm comp},)`.
        function (:obj:`str`, optional):
            Type of function used to fit the data.
        lower (:obj:`float`, optional):
            Number of standard deviations used for rejecting data
            **below** the mean residual. If None, no rejection is
            *performed. See
            :func:`utils.robust_polyfit_djs`.
        upper (:obj:`float`, optional):
            Number of standard deviations used for rejecting data
            **above** the mean residual. If None, no rejection is
            *performed. See
            :func:`utils.robust_polyfit_djs`.
        minx, maxx (:obj:`float`, optional):
            Minimum and maximum values used to rescale the
            independent axis data. If None, the minimum and maximum
            values of `coo` are used. See
            :func:`utils.robust_polyfit_djs`.
        maxrej (:obj:`int`, optional):
            Maximum number of points to reject during fit iterations.
            See :func:`utils.robust_polyfit_djs`.
        maxiter (:obj:`int`, optional):
            Maximum number of rejection iterations allows. To force
            no rejection iterations, set to 0.
        coo (`numpy.ndarray`_, optional):
            Floating-point array with the independent coordinates to
            use when fitting the PCA coefficients. If None, simply
            uses a running number. Shape must be :math:`(N_{\rm
            vec},)`.
        debug (:obj:`bool`, optional):
            Show plots useful for debugging.

    Returns:
        Returns four objects:
            - A boolean `numpy.ndarray`_ masking data (`coeff`) that
            were rejected during the polynomial fitting. Shape is the
            same as the input `coeff`.
            - A `list` of `numpy.ndarray`_ objects (or a single
            `numpy.ndarray`_), one per PCA component where the length
            of the 1D array is the number of coefficients fit to the
            PCA-component coefficients. The number of function
            coefficients is typically :math:`N_{\rm coeff} = o+1`.
            - The minimum and maximum coordinate values used to
            rescale the abscissa during the fitting.
    """
    # Check the input
    #   - Get the shape of the input data to fit
    _coeff = np.atleast_2d(coeff)
    if _coeff.ndim != 2:
        raise ValueError('Array with coefficiencts cannot be more than 2D')
    nvec, npca = _coeff.shape
    #   - Set the abscissa of the data if not provided and check its
    #   shape
    if coo is None:
        coo = np.arange(nvec, dtype=float)
    if coo.size != nvec:
        raise ValueError('Vector coordinates have incorrect shape.')
    #   - Check the order of the functions to fit
    _order = np.atleast_1d(order)
    if _order.size == 1:
        _order = np.full(npca, order, dtype=int)
    if _order.size != npca:
        raise ValueError('Function order must be a single number or one number per PCA component.')
    #   - Force the values of minx and maxx if they're not provided directly
    if minx is None:
        minx = np.amin(coo)
    if maxx is None:
        maxx = np.amax(coo)

    # Instantiate the output
    coeff_used = np.ones(_coeff.shape, dtype=bool)
    fit_coeff = [None]*npca

    # Fit the coefficients of each PCA component so that they can be
    # interpolated to other coordinates.
    for i in range(npca):
        coeff_used[:,i], fit_coeff[i] \
                = utils.robust_polyfit_djs(coo, _coeff[:,i], _order[i], function=function,
                                           maxiter=maxiter, lower=lower, upper=upper,
                                           maxrej=maxrej, sticky=False, minx=minx, maxx=maxx)
        if debug:
            # Visually check the fits
            xvec = np.linspace(np.amin(coo), np.amax(coo), num=100)
            rejected = np.invert(coeff_used[:,i])
            plt.scatter(coo, _coeff[:,i], marker='.', color='k', s=100, facecolor='none',
                        label='pca coeff')
            if np.any(rejected):
                plt.scatter(coo[rejected], _coeff[rejected,i], marker='x', color='C3', s=80, 
                            label='robust_polyfit_djs rejected')
            plt.plot(xvec, utils.func_val(fit_coeff[i], xvec, function, minx=minx, maxx=maxx),
                     linestyle='--', color='C0',
                     label='Polynomial fit of order={0}'.format(_order[i]))
            plt.xlabel('Trace Coordinate', fontsize=14)
            plt.ylabel('PCA Coefficient', fontsize=14)
            plt.title('PCA Fit for Dimension #{0}/{1}'.format(i+1, npca))
            plt.legend()
            plt.show()

    # Return arrays that match the shape of the input data
    if coeff.ndim == 1:
        return np.invert(coeff_used)[0], fit_coeff[0], minx, maxx
    return np.invert(coeff_used), fit_coeff, minx, maxx


def pca_predict(x, pca_coeff_fits, pca_components, pca_mean, mean, function='legendre'):
    r"""
    Use a model of the PCA coefficients to predict vectors at the
    specified coordinates.

    Args:
        x (:obj:`float`, `numpy.ndarray`_):
            One or more trace coordinates at which to sample the PCA
            coefficients and produce the PCA model.
        pca_coeff_fits (:obj:`list`, `numpy.ndarray`_): 
            A `list` of `numpy.ndarray`_ objects (or a single
            `numpy.ndarray`_), one per PCA component where the length
            of the 1D array is the number of coefficients fit to the
            PCA-component coefficients.
        pca_mean (`numpy.ndarray`_):
            The mean offset of the PCA decomposotion for each pixel.
            Shape is :math:`(N_{\rm pix},)`.
        mean (`numpy.ndarray`_):
            The mean offset applied to each vector before the PCA.
            Shape is :math:`(N_{\rm vec},)`.
    
    Returns:
        `numpy.ndarray`_: PCA constructed vectors, one per position
        `x`. Shape is either :math:`(N_{\rm pix},)` or :math:`(N_{\rm
        x},N_{\rm pix})`, depending on the input shape/type of `x`.
    """
    _x = np.atleast_1d(x)
    if _x.ndim != 1:
        raise ValueError('Coordinates for predicted vectors must be no more than 1D.')
    # Calculate the coefficients using the best fitting function
    npca = pca_components.shape[0]
    c = np.zeros((_x.size, npca), dtype=float)
    for i in range(npca):
        c[:,i] = utils.func_val(pca_coeff_fits[i], _x, function)
    # Calculate the predicted vectors and return them
    vectors = np.dot(c, pca_components) + pca_mean[None,:] + mean[:,None]
    return vectors if _x.size > 1 else vectors[0,:]

