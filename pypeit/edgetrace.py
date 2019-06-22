"""
Implements edge tracing.

The primary purpose of this module is to provide the classes/methods
used to trace slit edges. To start tracing, you first need a trace
image. For example, given a `pypeit_file`, you can generate the trace
image as follows::

    # Get the path to the pypeit file
    redux_path = os.path.split(pypeit_file)[0]

    # Instantiate the reduction
    rdx = PypeIt(pypeit_file, redux_path=redux_path, show=True)

    # Select a dectector (if necessary)
    rdx.par['rdx']['detnum'] = 3
    det = 3

    # Find the first science frame
    frame_indx = numpy.arange(len(rdx.fitstbl))
    in_grp = rdx.fitstbl.find_calib_group(0)
    is_science = rdx.fitstbl.find_frames('science')
    grp_science = frame_indx[is_science & in_grp]
    combid = numpy.unique(rdx.fitstbl['comb_id'][grp_science])[0]
    calib_ID = int(rdx.fitstbl['calib'][grp_science][0])

    # Get all the frames and background frames for this combination
    # group
    frames = numpy.where(rdx.fitstbl['comb_id'] == combid)[0]
    bg_frames = numpy.where(rdx.fitstbl['bkg_id'] == combid)[0]

    # Run the calibrations
    rdx.caliBrate.set_config(frames[0], 3, rdx.par['calibrations'])

    # Run the steps needed for creating the trace image
    rdx.caliBrate.get_bpm()
    rdx.caliBrate.get_bias()

    # Get the trace files
    trace_rows = rdx.fitstbl.find_frames('trace', calib_ID=calib_ID, index=True)
    trace_image_files = rdx.fitstbl.frame_paths(trace_rows)
    trace_master_key = rdx.fitstbl.master_key(trace_rows[0], det=det)
    binning = rdx.fitstbl['binning'][trace_rows[0]]

    # Build the trace image
    traceImage = traceimage.TraceImage(rdx.spectrograph, files=trace_image_files, det=det,
                                       par=rdx.caliBrate.par['traceframe'])
    traceImage.process(bias_subtract=rdx.caliBrate.msbias, trim=rdx.caliBrate.par['trim'],
                       apply_gain=True)

With the trace image, you then instantiate the edge tracing object as follows::

    # Instantiate edge tracing
    return edges.EdgeTraceSet(rdx.spectrograph, rdx.caliBrate.par['slitedges'],
                              master_key=trace_master_key, master_dir=rdx.caliBrate.master_dir,
                              qa_path=rdx.caliBrate.qa_path, img=traceImage,
                              bpm=rdx.caliBrate.msbpm, det=det)

.. todo::
    - Should simplify example for how to generate a trace image

.. _numpy.ndarray: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
"""
import os
import time
import inspect
from collections import OrderedDict

import numpy as np

from scipy import ndimage

import matplotlib
from matplotlib import pyplot as plt
from matplotlib import ticker, rc

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy import table

from pypeit.bitmask import BitMask
from pypeit.par.pypeitpar import EdgeTracePar
from pypeit.core import parse, pydl, procimg, pca, trace

from pypeit import msgs
from pypeit import utils
from pypeit import sampling
from pypeit import ginga
from pypeit import masterframe
from pypeit import io
from pypeit.traceimage import TraceImage
from pypeit.spectrographs.util import load_spectrograph
from pypeit.spectrographs import slitmask

class EdgeTraceBitMask(BitMask):
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
              ('OUTSIDEAPERTURE', 'Recentering yielded a centroid outside the moment aperture'),
                   ('EDGEBUFFER', 'Recentering yielded a centroid too close to the detector edge'),
                ('DISCONTINUOUS', 'Pixel included in a trace but part of a diseontinuous segment'),
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
        super(EdgeTraceBitMask, self).__init__(list(mask.keys()), descr=list(mask.values()))

    @property
    def bad_flags(self):
        """
        List the flags that mean the trace is bad.
        """
        return list(set(self.bits.keys()) - set(self.insert_flags))

    @property
    def insert_flags(self):
        """
        List of flags used to mark traces inserted for various
        reasons.
        """
        return ['USERINSERT', 'SYNCINSERT', 'MASKINSERT'] 


class EdgeTracePCA:
    r"""
    Class to build and interact with PCA model of edge traces.

    This is primarily a container class for the results of
    :func:`pypeit.core.pca.pca_decomposition`,
    :func:`pypeit.core.pca.fit_pca_coefficients`, and
    :func:`pypeit.core.pca.pca_predict`.

    Args:
        trace_cen (`numpy.ndarray`_):
            A floating-point array with the spatial location of each
            each trace edge. Shape is :math:`(N_{\rm spec}, N_{\rm
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
                = pca.fit_pca_coefficients(self.pca_coeffs, order, function=self.function,
                                           lower=lower, upper=upper, minx=minx, maxx=maxx,
                                           maxrej=maxrej, maxiter=maxiter, coo=self.trace_coo,
                                           debug=debug)

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
        return pca.pca_predict(x, self.fit_coeff, self.pca_components, self.pca_mean, x,
                               function=self.function).T



class EdgeTraceSet(masterframe.MasterFrame):
    r"""
    Core class that identifies, traces, and pairs edges in an image
    to define the slit apertures.

    

    TODO: MORE

    img, bpm, and det should be considered mutually exclusive compare to load

    load takes precedence.  I.e., if both img and load are provided, img is ignored!

    if img is provided, the initialization will also run
    :func:`initial_trace` or :func:`auto_trace`, depending on the
    value of `auto`, *and* save the output.

    Nominal run:
        - initial_trace
        - centroid_refine
        - fit_refine (calls fit_trace)
        - build_pca
        - peak_refine (calls peak_trace, which uses both pca_trace and fit_trace)

    Final trace is based on a run of fit_refine that pass through the detected peaks

    design and object data are only available if the spectrograph
    class has a get_slitmask function, and that the slit mask data
    includes the object information.

    TODO: Write a script that uses the empty function to document the
    data model for these two tables. And for EdgeTraceSet as a whole?

    TODO: Describe I/O interface

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
        img (`numpy.ndarray`_, :class:`pypeit.traceimage.TraceImage`, optional):
            Two-dimensional image used to trace slit edges. If a
            :class:`pypeit.traceimage.TraceImage` is provided, the
            raw files used to construct the image are saved.
        bpm (`numpy.ndarray`_, optional):
            Bad-pixel boolean mask for the trace image. Must have the
            same shape as `img`. If None, all pixels are assumed to
            be valid.
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
            If a trace image is provided (`img`), run
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
            `img` in :func:`initial_trace` or :func:`auto_trace` is a
            :class:`pypeit.traceimage.TraceImage` object.
        img (`numpy.ndarray`_):
            Image data used to trace slit edges.
        bpm (`numpy.ndarray`_):
            Boolean bad-pixel mask for the trace image; must have the
            same shape as `img`.
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
            trace centroid; see :class:`EdgeTraceBitMask`.
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
#    master_type = 'Trace'   # For MasterFrame base
    master_type = 'Edges'
    bitmask = EdgeTraceBitMask()    # Object used to define and toggle tracing mask bits
    def __init__(self, spectrograph, par, master_key=None, master_dir=None, qa_path=None,
                 img=None, bpm=None, det=1, binning=None, auto=False, save=False, load=False):

        # TODO: It's possible for the master key and the detector
        # number to be inconsistent...
        masterframe.MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
                                         master_key=master_key)

        self.spectrograph = spectrograph    # Spectrograph used to take the data
        self.par = par                      # Parameters used for slit edge tracing

        self.files = None               # Files used to construct the trace image
        self.img = None                 # The image used to find the slit edges
        self.bpm = None                 # Mask for the trace image
        self.det = None                 # Detector used for the trace image
        self.binning = None             # Detector ordered spectral then spatial
        self.sobel_sig = None           # Sobel filtered image used to detect edges
        self.sobel_sig_left = None      # Sobel filtered image used to trace left edges
        self.sobel_sig_right = None     # Sobel filtered image used to trace right edges
        # TODO: Need a separate mask for the sobel image?
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

        if img is not None and load:
            raise ValueError('Arguments img and load are mutually exclusive.  Choose to '
                             'either trace a new image or load a previous trace.')

        if load:
            # Attempt to load an existing master frame
            self.load()
        elif img is not None:
            # Provided a trace image so instantiate the object.
            if auto:
                self.auto_trace(img, bpm=bpm, det=det, binning=binning, save=save)
            else:
                self.initial_trace(img, bpm=bpm, det=det, binning=binning, save=save)

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
                    table.Column(name='TRACELPIX', dtype=float, length=length,
                                 description='Spatial pixel coordinate for left edge'),
                    table.Column(name='TRACERPIX', dtype=float, length=length,
                                 description='Spatial pixel coordinate for right edge'),
                    table.Column(name='SLITID', dtype=int, length=length,
                                 description='Slit ID Number'),
                    table.Column(name='SLITLFOC', dtype=float, length=length,
                                 description='Left edge of the slit in mm at the focal plane'),
                    table.Column(name='SLITRFOC', dtype=float, length=length,
                                 description='Right edge of the slit in mm at the focal plane'),
                    table.Column(name='SLITRA', dtype=float, length=length,
                                 description='Right ascension of the slit center (deg)'),
                    table.Column(name='SLITDEC', dtype=float, length=length,
                                 description='Declination of the slit center (deg)'),
                    table.Column(name='SLITLEN', dtype=float, length=length,
                                 description='Slit length (arcsec)'),
                    table.Column(name='SLITWID', dtype=float, length=length,
                                 description='Slit width (arcsec)'),
                    table.Column(name='SLITPA', dtype=float, length=length,
                                 description='Slit position angle onsky (deg from N through E)'),
                    table.Column(name='ALIGN', dtype=np.int16, length=length,
                                 description='Slit used for alignment, not target observations.')
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
                    table.Column(name='OBJRA', dtype=float, length=length,
                                 description='Right ascension of the object (deg)'),
                    table.Column(name='OBJDEC', dtype=float, length=length,
                                 description='Declination of the object (deg)'),
                    table.Column(name='SLITID', dtype=int, length=length,
                                 description='Slit ID Number'),
                    table.Column(name='SLITINDX', dtype=int, length=length,
                                 description='Row index of relevant slit in the design table')
                           ])

    def rectify(self, flux, bpm=None, extract_width=None, mask_threshold=0.5):
        r""""
        Rectify the provided image based on the current edge trace
        PCA model.

        The is primarily a wrapper for
        :func:`pypeit.sampling.rectify_image`; see its documentation
        for more detail.

        Args:
            flux (`numpy.ndarray`_):
                The 2D image to rectify. Its shape should match the
                image used to construct the edge traces:
                :math:`(N_{\rm spec}, N_{\rm spat})`.
            bpm (`numpy.ndarray`_, optional):
                Boolean bad-pixel mask for pixels to ignore in input
                image. If None, no pixels are masked in the
                rectification. If provided, shape must match `flux`.
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
             image and its boolean bad-pixel mask.
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
        return sampling.rectify_image(flux, self.pca.predict(ocol), bpm=bpm, ocol=ocol,
                                      max_ocol=self.nspat-1, extract_width=extract_width,
                                      mask_threshold=mask_threshold)

    def auto_trace(self, img, bpm=None, det=1, binning=None, save=False, debug=False):
        r"""
        Execute a fixed series of methods to automatically identify
        and trace slit edges.

        Args:
            img (`numpy.ndarray`_, :class:`pypeit.traceimage.TraceImage`):
                2D image used to trace slit edges. If a
                :class:`pypeit.traceimage.TraceImage` is provided,
                the raw files used to construct the image and on-chip
                binning are saved; the latter overrides any directly
                provided `binning`. The array should have shape
                :math:`(N_{\rm spec},N_{\rm spat})`; i.e., spectra
                are ordered along columns.
            bpm (`numpy.ndarray`_, optional):
                Bad-pixel mask for the trace image. Must have the
                same shape as `img`. If None, all pixels are assumed
                to be valid.
            det (:obj:`int`, optional):
                The 1-indexed detector number that provided the trace
                image. This is *only* used to determine whether or
                not bad columns in the image are actually along
                columns or along rows, as determined by
                :attr:`spectrograph` and the result of a call to
                :func:`pypeit.spectrograph.Spectrograph.raw_is_transposed`.
            binning (:obj:`str`, optional):
                On-detector binning of the data ordered spectral then
                spatial with format, e.g., `2,1`. Ignored if `img` is
                an instance of :class:`pypeit.traceimage.TraceImage`.
            save (:obj:`bool`, optional):
                Save the result to the master frame.
            debug (:obj:`bool`, optional):
                Run in debug mode.
        """
        self.initial_trace(img, bpm=bpm, det=det, binning=binning, save=False)
        self.centroid_refine()
        self.fit_refine(debug=debug)
        self.pca_refine(debug=debug)
        self.peak_refine(rebuild_pca=True, debug=debug)
        self.sync()
        # TODO: Add mask_refine() when it's ready
        self.log += [inspect.stack()[0][3]]
        if save:
            self.save()

    def initial_trace(self, img, bpm=None, det=1, binning='1,1', save=False):
        r"""
        Initialize the object for tracing a new image.

        This effectively reinstantiates the object and must be the
        first method called for tracing an image.  The algorithm:
            - Lightly boxcar smooths the trace image spectrally.
            - Replaces bad pixel columns, if a mask is provided.
            - Applies a Sobel filter to the trace image along columns
              to detect slit edges using steep positive gradients
              (left edges) and steep negative gradients (right
              edges). See
              :func:`pypeit.core.trace.detect_slit_edges`.
            - Follows the detected left and right edges along
              spectrally adjacent pixels to identify coherent traces.
              See :func:`pypeit.core.trace.identify_traces`.
            - Performs basic handling of orphaned left or right
              edges. See
              :func:`pypeit.core.trace.handle_orphan_edges`.
            - Initializes the attributes that provide the trace
              position for each spectral pixel based on these
              results.

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are `filt_iter`,
        `sobel_mode`, `edge_thresh`, `minimum_spec_length`,
        `follow_span`, and `valid_flux_thresh`.

        The results of this are, by default, saved to the master
        frame; see `save` argument and :func:`save`.

        Args:
            img (`numpy.ndarray`_, :class:`pypeit.traceimage.TraceImage`):
                2D image used to trace slit edges. If a
                :class:`pypeit.traceimage.TraceImage` is provided,
                the raw files used to construct the image and on-chip
                binning are saved; the latter overrides any directly
                provided `binning`. The array should have shape
                :math:`(N_{\rm spec},N_{\rm spat})`; i.e., spectra
                are ordered along columns.
            bpm (`numpy.ndarray`_, optional):
                Bad-pixel mask for the trace image. Must have the
                same shape as `img`. If None, all pixels are assumed
                to be valid.
            det (:obj:`int`, optional):
                The 1-indexed detector number that provided the trace
                image. This is *only* used to determine whether or
                not bad columns in the image are actually along
                columns or along rows, as determined by
                :attr:`spectrograph` and the result of a call to
                :func:`pypeit.spectrograph.Spectrograph.raw_is_transposed`.
            binning (:obj:`str`, optional):
                On-detector binning of the data ordered spectral then
                spatial with format, e.g., `2,1`. Ignored if `img` is
                an instance of :class:`pypeit.traceimage.TraceImage`.
            save (:obj:`bool`, optional):
                Save the result to the master frame.
        """
        # TODO: Add debugging argument and hooks
        # Parse the input based on its type
        if isinstance(img, TraceImage):
            self.files = img.file_list
            _img = img.image
            self.binning = img.binning
            # TODO: does TraceImage have a mask?
        else:
            _img = img
            self.binning = binning

        # TODO: keep the TraceImage object instead of deconstructing
        # it?  For direct input, use a base PypeItImage object

        # Check the input
        if _img.ndim != 2:
            raise ValueError('Trace image must be 2D.')
        self.img = _img
        self.nspec, self.nspat = self.img.shape
        self.bpm = np.zeros((self.nspec, self.nspat), dtype=bool) if bpm is None else bpm
        if self.bpm.shape != self.img.shape:
            raise ValueError('Mask is not the same shape as the trace image.')
        self.det = det

        # Lightly smooth the image before using it to trace edges
        _img = ndimage.uniform_filter(self.img, size=(3, 1), mode='mirror')

        # Replace bad-pixel columns if they exist
        # TODO: Do this before passing the image to this function?
        # Instead of just replacing columns, replace all bad pixels...
        if np.any(self.bpm):
            # Do we need to replace bad *rows* instead of bad columns?
            flip = self.spectrograph.raw_is_transposed(det=self.det)
            axis = 1 if flip else 0

            # Replace bad columns that cover more than half the image
            bad_cols = np.sum(self.bpm, axis=axis) > (self.bpm.shape[axis]//2)
            if flip:
                # Deal with the transposes
                _img = procimg.replace_columns(_img.T, bad_cols, copy=True,
                                               replace_with='linear').T
            else:
                _img = procimg.replace_columns(_img, bad_cols, copy=True,
                                               replace_with='linear')

        # Filter the trace image and use the filtered image to detect
        # slit edges
        # TODO: Decide if mask should be passed to this or not,
        # currently not because of issues when masked pixels happen to
        # land in slit gaps.
        self.sobel_sig, edge_img \
                = trace.detect_slit_edges(_img, median_iterations=self.par['filt_iter'],
                                          sobel_mode=self.par['sobel_mode'],
                                          sigdetect=self.par['edge_thresh'])

        # Empty out the images prepared for left and right tracing
        # until they're needed.
        self.sobel_sig_left = None
        self.sobel_sig_right = None

        # Identify traces by following the detected edges in adjacent
        # spectral pixels.
        minimum_spec_pix = self.nspec * self.par['minimum_spec_length']
        trace_id_img = trace.identify_traces(edge_img, follow_span=self.par['follow_span'],
                                             minimum_spec_length=minimum_spec_pix)

        # Update the traces by handling single orphan edges and/or
        # traces without any left or right edges.
        flux_valid = np.median(_img) > self.par['valid_flux_thresh']
        trace_id_img = trace.handle_orphan_edge(trace_id_img, self.sobel_sig, bpm=self.bpm,
                                                flux_valid=flux_valid, copy=True)

        # Set the ID image to a MaskedArray to ease some subsequent
        # calculations; pixels without a detected edge are masked.
        trace_id_img = np.ma.MaskedArray(trace_id_img, mask=trace_id_img == 0)

        # Find the set of trace IDs; left traces are negative, right
        # traces are positive
        self.traceid = np.unique(trace_id_img.compressed())

        # Initialize the mask bits array for the trace coordinates and
        # just begin by setting them all as having no edge.
        self.spat_msk = self.bitmask.turn_on(np.zeros((self.nspec, self.ntrace),
                                                dtype=self.bitmask.minimum_dtype()), 'NOEDGE')

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

        # No design or object data yet
        self.design = None
        self.objects = None

        # Restart the log
        self.log = [inspect.stack()[0][3]]

        # Save if requested
        if save:
            self.save()

    def save(self, outfile=None, overwrite=True, checksum=True):
        """
        Save the trace object for later re-instantiation.

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
            msgs.error('Master file exists: {0}'.format(_outfile) + msgs.newline()
                       + 'Set overwrite=True to overwrite it.')

        # Report name before starting to write it
        msgs.info('Writing master frame to {0}'.format(_outfile))

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

        # If registered against slit-mask design data, include the
        # source file and fit parameters in the appropriate header
        designhdr = None
        if self.design is not None:
            designhdr = fits.Header()
            designhdr['MASKFILE'] = (self.design.meta['MASKFILE'],
                                     'File that provided slit-mask design data')
            designhdr['MASKOFF'] = (self.design.meta['MASKOFF'],
                                    'Best-fitting slit-mask design offset (pix)')
            designhdr['MASKSCL'] = (self.design.meta['MASKSCL'],
                                    'Best-fitting slit-mask design scale (pix per mm)')

        # Determine if the file should be compressed
        compress = False
        if _outfile.split('.')[-1] == 'gz':
            _outfile = _outfile[:_outfile.rfind('.')] 
            compress = True

        # Build the list of extensions to write
        # TODO: Separately adding the design and object data is
        # necessary because of a bug in the behavior of empty binary
        # tables in gzipped files in astropy.io.fits (astropy==3.1.2).
        # This has since been fixed, but I need to check if it works in
        # the most recent release. If it is, building the hdu can be
        # done in one go, where the binary tables will just be empty if
        # no design/object data is avialable.
        hdu = fits.HDUList([fits.PrimaryHDU(header=prihdr),
                            fits.ImageHDU(data=self.img, name='TRACEIMG'),
                            fits.ImageHDU(data=self.bpm.astype(np.int16), name='TRACEBPM'),
                            fits.ImageHDU(data=self.sobel_sig, name='SOBELSIG'),
                            fits.ImageHDU(data=self.traceid, name='TRACEID'),
                            fits.ImageHDU(data=self.spat_cen, name='CENTER'),
                            fits.ImageHDU(data=self.spat_err, name='CENTER_ERR'),
                            fits.ImageHDU(header=mskhdr, data=self.spat_msk, name='CENTER_MASK'),
                            fits.ImageHDU(header=fithdr, data=self.spat_fit, name='CENTER_FIT')
                          ])
        if self.design is not None: 
            hdu += [fits.BinTableHDU(header=designhdr, data=self.design, name='DESIGN')]
        if self.objects is not None: 
            hdu += [fits.BinTableHDU(data=self.objects, name='OBJECTS')]
        # Write the fits file; note not everything is written. Some
        # arrays are reconstruced by the load function.
        hdu.writeto(_outfile, overwrite=True, checksum=checksum)

        # Compress the file if the output filename has a '.gz'
        # extension; this is slow but still faster than if you have
        # astropy.io.fits do it directly
        if compress:
            msgs.info('Compressing file to: {0}.gz'.format(_outfile))
            io.compress_file(_outfile, overwrite=overwrite)

    def load(self, validate=True, rebuild_pca=True):
        """
        Load and reinitialize the trace data.

        Data is read from :attr:`file_path` and used to overwrite any
        internal data. Specific comparisons of the saved data are
        performed to ensure the file is consistent with having been
        written by a consistent version of the code; see
        :func:`_reinit`.

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
        self.img = hdu['TRACEIMG'].data
        self.nspec, self.nspat = self.img.shape
        self.bpm = hdu['TRACEBPM'].data.astype(bool)
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

        # Get the design and object data if they exist
        ext = [h.name for h in hdu]
        if 'DESIGN' in ext:
            self.design = table.Table(hdu['DESIGN'].data)
            self.design.meta['MASKFILE'] = hdu['DESIGN'].header['MASKFILE']
            self.design.meta['MASKOFF'] = hdu['DESIGN'].header['MASKOFF']
            self.design.meta['MASKSCL'] = hdu['DESIGN'].header['MASKSCL']
        if 'OBJECTS' in ext:
            self.objects = table.Table(hdu['OBJECTS'].data)

        # Set the integer pixel values
        self.spat_img = np.round(self.spat_cen if self.spat_fit is None
                                 else self.spat_fit).astype(int)

        # Rebuild the PCA if it existed previously and requested
        self.pca_type = None if hdu[0].header['PCATYPE'] == 'None' else hdu[0].header['PCATYPE']
        self._reset_pca(rebuild_pca and self.pca_type is not None)

        self.log = io.parse_hdr_key_group(hdu[0].header, prefix='LOG')

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
            raise ValueError('This master frame was generated using different parameter values!')

    def flagged_bits(self):
        """
        List the flags that have been tripped.

        Returns:
            list: List of the unique flags tripped by any trace
            datum.
        """
        return np.unique(np.concatenate([self.bitmask.flagged_bits(b) 
                                            for b in np.unique(self.spat_msk)])).tolist()

    def show(self, trace=None, include_error=False, thin=1, in_ginga=False, include_img=False,
             img_buffer=100, flag=None):
        """
        Show a scatter plot of the current trace data and fit, if
        it's available.

        Args:
            trace (:obj:`int`, array-like, optional):
                Trace ID number to show. Will show both positive and
                negative numbers. If None, all traces are plotted.
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
            img_buffer (:obj:`int`, optional):
                Buffer for plot window with respect to trace image
                size.
            flag (:obj:`str`, :obj:`list`, optional):
                One or more flags to use when *selecting* masked data
                to include in the plot; unmasked data are always
                plotted. If None, no masked data is plotted. If
                'any', trace locations masked for any reason are
                plotted. Otherwise, only trace locations flagged by
                the specified bits are plotted.
        """
        if in_ginga:
            # Currently can only show in ginga if the edges are
            # synchronized into slits
            if not self.is_synced:
                raise ValueError('To show in ginga, slit edges must be left-right synchronized.')
            ginga.connect_to_ginga(raise_err=True, allow_new=True)
            # Show Image
            viewer, ch = ginga.show_image(self.img, 'Trace Image')
            right = self.traceid > 0
            trace = self.spat_cen if self.spat_fit is None else self.spat_fit
            ginga.show_slits(viewer, ch, trace[:,np.invert(right)], trace[:,right],
                             self.traceid[right], pstep=thin)
            return

        # TODO: color/type-code the masked trace coordinates...

        # Build the data to plot. Any masked data is excluded; masked
        # data to be plotted are held in a separate array. This means
        # that errors and fits are currently never plotted for masked
        # data.
        _flag = None if flag in ['any', None] else np.atleast_1d(flag)
        cen = np.ma.MaskedArray(self.spat_cen, mask=self.bitmask.flagged(self.spat_msk))
        msk = None if flag is None \
                else np.ma.MaskedArray(self.spat_cen,
                                       mask=np.invert(self.bitmask.flagged(self.spat_msk,
                                                                           flag=_flag)))
        err = np.ma.MaskedArray(self.spat_err, mask=np.ma.getmaskarray(cen).copy())
        fit = None if self.spat_fit is None \
                else np.ma.MaskedArray(self.spat_fit, mask=np.ma.getmaskarray(cen).copy())

        # Show the traced image
        if include_img:
            plt.imshow(self.img.T, origin='lower', interpolation='nearest', aspect='auto')

        # Select the traces to show
        trace_indx = np.ones(self.ntrace, dtype=bool) if trace is None \
                        else np.isin(np.absolute(self.traceid), np.atleast_1d(trace))

        # Spectral position
        spec = np.tile(np.arange(self.nspec), (self.ntrace,1)).T

        if include_error:
            # Show the errors
            plt.errorbar(spec[::thin,trace_indx].ravel(), cen[::thin,trace_indx].ravel(),
                         yerr=err[::thin,trace_indx].ravel(), fmt='none', ecolor='k',
                         elinewidth=0.5, alpha=0.3, capthick=0, zorder=3)

        # Plot the left trace points
        left = (self.traceid < 0) & trace_indx
        plt.scatter(spec[::thin,left], cen[::thin,left], marker='.', color='k', s=30,
                    lw=0, zorder=4, label='left edge measurements', alpha=0.8)
        if msk is not None:
            plt.scatter(spec[::thin,left], msk[::thin,left], marker='x', color='k', s=20,
                        lw=0.5, zorder=5, label='masked left edges', alpha=0.8)

        # Plot the right trace points
        right = (self.traceid > 0) & trace_indx
        plt.scatter(spec[::thin,right], cen[::thin,right], marker='.', color='0.7',
                    s=30, lw=0, zorder=4, label='right edge measurements', alpha=0.8)
        if msk is not None:
            plt.scatter(spec[::thin,right], msk[::thin,right], marker='x', color='0.7', s=20,
                        lw=0.5, zorder=5, label='masked right edges', alpha=0.8)

        if self.spat_fit is None:
            # No fits, so we're done
            plt.legend()
            plt.show()
            return

        # Plot the trace fits
        for i in range(self.ntrace):
            # If statement structure primarily for the labels. Only
            # difference between left and right is the color.
            if left[i]:
                left_line = plt.plot(spec[::thin,i], self.spat_fit[::thin,i], color='C3', lw=0.5,
                                     zorder=6)
            elif right[i]:
                right_line = plt.plot(spec[::thin,i], self.spat_fit[::thin,i], color='C1', lw=0.5,
                                      zorder=6)

        # Limits and labels
        plt.xlim(-img_buffer, self.nspec+img_buffer)
        plt.ylim(-img_buffer, self.nspat+img_buffer)
        left_line[0].set_label('left edge fit')
        right_line[0].set_label('right edge fit')
        plt.xlabel('Spectral pixel index')
        plt.ylabel('Spatial pixel index')
        plt.legend()
        plt.show()

    def qa_plot(self, fileroot=None, min_spat=20):
        """
        Build a series of QA plots showing the edge traces.

        .. warning::
            This method is slow.

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
        img_zlim = utils.growth_lim(self.img, 0.95, fac=1.05)
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
                                                  flag=self.bitmask.bad_flags))
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
            ax.imshow(self.img.T, origin='lower', interpolation='nearest', vmin=img_zlim[0],
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

        This is primarily a wrapper for
        :func:`pypeit.core.trace.prepare_sobel_for_trace`; the boxcar
        smoothing is always 5 pixels. Barring a instantiation of the
        object, this calculation is only done once per side by "lazy
        loading" :attr:`sobel_sig_left` or :attr:`sobel_sig_right`.

        Args:
            side (:obj:`str`):
                The side to return; must be 'left' or 'right'.
    
        Returns:
            `numpy.ndarray`_: The manipulated Sobel image relevant to
            tracing the specified edge side.
        """
        # TODO: Add boxcar to EdgeTracePar?
        boxcar = 5
        if side == 'left':
            if self.sobel_sig_left is None:
                self.sobel_sig_left = trace.prepare_sobel_for_trace(self.sobel_sig, boxcar=boxcar,
                                                                    side='left')
            return self.sobel_sig_left
        if side == 'right':
            if self.sobel_sig_right is None:
                self.sobel_sig_right = trace.prepare_sobel_for_trace(self.sobel_sig, boxcar=boxcar,
                                                                     side='right')
            return self.sobel_sig_right
        raise ValueError('Side must be left or right.')

    def centroid_refine(self, follow=True, start_row=None, continuous=False, check=True,
                        use_fit=False):
        """
        Refine the edge positions using a moment analysis and assess
        the results.

        The method runs in two primary modes, depending on the
        value of `follow`:
            - If `follow=False`, the function simply executes
              :func:`pypeit.core.trace.masked_centroid` to recenter
              *all* the locations for the left and right edges at each
              spectral row independently. The maximum shift between
              any input and output centroid is `max_shift_abs` from
              :attr:`par`.
            - If `follow=True`, the method uses
              :func:`pypeit.core.trace.follow_centroid` to refine the
              centroids of the currently identified traces. Each edge
              is recentered starting from `start_row` and then moving
              to higher and then to lower spectral rows. In this
              case, the center of the aperture used for each row is
              the first moment (centroid) of the trace from returned
              for the preceding row. The maximum shift between the
              input and output center for the first row analyzed is
              set by `max_shift_abs`. For each subsequent row, the
              maximum shift is set by the `max_shift_adj` value in
              :attr:`par`.

        In either case, :func:`pypeit.core.trace.masked_centroid` is
        run with uniform weighting and an aperture width set to be
        twice `fwhm_uniform` from :attr:`par`. Other used parameters
        from :attr:`par` (:class:`pypeit.par.pypeitpar.EdgeTracePar`)
        are `max_shift_abs`, `max_shift_adj`, `max_spat_error`, and
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
            follow (:obj:`bool`, optional):
                Perform the centroiding first at a single row (see
                `start_row`) and then move to higher and lower
                spectral rows in series. See
                :func:`pypeit.core.trace.follow_centroid`.
            start_row (:obj:`int`, optional):
                The starting row when following the trace between
                adjacent rows; see `follow`. If None, the starting
                row is set by the finding the row that crosses the
                most unmasked trace positions; see
                :func:`pypeit.core.trace.most_common_trace_row`.
                Value is ignored if `follow=False`.
            continuous (:obj:`bool`, optional):
                Keep only the continuous part of the traces from the
                starting row.
            check (:obj:`bool`, optional):
                Run :func:`check_traces` to flag bad traces.
            use_fit (:obj:`bool`, optional):
                Use the fitted traces as the starting point for the
                refinement. Otherwise, uses :`spat_img`. If True and
                :attr:`spat_fit` is None, the method will raise an
                exception.
        """
        # Check input
        if use_fit and self.spat_fit is None:
            raise ValueError('No fit data available.')

        # Parse parameters and report
        width = 2 * self.par['fwhm_uniform']
        maxshift_start = self.par['max_shift_abs']
        maxshift_follow = self.par['max_shift_adj']
        maxerror = self.par['max_spat_error']

        msgs.info('Width of window for centroiding the edges: {0:.1f}'.format(width))
        msgs.info('Max shift between spectrally adjacent pixels: {0:.2f}'.format(maxshift_follow))
        msgs.info('Max centroid error: {0:.2f}'.format(maxerror))
    
        # To improve performance, generate bogus ivar and mask once
        # here so that they don't have to be generated multiple times.
        # TODO: Keep these as work space as class attributes?
        ivar = np.ones_like(self.sobel_sig, dtype=float)
        _bpm = np.zeros_like(self.sobel_sig, dtype=bool) if self.bpm is None else self.bpm
        fwgt = np.ones_like(self.sobel_sig, dtype=float)

        # Book-keeping objects to keep track of which traces to remove
        rmtrace = np.zeros(self.ntrace, dtype=bool)

        # To hold the refined traces and mask
        spat = self.spat_fit if use_fit else self.spat_img
        cen = np.zeros_like(self.spat_cen)
        err = np.zeros_like(self.spat_err)
        msk = np.zeros_like(self.spat_msk, dtype=self.bitmask.minimum_dtype())

        # Refine left then right
        for side in ['left', 'right']:
            
            # Get the image relevant to tracing this side
            _sobel_sig = self._side_dependent_sobel(side)

            # Identify the traces on the correct side: Traces on the
            # left side are negative.
            indx = self.traceid < 0 if side == 'left' else self.traceid > 0
            msgs.info('Number to retrace: {0}'.format(np.sum(indx)))

            if follow:
                # Find the most common row index
                bpm = self.bitmask.flagged(self.spat_msk, flag=self.bitmask.bad_flags)
                _start_row = trace.most_common_trace_row(bpm) if start_row is None else start_row
                msgs.info('Starting row is: {0}'.format(_start_row))

                # Trace starting from this row
                msgs.info('Sequentially tracing first moment of Sobel-filtered image to higher '
                        'and lower spectral positions.')
                cen[:,indx], err[:,indx], msk[:,indx] \
                        = trace.follow_centroid(_sobel_sig, _start_row, spat[_start_row,indx],
                                                ivar=ivar, bpm=_bpm, fwgt=fwgt, width=width,
                                                maxshift_start=maxshift_start,
                                                maxshift_follow=maxshift_follow, maxerror=maxerror,
                                                continuous=continuous, bitmask=self.bitmask)
            else:
                cen[:,indx], err[:,indx], msk[:,indx] \
                        = trace.masked_centroid(_sobel_sig, spat[:,indx], width, ivar=ivar,
                                                bpm=_bpm, fwgt=fwgt, maxshift=maxshift_start,
                                                maxerror=maxerror, bitmask=self.bitmask,
                                                fill='bound')

            # Check the traces
            mincol = None if side == 'left' else 0
            maxcol = _sobel_sig.shape[1]-1 if side == 'left' else None
            good, bad = self.check_traces(cen, err, msk, subset=indx, mincol=mincol, maxcol=maxcol)

            # Save the results and update the book-keeping
            self.spat_cen[:,good] = cen[:,good]
            self.spat_err[:,good] = err[:,good]
            self.spat_msk[:,good | bad] |= msk[:,good | bad]
            rmtrace[bad] = True

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
              tolerance (`match_tol` in :attr:`par`).
            - Remove traces that do not cover at least some fraction of
              the detector (`minimum_spec_length` in :attr:`par`).
            - Remove traces that are at a minimum or maximum column
              (typically the edge of the detector).

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are `match_tol`
        and `minimum_spec_length`.

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
        # Indices of traces to check
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
                spat_bpm = self.bitmask.flagged(self.spat_msk, flag=self.bitmask.bad_flags)
                spat_img = np.ma.MaskedArray(self.spat_img, mask=spat_bpm)
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
        col = np.round(cen).astype(int)
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

    @property
    def is_synced(self):
        """
        Confirm slit edges are synced.
        """
        side = np.clip(self.traceid, -1, 1)
        return side.size % 2 == 0 and np.all(side[1:] + side[:-1] == 0)

    def check_synced(self, rebuild_pca=False):
        """
        Quality check and masking of the synchronized edges.

        Before executing this method, the slit edges must be
        synchronized (see :func:`sync`) and ordered spatially
        in left-right pairs (see :func:`spatial_sort`). The former is
        checked explicitly.

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are
        `minimum_slit_length` and `length_range`.

        Checks are:
            - Any trace falling off the edge of the detector is
              masked (see :class:`EdgeTraceBitMask`). This is the
              only check performed by default (i.e., when no keyword
              arguments are provided).
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
                If the pca exists and the masking or number of traces
                changes because of the performed checks, redo the PCA
                decomposition using the new traces and trace masks
                and the previous parameter set.
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
        if not self.is_synced:
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

        # TODO: Check that slit edges meet a minimum slit gap?

        if self.par['sync_clip']:
            # Remove traces that have been fully flagged as bad
            rmtrace = np.all(self.bitmask.flagged(self.spat_msk, flag=self.bitmask.bad_flags),
                             axis=0)
            if np.sum(rmtrace) == self.ntrace:
                msgs.error('All slit edges are fully masked!')
            self.remove_traces(rmtrace, rebuild_pca=rebuild_pca, sync_rm='both')
        elif new_masks:
            # Reset the PCA if new masks are applied
            self._reset_pca(rebuild_pca and self.pca is not None)

    def remove_traces(self, indx, resort=True, rebuild_pca=False, sync_rm='ignore'):
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
            sync_rm (:obj:`str`, optional):
                If the traces are left-right synchronized (see
                :func:`sync` and :func:`is_synced`), use this method
                to deal with edges paired with those to be removed.
                Methods are:
                    - 'ignore': Just remove the selected traces and
                      ignore the synchronization.
                    - 'both': If at least one of the traces in a pair
                      is selected for removal, remove both.
                    - 'neither': If only one of the traces in a pair
                      is selected for removal, remove neither.
        """
        # Make sure there are traces to remove
        if not np.any(indx):
            msgs.warn('No traces removed.')
            return

        # Deal with removing traces when they are left-right
        # synchronized
        _indx = np.atleast_1d(indx).copy()
        if self.is_synced and sync_rm != 'ignore':
            if sync_rm == 'both':
                _indx = np.repeat(np.any(_indx.reshape(-1,2), axis=1), 2)
                msgs.info('Removing both traces if either selected in synchonized pair.')
            elif sync_rm == 'neither':
                _indx = np.repeat(np.all(_indx.reshape(-1,2), axis=1), 2)
                msgs.info('Removing both traces only if both selected in synchonized pair.')
            else:
                raise ValueError('Unknown sync removal keyword: {0}'.format(sync_rm))
            
        msgs.info('Removing {0} edge traces.'.format(np.sum(_indx)))

        # Reset the trace data
        keep = np.invert(_indx)
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

    def spatial_sort(self, use_mean=False, use_fit=True):
        """
        Sort the traces spatially.

        The coordinates used for the sorting is either the measured
        centroids or the fitted parameterization (see `use_fit`). The
        fiducial coordinates that are sorted are either the mean of
        the unmasked coordinates over all spectral rows or the
        unmasked coordinates at a specified reference row (see
        `use_mean`).

        The trace IDs are also reassigned to be sorted spatially;
        i.e., the trace IDs for three synced slits would be `[-1, 1,
        -2, 2, -3, 3]`.
        
        All attributes are edited in-place.

        Args:
            use_mean (:obj:`bool`, optional):
                Sort according to the mean of the masked spatial
                positions. If False, the spatial position at a
                reference row is used, where the reference row is
                either the same as used by the PCA (if available) or
                the result of
                :func:`pypeit.core.trace.most_common_trace_row` using
                the current trace mask.
            use_fit (:obj:`bool`, optional):
                Sort according to the fit positions instead of the
                measured positions. Otherwise, only use the fit
                positions if they're available and the measured
                location is masked.
        """
        # Check input
        if use_fit and self.spat_fit is None:
            msgs.warn('Fit data is not available; cannot use it for spatially sorting the edges.')

        # Set up the coordinates to use
        bpm = self.bitmask.flagged(self.spat_msk, self.bitmask.bad_flags)
        cen = self.spat_fit.copy() if self.spat_fit is not None and use_fit \
                    else self.spat_cen.copy()
        # Replace masked values with the fit if it is available
        if self.spat_fit is not None and not use_fit:
                cen[bpm] = self.spat_fit[bpm]

        # Get the sorted indices
        if use_mean:
            # Sort the traces by their spatial position (always use
            # measured positions even if fit positions are available)
            srt = np.argsort(np.mean(cen, axis=0))
        else:
            # Sort according to the spatial position in one row
            reference_row = trace.most_common_trace_row(bpm) if self.pca is None \
                                else self.pca.reference_row
            msgs.info('Re-sorting edges based on where they cross row {0}'.format(reference_row))
            srt = np.argsort(cen[reference_row,:])

        # Resort the arrays
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
        the relevant attributes to `None` when `rebuild` is False).
        """
        if rebuild:
            # Rebuild the PCA using the previous parameters
            return self.build_pca(use_center=self.pca_type == 'center') 
        # Remove the existing PCA
        self.pca_type = None
        self.pca = None

    def current_trace_locations(self):
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
        Iteratively re-measure and fit a functional form to the edge
        location data.

        Primarily a wrapper for :func:`pypeit.core.trace.fit_trace`,
        run once per edge side (left and right).

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are
        `fit_function`, `fit_order`, `fwhm_uniform`, `fwhm_gaussian`,
        `fit_maxdev`, `fit_maxiter`, and `fit_niter`.

        Args:
            weighting (:obj:`str`, optional):
                The weighting to apply to the position within each
                integration window (see
                :func:`pypeit.core.trace.fit_trace`).
            debug (:obj:`bool`, optional):
                Run in debug mode.
            idx (`numpy.ndarray`_, optional):
                Array of strings with the IDs for each object. Used
                only if `debug` is true for the plotting. Default is
                just a running number.
        """
        # Parse parameters and report
        maxshift = self.par['max_shift_abs']
        maxerror = self.par['max_spat_error']
        function = self.par['fit_function']
        order = self.par['fit_order']
        fwhm = self.par['fwhm_uniform'] if weighting == 'uniform' else self.par['fwhm_gaussian']
        maxdev = self.par['fit_maxdev']
        maxiter = self.par['fit_maxiter']
        niter = self.par['fit_niter']
        xmin = 0.
        xmax = self.nspec-1.

        msgs.info('Max shift btwn input and remeasured edge centroids: {0:.2f}'.format(maxshift))
        msgs.info('Max centroid error: {0:.2f}'.format(maxerror))
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
        bpm = np.zeros_like(self.sobel_sig, dtype=bool) if self.bpm is None else self.bpm

        # Initialize arrays
        fit = np.zeros_like(self.spat_cen, dtype=float)
        cen = np.zeros_like(self.spat_cen, dtype=float)
        err = np.zeros_like(self.spat_cen, dtype=float)
        msk = np.zeros_like(self.spat_cen, dtype=self.bitmask.minimum_dtype())

        # Flag bad traces; excludes inserted traces
        spat_bpm = self.bitmask.flagged(self.spat_msk, flag=self.bitmask.bad_flags)

        # Fit both sides
        for side in ['left', 'right']:
            # Get the image relevant to tracing this side
            _sobel_sig = self._side_dependent_sobel(side)
            # Select traces on this side
            this_side = self.traceid < 0 if side == 'left' else self.traceid > 0
            # Perform the fit
            fit[:,this_side], cen[:,this_side], err[:,this_side], msk[:,this_side] \
                    = trace.fit_trace(_sobel_sig, self.spat_cen[:,this_side], order, ivar=ivar,
                                      bpm=bpm, trace_bpm=spat_bpm[:,this_side],
                                      weighting=weighting, fwhm=fwhm, maxshift=maxshift,
                                      maxerror=maxerror, function=function, maxdev=maxdev,
                                      maxiter=maxiter, niter=niter, bitmask=self.bitmask,
                                      debug=debug, idx=idx, xmin=xmin, xmax=xmax)

        # Save the results of the edge measurements
        self.spat_cen = cen
        self.spat_err = err
        # Add the new masking information; this merges the new and
        # existing mask because the existing mask was used by fit_trace
        # to ignore input trace data
        # TODO: Does applying this mask here mess up fits later?
        self.spat_msk |= msk
        # Save the model fits
        self.spat_fit = fit
        self.spat_fit_type = '{0} : order={1}'.format(function, order)
        # Set the pixelated trace data based on the fit, instead of the
        # measured centroids
        self.spat_img = np.round(self.spat_fit).astype(int)
        # Append function execution to log
        self.log += [inspect.stack()[0][3]]

    def build_pca(self, use_center=False, debug=False):
        """
        Build a PCA model of the current edge data.

        Primarily a wrapper that instantiates :attr:`pca`, which has
        type :class:`EdgeTracePCA`. After executing this, traces can
        be predicted using the pca by calling
        `self.pca.predict(spat)`; see :func:`EdgeTracePCA.predict`.

        If no parametrized function has been fit to the trace data or
        if specifically requested (see `use_center`), the PCA is
        based on the measured trace centroids; othwerwise, the PCA
        uses the parametrized trace fits.

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are `pca_n`,
        `pca_var_percent`, `pca_function`, `pca_order`, `pca_sigrej`,
        `pca_maxrej`, and `pca_maxiter`.

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

        # When constructing the PCA, ignore bad trace measurements
        # *and* any traces inserted by hand.
        bpm = self.bitmask.flagged(self.spat_msk)

        # Only use traces that are unmasked for at least some fraction of the detector
        # TODO: Make this minimum length check a paramter
        minimum_pca_length=0.9
        # TODO: Is there a way to propagate the mask to the PCA?
        # TODO: Keep a separate mask specifically for the fit data?
        use_trace = np.sum(np.invert(bpm), axis=0)/self.nspec > minimum_pca_length
        if np.sum(use_trace) < 2:
            msgs.error('Insufficient traces for PCA decomposition.')
        trace_inp = self.spat_cen[:,use_trace] if self.spat_fit is None or use_center \
                        else self.spat_fit[:,use_trace]
        msgs.info('Using {0}/{1} traces in the PCA analysis (omitting short traces).'.format(
                  np.sum(use_trace), self.ntrace))

        # Instantiate the PCA
        self.pca = EdgeTracePCA(trace_inp, npca=npca, pca_explained_var=pca_explained_var,
                                reference_row=trace.most_common_trace_row(bpm[:,use_trace]))

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

        # TODO: Compare with the fit data. Remove traces where the mean
        # offset between the PCA prediction and the measured centroids
        # are larger than some threshold

        # Log what was done
        self.log += [inspect.stack()[0][3]]

    def peak_refine(self, rebuild_pca=False, debug=False):
        """
        Refine the trace by isolating peaks and troughs in the
        Sobel-filtered image.

        This function *requires* that the PCA model exists; see
        :func:`build_pca` or :func:`pca_refine`. It is also primarily
        a wrapper for :func:`pypeit.core.trace.peak_trace`.

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are
        `edge_thresh`, `smash_range`, `edge_detect_clip`,
        trace_median_frac`, `trace_thresh`, `fit_function`,
        `fit_order`, `fwhm_uniform`, `fwhm_uniform`,
        `niter_gaussian`, `niter_gaussian`, `fit_maxdev`, and
        `fit_maxiter`.

        Args:
            rebuild_pca (:obj:`bool`, optional):
                This method fundamentally resets the trace data,
                meaning that the PCA is no longer valid. Use this
                boolean to have the method rebuild the PCA based on
                the refined traces. Note that the PCA is *not* then
                used to reset the fitted trace data; i.e.,
                :attr:`spat_fit` remains based on the output of
                :func:`pypeit.core.trace.peak_trace`.
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
        msgs.info('Image fraction for trace mask filter: {0}'.format(trace_median_frac))
        msgs.info('Threshold for trace masking: {0}'.format(trace_thresh))
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
        bpm = np.zeros_like(self.sobel_sig, dtype=bool) if self.bpm is None else self.bpm

        # Get the image relevant to tracing
        _sobel_sig = trace.prepare_sobel_for_trace(self.sobel_sig, boxcar=5, side=None)

        # TODO: Need to fix how the rest of this handles bad_trace

        # Find and trace both peaks and troughs in the image. The input
        # trace data (`trace` argument) is the PCA prediction of the
        # trace that passes through each spatial position at the
        # reference spectral pixel.
        fit, cen, err, msk, nleft \
                = trace.peak_trace(_sobel_sig, ivar=ivar, bpm=bpm,
                                   trace_map=self.pca.predict(np.arange(self.nspat)),
                                   smash_range=smash_range, peak_thresh=peak_thresh,
                                   peak_clip=peak_clip, trough=True,
                                   trace_median_frac=trace_median_frac, trace_thresh=trace_thresh,
                                   fwhm_uniform=fwhm_uniform, fwhm_gaussian=fwhm_gaussian,
                                   function=function, order=order, maxdev=maxdev, maxiter=maxiter,
                                   niter_uniform=niter_uniform, niter_gaussian=niter_gaussian,
                                   bitmask=self.bitmask, debug=debug)

        # Assess the output
        ntrace = fit.shape[1]
        if ntrace < self.ntrace:
            msgs.warn('Found fewer traces using peak finding than originally available.  '
                      'May want to reset peak threshold.')

#        # TODO: Check the traces? NOTE: This currently only identifies short
#        # traces.
#        good, bad = self.check_traces(trace_cen, trace_err, bad_trace)
#        left = np.zeros(trace_fit.shape[1])

        # Reset the trace data
        self.traceid = np.zeros(ntrace, dtype=int)
        self.traceid[:nleft] = -1-np.arange(nleft)
        self.traceid[nleft:] = 1+np.arange(ntrace-nleft)
        self.spat_fit = fit
        self.spat_fit_type = '{0} : order={1}'.format(function, order)
        self.spat_cen = cen
        self.spat_err = err
        # The mask is entirely new and shouldn't be merged with the
        # existing mask.
        self.spat_msk = msk
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

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are
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

        Returns:
            `numpy.ndarray`_: Reference positions for all edge
            traces, both for the existing and new traces.
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
        bpm = self.bitmask.flagged(self.spat_msk, flag=self.bitmask.bad_flags)
        reference_row = trace.most_common_trace_row(bpm) if self.pca is None \
                            else self.pca.reference_row

        # Check that the trace data are sorted at this row
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
        # NOTE: Slit length should already be masked; np.ma.diff
        # returns a masked value when one of the slit edges is masked.
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

        # TODO: Nothing should now be masked. Get rid of this once
        # satisfied that the coding is correct.
        if np.any(trace_ref.mask):
            raise ValueError('Coding error: this should not happen')
        trace_ref = trace_ref.data.ravel()

        # Check that the predicted reference positions don't cause slit
        # overlaps
        indx = np.where(add_edge[1:-1])[0]
        if len(indx) > 0:
            indx += 1
            msgs.warn('Found {0} overlapping slit edges.  '.format(len(indx))
                      + 'Offsetting to minimum gap of {0} pixel(s).'.format(min_slit_gap))
            # Predicted below an existing edge (will never add paired edges)
            too_lo = trace_ref[indx] < trace_ref[indx-1]
            trace_ref[indx[too_lo]] = trace_ref[indx[too_lo]-1]+min_slit_gap
            # Predicted above an existing edge
            too_hi = trace_ref[indx] > trace_ref[indx+1]
            trace_ref[indx[too_hi]] = trace_ref[indx[too_hi]+1]-min_slit_gap

        return trace_ref

    def nudge_traces(self, trace):
        r"""
        Nudge traces away from the detector edge.

        Traces are shifted spatially, up to a maximum value set by
        `max_nudge`, to be no closer than a minimum of `det_buffer`
        pixels from the detector edges. Both parameters are pulled
        from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`). No limit is
        imposed on the size of the shift if `max_nudge` is None.

        .. warning::
            A failure mode that is not dealt with is when multiple
            traces fall off the detector and are nudged to nearly the
            same location.

        Args:
            trace (`numpy.ndarray`_):
                Array with trace locations to adjust. Must be 2D with
                shape :math:`(N_{\rm spec}, N_{\rm trace})`.

        Returns:
            `numpy.ndarray`_: The nudged traces.
        """
        # Check input
        if self.par['max_nudge'] is not None and self.par['max_nudge'] <= 0:
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

        # NOTE: Should never happen, but this makes a compromise if a
        # trace crosses both the left and right spatial edge of the
        # detector.
        offset = np.clip(_buffer - np.amin(trace, axis=0), 0, self.par['max_nudge']) \
                    + np.clip(self.nspat - 1 - _buffer - np.amax(trace, axis=0),
                              None if self.par['max_nudge'] is None else -self.par['max_nudge'], 0)
        # Offset and return the traces
        return trace + offset[None,:]

    def sync(self, rebuild_pca=True, debug=False):
        """
        Match left and right edge traces to construct slit edge pairs.

        First, the method ensures that the edge traces are sorted
        spatiall; see :func:`spatial_sort`. Then it determines where
        traces need to be inserted to create a left-right pair (using
        :func:`_get_insert_locations`).

        The next steps are determine the reference positions where
        the traces should be inserted and the shape the trace should
        take. The former is determined by the `sync_center` parameter
        in :attr:`par`, and the latter is determined by `sync_predict`.

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

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are `det_buffer`
        and `sync_predict`.

        Args:
            rebuild_pca (:obj:`bool`, optional):
                If the pca exists and traces are removed (see
                :func:`check_synced`), rebuild the PCA using the new
                traces and the previous parameter set. Note that
                inserted traces are *not* included in the PCA
                decomposition.
            debug (:obj:`bool`, optional):
                Run in debug mode.
        """
        # Check input
        if self.par['sync_predict'] not in ['pca', 'nearest']:
            raise ValueError('Unknown trace mode: {0}'.format(self.par['sync_predict']))
        if self.par['sync_predict'] == 'pca' and self.pca is None:
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
            if self.par['sync_predict'] == 'pca':
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
        if self.par['sync_predict'] == 'pca':
            trace_add = self.pca.predict(trace_ref[add_edge])
        elif self.par['sync_predict'] == 'nearest':
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

        # Insert the new traces
        self.insert_traces(side[add_edge], trace_add, loc=add_indx[add_edge], mode='sync')

        if debug:
            msgs.info('Show instance includes inserted traces but before checking the sync.')
            self.show(thin=10, include_img=True)

        # Check the full synchronized list and log completion of the
        # method
        self.check_synced(rebuild_pca=rebuild_pca)
        self.log += [inspect.stack()[0][3]]

    def insert_traces(self, side, trace, loc=None, mode='user', resort=True):
        r"""
        Insert/append a set of edge traces.

        New traces to add are first nudged away from the detector
        edge (see :func:`nudge_traces`) according to parameters
        `max_nudge` and `det_buffer` from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`). They are then
        inserted or appended to the existing traces and masked
        according to the provided `mode`. The traces are added to
        *both* the measured centroid list and the fitted model data.
        Then the full list of traces can be resorted spatially,
        according to the provided `resort`.

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

        msgs.info('Inserting {0} new traces.'.format(ntrace))

        # Nudge the traces
        _trace = self.nudge_traces(_trace)

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

    def fully_masked_traces(self, flag=None):
        """
        Return a boolean array selecting traces that are fully masked
        by the provided flag(s).
        
        Args:
            flag (:obj:`str`, :obj:`list`, optional):
                The bit mask flags to select. If None, any flags are
                used. See :func:`pypeit.bitmask.Bitmask.flagged`.

        Returns:
            `numpy.ndarray`_: Boolean array selecting traces that are
            flagged at all spectral pixels.
        """
        return np.all(self.bitmask.flagged(self.spat_msk, flag=flag), axis=0)
    
    def mask_refine(self, design_file=None, allow_resync=False, debug=False):
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

        TODO: Traces don't need to be synchronized...

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

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are
        `mask_reg_maxiter`, `mask_reg_maxsep`, `mask_reg_sigrej`, and
        `ignore_alignment`.

        Args:
            design_file (:obj:`str`, optional):
                A file with the mask design data. If None, the method
                will use the first file in :attr:`files`; if
                :attr:`files` is also None, the method will raise an
                exception.
            debug (:obj:`bool`, optional):
                Run in debug mode.
        """
        # Still not done with this function...
        raise NotImplementedError()

        # The PCA decomposition must have already been determined
        if self.pca is None:
            msgs.error('Must first run the PCA analysis for the traces; run build_pca.')

        # Get the file to use when parsing the mask design information
        _design_file = (None if self.files is None else self.files[0]) if design_file is None \
                            else design_file
        if _design_file is None or not os.path.isfile(_design_file):
            msgs.error('Slit-mask design file not found or none provided.')

        # Get the paramters to use
        maxiter = self.par['mask_reg_maxiter']
        maxsep = self.par['mask_reg_maxsep']
        sigma = self.par['mask_reg_sigrej']
        ignore_alignment = self.par['ignore_alignment']

        # TODO: Set allow_resync and design_file to be a parameters, as
        # well?

        # Read the design data
        msgs.info('Reading slit-mask design information from: {0}'.format(_design_file))
        if self.spectrograph.get_slitmask(_design_file) is None:
            msgs.error('Unable to read design file or no slit-mask design reader '
                       'defined for {0}.'.format(self.spectrograph.spectrograph))

        # Match both left and right edges simultaneously
        x_design = np.array([self.spectrograph.slitmask.bottom[:,0],
                             self.spectrograph.slitmask.top[:,0]]).T.ravel()
        x_det = self.spat_fit[self.pca.reference_row,:]

        # Mask traces that are fully masked, except if they were
        # specifically inserted in a previous step
        x_det_bpm = self.fully_masked_traces(flag=self.bitmask.bad_flags) \
                        & np.invert(self.fully_masked_traces(flag=self.bitmask.insert_flags))

#        x_design = np.amin(self.spectrograph.slitmask.corners[:,:,0], axis=1)
#        side = self.traceid < 0
#        x_det = self.spat_fit[self.pca.reference_row,side]

#        x_design = np.amax(self.spectrograph.slitmask.corners[:,:,0], axis=1)
#        side = self.traceid > 0
#        x_det = self.spat_fit[self.pca.reference_row,side]

        # Estimate the scale in pixels/mm as the telescope platescale
        # in arcsec/mm divided by the detector platescale in
        # arcsec/pixel
        pix_per_mm = self.spectrograph.telescope.platescale() \
                        / self.spectrograph.detector[self.det-1]['platescale']

        # If the traces are synchronized, use the estimated scale to
        # first mask edges that yeild slits that are too small relative
        # to the range of slit lengths in the mask file.
        if self.is_synced:
            slit_len_det = np.diff(x_det.reshape(-1,2), axis=1).ravel()
            slit_len_mask = np.diff(x_design.reshape(-1,2), axis=1).ravel()*pix_per_mm
            indx = (slit_len_det < np.amin(slit_len_mask)/1.1) \
                        | (slit_len_det > np.amax(slit_len_mask)*1.1)
            if np.any(indx):
                msgs.info('Removing {0} edges that form (an) '.format(np.sum(indx)*2)
                          + 'errantly small or large slit(s) compared to the mask design data.')
                x_det_bpm[np.repeat(indx,2)] = True

        # Initial guess for the offset
        try:
            # Try using the spectrograph detector map
            self.spectrograph.get_detector_map()
            # Set the offset based on the location of this detector
            offset = self.spectrograph.detector_map.image_coordinates(
                            self.spectrograph.detector_map.npix[0]/2,
                            self.spectrograph.detector_map.npix[1]/2, detector=self.det,
                            in_mm=False)[0][0] - self.spectrograph.detector_map.npix[0]/2
            # Set the bounds to some nominal fraction of the detector
            # size and pix/mm scale; allow for a +/- 10% deviation in
            # the pixel scale
            # TODO: Is 10% generally enough (for any instrument)? Make
            # this a (spectograph-specific) parameter?
            offset_rng = [offset-0.1*self.spectrograph.detector_map.npix[0],
                          offset+0.1*self.spectrograph.detector_map.npix[0]]
        except:
            # No detector map
            msgs.warn('No detector map available for {0}'.format(self.spectrograph.spectrograph)
                      + '; attempting to match to slit-mask design anyway.')
            # Set the guess offset such that two sets of coordinates
            # are offset to their mean
            offset = np.mean(x_det) - np.mean(pix_per_mm * x_design)
            # Set the offset range
            offset_rng = [offset-np.absolute(np.amin(x_det)-np.amin(pix_per_mm*x_design))*1.1,
                          offset+np.absolute(np.amax(pix_per_mm*x_design)-np.amax(x_det))*1.1]

#        import pdb
#        pdb.set_trace()
#
#        slitmask.xc_trace(x_det, x_design, pix_per_mm)
#
#        pdb.set_trace()

        # The solution can be highly dependent on the initial guess for
        # the offset, so do an initial grid search to get close to the
        # solution.
        msgs.info('Running a grid search to try to find the best starting offset.')
        # Step by 2 pixels
        off = np.arange(offset_rng[0], offset_rng[1], 2).astype(float)
        rms = np.zeros_like(off, dtype=float)
        scl = np.zeros_like(off, dtype=float)
        par = np.array([0, pix_per_mm])
        bounds = np.array([offset_rng, [pix_per_mm/1.1, pix_per_mm*1.1]])
        register = slitmask.SlitRegister(x_det, x_design, trace_mask=x_det_bpm)

        # NOTE: The commented approach below gets the RMS at each
        # offset point just using the estimated scale. This is faster
        # than the approach taken, but results are sensitive to the
        # accuracy of the estimated scale, which can lead to problems
        # in corner cases.
#        for i in range(off.size):
#            print('Grid point: {0}/{1}'.format(i+1, off.size), end='\r')
#            par[0] = off[i]
#            register.par = par
#            minsep = register.match(unique=True)[1]
#            rms[i] = sigma_clipped_stats(minsep, sigma=5)[2]
#        print('Grid point: {0}/{0}'.format(off.size))

        # For each grid point, keep the offset fixed and find the best
        # scale. No rejection iterations are performed.
        for i in range(off.size):
            print('Grid point: {0}/{1}'.format(i+1, off.size), end='\r')
            par[0] = off[i]
            register.find_best_match(guess=par, fix=[True,False], bounds=bounds, penalty=False)
            minsep = register.match(unique=True)[1]
            scl[i] = register.par[1]
            rms[i] = sigma_clipped_stats(minsep, sigma=5)[2]
        print('Grid point: {0}/{0}'.format(off.size))

        # Use the grid point with the best RMS
        minindx = np.argmin(rms)
        offset = off[minindx]
        best_rms = rms[minindx]
        msgs.info('Minimum RMS ({0:.2f}) found with offset = {1:.2f}'.format(best_rms, offset))
        if debug:
            # Plot the result
            ax1 = plt.subplot(211)
            ax1.scatter(off, rms, color='k', marker='.', s=100, lw=0, zorder=0)
            ax1.scatter(offset, best_rms, color='C3', marker='x', s=50, zorder=1)
            ax1.set_xlabel('Trace Offset (pix)')
            ax1.set_ylabel('RMS (det-mask; pix)')
            ax1.set_title('Grid search for initial offset')
            ax2 = plt.subplot(212, sharex=ax1)
            ax2.scatter(off, scl, color='k', marker='.', s=100, lw=0, zorder=0)
            ax2.set_ylabel('Best-fit scale')
            plt.show()

        # Do the final fit with some rejection iterations
        register.find_best_match(guess=[offset, pix_per_mm], bounds=bounds, penalty=False,
                                 maxiter=maxiter, maxsep=maxsep, sigma=sigma, debug=debug)

        if debug:
            register.show(minmax=[0, self.nspat], synced=True)

        # Find the missing, bad, and masked traces
        missing, bad = register.trace_mismatch(minmax=[0, self.nspat], synced=True)
#        masked_by_registration = np.where(register.trace_mask & np.invert(x_det_bpm))[0]
#        bad = np.append(bad, masked_by_registration)
        bad = np.append(bad, np.where(register.trace_mask | x_det_bpm)[0])

        # Ignore missing alignment boxes
        if ignore_alignment:
            missing = missing[np.invert(self.spectrograph.slitmask.alignment_slit[missing//2])]
            found_alignment_slits = register.match_index[
                            self.spectrograph.slitmask.alignment_slit[register.match_index//2]]
            bad = np.append(bad, found_alignment_slits)

        # Report
        msgs.info('Best-fitting offset and scale for mask coordinates: {0:.2f} {1:.2f}'.format(
                    *register.par))
        msgs.info('Traces will {0} alignment slits'.format('exclude' if ignore_alignment
                                                             else 'include'))
        msgs.info('Number of missing mask traces to insert: {0}'.format(len(missing)))
        msgs.info('Number of bad or alignment traces to remove: {0}'.format(len(bad)))

        if self.is_synced and (len(missing) - len(bad)) % 2 != 0:
            if allow_resync:
                msgs.warning('Difference in added and removed traces is odd; will resync traces.')
            else:
                msgs.error('Difference in added and removed traces desyncronizes traces.')

        if len(bad) > 0:
            # Remove the bad traces and rebuild the pca
            rmtrace = np.zeros(self.ntrace, dtype=bool)
            rmtrace[bad] = True
            self.remove_traces(rmtrace, rebuild_pca=True)

        if len(missing) > 0:
            # Even indices are lefts, odd indices are rights
            side = missing % 2 * 2 - 1
            # Predict the traces using the PCA
            missing_traces = self.pca.predict(register.match_coo[missing])
            # Insert them
            self.insert_traces(side, missing_traces, mode='mask')

#        import pdb
#        pdb.set_trace()

        if len(bad) > 0 or len(missing) > 0:
            # Traces were removed and/or inserted, resync or recheck that the edges are synced.
            if (len(missing) - len(bad)) % 2 != 0 and allow_resync:
                self.sync(rebuild_pca=True)
            else:
                self.check_synced(rebuild_pca=True)

            # Reset the match after removing/inserting traces
            x_det = self.spat_fit[self.pca.reference_row,:]
            x_det_bpm = self.fully_masked_traces(flag=self.bitmask.bad_flags) \
                            & np.invert(self.fully_masked_traces(flag=self.bitmask.insert_flags))
            register = slitmask.SlitRegister(x_det, x_design, trace_mask=x_det_bpm,
                                             guess=[offset, pix_per_mm], bounds=bounds,
                                             penalty=False, maxiter=maxiter, maxsep=maxsep,
                                             sigma=sigma, debug=debug, fit=True)

            # TODO: This fit should *never* result in missing or bad
            # traces! Keep this for a while until we feel like we've
            # vetted the code well enough.
            missing, bad = register.trace_mismatch(minmax=[0, self.nspat], synced=True)
            if len(missing) != 0 or len(bad) != 0:
                 msgs.error('CODING ERROR: Should never find missing or bad traces in re-fit!')

        # Fill the slit-design and object tables
        self._fill_design_table(register, _design_file)
        self._fill_objects_table(register)

    def _fill_design_table(self, register, design_file):
        """
        Fill :attr:`design` based on the results of the design
        registration.

        Args:
            register (:class:`pypeit.spectrographs.slitmask.SlitRegister`):
                Object holding the result of the registration.
            design_file (:obj:`str`):
                File that provided the slit-mask design data.
        """
        # Index for the slit in the design data
        slit_index = register.match_index[register.match_index % 2 == 0]//2
        # Number of slits
        nslits = len(slit_index)
        # Instantiate as an empty table
        self.design = EdgeTraceSet.empty_design_table(rows=nslits)
        # Save the fit parameters and the source file as table metadata
        self.design.meta['MASKFILE'] = design_file
        self.design.meta['MASKOFF'] = register.par[0]
        self.design.meta['MASKSCL'] = register.par[1]
        # Fill the columns
        self.design['TRACEID'] = np.arange(nslits, dtype=self.design['TRACEID'].dtype)
        self.design['TRACESROW'] = np.full(nslits, self.pca.reference_row,
                                           dtype=self.design['TRACESROW'].dtype)
        self.design['TRACELPIX'] = self.spat_fit[self.pca.reference_row,self.traceid<0].astype(
                                        dtype=self.design['TRACELPIX'].dtype)
        self.design['TRACERPIX'] = self.spat_fit[self.pca.reference_row,self.traceid>0].astype(
                                        dtype=self.design['TRACERPIX'].dtype)
        self.design['SLITID'] = self.spectrograph.slitmask.slitid[slit_index].astype(
                                        dtype=self.design['SLITID'].dtype)
        self.design['SLITLFOC'] = register.mask_spat[register.match_index][self.traceid<0].astype(
                                        dtype=self.design['SLITLFOC'].dtype)
        self.design['SLITRFOC'] = register.mask_spat[register.match_index][self.traceid>0].astype(
                                        dtype=self.design['SLITRFOC'].dtype)
        if self.spectrograph.slitmask.onsky is not None:
            for i,key in enumerate(['SLITRA', 'SLITDEC', 'SLITLEN', 'SLITWID', 'SLITPA']):
                self.design[key] = self.spectrograph.slitmask.onsky[slit_index,i].astype(
                                        dtype=self.design[key].dtype)
        self.design['ALIGN'] = self.spectrograph.slitmask.alignment_slit[slit_index].astype(
                                        dtype=self.design['ALIGN'].dtype)

    def _fill_objects_table(self, register):
        """
        Fill :attr:`objects` based on the result of the design
        registration.

        Args:
            register (:class:`pypeit.spectrographs.slitmask.SlitRegister`):
                Object holding the result of the registration.
        """
        if self.spectrograph.slitmask.objects is None:
            # No object data available in slit mask design object
            self.objects = None
            return

        # Index for the slit in the design data
        slit_index = register.match_index[register.match_index % 2 == 0]//2
        # The index in the objects table are found by mapping the slit
        # index of each object in the design file to the slit index
        # included in the registration
        obj_index = utils.index_of_x_eq_y(self.spectrograph.slitmask.slitindx, slit_index,
                                          strict=True)
        # Number of objects
        nobj = len(obj_index)
        # Instantiate an empty table
        self.objects = EdgeTraceSet.empty_objects_table(rows=nobj)
        # Fill the columns
        for i,key in enumerate(['SLITID', 'OBJID', 'OBJRA', 'OBJDEC']):
                self.objects[key] = self.spectrograph.slitmask.objects[obj_index,i].astype(
                                        dtype=self.objects[key].dtype)
        # SLITINDX is the index of the slit in the `design` table, not
        # in the original slit-mask design data
        self.objects['SLITINDX'] = utils.index_of_x_eq_y(self.objects['SLITID'],
                                                         self.design['SLITID'], strict=True)


