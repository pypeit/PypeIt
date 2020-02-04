"""
The primary purpose of this module is to provide the classes/methods
used to trace slit edges.

For a command-line script that executes the automatic tracing, use
`pypeit_trace_edges`. As always, for a list of the script options,
run:

.. code-block:: bash

    $ pypeit_trace_edges -h

With a :ref:`pypeit_file`, a typical execution of the script would be:

.. code-block:: bash

    $ pypeit_trace_edges -f my_pypeit_file.pypeit

To show the trace results after completing each stage and/or to run
in debugging mode, use the `--show` and/or `--debug` options:

.. code-block:: bash

    $ pypeit_trace_edges -f my_pypeit_file.pypeit --debug --show

Programmatically, if you have a :ref:`pypeit_file` and a path for the reductions
(`redux_path`), an example of how to trace the slits in a single
detector is as follows::

    # Imports
    from pypeit.pypeit import PypeIt
    from pypeit import traceimage, edgetrace

    # Instantiate the PypeIt class to perform the necessary setup
    rdx = PypeIt(pypeit_file, redux_path=redux_path)

    # Find the trace frames files for a specific calibration group
    group = 0
    tbl_rows = rdx.fitstbl.find_frames('trace', calib_ID=group, index=True)
    files = rdx.fitstbl.frame_paths(tbl_rows)

    # Select a detector to trace
    det = 1

    # Setup the output paths for the trace file; these can be anything but
    # the defaults are below
    master_dir = rdx.par['calibrations']['caldir']
    master_key = rdx.fitstbl.master_key(tbl_rows[0], det=det)

    # Skip the bias subtraction, if reasonable; see
    # pypeit.biasframe.BiasFrame to construct a bias to subtract from
    # the TraceImage
    rdx.par['calibrations']['traceframe']['process']['bias'] = 'skip'

    # Construct the TraceImage
    traceImage = traceimage.TraceImage(rdx.spectrograph, files=files, det=det,
                                       par=rdx.par['calibrations']['traceframe'])
    traceImage.build_image()

    # Then run the edge tracing.  This performs the automatic tracing.
    edges = edgetrace.EdgeTraceSet(rdx.spectrograph, rdx.par['calibrations']['slitedges'],
                                   master_key=master_key, master_dir=master_dir, img=traceImage,
                                   det=det, auto=True)
    # You can look at the results using the show method:
    edges.show(thin=10, include_img=True, idlabel=True)
    # Or in ginga viewer
    edges.show(thin=10, in_ginga=True)
    # And you can save the results to a file
    edges.save()

If you want to instead start without a pypeit file, you could do the
following for, e.g., a single unbinned Keck DEIMOS flat-field
exposure in a fits file called `trace_file`::

    import os
    from pypeit import traceimage, edgetrace
    from pypeit.spectrographs.util import load_spectrograph

    spec = load_spectrograph('keck_deimos')
    par = spec.default_pypeit_par()
    par['calibrations']['traceframe']['process']['bias'] = 'skip'
    # Make any desired changes to the parameters here
    det = 3
    master_dir = os.path.split(os.path.abspath(trace_file))[0]
    master_key = 'test_trace_{0}'.format(det)

    traceImage = traceimage.TraceImage(spec, files=[trace_file], det=det,
                                       par=par['calibrations']['traceframe'])
    traceImage.build_image()

    edges = edgetrace.EdgeTraceSet(spec, par['calibrations']['slitedges'], master_key=master_key,
                                   master_dir=master_dir, img=traceImage, det=det, auto=True)
    edges.save()

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import time
import inspect
from collections import OrderedDict

from IPython import embed

import numpy as np

from scipy import ndimage

import matplotlib
from matplotlib import pyplot as plt
from matplotlib import ticker, rc

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy import table

from pypeit import msgs
from pypeit import utils
from pypeit import sampling
from pypeit import ginga
from pypeit import masterframe
from pypeit import io
from pypeit import slittrace
# TODO: Commented until EdgeTraceSet becomes a DataContainer
#from pypeit.datamodel import DataContainer
from pypeit.bitmask import BitMask
from pypeit.par.pypeitpar import EdgeTracePar
from pypeit.core import parse, pydl, procimg, pca, trace
from pypeit.traceimage import TraceImage
from pypeit.tracepca import TracePCA
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
                       ('NOEDGE', 'No edge found/input for this trace at this spatial column.'),
                    ('MATHERROR', 'A math error occurred during the calculation (e.g., div by 0)'),
                  ('MOMENTERROR', 'Recentering moment calculation had a large error'),
                   ('LARGESHIFT', 'Recentering resulted in a large shift'),
              ('OUTSIDEAPERTURE', 'Recentering yielded a centroid outside the moment aperture'),
                   ('EDGEBUFFER', 'Recentering yielded a centroid too close to the detector edge'),
                ('DISCONTINUOUS', 'Pixel included in a trace but part of a discontinuous segment'),
                    ('DUPLICATE', 'Trace is a duplicate based on trace matching tolerance'),
                   ('SHORTRANGE', 'Trace does not meet the minimum spectral range criterion'),
#                ('SHORTDETRANGE', 'Trace length does not meet trace detection threshold.'),
#                ('SHORTFITRANGE', 'Trace length does not meet fitting/PCA threshold.'),
                       ('HITMIN', 'Trace crosses the minimum allowed spatial column'),
                       ('HITMAX', 'Trace crosses the maximum allowed spatial column'),
                  ('OFFDETECTOR', 'Trace lands off, or within `det_buffer` of, the detector edge'),
                   ('USERINSERT', 'Trace was inserted as requested by user'),
                   ('SYNCINSERT', 'Trace was inserted during left and right edge sync'),
                    ('SYNCERROR', 'Trace synchronization error, likely due to edge effects'),
                   ('MASKINSERT', 'Trace was inserted based on drilled slit-mask locations'),
                 ('ORPHANINSERT', 'Trace was inserted to match an orphaned edge'),
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
        return ['USERINSERT', 'SYNCINSERT', 'MASKINSERT', 'ORPHANINSERT']

    @property
    def exclude_flags(self):
        # We will not exclude SYNCINSERT slits edges from the masking
        #   These can easily crop up for star boxes and at the edge of the detector
        return ['USERINSERT', 'MASKINSERT', 'ORPHANINSERT']


class EdgeTraceSet(masterframe.MasterFrame):
    r"""
    Core class that identifies, traces, and pairs edges in an image
    to define the slit apertures.

    The instantiation of the object can either be used to produce an
    empty placeholder that you then use multiple times to trace
    different images, or you can have the tracing begin immediately
    upon instantiation. For the latter, you must provide (at minimum)
    the `img` to trace, and the initialization then run
    :func:`initial_trace` or :func:`auto_trace`, depending on the
    value of `auto`. To automatically have the instantiation save the
    results, set `save=True` on instantiation.

    To load an existing master file with the result of a trace, you
    can either use the :attr:`from_file` method::

        edges = EdgeTraceSet.from_file(file)

    or instantiate the object as would have been done initially and
    then check if the trace exists and load it::
    
        edges = EdgeTraceSet(spec, par)
        if edges.exists:
            edges.load()

    or you can attempt to load directly on instantiation::

        edges = EdgeTraceSet(spec, par, load=True)

    In the latter case, note that the `load` argument takes
    precedence and an exception is raised if one provides both `img`
    and sets `load=True`.

    Most commonly, one will use the automatic tracing routine to
    trace the slit edges; see the description of the steps used
    during auto-tracing in the docs for :func:`auto_trace`.

    The success of the tracing critically depends on the parameters
    used. The defaults are tuned for each spectrograph based on
    testing using data in the pypeit development suite. See
    :ref:`pypeitpar` for the full documentation of the
    :class:`pypeit.par.pypeitpar.EdgeTracePar` parameters. Note that
    the :class:`pypeit.par.pypeitpar.TraceSlitsPar` parameter group
    has been deprecated.

    Finally, note that the :attr:`design` and :attr:`object` data are
    currently empty, as part of a development path for matching slits
    traced on the detector to slits expected from provided metadata.
    Once finished these objects will only contain data for
    spectrograph output files that provide the relevant metadata.
   
    .. todo:
        - Include a method/check that determines if traces cross one
        another anywhere.
        - Provide some guidance for the parameters to use.

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The object that sets the instrument used to take the
            observations. Used to set :attr:`spectrograph`.
        par (:class:`pypeit.par.pypeitpar.EdgeTracePar`):
            The parameters used to guide slit tracing. Used to set
            :attr:`par`.
        master_key (:obj:`str`, optional):
            The string identifier for the instrument configuration.  See
            :class:`pypeit.masterframe.MasterFrame`.
        master_dir (:obj:`str`, optional):
            Path to master frames.
        qa_path (:obj:`str`, optional):
            Directory for QA output. If None, no QA plots are
            provided.
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
            image.  Cannot be `None`.
        binning (`str`, optional):
            Comma-separated binning along the spectral and spatial
            directions following the PypeIt convention (e.g., '2,1').
            This is used to set the pixel scale of the image in
            arcsec per pixel, as needed for some assessments of the
            edge traces.
        auto (:obj:`bool`, optional):
            If a trace image is provided (`img`), run
            :func:`auto_trace` instead of :func:`initial_trace`.
        debug (:obj:`bool`, optional):
            Run in debug mode.
        show_stages (:obj:`bool`, optional):
            After ever stage of the auto trace prescription
            (`auto=True`), show the traces against the image using
            :func:`show`.
        save (:obj:`bool`, optional):
            Save the result to the master frame.
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
            The list of raw files used to construct the trace image
            (:attr:`img`). Only defined if argument `img` in
            :func:`initial_trace` or :func:`auto_trace` is a
            :class:`pypeit.traceimage.TraceImage` object.
        img (`numpy.ndarray`_):
            See argument list.
        bpm (`numpy.ndarray`_):
            See argument list.
        det (:obj:`int`):
            See argument list.
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
        binning (:obj:`str`, optional):
            On-detector binning of the data ordered spectral then
            spatial with format, e.g., `2,1`. Ignored if `img` is
            an instance of :class:`pypeit.traceimage.TraceImage`.
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
            Error in slit edge locations; measurements without errors
            have their errors set to -1.
        spat_msk (`numpy.ndarray`_):
            An integer array with the mask bits assigned to each
            trace centroid; see :class:`EdgeTraceBitMask`.
        spat_fit (`numpy.ndarray`_):
            A model fit to the `spat_cen` data.
        spat_fit_type (:obj:`str`):
            An informational string identifier for the type of model
            used to fit the trace data.
        pca (:obj:`list`, :class:`pypeit.tracepca.TracePCA`):
            Result of a PCA decomposition of the edge traces, used to
            predict new traces. This can either be a single
            :class:`pypeit.tracepca.TracePCA` object or a list of two
            :class:`pypeit.tracepca.TracePCA` objects if the PCA
            decomposition is peformed for the left (`pca[0]`) and
            right (`pca[1]`) traces separately.
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
            Directory for QA output. If None, no QA plots are
            provided.
        log (:obj:`list`):
            A list of strings indicating the main methods applied
            when tracing.
    """
    master_type = 'Edges'
    bitmask = EdgeTraceBitMask()    # Object used to define and toggle tracing mask bits
    def __init__(self, spectrograph, par, master_key=None, master_dir=None, qa_path=None,
                 img=None, bpm=None, det=1, binning=None, auto=False, debug=False,
                 show_stages=False, save=False, load=False):

        # TODO: It's possible for the master key and the detector
        # number to be inconsistent...
        masterframe.MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
                                         master_key=master_key, file_format='fits.gz')

        # TODO: Should make file_format a class attribute

        # TODO: Add type-checking for spectrograph and par
        self.spectrograph = spectrograph    # Spectrograph used to take the data
        self.par = par                      # Parameters used for slit edge tracing

        self.files = None               # Files used to construct the trace image
        self.img = None                 # The image used to find the slit edges
        self.bpm = None                 # Mask for the trace image
        self.det = None                 # Detector used for the trace image
        self.sobel_sig = None           # Sobel filtered image used to detect edges
        self.sobel_sig_left = None      # Sobel filtered image used to trace left edges
        self.sobel_sig_right = None     # Sobel filtered image used to trace right edges
        # TODO: Need a separate mask for the sobel image?
        self.nspec = None               # The shape of the trace image is (nspec,nspat)
        self.nspat = None
        self.binning = None             # Detector ordered spectral then spatial

        self.traceid = None             # The ID numbers for each trace
        self.spat_img = None            # (Integer) Pixel nearest the slit edge for each trace
        self.spat_cen = None            # (Floating-point) Spatial coordinate of the slit edges
                                        # for each spectral pixel
        self.spat_err = None            # Error in the slit edge spatial coordinate
        self.spat_msk = None            # Mask for the slit edge position for each spectral pixel

        self.spat_fit = None            # The result of modeling the slit edge positions
        self.spat_fit_type = None       # The type of fitting performed
        
        self.pca = None                 # One or two TracePCA objects with PCA decomposition
        self.pca_type = None            # Measurements used to construct the PCA (center or fit)

        self.design = None              # Table that collates slit-mask design data matched to
                                        # the edge traces
        self.objects = None             # Table that collates object information, if available
                                        # in the slit-mask design, matched to the `design` table.

        self.qa_path = qa_path          # Directory for QA output

        self.log = None                 # Log of methods applied

        if img is not None and load:
            msgs.error('Arguments img and load are mutually exclusive.  Choose to either trace '
                       'a new image or load a previous trace.')

        if load:
            # Attempt to load an existing master frame
            self.load()
        elif img is not None:
            # Provided a trace image so instantiate the object.
            if auto:
                self.auto_trace(img, bpm=bpm, det=det, binning=binning, save=save, debug=debug,
                                show_stages=show_stages)
            else:
                # JFH Is this option every used?
                self.initial_trace(img, bpm=bpm, det=det, binning=binning, save=save)

    def _reinit_trace_data(self):
        """
        Convenience method to set all attributes related to trace data to `None`.
        """
        self.traceid = None
        self.spat_img = None
        self.spat_cen = None
        self.spat_err = None
        self.spat_msk = None
        self.spat_fit = None
        self.spat_fit_type = None
        self.pca = None
        self.pca_type = None
        self.design = None
        self.objects = None

    @property
    def ntrace(self):
        """
        The number of edges (left and right) traced.
        """
        return 0 if self.traceid is None else self.traceid.size

    # ADDED in rmtdict
    @property
    def nslits(self):
        if self.is_synced:
            return self.ntrace//2
        # TODO: Maybe this should only throw a warning
        msgs.error('Number of slits undefined because edges are not left-right synchronized.')

    @staticmethod
    def empty_design_table(rows=None):
        """
        Construct an empty `design` table.

        Args:
            rows (:obj:`int`, optional):
                Number of table rows for each table column, expected
                to be the number of matched slits. If None, the table
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
                                 description='Slit used for alignment (1-yes; 0-no), not target '
                                             'observations.')
                           ])

    @staticmethod
    def empty_objects_table(rows=None):
        """
        Construct an empty `objects` table.

        Args:
            rows (:obj:`int`, optional):
                Number of table rows for each table column, expected
                to be the number of objects. If None, the table has
                empty columns.

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

    def rectify(self, flux, bpm=None, extract_width=None, mask_threshold=0.5, side='left'):
        r""""
        Rectify the provided image based on the current edge trace
        PCA model.

        The is primarily a wrapper for
        :func:`pypeit.sampling.rectify_image`; see its documentation
        for more detail.

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are
        `left_right_pca`.

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
                None, the image rectification is performed using
                :class:`pypeit.sampling.Resample` along each row
                (spectral position).
            mask_threshold (:obj:`float`, optional):
                Either due to `mask` or the bounds of the provided
                `flux`, pixels in the rectified image may not be fully
                covered by valid pixels in `flux`. Pixels in the
                output image with less than this fractional coverage
                by input pixels are flagged in the output.
            side (:obj:`str`, optional):
                If the PCA decomposition was performed for the left
                and right traces separately, this selects the PCA to
                use when constructing the rectification coordinate
                system. This is ignored if both were used in the PCA.

        Returns:
             Two `numpy.ndarray`_ objects are returned both with
             shape :math:`(N_{\rm spec}, N_{\rm spat})`, the rectified
             image and its boolean bad-pixel mask.
        """
        if self.pca is None:
            msgs.error('Must first run the PCA analysis for the traces; run build_pca.')
        pca = self.pca[0 if side == 'left' else 1] if self.par['left_right_pca'] else self.pca

        # Get the traces that cross the reference spatial position at
        # the first and last pixels of the image
        first_last_trace = pca.predict(np.array([0,self.nspat-1]))
        # Use these two traces to define the spatial pixel coordinates
        # to sample
        start = np.ceil(np.amax(np.amin(first_last_trace, axis=1))).astype(int)
        buffer = self.nspat - np.floor(np.amin(np.amax(first_last_trace, axis=1))).astype(int) \
                    + start
        # TODO: This has its limitations if the PCA is highly non-linear
        # Rectify the image
        ocol = np.arange(self.nspat+buffer)-start
        return sampling.rectify_image(flux, pca.predict(ocol), bpm=bpm, ocol=ocol,
                                      max_ocol=self.nspat-1, extract_width=extract_width,
                                      mask_threshold=mask_threshold)

    def auto_trace(self, img, bpm=None, det=1, binning=None, save=False, debug=False,
                   show_stages=False):
        r"""
        Execute a fixed series of methods to automatically identify
        and trace slit edges.

        The current algorithm is:
            - Detect and follow slit edges using :func:`initial_trace`.
            - Refine the measured centroids of the edge locations
              using :func:`centroid_refine`.
            - Fit the measured centroids with a polynomial using
              :func:`fit_refine`, which is basically a wrapper for
              :func:`pypeit.core.trace.fit_trace`. If a PCA
              decomposition of the traces is not possible because
              there are too few traces, the next two steps are
              skipped (skipping down to edge synchronization), and
              *the final measured slit edge traces* are the result of
              this step.
            - Construct a PCA decomposition of the fitted forms for
              the traces using :func:`pca_refine`. Traces are either
              decomposed as a group or split into separate left and
              right decompositions.
            - Use the PCA decomposition to rectify the trace image
              such that the slit edges should be aligned, collapse
              the image to get a high-signal-to-noise detection of
              each slit edge, and use these locations to predict,
              remeasure, and refit the slit edges; see
              :func:`peak_refine`. *The final measured slit edge
              traces* are based on the result of :func:`peak_refine`.
            - Synchronize the left and right traces into pairs that
              define slit apertures using :func:`sync`.
            - Use :func:`add_user_traces` and :func:`rm_user_traces`
              to add and remove traces as defined by the
              user-provided lists in the :attr:`par`.
            - If :attr:`qa_path` is not None, the QA plots are
              constructed and written; see :func:`qa_plot`.
            - Use :func:`save` to save the results, if requested.

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
                image.  Cannot be `None`.
            binning (:obj:`str`, optional):
                On-detector binning of the data ordered spectral then
                spatial with format, e.g., `2,1`. Ignored if `img` is
                an instance of :class:`pypeit.traceimage.TraceImage`.
            save (:obj:`bool`, optional):
                Save the result to the master frame.
            debug (:obj:`bool`, optional):
                Run in debug mode.
            show_stages (:obj:`bool`, optional):
                After ever stage of the auto trace, execute
                :func:`show` to show the results. These calls to show
                are always::

                    self.show(thin=10, include_img=True, idlabel=True)

        """
        # Perform the initial edge detection and trace identification
        self.initial_trace(img, bpm=bpm, det=det, binning=binning, save=False)
        if show_stages:
            self.show(thin=10, include_img=True, idlabel=True)

        # Initial trace can result in no edges found
        if not self.is_empty:
            # Refine the locations of the trace using centroids of the
            # features in the Sobel-filtered image.
            self.centroid_refine()
            if show_stages:
                self.show(thin=10, include_img=True, idlabel=True)

        # Initial trace can result in no edges found, or centroid
        # refinement could have removed all traces (via `check_traces`)
        if not self.is_empty:
            # Fit the trace locations with a polynomial
            self.fit_refine(debug=debug)
            # Use the fits to determine if there are any discontinous
            # trace centroid measurements that are actually components
            # of the same slit edge
            self.merge_traces(debug=debug)
            if show_stages:
                self.show(thin=10, include_img=True, idlabel=True)

        # Check if the PCA decomposition is possible; this should catch
        # long slits
        if self.par['auto_pca'] and self.can_pca():
            # Use a PCA decomposition to parameterize the trace
            # functional forms
            self.pca_refine(debug=debug)
            if show_stages:
                self.show(thin=10, include_img=True, idlabel=True)

            # Use the results of the PCA decomposition to rectify and
            # detect peaks/troughs in the spectrally collapsed
            # Sobel-filtered image, then use those peaks to further
            # refine the edge traces
            self.peak_refine(rebuild_pca=True, debug=debug)
            if show_stages:
                self.show(thin=10, include_img=True, idlabel=True)

        elif self.par['sync_predict'] == 'pca':
            # TODO: This causes the code to fault. Maybe there's a way
            # to catch this earlier on?
            msgs.error('Sync predict cannot use PCA because too few edges were found.  If you are '
                       'reducing multislit or echelle data, you may need a better trace image or '
                       'change the mode used to predict traces (see below).  If you are reducing '
                       'longslit data, make sure to set the sync_predict parameter to nearest: '
                       + msgs.newline() +
                       '    [calibrations]' + msgs.newline() +
                       '        [[slitedges]]' + msgs.newline() +
                       '            sync_predict = nearest')
#            self.par['sync_predict'] = 'nearest'

            # NOTE: If the PCA decomposition is possible, the
            # subsequent call to trace.peak_trace (called by
            # peak_refine) removes all the existing traces and replaces
            # them with traces at the peak locations. Those traces are
            # then fit, regardless of the length of the centroid
            # measurements. So coming out of peak_refine there are no
            # traces that are fully masked, which is not true if that
            # block of code isn't run. That means for the left-right
            # synchronization to work correctly, we have to remove
            # fully masked traces. This is done inside sync().

        # Left-right synchronize the traces
        self.sync()
        if show_stages:
            self.show(thin=10, include_img=True, idlabel=True)

        # First manually remove some traces, just in case a user
        # wishes to manually place a trace nearby a trace that
        # was automatically identified. One problem with adding
        # slits first is that we may have to sync the slits again.
        if self.par['rm_slits'] is not None:
            self.rm_user_traces(trace.parse_user_slits(self.par['rm_slits'], self.det, rm=True))

        # Add user traces
        if self.par['add_slits'] is not None:
            self.add_user_traces(trace.parse_user_slits(self.par['add_slits'], self.det))

        # TODO: Add a parameter and an if statement that will allow for
        # this.
        # `peak_refine` ends with the traces being described by a
        # polynomial. Instead finish by reconstructing the trace models
        # using the PCA decomposition
#        self.pca_refine(debug=debug)
#        if show_stages:
#            self.show(thin=10, include_img=True, idlabel=True)
            
        # TODO: Add mask_refine() when it's ready

        # Write the qa plots
        # TODO: Should maybe have a keyword argument that will allow
        # this to be skipped, even if the path is not None.
        if self.qa_path is not None:
            self.qa_plot(fileroot=self.file_name.split('.')[0])

        # Add this to the log
        self.log += [inspect.stack()[0][3]]
        if save:
            # Save the object to a file
            self.save()

    def initial_trace(self, img, bpm=None, det=1, binning='1,1', save=False):
        r"""
        Initialize the object for tracing a new image.

        This effectively reinstantiates the object and must be the
        first method called for tracing an image.  The algorithm:

            - Lightly boxcar smooths the trace image spectrally.
            - Replaces bad pixel columns, if a mask is provided.
            - Applies a Sobel filter to the trace image along the
              spatial axis (columns) to detect slit edges using steep
              positive gradients (left edges) and steep negative
              gradients (right edges). See
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
        `sobel_mode`, `edge_thresh`, `follow_span`,
        `valid_flux_thresh`, and `det_buffer`.

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
                image.  Cannot be `None`.
            binning (:obj:`str`, optional):
                On-detector binning of the data ordered spectral then
                spatial with format, e.g., `2,1`. Ignored if `img` is
                an instance of :class:`pypeit.traceimage.TraceImage`.
            save (:obj:`bool`, optional):
                Save the result to the master frame.
        """
        msgs.info('-'*50)
        msgs.info('{0:^50}'.format('Initialize Edge Tracing'))
        msgs.info('-'*50)

        # TODO: Add debugging argument and hooks
        # Parse the input based on its type
        if isinstance(img, TraceImage):
            self.files = img.file_list
            _img = img.pypeitImage.image
            self.binning = img.pypeitImage.binning
            # TODO: does TraceImage have a mask?
        else:
            _img = img
            self.binning = binning

        # TODO: keep the TraceImage object instead of deconstructing
        # it?  For direct input, use a base PypeItImage object

        # Check the input
        if _img.ndim != 2:
            msgs.error('Trace image must be 2D.')
        self.img = _img
        self.nspec, self.nspat = self.img.shape
        self.bpm = np.zeros((self.nspec, self.nspat), dtype=bool) if bpm is None else bpm
        if self.bpm.shape != self.img.shape:
            msgs.error('Mask is not the same shape as the trace image.')
        if det is None or det < 1:
            msgs.error('Detector must be an integer that is >=1.')
        self.det = det

        # Lightly smooth the image before using it to trace edges
        # TODO: Make this filter size a parameter?
        # NOTE: This was _make_binarr()
        _img = ndimage.uniform_filter(self.img, size=(3, 1), mode='mirror')

        # Replace bad-pixel columns if they exist
        # TODO: Do this before passing the image to this function?
        # Instead of just replacing columns, replace all bad pixels...
        # NOTE: This was previously done at the beginning of
        # edgearr_from_binarr
        if np.any(self.bpm):

            # TODO: For tracing, we really only care about bad spec
            # lines, i.e. columns in the PypeIt frame. And this is how
            # BPM is oriented. I am now setting it to deal with both
            # options. Do we need to replace bad *rows* instead of bad
            # columns?

            #flip = self.spectrograph.raw_is_transposed(det=self.det)
            #axis = 1 if flip else 0

            # Replace bad columns that cover more than half the image
            for flip, axis in zip([False,True], [0,1]):
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
                = trace.detect_slit_edges(_img, bpm=self.bpm, median_iterations=self.par['filt_iter'],
                                          sobel_mode=self.par['sobel_mode'],
                                          sigdetect=self.par['edge_thresh'])
        # Empty out the images prepared for left and right tracing
        # until they're needed.
        self.sobel_sig_left = None
        self.sobel_sig_right = None

        # Identify traces by following the detected edges in adjacent
        # spectral pixels.
        # TODO: Use the user parameter instead of the hard-coded
        # number?
#        minimum_spec_length = self.nspec * self.par['det_min_spec_length']
        minimum_spec_length = 50
        _trace_id_img = trace.identify_traces(edge_img, follow_span=self.par['follow_span'],
                                              minimum_spec_length=minimum_spec_length)

        ################################################################
        #
        # TODO: Make this an option
#        # Update the traces by handling single orphan edges and/or
#        # traces without any left or right edges.
#        # NOTE: This combines the functionality of
#        # edgearr_add_left_right and edgearr_final_left_right
#
#        # TODO: Why is this needed? Comments say that this is primarily
#        # in place for LRISb, but can this be handled better with sync()?
#
#        flux_valid = np.median(_img) > self.par['valid_flux_thresh']
#        trace_id_img = trace.handle_orphan_edges(_trace_id_img, self.sobel_sig, bpm=self.bpm,
#                                                  flux_valid=flux_valid,
#                                                  buffer=self.par['det_buffer'], copy=True)
#
#        # Flag any inserted edges
#        # TODO: handle_orphan_edge can delete or insert edges; deleted
#        # edges are ignored hereafter, inserted edges are flagged
#        inserted_edge = (trace_id_img != 0) & (trace_id_img != _trace_id_img)
        trace_id_img = _trace_id_img
        inserted_edge = np.zeros(trace_id_img.shape, dtype=bool)
        ################################################################

        # Check that edges were found
        if np.all(trace_id_img == 0):
            msgs.warn('No edges found!  Trace data will be empty.')
            self._reinit_trace_data()
            self.log = [inspect.stack()[0][3]]
            if save:
                self.save()
            return

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

            # Flag any insert traces
            row, col = np.where(np.invert(trace_id_img.mask) & inserted_edge
                                    & (trace_id_img.data == self.traceid[i]))
            if len(row) > 0:
                self.spat_msk[row,i] = self.bitmask.turn_on(self.spat_msk[row,i], 'ORPHANINSERT')

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

    def save(self, outfile=None, overwrite=True, checksum=True, float_dtype='float32'):
        """
        Save the trace object for later re-instantiation.

        Args:
            outfile (:obj:`str`, optional):
                Name for the output file.  Defaults to
                :attr:`master_file_path`.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file.
            checksum (:obj:`bool`, optional):
                Passed to `astropy.io.fits.HDUList.writeto`_ to add
                the DATASUM and CHECKSUM keywords fits header(s).
            float_dtype (:obj:`str`, optional):
                Convert floating-point data to this data type before
                writing.  Default is 32-bit precision.
        """

        _outfile = self.master_file_path if outfile is None else outfile
        # Check if it exists
        if os.path.exists(_outfile) and not overwrite:
            msgs.error('Master file exists: {0}'.format(_outfile) + msgs.newline()
                       + 'Set overwrite=True to overwrite it.')

        # Report name before starting to write it
        msgs.info('Writing master frame to {0}'.format(_outfile))

        # Build the primary header
        #   - Initialize with basic metadata
        prihdr = self.build_master_header()
        #   - Add the qa path
        prihdr['QADIR'] = (self.qa_path, 'PypeIt: QA directory')
        #   - Add metadata specific to this class -- This is added as PYP_SPEC in the build_master_header() call
        #prihdr['SPECT'] = (self.spectrograph.spectrograph, 'PypeIt: Spectrograph Name')
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
        if self.log is not None and len(self.log) > 0:
            ndig = int(np.log10(len(self.log)))+1
            for i,m in enumerate(self.log):
                prihdr['LOG{0}'.format(str(i+1).zfill(ndig))] \
                        = (m, '{0}: Completed method'.format(self.__class__.__name__))
        #   - PCA type, used for rebuilding the PCA when loading
        prihdr['PCATYPE'] = ('None' if self.pca is None else self.pca_type,
                             'PypeIt: Edge trace PCA type')
        #   - The dispersion name is needed to setup some spectrographs
        if self.spectrograph.dispname is not None:
            prihdr['DISPNAME'] = (self.spectrograph.dispname, 'Spectrograph disperser')
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
                            fits.ImageHDU(data=self.img.astype(float_dtype), name='TRACEIMG'),
                            fits.ImageHDU(data=self.bpm.astype(np.int16), name='TRACEBPM'),
                            fits.ImageHDU(data=self.sobel_sig.astype(float_dtype),
                                          name='SOBELSIG'),
                            fits.ImageHDU(data=self.traceid, name='TRACEID'),
                            fits.ImageHDU(data=self.spat_cen.astype(float_dtype), name='CENTER'),
                            fits.ImageHDU(data=self.spat_err.astype(float_dtype),
                                          name='CENTER_ERR'),
                            fits.ImageHDU(header=mskhdr, data=self.spat_msk, name='CENTER_MASK'),
                            fits.ImageHDU(header=fithdr, data=self.spat_fit.astype(float_dtype),
                                          name='CENTER_FIT')])
        if self.pca is not None:
            if self.par['left_right_pca']:
                hdu += [self.pca[0].to_hdu(name='LPCA'), self.pca[1].to_hdu(name='RPCA')]
            else:
                hdu += [self.pca.to_hdu()]
        # TODO: These things are going to go in slittrace.SlitTraceSet
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

    def load(self, validate=True, rebuild_pca=False):
        """
        Load and reinitialize the trace data.

        Data is read from :attr:`master_file_path` and used to
        overwrite any internal data. Specific comparisons of the
        saved data are performed to ensure the file is consistent
        with having been written by a consistent version of the code;
        see :func:`_reinit`.

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
                Rebuild the PCA decomposition of the traces based on
                the loaded trace data, but only if the PCA had been
                originally produced for the saved object (as
                indicated by the fits header). Otherwise, any
                existing saved PCA data will be used to construct the
                PCA object(s). If the header indicates that the PCA
                was not originally performed, this is ignored.
        Raises:
            FileNotFoundError:
                Raised if no data has been written for this master
                frame.
            ValueError:
                Raised if validation of the data fails (actually
                raised by :func:`_reinit`).
        """
        # Check the file exists
        if not self.exists:
            msgs.error('File does not exit: {0}'.format(self.master_file_path))
        with fits.open(self.master_file_path) as hdu:
            # Re-initialize and validate
            self._reinit(hdu, validate=validate, rebuild_pca=rebuild_pca)

    @classmethod
    def from_file(cls, filename, rebuild_pca=False):
        """
        Instantiate using data from a file.

        To reload data that has been saved for an existing
        instantiation, use :func:`load`.

        Args:
            filename (:obj:`str`):
                Fits file produced by :func:`save`.
            rebuild_pca (:obj:`bool`, optional):
                Rebuild the PCA decomposition of the traces based on
                the loaded trace data, but only if the PCA had been
                originally produced for the saved object (as
                indicated by the fits header). Otherwise, any
                existing saved PCA data will be used to construct the
                PCA object(s). If the header indicates that the PCA
                was not originally performed, this is ignored.
        """

        # TODO: Consolidate this with items_from_master_file in
        # masterframe.py?

        # Check the file exists
        if not os.path.isfile(filename):
            msgs.error('File does not exit: {0}'.format(filename))
        msgs.info('Loading EdgeTraceSet data from: {0}'.format(filename))
        with fits.open(filename) as hdu:
            this = cls(load_spectrograph(hdu[0].header['PYP_SPEC']),
                       EdgeTracePar.from_header(hdu[0].header),
                       master_key=hdu[0].header['MSTRKEY'], master_dir=hdu[0].header['MSTRDIR'],
                       qa_path=hdu[0].header['QADIR'])

            # Re-initialize and validate
            this._reinit(hdu, rebuild_pca=rebuild_pca)
        return this

    def _reinit(self, hdu, validate=True, rebuild_pca=False):
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
                Rebuild the PCA decomposition of the traces based on
                the loaded trace data, but only if the PCA had been
                originally produced for the saved object (as
                indicated by the fits header). Otherwise, any
                existing saved PCA data will be used to construct the
                PCA object(s). If the header indicates that the PCA
                was not originally performed, this is ignored.
        """
        # Read and assign data from the fits file
        self.files = io.parse_hdr_key_group(hdu[0].header, prefix='RAW')
        if len(self.files) == 0:
            self.files = None

        # TODO: These now case back to float64, regardless of how they
        # were written (typically float32). Not sure this is a great
        # idea, but I put it here so that get_slits always casts to
        # float64.  Not exactly sure what make sense...
        self.img = hdu['TRACEIMG'].data.astype(float)
        self.nspec, self.nspat = self.img.shape
        self.bpm = hdu['TRACEBPM'].data.astype(bool)
        self.det = hdu[0].header['DET']
        self.binning = hdu[0].header['BINNING']
        self.sobel_sig = hdu['SOBELSIG'].data.astype(float)
        self.traceid = hdu['TRACEID'].data
        self.spat_cen = hdu['CENTER'].data.astype(float)
        self.spat_err = hdu['CENTER_ERR'].data.astype(float)
        self.spat_msk = hdu['CENTER_MASK'].data
        self.spat_fit = hdu['CENTER_FIT'].data.astype(float)
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
        self.spat_img = None if self.traceid is None \
                            else (np.round(self.spat_cen if self.spat_fit is None
                                    else self.spat_fit).astype(int))

        # Read or rebuild the PCA if it existed previously
        self.pca_type = None if hdu[0].header['PCATYPE'] == 'None' else hdu[0].header['PCATYPE']
        if self.pca_type is None:
            self.pca = None
        elif rebuild_pca:
            if self.can_pca():
                self._reset_pca(True)
            else:
                # TODO: Should this throw a warning instead?
                msgs.error('Traces do not meet necessary criteria for the PCA decomposition.')
        else:
            self.pca = [ TracePCA.from_hdu(hdu[ext]) for ext in ['LPCA', 'RPCA']] \
                            if self.par['left_right_pca'] else TracePCA.from_hdu(hdu['PCA'])

        # Read the disperser if possible
        if self.spectrograph.dispname is None:
            try:
                self.spectrograph.dispname = hdu[0].header['DISPNAME']
            except:
                # Assume the spectrgrograph *NEVER* defines the disperser
                pass

        # Read the log
        self.log = io.parse_hdr_key_group(hdu[0].header, prefix='LOG')

        # Reinit the left and right sobel images
        self.sobel_sig_left = None
        self.sobel_sig_right = None

        # Finished, if not validating
        if not validate:
            return

        # Check the package versions used to create the file
        if not io.header_version_check(hdu['PRIMARY'].header):
            msgs.warn('This file was written with different package versions.  You may need to '
                      'redo the reductions to use the file within the current environment!')

        # Test the bitmask has the same keys and key values
        hdr_bitmask = BitMask.from_header(hdu['CENTER_MASK'].header)
        if hdr_bitmask.bits != self.bitmask.bits:
            msgs.warn('The bitmask in this fits file appear to be out of date!  Will continue '
                      'by using old bitmask but errors may occur.  You should recreate this '
                      'master frame.')
            self.bitmask = hdr_bitmask

        # Test the spectrograph is the same
        if self.spectrograph.spectrograph != hdu[0].header['PYP_SPEC']:
            msgs.error('Data used for this master frame was from a different spectrograph!')

        # Check that the disperser is correct. NOTE: If disperser was
        # None when this method was called, this will automatically be
        # true.
        if self.spectrograph.dispname is not None \
                and self.spectrograph.dispname != hdu[0].header['DISPNAME']:
            msgs.error('Current spectrograph disperser does not match file data.')

        # Test the parameters used are the same
        par = EdgeTracePar.from_header(hdu[0].header)
        if self.par.data != par.data:
            # TODO: The above inequality works for non-nested ParSets,
            # but will need to be more careful for nested ones, or just
            # avoid writing nested ParSets to headers...
            msgs.error('This master frame was generated using different parameter values!')

    @property
    def flagged_bits(self):
        """
        List the flags that have been tripped.

        Returns:
            list: List of the unique flags tripped by any trace
            datum.
        """
        return [] if self.spat_msk is None \
                    else np.unique(np.concatenate([self.bitmask.flagged_bits(b) 
                                                    for b in np.unique(self.spat_msk)])).tolist()

    # TODO: Break the ginga commands out into a separate method?
    ## TODO It is confusing that this show routine shows images flipped from the PypeIt convention.
    ## It should be rewritten to show images with spectral direction vertical like all our other QA.
    def show(self, include_error=False, thin=1, in_ginga=False, include_img=False,
             include_sobel=False, img_buffer=100, flag=None, idlabel=False, slits=None,
             original=False):
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
                Display the trace data in a ginga. If a ginga window
                is not open, this instantiates one; otherwise, the
                existing ginga window is cleared. The trace image is
                shown in one window and the sobel image is shown in a
                second. The edge traces themselves are shown in both
                windows. Any masking is ignored except that any
                traces that are fully masked (see
                :func:`fully_masked_traces`) are *not* shown. If a
                SlitTraceSet object is *not* provided, the data shown
                is the modeled results for the trace if it exists,
                and the measured trace locations otherwise.
            include_img (:obj:`bool`, optional):
                Overlay the trace data on the trace image; mutually
                exclusive with `include_sobel`.
            include_sobel (:obj:`bool`, optional):
                Overlay the trace data on the Sobel-filtered trace
                image; mutually exclusive with `include_img`.
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
            idlabel (:obj:`bool`, optional):
                Label each trace by their ID numbers.
            slits (:class:`pypeit.slittrace.SlitTraceSet`, optional):
                Slits to plot instead of those kept internally. Note
                that, if this is provided, the "modeled" and
                "measured" slit locations are identical.
            original (:obj:`bool`, optional):
                When ``slits`` are provided and tweaked slits are
                available, this selects which traces to show. If
                True, show the original slit traces; if False, show
                the tweaked ones.
        """
        if include_img and include_sobel:
            msgs.error('Cannot show both the trace image and the filtered version.')

        # Build the slit edge data to plot.
        if slits is None:
            # Use the internals. Any masked data is excluded; masked
            # data to be plotted are held in a separate array. This
            # means that errors and fits are currently never plotted
            # for masked data.
            _flag = None if flag in ['any', None] else np.atleast_1d(flag)
            cen = np.ma.MaskedArray(self.spat_cen, mask=self.bitmask.flagged(self.spat_msk))
            fit = self.spat_fit
            err = np.ma.MaskedArray(self.spat_err, mask=np.ma.getmaskarray(cen).copy())
            msk = None if flag is None \
                    else np.ma.MaskedArray(self.spat_cen, mask=np.invert(
                                           self.bitmask.flagged(self.spat_msk, flag=_flag)))
            is_left = self.is_left
            is_right = self.is_right
            _include_error = include_error
            gpm = np.invert(self.fully_masked_traces(flag=self.bitmask.bad_flags,
                                                     exclude=self.bitmask.exclude_flags))
            traceid = self.traceid
            nslits = np.amax(traceid)   # Only used if synced is True
            synced = self.is_synced
        else:
            # Use the provided SlitTraceSet
            _include_error = False
            if include_error:
                msgs.warn('SlitTraceSet object has no errors.')
            left, right = slits.select_edges(original=original)
            cen = np.hstack((left,right))
            fit = cen
            msk = None
            nslits = slits.nslits
            is_left = np.ones(2*nslits, dtype=bool)
            is_left[nslits:] = False
            is_right = np.invert(is_left)
            gpm = np.ones(2*nslits, dtype=bool)
            traceid = np.concatenate((-np.arange(nslits), np.arange(nslits)))
            synced = True

        if in_ginga:
            # Set up the appropriate keyword arguments for the IDs
            id_kwargs = {'slit_ids': np.arange(nslits)+1} if synced \
                            else {'left_ids': traceid[gpm & is_left],
                                  'right_ids': traceid[gpm & is_right]}
            _trc = cen if fit is None else fit

            # Connect to or instantiate ginga window
            ginga.connect_to_ginga(raise_err=True, allow_new=True)
            # Clear the viewer and show the trace image
            trace_viewer, trace_ch = ginga.show_image(self.img, chname='Trace Image', clear=True)
            if not self.is_empty:
                ginga.show_slits(trace_viewer, trace_ch, _trc[:,gpm & is_left],
                                 _trc[:,gpm & is_right], pstep=thin, synced=synced, **id_kwargs)

            # Show the Sobel sigma image (do *not* clear)
            sobel_viewer, sobel_ch = ginga.show_image(self.sobel_sig, chname='Sobel Filtered')
            if not self.is_empty:
                ginga.show_slits(sobel_viewer, sobel_ch, _trc[:,gpm & is_left],
                                 _trc[:,gpm & is_right], pstep=thin, synced=synced, **id_kwargs)
            return

        # Show the traced image
        if include_img:
            img_zlim = utils.growth_lim(self.img, 0.95, fac=1.05)
            plt.imshow(self.img, origin='lower', interpolation='nearest', aspect='auto',
                       vmin=img_zlim[0], vmax=img_zlim[1])
        elif include_sobel:
            sob_zlim = utils.growth_lim(self.sobel_sig, 0.95, fac=1.05)
            plt.imshow(self.sobel_sig, origin='lower', interpolation='nearest', aspect='auto',
                       vmin=sob_zlim[0], vmax=sob_zlim[1])

        if self.is_empty:
            msgs.info('No traces defined.')
            plt.xlim(-img_buffer, self.nspat+img_buffer)
            plt.ylim(-img_buffer, self.nspec+img_buffer)
            plt.ylabel('Spectral pixel index')
            plt.xlabel('Spatial pixel index')
            plt.show()
            return

        # Spectral position
        spec = np.tile(np.arange(self.nspec), (self.ntrace,1)).T

        if _include_error:
            indx = is_left | is_right
            # Show the errors
            plt.errorbar(cen[::thin,indx].ravel(), spec[::thin,indx].ravel(),
                         xerr=err[::thin,indx].ravel(), fmt='none', ecolor='k', elinewidth=0.5,
                         alpha=0.3, capthick=0, zorder=3)
        # Plot the left trace points
        plt.scatter(cen[::thin,is_left], spec[::thin,is_left], marker='.', color='k', s=30,
                    lw=0, zorder=4, label='left edge measurements', alpha=0.8)
        if msk is not None:
            plt.scatter(msk[::thin,is_left], spec[::thin,is_left], marker='x', color='C3', s=20,
                        lw=0.5, zorder=5, label='masked left edges', alpha=0.8)

        # Plot the right trace points
        plt.scatter(cen[::thin,is_right], spec[::thin,is_right], marker='.', color='0.7',
                    s=30, lw=0, zorder=4, label='right edge measurements', alpha=0.8)
        if msk is not None:
            plt.scatter(msk[::thin,is_right], spec[::thin,is_right], marker='x', color='C1', s=20,
                        lw=0.5, zorder=5, label='masked right edges', alpha=0.8)

        if idlabel:
            # Label the traces by their ID number
            for i in range(self.ntrace):
                indx = np.invert(cen.mask[:,i])
                if not np.any(indx):
                    continue
                _spec = spec[:,i][indx]
                _cen = cen[:,i][indx]
                plt.text(_cen[_spec.size//2], _spec[_spec.size//2], str(self.traceid[i]),
                         color='k', fontsize=16, alpha=0.7, zorder=10)

        if fit is None:
            # No fits, so we're done
            plt.legend()
            plt.show()
            return

        # Plot the trace fits
        for i in range(self.ntrace):
            if not gpm[i]:
                continue
            # If statement structure primarily for the labels. Only
            # difference between left and right is the color.
            if is_left[i]:
                left_line = plt.plot(fit[::thin,i], spec[::thin,i], color='C3', lw=1, zorder=6)
            elif is_right[i]:
                right_line = plt.plot(fit[::thin,i], spec[::thin,i], color='C1', lw=1, zorder=6)

            if idlabel and np.all(cen.mask[:,i]):
                plt.text(spec[self.nspec//2,i], fit[self.nspec//2,i], str(traceid[i]), color='k',
                         fontsize=16, alpha=0.7, zorder=10)

        # Limits and labels
        plt.xlim(-img_buffer, self.nspat+img_buffer)
        plt.ylim(-img_buffer, self.nspec+img_buffer)
        if np.any(is_left & gpm):
            left_line[0].set_label('left edge fit')
        if np.any(is_right & gpm):
            right_line[0].set_label('right edge fit')
        plt.ylabel('Spectral pixel index')
        plt.xlabel('Spatial pixel index')
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
        if self.is_empty:
            msgs.error('No traces for QA plot.')

        # Restore matplotlib defaults
        # TODO: Is this going to screw up later plots?
        matplotlib.rcParams.update(matplotlib.rcParamsDefault)

        # Set font size
        rc('font', size=8)

        # Spectral pixel coordinate vector and global plot limits
        spec = np.arange(self.nspec)
        ylim = [-1,self.nspec]
        img_zlim = utils.growth_lim(self.img, 0.95, fac=1.05)
        sob_zlim = utils.growth_lim(self.sobel_sig, 0.95, fac=1.05)

        # Set figure
        w,h = plt.figaspect(1)
        fig = plt.figure(figsize=(1.5*w,1.5*h))

        # Grid for plots
        n = np.array([4,3])
        buff = np.array([0.05, 0.02])
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
            ax_x0 = strt[0]+ii*(buff[0]+delt[0])
            ax_y = strt[1]+(n[1]-jj-1)*(buff[1]+delt[1])

            # Spatial pixel plot limits for this trace
            indx = np.invert(self.bitmask.flagged(self.spat_msk[:,i], flag=self.bitmask.bad_flags))
            xlim = utils.growth_lim(self.spat_cen[indx,i], 1.0, fac=2.0)
            if min_spat is not None and np.diff(xlim) < min_spat:
                xlim = np.sum(xlim)/2 + np.array([-1,1])*min_spat/2

            # Plot the trace image and the fit (if it exists)
            ax = fig.add_axes([ax_x0, ax_y, delt[0]/2, delt[1]])
            ax.minorticks_on()
            ax.tick_params(which='major', length=10, direction='in', top=True, right=True)
            ax.tick_params(which='minor', length=5, direction='in', top=True, right=True)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.imshow(self.img, origin='lower', interpolation='nearest', vmin=img_zlim[0],
                      vmax=img_zlim[1], aspect='auto')
            if self.spat_fit is not None:
                ax.plot(self.spat_fit[:,i], spec, color='C3' if self.traceid[i] < 0 else 'C1')
            ax.text(0.07, 0.93,'{0}'.format(self.traceid[i]), ha='left', va='center',
                    transform=ax.transAxes, fontsize=12, bbox=dict(facecolor='white', alpha=0.3))
            if ii == 0:
                ax.text(-0.55, 0.5, 'Spectral Coordinate (pix)', ha='center', va='center',
                        transform=ax.transAxes, rotation='vertical')

            # Plot the filtered image and the fit (if it exists)
            ax = fig.add_axes([ax_x0 + delt[0]/2, ax_y, delt[0]/2, delt[1]])
            ax.minorticks_on()
            ax.tick_params(which='major', length=10, direction='in', top=True, right=True)
            ax.tick_params(which='minor', length=5, direction='in', top=True, right=True)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.yaxis.set_major_formatter(ticker.NullFormatter())
            ax.imshow(self.sobel_sig, origin='lower', interpolation='nearest', vmin=sob_zlim[0],
                      vmax=sob_zlim[1], aspect='auto')
            if self.spat_fit is not None:
                ax.plot(self.spat_fit[:,i], spec, color='C3' if self.traceid[i] < 0 else 'C1')
            if jj == n[1]-1:
                ax.text(0.0, -0.1, 'Spatial Coordinate (pix)', ha='center', va='center',
                        transform=ax.transAxes)

            # Prepare for the next trace plot
            j += 1
            if j == np.prod(n) or i == self.ntrace-1:
                j = 0
                if fileroot is None:
                    plt.show()
                else:
                    page += 1
                    ofile = os.path.join(self.qa_path, 'PNGs',
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
                self.sobel_sig_left = trace.prepare_sobel_for_trace(self.sobel_sig, bpm=self.bpm, boxcar=boxcar,
                                                                    side='left')
            return self.sobel_sig_left
        if side == 'right':
            if self.sobel_sig_right is None:
                self.sobel_sig_right = trace.prepare_sobel_for_trace(self.sobel_sig, bpm=self.bpm, boxcar=boxcar,
                                                                     side='right')
            return self.sobel_sig_right
        msgs.error('Side must be left or right.')

    def centroid_refine(self, follow=True, start_indx=None, continuous=False, check=True,
                        use_fit=False):
        """
        Refine the edge positions using a moment analysis and assess
        the results.

        The method runs in two primary modes, depending on the
        value of `follow`:

        When `follow=False`, the function simply executes
        :func:`pypeit.core.trace.masked_centroid` to recenter
        *all* the locations for the left and right edges at each
        spectral position (row) independently. The maximum shift
        between any input and output centroid is `max_shift_abs` from
        :attr:`par`.

        When `follow=True`, the method uses
        :func:`pypeit.core.trace.follow_centroid` to recenter each
        edge starting from the `start_indx` spectral index position
        (row) and then moving to higher and lower spectral rows. In
        this case, the center of the aperture used for each spectral
        row is the centroid of the trace measured for the preceding
        spectral row. The maximum shift between the input and output
        center for the first position analyzed is set by
        `max_shift_abs` and the maximum shift for each subsequent
        spectral row is set by `max_shift_adj`; both parameters are
        provided in :attr:`par`. In this approach, it's typically
        best to let the method determine the starting spectral row
        instead of providing it directly. If left to its own devices,
        it will iterate through all the traces collecting those that
        all cross a specific spectral row into groups that it can
        follow simultaneously. If a starting specral row is provided
        directly, all traces must cross that row.

        Regardless of the value of `follow`,
        :func:`pypeit.core.trace.masked_centroid` is run with uniform
        weighting and an aperture width set to be twice
        `fwhm_uniform` from :attr:`par`. Other used parameters from
        :attr:`par` (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are
        `max_shift_abs`, `max_shift_adj`, `max_spat_error`, and
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
                `start_indx`) and then move to higher and lower
                spectral rows in series. See
                :func:`pypeit.core.trace.follow_centroid`.
            start_indx (:obj:`int`, optional):
                The index of the starting spectral row when following
                the trace between adjacent rows; see `follow`. If
                None, the starting spectral row is set by finding the
                spectral row that crosses the most unmasked trace
                positions; see
                :func:`pypeit.core.trace.most_common_trace_row`.
                Value is ignored if `follow=False`.
            continuous (:obj:`bool`, optional):
                Keep only the continuous part of the traces from the
                starting spectral row.
            check (:obj:`bool`, optional):
                Run :func:`check_traces` to flag bad traces.
            use_fit (:obj:`bool`, optional):
                Use the fitted traces as the starting point for the
                refinement. Otherwise, uses :`spat_img`. If True and
                :attr:`spat_fit` is None, the method will raise an
                exception.
        """
        # Check that there are traces to refine!
        if self.is_empty:
            # TODO: Continue to have this fault?
            msgs.error('No traces are defined.')

        # Check input
        if use_fit and self.spat_fit is None:
            msgs.error('No fit data available.')

        # Parse parameters and report
        width = 2 * self.par['fwhm_uniform']
        maxshift_start = self.par['max_shift_abs']
        maxshift_follow = self.par['max_shift_adj']
        maxerror = self.par['max_spat_error']
        minimum_spec_length = self.par['det_min_spec_length']*self.nspec

        msgs.info('-'*50)
        msgs.info('{0:^50}'.format('Edge Centroid Refinement'))
        msgs.info('-'*50)
        msgs.info('Width of window for centroiding the edges: {0:.1f}'.format(width))
        msgs.info('Max shift between spectrally adjacent pixels: {0:.2f}'.format(maxshift_follow))
        msgs.info('Max centroid error: {0}'.format(maxerror))
        msgs.info('Minimum spectral pixels for a valid trace: {0}'.format(minimum_spec_length))
    
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

        # Ignore inserted traces
        inserted = self.fully_masked_traces(flag=self.bitmask.insert_flags)
        if np.any(inserted):
            cen[:,inserted] = self.spat_cen[:,inserted]
            err[:,inserted] = self.spat_err[:,inserted]
            msk[:,inserted] = self.spat_msk[:,inserted]

        # Refine left then right
        for side in ['left', 'right']:
            
            # Get the image relevant to tracing this side
            _sobel_sig = self._side_dependent_sobel(side)

            # Find the traces to refine, must be on the correct side
            # and must not have been inserted
            indx = (self.is_left if side == 'left' else self.is_right) & np.invert(inserted)
            if not np.any(indx):
                continue

            msgs.info('Found {0} {1} edge trace(s) to refine'.format(np.sum(indx), side))

            if follow:
                # Find the bad trace positions
                bpm = self.bitmask.flagged(self.spat_msk, flag=self.bitmask.bad_flags)
                untraced = indx.copy()
                while np.any(untraced):
                    # Get the starting row
                    _start_indx = trace.most_common_trace_row(bpm[:,untraced], valid_frac=1.) \
                                        if start_indx is None else start_indx
                    # Select the edges to follow
                    to_trace = untraced & np.invert(bpm[_start_indx,:])
                    if not np.any(to_trace):
                        # Something has gone wrong
                        # TODO: Get rid of this when convinced it won't
                        # get tripped...
                        msgs.error('Traces remain but could not select good starting position.')

                    ## TODO row and column should not be used here in the output. Adopt the PypeIt convention spec, spat
                    msgs.info('Following {0} {1} edge(s) '.format(np.sum(to_trace), side)
                              + 'from row {0}; '.format(_start_indx)
                              + '{0} trace(s) remain.'.format(np.sum(untraced)-np.sum(to_trace)))
                    # Follow the centroid of the Sobel-filtered image
                    cen[:,to_trace], err[:,to_trace], msk[:,to_trace] \
                            = trace.follow_centroid(_sobel_sig, _start_indx,
                                                    spat[_start_indx,to_trace],
                                                    ivar=ivar, bpm=_bpm, fwgt=fwgt, width=width,
                                                    maxshift_start=maxshift_start,
                                                    maxshift_follow=maxshift_follow,
                                                    maxerror=maxerror, continuous=continuous,
                                                    bitmask=self.bitmask)
                    # Update untraced
                    untraced[to_trace] = False
            else:
                cen[:,indx], err[:,indx], msk[:,indx] \
                        = trace.masked_centroid(_sobel_sig, spat[:,indx], width, ivar=ivar,
                                                bpm=_bpm, fwgt=fwgt, maxshift=maxshift_start,
                                                maxerror=maxerror, bitmask=self.bitmask,
                                                fill='bound')

            # Include previous mask in result
            # TODO: This means if the edge was not detected by
            # detect_slit_edges, it won't get used here. Need to check
            # how this meshes with the previous behavior.
            msk[:,indx] |= self.spat_msk[:,indx]

            # Check the traces
            min_spatial = None if side == 'left' else 0
            max_spatial = _sobel_sig.shape[1]-1 if side == 'left' else None
            good, bad = self.check_traces(cen=cen, msk=msk, subset=indx, min_spatial=min_spatial,
                                          max_spatial=max_spatial,
                                          minimum_spec_length=minimum_spec_length)

            # Save the results: only keep the good centroids, but merge all the masks
            self.spat_cen[:,good] = cen[:,good]
            self.spat_err[:,good] = err[:,good]
            self.spat_msk[:,indx] |= msk[:,indx]
            rmtrace[bad] = True

        # Update the image coordinates
        self.spat_img = np.round(self.spat_cen).astype(int)

        # Erase any previous fitting, PCA, and slit-mask-match results
        # TODO: Should probably get sectioned off as a method because a
        # few methods do this.  Add an option to _reinit_trace_data?
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

    def check_traces(self, cen=None, msk=None, subset=None, min_spatial=None, max_spatial=None,
                     minimum_spec_length=None):
        r"""
        Validate new trace data to be added.

        Steps are:
            - Remove duplicates based on the provided matching
              tolerance (`match_tol` in :attr:`par`). Beware that
              this should be done by providing `cen` directly.
            - Remove traces that do not cover at least some fraction
              of the detector (see `minimum_spec_length`).
            - Remove traces that are at a minimum or maximum spatial
              position (column) (typically the edge of the detector;
              see `min_spatial` and `max_spatial`).

        The only used parameter from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) is `match_tol`.

        .. warning::
            `msk` is edited in-place

        Args:
            cen (`numpy.ndarray`_, optional):
                The adjusted center of the refined traces. Shape is
                :math:`(N_{\rm spec}, N_{\rm refine},)`. If None, use
                :attr:`spat_cen`.
            msk (`numpy.ndarray`_, optional):
                The mask bits for the adjusted center of the refined
                traces. Shape is :math:`(N_{\rm spec}, N_{\rm
                refine},)`. If None, use :attr:`spat_msk`. This is
                edited in-place!
            subset (`numpy.ndarray`_, optional):
                Boolean array selecting the traces to compare. Shape
                is :math:`(N_{\rm trace},)`, with :math:`N_{\rm
                refine}` True values. It is expected that all the
                traces selected by a subset must be from the same
                slit side (left or right). If None, no repeat traces
                can be identified. 
            min_spatial (:obj:`int`, optional):
                Clip traces that hit this minimum spatial index
                (column) at the center spectral row
                (`self.nspec//2`). If None, no traces clipped.
            max_spatial (:obj:`int`, optional):
                Clip traces that hit this maximum spatial index
                (column) at the center spectral row
                (`self.nspec//2`). If None, no traces clipped.
            minimum_spec_length (:obj:`float`, optional):
                The minimum number of spectral rows in an edge trace.

        Returns:
            Returns two boolean arrays selecting the good and bad
            traces.  Shapes are :math:`(N_{\rm trace},)`.
        """
        if self.is_empty:
            msgs.warn('No traces to check.')

        # Indices of traces to check
        indx = np.ones(self.ntrace, dtype=bool) if subset is None else subset
        # Data to compare
        _cen = self.spat_cen if cen is None else cen
        _msk = self.spat_msk if msk is None else msk
        # Ignore inserted traces when flagging
        _bpm = self.bitmask.flagged(_msk, self.bitmask.bad_flags)

        # Find repeat traces; comparison of traces must include
        # unmasked trace data, be traces of the same edge (left or
        # right), and be within the provided matching tolerance
        repeat = np.zeros_like(indx, dtype=bool)
        if subset is not None:
            msgs.info('Tolerance for finding repeat traces: {0:.1f}'.format(self.par['match_tol']))
            side = -1 if np.all(self.traceid[indx] < 0) else 1
            compare = (side*self.traceid > 0) & np.invert(indx)
            if np.any(compare):
                # Use masked arrays to ease exclusion of masked data
                _col = np.ma.MaskedArray(np.round(_cen).astype(int), mask=_bpm)
                spat_bpm = self.bitmask.flagged(self.spat_msk, flag=self.bitmask.bad_flags)
                spat_img = np.ma.MaskedArray(self.spat_img, mask=spat_bpm)
                mindiff = np.ma.amin(np.absolute(_col[:,indx,None]-spat_img[:,None,compare]),
                                    axis=(0,2))
                # TODO: This tolerance uses the integer image
                # coordinates, not the floating-point centroid
                # coordinates....
                repeat[indx] =  (mindiff.data < self.par['match_tol']) & np.invert(mindiff.mask)
                if np.any(repeat):
                    _msk[:,repeat] = self.bitmask.turn_on(_msk[:,repeat], 'DUPLICATE')
                    msgs.info('Found {0} repeat trace(s).'.format(np.sum(repeat)))

        # Find spectrally short traces
        short = np.zeros_like(indx, dtype=bool)
        if minimum_spec_length is not None:
            msgs.info('Minimum spectral length of any trace (pixels): {0:.2f}'.format(
                      minimum_spec_length))
            short[indx] = np.sum(np.invert(_bpm[:,indx]), axis=0) < minimum_spec_length
            if np.any(short):
                _msk[:,short] = self.bitmask.turn_on(_msk[:,short], 'SHORTRANGE')
                msgs.info('Found {0} short trace(s).'.format(np.sum(short)))

        # Find traces that are at the minimum column at the center
        # spectral row
        # TODO: Why only the center row?
        col = np.round(_cen).astype(int)
        hit_min = np.zeros_like(indx, dtype=bool)
        if min_spatial is not None:
            hit_min[indx] = (col[self.nspec//2,indx] <= min_spatial) \
                                & np.invert(_bpm[self.nspec//2,indx])
            if np.any(hit_min):
                _msk[:,hit_min] = self.bitmask.turn_on(_msk[:,hit_min], 'HITMIN')
                msgs.info('{0} trace(s) hit the minimum centroid value.'.format(np.sum(hit_min)))
         
        # Find traces that are at the maximum column at the center
        # spectral row
        # TODO: Why only the center row?
        hit_max = np.zeros_like(indx, dtype=bool)
        if max_spatial is not None:
            hit_max[indx] = (col[self.nspec//2,indx] >= max_spatial) \
                                & np.invert(_bpm[self.nspec//2,indx])
            if np.any(hit_max):
                _msk[:,hit_max] = self.bitmask.turn_on(_msk[:,hit_max], 'HITMAX')
                msgs.info('{0} trace(s) hit the maximum centroid value.'.format(np.sum(hit_max)))

        # Find traces, or trace regions, that fall off the detector
        off_detector = np.zeros_like(indx, dtype=bool)
        if self.par['det_buffer'] > 0:
            pix_off_detector = np.zeros(_cen.shape, dtype=bool)
            pix_off_detector[:,indx] = (_cen[:,indx] < self.par['det_buffer']) \
                                            | (_cen[:,indx] > self.nspat-self.par['det_buffer'])
            off_detector[indx] = np.all(pix_off_detector[:,indx], axis=0)
            _msk[pix_off_detector] = self.bitmask.turn_on(_msk[pix_off_detector], 'OFFDETECTOR')

        # TODO: Check that traces (trace fits?) don't cross?

        # Good traces
        bad = indx & (repeat | short | hit_min | hit_max | off_detector)
        msgs.info('Identified {0} bad trace(s) in all.'.format(np.sum(bad)))
        good = indx & np.invert(bad)
        return good, bad

    def merge_traces(self, merge_frac=0.5, refit=True, debug=False):
        """
        Merge traces based on their spatial separation.

        Traces are merged if:
            - they are from the same slit side (left or right)
            - the fitted trace locations are separated by less than
              `match_tol` (from :attr:`par`) for more than
              `merge_frac` of the spectral range of the detector.

        Since the merging is based on the *fitted* trace location,
        the traces must have been fit beforehand; see
        :func:fit_refine`.

        When there are traces found that should be merged, the
        unmasked centroid measurements from the shorter trace(s) are
        added to the longest trace. The traces are then resorted and
        given new trace IDs.

        .. warning::
            - This automatically removes any existing PCA
              decomposition
            - If functional fits exist and `refit` is False, *all*
              fits are removed.
            - If traces are merged, a warning is issued that the
              traces may no longer be left-right synchronized.

        Args:
            merge_frac (:obj:`float`, optional):
                The fraction of the spectral length of the detector
                that should be less than `match_tol` in :attr:`par`
                used to find traces to merge.
            refit (:obj:`bool`, optional):
                If fitted models exist and traces are merged, run
                :func:`fit_refine` with its default arguments after
                the traces have been merged to refit the trace data.
            debug (:obj:`bool`, optional):
                Only passed to :func:`fit_refine` if `refit` is True.
        """
        if self.is_empty:
            msgs.warn('No traces to merge.')
        if self.spat_fit is None:
            msgs.error('Trace merging requires model fits to the trace location; run fit_refine.')
        _refit = refit
        if refit and self.spat_fit is None:
            msgs.warn('No previous fits existed, so fitting will not be redone.')
            _refit = False

        # Construct the bad pixel mask depending whether we matching
        # models or measurements
        bpm = np.tile(self.fully_masked_traces(flag=self.bitmask.bad_flags),(self.nspec,1))

        # The center values used to test for matches; can either be the
        # modeled or measured data.
        cen_match = np.ma.MaskedArray(self.spat_fit, mask=bpm)
        # The center values used to construct the merged trace; always
        # the measured values
        cen_merge = np.ma.MaskedArray(self.spat_cen, mask=self.spat_msk.astype(bool))
        rmtrace = np.zeros(self.ntrace, dtype='bool')

        for side in ['left', 'right']:
            # Match traces on the same side and 
            indx = np.where((self.is_left if side == 'left' else self.is_right)
                                & np.invert(np.all(cen_match.mask, axis=0)))[0]
            if indx.size == 0:
                continue

            for i in range(indx.size - 1):
                if rmtrace[indx[i]]:
                    continue
                merge = np.ma.sum(np.absolute(cen_match[:,indx[i],None] - cen_match[:,indx[i+1:]])
                                     < self.par['match_tol'], axis=0).filled(0) \
                            > self.nspec*merge_frac
                if not np.any(merge):
                    continue
                rmtrace[indx[i+1:]] = merge
                merge = np.append(indx[i], indx[i+1:][merge])
                msgs.info('Merging traces: {0}'.format(self.traceid[merge]))
                merged_trace = np.ma.mean(cen_merge[:,merge], axis=1)
                gpm = np.invert(np.ma.getmaskarray(merged_trace))
                self.spat_cen[gpm,indx[i]] = merged_trace[gpm]
                self.spat_msk[gpm,indx[i]] = 0

        if not np.any(rmtrace):
            msgs.info('No traces merged.')
            return

        # Remove traces and resort them
        self.remove_traces(rmtrace)

        # Erase any previous fitting, PCA, and slit-mask-match results
        # TODO: Should probably get sectioned off as a method because a
        # few methods do this.
        self.spat_fit_type = None
        self.spat_fit = None
        self.pca_type = None
        self.pca = None
        self.design = None
        self.objects = None

        if _refit:
            self.fit_refine(debug=debug)

    @property
    def is_left(self):
        """Boolean array selecting the left traces."""
        if self.is_empty:
            msgs.error('No traces!')
        return self.traceid < 0
    
    @property
    def is_right(self):
        """Boolean array selecting the right traces."""
        if self.is_empty:
            msgs.error('No traces!')
        return self.traceid > 0

    @property
    def is_empty(self):
        """Flag that object has no trace data."""
        return self.traceid is None
        
    @property
    def is_synced(self):
        """
        Confirm slit edges are synced.
        """
        if self.is_empty:
            return False
        # Do it
        gpm = np.invert(self.fully_masked_traces(flag=self.bitmask.bad_flags,
                                                 exclude=self.bitmask.exclude_flags))
        side = np.clip(self.traceid[gpm], -1, 1)
        return side[0] == -1 and side.size % 2 == 0 and np.all(side[1:] + side[:-1] == 0)

    def check_synced(self, rebuild_pca=False):
        """
        Quality check and masking of the synchronized edges.

        Before executing this method, the slit edges must be
        synchronized (see :func:`sync`) and ordered spatially in
        left-right pairs (see :func:`spatial_sort`). The former is
        checked explicitly. Any traces fully masked as bad (see
        :func:`clean_traces`) are removed, along with its
        synchronized partner.

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are
        `minimum_slit_gap`, `minimum_slit_length`, and
        `length_range`.

        Checks are:
            - Any trace falling off the edge of the detector is
              masked (see :class:`EdgeTraceBitMask`). This is the
              only check performed by default (i.e., when no keyword
              arguments are provided).
            - Traces that form slit gaps (the median difference
              between the right and left traces of adjacent slits)
              that are below an absolute tolerance are removed and
              the two relevant slits are merged. This is done before
              the checks of the slit length below such that the
              merged slit is assessed in any expected slit length
              constraints.
            - Traces that form a slit with a length (the median
              difference between the left and right edges) below an
              absolute tolerance (i.e., `right-left < atol`) are
              masked. The absolute tolerance is set using the
              platescale provided by the spectrograph class, the
              spatial binning (from :attr:`binning`), and the minimum
              slit length in arcsec (`minimum_slit_length` in
              :attr:`par`).
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
        if self.is_empty:
            msgs.warn('No traces to check.')

        # Remove any fully masked traces and its synced counterpart;
        # force the removal of traces marked as SYNCERROR, even if
        # those traces were inserted.
        self.clean_traces(sync_rm='both', force_flag='SYNCERROR')

        # Use the fit data if available
        trace_cen = self.spat_cen if self.spat_fit is None else self.spat_fit

        # Flag trace locations falling off the detector. This general
        # check that is independent of whether or not the traces are
        # synced. However, this could mask traces based on the *fitted*
        # position, instead of the measured centroid position.
        # TODO: Keep a separate mask for the fit positions?
        indx = (trace_cen < 0) | (trace_cen >= self.nspat)
        if np.any(indx):
            new_masks = True
            self.spat_msk[indx] = self.bitmask.turn_on(self.spat_msk[indx], 'OFFDETECTOR')

        # Check the slits are synced
        if not self.is_synced:
            msgs.error('Edge traces are not yet (or improperly) synced.  Either sync() failed '
                       'or has not yet been executed.')

        # Parse parameters and report
        gap_atol = None
        length_atol = None
        length_rtol = self.par['length_range']
        if self.par['minimum_slit_length'] is not None or self.par['minimum_slit_gap'] is not None:
            platescale = parse.parse_binning(self.binning)[1] \
                            * self.spectrograph.detector[self.det-1]['platescale']
            msgs.info('Binning: {0}'.format(self.binning))
            msgs.info('Platescale per binned pixel: {0}'.format(platescale))
            if self.par['minimum_slit_length'] is not None:
                length_atol = self.par['minimum_slit_length']/platescale
            if self.par['minimum_slit_gap'] is not None:
                gap_atol = self.par['minimum_slit_gap']/platescale

        msgs.info('Minimum slit gap (binned pixels): {0}'.format(gap_atol))
        msgs.info('Minimum slit length (binned pixels): {0}'.format(length_atol))
        msgs.info('Range in slit length not limited' if length_rtol is None else
                  'Range in slit length limited to +/-{0:.1f}%'.format(length_rtol*100))

        # TODO: Should here and below only use the unmasked parts of
        # the trace for the slit gap and length computations?

        # Keep track of whether or not any new masks were applied
        new_masks = False

        if gap_atol is not None:
            # Calculate the slit gaps
            slit_gap = np.median(np.diff(trace_cen[:,1:], axis=1)[:,::2], axis=0)
            indx = slit_gap < gap_atol
            if np.any(indx):
                # TODO: Allow for these traces to be flagged instead of just removed?
                msgs.info('Found {0} slit(s) with gaps below {1} arcsec ({2:.2f} pixels).'.format(
                            np.sum(indx), self.par['minimum_slit_gap'], gap_atol))
                rmtrace = np.concatenate(([False],np.repeat(indx,2),[False]))
                self.remove_traces(rmtrace, rebuild_pca=rebuild_pca)
                # TODO: This should never happen, but keep this around
                # until we're sure it doesn't.
                if self.is_empty:
                    msgs.error('Coding error: Removing gaps removed all traces.')
                # Reset the trace center data to use
                trace_cen = self.spat_cen if self.spat_fit is None else self.spat_fit

        # Calculate the slit length and gap
        slit_length = np.median(np.squeeze(np.diff(trace_cen.reshape(self.nspec,-1,2), axis=-1)),
                                axis=0)
        if length_atol is not None:
            # Find any short slits (flag both edges of the slit)
            indx = np.repeat(slit_length < length_atol, 2)
            if np.sum(indx) == self.ntrace:
                if self.ntrace == 2:
                    # TODO: This is a kludge that was necessary to
                    # address the features seen in the Keck_LRIS_blue
                    # long-slit data in the dev suite. How often does
                    # something like this happen?
                    # TODO: All of this should get pulled out into a
                    # corner-case method.
                    msgs.warn('The single slit found has been rejected because it is too short.'
                              '  If this was by mistake, re-run pypeit with a smaller '
                              '`minimum_slit_length` parameter.  Otherwise, we assume this is a '
                              'long-slit with one edge off the detector and with the current '
                              'slit edges errantly isolating some feature in the data.')
                    # TODO: May want to limit the number of columns included in this calculation.
                    if np.mean(self.img[:,int(np.ceil(np.max(trace_cen[:,1]))):]) \
                            > np.mean(self.img[:,:int(np.floor(np.min(trace_cen[:,0])))]):
                        msgs.warn('The mean of the trace image to the right of the right trace '
                                  'is larger than it is to the left of the left trace. Removing '
                                  'the right trace and re-synchronizing.')
                        self.remove_traces(np.array([False,True]))
                    else:
                        msgs.warn('The mean of the trace image to the left of the left trace '
                                  'is larger than it is to the right of the right trace. Removing '
                                  'the right trace and re-synchronizing.')
                        self.remove_traces(np.array([True,False]))
                    # TODO: I *really* don't like this because it has
                    # the potential to yield an infinite loop, but it's
                    # also the simplest approach.
                    self.sync(rebuild_pca=False)
                    return
                msgs.warn('All slits are too short!')
            if np.any(indx):
                new_masks = True
                msgs.info('Rejecting {0} slits that are too short.'.format(np.sum(indx)))
                self.spat_msk[:,indx] = self.bitmask.turn_on(self.spat_msk[:,indx], 'SHORTSLIT')


        if length_rtol is not None:
            # Find slits that are not within the provided fraction of
            # the median length
            indx = np.repeat(np.absolute(np.log(slit_length/np.median(slit_length)))
                             > np.log(1+length_rtol), 2)
            if np.any(indx):
                new_masks = True
                msgs.info('Rejecting {0} abnormally long or short slits.'.format(np.sum(indx)))
                self.spat_msk[:,indx] = self.bitmask.turn_on(self.spat_msk[:,indx], 'ABNORMALSLIT')

        # TODO: Check that slit edges meet a minimum slit gap?

        if self.par['sync_clip']:
            # Remove traces that have been fully flagged as bad
            rmtrace = self.fully_masked_traces(flag=self.bitmask.bad_flags)
            self.remove_traces(rmtrace, rebuild_pca=rebuild_pca, sync_rm='both')
            if self.is_empty:
                msgs.warn('Assuming a single long-slit and continuing.')
                self.bound_detector()
        elif new_masks:
            # Reset the PCA if new masks are applied
            self._reset_pca(rebuild_pca and self.pca is not None and self.can_pca())

    def rm_user_traces(self, rm_traces):
        """
        Parse the user input traces to remove

        Args:
            rm_user_traces (list):
              y_spec, x_spat pairs

        Returns:

        """
        if not self.is_synced:
            msgs.error("This method should not be run until after the slits are synced")
        # Setup
        lefts = self.spat_fit[:, self.is_left]
        rights = self.spat_fit[:, self.is_right]
        indx = np.zeros(self.ntrace, dtype=bool)
        # Loop me
        for rm_trace in rm_traces:
            # Deconstruct
            y_spec, xcen = rm_trace
            #
            lefty = lefts[y_spec,:]
            righty = rights[y_spec,:]
            # Match?
            bad_slit = (lefty < xcen) & (righty > xcen)
            if np.any(bad_slit):
                # Double check
                if np.sum(bad_slit) != 1:
                    msgs.error("Something went horribly wrong in edge tracing")
                #
                idx = np.where(bad_slit)[0][0]
                indx[2*idx:2*idx+2] = True
                msgs.info("Removing user-supplied slit at {},{}".format(xcen, y_spec))
        # Remove
        self.remove_traces(indx, sync_rm='both')

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

                    - ``'ignore'``: Just remove the selected traces and
                      ignore the synchronization.
                    - ``'both'``: If at least one of the traces in a
                      pair is selected for removal, remove both.
                    - ``'neither'``: If only one of the traces in a pair
                      is selected for removal, remove neither.

        """
        # Make sure there are traces to remove
        if not np.any(indx):
            msgs.warn('No trace to remove.')
            return
        # Check the input
        _indx = np.atleast_1d(indx).copy()
        if _indx.size != self.ntrace:
            msgs.error('Boolean array selecting traces to remove has incorrect length.')

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
                msgs.error('Unknown sync removal keyword: {0}'.format(sync_rm))

        if np.all(_indx):
            msgs.warn('All traces removed!')
            self._reinit_trace_data()
            return
            
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
        self._reset_pca(rebuild_pca and self.pca is not None and self.can_pca())

    def clean_traces(self, sync_rm='ignore', force_flag=None):
        """
        Remove any traces that are fully masked as bad.

        For removal, the traces must be fully masked and
        *none* of its pixels can be masked by one of the insertion
        flags; see :func:`fully_masked_traces`. 

        To force removal of traces with certain flags, regardless of
        their insertion status, use `force_flag`. See
        :class:`EdgeTraceBitMask` for list of flags.

        Args:
            sync_rm (:obj:`str`, optional):
                If the traces are left-right synchronized (see
                :func:`sync` and :func:`is_synced`), use this method
                to deal with edges paired with those to be removed.
                See :func:`remove_traces`.
            force_flag (:obj:`str`, :obj:`list`, optional):
                Force removal of traces fully masked with these
                flags, even if they are also flagged as having been
                inserted.
        """
        if self.is_empty:
            msgs.warn('No traces to clean.')
            return

        # Traces to remove
        rmtrace = self.fully_masked_traces(flag=self.bitmask.bad_flags,
                                           exclude=self.bitmask.insert_flags)
        if force_flag is not None:
            rmtrace |= self.fully_masked_traces(flag=force_flag)

        if np.any(rmtrace):
            # The removed traces should not have been included in the
            # PCA decomposition to begin with, but this call to remove
            # traces has to "rebuild" the PCA because it will remove it
            # otherwise.
            # TODO: Fix remove_traces to check if traces used in the
            # construction of the PCA, and only then use the
            # rebuild_pca flag to decide if it should rebuild.
            self.remove_traces(rmtrace, rebuild_pca=True, sync_rm=sync_rm)

    def spatial_sort(self, use_mean=False, use_fit=True):
        """
        Sort the traces spatially.

        The coordinates used for the sorting is either the measured
        centroids or the fitted parameterization (see `use_fit`). The
        fiducial coordinates that are sorted are either the mean of
        the unmasked coordinates over all spectral rows or the
        unmasked coordinates at a specified reference spectral row
        (see `use_mean`).

        The trace IDs are also reassigned to be sorted spatially;
        i.e., the trace IDs for three synced slits would be `[-1, 1,
        -2, 2, -3, 3]`.
        
        All attributes are edited in-place.

        Args:
            use_mean (:obj:`bool`, optional):
                Sort according to the mean of the masked spatial
                positions. If False, the spatial position at a
                reference spectral row is used, where the reference
                spectral row is either the same as used by the PCA
                (if available) or the result of
                :func:`pypeit.core.trace.most_common_trace_row` using
                the current trace mask.
            use_fit (:obj:`bool`, optional):
                Sort according to the fit positions instead of the
                measured positions. Otherwise, only use the fit
                positions if they're available and the measured
                location is masked.
        """
        if self.is_empty:
            msgs.error('No traces to sort.')

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
                                else (self.pca[0].reference_row if self.par['left_right_pca']
                                    else self.pca.reference_row)
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
        if self.is_empty:
            return edge_img
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
        # Check that there are traces to refine!
        if self.is_empty:
            msgs.error('No traces to refine!')

        # Parse parameters and report
        maxshift = self.par['max_shift_abs']
        maxerror = self.par['max_spat_error']
        function = self.par['fit_function']
        order = self.par['fit_order']
        fwhm = self.par['fwhm_uniform'] if weighting == 'uniform' else self.par['fwhm_gaussian']
        maxdev = self.par['fit_maxdev']
        maxiter = self.par['fit_maxiter']
        niter = self.par['fit_niter']
        minimum_spec_length = self.par['fit_min_spec_length']*self.nspec
        xmin = 0.
        xmax = self.nspec-1.

        msgs.info('-'*50)
        msgs.info('{0:^50}'.format('Fitting Polynomial to Edge Trace'))
        msgs.info('-'*50)
        msgs.info('Max shift btwn input and remeasured edge centroids: {0:.2f}'.format(maxshift))
        msgs.info('Max centroid error: {0}'.format(maxerror))
        msgs.info('Trace fitting function: {0}'.format(function))
        msgs.info('Trace fitting order: {0}'.format(order))
        msgs.info('Weighting for remeasuring edge centroids: {0}'.format(weighting))
        msgs.info('FWHM parameter for remeasuring edge centroids: {0:.1f}'.format(fwhm))
        msgs.info('Maximum deviation for fitted data: {0:.1f}'.format(maxdev))
        msgs.info('Maximum number of rejection iterations: {0}'.format(maxiter))
        msgs.info('Number of remeasuring and refitting iterations: {0}'.format(niter))

#        masked_cen = np.ma.MaskedArray(self.spat_cen, mask=self.spat_msk.astype(bool))
#        indx = np.argmax(np.sum(np.invert(masked_cen.mask),axis=0))
#        diff_cen = masked_cen - masked_cen[:,indx,None]
#        meanoffset = np.ma.mean(diff_cen, axis=0)
#        offset_cen = masked_cen - meanoffset[None,:]
#        x = np.arange(self.nspec)
#        for i in range(self.ntrace):
#            plt.scatter(x, offset_cen[:,i])
#        plt.show()

        # Check the traces to make sure they meet the minimum length.
        # This modifies self.spat_msk directly.
        self.check_traces(minimum_spec_length=minimum_spec_length)

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

        # Flag bad traces; this explicitly does *not* exclude inserted traces
        spat_bpm = self.bitmask.flagged(self.spat_msk, flag=self.bitmask.bad_flags)

        # Fit both sides
        for side in ['left', 'right']:
            # Get the image relevant to tracing this side
            _sobel_sig = self._side_dependent_sobel(side)
            # Select traces on this side and that are not fully masked
            indx = (self.is_left if side == 'left' else self.is_right) \
                        & np.invert(np.all(spat_bpm, axis=0))
            if not np.any(indx):
                continue

            # Perform the fit
            fit[:,indx], cen[:,indx], err[:,indx], msk[:,indx], _ \
                    = trace.fit_trace(_sobel_sig, self.spat_cen[:,indx], order, ivar=ivar,
                                      bpm=bpm, trace_bpm=spat_bpm[:,indx],
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

    def can_pca(self):
        """
        Determine if traces are suitable for PCA decomposition.

        The criterion is that a minimum number of traces
        (``pca_min_edges``) must cover more than the fraction of the
        full spectral range specified by `fit_min_spec_length` in
        :attr:`par`. Traces that are inserted are ignored.

        If the PCA decomposition will be performed on the left and
        right traces separately, the function will return `False` if
        there are fewer than the minimum left *or* right edge traces.

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are
        ``fit_min_spec_length``, ``left_right_pca``, and ``pca_min_edges``.

        .. warning::
            This function calls :func:`check_trace` using
            `fit_min_spec_length` to flag short traces, meaning that
            :attr:`spat_msk` will can be altered by this call.

        Returns:
            :obj:`bool`: Flag that traces meet criterion for PCA
            decomposition.
        """
        # Check that there are traces to refine!
        if self.is_empty:
            return False

        # Set and report the minimum length needed for the PCA in
        # pixels
        minimum_spec_length = self.par['fit_min_spec_length']*self.nspec
        msgs.info('Minium length of traces to include in the PCA: {0}'.format(minimum_spec_length))

        # This call to check_traces will flag any trace with a length
        # below minimum_spec_length as SHORTRANGE
        # TODO: This means that SHORTRANGE can actually have two
        # different meanings. Instead use one bit for "too short for
        # detection" and a separate one for "too short to fit"?
        self.check_traces(minimum_spec_length=minimum_spec_length)

        # Find the valid traces
        # NOTE: Because of the run of check_traces above, short traces
        # are fully flagged meaning that we can just check if the
        # length of the trace is larger than 0.
        good = np.sum(np.invert(self.bitmask.flagged(self.spat_msk)), axis=0) > 0

        # Returned value depends on whether or not the left and right
        # traces are done separately
        return np.sum(good[self.is_left]) > self.par['pca_min_edges'] \
                    and np.sum(good[self.is_right]) > self.par['pca_min_edges'] \
                    if self.par['left_right_pca'] else np.sum(good) > self.par['pca_min_edges']

    def predict_traces(self, spat_cen, side=None):
        """
        Use the PCA decomposition to predict traces.

        The PCA decomposition must be available via :attr:`pca`; see
        :func:`build_pca`. This is a convenience method to handle the
        PCA predictions given that left and right traces can be
        decomposed separately or simultaneously.

        Args:
            spat_cen (:obj:`float`, `numpy.ndarray`):
                A single value or 1D array with the spatial location
                (column) for 1 or more traces to predict. The
                predicted traces will pass through these spatial
                positions (columns) and the reference spectral row
                set for the PCA decomposition; see :func:`build_pca`.
            side (:obj:`float`, `numpy.ndarray`, optional):
                A single value or 1D integer array indicating the
                edge side to be predicted; -1 for left and 1 for
                right. Must be the same length as `spat_cen`. This is
                only used if the PCA is side-dependent, and *must* be
                provided in the case that it is (see `left_right_pca`
                in :attr:`par`).

        Returns:
            `numpy.ndarray`: A 1D or 2D array of size :attr:`nspec`
            by the length of the position array provided. If a single
            coordinate is provided, a single trace vector is
            returned.
        """
        if self.pca is None:
            msgs.error('Must first run the PCA analysis fo the traces; run build_pca.')

        _spat_cen = np.atleast_1d(spat_cen)
        _side = np.atleast_1d(side)
        if _spat_cen.size != _side.size:
            msgs.error('Spatial locations and side integers must have the same shape.')

        if self.par['left_right_pca']:
            trace_add = np.zeros((self.nspec,0), dtype='float')
            for s,p in zip([-1,1], self.pca):
                indx = _side == s
                if not np.any(indx):
                    continue
                trace_add = np.hstack((trace_add, p.predict(np.atleast_1d(_spat_cen[indx]))))
        else:
            trace_add = self.pca.predict(_spat_cen)

        return trace_add if isinstance(spat_cen, np.ndarray) else trace_add.ravel()

    def build_pca(self, use_center=False, debug=False):
        """
        Build a PCA model of the current edge data.

        Primarily a wrapper that instantiates :attr:`pca`. If left
        and right traces are analyzed separately, :attr:`pca` will be
        a list of two :class:`pypeit.tracepca.TracePCA` objects;
        otherwise :attr:`pca` is a single
        :class:`pypeit.tracepca.TracePCA` object. After executing
        this, traces can be predicted by calling
        `self.pca.predict(spat)` (see
        :func:`pypeit.tracepca.TracePCA.predict`) if both sides are
        analyzed simultaneously; otherwise, left and tright traces
        can be predicted using `self.pca[0].predict(spat)`. and
        `self.pca[1].predict(spat)`, respectively. See
        :func:`predict_traces`.

        If no parametrized function has been fit to the trace data or
        if specifically requested (see `use_center`), the PCA is
        based on the measured trace centroids (:attr:`spat_cen`);
        othwerwise, the PCA uses the parametrized trace fits
        (:attr:`spat_fit`).

        The reference spectral row used for the decomposition (see
        :class:`pypeit.tracepca.TracePCA`) is set by
        :func:`pypeit.core.trace.most_common_trace_row` using the
        existing mask. If treating left and right traces separately,
        the reference spectral row is the same for both PCA
        decompositions.

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are
        `fit_min_spec_length`, `left_right_pca`, `pca_n`,
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
        if self.is_empty:
            msgs.error('No traces exist.')

        # Parse parameters and report
        left_right_pca = self.par['left_right_pca']
        npca = self.par['pca_n']
        pca_explained_var = self.par['pca_var_percent']
        function = self.par['pca_function']
        order = self.par['pca_order']
        lower, upper = self.par['pca_sigrej'] if hasattr(self.par['pca_sigrej'], '__len__') \
                        else (self.par['pca_sigrej'],)*2
        maxrej = self.par['pca_maxrej']
        maxiter = self.par['pca_maxiter']

        msgs.info('-'*50)
        msgs.info('{0:^50}'.format('Constructing PCA interpolator'))
        msgs.info('-'*50)
        msgs.info('PCA composition of the left and right traces is done {0}.'.format(
                    'separately' if left_right_pca else 'simultaneously'))
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

        # Check if the PCA decomposition can be performed
        if not self.can_pca():
            msgs.error('Traces do not meet necessary criteria for the PCA decomposition.')

        # Set the data used to construct the PCA
        self.pca_type = 'center' if self.spat_fit is None or use_center else 'fit'

        # When constructing the PCA, ignore bad trace measurements
        # *and* any traces inserted by hand.
        bpm = self.bitmask.flagged(self.spat_msk)

        # TODO: Is there a way to propagate the mask to the PCA?
        # TODO: Keep a separate mask specifically for the fit data? e.g., spat_fit_msk

        # The call to can_pca means that short traces are fully masked
        # and that valid traces will be any trace with unmasked pixels.
        use_trace = np.sum(np.invert(bpm), axis=0) > 0

        # Set the reference row so that, regardless of whether the PCA
        # is for the left, right, or all traces, the reference row is
        # always the same.
        reference_row = trace.most_common_trace_row(bpm[:,use_trace])
#        reference_row = self.npec//2

        # Setup the list of traces to use in a single object so that we
        # can loop though them, regardless of whether we're performing
        # the PCA for left and right traces separately
        pcaindx = [use_trace]
        if left_right_pca:
            pcaindx = [None, None]
            for i, side in enumerate(['left', 'right']):
                pcaindx[i] = (self.is_left if side == 'left' else self.is_right) & use_trace
                msgs.info('Using {0}/{1} of the {2} traces in the PCA analysis.'.format(
                                np.sum(pcaindx[i]), self.ntrace, side))

        # Run the PCA decomposition and construct its interpolator
        self.pca = [None]*len(pcaindx)
        for i,indx in enumerate(pcaindx):
            # Grab the trace data. This uses all the data, even if some
            # of it is masked.
            trace_inp = self.spat_cen[:,indx] if self.spat_fit is None or use_center \
                            else self.spat_fit[:,indx]

            # Instantiate the PCA
            self.pca[i] = TracePCA(trace_inp, npca=npca, pca_explained_var=pca_explained_var,
                                   reference_row=reference_row)

            # Set the order of the function fit to the PCA
            # coefficiencts: Order is set to cascade down to lower
            # order for components that account for a smaller
            # percentage of the variance.
            _order = np.clip(order - np.arange(self.pca[i].npca), 1, None).astype(int)
            msgs.info('Order of function fit to each component: {0}'.format(_order))

            # Apply a 10% relative error to each coefficient. This
            # performs better than use_mad, since larger coefficients
            # will always be considered inliers, if the coefficients
            # vary rapidly with order as they sometimes do.
            #ivar = utils.inverse(np.square(np.fmax(0.1*np.abs(self.pca[i].pca_coeffs), 0.1)))
            #ivar = None

            # TODO: Instead, weight by the mean/median value of
            # sobel_sig along each trace.

            # Run the fit
            self.pca[i].build_interpolator(_order, function=function, lower=lower,
                                           upper=upper, minx=0., maxx=self.nspat-1., maxrej=maxrej,
                                           maxiter=maxiter, debug=debug)

            # TODO: Use the rejected pca coefficiencts to reject traces?

        # If left and right use the same PCA, get rid of the list
        if not left_right_pca:
            self.pca = self.pca[0]

    def pca_refine(self, use_center=False, debug=False, force=False):
        """
        Use a PCA decomposition to refine the traces.

        If no parametrized function has been fit to the trace data or
        if specifically requested (see `use_center`), the PCA is
        based on the measured trace centroids; othwerwise, the PCA
        uses the parametrized trace fits.

        If needed or forced to, this first executes :func:`build_pca`
        and then uses :func:`predict_traces` to use the PCA to reset
        the trace data.

        Only used parameter from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) is
        `left_right_pca`.

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
        if self.is_empty:
            msgs.error('No traces to refine!')

        # Perform the PCA decomposition if necessary
        _pca_type = 'center' if use_center or self.spat_fit is None else 'fit'
        if force or self.pca is None or self.pca_type != _pca_type:
            self.build_pca(use_center=use_center, debug=debug)

        self.spat_fit_type = 'pca'

        # Predict the traces using the PCA
        reference_row = self.pca[0].reference_row if self.par['left_right_pca'] \
                            else self.pca.reference_row
        trace_ref = self.spat_cen[reference_row,:] if self.pca_type == 'center' \
                            else self.spat_fit[reference_row,:]
        side = self.is_right.astype(int)*2-1
        self.spat_fit = self.predict_traces(trace_ref, side)

        # TODO: Compare with the fit data. Remove traces where the mean
        # offset between the PCA prediction and the measured centroids
        # are larger than some threshold?

        # Log what was done
        self.log += [inspect.stack()[0][3]]

    def peak_refine(self, rebuild_pca=False, debug=False):
        """
        Refine the trace by isolating peaks and troughs in the
        Sobel-filtered image.

        This function *requires* that the PCA model exists; see
        :func:`build_pca` or :func:`pca_refine`. It is also primarily
        a wrapper for :func:`pypeit.core.trace.peak_trace`.

        If the left and right traces have separate PCA
        decompositions, this function makes one call to
        :func:`pypeit.core.trace.peak_trace` for each side.
        Otherwise, a single call is made to
        :func:`pypeit.core.trace.peak_trace` where both the peak and
        troughs in :attr:`sobel_sig` are detected and traced.

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are
        `left_right_pca`, `edge_thresh`, `smash_range`,
        `edge_detect_clip`, `trace_median_frac`, `trace_thresh`,
        `fit_function`, `fit_order`, `fwhm_uniform`, `fwhm_uniform`,
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
        # Check that there are traces to refine!
        if self.is_empty:
            msgs.error('No traces are defined.')

        if self.pca is None:
            msgs.error('Must first run the PCA analysis fo the traces; run build_pca.')

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

        msgs.info('-'*50)
        msgs.info('{0:^50}'.format('Refining traces using collapsed Sobel image'))
        msgs.info('-'*50)
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

        # Generate bogus ivar and mask once here so that they don't
        # have to be generated multiple times.
        # TODO: Keep these as work space as class attributes? so that
        # they don't need to be reinstantiated.
        ivar = np.ones_like(self.sobel_sig, dtype=float)
        bpm = np.zeros_like(self.sobel_sig, dtype=bool) if self.bpm is None else self.bpm

        # Treatment is different if the PCA was done for all traces or
        # separately for left and right traces
        if self.par['left_right_pca']:
            # Initialize the arrays holding the results for both sides
            fit = np.zeros((self.nspec,0), dtype='float')
            cen = np.zeros((self.nspec,0), dtype='float')
            err = np.zeros((self.nspec,0), dtype='float')
            msk = np.zeros((self.nspec,0), dtype=self.bitmask.minimum_dtype())

            # Iterate through each side
            for i,side in enumerate(['left', 'right']):
                # Get the image relevant to tracing
                _sobel_sig = trace.prepare_sobel_for_trace(self.sobel_sig, bpm=self.bpm, boxcar=5, side=side)

                _fit, _cen, _err, _msk, nside \
                        = trace.peak_trace(_sobel_sig, ivar=ivar, bpm=bpm,
                                           trace_map=self.pca[i].predict(np.arange(self.nspat)),
                                           smash_range=smash_range, peak_thresh=peak_thresh,
                                           peak_clip=peak_clip, trace_median_frac=trace_median_frac,
                                           trace_thresh=trace_thresh, fwhm_uniform=fwhm_uniform,
                                           fwhm_gaussian=fwhm_gaussian, function=function,
                                           order=order, maxdev=maxdev, maxiter=maxiter,
                                           niter_uniform=niter_uniform,
                                           niter_gaussian=niter_gaussian, bitmask=self.bitmask,
                                           debug=debug)
                fit = np.hstack((fit,_fit))
                cen = np.hstack((cen,_cen))
                err = np.hstack((err,_err))
                msk = np.hstack((msk,_msk))
                if side == 'left':
                    nleft = nside
        else:
            # Get the image relevant to tracing
            _sobel_sig = trace.prepare_sobel_for_trace(self.sobel_sig, bpm=self.bpm, boxcar=5, side=None)

            # Find and trace both peaks and troughs in the image. The
            # input trace data (`trace` argument) is the PCA prediction
            # of the trace that passes through each spatial position at
            # the reference spectral pixel.
            fit, cen, err, msk, nleft \
                    = trace.peak_trace(_sobel_sig, ivar=ivar, bpm=bpm,
                                       trace_map=self.pca.predict(np.arange(self.nspat)),
                                       smash_range=smash_range, peak_thresh=peak_thresh,
                                       peak_clip=peak_clip, trough=True,
                                       trace_median_frac=trace_median_frac,
                                       trace_thresh=trace_thresh, fwhm_uniform=fwhm_uniform,
                                       fwhm_gaussian=fwhm_gaussian, function=function, order=order,
                                       maxdev=maxdev, maxiter=maxiter, niter_uniform=niter_uniform,
                                       niter_gaussian=niter_gaussian, bitmask=self.bitmask,
                                       debug=debug)

        # Assess the output
        ntrace = fit.shape[1]
        if ntrace < self.ntrace:
            msgs.warn('Found fewer traces using peak finding than originally available.  '
                      'May want to reset peak threshold.')

#        # TODO: Check the traces? NOTE: This currently only identifies short
#        # traces.
#        good, bad = self.check_traces(cen, msk)
#        left = np.zeros(fit.shape[1])

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
        self._reset_pca(rebuild_pca and self.can_pca())
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

    def _get_reference_locations(self, trace_cen, add_edge):
        """
        Insert the reference locations for the traces to add.

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are
        `sync_center`, `sync_to_edge`, and `gap_offset`. See
        :func:`sync`.

        Args:
            trace_cen (`numpy.ndarray`_):
                Trace data to use for determining new edge locations.
            add_edge (`numpy.ndarray`_):
                Boolean array indicating that a trace in the new
                array is an added trace. The number of *False*
                entries in `add_edge` should match the length of the
                2nd axis of `trace_cen`.

        Returns:
            `numpy.ndarray`_: Reference positions for all edge
            traces, both for the existing and new traces.
        """
        # Parse parameters and report
        center_mode = self.par['sync_center']
        to_edge = self.par['sync_to_edge']
        gap_offset = self.par['gap_offset']

        msgs.info('Mode used to set spatial position of new traces: {0}'.format(center_mode))
        msgs.info('For first left and last right, set trace to the edge: {0}'.format(to_edge))
        if center_mode == 'gap':
            msgs.info('Gap offset for adjacent slits: {0}'.format(gap_offset))

        # Get the reference row for the placement calculation; allow
        # the use of inserted traces.
        bpm = self.bitmask.flagged(self.spat_msk, flag=self.bitmask.bad_flags)
        reference_row = trace.most_common_trace_row(bpm) if self.pca is None \
                            else (self.pca[0].reference_row if self.par['left_right_pca']
                                    else self.pca.reference_row)

        # Check that the trace data are sorted at this spectral row
        if not np.array_equal(np.arange(trace_cen.shape[1]),
                              np.argsort(trace_cen[reference_row,:])):
            msgs.error('Trace data must be spatially sorted.')

        # Build a masked array with the trace positions at that
        # spectral row, masked where new traces are supposed to go.
        trace_ref = np.ma.masked_all(add_edge.size)
        trace_ref[np.invert(add_edge)] = trace_cen[reference_row,:]
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
        if center_mode in ['median', 'gap']:
            offset = np.full(nslits, np.ma.median(slit_length), dtype=float)
        elif center_mode == 'nearest':
            # Find the index of the nearest slit with both existing
            # edges (i.e. has an unmasked slit length)
            nearest = utils.nearest_unmasked(slit_center, use_indices=True)
            # The offset is the slit length of the nearest valid slit
            offset = slit_length.data[nearest]
        else:
            msgs.error('Unknown trace centering mode: {0}'.format(center_mode))

        # Set the new edge trace reference locations
        for slit in range(nslits):
            if not slit_length.mask[slit]:
                # Both slit edges already defined
                continue
            if trace_ref.mask[slit,0]:
                # Add the left edge
                if slit > 0 and center_mode == 'gap':
                    trace_ref[slit,0] = trace_ref[slit-1,1] + gap_offset
                else:     
                    trace_ref[slit,0] = 0 if slit == 0 and to_edge \
                                            else trace_ref[slit,1] - offset[slit]
                continue
            # Add the right edge
            if slit < nslits-1 and center_mode == 'gap':
                trace_ref[slit,1] = trace_ref[slit+1,0] - gap_offset
            else:
                trace_ref[slit,1] = self.nspat - 1 if slit == nslits-1 and to_edge \
                                        else trace_ref[slit,0] + offset[slit]

        # TODO: Nothing should now be masked. Get rid of this once
        # satisfied that the coding is correct.
        if np.any(trace_ref.mask):
            msgs.error('Coding error: this should not happen')
        trace_ref = trace_ref.data.ravel()

        # Check that the predicted reference positions don't cause slit
        # overlaps
        indx = np.where(add_edge[1:-1])[0]
        if len(indx) > 0:
            indx += 1
            # Predicted below an existing edge (will never add paired edges)
            too_lo = trace_ref[indx] < trace_ref[indx-1]
            trace_ref[indx[too_lo]] = trace_ref[indx[too_lo]-1] + gap_offset
            noffset = np.sum(too_lo)
            # Predicted above an existing edge
            too_hi = trace_ref[indx] > trace_ref[indx+1]
            trace_ref[indx[too_hi]] = trace_ref[indx[too_hi]+1] - gap_offset
            noffset += np.sum(too_hi)
            if noffset > 0:
                msgs.warn('Reference locations for {0} slit edges adjusted '.format(noffset)
                          + 'to have a slit gap of {0} pixel(s).'.format(gap_offset))

        return trace_ref

    def nudge_traces(self, trace_cen):
        r"""
        Nudge traces away from the detector edge.

        Traces are shifted spatially, up to a maximum value
        (`max_nudge`), to be no closer than a minimum number
        (`det_buffer`) pixels from the detector edges. Both
        parameters are pulled from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`). No limit is
        imposed on the size of the shift if `max_nudge` is None.

        .. warning::
            A failure mode that is not dealt with is when multiple
            traces fall off the detector and are nudged to nearly the
            same location.

        Args:
            trace_cen (`numpy.ndarray`_):
                Array with trace locations to adjust. Must be 2D with
                shape :math:`(N_{\rm spec}, N_{\rm trace})`.

        Returns:
            `numpy.ndarray`_: The nudged traces.
        """
        # Check input
        if self.par['max_nudge'] is not None and self.par['max_nudge'] <= 0:
            # Nothing to do
            return trace_cen
        if trace_cen.shape[0] != self.nspec:
            msgs.error('Traces have incorrect length.')
        _buffer = self.par['det_buffer']
        if _buffer < 0:
            msgs.warn('Buffer must be greater than 0; ignoring.')
            _buffer = 0

        msgs.info('Nudging traces, by at most {0} pixel(s)'.format(self.par['max_nudge'])
                  + ', to be no closer than {0} pixel(s) from the detector edge.'.format(_buffer))

        # NOTE: Should never happen, but this makes a compromise if a
        # trace crosses both the left and right spatial edge of the
        # detector.
        offset = np.clip(_buffer - np.amin(trace_cen, axis=0), 0, self.par['max_nudge']) \
                    + np.clip(self.nspat - 1 - _buffer - np.amax(trace_cen, axis=0),
                              None if self.par['max_nudge'] is None else -self.par['max_nudge'], 0)

        # TODO: This nudging can actually make the slit length
        # negative. This isn't handled in this method because it's only
        # meant to nudge the traces away from the detector edge,
        # independent of what edge it is or its counterpart. For now,
        # check_synced handles what to do if the slit length is
        # negative.

        # Offset and return the traces
        return trace_cen + offset[None,:]

    def sync(self, rebuild_pca=True, debug=False):
        """
        Match left and right edge traces to construct slit edge pairs.

        First, the method removes any fully masked traces (see
        :func:`clean_traces`). If this leaves the object empty, two
        traces are added at the edge of the detector; otherwise, the
        method ensures that the edge traces are sorted spatially (see
        :func:`spatial_sort`) and determines where traces need to be
        inserted to create a full set of left-right pairs (see
        :func:`_get_insert_locations`).

        The next steps are to determine the reference positions where
        the traces should be inserted and the shape the trace should
        take. The former is determined by the `sync_center` parameter
        in :attr:`par` and the latter is determined by `sync_predict`.

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
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are
        `det_buffer`, `left_right_pca`, and `sync_predict`.

        .. warning::
            Synchronizing the left and right edges requires that
            traces that are fully masked as bad must be removed (this
            does not include traces that are "masked" as having been
            deliberately inserted), and these traces are removed
            regardless of the user-specified `clip` in :attr:`par`.

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
        # Remove any fully masked traces
        self.clean_traces()

        # Make sure there are still traces left
        if self.is_empty:
            msgs.warn('No traces left!  Left and right edges placed at detector boundaries.')
            self.bound_detector()

        # Make sure that the traces are sorted spatially
        # TODO: This should be the convention of the class and should
        # *always* be true; instead check for this and raise an error
        # if it's not?
        self.spatial_sort()

        # If the traces are already synced, check them and log the
        # function as completed
        if self.is_synced:
            self.check_synced(rebuild_pca=rebuild_pca and self.pca is not None)
            self.log += [inspect.stack()[0][3]]
            return

        # Edges are currently not synced, so check the input
        if self.par['sync_predict'] not in ['pca', 'nearest']:
            msgs.error('Unknown trace mode: {0}'.format(self.par['sync_predict']))
        if self.par['sync_predict'] == 'pca' and self.pca is None:
            msgs.error('The PCA decomposition does not exist.  Either run self.build_pca or use '
                       'a different trace_mode.')

        # Find the edges to add, what side they're on, and where to
        # insert them into the existing trace array
        side, add_edge, add_indx = self._get_insert_locations()
        if not np.any(add_edge):
            # No edges to add
            return

        # Report
        msgs.info('-'*50)
        msgs.info('{0:^50}'.format('Synchronizing left and right traces'))
        msgs.info('-'*50)
        msgs.info('Found {0} left and {1} right trace(s) to add.'.format(
                    np.sum((side == -1) & add_edge), np.sum((side == 1) & add_edge)))

        # Allow the edges to be synced, even if a fit hasn't been done yet
        trace_cen = self.spat_cen if self.spat_fit is None else self.spat_fit

        # Instantiate the traces to add
        trace_add = np.zeros((self.nspec, np.sum(add_edge)), dtype=float)

        # If there was only one edge, just add the other one
        if side.size == 2:
            msgs.warn('Only one edge traced.  Ignoring center_mode and adding edge at the '
                      'opposite edge of the detector.')
            msgs.info('Detector edge buffer: {0}'.format(self.par['det_buffer']))
            # TODO: PCA would have failed because there needs to be at
            # least two traces. Get rid of this test once satisfied
            # that this exception is never raised...
            if self.par['sync_predict'] == 'pca':
                msgs.error('Coding error: this should not happen.')
            # Set the offset to add to the existing trace
            offset = self.par['det_buffer'] - np.amin(trace_cen[:,0]) if add_edge[0] \
                        else self.nspat - np.amax(trace_cen[:,0]) - self.par['det_buffer']
            # Construct the trace to add and insert it
            trace_add[:,0] = trace_cen[:,0] + offset
            self.insert_traces(side[add_edge], trace_add, loc=add_indx[add_edge], mode='sync')
            return

        # Get the reference locations for the new edges
        trace_ref = self._get_reference_locations(trace_cen, add_edge)

        # Predict the traces either using the PCA or using the nearest slit edge
        if self.par['sync_predict'] == 'pca':
            trace_add = self.predict_traces(trace_ref[add_edge], side[add_edge])
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
            # reference trace position to construct the new traces.
            trace_add = trace_cen[:,indx[nearest[add_edge]]] + trace_ref[add_edge] \
                            - trace_ref[nearest[add_edge]]

        # Insert the new traces and resort them spatially
        self.insert_traces(side[add_edge], trace_add, loc=add_indx[add_edge], mode='sync')

        # The sorted edges should now be arranged correctly. If not, it
        # should be because the traces were nudged away from the
        # detector edge and caused "negative" slit lengths...
        side = np.clip(self.traceid, -1, 1)
        indx = np.zeros(side.size, dtype=bool)
        indx[::2] = side[::2] != -1
        indx[1::2] = side[1::2] != 1
        if np.all(indx):
            msgs.error('Catastrophic error in left-right synchronization.  Edge order incorrect.')
        if np.any(indx):
            msgs.warn('Synchronized traces are not properly ordered, likely because they '
                      'have been placed close to the detector edges. Flagging '
                      '{0} traces that are not properly sorted for removal.'.format(np.sum(indx)))
            # Mask the traces as due to a synchronization error
            # NOTE: These are only masked here so that they can be
            # plotted if debug is True. Because the full traces are
            # masked, they're immediately removed by check_synced.
            self.spat_msk[:,indx] = self.bitmask.turn_on(self.spat_msk[:,indx], 'SYNCERROR')

        if debug:
            msgs.info('Show instance includes inserted traces but before checking the sync.')
            self.show(thin=10, include_img=True, idlabel=True, flag='any')

        # Check the full synchronized list and log completion of the
        # method
        self.check_synced(rebuild_pca=rebuild_pca)
        self.log += [inspect.stack()[0][3]]

    def add_user_traces(self, user_traces):
        """
        Add user-defined slit(s)

        Args:
            user_slits (list):

        """
        sides = []
        new_traces = np.zeros((self.nspec, len(user_traces)*2))
        # Add user input slits
        for kk, new_slit in enumerate(user_traces):
            # Parse
            y_spec, x_spat0, x_spat1 = new_slit
            msgs.info("Adding new slits at x0, x1 (left, right)".format(x_spat0, x_spat1))
            #
            # TODO -- Use existing traces (ideally the PCA) not just vertical lines!
            sides.append(-1)
            sides.append(1)
            # Trace cen
            new_traces[:,kk*2] = x_spat0
            new_traces[:,kk*2+1] = x_spat1
        # Insert
        self.insert_traces(np.array(sides), new_traces, mode='user')
        # Sync
        self.check_synced(rebuild_pca=False)

    def insert_traces(self, side, trace_cen, loc=None, mode='user', resort=True):
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

        .. warning::

            If the traces can be nudged away from the detector edge,
            the offset can, e.g., place an inserted left edge to the
            right of its associated right edge. This possibility is
            currently *not* handled by this function.

        Args:
            side (:obj:`int`, `numpy.ndarray`_):
                Side for each trace to be added: -1 for left, 1 for
                right. Shape is :math:`(N_{\rm new},)`.
            trace_cen (`numpy.ndarray`_):
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

                    - ``None``: Traces are simply inserted without
                      flagging.
                    - ``'user'``: Traces are the result of a user
                      request.
                    - ``'sync'``: Traces were generated by synchronizing
                      left and right traces.
                    - ``'mask'``: Traces were generated based on the
                      expected slit positions from mask design data.

            resort (:obj:`bool`, optional):
                Resort the traces in the spatial dimension; see
                :func:`spatial_sort`.

        """
        # Check input
        _side = np.atleast_1d(side)
        ntrace = _side.size
        _trace_cen = trace_cen.reshape(-1,1) if trace_cen.ndim == 1 else trace_cen
        if _trace_cen.shape[1] != ntrace:
            msgs.error('Number of sides does not match the number of traces to insert.')
        if loc is None:
            # Insertion locations not provided so append
            loc = np.full(ntrace, self.ntrace, dtype=int)
        if loc.size != ntrace:
            msgs.error('Number of sides does not match the number of insertion locations.')

        msgs.info('Inserting {0} new traces.'.format(ntrace))

        # Nudge the traces
        _trace_cen = self.nudge_traces(_trace_cen)

        # Set the mask
        mask = np.zeros(_trace_cen.shape, dtype=self.bitmask.minimum_dtype())
        # Flag the traces pixels that fall off the detector
        indx = (_trace_cen < 0) | (_trace_cen >= self.nspat)
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
            last_left = 0 if self.traceid is None or np.sum(self.is_left) == 0 \
                            else np.amin(self.traceid)
            _traceid[indx] = last_left - 1 - np.arange(np.sum(indx))
        indx = _side > 0
        if np.any(indx):
            last_right = 0 if self.traceid is None or np.sum(self.is_right) == 0 \
                            else np.amax(self.traceid)
            _traceid[indx] = last_right + 1 + np.arange(np.sum(indx))

        if self.is_empty:
            # No traces exist, so set the internals directly
            self.traceid = _traceid
            self.spat_img = np.round(_trace_cen).astype(int)
            self.spat_cen = _trace_cen
            self.spat_err = np.zeros(_trace_cen.shape, dtype=float)
            self.spat_msk = mask
            self.spat_fit = _trace_cen
            return

        # Add the new traces. The new traces are added to both the
        # fitted list and the center list!
        self.traceid = np.insert(self.traceid, loc, _traceid)
        self.spat_img = np.insert(self.spat_img, loc, np.round(_trace_cen).astype(int), axis=1)
        self.spat_cen = np.insert(self.spat_cen, loc, _trace_cen, axis=1)
        self.spat_err = np.insert(self.spat_err, loc,
                                  np.zeros(_trace_cen.shape, dtype=float), axis=1)
        self.spat_msk = np.insert(self.spat_msk, loc, mask, axis=1)
        self.spat_fit = np.insert(self.spat_fit, loc, _trace_cen, axis=1)

        if resort:
            self.spatial_sort()

    def bound_detector(self):
        """
        Insert traces at both detector boundaries.

        Accounting for the requested detector buffer, a left and
        right trace are placed at the detector edges. Traces are
        masked as user-inserted.

        Only used parameter from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) is
        `det_buffer`.
        """
        self.insert_traces(np.array([-1,1]),
                           np.array([np.full(self.nspec, self.par['det_buffer'], dtype='float'),
                                     np.full(self.nspec, self.nspat-1-self.par['det_buffer'],
                                             dtype='float')]).T)

    def fully_masked_traces(self, flag=None, exclude=None):
        """
        Find fully masked edge traces.

        Traces are identified as masked by one or more of the flags
        in `flag` and explicitly not flagged by any flag in `exclude`
        over the fully spectral range of the detector.
        
        Args:
            flag (:obj:`str`, :obj:`list`, optional):
                The bit mask flags to select. If None, any flags are
                used. See :func:`pypeit.bitmask.Bitmask.flagged`.
            exclude (:obj:`str`, :obj:`list`, optional):
                A set of flags to explicitly exclude from
                consideration as a masked trace. I.e., if any
                spectral pixel in the trace is flagged with one of
                these flags, it will not be considered a fully masked
                trace. This is typically used to exclude inserted
                traces from being considered as a bad trace.

        Returns:
            `numpy.ndarray`_: Boolean array selecting traces that are
            flagged at all spectral pixels.
        """
        if self.is_empty:
            return None
        bpm = np.all(self.bitmask.flagged(self.spat_msk, flag=flag), axis=0)
        if exclude is not None:
            bpm &= np.invert(np.any(self.bitmask.flagged(self.spat_msk, flag=exclude), axis=0))
        return bpm
    
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
        `left_right_pca`, `mask_reg_maxiter`, `mask_reg_maxsep`,
        `mask_reg_sigrej`, and `ignore_alignment`.

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

        # Check that there are traces to refine!
        if self.is_empty:
            msgs.error('No traces to refine.')

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
        reference_row = self.pca[0].reference_row if self.par['left_right_pca'] \
                            else self.pca.reference_row
        x_det = self.spat_fit[reference_row,:]

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
            raise NotImplementedError()
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
            missing_traces = self.predict_traces(register.match_coo[missing], side)
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
            reference_row = self.pca[0].reference_row if self.par['left_right_pca'] \
                            else self.pca.reference_row
            # Reset the match after removing/inserting traces
            x_det = self.spat_fit[reference_row,:]
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
        # Reference row
        reference_row = self.pca[0].reference_row if self.par['left_right_pca'] \
                            else self.pca.reference_row
        # Instantiate as an empty table
        self.design = EdgeTraceSet.empty_design_table(rows=nslits)
        # Save the fit parameters and the source file as table metadata
        self.design.meta['MASKFILE'] = design_file
        self.design.meta['MASKOFF'] = register.par[0]
        self.design.meta['MASKSCL'] = register.par[1]
        # Fill the columns
        self.design['TRACEID'] = np.arange(nslits, dtype=self.design['TRACEID'].dtype)
        self.design['TRACESROW'] = np.full(nslits, reference_row,
                                           dtype=self.design['TRACESROW'].dtype)
        self.design['TRACELPIX'] = self.spat_fit[reference_row,self.traceid<0].astype(
                                        dtype=self.design['TRACELPIX'].dtype)
        self.design['TRACERPIX'] = self.spat_fit[reference_row,self.traceid>0].astype(
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

    def slit_spatial_center(self, normalized=True, spec=None, resort=False, use_center=False):
        """
        Return the spatial coordinate of the center of each slit.

        The slit edges must be left-right synchronized.

        Args:
            normalized (:obj:`bool`, optional):
                Return coordinates normalized by the size of the
                detector.
            spec (:obj:`int`, optional):
                Spectral position (row) at which to return the
                spatial position. If ``None``, set at the central row
                (i.e., ``self.nspat//2``)
            resort (:obj:`bool`, optional):
                Ensure that the traces are sorted. If ``False``, they
                are assumed to be.
            use_center (:obj:`bool`, optional):
                Use the measured centroids to define the slit edges
                even if the slit edges have been otherwise modeled.

        Returns:
            `numpy.ndarray`_: Spatial coordinates of the slit centers
            in pixels or in fractions of the detector.
        """
        if not self.is_synced:
            msgs.error('EdgeTraceSet must be synced to compute slit centers.')

        # TODO: Use reference_row by default? Except that it's only
        # defined if the PCA is defined.
        _spec = self.nspat//2 if spec is None else spec

        # Resort if requested
        if resort:
            self.spatial_sort()

        # Synced, spatially sorted traces are always ordered in left,
        # right pairs
        trace_cen = self.spat_cen if self.spat_fit is None or use_center else self.spat_fit
        slit_c = np.mean(trace_cen[_spec,:].reshape(-1,2), axis=1)
        return slit_c/self.nspat if normalized else slit_c

    def get_slits(self):
        """
        Use the data to instatiate the relevant
        :class:`pypeit.slittrace.SlitTraceSet` object.

        The traced edges must have first been organized into slits;
        see :func:`sync`.

        The :class:`pypeit.slittrace.SlitTraceSet` object will use
        the same :attr:`master_key`, :attr:`master_dir`, and
        :attr:`reuse_masters` as this parent :class:`EdgeTraceSet`
        object.
        """
        if not self.is_synced:
            msgs.error('Edges must be synced to construct SlitTraceSet object.')

        # Select the good traces
        gpm = np.invert(self.fully_masked_traces(flag=self.bitmask.bad_flags,
                                                 exclude=self.bitmask.exclude_flags))

        # Parse the data for the SlitTraceSet
        left = self.spat_fit[:,gpm & self.is_left]
        right = self.spat_fit[:,gpm & self.is_right]
        binspec, binspat = parse.parse_binning(self.binning)
        slitspat = slittrace.slit_spat_pos(left, right, self.nspat)
        specmin, specmax = self.spectrograph.slit_minmax(slitspat, binspectral=binspec)

        # Instantiate and return
        return slittrace.SlitTraceSet(left=left, right=right, nspat=self.nspat,
                                      spectrograph=self.spectrograph.spectrograph, specmin=specmin,
                                      specmax=specmax, binspec=binspec, binspat=binspat,
                                      pad=self.par['pad'], master_key=self.master_key,
                                      master_dir=self.master_dir, reuse=self.reuse_masters)

