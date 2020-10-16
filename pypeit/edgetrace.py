# TODO: This docstring needs to be updated!!
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
    # TODO -- Update this doc
    traceImage = pypeit.images.buildcalibration.TraceImage(rdx.spectrograph, files=files, det=det,
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

    # TODO - update this
    traceImage = traceimage.TraceImage(spec, files=[trace_file], det=det,
                                       par=par['calibrations']['traceframe'])
    traceImage.build_image()

    edges = edgetrace.EdgeTraceSet(spec, par['calibrations']['slitedges'], master_key=master_key,
                                   master_dir=master_dir, img=traceImage, det=det, auto=True)
    edges.save()

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

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
from pypeit import io
from pypeit import utils
from pypeit import sampling
from pypeit import masterframe
from pypeit import slittrace
from pypeit.datamodel import DataContainer
from pypeit.bitmask import BitMask
from pypeit.display import display
from pypeit.par.pypeitpar import EdgeTracePar
from pypeit.core import parse, pydl, procimg, pca, trace, slitdesign_matching
from pypeit.images.buildimage import TraceImage
from pypeit.tracepca import TracePCA
from pypeit.spectrographs.spectrograph import Spectrograph
from pypeit.spectrographs.util import load_spectrograph


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
                       ('HITMIN', 'Trace crosses the minimum allowed spatial column'),
                       ('HITMAX', 'Trace crosses the maximum allowed spatial column'),
                  ('OFFDETECTOR', 'Trace lands off, or within `det_buffer` of, the detector edge'),
                   ('USERINSERT', 'Trace was inserted as requested by user'),
                   ('SYNCINSERT', 'Trace was inserted during left and right edge sync'),
                    ('SYNCERROR', 'Trace synchronization error, likely due to edge effects'),
                   ('MASKINSERT', 'Trace was inserted based on drilled slit-mask locations'),
                 ('ORPHANINSERT', 'Trace was inserted to match an orphaned edge'),
                    ('SHORTSLIT', 'Slit formed by left and right edge is too short to be valid'),
                      ('BOXSLIT', 'Slit formed by left and right edge is valid (large enough '
                                  'to be a valid slit), but too short to be a science slit'),
                 ('ABNORMALSLIT', 'Slit formed by left and right edge has abnormal length'),
                   ('USERRMSLIT', 'Slit removed by user'),
                      ('NOORDER', 'Unable to associate this trace with an echelle order (echelle '
                                  ' spectrographs only)'),
                ('ORDERMISMATCH', 'Slit traces are not well matched to any echelle order (echelle '
                                  ' spectrographs only)')])
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
    def order_flags(self):
        """
        List of flags related to the echelle order.
        """
        return ['NOORDER', 'ORDERMISMATCH']

    @property
    def exclude_flags(self):
        """
        List of flags to exclude when finding bad trace data.
        """
        return self.insert_flags + self.order_flags + ['BOXSLIT']


class EdgeTraceSet(DataContainer):
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

    To load an existing master file with the result of a trace, 
    use the :attr:`from_file` method::

        edges = EdgeTraceSet.from_file(file)

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
        traceimg (:class:`pypeit.images.buildcalibration.TraceImage`):
            Two-dimensional image used to trace slit edges. If a
            :class:`pypeit.images.buildcalibration.TraceImage` is provided, the
            raw files used to construct the image are saved.
            This includes detector info (e.g. binning, platescale)
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
        bpm (`numpy.ndarray`_, optional):
            Bad-pixel boolean mask for the trace image. Must have the
            same shape as `img`. If None, all pixels are assumed to
            be valid.
        auto (:obj:`bool`, optional):
            Find the edge traces using :func:`auto_trace`. If False,
            the trace data will only be the result of running
            :func:`initial_trace`.
        debug (:obj:`bool`, optional):
            Run in debug mode.
        show_stages (:obj:`bool`, optional):
            After ever stage of the auto trace prescription
            (`auto=True`), show the traces against the image using
            :func:`show`.

    Attributes:
        traceimg
            (:class:`pypeit.images.buildcalibration.TraceImage`):
            See argument list.
        spectrograph
            (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            See argument list.
        par (:class:`pypeit.par.pypeitpar.EdgeTracePar`):
            See argument list.
        files (:obj:`list`):
            The list of raw files used to construct the trace image
            (:attr:`img`). Only defined if argument `img` in
            :func:`initial_trace` or :func:`auto_trace` is a
            :class:`pypeit.images.buildcalibration.TraceImage` object.
        img (`numpy.ndarray`_):
            Convenience for now.
        bpm (`numpy.ndarray`_):
            See argument list.
        det (:obj:`int`):
            See argument list.
        sobelsig (`numpy.ndarray`_)):
            Sobel-filtered image used to detect left and right edges
            of slits.
        sobelsig_left (`numpy.ndarray`_):
            Lazy-loaded version of `sobelsig` that clips the
            features related to right edges. Only kept for
            convenience.
        sobelsig_right (`numpy.ndarray`_):
            Lazy-loaded version of `sobelsig` that clips the
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
        edge_img (`numpy.ndarray`_):
            An integer array with the spatial pixel nearest to each
            trace edge. This is identically::

                self.edge_img = np.round(self.edge_cen
                                         if self.edge_fit is None
                                         else self.edge_fit).astype(int)

        edge_cen (`numpy.ndarray`_):
            A floating-point array with the location of the slit edge
            for each spectral pixel *as measured* from the trace
            image. Shape is :math:`(N_{\rm spec},N_{\rm trace})`.
        edge_err (`numpy.ndarray`_):
            Error in slit edge locations; measurements without errors
            have their errors set to -1.
        edge_msk (`numpy.ndarray`_):
            An integer array with the mask bits assigned to each
            trace centroid; see :class:`EdgeTraceBitMask`.
        edge_fit (`numpy.ndarray`_):
            A model fit to the `edge_cen` data.
        fittype (:obj:`str`):
            An informational string identifier for the type of model
            used to fit the trace data.
        pca (:obj:`list`, :class:`pypeit.tracepca.TracePCA`):
            Result of a PCA decomposition of the edge traces, used to
            predict new traces. This can either be a single
            :class:`pypeit.tracepca.TracePCA` object or a list of two
            :class:`pypeit.tracepca.TracePCA` objects if the PCA
            decomposition is peformed for the left (`pca[0]`) and
            right (`pca[1]`) traces separately.
        pcatype (:obj:`str`)
            An informational string indicating which data were used
            in the PCA decomposition, 'center' for `edge_cen` or
            'fit' for `edge_fit`.
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
    """Root name for the master frame file."""

    master_file_format = 'fits.gz'
    """Master frame file format."""

    bitmask = EdgeTraceBitMask()
    """BitMask instance."""

    version = '1.0.0'
    """DataContainer datamodel version."""

    output_float_dtype = np.float32
    """Regardless of datamodel, output floating-point data have this fixed bit size."""

    datamodel = {'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 'dispname': dict(otype=str,
                                  descr='Spectrograph disperser name.  Primarily needed for '
                                        'reloading an existing MasterEdge file.'),
                 'traceimg': dict(otype=TraceImage,
                                   descr='Image used to construct the edge traces.'),
                 'nspec': dict(otype=int, descr='Image pixels in the spectral direction.'),
                 'nspat': dict(otype=int, descr='Image pixels in the spatial direction.'),
                 'tracebpm': dict(otype=np.ndarray, atype=(bool, np.bool_),
                                  descr='Bad-pixel mask for trace image'),
                 'sobelsig': dict(otype=np.ndarray, atype=(float, np.floating),
                                  descr='Sobel-filtered image used to detect edges'),
                 'traceid': dict(otype=np.ndarray, atype=(int, np.integer),
                                 descr='ID number for the edge traces.  Negative and positive '
                                       'IDs are for, respectively, left and right edges.'),
                 'orderid': dict(otype=np.ndarray, atype=(int, np.integer),
                                 descr='For echelle spectrographs, this is the order ID number '
                                       'for the edge traces.  Negative and positive IDs are for, '
                                       'respectively, left and right edges.'),
                 'edge_cen': dict(otype=np.ndarray, atype=(float, np.floating),
                                  descr='(Floating-point) Measured spatial coordinate of the '
                                        'edge traces for each spectral pixel.  Shape is '
                                        '(Nspec,Ntrace).'),
                 'edge_err': dict(otype=np.ndarray, atype=(float, np.floating),
                                  descr='Error in the measured spatial coordinate edge traces.'),
                 'edge_msk': dict(otype=np.ndarray, atype=bitmask.minimum_dtype(),
                                  descr='Bitmask for the edge trace positions.'),
                 'edge_fit': dict(otype=np.ndarray, atype=(float, np.floating),
                                  descr='The best-fit model result for the trace edge.'),
                 'fittype': dict(otype=str,
                                 descr='An informational string identifying the type of model '
                                       'used to fit the trace data.  Either ``pca`` for a PCA '
                                       'decomposition or the polynomial function type and order'),
                 'pca': dict(otype=TracePCA,
                             descr='The PCA decomposition of all edge traces.  Ignored if the '
                                   'left_right_pca parameter is True.'),
                 'left_pca': dict(otype=TracePCA,
                                  descr='The PCA decomposition of the left-edge traces.  Ignored '
                                        'if the left_right_pca parameter is False.'),
                 'right_pca': dict(otype=TracePCA,
                                   descr='The PCA decomposition of the right-edge traces.  '
                                         'Ignored if the left_right_pca parameter is False.'),
                 'pcatype': dict(otype=str,
                                 descr='String identifier for the measurements used to construct '
                                       'the PCA (center or fit)')}
    """DataContainer datamodel."""

    def __init__(self, traceimg, spectrograph, par, bpm=None, qa_path=None, auto=False,
                 debug=False, show_stages=False):

        # Instantiate as an empty DataContainer
        super(EdgeTraceSet, self).__init__()

        # TODO: TraceImage has a bad-pixel mask. Why isn't the
        # traceimg.bpm the same as the keyword argument bpm? Currently,
        # the traceimg.bpm array is not used. Should this be combined
        # with the provided bpm?

        # Check input types
        if not isinstance(traceimg, TraceImage):
            msgs.error('Input traceimg must be a TraceImage object.')
        if not isinstance(spectrograph, Spectrograph):
            msgs.error('Input spectrograph must be a Spectrograph object.')
        if not isinstance(par, EdgeTracePar):
            msgs.error('Input par must be an EdgeTracePar object.')

        # TODO:
        #   - Change spectrograph.spectrograph to spectrograph.name
        #   - Change self.PYP_SPEC to self.specname
        self.traceimg = traceimg                        # Input TraceImage
        self.nspec, self.nspat = self.traceimg.shape    # The shape of the trace image
        self.spectrograph = spectrograph                # Spectrograph used to take the data
        self.PYP_SPEC = spectrograph.spectrograph       # For the Header.  Will be in datamodel
        self.dispname = spectrograph.dispname           # Spectrograph disperser
        self.par = par                                  # Parameters used for slit edge tracing
        self.qa_path = qa_path                          # Directory for QA plots

        # NOTE: This means that, no matter what, every instance of
        # EdgeTraceSet should have a sobelsig attribute that is *not*
        # None.
        if auto:
            # Run the automatic tracing
            self.auto_trace(bpm=bpm, debug=debug, show_stages=show_stages)
        else:
            # Only get the initial trace
            self.initial_trace(bpm=bpm)

    def _init_internals(self):
        """Add any attributes that are *not* part of the datamodel."""
        self.spectrograph = None        # Spectrograph instance
        self.par = None                 # EdgeTracePar instance
        self.qa_path = None             # Path for the QA plots
        self.edge_img = None            # Array with the spatial pixel nearest to each trace edge.
        self.sobelsig_left = None      # Sobel filtered image used to trace left edges
        self.sobelsig_right = None     # Sobel filtered image used to trace right edges
        self.design = None              # Table that collates slit-mask design data matched to
                                        # the edge traces
        self.objects = None             # Table that collates object information, if available
                                        # in the slit-mask design, matched to the `design` table.
        self.log = None                 # Log of methods applied
        self.master_key = None          # Calibration key for master frame
        self.master_dir = None          # Directory for Master frames
        self.maskdef_id = None          # Slit ID number from slit-mask design to record in SlitTraceSet

    def _reinit_trace_data(self):
        """
        Convenience method to set all attributes related to trace data to `None`.
        """
        self.traceid = None
        self.orderid = None
        self.edge_img = None
        self.edge_cen = None
        self.edge_err = None
        self.edge_msk = None
        self.edge_fit = None
        self.fittype = None
        self.pca = None
        self.left_pca = None
        self.right_pca = None
        self.pcatype = None
        self.design = None
        self.objects = None
        self.maskdef_id = None

    @property
    def ntrace(self):
        """
        The number of edges (left and right) traced.
        """
        return 0 if self.traceid is None else self.traceid.size

    @property
    def nslits(self):
        if self.is_synced:
            return self.ntrace//2
        msgs.error('Number of slits undefined because edges are not left-right synchronized.')

    # TODO: Add self.design to the data model when we're ready to match
    # to the slit-mask design data.
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
                    # table.Column(name='SLITLFOC', dtype=float, length=length,
                    #              description='Left edge of the slit in mm at the focal plane'),
                    # table.Column(name='SLITRFOC', dtype=float, length=length,
                    #              description='Right edge of the slit in mm at the focal plane'),
                    table.Column(name='SLITLOPT', dtype=float, length=length,
                                description='Left edge of the slit in pixel from optical model'),
                    table.Column(name='SLITROPT', dtype=float, length=length,
                                description='Right edge of the slit in pixel from optical model'),
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
        if self.pcatype is None:
            msgs.error('Must first run the PCA analysis for the traces; run build_pca.')
        pca = (self.left_pca if side == 'left' else self.right_pca) \
                    if self.par['left_right_pca'] else self.pca

        # Get the traces that cross the reference spatial position at
        # the first and last pixels of the image
        first_last_trace = pca.predict(np.array([0,self.nspat-1]))
        # Use these two traces to define the spatial pixel coordinates
        # to sample
        start = np.ceil(np.amax(np.amin(first_last_trace, axis=1))).astype(int)
        buffer = self.nspat - np.floor(np.amin(np.amax(first_last_trace, axis=1))).astype(int) \
                    + start
        # Rectify the image
        # TODO: This has its limitations if the PCA is highly non-linear.
        ocol = np.arange(self.nspat+buffer)-start
        return sampling.rectify_image(flux, pca.predict(ocol), bpm=bpm, ocol=ocol,
                                      max_ocol=self.nspat-1, extract_width=extract_width,
                                      mask_threshold=mask_threshold)

    def auto_trace(self, bpm=None, debug=False, show_stages=False):
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
            - Use :func:`maskdesign_matching` to match the slit
              edge traces found with the ones predicted by the
              slit-mask design.
            - Use :func:`add_user_traces` and :func:`rm_user_traces`
              to add and remove traces as defined by the
              user-provided lists in the :attr:`par`.

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are
        `use_maskdesign`.

        Args:
            bpm (`numpy.ndarray`_, optional):
                Bad-pixel mask for the trace image. Must have the
                same shape as `img`. If None, all pixels are assumed
                to be valid.
            debug (:obj:`bool`, optional):
                Run in debug mode.
            show_stages (:obj:`bool`, optional):
                After ever stage of the auto trace, execute
                :func:`show` to show the results. These calls to show
                are always::

                    self.show()

        """
        # Perform the initial edge detection and trace identification
        self.initial_trace(bpm=bpm)
        if show_stages:
            self.show(title='Initial identification and tracing of slit edges')

        # Initial trace can result in no edges found
        if not self.is_empty:
            # Refine the locations of the trace using centroids of the
            # features in the Sobel-filtered image.
            self.centroid_refine()
            if show_stages:
                self.show(title='Refinement of edge locations')

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
                self.show(title='Polynomial fit to the edge locations')

        # Check if the PCA decomposition is possible; this should catch
        # long slits
        if self.par['auto_pca'] and self.can_pca():
            # Use a PCA decomposition to parameterize the trace
            # functional forms
            self.pca_refine(debug=debug)
            if show_stages:
                self.show(title='PCA refinement of the trace models')

            # Use the results of the PCA decomposition to rectify and
            # detect peaks/troughs in the spectrally collapsed
            # Sobel-filtered image, then use those peaks to further
            # refine the edge traces
            self.peak_refine(rebuild_pca=True, debug=debug)
            if show_stages:
                self.show(title='Result after re-identifying slit edges from a spectrally '
                                'collapsed image.')

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
            self.show(title='After synchronizing left-right traces into slits')

        # [DP] `maskdesign_matching` for now is only matching the traces found with the ones predicted by
        # the slit-mask design. If we include the piece of code in which `maskdesign_matching` recovers
        # missing traces, we should move the following 2 lines before `self.sync()`
        if self.par['use_maskdesign']:
            msgs.info('-' * 50)
            msgs.info('{0:^50}'.format('Matching traces to the slit-mask design'))
            msgs.info('-' * 50)
            self.maskdesign_matching(debug=debug)

        # First manually remove some traces, just in case a user
        # wishes to manually place a trace nearby a trace that
        # was automatically identified. One problem with adding
        # slits first is that we may have to sync the slits again.
        ad_rm = False
        if self.par['rm_slits'] is not None:
            rm_user_slits = trace.parse_user_slits(self.par['rm_slits'],
                                                   self.traceimg.detector.det, rm=True)
            if rm_user_slits is not None:
                ad_rm = True
                self.rm_user_traces(rm_user_slits)

        # Add user traces
        if self.par['add_slits'] is not None:
            add_user_slits = trace.parse_user_slits(self.par['add_slits'],
                                                    self.traceimg.detector.det)
            if add_user_slits is not None:
                ad_rm = True
                self.add_user_traces(add_user_slits)
        if show_stages and ad_rm:
            self.show(title='After user-dictated adding/removing slits')

        # TODO: Add a parameter and an if statement that will allow for
        # this.
        # `peak_refine` ends with the traces being described by a
        # polynomial. Instead finish by reconstructing the trace models
        # using the PCA decomposition
#        self.pca_refine(debug=debug)
#        if show_stages:
#            self.show()
            
        # TODO: Add mask_refine() when it's ready

        # Add this to the log
        self.log += [inspect.stack()[0][3]]

    def initial_trace(self, bpm=None):
        r"""
        Perform the initial trace of the image.

        This effectively reinstantiates the object and must be the
        first method called for tracing an image.  The algorithm:

            - Lightly boxcar smooths the trace image spectrally.
            - Replaces pixel columns and rows that are substantially
              masked, if a bad-pixel mask is provided (see ``bpm``).
            - Applies a Sobel filter to the trace image along the
              spatial axis (columns) to detect slit edges using steep
              positive gradients (left edges) and steep negative
              gradients (right edges). See
              :func:`pypeit.core.trace.detect_slit_edges`.
            - Follows the detected left and right edges along
              spectrally adjacent pixels to identify coherent traces.
              See :func:`pypeit.core.trace.identify_traces`.
            - Initializes the attributes that provide the trace
              position for each spectral pixel based on these
              results.

        Used parameters from :attr:`par`
        (:class:`~pypeit.par.pypeitpar.EdgeTracePar`) are
        ``filt_iter``, ``sobel_mode``, ``edge_thresh``, and
        ``follow_span``.

        Args:
            bpm (`numpy.ndarray`_, optional):
                Bad-pixel mask for the trace image. Must have the
                same shape as `img`. If None, all pixels are assumed
                to be valid.
        """
        msgs.info('-'*50)
        msgs.info('{0:^50}'.format('Initialize Edge Tracing'))
        msgs.info('-'*50)

        # TODO: Put in some checks that makes sure the relevant image
        # attributes are not None?

        self.tracebpm = np.zeros((self.nspec, self.nspat), dtype=bool) \
                            if bpm is None else bpm.astype(bool)
        if self.tracebpm.shape != self.traceimg.shape:
            msgs.error('Mask is not the same shape as the trace image.')

        # Lightly smooth the image before using it to trace edges
        # TODO: Make this filter size a parameter?
        _img = ndimage.uniform_filter(self.traceimg.image, size=(3, 1), mode='mirror')

        # Replace bad-pixel columns if they exist
        # TODO: Do this before passing the image to this function?
        # Instead of just replacing columns, replace all bad pixels...
        if np.any(self.tracebpm):

            # TODO: For tracing, we really only care about bad spec
            # lines, i.e. columns in the PypeIt frame. And this is how
            # BPM is oriented. I am now setting it to deal with both
            # options. Do we need to replace bad *rows* instead of bad
            # columns?

            # Replace bad columns and rows that cover more than half
            # the image
            for flip, axis in zip([False,True], [0,1]):
                bad_cols = np.sum(self.tracebpm, axis=axis) > (self.tracebpm.shape[axis]//2)
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
        self.sobelsig, edge_img \
                = trace.detect_slit_edges(_img, bpm=self.tracebpm,
                                          median_iterations=self.par['filt_iter'],
                                          sobel_mode=self.par['sobel_mode'],
                                          sigdetect=self.par['edge_thresh'])
        # Empty out the images prepared for left and right tracing
        # until they're needed.
        self.sobelsig_left = None
        self.sobelsig_right = None

        # Identify traces by following the detected edges in adjacent
        # spectral pixels.
        # TODO: Use the user parameter instead of the hard-coded
        # number?
#        minimum_spec_length = self.nspec * self.par['det_min_spec_length']
        minimum_spec_length = 50
        trace_id_img = trace.identify_traces(edge_img, follow_span=self.par['follow_span'],
                                             minimum_spec_length=minimum_spec_length)
        # TODO: Remove this
        inserted_edge = np.zeros(trace_id_img.shape, dtype=bool)

        # Check that edges were found
        if np.all(trace_id_img == 0):
            msgs.warn('No edges found!  Trace data will be empty.')
            self._reinit_trace_data()
            self.log = [inspect.stack()[0][3]]
            return

        # Set the ID image to a MaskedArray to ease some subsequent
        # calculations; pixels without a detected edge are masked.
        trace_id_img = np.ma.MaskedArray(trace_id_img, mask=trace_id_img == 0)

        # Find the set of trace IDs; left traces are negative, right
        # traces are positive
        self.traceid = np.unique(trace_id_img.compressed())

        # Initialize the mask bits array for the trace coordinates and
        # just begin by setting them all as having no edge.
        self.edge_msk = self.bitmask.turn_on(np.zeros((self.nspec, self.ntrace),
                                                dtype=self.bitmask.minimum_dtype()), 'NOEDGE')

        # Save the input trace edges and remove the mask for valid edge
        # coordinates
        self.edge_img = np.zeros((self.nspec, self.ntrace), dtype=int)
        for i in range(self.ntrace):
            row, col = np.where(np.invert(trace_id_img.mask)
                                    & (trace_id_img.data == self.traceid[i]))
            self.edge_img[row,i] = col
            self.edge_msk[row,i] = 0            # Turn-off the mask

            # Flag any insert traces
            row, col = np.where(np.invert(trace_id_img.mask) & inserted_edge
                                    & (trace_id_img.data == self.traceid[i]))
            if len(row) > 0:
                self.edge_msk[row,i] = self.bitmask.turn_on(self.edge_msk[row,i], 'ORPHANINSERT')

        # Instantiate objects to store the floating-point trace
        # centroids and errors.
        self.edge_cen = self.edge_img.astype(float)   # This makes a copy
        self.edge_err = np.zeros((self.nspec, self.ntrace), dtype=float)

        # No fitting has been done yet
        self.fittype = None
        self.edge_fit = None

        # No PCA has been constructed yet
        self.pca = None
        self.left_pca = None
        self.right_pca = None
        self.pcatype = None

        # No design or object data yet
        # TODO: Move these to the datamodel when we have the slit-mask
        # design matching working.
        self.design = None
        self.objects = None

        # Restart the log
        self.log = [inspect.stack()[0][3]]

    def _base_header(self, hdr=None):
        """
        Construct the baseline header for all HDU extensions.

        This appends the :class:`EdgeTracePar` and
        :class:`EdgeTraceBitMask` data to the headers of all HDU
        extensions. This is overkill, but avoids overriding
        :func:`pypeit.datamodel.DataContainer.to_hdu` by just
        including the data in all HDU extensions.

        Args:
            hdr (`astropy.io.fits.Header`_, optional):
                Baseline header to add to all returned HDUs. If None,
                set by :func:`pypeit.io.initialize_header()`.

        Returns:
            `astropy.io.fits.Header`_: Header object to include in
            all HDU extensions.
        """
        _hdr = super(EdgeTraceSet, self)._base_header(hdr=hdr)
        _hdr['QAPATH'] = 'None' if self.qa_path is None else self.qa_path
        self.par.to_header(_hdr)
        self.bitmask.to_header(_hdr)
        return _hdr

    def _bundle(self):
        """
        Bundle the object for writing using
        :func:`~pypeit.datamodel.DataContainer.to_hdu`.
        """

        # TODO: This will need to be edited if/when astropy.table.Table
        # objects with the slit matching results get added to the
        # datamodel.

        # All of the datamodel elements that can be written as a header
        # are written to the SOBELSIG extension.
        d = []
        for key in self.keys():
            if self[key] is None:
                continue
            if isinstance(self[key], DataContainer):
                d += [{key: self[key]}]
                continue
            if self.datamodel[key]['otype'] == np.ndarray:
                # TODO: We might want a general solution in
                # DataContainer that knows to convert (and revert)
                # boolean arrays to np.uint8 types. Does anyone know a
                # way around this? It's annoying that ImageHDUs cannot
                # have boolean type arrays.
                if issubclass(self[key].dtype.type, (bool, np.bool_)):
                    d += [{key: self[key].astype(np.uint8)}]
                elif issubclass(self[key].dtype.type, (float, np.floating)):
                    d += [{key: self[key].astype(self.output_float_dtype)}]
                else:
                    d += [{key: self[key]}]

        # Find the element with sobelsig and add all the non-array,
        # non-DataContainer elements to that extension. This is done
        # after the first loop to just to make sure that the order of
        # extensions matches the order in the datamodel (although
        # that's not really necessary).
        for _d in d:
            if list(_d.keys()) != ['sobelsig']:
                continue
            for key in self.keys():
                if self[key] is None or self.datamodel[key]['otype'] == np.ndarray \
                        or isinstance(self[key], DataContainer):
                    continue
                _d[key] = self[key]

        return d

    def to_hdu(self, **kwargs):
        """
        Construct the HDUList to write.

        This overrides :func:`pypeit.datamodel.DataContainer.to_hdu`
        and has the same calling sequence (i.e., ``kwargs`` are
        passed directly to the base class method). This is needed to
        deal with the multiple :class:`~pypeit.tracepca.TracePCA`
        objects included in the output.

        Returns:
            :obj:`list`, `astropy.io.fits.HDUList`_: A list of HDUs,
            where the type depends on the value of ``add_primary``.
        """
        if self.pcatype is None or self.par['left_right_pca'] is False:
            # Do not need to change the default behavior if the PCA
            # doesn't exist or there is only one PCA for both left and
            # right edges.
            return super(EdgeTraceSet, self).to_hdu(**kwargs)

        # TODO: We need a better solution for multiple levels of nested
        # DataContainers. Here the commpication is that we're writing
        # many DataContainers of a single derived class (TracePCA) to
        # the same file.

        # Temporarily erase the left and right pca so
        # they're not written
        _left_pca, _right_pca = self.left_pca, self.right_pca
        self.left_pca, self.right_pca = None, None

        # Run the default (with add_primary = False)
        hdu = super(EdgeTraceSet, self).to_hdu(**kwargs)

        # Reset them
        self.left_pca, self.right_pca = _left_pca, _right_pca

        # Add in the left and right PCAs with the correct extension
        # names. Can only append on HDU at a time ...
        for h in self.left_pca.to_hdu(add_primary=False, hdu_prefix='LEFT_'):
            hdu.append(h)
        for h in self.right_pca.to_hdu(add_primary=False, hdu_prefix='RIGHT_'):
            hdu.append(h)
        # Finally return the full list
        return hdu

    @classmethod
    def from_hdu(cls, hdu, hdu_prefix=None, chk_version=True):
        """
        Instantiate the object from an HDU extension.

        This overrides the base-class method. Overriding this method
        is preferrable to overriding the ``_parse`` method because it
        makes it easier to deal with the muliple
        :class:`~pypeit.datamodel.DataContainer` objects contained by
        :class:`EdgeTraceSet`.

        Args:
            hdu (`astropy.io.fits.HDUList`_, `astropy.io.fits.ImageHDU`_, `astropy.io.fits.BinTableHDU`_):
                The HDU(s) with the data to use for instantiation.
            hdu_prefix (:obj:`str`, optional):
                Maintained for consistency with the base class but is
                not used by this method.
            chk_version (:obj:`bool`, optional):
                If True, raise an error if the datamodel version or
                type check failed. If False, throw a warning only.
        """
        # Run the default parser to get most of the data. This won't
        # parse traceimg because it's not a single-extension
        # DataContainer. It *will* parse pca, left_pca, and right_pca,
        # if they exist, but not their model components.
        d, version_passed, type_passed, parsed_hdus = super(EdgeTraceSet, cls)._parse(hdu)
        if not type_passed:
            msgs.error('The HDU(s) cannot be parsed by a {0} object!'.format(cls.__name__))
        if not version_passed:
            _f = msgs.error if chk_version else msgs.warn
            _f('Current version of {0} object in code (v{1})'.format(cls.__name__, cls.version)
               + ' does not match version used to write your HDU(s)!')

        # Instantiate the TraceImage from the header
        d['traceimg'] = TraceImage.from_hdu(hdu, chk_version=chk_version)

        # Check if there should be any PCAs
        parsed_pcas = np.any(['PCA' in h for h in parsed_hdus]) 
        if d['pcatype'] is not None and not parsed_pcas:
            msgs.error('CODING ERROR: Expect to parse PCA headers if pcatype is present.')

        # Instantiate the TracePCAs using the appropriate hdus.
        if d['pcatype'] is not None:
            if 'PCA' in parsed_hdus:
                d['pca'] = TracePCA.from_hdu(hdu, hdu_prefix=None)
            if 'LEFT_PCA' in parsed_hdus:
                d['left_pca'] = TracePCA.from_hdu(hdu, hdu_prefix='LEFT_')
            if 'RIGHT_PCA' in parsed_hdus:
                d['right_pca'] = TracePCA.from_hdu(hdu, hdu_prefix='RIGHT_')

        # Convert data types
        for key in d.keys():
            if d[key] is not None and cls.datamodel[key]['otype'] == np.ndarray:
                d[key] = d[key].astype(cls.datamodel[key]['atype'][0] 
                                        if hasattr(cls.datamodel[key]['atype'], '__len__') 
                                        else cls.datamodel[key]['atype'])

        # Instantiate
        self = super(EdgeTraceSet, cls).from_dict(d=d)

        # Set the integer pixel values
        self.edge_img = None if self.traceid is None \
                            else (np.round(self.edge_cen if self.edge_fit is None
                                    else self.edge_fit).astype(int))

        # Reinstantiate elements that are not part of the datamodel
        # NOTE: The SOBELSIG extension is used here as the reference
        # header because some headers are overwritten by the TraceImage
        # or TracePCA DataContainers
        self.spectrograph = load_spectrograph(hdu['SOBELSIG'].header['PYP_SPEC'])
        self.spectrograph.dispname = self.dispname
        self.par = EdgeTracePar.from_header(hdu['SOBELSIG'].header)
        self.qa_path = hdu['SOBELSIG'].header['QAPATH']

        # Check the bitmasks
        hdr_bitmask = BitMask.from_header(hdu['SOBELSIG'].header)
        if chk_version and hdr_bitmask.bits != self.bitmask.bits:
            msgs.error('The bitmask in this fits file appear to be out of date!  Recreate this '
                       'master frame either by rerunning run_pypeit or pypeit_trace_edges.')

        return self

    @property
    def flagged_bits(self):
        """
        List the flags that have been tripped.

        Returns:
            list: List of the unique flags tripped by any trace
            datum.
        """
        return [] if self.edge_msk is None \
                    else np.unique(np.concatenate([self.bitmask.flagged_bits(b) 
                                                    for b in np.unique(self.edge_msk)])).tolist()

    # TODO: Break the ginga commands out into a separate method?
    def show(self, include_error=False, thin=10, in_ginga=False, include_img=True,
             include_sobel=False, img_buffer=100, flag=None, idlabel=True, slits=None,
             original=False, title=None):
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
                Label each trace by their ID numbers. If the data is
                from an echelle spectrograph and the echelle order
                has been matched to the synchronized slit-trace
                pairs, the ID labels are based on the order matching.
            slits (:class:`~pypeit.slittrace.SlitTraceSet`, optional):
                Slits to plot instead of those kept internally. Note
                that, if this is provided, the "modeled" and
                "measured" slit locations are identical.
            original (:obj:`bool`, optional):
                When ``slits`` are provided and tweaked slits are
                available, this selects which traces to show. If
                True, show the original slit traces; if False, show
                the tweaked ones.
            title (:obj:`str`, optional):
                Title for the matplotlib plot (ignored if shown in
                ginga).  If None, plot is not given a title.
        """
        if include_img and include_sobel:
            msgs.error('Cannot show both the trace image and the filtered version.')

        # TODO: Clean and consolidate the objects needed for either the
        # ginga or matplotlib methods so that this isn't as onerous.

        # TODO: Currently barfs if object is empty!

        # Build the slit edge data to plot.
        if slits is None:
            # Use the internals. Any masked data is excluded; masked
            # data to be plotted are held in a separate array. This
            # means that errors and fits are currently never plotted
            # for masked data.
            _flag = None if flag in ['any', None] else np.atleast_1d(flag)
            cen = np.ma.MaskedArray(self.edge_cen, mask=self.bitmask.flagged(self.edge_msk))
            fit = self.edge_fit
            err = np.ma.MaskedArray(self.edge_err, mask=np.ma.getmaskarray(cen).copy())
            msk = None if flag is None \
                    else np.ma.MaskedArray(self.edge_cen, mask=np.invert(
                                           self.bitmask.flagged(self.edge_msk, flag=_flag)))
            is_left = self.is_left
            is_right = self.is_right
            _include_error = include_error
            # This selects all unmasked traces, traces that were
            # inserted, and traces that have been identified as coming
            # from box slits.
            gpm = np.logical_not(self.fully_masked_traces(flag=self.bitmask.bad_flags,
                                                          exclude=self.bitmask.exclude_flags))
            traceid = self.traceid if self.orderid is None else self.orderid
            nslits = np.amax(traceid)   # Only used if synced is True
            synced = self.is_synced
            slit_ids = None
            if synced:
                _trc = cen if fit is None else fit
                half = _trc.shape[0] // 2
                slit_ids = ((_trc[half, gpm & is_left]
                             + _trc[half, gpm & is_right]) / 2.).astype(int)
        else:
            # Use the provided SlitTraceSet
            _include_error = False
            if include_error:
                msgs.warn('SlitTraceSet object has no errors.')
            left, right, _ = slits.select_edges() # original=original)
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
            # JFH Trying to make this work for the case where slits are
            # known, but this functionality appears to be never used? I
            # think slits needs its own show method.
            slit_ids = slits.slitord_id

        # TODO: The above should set everything (self shouldn't be used
        # below). We need a single show method that both EdgeTraceSet
        # and SlitTraceSet can use.

        if in_ginga:
            # Set up the appropriate keyword arguments for the IDs
            id_kwargs = {'slit_ids': slit_ids} if synced \
                            else {'left_ids': traceid[gpm & is_left],
                                  'right_ids': traceid[gpm & is_right]}
            _trc = cen if fit is None else fit

            # Connect to or instantiate ginga window
            display.connect_to_ginga(raise_err=True, allow_new=True)
            # Clear the viewer and show the trace image
            trace_viewer, trace_ch = display.show_image(self.traceimg.image, chname='Trace Image',
                                                        clear=True)
            if not self.is_empty:
                display.show_slits(trace_viewer, trace_ch, _trc[:,gpm & is_left],
                                   _trc[:,gpm & is_right], pstep=thin, synced=synced, **id_kwargs)

            # Show the Sobel sigma image (do *not* clear)
            sobel_viewer, sobel_ch = display.show_image(self.sobelsig, chname='Sobel Filtered')
            if not self.is_empty:
                display.show_slits(sobel_viewer, sobel_ch, _trc[:,gpm & is_left],
                                   _trc[:,gpm & is_right], pstep=thin, synced=synced, **id_kwargs)
            return

        # Show the traced image
        if include_img:
            img_zlim = utils.growth_lim(self.traceimg.image, 0.95, fac=1.05)
            plt.imshow(self.traceimg.image, origin='lower', interpolation='nearest',
                       aspect='auto', vmin=img_zlim[0], vmax=img_zlim[1])
        elif include_sobel:
            sob_zlim = utils.growth_lim(self.sobelsig, 0.95, fac=1.05)
            plt.imshow(self.sobelsig, origin='lower', interpolation='nearest', aspect='auto',
                       vmin=sob_zlim[0], vmax=sob_zlim[1])

        if title is not None:
            plt.title(title)
        plt.ylabel('Spectral pixel index')
        plt.xlabel('Spatial pixel index')
        plt.xlim(-img_buffer, self.nspat+img_buffer)
        plt.ylim(-img_buffer, self.nspec+img_buffer)

        if self.is_empty:
            msgs.info('No traces defined.')
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
                indx = np.logical_not(np.ma.getmaskarray(cen)[:,i])
                if not np.any(indx):
                    continue
                _spec = spec[:,i][indx]
                _cen = cen[:,i][indx]
                plt.text(_cen[_spec.size//2], _spec[_spec.size//2], str(traceid[i]),
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
                plt.text(fit[self.nspec//2,i], spec[self.nspec//2,i], str(traceid[i]), color='k',
                         fontsize=16, alpha=0.7, zorder=10)

        # Limits and labels
        if np.any(is_left & gpm):
            left_line[0].set_label('left edge fit')
        if np.any(is_right & gpm):
            right_line[0].set_label('right edge fit')
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
        img_zlim = utils.growth_lim(self.traceimg.image, 0.95, fac=1.05)
        sob_zlim = utils.growth_lim(self.sobelsig, 0.95, fac=1.05)

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
            indx = np.logical_not(self.bitmask.flagged(self.edge_msk[:,i],
                                                       flag=self.bitmask.bad_flags))
            xlim = utils.growth_lim(self.edge_cen[indx,i], 1.0, fac=2.0)
            if min_spat is not None and np.diff(xlim) < min_spat:
                xlim = np.sum(xlim)/2 + np.array([-1,1])*min_spat/2

            # Plot the trace image and the fit (if it exists)
            ax = fig.add_axes([ax_x0, ax_y, delt[0]/2, delt[1]])
            ax.minorticks_on()
            ax.tick_params(which='major', length=10, direction='in', top=True, right=True)
            ax.tick_params(which='minor', length=5, direction='in', top=True, right=True)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.imshow(self.traceimg.image, origin='lower', interpolation='nearest',
                      vmin=img_zlim[0], vmax=img_zlim[1], aspect='auto')
            if self.edge_fit is not None:
                ax.plot(self.edge_fit[:,i], spec, color='C3' if self.traceid[i] < 0 else 'C1')
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
            ax.imshow(self.sobelsig, origin='lower', interpolation='nearest', vmin=sob_zlim[0],
                      vmax=sob_zlim[1], aspect='auto')
            if self.edge_fit is not None:
                ax.plot(self.edge_fit[:,i], spec, color='C3' if self.traceid[i] < 0 else 'C1')
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
        :func:`~pypeit.core.trace.prepare_sobel_for_trace`; the
        boxcar smoothing is always 5 pixels. Barring a instantiation
        of the object, this calculation is only done once per side by
        "lazy loading" :attr:`sobelsig_left` or
        :attr:`sobelsig_right`.

        Args:
            side (:obj:`str`):
                The side to return; must be ``'left'`` or ``'right'``.
    
        Returns:
            `numpy.ndarray`_: The manipulated Sobel image relevant to
            tracing the specified edge side.

        Raises:
            PypeItError:
                Raised if ``side`` is not ``'left'`` or ``'right'``.
        """
        # TODO: Add boxcar to EdgeTracePar?
        boxcar = 5
        if side == 'left':
            if self.sobelsig_left is None:
                self.sobelsig_left = trace.prepare_sobel_for_trace(self.sobelsig,
                                                                   bpm=self.tracebpm,
                                                                   boxcar=boxcar, side='left')
            return self.sobelsig_left
        if side == 'right':
            if self.sobelsig_right is None:
                self.sobelsig_right = trace.prepare_sobel_for_trace(self.sobelsig,
                                                                    bpm=self.tracebpm,
                                                                    boxcar=boxcar, side='right')
            return self.sobelsig_right
        msgs.error('Side must be left or right.')

    def centroid_refine(self, follow=True, start_indx=None, continuous=False, use_fit=False):
        """
        Refine the edge positions using a moment analysis and assess
        the results.

        The starting point for all trace positions is
        :attr:`edge_cen`, unless ``use_fit`` is True.

        The method runs in two primary modes, depending on the
        value of ``follow``:

        When ``follow=False``, the function simply executes
        :func:`~pypeit.core.trace.masked_centroid` to recenter
        *all* the locations for the left and right edges at each
        spectral position (row) independently. The maximum shift
        between any input and output centroid is ``max_shift_abs``
        from :attr:`par`.

        When ``follow=True``, the method uses
        :func:`~pypeit.core.trace.follow_centroid` to recenter each
        edge starting from the ``start_indx`` spectral index position
        (row) and then follows the centroid to higher and lower
        spectral rows. In this case, the center of the aperture used
        for each spectral row is the centroid of the trace measured
        for the preceding spectral row. The maximum shift between the
        input and output center for the first position analyzed is
        set by ``max_shift_abs`` and the maximum shift for each
        subsequent spectral row is set by ``max_shift_adj``; both
        parameters are provided in :attr:`par`. In this approach,
        it's typically best to let the method determine the starting
        spectral row instead of providing it directly. If left to its
        own devices, it will iterate through all the traces
        collecting those that all cross a specific spectral row into
        groups that it can follow simultaneously. If a starting
        specral row is provided directly, all traces must cross that
        row.

        Regardless of the value of ``follow``,
        :func:`~pypeit.core.trace.masked_centroid` is run with
        uniform weighting and an aperture width set to be twice
        ``fwhm_uniform`` from :attr:`par`.

        During the refinement, :func:`check_traces` after each new
        trace is refined, indentifying repeat and errorneous traces.
        (The minimum length of the trace passed to
        :func:`check_traces` is set by ``det_min_spec_length`` in
        :attr:`par`.) If ``clip`` in :attr:`par` is True, all bad
        traces are removed (see :func:`remove_traces`).

        Other used parameters from :attr:`par`
        (:class:`~pypeit.par.pypeitpar.EdgeTracePar`) are
        ``max_shift_abs``, ``max_shift_adj``, and ``max_spat_error``.

        .. warning::

            - This function modifies the internal trace arrays **in
              place**.
            - Because this changes :attr:`edge_cen` and
              :attr:`edge_err`, any model fitting of these data are
              erased by this function! I.e., :attr:`edge_fit` and
              :attr:`fittype` are set to None.
            - This *always* removes the PCA if it exists.

        Args:
            follow (:obj:`bool`, optional):
                Perform the centroiding first at a single row (see
                ``start_indx``) and then move to higher and lower
                spectral rows in series. See
                :func:`~pypeit.core.trace.follow_centroid`.
            start_indx (:obj:`int`, optional):
                The index of the starting spectral row when following
                the trace between adjacent rows; see ``follow``. If
                None, the starting spectral row is set by finding the
                spectral row that crosses the most unmasked trace
                positions; see
                :func:`~pypeit.core.trace.most_common_trace_row`.
                Value is ignored if ``follow=False``.
            continuous (:obj:`bool`, optional):
                Keep only the continuous part of the traces from the
                starting spectral row. See
                :func:`~pypeit.core.trace.follow_centroid`.
            use_fit (:obj:`bool`, optional):
                Use the fitted traces as the starting point for the
                refinement. Otherwise, uses :`edge_cen`. If True and
                :attr:`edge_fit` is None, the method will raise an
                exception.
        """
        # Check that there are traces to refine!
        if self.is_empty:
            # TODO: Continue to have this fault?
            msgs.error('No traces are defined.')

        # Check input
        if use_fit and self.edge_fit is None:
            msgs.error('No fit data available.')

        # Parse parameters and report
        width = 2 * self.par['fwhm_uniform']
        maxshift_start = self.par['max_shift_abs']
        maxshift_follow = self.par['max_shift_adj']
        maxerror = self.par['max_spat_error']
        minimum_spec_length = self.par['det_min_spec_length']*self.nspec

        # Report
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
        ivar = np.ones_like(self.sobelsig, dtype=float)
        _bpm = np.zeros_like(self.sobelsig, dtype=bool) \
                    if self.tracebpm is None else self.tracebpm
        fwgt = np.ones_like(self.sobelsig, dtype=float)

        # Book-keeping objects to keep track of which traces to remove
        rmtrace = np.zeros(self.ntrace, dtype=bool)

        # To hold the refined traces and mask
        spat = self.edge_fit if use_fit else self.edge_img
        cen = np.zeros_like(self.edge_cen)
        err = np.zeros_like(self.edge_err)
        msk = np.zeros_like(self.edge_msk, dtype=self.bitmask.minimum_dtype())

        # Ignore inserted traces
        inserted = self.fully_masked_traces(flag=self.bitmask.insert_flags)
        if np.any(inserted):
            cen[:,inserted] = self.edge_cen[:,inserted]
            err[:,inserted] = self.edge_err[:,inserted]
            msk[:,inserted] = self.edge_msk[:,inserted]

        # Refine left then right
        for side in ['left', 'right']:
            
            # Get the image relevant to tracing this side
            _sobelsig = self._side_dependent_sobel(side)

            # Find the traces to refine, must be on the correct side
            # and must not have been inserted
            indx = (self.is_left if side == 'left' else self.is_right) & np.logical_not(inserted)
            if not np.any(indx):
                continue

            msgs.info('Found {0} {1} edge trace(s) to refine'.format(np.sum(indx), side))

            if follow:
                # Find the bad trace positions
                bpm = self.bitmask.flagged(self.edge_msk, flag=self.bitmask.bad_flags)
                untraced = indx.copy()
                while np.any(untraced):
                    # Get the starting row
                    _start_indx = trace.most_common_trace_row(bpm[:,untraced], valid_frac=1.) \
                                        if start_indx is None else start_indx
                    # Select the edges to follow
                    to_trace = untraced & np.logical_not(bpm[_start_indx,:])
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
                            = trace.follow_centroid(_sobelsig, _start_indx,
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
                        = trace.masked_centroid(_sobelsig, spat[:,indx], width, ivar=ivar,
                                                bpm=_bpm, fwgt=fwgt, maxshift=maxshift_start,
                                                maxerror=maxerror, bitmask=self.bitmask,
                                                fill='bound')

            # Include previous mask in result
            # TODO: This means if the edge was not detected by
            # detect_slit_edges, it won't get used here. Need to check
            # how this meshes with the previous behavior.
            msk[:,indx] |= self.edge_msk[:,indx]

            # Check the traces
            min_spatial = None if side == 'left' else 0
            max_spatial = _sobelsig.shape[1]-1 if side == 'left' else None
            good, bad = self.check_traces(cen=cen, msk=msk, subset=indx, min_spatial=min_spatial,
                                          max_spatial=max_spatial,
                                          minimum_spec_length=minimum_spec_length)

            # Save the results: only keep the good centroids, but merge all the masks
            self.edge_cen[:,good] = cen[:,good]
            self.edge_err[:,good] = err[:,good]
            self.edge_msk[:,indx] |= msk[:,indx]
            rmtrace[bad] = True

        # Update the image coordinates
        self.edge_img = np.round(self.edge_cen).astype(int)

        # Erase any previous fitting, PCA, and slit-mask-match results
        # TODO: Should probably get sectioned off as a method because a
        # few methods do this.  Add an option to _reinit_trace_data?
        self.fittype = None
        self.edge_fit = None
        self.pcatype = None
        self.pca = None
        self.left_pca = None
        self.right_pca = None
        self.design = None
        self.objects = None

        # Remove bad traces and re-order the trace IDs
        if self.par['clip']:
            self.remove_traces(rmtrace)

        # Add to the log
        self.log += [inspect.stack()[0][3]]

    def trace_pixels_off_detector(self, cen=None):
        r"""
        Determine if the traces land outside of the detector and the
        defined edge buffer (``det_buffer`` in :attr:`par`).

        Args:
            cen (`numpy.ndarray`_, optional):
                Array with the trace locations. Shape is
                :math:`(N_{\rm spec}, N_{\rm trace})`. If None, uses
                :attr:`edge_cen`.
    
        Returns:
            `numpy.ndarray`_: Boolean array with shape :math:`(N_{\rm
            spec}, N_{\rm trace})` selecting pixels that are (True)
            or are not (False) "off" the detector.

        Raises:
            PypeItError:
                Raised if ``cen`` is not provided and
                :attr:`edge_cen` is None.
        """
        buff = 0 if self.par['det_buffer'] is None else self.par['det_buffer']
        if buff < 0:
            msgs.warn('Detector buffer must be >=0 (input was {0}).  Setting buffer to 0.'.format(
                      self.par['det_buffer']))
            buff = 0
        if cen is None:
            cen = self.edge_cen
        if cen is None:
            msgs.error('No trace locations!')
        return (cen < buff) | (cen > self.nspat - buff)

    def check_traces(self, cen=None, msk=None, subset=None, min_spatial=None, max_spatial=None,
                     minimum_spec_length=None):
        r"""
        Validate new or existing trace data.

        This function *only* changes the mask status of the traces and
        returns flags for traces that are good or bad.  Traces are
        flagged as bad if they meet any of the following conditions:

            - If and only if ``subset`` is provided, the selected subset
              of traces are checked against all other traces to find
              repeats.  Repeats are identified as those tracing the same
              side of a slit and having a minimum spatial separation of
              less than the provided ``match_tol`` parameter in
              :attr:`par`.  Beware that this should be done by providing
              ``cen`` directly.
            - Traces that do not extend for at least some fraction
              of the detector (see ``minimum_spec_length``).
            - Traces with spatial positions at the central spectral
              channel that are above or below the provided threshold
              (see ``min_spatial`` and ``max_spatial``).
            - Traces that fully land off the edge of the detector; see
              :func:`trace_pixels_off_detector`.

        The only used parameter from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) is ``match_tol``.

        .. warning::

            :attr:`msk` is edited in-place

        Args:
            cen (`numpy.ndarray`_, optional):
                The adjusted center of the refined traces. Shape is
                :math:`(N_{\rm spec}, N_{\rm refine},)`. If None, use
                :attr:`edge_cen`.
            msk (`numpy.ndarray`_, optional):
                The mask bits for the adjusted center of the refined
                traces. Shape is :math:`(N_{\rm spec}, N_{\rm
                refine},)`. If None, use :attr:`edge_msk`. This is
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
            :obj:`tuple`: Returns two boolean arrays selecting the good
            and bad traces.  Both are returned because neither will be
            true for traces not selected by ``subset`` if it is
            provided.  The shape of both arrays is :math:`(N_{\rm
            trace},)`.

        """
        if self.is_empty:
            msgs.warn('No traces to check.')

        # Indices of traces to check
        indx = np.ones(self.ntrace, dtype=bool) if subset is None else subset
        # Data to compare
        _cen = self.edge_cen if cen is None else cen
        _msk = self.edge_msk if msk is None else msk
        # Ignore inserted traces when flagging
        _bpm = self.bitmask.flagged(_msk, flag=self.bitmask.bad_flags)

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
                edge_bpm = self.bitmask.flagged(self.edge_msk, flag=self.bitmask.bad_flags)
                edge_img = np.ma.MaskedArray(self.edge_img, mask=edge_bpm)
                mindiff = np.ma.amin(np.absolute(_col[:,indx,None]-edge_img[:,None,compare]),
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
        pix_off_detector = np.zeros(_cen.shape, dtype=bool)
        pix_off_detector[:,indx] = self.trace_pixels_off_detector(cen=_cen[:,indx])
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
            - the *fitted* trace locations of two different traces
              are separated by less than ``match_tol`` (from
              :attr:`par`) for more than ``merge_frac`` of the
              spectral range of the detector.

        Since the merging is based on the *fitted* trace location,
        the traces must have been fit beforehand; see
        :func:`fit_refine`.

        When there are traces found that should be merged, the
        unmasked centroid measurements from the shorter trace(s) are
        added to the longest trace (meaning that the shorter trace is
        effectively removed). The traces are then resorted and given
        new trace IDs.

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
                that should be less than ``match_tol`` in :attr:`par`
                used to find traces to merge.
            refit (:obj:`bool`, optional):
                If fitted models exist and traces are merged, run
                :func:`fit_refine` with its default arguments after
                the traces have been merged to refit the trace data.
            debug (:obj:`bool`, optional):
                Only passed to :func:`fit_refine` if ``refit`` is
                True.
        """
        if self.is_empty:
            msgs.warn('No traces to merge.')
        if self.edge_fit is None:
            msgs.error('Trace merging requires model fits to the trace location; run fit_refine.')
        _refit = refit
        if refit and self.edge_fit is None:
            msgs.warn('No previous fits existed, so fitting will not be redone.')
            _refit = False

        # Construct the bad pixel mask depending whether we matching
        # models or measurements
        bpm = np.tile(self.fully_masked_traces(flag=self.bitmask.bad_flags),(self.nspec,1))

        # The center values used to test for matches; can either be the
        # modeled or measured data.
        cen_match = np.ma.MaskedArray(self.edge_fit, mask=bpm)
        # The center values used to construct the merged trace; always
        # the measured values
        cen_merge = np.ma.MaskedArray(self.edge_cen, mask=self.edge_msk.astype(bool))
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
                self.edge_cen[gpm,indx[i]] = merged_trace[gpm]
                self.edge_msk[gpm,indx[i]] = 0

        if not np.any(rmtrace):
            msgs.info('No traces merged.')
            return

        # Remove traces and resort them
        self.remove_traces(rmtrace)

        # Erase any previous fitting, PCA, and slit-mask-match results
        # TODO: Should probably get sectioned off as a method because a
        # few methods do this.
        self.fittype = None
        self.edge_fit = None
        self.pcatype = None
        self.pca = None
        self.left_pca = None
        self.right_pca = None
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

    def good_traces(self, include_box=False, good_orders=False):
        r"""
        Select traces that are *not* fully masked.

        This is just a wrapper for :func:`fully_masked_traces`.

        Args:
            include_box (:obj:`bool`, optional):
                Include any identified box slits in the list of good
                traces.
            good_orders (:obj:`bool`, optional):
                Only include those traces that are well-matched to
                echelle orders. This is effectively ignored by
                multi-slit spectrographs because none of the order
                flags should have been tripped for them.

        Returns:
            `numpy.ndarray`_: Boolean array with shape :math:`(N_{\rm
            trace},)` flagging good traces.
        """
        if self.spectrograph.pypeline == 'Echelle' and good_orders and self.orderid is None:
            msgs.warn('Orders undefined! Selecting all traces. To select good orders only, first '
                      'run match_order().')
        bad_flags = self.bitmask.bad_flags
        exclude = self.bitmask.insert_flags
        if include_box:
            exclude += ['BOXSLIT']
        if good_orders:
            bad_flags += self.bitmask.order_flags
        else:
            exclude += self.bitmask.order_flags
        return np.logical_not(self.fully_masked_traces(flag=bad_flags, exclude=exclude))
        
    @property
    def is_synced(self):
        """
        Boolean that left-right traces are synchronized into slits.

        For this to be true:

            - **all** traces must be valid; i.e.,
              ``np.all(self.good_traces(include_box=True))`` must be
              True.

            - The sorted traces must be ordered
              left-right-left-right-...

        """
        if self.is_empty:
            return False
        if not np.all(self.good_traces(include_box=True)):
            return False
        side = np.clip(self.traceid, -1, 1)
        if len(side) == 0:
            return False
        return side[0] == -1 and side.size % 2 == 0 and np.all(side[1:] + side[:-1] == 0)

    def _masked_single_slit(self, trace_cen):
        """
        Handle masking a single slit as too short.

        This is largely a kludge that was necessary to address the
        features seen in the Keck_LRIS_blue long-slit data in the dev
        suite.
        
        .. todo::

            We need to understand how often a single slit is masked
            and if this is always the right way to deal with it.

        Args:
            trace_cen (`numpy.ndarray`_):
                Trace data to use for determining new edge locations.
        """
        if self.ntrace != 2:
            raise ValueError('Coding error:  Should only get here if there are two traces.')

        msgs.warn('The single slit found has been rejected because it is too short.  If this '
                  'was by mistake, re-run pypeit with a smaller `minimum_slit_length` parameter.'
                  '  Otherwise, we assume this is a long-slit with one edge off the detector '
                  'and with the current slit edges errantly isolating some feature in the data.')

        # TODO: May want to limit the number of columns included in this calculation.
        if np.mean(self.traceimg.image[:,int(np.ceil(np.max(trace_cen[:,1]))):]) \
                > np.mean(self.traceimg.image[:,:int(np.floor(np.min(trace_cen[:,0])))]):
            msgs.warn('The mean of the trace image to the right of the right trace is larger '
                      'than it is to the left of the left trace. Removing the right trace and '
                      're-synchronizing.')
            self.remove_traces(np.array([False,True]))
        else:
            msgs.warn('The mean of the trace image to the left of the left trace is larger than '
                      'it is to the right of the right trace. Removing the right trace and '
                      're-synchronizing.')
            self.remove_traces(np.array([True,False]))

        self.sync(rebuild_pca=False)

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
        `minimum_slit_gap`, `minimum_slit_length`,
        `minimum_slit_length_sci`, and `length_range`.

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
              masked as SHORTSLIT (see :class:`EdgeTraceBitMask`).
              The absolute tolerance is set using the platescale
              provided by the spectrograph class, the spatial binning
              (from :attr:`binning`), and the minimum slit length in
              arcsec (`minimum_slit_length` in :attr:`par`).
            - Traces that form a slit with a length (the median
              difference between the left and right edges) below an
              absolute tolerance (i.e., ``right-left < atol``) **for
              a science slit** are masked as either BOXSLIT or
              SHORTSLIT (see :class:`EdgeTraceBitMask`), depending on
              the value of ``minimum_slit_length``: if
              ``minimum_slit_length`` is None, all are flagged as
              BOXSLIT; otherwise, BOXSLITs are those that are larger
              than ``minimum_slit_length`` but smaller than
              ``minimum_slit_length_sci``. The absolute tolerance is
              set using the platescale provided by the spectrograph
              class, the spatial binning (from :attr:`binning`), and
              the minimum slit length in arcsec for the science slits
              (`minimum_slit_length_sci` in :attr:`par`).
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

        # Decide if the PCA should be rebuilt
        _rebuild_pca = rebuild_pca and self.pcatype is not None and self.can_pca()

        # Remove any fully masked traces and its synced counterpart;
        # force the removal of traces marked as SYNCERROR, even if
        # those traces were inserted.
        self.clean_traces(force_flag='SYNCERROR', rebuild_pca=_rebuild_pca, sync_mode='both',
                          assume_synced=True)

        # Use the fit data if available
        trace_cen = self.edge_cen if self.edge_fit is None else self.edge_fit

        # Flag trace locations falling off the detector. This general
        # check that is independent of whether or not the traces are
        # synced. However, this could mask traces based on the *fitted*
        # position, instead of the measured centroid position.
        # TODO: Keep a separate mask for the fit positions?
        indx = self.trace_pixels_off_detector(cen=trace_cen)
        if np.any(indx):
            self.edge_msk[indx] = self.bitmask.turn_on(self.edge_msk[indx], 'OFFDETECTOR')

        # Check the slits are synced
        if not self.is_synced:
            msgs.error('Edge traces are not yet (or improperly) synced.  Either sync() failed '
                       'or has not yet been executed.')

        # Parse parameters and report
        gap_atol = None
        length_atol = None
        length_atol_sci = None
        length_rtol = self.par['length_range']
        if self.par['minimum_slit_length'] is not None \
                or self.par['minimum_slit_length_sci'] is not None \
                or self.par['minimum_slit_gap'] is not None:
            platescale = parse.parse_binning(self.traceimg.detector.binning)[1] \
                            * self.traceimg.detector['platescale']
            msgs.info('Binning: {0}'.format(self.traceimg.detector.binning))
            msgs.info('Platescale per binned pixel: {0}'.format(platescale))
            if self.par['minimum_slit_length'] is not None:
                length_atol = self.par['minimum_slit_length']/platescale
            if self.par['minimum_slit_length_sci'] is not None:
                length_atol_sci = self.par['minimum_slit_length_sci']/platescale
            if self.par['minimum_slit_gap'] is not None:
                gap_atol = self.par['minimum_slit_gap']/platescale

        msgs.info('Minimum slit gap (binned pixels): {0}'.format(gap_atol))
        msgs.info('Minimum slit length (binned pixels): {0}'.format(length_atol))
        msgs.info('Minimum science slit length (binned pixels): {0}'.format(length_atol_sci))
        msgs.info('Range in slit length not limited' if length_rtol is None else
                  'Range in slit length limited to +/-{0:.1f}%'.format(length_rtol*100))

        # TODO: Should here and below only use the unmasked parts of
        # the trace for the slit gap and length computations?

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
                trace_cen = self.edge_cen if self.edge_fit is None else self.edge_fit

        # Calculate the slit length and gap
        slit_length = np.median(np.squeeze(np.diff(trace_cen.reshape(self.nspec,-1,2), axis=-1)),
                                axis=0)
        if length_atol is not None:
            # Find any short slits (flag both edges of the slit)
            indx = np.repeat(slit_length < length_atol, 2)
            if np.sum(indx) == self.ntrace:
                if self.ntrace == 2:
                    # TODO: I *really* don't like this because it has
                    # the potential to yield an infinite loop, but it's
                    # also the simplest approach.
                    return self._masked_single_slit(trace_cen)
                msgs.warn('All slits are too short!')
            if np.any(indx):
                msgs.info('Rejecting {0} slits that are too short.'.format(np.sum(indx)))
                self.edge_msk[:,indx] = self.bitmask.turn_on(self.edge_msk[:,indx], 'SHORTSLIT')

        if length_atol_sci is not None:
            # Find any box slits (flag both edges of the slit)
            indx = slit_length < length_atol_sci
            if length_atol is not None:
                indx &= (slit_length > length_atol)
            indx = np.repeat(indx, 2)
            if np.sum(indx) == self.ntrace:
                if self.ntrace == 2:
                    # TODO: I *really* don't like this because it has
                    # the potential to yield an infinite loop, but it's
                    # also the simplest approach.
                    return self._masked_single_slit(trace_cen)
                msgs.warn('All slits are too short!')
            if np.any(indx):
                msgs.info('Identifying/Rejecting {0} slits as box slits.'.format(np.sum(indx)))
                self.edge_msk[:,indx] = self.bitmask.turn_on(self.edge_msk[:,indx], 'BOXSLIT')

        if length_rtol is not None:
            # Find slits that are not within the provided fraction of
            # the median length
            indx = np.repeat(np.absolute(np.log(slit_length/np.median(slit_length)))
                             > np.log(1+length_rtol), 2)
            if np.any(indx):
                msgs.info('Rejecting {0} abnormally long or short slits.'.format(np.sum(indx)))
                self.edge_msk[:,indx] = self.bitmask.turn_on(self.edge_msk[:,indx], 'ABNORMALSLIT')

        # TODO: Check that slit edges meet a minimum slit gap?

        # Find all traces to remove
        rmtrace = self.fully_masked_traces(flag=self.bitmask.bad_flags, exclude=self.bitmask.exclude_flags)
        # Make sure to also remove the synced one
        rmtrace = self.synced_selection(rmtrace, mode='both', assume_synced=True)
        # Remove 'em
        self.remove_traces(rmtrace, rebuild_pca=_rebuild_pca)
        if self.is_empty:
            msgs.warn('Assuming a single long-slit and continuing.')
            self.bound_detector()

    def rm_user_traces(self, rm_traces):
        """
        Parse the user input traces to remove

        Args:
            rm_user_traces (list):
              y_spec, x_spat pairs

        Returns:

        """
        if not self.is_synced:
            msgs.error('Trace removal should only be executed after traces have been '
                       'synchronized into left-right slit pairs; run sync()')
        # Setup
        lefts = self.edge_fit[:, self.is_left]
        rights = self.edge_fit[:, self.is_right]
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
                # Mask
                self.bitmask.turn_on(self.edge_msk[:,indx], 'USERRMSLIT')

        # Syncronize the selection
        indx = self.synced_selection(indx, mode='both', assume_synced=True)
        # TODO: Add rebuild_pca ?
        # Remove
        self.remove_traces(indx)

    # TODO -- Add an option to distinguish between an actual remove and a flagging
    def remove_traces(self, indx, resort=True, rebuild_pca=False):
        r"""
        Remove a set of traces.

        This method directly alters the attributes containing the
        trace data.

        .. warning::

            If the traces are left-right synchronized, you should
            decide how you want to treat the traces as pairs. See
            :func:`synced_selection`. For example, to remove both
            traces in a pair if either are selected, run::

                edges.remove_traces(edges.synced_selection(indx, mode='both'))

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
            msgs.warn('No trace to remove.')
            return

        if np.all(indx):
            msgs.warn('All traces removed!')
            self._reinit_trace_data()
            return
            
        msgs.info('Removing {0} edge traces.'.format(np.sum(indx)))

        # Reset the trace data
        keep = np.logical_not(indx)
        self.edge_img = self.edge_img[:,keep]
        self.edge_cen = self.edge_cen[:,keep]
        self.edge_err = self.edge_err[:,keep]
        self.edge_msk = self.edge_msk[:,keep]
        if self.edge_fit is not None:
            self.edge_fit = self.edge_fit[:,keep]
        self.traceid = self.traceid[keep]

        if resort:
            # Resort by the spatial dimension
            self.spatial_sort()

        # Reset the PCA
        self._reset_pca(rebuild_pca and self.pcatype is not None and self.can_pca())

    def synced_selection(self, indx, mode='ignore', assume_synced=False):
        r"""
        Coordinate selection of traces based on their synchronization
        status.

        Args:
            indx (`numpy.ndarray`_):
                Boolean vector selecting a set of traces. The shape
                should be :math:`(N_{\rm trace},)`.
            mode (:obj:`str`, optional):
                If the traces are left-right synchronized (see
                :func:`sync` and :func:`is_synced`), this method
                dictates how the selection of traces should be
                changed to deal with edges that have been
                synchronized into left-right pairs. Methods are:

                    - ``'ignore'``: Ignore the synchronization and
                      just return the input.
                    - ``'both'``: If one of the traces in a pair is
                      selected, select both.
                    - ``'neither'``: If only one of the traces in a
                      pair is selected, *deselect* both.
            assume_synced (:obj:`bool`, optional):
                Assume the set of traces is synced when identifying
                which traces to select/de-select if mode is either
                ``'both'`` or ``'neither'``. This essentially ignores
                the sync status of the traces and will fault if the
                number of traces is not even. If False, the sync
                status is checked using :func:`is_synced`.

        Raises:
            PypeItError:
                Raised if the input ``indx`` doesn't have the correct
                length, if mode is *not* ``'ignore'`` and the edges
                have not been left-right synchronized and
                ``assume_synced`` is False, or if the provided
                ``mode`` is not valid.

        Returns:
            `numpy.ndarray`_: The updated boolean vector selecting
            traces based on the input mode and the synchronization
            status.
        """
        if mode == 'ignore':
            return indx

        if indx.size != self.ntrace:
            msgs.error('Boolean array selecting traces to remove has incorrect length.')

        if not assume_synced and not self.is_synced:
            msgs.error('To synchronize the trace selection, it is expected that the traces have '
                       'been left-right synchronized.  Either run sync() to sychronize, ignore '
                       'the synchronization (which may raise an exception) by setting '
                       'assume_synced=True.')
        if mode == 'both':
            return np.repeat(np.any(indx.reshape(-1,2), axis=1), 2)
        elif mode == 'neither':
            return np.repeat(np.all(indx.reshape(-1,2), axis=1), 2)
        msgs.error('Unknown synchronized trace selection mode: {0}'.format(mode))

    def clean_traces(self, force_flag=None, rebuild_pca=True, sync_mode='ignore',
                     assume_synced=False):
        """
        Remove any traces that are fully masked as bad.

        Traces selected for removal must be fully masked; see
        :func:`fully_masked_traces`.

        By default, flags used to select bad trace measurements are
        those provided by :func:`EdgeTraceBitMask.bad_flags`, and
        those flags excluded from designating a trace as bad are
        provided by :func:`EdgeTraceBitMask.exclude_flags`. To force
        removal of traces with certain flags, regardless of these two
        groups, use ``force_flag``. See :class:`EdgeTraceBitMask` for
        list of flags.

        Args:
            force_flag (:obj:`str`, :obj:`list`, optional):
                Force inclusion of these flags in assessing whether
                or not a trace is fully masked.
            rebuild_pca (:obj:`bool`, optional):
                Rebuild the PCA decomposition of the traces based on
                the loaded trace data. Passed directly to
                :func:`remove_traces`, which is only called if traces
                are in fact removed.
            sync_mode (:obj:`str`, optional):
                If the traces are left-right synchronized (see
                :func:`sync` and :func:`is_synced`), use this method
                to deal with edges paired with those to be removed.
                See :func:`synced_selection`.
            assume_synced (:obj:`bool`, optional):
                Assume the set of traces is synced. See
                :func:`synced_selection`.
        """
        if self.is_empty:
            msgs.warn('No traces to clean.')
            return

        # Traces to remove
        rmtrace = self.fully_masked_traces(flag=self.bitmask.bad_flags,
                                           exclude=self.bitmask.exclude_flags)
        if force_flag is not None:
            rmtrace |= self.fully_masked_traces(flag=force_flag)

        if np.any(rmtrace):
            # The removed traces should not have been included in the
            # PCA decomposition to begin with, but this call to remove
            # traces has to "rebuild" the PCA because it will remove it
            # otherwise. Note that `remove_traces` only rebuilds the
            # PCA if the parameter passed to it is true AND :attr:`pca`
            # previously existed AND the traces can be pca'd to begin
            # with.
            rmtrace = self.synced_selection(rmtrace, mode=sync_mode, assume_synced=assume_synced)
            self.remove_traces(rmtrace, rebuild_pca=rebuild_pca)

    def spatial_sort(self, use_mean=False, use_fit=True):
        """
        Sort the traces spatially.

        The coordinates used for the sorting are either the measured
        centroids or the fitted parameterization (see ``use_fit``). The
        fiducial coordinates that are sorted are either the mean of
        the unmasked coordinates over all spectral rows or the
        unmasked coordinates at a specified reference spectral row
        (see ``use_mean``).

        The trace IDs are also reassigned to be sorted spatially;
        i.e., the trace IDs for three synced slits would be ``[-1, 1,
        -2, 2, -3, 3]``.
        
        All attributes are edited in-place.

        Args:
            use_mean (:obj:`bool`, optional):
                Sort according to the mean of the masked spatial
                positions. If False, the spatial position at a
                reference spectral row is used, where the reference
                spectral row is either the same as used by the PCA
                (if available) or the result of
                :func:`~pypeit.core.trace.most_common_trace_row`
                using the current trace mask.
            use_fit (:obj:`bool`, optional):
                Sort according to the fit positions instead of the
                measured positions. Otherwise, only use the fit
                positions if they're available and the measured
                location is masked.
        """
        if self.is_empty:
            msgs.error('No traces to sort.')

        # Check input
        if use_fit and self.edge_fit is None:
            msgs.warn('Fit data is not available; cannot use it for spatially sorting the edges.')

        # Set up the coordinates to use
        bpm = self.bitmask.flagged(self.edge_msk, self.bitmask.bad_flags)
        cen = self.edge_fit.copy() if self.edge_fit is not None and use_fit \
                    else self.edge_cen.copy()
        # Replace masked values with the fit if it is available
        if self.edge_fit is not None and not use_fit:
                cen[bpm] = self.edge_fit[bpm]

        # Get the sorted indices
        if use_mean:
            # Sort the traces by their spatial position (always use
            # measured positions even if fit positions are available)
            srt = np.argsort(np.mean(cen, axis=0))
        else:
            # Sort according to the spatial position in one row
            reference_row = trace.most_common_trace_row(bpm) if self.pcatype is None \
                                else (self.left_pca.reference_row if self.par['left_right_pca']
                                    else self.pca.reference_row)
            msgs.info('Re-sorting edges based on where they cross row {0}'.format(reference_row))
            srt = np.argsort(cen[reference_row,:])

        # Resort the arrays
        self.traceid = self.traceid[srt]
        self.edge_img = self.edge_img[:,srt]
        self.edge_cen = self.edge_cen[:,srt]
        self.edge_err = self.edge_err[:,srt]
        self.edge_msk = self.edge_msk[:,srt]
        if self.edge_fit is not None:
            self.edge_fit = self.edge_fit[:,srt]

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
            return self.build_pca(use_center=self.pcatype == 'center') 
        # Remove the existing PCA
        self.pcatype = None
        self.pca = None
        self.left_pca = None
        self.right_pca = None

    def current_trace_locations(self):
        """
        Return an image with the trace IDs at the locations of each
        edge in the original image.
        """
        edge_img = np.zeros((self.nspec, self.nspat), dtype=int)
        if self.is_empty:
            return edge_img
        i = np.tile(np.arange(self.nspec), (self.ntrace,1)).T.ravel()
        edge_img[i, self.edge_img.ravel()] = np.tile(self.traceid, (self.nspec,1)).ravel()
        return edge_img

    def fit_refine(self, weighting='uniform', debug=False, idx=None):
        """
        Iteratively re-measure and fit a functional form to the edge
        locations.

        Before fitting the traces, :func:`check_traces` is used to
        flag traces that do not meet a minimum length required for
        fitting (set by ``fit_min_spec_length`` in :attr:`par`).

        After this, the function is primarily a wrapper for
        :func:`~pypeit.core.trace.fit_trace`, run once per edge side
        (left and right). Both the measured centers and the fitted
        model are modified by :func:`~pypeit.core.trace.fit_trace`,
        such that :attr:`edge_cen`, :attr:`edge_err`,
        :attr:`edge_msk`, :attr:`edge_fit`, :attr:`edge_fit_type`, and
        :attr:`edge_img` are all modified by this method. Note that
        only traces that are *not* fully flagged are fit.

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are
        ``max_shift_abs``, ``max_spat_error``, ``fit_function``,
        ``fit_order``, ``fwhm_uniform``, ``fwhm_gaussian``,
        ``fit_maxdev``, ``fit_maxiter``, ``fit_niter``, and
        ``fit_min_spec_length``.

        Args:
            weighting (:obj:`str`, optional):
                The weighting to apply to the position within each
                integration window (see
                :func:`~pypeit.core.trace.fit_trace`).
            debug (:obj:`bool`, optional):
                Run :func:`~pypeit.core.trace.fit_trace` in debug
                mode.
            idx (`numpy.ndarray`_, optional):
                Array of strings with the IDs for each object. Used
                only if ``debug`` is true for the plotting (see
                :func:`~pypeit.core.trace.fit_trace`).
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

        # Check the traces to make sure they meet the minimum length.
        # This modifies self.edge_msk directly.
        self.check_traces(minimum_spec_length=minimum_spec_length)

        # Generate bogus ivar and mask once here so that they don't
        # have to be generated multiple times.
        # TODO: Keep these as work space as class attributes?
        ivar = np.ones_like(self.sobelsig, dtype=float)
        bpm = np.zeros_like(self.sobelsig, dtype=bool) if self.tracebpm is None else self.tracebpm

        # Initialize arrays
        fit = np.zeros_like(self.edge_cen, dtype=float)
        cen = np.zeros_like(self.edge_cen, dtype=float)
        err = np.zeros_like(self.edge_cen, dtype=float)
        msk = np.zeros_like(self.edge_cen, dtype=self.bitmask.minimum_dtype())

        # Flag bad traces; this explicitly does *not* exclude inserted traces
        edge_bpm = self.bitmask.flagged(self.edge_msk, flag=self.bitmask.bad_flags)

        # Fit both sides
        for side in ['left', 'right']:
            # Get the image relevant to tracing this side
            _sobelsig = self._side_dependent_sobel(side)
            # Select traces on this side and that are not fully masked
            indx = (self.is_left if side == 'left' else self.is_right) \
                        & np.invert(np.all(edge_bpm, axis=0))
            if not np.any(indx):
                continue

            # Perform the fit
            fit[:,indx], cen[:,indx], err[:,indx], msk[:,indx], _ \
                    = trace.fit_trace(_sobelsig, self.edge_cen[:,indx], order, ivar=ivar,
                                      bpm=bpm, trace_bpm=edge_bpm[:,indx],
                                      weighting=weighting, fwhm=fwhm, maxshift=maxshift,
                                      maxerror=maxerror, function=function, maxdev=maxdev,
                                      maxiter=maxiter, niter=niter, bitmask=self.bitmask,
                                      debug=debug, idx=idx, xmin=xmin, xmax=xmax)

        # Save the results of the edge measurements
        self.edge_cen = cen
        self.edge_err = err
        # Add the new masking information; this merges the new and
        # existing mask because the existing mask was used by fit_trace
        # to ignore input trace data
        self.edge_msk |= msk
        # Save the model fits
        self.edge_fit = fit
        self.fittype = '{0} : order={1}'.format(function, order)
        # Set the pixelated trace data based on the fit, instead of the
        # measured centroids
        self.edge_img = np.round(self.edge_fit).astype(int)
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
            :attr:`edge_msk` will can be altered by this call.

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
        good = np.sum(np.invert(self.bitmask.flagged(self.edge_msk)), axis=0) > 0

        # Returned value depends on whether or not the left and right
        # traces are done separately
        return np.sum(good[self.is_left]) > self.par['pca_min_edges'] \
                    and np.sum(good[self.is_right]) > self.par['pca_min_edges'] \
                    if self.par['left_right_pca'] else np.sum(good) > self.par['pca_min_edges']

    def predict_traces(self, edge_cen, side=None):
        """
        Use the PCA decomposition to predict traces.

        The PCA decomposition must be available via :attr:`pca` or
        :attr:`left_pca` and :attr:`right_pca`; see
        :func:`build_pca`. This is a convenience method to handle the
        PCA predictions given that left and right traces can be
        decomposed separately or simultaneously.

        Args:
            edge_cen (:obj:`float`, `numpy.ndarray`):
                A single value or 1D array with the spatial location
                (column) for 1 or more traces to predict. The
                predicted traces will pass through these spatial
                positions (columns) and the reference spectral row
                set for the PCA decomposition; see :func:`build_pca`.
            side (:obj:`float`, `numpy.ndarray`, optional):
                A single value or 1D integer array indicating the
                edge side to be predicted; -1 for left and 1 for
                right. Must be the same length as `edge_cen`. This is
                only used if the PCA is side-dependent, and *must* be
                provided in the case that it is (see `left_right_pca`
                in :attr:`par`).

        Returns:
            `numpy.ndarray`: A 1D or 2D array of size :attr:`nspec`
            by the length of the position array provided. If a single
            coordinate is provided, a single trace vector is
            returned.
        """
        if self.pcatype is None:
            msgs.error('Must first run the PCA analysis fo the traces; run build_pca.')

        _edge_cen = np.atleast_1d(edge_cen)
        _side = np.atleast_1d(side)
        if _edge_cen.size != _side.size:
            msgs.error('Spatial locations and side integers must have the same shape.')

        if self.par['left_right_pca']:
            trace_add = np.zeros((self.nspec,0), dtype='float')
            for s,p in zip([-1,1], [self.left_pca,self.right_pca]):
                indx = _side == s
                if not np.any(indx):
                    continue
                trace_add = np.hstack((trace_add, p.predict(np.atleast_1d(_edge_cen[indx]))))
        else:
            trace_add = self.pca.predict(_edge_cen)

        return trace_add if isinstance(edge_cen, np.ndarray) else trace_add.ravel()

    def build_pca(self, use_center=False, debug=False):
        """
        Build a PCA model of the current edge data.

        Primarily a wrapper that instantiates :attr:`pca`, or
        :attr:`left_pca` and :attr:`right_pca` if left and right
        traces are analyzed separately. All of these objects will be
        a :class:`pypeit.tracepca.TracePCA`, if instantiated. After
        executing this method, traces can be predicted by
        :func:`pypeit.tracepca.TracePCA.predict` for the relevant PCA
        decomposition; see :func:`predict_traces`.

        If no parametrized function has been fit to the trace data or
        if specifically requested (see `use_center`), the PCA is
        based on the measured trace centroids (:attr:`edge_cen`);
        othwerwise, the PCA uses the parametrized trace fits
        (:attr:`edge_fit`).

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
        if self.pcatype is not None:
            msgs.warn('PCA model already exists and will be overwritten.')
        if self.edge_fit is None and not use_center:
            msgs.warn('No trace fits exits.  PCA based on trace centroid measurements.')

        # Check if the PCA decomposition can be performed
        if not self.can_pca():
            msgs.error('Traces do not meet necessary criteria for the PCA decomposition.')

        # Set the data used to construct the PCA
        self.pcatype = 'center' if self.edge_fit is None or use_center else 'fit'

        # When constructing the PCA, ignore bad trace measurements
        # *and* any traces inserted by hand.
        bpm = self.bitmask.flagged(self.edge_msk)

        # TODO: Is there a way to propagate the mask to the PCA?
        # TODO: Keep a separate mask specifically for the fit data? e.g., edge_fit_msk

        # The call to can_pca means that short traces are fully masked
        # and that valid traces will be any trace with unmasked pixels.
        use_trace = np.sum(np.invert(bpm), axis=0) > 0

        # Set the reference row so that, regardless of whether the PCA
        # is for the left, right, or all traces, the reference row is
        # always the same.
        reference_row = trace.most_common_trace_row(bpm[:,use_trace])

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
        _pca = [None]*len(pcaindx)
        for i,indx in enumerate(pcaindx):
            # Grab the trace data. This uses all the data, even if some
            # of it is masked.
            trace_inp = self.edge_cen[:,indx] if self.edge_fit is None or use_center \
                            else self.edge_fit[:,indx]

            # Instantiate the PCA
            _pca[i] = TracePCA(trace_inp, npca=npca, pca_explained_var=pca_explained_var,
                                   reference_row=reference_row)

            # Set the order of the function fit to the PCA
            # coefficiencts: Order is set to cascade down to lower
            # order for components that account for a smaller
            # percentage of the variance.
            _order = np.clip(order - np.arange(_pca[i].npca), 1, None).astype(int)
            msgs.info('Order of function fit to each component: {0}'.format(_order))

            # Apply a 10% relative error to each coefficient. This
            # performs better than use_mad, since larger coefficients
            # will always be considered inliers, if the coefficients
            # vary rapidly with order as they sometimes do.
            #ivar = utils.inverse(np.square(np.fmax(0.1*np.abs(_pca[i].pca_coeffs), 0.1)))
            #ivar = None

            # TODO: Instead, weight by the mean/median value of
            # sobel_sig along each trace.

            # Run the fit
            _pca[i].build_interpolator(_order, function=function, lower=lower, upper=upper,
                                       minx=0., maxx=self.nspat-1., maxrej=maxrej, maxiter=maxiter,
                                       debug=debug)

            # TODO: Use the rejected pca coefficiencts to reject traces?

        # Save the result
        if left_right_pca:
            self.left_pca, self.right_pca = _pca
        else:
            self.pca = _pca[0]

    def pca_refine(self, use_center=False, debug=False, force=False):
        """
        Use a PCA decomposition to refine the traces.

        If no parametrized function has been fit to the trace data or
        if specifically requested (see ``use_center``), the PCA is
        based on the measured trace centroids (:attr:`edge_cen`);
        othwerwise, the PCA uses the parametrized trace fits
        (:attr:`edge_fit`).

        If needed or forced to, this first executes :func:`build_pca`
        and then uses :func:`predict_traces` to use the PCA to reset
        the trace data.

        Only used parameter from :attr:`par`
        (:class:`~pypeit.par.pypeitpar.EdgeTracePar`) is
        ``left_right_pca``.

        Args:
            use_center (:obj:`bool`, optional):
                Use the center measurements for the PCA decomposition
                instead of the functional fit to those data. This is
                only relevant if both are available. If no fits have
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
        _pcatype = 'center' if use_center or self.edge_fit is None else 'fit'
        if force or self.pcatype is None or self.pcatype != _pcatype:
            self.build_pca(use_center=use_center, debug=debug)

        self.fittype = 'pca'

        # Predict the traces using the PCA
        reference_row = self.left_pca.reference_row if self.par['left_right_pca'] \
                            else self.pca.reference_row
        trace_ref = self.edge_cen[reference_row,:] if self.pcatype == 'center' \
                            else self.edge_fit[reference_row,:]
        side = self.is_right.astype(int)*2-1
        self.edge_fit = self.predict_traces(trace_ref, side)

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
        a wrapper for :func:`~pypeit.core.trace.peak_trace`. See the
        documentation of that function for the explanation of the
        algorithm.

        If the left and right traces have separate PCA
        decompositions, this function makes one call to
        :func:`~pypeit.core.trace.peak_trace` for each side.
        Otherwise, a single call is made to
        :func:`~pypeit.core.trace.peak_trace` where both the peak and
        troughs in :attr:`sobelsig` are detected and traced.

        Note that this effectively reinstantiates much of the object
        attributes, including :attr:`traceid` :attr:`edge_cen`
        :attr:`edge_err` :attr:`edge_msk` :attr:`edge_img`
        :attr:`edge_fit`, and :attr:`fittype`.

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are
        ``left_right_pca``, ``edge_thresh``, ``smash_range``,
        ``edge_detect_clip``, ``trace_median_frac``, ``trace_thresh``,
        ``fit_function``, ``fit_order``, ``fwhm_uniform``, ``fwhm_uniform``,
        ``niter_gaussian``, ``niter_gaussian``, ``fit_maxdev``, and
        ``fit_maxiter``.

        Args:
            rebuild_pca (:obj:`bool`, optional):
                This method fundamentally resets the trace data,
                meaning that the PCA is no longer valid. Use this
                boolean to have the method rebuild the PCA based on
                the refined traces. Note that the PCA is *not* then
                used to reset the fitted trace data; i.e.,
                :attr:`edge_fit` remains based on the output of
                :func:`~pypeit.core.trace.peak_trace`.
            debug (:obj:`bool`, optional):
                Run in debug mode.

        Raises:
            PypeItError:
                Raised if :attr:`pca` is not defined.

        """
        # Check that there are traces to refine!
        if self.is_empty:
            msgs.error('No traces are defined.')

        if self.pcatype is None:
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
        ivar = np.ones_like(self.sobelsig, dtype=float)
        bpm = np.zeros_like(self.sobelsig, dtype=bool) if self.tracebpm is None else self.tracebpm

        # Treatment is different if the PCA was done for all traces or
        # separately for left and right traces
        if self.par['left_right_pca']:
            # Initialize the arrays holding the results for both sides
            fit = np.zeros((self.nspec,0), dtype='float')
            cen = np.zeros((self.nspec,0), dtype='float')
            err = np.zeros((self.nspec,0), dtype='float')
            msk = np.zeros((self.nspec,0), dtype=self.bitmask.minimum_dtype())

            # Iterate through each side
            for side in ['left', 'right']:
                # Get the image relevant to tracing
#                _sobelsig = trace.prepare_sobel_for_trace(self.sobelsig, bpm=self.bpm, boxcar=5,
#                                                          side=side)
                _sobelsig = self._side_dependent_sobel(side)
                _pca = self.left_pca if side == 'left' else self.right_pca

                _fit, _cen, _err, _msk, nside \
                        = trace.peak_trace(_sobelsig, ivar=ivar, bpm=bpm,
                                           trace_map=_pca.predict(np.arange(self.nspat)),
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
            _sobelsig = trace.prepare_sobel_for_trace(self.sobelsig, bpm=self.tracebpm, boxcar=5,
                                                      side=None)

            # Find and trace both peaks and troughs in the image. The
            # input trace data (`trace` argument) is the PCA prediction
            # of the trace that passes through each spatial position at
            # the reference spectral pixel.
            fit, cen, err, msk, nleft \
                    = trace.peak_trace(_sobelsig, ivar=ivar, bpm=bpm,
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

        # Reset the trace data
        self.traceid = np.zeros(ntrace, dtype=int)
        self.traceid[:nleft] = -1-np.arange(nleft)
        self.traceid[nleft:] = 1+np.arange(ntrace-nleft)
        self.edge_fit = fit
        self.fittype = '{0} : order={1}'.format(function, order)
        self.edge_cen = cen
        self.edge_err = err
        # The mask is entirely new and shouldn't be merged with the
        # existing mask.
        self.edge_msk = msk
        self.edge_img = np.round(self.edge_fit).astype(int)

        # Spatially sort the traces
        self.spatial_sort()
        # Reset the PCA
        self._reset_pca(rebuild_pca and self.can_pca())
        self.log += [inspect.stack()[0][3]]

    # TODO: Make this a core function?
    def _get_insert_locations(self):
        """
        Find where edges need to be inserted.

        This only determines where the left-right ordering of the
        traces implies that a trace needs to be inserted. Where the
        trace is inserted and with what shape is determined by
        :func:`_get_reference_locations` and :func:`sync`,
        respectively.

        Returns:
            :obj:`tuple`: Three `numpy.ndarray`_ objects are
            returned:

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
        Determine the reference locations for traces to add during
        the left-right synchronization.

        The positions of new traces are determined either by the
        median slit length of the existing left-right pairs, or based
        on the slit length of the nearest left-right pair. This
        function only determines the positions for the new traces at
        the reference spetral location (row). The shape is determined
        by :func:`sync`.

        Used parameters from :attr:`par`
        (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are
        ``sync_center``, ``sync_to_edge``, and ``gap_offset``.

        Args:
            trace_cen (`numpy.ndarray`_):
                Trace data to use for determining new edge locations.
            add_edge (`numpy.ndarray`_):
                Boolean array indicating that a trace in the new
                array is an added trace. The number of False entries
                in ``add_edge`` should match the length of the 2nd
                axis of ``trace_cen``.

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
        bpm = self.bitmask.flagged(self.edge_msk, flag=self.bitmask.bad_flags)
        reference_row = trace.most_common_trace_row(bpm) if self.pcatype is None \
                            else (self.left_pca.reference_row if self.par['left_right_pca']
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

        if np.all(np.logical_not(np.absolute(offset) > 0)):
            # No offsets necessary
            return trace_cen

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
        Match left and right edge traces to construct slit edge
        pairs.

        Synchronization of the slits proceeds as follows:

            - Any fully masked traces are removed using
              :func:`clean_traces`.
            - If that operation removes *all* traces (see
              :func:`is_empty`), two traces are added at the edge of
              the detector.
            - The spatial order of the traces is sorted (see
              :func:`spatial_sort`).
            - At this point, if the traces are already synced (see
              :func:`is_synced`), the synchronization is checked using
              :func:`check_synced` and the method finishes.
            - If the traces need to be synced, the method determines
              which traces need an inserted trace to make a pair (see
              :func:`_get_insert_locations`).
            - The reference locations for the new traces are
              determined by :func:`_get_reference_locations`.
            - The shape of the inserted traces is set according to
              the ``sync_predict`` parameter in :attr:`par`. The
              shape is either predicted by the PCA decomposition or
              taken to be exactly the same shape as the nearest left
              or right edge.
            - These new traces are then inserted using
              :func:`insert_traces` and flagged as being inserted by
              the synchronization process.
            - Assuming there wasn't an error in the insertion scheme,
              the synchronized traces are then checked by
              :func:`check_synced` and the method finishes.

        Before the last step above, the synchronization is checked.
        Certain corner cases can lead to catastrophic errors in where
        the inserted traces are placed such that the left-right
        ordering of the synchronized traces is incorrect and
        exception is raised.
            
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
        # Remove any fully masked traces. Keeps any inserted or
        # box-slit traces.
        self.clean_traces(rebuild_pca=rebuild_pca)

        # Make sure there are still traces left
        if self.is_empty:
            msgs.warn('No traces left!  Left and right edges placed at detector boundaries.')
            self.bound_detector()

        # Make sure that the traces are sorted spatially
        self.spatial_sort()

        # If the traces are already synced, check them and log the
        # function as completed
        if self.is_synced:
            self.check_synced(rebuild_pca=rebuild_pca and self.pcatype is not None)
            self.log += [inspect.stack()[0][3]]
            return

        # Edges are currently not synced, so check the input
        if self.par['sync_predict'] not in ['pca', 'nearest']:
            msgs.error('Unknown trace mode: {0}'.format(self.par['sync_predict']))
        if self.par['sync_predict'] == 'pca' and self.pcatype is None:
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
        trace_cen = self.edge_cen if self.edge_fit is None else self.edge_fit

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
            self.edge_msk[:,indx] = self.bitmask.turn_on(self.edge_msk[:,indx], 'SYNCERROR')

        if debug:
            msgs.info('Show instance includes inserted traces but before checking the sync.')
            self.show(flag='any')

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

    def insert_traces(self, side, trace_cen, loc=None, mode='user', resort=True, nudge=True):
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
            nudge (:obj:`bool`, optional):
                Allow the traces to be nudged away from the detector
                edge according to :attr:`par` and
                :func:`nudge_traces`.

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
        if nudge:
            _trace_cen = self.nudge_traces(_trace_cen)

        # Set the mask
        mask = np.zeros(_trace_cen.shape, dtype=self.bitmask.minimum_dtype())
        # Flag the traces pixels that fall off the detector
        indx = self.trace_pixels_off_detector(cen=_trace_cen)
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
            self.edge_img = np.round(_trace_cen).astype(int)
            self.edge_cen = _trace_cen
            self.edge_err = np.zeros(_trace_cen.shape, dtype=float)
            self.edge_msk = mask
            self.edge_fit = _trace_cen
            return

        # Add the new traces. The new traces are added to both the
        # fitted list and the center list!
        self.traceid = np.insert(self.traceid, loc, _traceid)
        self.edge_img = np.insert(self.edge_img, loc, np.round(_trace_cen).astype(int), axis=1)
        self.edge_cen = np.insert(self.edge_cen, loc, _trace_cen, axis=1)
        self.edge_err = np.insert(self.edge_err, loc,
                                  np.zeros(_trace_cen.shape, dtype=float), axis=1)
        self.edge_msk = np.insert(self.edge_msk, loc, mask, axis=1)
        self.edge_fit = np.insert(self.edge_fit, loc, _trace_cen, axis=1)

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
        bpm = np.all(self.bitmask.flagged(self.edge_msk, flag=flag), axis=0)
        if exclude is not None:
            bpm &= np.logical_not(np.any(self.bitmask.flagged(self.edge_msk, flag=exclude), axis=0))
        return bpm
    
    def maskdesign_matching(self, debug=False):
        """
        Match slit info from the mask design data to the traced slits.

        Use of this method requires:
            - a PCA decomposition is available,
            - :attr:`spectrograph` has a viable `get_slitmask` method
              to read slit mask design data. This data can be pulled
              from one of the files used to construct the trace image.
            - :attr:`spectrograph` has a viable `get_grating` method
              which provides the grating info to recover the optical model.
            - :attr:`spectrograph` has a viable `get_amapbmap` method
              which provides pre- and post-grating maps of the detector
              used convert the mask design data from mm to pixels.

        The method use a collection of scripts in pypeit.core.slitdesign_matching
        which are taken from DEEP2 IDL-based pipeline for DEIMOS data.


        Args:
            debug (:obj:`bool`, optional):
                Run in debug mode.
        """

        # Check that there are traces to match!
        if self.is_empty:
            msgs.error('No traces to match.')

        # The PCA decomposition must have already been determined
        if self.pcatype is None:
            msgs.error('Must first run the PCA analysis for the traces; run build_pca.')

        # `traceimg` must have knowledge of the flat frame that built it
        if self.spectrograph.get_slitmask(self.traceimg.files[0]) is None:
            msgs.error('Unable to read slitmask design info')
        if self.spectrograph.get_grating(self.traceimg.files[0]) is None:
            msgs.error('Unable to read grating info')
        if self.spectrograph.get_amapbmap(self.traceimg.files[0]) is None:
            msgs.error('Unable to read amap and bmap')

        # Match left and right edges separately
        # Sort slits in mm from the slit-mask design
        sortindx = np.argsort(self.spectrograph.slitmask.center[:, 0])

        # Left (bottom) and right (top) traces in pixels from optical model (image plane and detector)
        # bottom
        omodel_bcoo = self.spectrograph.mask_to_pixel_coordinates(x=self.spectrograph.slitmask.bottom[:, 0],
                                                                 y=self.spectrograph.slitmask.bottom[:, 1])
        bedge_img, ccd_b, bedge_pix = omodel_bcoo[0], omodel_bcoo[2], omodel_bcoo[3]

        # top
        omodel_tcoo = self.spectrograph.mask_to_pixel_coordinates(x=self.spectrograph.slitmask.top[:, 0],
                                                                  y=self.spectrograph.slitmask.top[:, 1])
        tedge_img, ccd_t, tedge_pix = omodel_tcoo[0], omodel_tcoo[2], omodel_tcoo[3]

        # Per each slit we take the median value of the traces over the wavelength direction. These medians will be used
        # for the cross-correlation with the traces found in the images.
        ccdnum=self.traceimg.detector.det
        omodel_bspat = np.zeros(self.spectrograph.slitmask.nslits)
        omodel_tspat = np.zeros(self.spectrograph.slitmask.nslits)

        for i in range(omodel_bspat.size):
            # We "flag" the left and right traces predicted by the optical model that are outside of the
            # current detector, by giving a value of -1.
            # bottom
            omodel_bspat[i] = -1 if bedge_pix[i, ccd_b[i, :] == ccdnum].shape[0] < 10 else np.median(
                                                                                    bedge_pix[i, ccd_b[i, :] == ccdnum])
            # top
            omodel_tspat[i] = -1 if tedge_pix[i, ccd_t[i, :] == ccdnum].shape[0] < 10 else np.median(
                                                                                    tedge_pix[i, ccd_t[i, :] == ccdnum])

            # If a left (or right) trace is outside of the detector, the corresponding right (or left) trace
            # is determined using the pixel position from the image plane.
            whgood = np.where(tedge_img[i, :] > -1e4)[0]
            npt_img = whgood.shape[0] // 2
            # This is hard-coded for DEIMOS, since it refers to the detectors configuration
            whgood = whgood[:npt_img] if ccdnum <= 4 else whgood[npt_img:]
            if omodel_bspat[i] == -1 and omodel_tspat[i] >= 0:
                omodel_bspat[i] = omodel_tspat[i] - np.median((tedge_img - bedge_img)[i, whgood])
            if omodel_tspat[i] == -1 and omodel_bspat[i] >= 0:
                omodel_tspat[i] = omodel_bspat[i] + np.median((tedge_img - bedge_img)[i, whgood])

            # If the `omodel_bspat` is greater than `omodel_tspat` we switch the order
            if omodel_bspat[i] > omodel_tspat[i]:
                invert_order = omodel_bspat[i]
                omodel_bspat[i] = omodel_tspat[i]
                omodel_tspat[i] = invert_order

        # This print a QA table with info on the slits (sorted from left to right) that fall in the current detector.
        # The only info provided here is `slitid`, which is called `dSlitId` in the DEIMOS design file. I had to remove
        # `slitindex` because not always matches `SlitName` from the DEIMOS design file.
        if not debug:
            num = 0
            msgs.info('*' * 18)
            msgs.info('{0:^6s} {1:^12s}'.format('N.', 'dSlitId'))
            msgs.info('{0:^6s} {1:^12s}'.format('-' * 5, '-' * 9))
            for i in range(sortindx.shape[0]):
                if omodel_bspat[sortindx][i] != -1 or omodel_tspat[sortindx][i] != -1:
                    msgs.info('{0:^6d} {1:^12d}'.format(num, self.spectrograph.slitmask.slitid[sortindx][i]))
                    num += 1
            msgs.info('*' * 18)

        # If instead we run this method in debug mode, we print more info useful for comparison, for example, with
        # the IDL-based pipeline.
        if debug:
            num = 0
            msgs.info('*' * 92)
            msgs.info('{0:^5s} {1:^10s} {2:^12s} {3:^12s} {4:^14s} {5:^16s} {6:^16s}'.format('N.',
                                                                                             'dSlitId', 'slitLen(mm)',
                                                                                             'slitWid(mm)',
                                                                                             'spat_cen(mm)',
                                                                                             'omodel_bottom(pix)',
                                                                                             'omodel_top(pix)'))
            msgs.info('{0:^5s} {1:^10s} {2:^12s} {3:^12s} {4:^14s} {5:^16s} {6:^14s}'.format('-' * 4, '-' * 9, '-' * 11,
                                                                                             '-' * 11, '-' * 13,
                                                                                             '-' * 18, '-' * 15))
            for i in range(sortindx.shape[0]):
                if omodel_bspat[sortindx][i] != -1 or omodel_tspat[sortindx][i] != -1:
                    msgs.info('{0:^5d}{1:^14d} {2:^9.3f} {3:^12.3f} {4:^14.3f}    {5:^16.2f} {6:^14.2f}'
                              .format(num, self.spectrograph.slitmask.slitid[sortindx][i],
                                         self.spectrograph.slitmask.length[sortindx][i],
                                         self.spectrograph.slitmask.width[sortindx][i],
                                         self.spectrograph.slitmask.center[:, 0][sortindx][i],
                                         omodel_bspat[sortindx][i], omodel_tspat[sortindx][i]))
                    num += 1
            msgs.info('*' * 92)

        reference_row = self.left_pca.reference_row if self.par['left_right_pca'] else self.pca.reference_row
        spat_bedge = self.edge_fit[reference_row, self.is_left]
        spat_tedge = self.edge_fit[reference_row, self.is_right]

        # # [DP] I am not sure how and if we need to incorporate the lines below.
        # # Mask traces that are fully masked, except if they were specifically inserted in a previous step
        # x_det_bpm = self.fully_masked_traces(flag=self.bitmask.bad_flags, exclude=self.bitmask.insert_flags)
        # spat_bedge_bpm = self.fully_masked_traces(flag=self.bitmask.bad_flags,
        #                                           exclude=self.bitmask.insert_flags)[self.is_left]
        # spat_tedge_bpm = self.fully_masked_traces(flag=self.bitmask.bad_flags,
        #                                           exclude=self.bitmask.insert_flags)[self.is_right]

        # It seems from the IDL pipeline that left and right edges from the optical model are occasionally switched
        wh = omodel_tspat != omodel_bspat
        switched = np.mean(omodel_tspat[wh] - omodel_bspat[wh]) < 0
        # Matching
        # `offsets_range` is the range of offsets in pixels allowed between the slit positions
        # predicted by the mask design and the traced slit positions.
        offsets_range = [-self.par['maskdesign_maxsep'], self.par['maskdesign_maxsep']]
        if not switched:
            # Bottom slit edge
            ind_b, dupl_b, coeff_b, sigres_b = \
                slitdesign_matching.slit_match(spat_bedge, omodel_bspat, step=self.par['maskdesign_step'],
                                               xlag_range=offsets_range, sigrej=self.par['maskdesign_sigrej'],
                                               print_matches=debug, edge='left')
            # Top slit edge
            ind_t, dupl_t, coeff_t, sigres_t = \
                slitdesign_matching.slit_match(spat_tedge, omodel_tspat, step=self.par['maskdesign_step'],
                                               xlag_range=offsets_range, sigrej=self.par['maskdesign_sigrej'],
                                               print_matches=debug, edge='right')
        else:
            # Bottom slit edge
            ind_b, dupl_b, coeff_b, sigres_b = \
                slitdesign_matching.slit_match(spat_bedge, omodel_tspat, step=self.par['maskdesign_step'],
                                               xlag_range=offsets_range, sigrej=self.par['maskdesign_sigrej'],
                                               print_matches=debug, edge='right')
            # Top slit edge
            ind_t, dupl_t, coeff_t, sigres_t = \
                slitdesign_matching.slit_match(spat_tedge, omodel_bspat, step=self.par['maskdesign_step'],
                                               xlag_range=offsets_range, sigrej=self.par['maskdesign_sigrej'],
                                               print_matches=debug, edge='left')

        if debug:
            plt.scatter(spat_bedge, omodel_bspat[ind_b], s=80, lw=2, marker='+', color='g', zorder=1,
                        label='Bottom edge: RMS={}'.format(round(sigres_b, 4)))
            plt.scatter(spat_tedge, omodel_tspat[ind_t], s=40, lw=1, marker='D', facecolors='none',
                        edgecolors='r', zorder=0, label='Top edge: RMS={}'.format(round(sigres_t, 4)))
            if np.any(dupl_b):
                plt.scatter(spat_bedge[dupl_b], spat_bedge[dupl_b], s=80, lw=2, marker='+', color='m',
                            zorder=1, label='Duplicate match (Bottom)')
            if np.any(dupl_t):
                plt.scatter(spat_tedge[dupl_t], spat_tedge[dupl_t], s=80, lw=1, marker='D', facecolors='none',
                            edgecolor='orange', zorder=1, label='Duplicate match (Top)')
            plt.plot(np.linspace(0, self.traceimg.shape[1]), np.linspace(0, self.traceimg.shape[1]), 'b-', zorder=-1)
            plt.xlabel('Edges from trace')
            plt.ylabel('Edges from model')
            plt.xlim(0, self.traceimg.shape[1] + 20)
            plt.ylim(0, self.traceimg.shape[1] + 20)
            plt.legend()
        msgs.info('SLIT_MATCH: RMS residuals for left and right edges: {}, {} pixels'.format(sigres_b, sigres_t))

        bot_edge_pred = coeff_b[0] + coeff_b[1]*omodel_bspat if not switched else coeff_b[0] + coeff_b[1]*omodel_tspat
        top_edge_pred = coeff_t[0] + coeff_t[1]*omodel_tspat if not switched else coeff_t[0] + coeff_t[1]*omodel_bspat

        # Find if there are missing traces.
        # Need exactly one occurrence of each index in "need"
        buffer = 20.
        need = ((top_edge_pred > buffer) & (bot_edge_pred < (self.traceimg.shape[1] - 1 - buffer))) & \
               ((omodel_bspat != -1) | (omodel_tspat != -1))

        # bottom edges
        needadd_b = need.copy()
        needadd_b[ind_b] = False
        needind_b = np.where(needadd_b)[0]  # edges we are missing

        # top edges
        needadd_t = need.copy()
        needadd_t[ind_t] = False
        needind_t = np.where(needadd_t)[0]  # edges we are missing

        if (needind_b.shape[0] > 0) | (needind_t.shape[0] > 0):
            msgs.warn('Missing edge traces: {} left and {} right'.format(needind_b.shape[0], needind_t.shape[0]))
            msgs.warn('Some of them may have been removed by setting the parameter `minimum_slit_length`')

        if debug:
            slitdesign_matching.plot_matches(self.edge_fit[:,self.is_left], ind_b, omodel_bspat, spat_bedge,
                                             reference_row, self.spectrograph.slitmask.slitindx, nspat=self.nspat,
                                             duplicates=dupl_b, missing=needind_b, edge='left')
            slitdesign_matching.plot_matches(self.edge_fit[:,self.is_right], ind_t, omodel_tspat, spat_tedge,
                                             reference_row, self.spectrograph.slitmask.slitindx, nspat=self.nspat,
                                             duplicates=dupl_t, missing=needind_t, edge='right')

        # [DP] The code below is to add traces that are predicted but not found. For now we leave it commented, as we
        # incorporate it in a second moment.

        # if needind_b.shape[0] > 0:
        #     ind_b = np.append(ind_b, needind_b)
        #     sortind_b = np.argsort(utils.index_of_x_eq_y(self.spectrograph.slitmask.slitindx[sortindx], ind_b, strict=True))
        #     ind_b = ind_b[sortind_b]
        #     lside = -np.ones(bot_edge_pred[needind_b].shape[0], dtype=int)
        #     missing_left_traces = self.predict_traces(bot_edge_pred[needind_b], side=lside)
        #     self.insert_traces(lside, missing_left_traces, mode='mask')
        #
        # if needind_t.shape[0] > 0:
        #     ind_t = np.append(ind_t, needind_t)
        #     sortind_t = np.argsort(utils.index_of_x_eq_y(self.spectrograph.slitmask.slitindx[sortindx], ind_t, strict=True))
        #     ind_t = ind_t[sortind_t]
        #     rside = np.ones(top_edge_pred[needind_t].shape[0], dtype=int)
        #     missing_right_traces = self.predict_traces(top_edge_pred[needind_t], side=rside)
        #     self.insert_traces(rside, missing_right_traces, mode='mask')
        #
        # if ((needind_b.shape[0] > 0)|(needind_t.shape[0] > 0))&(debug is True):
        #     slitdesign_matching.plot_matches(self.edge_fit[:, self.is_left], ind_b, omodel_bspat,
        #                  spat_bedge, reference_row, self.spectrograph.slitmask.slitindx, edge='left')
        #     slitdesign_matching.plot_matches(self.edge_fit[:, self.is_right], ind_t, omodel_tspat,
        #                  spat_tedge, reference_row, self.spectrograph.slitmask.slitindx, edge='right')

        # [DP] The index resulting from the matching are provided separately for left and right edges (ind_b, ind_t)
        # Therefore I can provide only one of the two here. If for a specific slit the right trace is found but not the
        # left one, that slit will not have an associated id.
        self.maskdef_id = self.spectrograph.slitmask.slitid[ind_b]
        # In the eventuality (hopefully remote) that some matches are duplicates, some traces will not have
        # a `maskdef_id` associated. We assign a value of `maskdef_id=-99` to those traces.
        if np.any(dupl_b):
            self.maskdef_id[dupl_b] = -99

        # [DP] The following lines create two attributes, `design` and `object`, which store
        # the matched slit-design and object information.
        # Traced edges MUST be synchronized and NO duplicate matches should exist.
        # TODO: pass the `design` and `object` to `EdgeTraceSet` and/or `SlitTraceSet` datamodel.
        # TODO: [DP] There is some inconsistency for the two scripts below if there are duplicate matches.
        #  For now I'll just put a condition, but it needs to be fixed.
        if not np.any(dupl_b):
            self._fill_design_table(ind_b, coeff_b, omodel_bspat, omodel_tspat)
            self._fill_objects_table(ind_b)

    def _fill_design_table(self, ind, coeff, omodel_bspat, omodel_tspat):
        """
        Fill :attr:`design` based on the results of the design
        registration.

        Args:
            ind (:obj:`int`):
                matched index for the slit-mask design data.
            coeff (:obj:`numpy.array`):
                Fit parameters of the cross-correlation between slit-mask design
                and traced edges.
            omodel_bspat, omodel_tspat (:obj:`float`):
                Left and right spatial position of the slit edges from optical model
        """
        # Number of slits
        nslits = len(ind)
        # Reference row
        reference_row = self.left_pca.reference_row if self.par['left_right_pca'] \
                            else self.pca.reference_row
        # Instantiate as an empty table
        self.design = EdgeTraceSet.empty_design_table(rows=nslits)
        # Save the fit parameters and the source file as table metadata
        self.design.meta['MASKFILE'] = None
        self.design.meta['MASKOFF'] = coeff[0]
        self.design.meta['MASKSCL'] = coeff[1]
        # Fill the columns
        self.design['TRACEID'] = np.arange(nslits, dtype=self.design['TRACEID'].dtype)
        self.design['TRACESROW'] = np.full(nslits, reference_row,
                                           dtype=self.design['TRACESROW'].dtype)
        self.design['TRACELPIX'] = self.edge_fit[reference_row,self.traceid<0].astype(
                                        dtype=self.design['TRACELPIX'].dtype)
        self.design['TRACERPIX'] = self.edge_fit[reference_row,self.traceid>0].astype(
                                        dtype=self.design['TRACERPIX'].dtype)
        self.design['SLITID'] = self.spectrograph.slitmask.slitid[ind].astype(
                                        dtype=self.design['SLITID'].dtype)
        self.design['SLITLOPT'] = omodel_bspat[ind].astype(dtype=self.design['SLITLOPT'].dtype)
        self.design['SLITROPT'] = omodel_tspat[ind].astype(dtype=self.design['SLITROPT'].dtype)
        if self.spectrograph.slitmask.onsky is not None:
            for i,key in enumerate(['SLITRA', 'SLITDEC', 'SLITLEN', 'SLITWID', 'SLITPA']):
                self.design[key] = self.spectrograph.slitmask.onsky[ind,i].astype(
                                        dtype=self.design[key].dtype)
        self.design['ALIGN'] = self.spectrograph.slitmask.alignment_slit[ind].astype(
                                        dtype=self.design['ALIGN'].dtype)

    def _fill_objects_table(self, ind):
        """
        Fill :attr:`objects` based on the result of the design
        registration.

        Args:
            ind (:obj:`int`):
                matched index for the slit-mask design data.
        """
        if self.spectrograph.slitmask.objects is None:
            # No object data available in slit mask design object
            self.objects = None
            return

        # The index in the objects table are found by mapping the slit
        # index of each object in the design file to the slit index
        # included in the registration
        obj_index = utils.index_of_x_eq_y(self.spectrograph.slitmask.slitindx, ind,
                                          strict=True)
        # Number of objects
        nobj = len(obj_index)
        # Instantiate an empty table
        self.objects = EdgeTraceSet.empty_objects_table(rows=nobj)
        # Fill the columns
        for i,key in enumerate(['SLITID', 'OBJID', 'OBJRA', 'OBJDEC']):
                self.objects[key] = self.spectrograph.slitmask.objects[obj_index,i].astype(
                                        dtype=self.objects[key].dtype)
        #TODO it would be good to add also the 'OBJECT' keyword from the slit-mask design.
        # 'OBJECT' is the name that the observer gives to each target and therefore more easily recognizable

        # SLITINDX is the index of the slit in the `design` table, not
        # in the original slit-mask design data
        self.objects['SLITINDX'] = utils.index_of_x_eq_y(self.objects['SLITID'],
                                                         self.design['SLITID'], strict=True)



# NOTE: I'd like us to keep this commented mask_refine function around
# for the time being.
        # def mask_refine(self, design_file=None, allow_resync=False, debug=False):
        #     """
        #     Use the mask design data to refine the edge trace positions.
        #
        #     Use of this method requires:
        #         - a PCA decomposition is available,
        #         - the traces are synchronized into left-right pairs, and
        #         - :attr:`spectrograph` has a viable `get_slitmask` method
        #           to read slit mask design data from a file. That file is
        #           either provided directly or pulled from one of the
        #           files used to construct the trace image; see
        #           `design_file`. The result of the `get_slitmask` method
        #           must provide a
        #           :class:`pypeit.spectrographs.slitmask.SlitMask` object
        #           with the slit-mask design data.
        #
        #     TODO: Traces don't need to be synchronized...
        #
        #     Also useful, but not required, is for :attr:`spectrograph` to
        #     have a viable `get_detector_map` method that provides a
        #     :class:`pypeit.spectrograph.opticalmodel.DetectorMap` object,
        #     which is used to provide a guess offset between the slit-mask
        #     focal-plane positions and the trace pixel positions. If no
        #     such `get_detector_method` exists, the guess offset is::
        #
        #         this
        #
        #     and the match between expected and traced slit positions may
        #     be unstable.
        #
        #     The method uses
        #     :class:`pypeit.spectrographs.slitmask.SlitRegister` to match
        #     the expected and traced position and identify both missing
        #     and erroneous trace locations. The former are used to add new
        #     traces and the latter are removed. The method also constructs
        #     the :attr:`design` and :attr:`objects` tables, depending on
        #     the data accessible via the
        #     :class:`pypeit.spectrographs.slitmask.SlitMask` instance.
        #
        #     Used parameters from :attr:`par`
        #     (:class:`pypeit.par.pypeitpar.EdgeTracePar`) are
        #     `left_right_pca`, `mask_reg_maxiter`, `mask_reg_maxsep`,
        #     `mask_reg_sigrej`, and `ignore_alignment`.
        #
        #     Args:
        #         design_file (:obj:`str`, optional):
        #             A file with the mask design data. If None, the method
        #             will use the first file in :attr:`files`; if
        #             :attr:`files` is also None, the method will raise an
        #             exception.
        #         debug (:obj:`bool`, optional):
        #             Run in debug mode.
        #     """
        #     # Still not done with this function...
        #     raise NotImplementedError()
        #
        #     # Check that there are traces to refine!
        #     if self.is_empty:
        #         msgs.error('No traces to refine.')
        #
        #     # The PCA decomposition must have already been determined
        #     if self.pcatype is None:
        #         msgs.error('Must first run the PCA analysis for the traces; run build_pca.')
        #
        #     # Get the file to use when parsing the mask design information
        #     _design_file = (None if self.traceimg.files is None else self.traceimg.files[0]) \
        #         if design_file is None else design_file
        #     if _design_file is None or not os.path.isfile(_design_file):
        #         msgs.error('Slit-mask design file not found or none provided.')
        #
        #     # Get the paramters to use
        #     maxiter = self.par['mask_reg_maxiter']
        #     maxsep = self.par['mask_reg_maxsep']
        #     sigma = self.par['mask_reg_sigrej']
        #     ignore_alignment = self.par['ignore_alignment']
        #
        #     # TODO: Set allow_resync and design_file to be a parameters, as
        #     # well?
        #
        #     # Read the design data
        #     msgs.info('Reading slit-mask design information from: {0}'.format(_design_file))
        #     if self.spectrograph.get_slitmask(_design_file) is None:
        #         msgs.error('Unable to read design file or no slit-mask design reader '
        #                    'defined for {0}.'.format(self.spectrograph.spectrograph))
        #
        #     # Match both left and right edges simultaneously
        #     x_design = np.array([self.spectrograph.slitmask.bottom[:, 0],
        #                          self.spectrograph.slitmask.top[:, 0]]).T.ravel()
        #     reference_row = self.left_pca.reference_row if self.par['left_right_pca'] \
        #         else self.pca.reference_row
        #     x_det = self.edge_fit[reference_row, :]
        #
        #     # Mask traces that are fully masked, except if they were
        #     # specifically inserted in a previous step
        #     # TODO: Should the BOXSLITS also be included here?
        #     x_det_bpm = self.fully_masked_traces(flag=self.bitmask.bad_flags,
        #                                          exclude=self.bitmask.insert_flags)
        #
        #     #        x_design = np.amin(self.spectrograph.slitmask.corners[:,:,0], axis=1)
        #     #        side = self.traceid < 0
        #     #        x_det = self.edge_fit[self.pca.reference_row,side]
        #
        #     #        x_design = np.amax(self.spectrograph.slitmask.corners[:,:,0], axis=1)
        #     #        side = self.traceid > 0
        #     #        x_det = self.edge_fit[self.pca.reference_row,side]
        #
        #     # Estimate the scale in pixels/mm as the telescope platescale
        #     # in arcsec/mm divided by the detector platescale in
        #     # arcsec/pixel
        #     pix_per_mm = self.spectrograph.telescope.platescale() \
        #                  / self.traceimg.detector['platescale']
        #     # / self.spectrograph.detector[self.det - 1]['platescale']
        #
        #     # If the traces are synchronized, use the estimated scale to
        #     # first mask edges that yeild slits that are too small relative
        #     # to the range of slit lengths in the mask file.
        #     if self.is_synced:
        #         slit_len_det = np.diff(x_det.reshape(-1, 2), axis=1).ravel()
        #         slit_len_mask = np.diff(x_design.reshape(-1, 2), axis=1).ravel() * pix_per_mm
        #         indx = (slit_len_det < np.amin(slit_len_mask) / 1.1) \
        #                | (slit_len_det > np.amax(slit_len_mask) * 1.1)
        #         if np.any(indx):
        #             msgs.info('Removing {0} edges that form (an) '.format(np.sum(indx) * 2)
        #                       + 'errantly small or large slit(s) compared to the mask design data.')
        #             x_det_bpm[np.repeat(indx, 2)] = True
        #
        #     # Initial guess for the offset
        #     try:
        #         raise NotImplementedError()
        #         # Try using the spectrograph detector map
        #         self.spectrograph.get_detector_map()
        #         # Set the offset based on the location of this detector
        #         offset = self.spectrograph.detector_map.image_coordinates(
        #             self.spectrograph.detector_map.npix[0] / 2,
        #             self.spectrograph.detector_map.npix[1] / 2,
        #             detector=self.traceimg.detector.det,
        #             in_mm=False)[0][0] - self.spectrograph.detector_map.npix[0] / 2
        #         # Set the bounds to some nominal fraction of the detector
        #         # size and pix/mm scale; allow for a +/- 10% deviation in
        #         # the pixel scale
        #         # TODO: Is 10% generally enough (for any instrument)? Make
        #         # this a (spectrograph-specific) parameter?
        #         offset_rng = [offset - 0.1 * self.spectrograph.detector_map.npix[0],
        #                       offset + 0.1 * self.spectrograph.detector_map.npix[0]]
        #     except:
        #         # No detector map
        #         msgs.warn('No detector map available for {0}'.format(self.spectrograph.spectrograph)
        #                   + '; attempting to match to slit-mask design anyway.')
        #         # Set the guess offset such that two sets of coordinates
        #         # are offset to their mean
        #         offset = np.mean(x_det) - np.mean(pix_per_mm * x_design)
        #         # Set the offset range
        #         offset_rng = [offset - np.absolute(np.amin(x_det) - np.amin(pix_per_mm * x_design)) * 1.1,
        #                       offset + np.absolute(np.amax(pix_per_mm * x_design) - np.amax(x_det)) * 1.1]
        #
        #     #        import pdb
        #     #        pdb.set_trace()
        #     #
        #     #        slitmask.xc_trace(x_det, x_design, pix_per_mm)
        #     #
        #     #        pdb.set_trace()
        #
        #     # The solution can be highly dependent on the initial guess for
        #     # the offset, so do an initial grid search to get close to the
        #     # solution.
        #     msgs.info('Running a grid search to try to find the best starting offset.')
        #     # Step by 2 pixels
        #     off = np.arange(offset_rng[0], offset_rng[1], 2).astype(float)
        #     rms = np.zeros_like(off, dtype=float)
        #     scl = np.zeros_like(off, dtype=float)
        #     par = np.array([0, pix_per_mm])
        #     bounds = np.array([offset_rng, [pix_per_mm / 1.1, pix_per_mm * 1.1]])
        #     register = slitmask.SlitRegister(x_det, x_design, trace_mask=x_det_bpm)
        #
        #     # NOTE: The commented approach below gets the RMS at each
        #     # offset point just using the estimated scale. This is faster
        #     # than the approach taken, but results are sensitive to the
        #     # accuracy of the estimated scale, which can lead to problems
        #     # in corner cases.
        #     #        for i in range(off.size):
        #     #            print('Grid point: {0}/{1}'.format(i+1, off.size), end='\r')
        #     #            par[0] = off[i]
        #     #            register.par = par
        #     #            minsep = register.match(unique=True)[1]
        #     #            rms[i] = sigma_clipped_stats(minsep, sigma=5)[2]
        #     #        print('Grid point: {0}/{0}'.format(off.size))
        #
        #     # For each grid point, keep the offset fixed and find the best
        #     # scale. No rejection iterations are performed.
        #     for i in range(off.size):
        #         print('Grid point: {0}/{1}'.format(i + 1, off.size), end='\r')
        #         par[0] = off[i]
        #         register.find_best_match(guess=par, fix=[True, False], bounds=bounds, penalty=False)
        #         minsep = register.match(unique=True)[1]
        #         scl[i] = register.par[1]
        #         rms[i] = sigma_clipped_stats(minsep, sigma=5)[2]
        #     print('Grid point: {0}/{0}'.format(off.size))
        #
        #     # Use the grid point with the best RMS
        #     minindx = np.argmin(rms)
        #     offset = off[minindx]
        #     best_rms = rms[minindx]
        #     msgs.info('Minimum RMS ({0:.2f}) found with offset = {1:.2f}'.format(best_rms, offset))
        #     if debug:
        #         # Plot the result
        #         ax1 = plt.subplot(211)
        #         ax1.scatter(off, rms, color='k', marker='.', s=100, lw=0, zorder=0)
        #         ax1.scatter(offset, best_rms, color='C3', marker='x', s=50, zorder=1)
        #         ax1.set_xlabel('Trace Offset (pix)')
        #         ax1.set_ylabel('RMS (det-mask; pix)')
        #         ax1.set_title('Grid search for initial offset')
        #         ax2 = plt.subplot(212, sharex=ax1)
        #         ax2.scatter(off, scl, color='k', marker='.', s=100, lw=0, zorder=0)
        #         ax2.set_ylabel('Best-fit scale')
        #         plt.show()
        #
        #     # Do the final fit with some rejection iterations
        #     register.find_best_match(guess=[offset, pix_per_mm], bounds=bounds, penalty=False,
        #                              maxiter=maxiter, maxsep=maxsep, sigma=sigma, debug=debug)
        #
        #     if debug:
        #         register.show(minmax=[0, self.nspat], synced=True)
        #
        #     # Find the missing, bad, and masked traces
        #     missing, bad = register.trace_mismatch(minmax=[0, self.nspat], synced=True)
        #     #        masked_by_registration = np.where(register.trace_mask & np.invert(x_det_bpm))[0]
        #     #        bad = np.append(bad, masked_by_registration)
        #     bad = np.append(bad, np.where(register.trace_mask | x_det_bpm)[0])
        #
        #     # Ignore missing alignment boxes
        #     if ignore_alignment:
        #         missing = missing[np.invert(self.spectrograph.slitmask.alignment_slit[missing // 2])]
        #         found_alignment_slits = register.match_index[
        #             self.spectrograph.slitmask.alignment_slit[register.match_index // 2]]
        #         bad = np.append(bad, found_alignment_slits)
        #
        #     # Report
        #     msgs.info('Best-fitting offset and scale for mask coordinates: {0:.2f} {1:.2f}'.format(
        #         *register.par))
        #     msgs.info('Traces will {0} alignment slits'.format('exclude' if ignore_alignment
        #                                                        else 'include'))
        #     msgs.info('Number of missing mask traces to insert: {0}'.format(len(missing)))
        #     msgs.info('Number of bad or alignment traces to remove: {0}'.format(len(bad)))
        #
        #     if self.is_synced and (len(missing) - len(bad)) % 2 != 0:
        #         if allow_resync:
        #             msgs.warning('Difference in added and removed traces is odd; will resync traces.')
        #         else:
        #             msgs.error('Difference in added and removed traces desyncronizes traces.')
        #
        #     if len(bad) > 0:
        #         # Remove the bad traces and rebuild the pca
        #         rmtrace = np.zeros(self.ntrace, dtype=bool)
        #         rmtrace[bad] = True
        #         self.remove_traces(rmtrace, rebuild_pca=True)
        #
        #     if len(missing) > 0:
        #         # Even indices are lefts, odd indices are rights
        #         side = missing % 2 * 2 - 1
        #         # Predict the traces using the PCA
        #         missing_traces = self.predict_traces(register.match_coo[missing], side)
        #         # Insert them
        #         self.insert_traces(side, missing_traces, mode='mask')
        #
        #     #        import pdb
        #     #        pdb.set_trace()
        #
        #     if len(bad) > 0 or len(missing) > 0:
        #         # Traces were removed and/or inserted, resync or recheck that the edges are synced.
        #         if (len(missing) - len(bad)) % 2 != 0 and allow_resync:
        #             self.sync(rebuild_pca=True)
        #         else:
        #             self.check_synced(rebuild_pca=True)
        #         reference_row = self.left_pca.reference_row if self.par['left_right_pca'] \
        #             else self.pca.reference_row
        #         # Reset the match after removing/inserting traces
        #         x_det = self.edge_fit[reference_row, :]
        #         # TODO: Should the BOXSLITS also be included here?
        #         x_det_bpm = self.fully_masked_traces(flag=self.bitmask.bad_flags,
        #                                              exclude=self.bitmask.insert_flags)
        #         register = slitmask.SlitRegister(x_det, x_design, trace_mask=x_det_bpm,
        #                                          guess=[offset, pix_per_mm], bounds=bounds,
        #                                          penalty=False, maxiter=maxiter, maxsep=maxsep,
        #                                          sigma=sigma, debug=debug, fit=True)
        #
        #         # TODO: This fit should *never* result in missing or bad
        #         # traces! Keep this for a while until we feel like we've
        #         # vetted the code well enough.
        #         missing, bad = register.trace_mismatch(minmax=[0, self.nspat], synced=True)
        #         if len(missing) != 0 or len(bad) != 0:
        #             msgs.error('CODING ERROR: Should never find missing or bad traces in re-fit!')
        #
        #     # Fill the slit-design and object tables
        #     self._fill_design_table(register, _design_file)
        #     self._fill_objects_table(register)

    def slit_spatial_center(self, normalized=True, spec=None, use_center=False, include_box=False):
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
            use_center (:obj:`bool`, optional):
                Use the measured centroids to define the slit edges
                even if the slit edges have been otherwise modeled.
            include_box (:obj:`bool`, optional):
                Include box slits in the calculated coordinates.

        Returns:
            `numpy.ma.MaskedArray`_: Spatial coordinates of the slit
            centers in pixels or in fractions of the detector. Masked
            locations are for bad/excluded slits.
        """
        if not self.is_synced:
            msgs.error('EdgeTraceSet must be synced to compute slit centers.')

        # Select the good traces
        gpm = self.good_traces(include_box=include_box)
        good_slit = np.all(gpm.reshape(-1,2), axis=1)

        # TODO: Use reference_row by default? Except that it's only
        # defined if the PCA is defined.
        _spec = self.nspec//2 if spec is None else spec

        # Synced, spatially sorted traces are always ordered in left,
        # right pairs
        trace_cen = self.edge_cen[:,gpm] if self.edge_fit is None or use_center \
                        else self.edge_fit[:,gpm]
        slit_c = np.ma.MaskedArray(np.zeros(self.nslits, dtype=float))
        slit_c[np.logical_not(good_slit)] = np.ma.masked
        slit_c[good_slit] = np.mean(trace_cen[_spec,:].reshape(-1,2), axis=1)
        return slit_c/self.nspat if normalized else slit_c

    def match_order(self):
        """
        Match synchronized slits to the expected echelle orders.

        For non-Echelle spectrographs (selected based on the
        ``pypeline`` attribute of :attr:`spectrograph`), this simply
        returns without doing anything.

        For Echelle spectrographs, this finds the best matching order
        for each left-right trace pair; the traces must have already
        been synchronized into left-right pairs. Currently, this is a
        very simple, non-optimized match:

            - The closest order from
              ``self.spectrograph.order_spat_pos`` is matched to each
              slit.
            - Any slit that is not matched to an order (only if the
              number of slits is larger than the number of expected
              orders) is flagged as ``NOORDER``.
            - If multiple slits are matched to the same order, the
              slit with the smallest on-detector separation is kept
              and the other match is ignored.
            - Any slit matched to an order with a separation above
              the provided tolerance is flagged as ORDERMISMATCH.
            - A warning is issued if the number of valid matches
              is not identical to the number of expected orders
              (``self.spectrograph.norders``). The warning includes
              the list of orders that were not identified.

        The match tolerance is et by the parameter ``order_match``.
        An offset can be applied to improve the match via the
        parameter ``order_offset``; i.e., this should minimize the
        difference between the expected order positions and
        ``self.slit_spatial_center() + self.par['order_offset']``.
        Both ``order_match`` and ``order_offset`` are given in
        fractions of the detector size along the spatial axis.

        The result of this method is to instantiate :attr:`orderid`.

        Raises:
            PypeItError:
                Raised if the number of orders or their spatial
                locations are not defined for an Echelle
                spectrograph.
        """
        if self.spectrograph.pypeline != 'Echelle':
            return

        if self.spectrograph.norders is None:
            msgs.error('Coding error: norders not defined for {0}!'.format(
                        self.spectrograph.__class__.__name__))
        if self.spectrograph.orders is None:
            msgs.error('Coding error: orders not defined for {0}!'.format(
                        self.spectrograph.__class__.__name__))
        if self.spectrograph.order_spat_pos is None:
            msgs.error('Coding error: order_spat_pos not defined for {0}!'.format(
                       self.spectrograph.__class__.__name__))

        offset = self.par['order_offset']
        if offset is None:
            offset = 0.0

        # This requires the slits to be synced! Masked elements in
        # slit_cen are for bad slits.
        slit_cen = self.slit_spatial_center()

        # Calculate the separation between the order and every
        sep = self.spectrograph.order_spat_pos[:,None] - slit_cen[None,:] - offset
        # Find the smallest offset for each order
        slit_indx = np.ma.MaskedArray(np.ma.argmin(np.absolute(sep), axis=1))

        # Minimum separation between the order and its matching slit;
        # keep the signed value for reporting, but used the absolute
        # value of the difference for vetting below.
        sep = sep[(np.arange(self.spectrograph.norders),slit_indx)]
        min_sep = np.absolute(sep)

        # Report
        msgs.info('Before vetting, the echelle order, matching left-right trace pair, and '
                  'matching separation are:')
        msgs.info(' {0:>6} {1:>4} {2:>6}'.format('ORDER', 'PAIR', 'SEP'))
        msgs.info(' {0} {1} {2}'.format('-'*6, '-'*4, '-'*6))
        for i in range(self.spectrograph.norders):
            msgs.info(' {0:>6} {1:>4} {2:6.3f}'.format(self.spectrograph.orders[i], i+1, sep[i]))
        msgs.info(' {0} {1} {2}'.format('-'*6, '-'*4, '-'*6))

        # Single slit matched to multiple orders
        uniq, cnts = np.unique(slit_indx.compressed(), return_counts=True)
        for u in uniq[cnts > 1]:
            indx = slit_indx == u
            # Keep the one with the smallest separation
            indx[np.argmin(min_sep[indx]) + np.where(indx)[0][1:]] = False
            # Disassociate the other order from any slit
            slit_indx[np.logical_not(indx) & (slit_indx == u)] = np.ma.masked

        # Flag and remove orders separated by more than the provided
        # threshold
        if self.par['order_match'] is not None:
            indx = (min_sep > self.par['order_match']) \
                        & np.logical_not(np.ma.getmaskarray(min_sep))
            if np.any(indx):
                # Flag the associated traces
                _indx = np.isin(np.absolute(self.traceid), (slit_indx[indx]).compressed()+1)
                self.edge_msk[:,_indx] = self.bitmask.turn_on(self.edge_msk[:,_indx],
                                                              'ORDERMISMATCH')
                # Disassociate these orders from any slit
                slit_indx[indx] = np.ma.masked

        # Unmatched slits
        indx = np.logical_not(np.isin(np.arange(self.nslits), slit_indx.compressed()))
        if np.any(indx):
            # This works because the traceids are sorted and synced
            indx = np.repeat(indx, 2)
            self.edge_msk[:,indx] = self.bitmask.turn_on(self.edge_msk[:,indx], 'NOORDER')

        # Warning that there are missing orders
        missed_order = np.ma.getmaskarray(slit_indx)
        if np.any(missed_order):
            msgs.warn('Did not find all orders!  Missing orders: {0}'.format(
                        ', '.join(self.spectrograph.orders[missed_order].astype(str))))

        # Instantiate the order ID; 0 means the order is unassigned
        self.orderid = np.zeros(self.nslits*2, dtype=int)
        found_orders = self.spectrograph.orders[np.logical_not(missed_order)]
        nfound = len(found_orders)
        indx = (2*slit_indx.compressed()[:,None] + np.tile(np.array([0,1]), (nfound,1))).ravel()
        self.orderid[indx] = (np.array([-1,1])[None,:]*found_orders[:,None]).ravel()

    def get_slits(self):
        """
        Use the data to instatiate the relevant
        :class:`~pypeit.slittrace.SlitTraceSet` object.

        The traced edges must have first been organized into slits;
        see :func:`sync`.

        The method automatically calls :func:`match_order` and will
        only include the "slits" that have been correctly matched to
        known orders.
        
        The :class:`~pypeit.slittrace.SlitTraceSet` object will use
        the same :attr:`master_key`, :attr:`master_dir`, and
        :attr:`reuse_masters` as this parent :class:`EdgeTraceSet`
        object.

        Returns:
            :class:`~pypeit.slittrace.SlitTraceSet`: Object holding
            the slit traces.
        """
        if not self.is_synced:
            msgs.error('Edges must be synced to construct SlitTraceSet object.')

        # For echelle spectrographs, match the left-right trace pairs
        # to echelle orders
        self.match_order()

        # Select the good traces, including only those correctly
        # matched to echelle orders for echelle spectrographs and
        # including any box slits for multi-slit observations
        gpm = self.good_traces(include_box=True, good_orders=True)
        gpm = self.synced_selection(gpm, mode='neither')

        # Initialize SlitTrace mask
        slit_bitmask = slittrace.SlitTraceBitMask()
        nslit = int(np.sum(gpm)) // 2
        slit_msk = np.zeros(nslit, dtype=slit_bitmask.minimum_dtype())
        # Loop on left edges
        for sidx, eidx in enumerate(np.where(gpm & self.is_left)[0]):
            # Loop on SlitTrace mask keys
            for key in slit_bitmask.keys():
                if key in self.bitmask.keys() and np.all(self.bitmask.flagged(
                        self.edge_msk[:, eidx], flag=key)):
                    slit_msk[sidx] = slit_bitmask.turn_on(slit_msk[sidx], key)

        # Parse the data for the SlitTraceSet
        left = self.edge_fit[:,gpm & self.is_left]
        right = self.edge_fit[:,gpm & self.is_right]
        binspec, binspat = parse.parse_binning(self.traceimg.detector.binning)
        ech_order = None if self.orderid is None else self.orderid[gpm][1::2]
        if self.spectrograph.spec_min_max is None or ech_order is None:
            specmin = np.asarray([-np.inf]*nslit)
            specmax = np.asarray([np.inf]*nslit)
        else:
            specmin, specmax = self.spectrograph.spec_min_max
            indx = utils.index_of_x_eq_y(self.spectrograph.orders, ech_order)
            specmin = specmin[indx]/binspec
            specmax = specmax[indx]/binspec

        # Instantiate and return
        return slittrace.SlitTraceSet(left, right, self.spectrograph.pypeline, nspat=self.nspat,
                                      PYP_SPEC=self.spectrograph.spectrograph, specmin=specmin,
                                      specmax=specmax, binspec=binspec, binspat=binspat,
                                      pad=self.par['pad'], mask_init=slit_msk,
                                      maskdef_id=self.maskdef_id,
                                      ech_order=ech_order)


