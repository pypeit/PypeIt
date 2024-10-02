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
detector is as follows:

.. code-block:: python

    # Imports
    from pypeit.pypeit import PypeIt
    from pypeit import edgetrace
    from pypeit.images import buildimage

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
    calib_dir = rdx.par['calibrations']['caldir']
    setup = rdx.fitstbl['setup'][tbl_rows[0]]
    calib_id = rdx.fitstbl['calib'][tbl_rows[0]]

    # Skip the bias subtraction, if reasonable; see
    # pypeit.biasframe.BiasFrame to construct a bias to subtract from
    # the TraceImage
    rdx.par['calibrations']['traceframe']['process']['bias'] = 'skip'

    # Construct the TraceImage
    traceImage = buildimage.buildimage_fromlist(rdx.spectrograph, det,
                                                rdx.par['calibrations']['traceframe'],
                                                files, calib_dir=self.calib_dir,
                                                setup=setup, calib_id=calib_id)

    # Then run the edge tracing.  This performs the automatic tracing.
    edges = edgetrace.EdgeTraceSet(traceImage, rdx.spectrograph,
                                   rdx.par['calibrations']['slitedges'], auto=True)
    # You can look at the results using the show method:
    edges.show()
    # Or in ginga viewer
    edges.show(in_ginga=True)
    # And you can save the results to a file
    edges.to_file()

If you want to instead start without a pypeit file, you could do the
following for, e.g., a single unbinned Keck DEIMOS flat-field
exposure in a fits file called `trace_file`:

.. code-block:: python

    import os
    from pypeit import edgetrace
    from pypeit.images import buildimage
    from pypeit.spectrographs.util import load_spectrograph

    spec = load_spectrograph('keck_deimos')
    par = spec.default_pypeit_par()
    par['calibrations']['traceframe']['process']['bias'] = 'skip'
    # Make any desired changes to the parameters here
    det = 3
    calib_dir = par['calibrations']['caldir']

    # Construct the TraceImage
    traceImage = buildimage.buildimage_fromlist(spec, det, par['calibrations']['traceframe'],
                                                [trace_file], calib_dir=self.calib_dir,
                                                setup='A', calib_id=1)

    edges = edgetrace.EdgeTraceSet(traceImage, spec, par['calibrations']['slitedges'], auto=True)
    edges.to_file()

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import os
import inspect
from pathlib import Path
from collections import OrderedDict

from IPython import embed

import numpy as np

from scipy import ndimage

import matplotlib
from matplotlib import pyplot as plt
from matplotlib import ticker, rc

from astropy import table

from pypeit import msgs
from pypeit import utils
from pypeit import sampling
from pypeit import slittrace
from pypeit.datamodel import DataContainer
from pypeit import calibframe
from pypeit.bitmask import BitMask
from pypeit.display import display
from pypeit.par.pypeitpar import EdgeTracePar
from pypeit.core import parse, procimg, trace, slitdesign_matching
from pypeit.core import fitting
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
           ('ABNORMALSLIT_SHORT', 'Slit formed by left and right edge is abnormally short'),
            ('ABNORMALSLIT_LONG', 'Slit formed by left and right edge is abnormally long'),
                   ('USERRMSLIT', 'Slit removed by user'),
                      ('NOORDER', '(DEPRECATED as of version 1.15.0; use ORDERMISMATCH).  '
                                  'Unable to associate this trace with an echelle order (echelle '
                                  'spectrographs only)'),
                ('ORDERMISMATCH', 'Slit traces are not well matched to any echelle order (echelle '
                                  'spectrographs only)'),
                  ('ORDERINSERT', 'Trace was inserted as the expected location of an echelle '
                                  'order missed by the automated tracing'),
            ('LARGELENGTHCHANGE', 'Large difference in the slit length as a function of '
                                  'wavelength.')])
        super().__init__(list(mask.keys()), descr=list(mask.values()))

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
        return ['USERINSERT', 'SYNCINSERT', 'MASKINSERT', 'ORPHANINSERT', 'ORDERINSERT']

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


class EdgeTraceSet(calibframe.CalibFrame):
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

    To load an existing calibration file with the result of a trace, use the
    :attr:`from_file` method::

        edges = EdgeTraceSet.from_file(file)

    Most commonly, one will use the automatic tracing routine to
    trace the slit edges; see the description of the steps used
    during auto-tracing in the docs for :func:`auto_trace`.

    The success of the tracing critically depends on the parameters
    used. The defaults are tuned for each spectrograph based on
    testing using data in the PypeIt development suite. See
    :ref:`pypeitpar` for the full documentation of the
    :class:`~pypeit.par.pypeitpar.EdgeTracePar` parameters.

    Finally, note that the :attr:`design` and :attr:`object` data are
    currently empty, as part of a development path for matching slits
    traced on the detector to slits expected from provided metadata.
    Once finished these objects will only contain data for
    spectrograph output files that provide the relevant metadata.

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_edgetraceset.rst
   
    .. todo:
        - Include a method/check that determines if traces cross one
        another anywhere.
        - Provide some guidance for the parameters to use.

    Args:
        traceimg (:class:`~pypeit.images.buildimage.TraceImage`):
            Two-dimensional image used to trace slit edges.  The object provides
            the image, the bad-pixel mask, the detector information, and (one
            of) the original raw file name when matching slits to a design file.
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The object that sets the instrument used to take the
            observations. Used to set :attr:`spectrograph`.
        par (:class:`~pypeit.par.pypeitpar.EdgeTracePar`):
            The parameters used to guide slit tracing. Used to set
            :attr:`par`.
        qa_path (:obj:`str`, `Path`_, optional):
            Directory for QA output. If None, no QA plots are
            provided.
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
            (:class:`~pypeit.images.buildimage.TraceImage`):
            See argument list.
        spectrograph
            (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            See argument list.
        par (:class:`~pypeit.par.pypeitpar.EdgeTracePar`):
            See argument list.
        files (:obj:`list`):
            The list of raw files used to construct the trace image
            (:attr:`img`). Only defined if argument `img` in
            :func:`initial_trace` or :func:`auto_trace` is a
            :class:`~pypeit.images.buildimage.TraceImage` object.
        img (`numpy.ndarray`_):
            Convenience for now.
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
        pca (:obj:`list`, :class:`~pypeit.tracepca.TracePCA`):
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
        objects (`astropy.table.Table`_):
            Collated object ID and coordinate information matched to
            the design table.
        qa_path (`Path`_):
            Directory for QA output. If None, no QA plots are
            provided.
        log (:obj:`list`):
            A list of strings indicating the main methods applied
            when tracing.
        maskdef_id (`numpy.ndarray`_):
            An integer array with the slitmask IDs assigned to each
            trace.
        omodel_bspat (`numpy.ndarray`_):
            A floating-point array with the location of the slit LEFT edge,
            averaged along the spectral direction, predicted by the optical model
            (before x-correlation with traced edges)
        omodel_tspat (`numpy.ndarray`_):
            A floating-point array with the location of the slit RIGHT edge,
            averaged along the spectral direction, predicted by the optical model
            (before x-correlation with traced edges)
        coeff_b (`numpy.ndarray`_):
            A floating-point array with the coefficients (offset, scale) of the
            x-correlation between LEFT edges predicted by the optical model and the
            ones traced on the image.
        coeff_t (`numpy.ndarray`_):
            A floating-point array with the coefficients (offset, scale) of the
            x-correlation between RIGHT edges predicted by the optical model and the
            ones traced on the image.
        maskfile (:obj:`str`):
            Full path to the file used to extract slit-mask information
    """

    calib_type = 'Edges'
    """Root name for the calibration frame file."""

    calib_file_format = 'fits.gz'
    """Calibration frame file format."""

    bitmask = EdgeTraceBitMask()
    """Bit interpreter for the edge tracing mask."""

    version = '1.0.1'
    """DataContainer datamodel version."""

    # TODO: Add this to the DataContainer base class?
    output_float_dtype = np.float32
    """Regardless of datamodel, output floating-point data have this fixed bit size."""

    internals = calibframe.CalibFrame.internals \
                 + ['spectrograph',     # Spectrograph instance
                    'par',              # EdgeTracePar instance
                    'qa_path',          # Path for the QA plots
                    'edge_img',         # Array with the spatial pixel nearest to each trace edge.
                    'sobelsig_left',    # Sobel filtered image used to trace left edges
                    'sobelsig_right',   # Sobel filtered image used to trace right edges
                    'design',           # Table that collates slit-mask design data matched to
                                        #   the edge traces
                    'objects',          # Table that collates object information, if available
                                        #   in the slit-mask design, matched to the `design` table.
                    'log',              # Log of methods applied
                    'omodel_bspat',     # Left edges predicted by the optical model
                                        #   (before x-correlation)
                    'omodel_tspat',     # Right edges predicted by the optical model
                                        #   (before x-correlation)
                    'cc_params_b',      # Parameters of the x-correlation between LEFT edges
                                        #   predicted by the slitmask design and the one traced
                                        #   on the image.
                    'cc_params_t',      # Parameters of the x-correlation between RIGHT edges
                                        #   predicted by the slitmask design and the one traced
                                        #   on the image.
                    'maskfile',         # File used to slurp in slit-mask design
                    'slitmask',         # SlitMask instance that hold info on slitmask design
                    'success']          # Flag that the automatic edge tracing was successful
    """
    Attributes kept separate from the datamodel.
    """

    datamodel = {'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 'dispname': dict(otype=str, descr='Spectrograph disperser name.'),
                 'traceimg': dict(otype=TraceImage,
                                   descr='Image used to construct the edge traces; see '
                                         ':class:`~pypeit.images.buildimage.TraceImage` and '
                                         ':class:`~pypeit.images.pypeitimage.PypeItImage`.'),
                 'nspec': dict(otype=int, descr='Image pixels in the spectral direction.'),
                 'nspat': dict(otype=int, descr='Image pixels in the spatial direction.'),
                 'tracebpm': dict(otype=np.ndarray, atype=(bool, np.bool_),
                                  descr='Bad-pixel mask for trace image'),
                 'sobelsig': dict(otype=np.ndarray, atype=(float, np.floating),
                                  descr='Sobel-filtered image used to detect edges'),
                 'traceid': dict(otype=np.ndarray, atype=(int, np.integer),
                                 descr='ID number for the edge traces.  Negative and positive '
                                       'IDs are for, respectively, left and right edges.'),
                 'maskdef_id': dict(otype=np.ndarray, atype=(int, np.integer),
                                 descr='slitmask ID number for the edge traces. '
                                       'IDs are for, respectively, left and right edges.  Only '
                                       'defined if mask-design metadata is available.'),
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
                             descr='The PCA decomposition of all edge traces.  Not defined if PCA '
                                   'separated between left and right traces (i.e., the '
                                   'left_right_pca parameter is True).  See '
                                   ':class:`~pypeit.tracepca.TracePCA`.'),
                 'left_pca': dict(otype=TracePCA,
                                  descr='The PCA decomposition of the left-edge traces.  Not '
                                        'defined if PCA performed on all traces, regardless of '
                                        'edge side (i.e., the left_right_pca parameter is '
                                        'False).  See :class:`~pypeit.tracepca.TracePCA`.'),
                 'right_pca': dict(otype=TracePCA,
                                   descr='The PCA decomposition of the right-edge traces.  Not '
                                         'defined if PCA performed on all traces, regardless of '
                                         'edge side (i.e., the left_right_pca parameter is '
                                         'False).  See :class:`~pypeit.tracepca.TracePCA`.'),
                 'pcatype': dict(otype=str,
                                 descr='String identifier for the measurements used to construct '
                                       'the PCA (center or fit)')}
    """DataContainer datamodel."""

    def __init__(self, traceimg, spectrograph, par, qa_path=None, auto=False, debug=False,
                 show_stages=False):

        # Instantiate as an empty DataContainer
        super().__init__()

        # Check input types
        if not isinstance(traceimg, TraceImage):
            msgs.error('Input traceimg must be a TraceImage object.')
        if not isinstance(spectrograph, Spectrograph):
            msgs.error('Input spectrograph must be a Spectrograph object.')
        if not isinstance(par, EdgeTracePar):
            msgs.error('Input par must be an EdgeTracePar object.')

        self.traceimg = traceimg                        # Input TraceImage
        self.nspec, self.nspat = self.traceimg.shape    # The shape of the trace image
        self.spectrograph = spectrograph                # Spectrograph used to take the data
        self.PYP_SPEC = spectrograph.name               # For the Header.  Will be in datamodel
        self.dispname = spectrograph.dispname           # Spectrograph disperser
        self.par = par                                  # Parameters used for slit edge tracing
        self.maskdef_id = None                          # Slit ID number from slit-mask design
                                                        # matched to traced slits
        # Directory for QA plots
        self.qa_path = None if qa_path is None else Path(qa_path).absolute()

        # Inherit the calibration frame attributes from the trace image:
        self.copy_calib_internals(self.traceimg)

        # NOTE: This means that, no matter what, every instance of
        # EdgeTraceSet should have a sobelsig attribute that is *not*
        # None.
        self.success = False
        if auto:
            # Run the automatic tracing
            self.auto_trace(debug=debug, show_stages=show_stages)
        else:
            # Only get the initial trace
            self.initial_trace()

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
        self.omodel_bspat = None
        self.omodel_tspat = None
        self.cc_params_b = None
        self.cc_params_t = None
        self.maskfile = None
        self.slitmask = None
        self.success = False

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
                    table.Column(name='SPAT_ID', dtype=int, length=length,
                                 description='ID Number assigned by the pypeline to each slit'),
                    table.Column(name='MASKDEF_ID', dtype=int, length=length,
                                 description='Slit ID Number from slit-mask design'),
                    # table.Column(name='SLITLFOC', dtype=float, length=length,
                    #              description='Left edge of the slit in mm at the focal plane'),
                    # table.Column(name='SLITRFOC', dtype=float, length=length,
                    #              description='Right edge of the slit in mm at the focal plane'),
                    table.Column(name='SLITLMASKDEF', dtype=float, length=length,
                                description='Left edge of the slit in pixel from slit-mask design before x-correlation'),
                    table.Column(name='SLITRMASKDEF', dtype=float, length=length,
                                description='Right edge of the slit in pixel from slit-mask design before x-correlation'),
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
        # TODO -- Rename TOPDIST -> LEFTDIST and BOTDIST -> RIGHTDIST
        return table.Table([
                    table.Column(name='OBJID', dtype=int, length=length,
                                 description='Object ID Number'),
                    table.Column(name='OBJRA', dtype=float, length=length,
                                 description='Right ascension of the object (deg)'),
                    table.Column(name='OBJDEC', dtype=float, length=length,
                                 description='Declination of the object (deg)'),
                    table.Column(name='OBJNAME', dtype='<U32', length=length,
                                 description='Object name assigned by the observer'),
                    table.Column(name='OBJMAG', dtype=float, length=length,
                                 description='Object magnitude provided by the observer'),
                    table.Column(name='OBJMAG_BAND', dtype='<U32', length=length,
                                 description='Band of the magnitude provided by the observer'),
                    table.Column(name='MASKDEF_ID', dtype=int, length=length,
                                 description='Slit ID Number from slit-mask design'),
                    table.Column(name='OBJ_TOPDIST', dtype=float, length=length,
                                 description='Projected distance (in arcsec) of the object from the '
                                             'left edge of the slit (in PypeIt orientation).'),
                    table.Column(name='OBJ_BOTDIST', dtype=float, length=length,
                                 description='Projected distance (in arcsec) of the object from the '
                                             'right edge of the slit (in PypeIt orientation)'),
                    table.Column(name='TRACEID', dtype=int, length=length,
                                 description='Row index that matches TRACEID in the design table')
                           ])

    def rectify(self, flux, bpm=None, extract_width=None, mask_threshold=0.5, side='left'):
        r""""
        Rectify the provided image based on the current edge trace
        PCA model.

        The is primarily a wrapper for
        :func:`pypeit.sampling.rectify_image`; see its documentation
        for more detail.

        Used parameters from :attr:`par`
        (:class:`~pypeit.par.pypeitpar.EdgeTracePar`) are
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
        _pca = (self.left_pca if side == 'left' else self.right_pca) \
                    if self.par['left_right_pca'] else self.pca

        # Get the traces that cross the reference spatial position at
        # the first and last pixels of the image
        first_last_trace = _pca.predict(np.array([0,self.nspat-1]))
        # Use these two traces to define the spatial pixel coordinates
        # to sample
        start = np.ceil(np.amax(np.amin(first_last_trace, axis=1))).astype(int)
        buffer = self.nspat - np.floor(np.amin(np.amax(first_last_trace, axis=1))).astype(int) \
                    + start
        # Rectify the image
        # TODO: This has its limitations if the PCA is highly non-linear.
        ocol = np.arange(self.nspat+buffer)-start
        return sampling.rectify_image(flux, _pca.predict(ocol), bpm=bpm, ocol=ocol,
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
        (:class:`~pypeit.par.pypeitpar.EdgeTracePar`) are
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

        # Match the traces found in the image with the ones predicted by
        # the slit-mask design. If not expected traces are found in the image, they
        # will be removed. If traces are missed, they will be added.
        if not self.is_empty and self.par['use_maskdesign']:
            msgs.info('-' * 50)
            msgs.info('{0:^50}'.format('Matching traces to the slit-mask design'))
            msgs.info('-' * 50)
            self.maskdesign_matching(debug=debug)
            if show_stages:
                self.show(title='After matching to slit-mask design metadata.')
            if np.all(self.bitmask.flagged(self.edge_msk, self.bitmask.bad_flags)):
                msgs.error('All traces masked!  Problem with mask-design matching, which may be '
                           'due to spurious edges.  Try changing the edge detection threshold '
                           '(edge_thresh) and troubleshooting the problem using the '
                           'pypeit_trace_edges script.')

        if self.par['auto_pca'] and not self.can_pca() and not self.is_empty and self.par['sync_predict'] == 'pca':
            # TODO: This causes the code to fault. Maybe there's a way
            # to catch this earlier on?
            msgs.warn('Sync predict cannot use PCA because too few edges were found.  If you are '
                       'reducing multislit or echelle data, you may need a better trace image or '
                       'change the mode used to predict traces (see below).  If you are reducing '
                       'longslit data, make sure to set the sync_predict parameter to nearest: '
                       + msgs.newline() +
                       '    [calibrations]' + msgs.newline() +
                       '        [[slitedges]]' + msgs.newline() +
                       '            sync_predict = nearest')
        #            self.par['sync_predict'] = 'nearest'
            self.success = False
        else:
            # Left-right synchronize the traces
            # TODO: If the object "is_empty" at this point, sync adds two edges at
            # the left and right edges of the detectors. This means that detectors
            # with no slits (e.g., an underfilled mask in DEIMOS) will be treated
            # like a long-slit observation. At best, that will lead to a lot of
            # wasted time in the reductions; at worst, it will just cause the code
            # to fault later on.
            self.success = self.sync()
            if not self.success:
                return
            if show_stages:
                self.show(title='After synchronizing left-right traces into slits')

        if not self.is_empty and self.par['add_missed_orders']:
            # Refine the order traces
            self.order_refine(debug=debug)
            # Check that the edges are still synced
            if not self.is_synced:
                msgs.error('Traces are no longer synced after adding in missed orders.')

#           KBW: Keep this code around for a while.  It is the old code that
#           resynced the edges just after adding in new orders.  Nominally, this
#           shouldn't be necessary, but the comment suggests this may be
#           necessary if orders are missed.  We should keep this around until
#           we're sure it's not needed.
#            # Check that the edges are still sinked (not overkill if orders are
#            # missed)
#            self.success = self.sync(debug=True)
#            if not self.success:
#                return

            if show_stages:
                self.show(title='After adding in missing orders')

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
                self.add_user_traces(add_user_slits, method=self.par['add_predict'])
        if show_stages and ad_rm:
            self.show(title='After user-dictated adding/removing slits')

        # TODO: If slits are added/removed, should the code check again if the
        # traces are synced?

        # TODO: Add a parameter and an if statement that will allow for
        # this.
        # `peak_refine` ends with the traces being described by a
        # polynomial. Instead finish by reconstructing the trace models
        # using the PCA decomposition
#        self.pca_refine(debug=debug)
#        if show_stages:
#            self.show()

        # Add this to the log
        # TODO: Is this log ever used? We should probably get rid of it...
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
        # TODO: Construct the default bpm using *any* pixel that is masked, not
        # just the ones flagged as BPM?
        self.tracebpm = self.traceimg.select_flag(flag='BPM') if bpm is None else bpm.astype(bool)
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

        # Only after the correction for the bad column, update the bmp with the
        # regions to exclude from slit tracing
        if self.par['exclude_regions'] is not None:
            reg_start, reg_end = self._parse_exclude_regions()
            for i in range(reg_start.size):
                self.tracebpm[:, reg_start[i]:reg_end[i]] = True

        self.sobelsig, edge_img \
                = trace.detect_slit_edges(_img, bpm=self.tracebpm,
                                          median_iterations=self.par['filt_iter'],
                                          sobel_mode=self.par['sobel_mode'],
                                          sigdetect=self.par['edge_thresh'],
                                          sobel_enhance=self.par['sobel_enhance'])
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
            row, col = np.where(np.logical_not(trace_id_img.mask)
                                    & (trace_id_img.data == self.traceid[i]))
            self.edge_img[row,i] = col
            self.edge_msk[row,i] = 0            # Turn-off the mask

            # Flag any insert traces
            row, col = np.where(np.logical_not(trace_id_img.mask) & inserted_edge
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

    def _parse_exclude_regions(self):
        """
        Parse the `exclude_regions` parset.

        Returns:
            :obj:`tuple`: Returns two arrays with the starting
            and ending pixels of the regions to exclude in this detector
        """
        if self.par['exclude_regions'] is None:
            msgs.error('No regions to exclude have been provided. '
                       'To do so, see parameter `exclude_regions` in `EdgeTracePar`')

        # create the arrays with det, starting pixels and ending pixels
        dets = np.zeros(len(self.par['exclude_regions']), dtype=int)
        reg_start = np.zeros(len(self.par['exclude_regions']), dtype=int)
        reg_end = np.zeros(len(self.par['exclude_regions']), dtype=int)
        # fill them
        for i, region in enumerate(self.par['exclude_regions']):
            dets[i], reg_start[i], reg_end[i] = [int(s) for s in region.split(':')]
        # select the current detector
        this_det = dets == self.spectrograph.ndet

        return reg_start[this_det], reg_end[this_det]

    def _base_header(self, hdr=None):
        """
        Construct the baseline header for all HDU extensions.

        This appends the :class:`~pypeit.par.pypeitpar.EdgeTracePar` and
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
        _hdr = super()._base_header(hdr=hdr)
        _hdr['QAPATH'] = 'None' if self.qa_path is None else str(self.qa_path)
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
            return super().to_hdu(**kwargs)

        # TODO: We need a better solution for multiple levels of nested
        # DataContainers. Here the commpication is that we're writing
        # many DataContainers of a single derived class (TracePCA) to
        # the same file.

        # Temporarily erase the left and right pca so
        # they're not written
        _left_pca, _right_pca = self.left_pca, self.right_pca
        self.left_pca, self.right_pca = None, None

        # Run the default (with add_primary = False)
        hdu = super().to_hdu(**kwargs)

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
        d, version_passed, type_passed, parsed_hdus = cls._parse(hdu)
        # Check
        cls._check_parsed(version_passed, type_passed, chk_version=chk_version)

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
        self = super().from_dict(d=d)

        # Calibration frame attributes
        # NOTE: If multiple HDUs are parsed, this assumes that the information
        # necessary to set all the calib internals is always in *every* header.
        # BEWARE!
        self.calib_keys_from_header(hdu[parsed_hdus[0]].header)

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
        self.qa_path = Path(hdu['SOBELSIG'].header['QAPATH']).absolute()

        # Check the bitmasks
        hdr_bitmask = BitMask.from_header(hdu['SOBELSIG'].header)
        if chk_version and hdr_bitmask.bits != self.bitmask.bits:
            msgs.error('The bitmask in this fits file appear to be out of date!  Recreate this '
                       'file by re-running the relevant script or set chk_version=False.',
                       cls='PypeItBitMaskError')

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
        if slits is None and not self.is_empty:
            # Use the internals. Any masked data is excluded; masked
            # data to be plotted are held in a separate array. This
            # means that errors and fits are currently never plotted
            # for masked data.
            _flag = None if flag in ['any', None] else np.atleast_1d(flag)
            cen = np.ma.MaskedArray(self.edge_cen, mask=self.bitmask.flagged(self.edge_msk))
            fit = self.edge_fit
            err = np.ma.MaskedArray(self.edge_err, mask=np.ma.getmaskarray(cen).copy())
            msk = None if flag is None \
                    else np.ma.MaskedArray(self.edge_cen, mask=np.logical_not(
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
            maskdef_ids = None
            if synced:
                _trc = cen if fit is None else fit
                half = _trc.shape[0] // 2
                slit_ids = ((_trc[half, gpm & is_left]
                             + _trc[half, gpm & is_right]) / 2.).astype(int)
                if self.maskdef_id is not None:
                    maskdef_ids = self.maskdef_id[gpm & self.is_left]
                    maskdef_ids[maskdef_ids == -99] = self.maskdef_id[gpm & self.is_right][maskdef_ids == -99]
        elif slits is not None:
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
            is_right = np.logical_not(is_left)
            gpm = np.ones(2*nslits, dtype=bool)
            traceid = np.concatenate((-np.arange(nslits), np.arange(nslits)))
            synced = True
            # JFH Trying to make this work for the case where slits are
            # known, but this functionality appears to be never used? I
            # think slits needs its own show method.
            slit_ids = slits.slitord_id
            maskdef_ids = slits.maskdef_id

        # TODO: The above should set everything (self shouldn't be used
        # below). We need a single show method that both EdgeTraceSet
        # and SlitTraceSet can use.

        if in_ginga:
            # Set up the appropriate keyword arguments for the IDs
            id_kwargs = {'slit_ids': slit_ids} if synced \
                            else {'left_ids': traceid[gpm & is_left],
                                  'right_ids': traceid[gpm & is_right]}
            id_kwargs['maskdef_ids'] = maskdef_ids
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
                    ofile = self.qa_path / 'PNGs', f'{fileroot}_{str(page).zfill(ndig)}.png'
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
#        _bpm = np.zeros_like(self.sobelsig, dtype=bool) \
        _bpm = self.traceimg.select_flag(flag='BPM') \
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
                bad_trace_pixels = self.bitmask.flagged(self.edge_msk, flag=self.bitmask.bad_flags)
                untraced = indx.copy()
                while np.any(untraced):
                    # Get the starting row
                    _start_indx = trace.most_common_trace_row(bad_trace_pixels[:,untraced],
                                                              valid_frac=1.) \
                                        if start_indx is None else start_indx
                    # Select the edges to follow
                    to_trace = untraced & np.logical_not(bad_trace_pixels[_start_indx,:])
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
        (:class:`~pypeit.par.pypeitpar.EdgeTracePar`) is ``match_tol``.

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
            compare = (side*self.traceid > 0) & np.logical_not(indx)
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
                repeat[indx] = (mindiff.data < self.par['match_tol']) & np.logical_not(mindiff.mask)
                if np.any(repeat):
                    _msk[:,repeat] = self.bitmask.turn_on(_msk[:,repeat], 'DUPLICATE')
                    msgs.info('Found {0} repeat trace(s).'.format(np.sum(repeat)))

        # Find spectrally short traces
        short = np.zeros_like(indx, dtype=bool)
        if minimum_spec_length is not None:
            msgs.info('Minimum spectral length of any trace (pixels): {0:.2f}'.format(
                      minimum_spec_length))
            short[indx] = np.sum(np.logical_not(_bpm[:,indx]), axis=0) < minimum_spec_length
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
                                & np.logical_not(_bpm[self.nspec//2,indx])
            if np.any(hit_min):
                _msk[:,hit_min] = self.bitmask.turn_on(_msk[:,hit_min], 'HITMIN')
                msgs.info('{0} trace(s) hit the minimum centroid value.'.format(np.sum(hit_min)))
         
        # Find traces that are at the maximum column at the center
        # spectral row
        # TODO: Why only the center row?
        hit_max = np.zeros_like(indx, dtype=bool)
        if max_spatial is not None:
            hit_max[indx] = (col[self.nspec//2,indx] >= max_spatial) \
                                & np.logical_not(_bpm[self.nspec//2,indx])
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
        good = indx & np.logical_not(bad)
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
                                & np.logical_not(np.all(cen_match.mask, axis=0)))[0]
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
                gpm = np.logical_not(np.ma.getmaskarray(merged_trace))
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

    def _flag_edges(self, trace_cen, indx, flg):
        """
        Convenience function for slit checking, which performs the same
        operations to the set of edges flagged for the many reasons iterated
        through in :func:`check_synced`.

        This modifies :attr:`edge_msk` directly.

        Args:
            trace_cen (`numpy.ndarray`_):
                The 2D array with the edge traces.
            indx (`numpy.ndarray`_):
                The boolean array selecting the edges to be flagged.
            flg (:obj:`str`):
                The bit flag to be assigned.
        """
        if np.sum(indx) == self.ntrace:
            if self.ntrace == 2:
                # TODO: I *really* don't like this because it has
                # the potential to yield an infinite loop, but it's
                # also the simplest approach.
                return self._masked_single_slit(trace_cen)
            msgs.warn('All slits have been flagged!')
        if np.any(indx):
            msgs.info(f'Flagging {np.sum(indx)//2} slits as {flg}!')
            self.edge_msk[:,indx] = self.bitmask.turn_on(self.edge_msk[:,indx], flg)

    def check_synced(self, rebuild_pca=False):
        """
        Quality check and masking of the synchronized edges.

        Before executing this method, the slit edges must be synchronized (see
        :func:`sync`) and ordered spatially in left-right pairs (see
        :func:`spatial_sort`); only the former is checked explicitly. Any traces
        fully masked as bad (see :func:`clean_traces`) are removed, along with
        its synchronized partner.

        Used parameters from :attr:`par`
        (:class:`~pypeit.par.pypeitpar.EdgeTracePar`) are ``minimum_slit_gap``,
        ``minimum_slit_length``, ``minimum_slit_length_sci``, and
        ``length_range``.

        Checks are:

            - Any trace falling off the edge of the detector is masked (see
              :class:`EdgeTraceBitMask`). This is the only check performed by
              default (i.e., when no keyword arguments are provided).
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

        Returns:
            :obj:`bool`: Flag that checking and cleaning traces maintained
            left-right syncronization.  If False, traces will need to be
            re-syncronized.
        """
        if self.is_empty:
            msgs.warn('No traces to check.')
            return

        # Decide if the PCA should be rebuilt
        _rebuild_pca = rebuild_pca and self.pcatype is not None and self.can_pca()
        if rebuild_pca and not _rebuild_pca:
            msgs.warn('Rebuilding the PCA was requested but is not possible.')

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
        dlength_rtol = self.par['dlength_range']
        dlength_atol = None
        length_atol = None
        length_atol_sci = None
        length_rtol = self.par['length_range']
        if self.par['minimum_slit_dlength'] is not None \
                or self.par['minimum_slit_length'] is not None \
                or self.par['minimum_slit_length_sci'] is not None \
                or self.par['minimum_slit_gap'] is not None:
            platescale = parse.parse_binning(self.traceimg.detector.binning)[1] \
                            * self.traceimg.detector['platescale']
            msgs.info('Binning: {0}'.format(self.traceimg.detector.binning))
            msgs.info('Platescale per binned pixel: {0}'.format(platescale))
            if self.par['minimum_slit_dlength'] is not None:
                dlength_atol = self.par['minimum_slit_dlength']/platescale
            if self.par['minimum_slit_length'] is not None:
                length_atol = self.par['minimum_slit_length']/platescale
            if self.par['minimum_slit_length_sci'] is not None:
                length_atol_sci = self.par['minimum_slit_length_sci']/platescale
            if self.par['minimum_slit_gap'] is not None:
                gap_atol = self.par['minimum_slit_gap']/platescale

        msgs.info('Minimum slit gap (binned pixels): {0}'.format(gap_atol))
        msgs.info('Minimum change in slit length (binned pixels): {0}'.format(dlength_atol))
        msgs.info('Range in the change in slit length not limited' if dlength_rtol is None else
                  f'Range in the change in slit length limited to +/-{dlength_rtol*100:.1f}%')
        msgs.info('Minimum slit length (binned pixels): {0}'.format(length_atol))
        msgs.info('Minimum science slit length (binned pixels): {0}'.format(length_atol_sci))
        msgs.info('Range in slit length not limited' if length_rtol is None else
                  f'Range in slit length limited to +/-{length_rtol*100:.1f}%')

        if length_rtol is None and self.par['overlap']:
            msgs.warn('Overlap keyword ignored!  Must set length_range to identify abnormally '
                      'short slits.')

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
                self.remove_traces(rmtrace, rebuild_pca=_rebuild_pca)
                # TODO: This should never happen, but keep this around
                # until we're sure it doesn't.
                if self.is_empty:
                    msgs.error('Coding error: Removing gaps removed all traces.')
                # Reset the trace center data to use
                trace_cen = self.edge_cen if self.edge_fit is None else self.edge_fit

        # Calculate the slit length and gap
        slit_length = np.squeeze(np.diff(trace_cen.reshape(self.nspec,-1,2), axis=-1))
        med_slit_length = np.median(slit_length, axis=0)

        if dlength_atol is not None:
            # Check the length of each slit against the median value, flagging
            # any above the provided *absolute* tolerance (both edges are
            # flagged)
            indx = np.repeat(np.any(np.absolute(slit_length-med_slit_length[None,:])
                                    > dlength_atol, axis=0), 2)
            self._flag_edges(trace_cen, indx, 'LARGELENGTHCHANGE')

        if dlength_rtol is not None:
            # For slits with larger than 0 lengths (i.e., valid slits), check
            # their length against the median value, flagging any exhibit a
            # *relative* change above the provided tolerance (both edges are
            # flagged)
            indx = np.all(slit_length > 0, axis=0)
            dlen = np.zeros_like(slit_length)
            dlen[:,indx] = np.log(np.absolute(slit_length[:,indx]/med_slit_length[None,indx]))
            indx = np.repeat(np.logical_not(indx) 
                                | np.any(np.absolute(dlen) > np.log(1+dlength_rtol), axis=0), 2)
            self._flag_edges(trace_cen, indx, 'LARGELENGTHCHANGE')

        if length_atol is not None:
            # Find any short slits (flag both edges of the slit)
            indx = np.repeat(med_slit_length < length_atol, 2)
            self._flag_edges(trace_cen, indx, 'SHORTSLIT')

        if length_atol_sci is not None:
            # Find any box slits (flag both edges of the slit)
            indx = med_slit_length < length_atol_sci
            if length_atol is not None:
                indx &= (med_slit_length > length_atol)
            indx = np.repeat(indx, 2)
            self._flag_edges(trace_cen, indx, 'BOXSLIT')

        # TODO: Might want to do this first, then do overlap, then redo all of
        # the above after rejiggering the overlap regions.
        if length_rtol is not None:
            # Find slits that are abnormally short
            short = np.repeat(np.log(med_slit_length/np.median(med_slit_length))
                              < np.log(1-length_rtol), 2)
            if np.any(short):
                msgs.info(f'Flagging {np.sum(short)} abnormally short slit edges.')
                self.edge_msk[:,short] \
                        = self.bitmask.turn_on(self.edge_msk[:,short], 'ABNORMALSLIT_SHORT')
            long = np.repeat(np.log(med_slit_length/np.median(med_slit_length))
                             > np.log(1+length_rtol), 2)
            if np.any(long):
                msgs.info(f'Flagging {np.sum(long)} abnormally long slit edges.')
                self.edge_msk[:,long] \
                        = self.bitmask.turn_on(self.edge_msk[:,long], 'ABNORMALSLIT_LONG')

        # TODO: Consider removing slits that have large length changes.  Like so:
#        # Remove traces that have significant changes in their spatial extent
#        # along the dispersion direction.
#        dl_flag = self.fully_masked_traces(flag='LARGELENGTHCHANGE')
#        if np.any(dl_flag):
#            msgs.info(f'Removing {np.sum(dl_flag)} traces because of large spatial extent '
#                      ' changes along the dispersion direction.')
#            self.remove_traces(dl_flag, rebuild_pca=_rebuild_pca)

        # Try to detect overlap between adjacent slits by finding abnormally
        # short slits.
        #   - Find abnormally short slits that *do not* include inserted edges;
        #     i.e., these must be *detected* edges, not inserted ones.
        #   - *Both* edges in the fit must be flagged because of this
        #     requirement that the trace not be inserted.  This means that we
        #     set mode='neither' when running synced_selection.  I also set
        #     assume_synced=True: the traces should be synced if the code has
        #     made it this far.  Any flags that would indicate otherwise will
        #     have been set by this function.
        short = self.fully_masked_traces(flag='ABNORMALSLIT_SHORT',
                                         exclude=self.bitmask.insert_flags)
        short = self.synced_selection(short, mode='neither', assume_synced=True)
        if self.par['overlap'] and np.any(short):
            msgs.info('Assuming slits flagged as abnormally short are actually due to '
                      'overlapping slit edges.')
            rmtrace = np.zeros(self.ntrace, dtype=bool)
            # Find sets of adjacent short slits and assume they all select
            # adjacent overlap regions.
            short_slits = utils.contiguous_true(short)
            sync_inserts = self.fully_masked_traces(flag='SYNCINSERT')
            for slc in short_slits:
                # Remove the edges just before and after this region of short
                # slits, if they inserted by the left-right syncing.
                rmtrace[max(0,slc.start-1)] = sync_inserts[max(0,slc.start-1)]
                rmtrace[min(self.ntrace-1, slc.stop)] = sync_inserts[min(self.ntrace-1, slc.stop)]
                # Flip the sign of the edges (i.e., turn lefts into rights and
                # vice versa).  This turns the overlap regions into slit gaps.
                self.traceid[slc] *= -1
                # Turn off the masking
                self.edge_msk[:,slc] \
                        = self.bitmask.turn_off(self.edge_msk[:,slc], 'ABNORMALSLIT_SHORT')

            # Remove the flagged traces, resort the edges, and rebuild the pca
            self.remove_traces(rmtrace, rebuild_pca=_rebuild_pca)

            # If this de-synchronizes the traces, we effectively have to start
            # the synchronization process over again, with the adjustments for
            # the "short" slits that are assumed to be overlap regions.
            if not self.is_synced:
                msgs.info('Checking/cleaning traces for overlap led to de-syncronization.')
                return False

        # TODO: Check that slit edges meet a minimum slit gap?

        # Find all traces to remove
        rmtrace = self.fully_masked_traces(flag=self.bitmask.bad_flags,
                                           exclude=self.bitmask.exclude_flags)
        # Make sure to also remove the synced one
        rmtrace = self.synced_selection(rmtrace, mode='both', assume_synced=True)
        # Remove 'em
        self.remove_traces(rmtrace, rebuild_pca=_rebuild_pca)
        if self.is_empty:
            msgs.warn('Assuming a single long-slit and continuing.')
            self.bound_detector()
        return True

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
                # JFH Print out in the spec:spat format that corresponds to the pypeit file
                msgs.info("Removing user-supplied slit at spec:spat {}:{}".format(y_spec, xcen))
                # Mask
                self.bitmask.turn_on(self.edge_msk[:,indx], 'USERRMSLIT')

        # Syncronize the selection
        indx = self.synced_selection(indx, mode='both', assume_synced=True)
        # TODO: Add rebuild_pca ?
        # Remove
        self.remove_traces(indx)

    # TODO:
    #   - Add an option to distinguish between an actual remove and a flagging
    #   - Allow traces to be removed but keep the PCA if it exists instead of
    #     rebuilding it?
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
                If the pca exists, rebuild it using the new traces and the
                previous parameter set.  If False, any existing PCA model is
                removed.
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
        if self.maskdef_id is not None:
            self.maskdef_id = self.maskdef_id[keep]

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
                       'been left-right synchronized.  Either run sync() to sychronize or ignore '
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
        if self.maskdef_id is not None:
            self.maskdef_id = self.maskdef_id[srt]

        # Reorder the trace numbers
        indx = self.traceid < 0
        self.traceid[indx] = -1-np.arange(np.sum(indx))
        indx = np.logical_not(indx)
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
        (:class:`~pypeit.par.pypeitpar.EdgeTracePar`) are
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
        bpm = self.traceimg.select_flag(flag='BPM') if self.tracebpm is None else self.tracebpm

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
                        & np.logical_not(np.all(edge_bpm, axis=0))
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
        (:class:`~pypeit.par.pypeitpar.EdgeTracePar`) are
        ``fit_min_spec_length``, ``left_right_pca``, and ``pca_min_edges``.

        .. warning::

            This function calls :func:`check_traces` using
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
        msgs.info('Minimum length of traces to include in the PCA: {0}'.format(minimum_spec_length))

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
        good = np.sum(np.logical_not(self.bitmask.flagged(self.edge_msk)), axis=0) > 0

        # Returned value depends on whether or not the left and right
        # traces are done separately
        return np.sum(good[self.is_left]) > self.par['pca_min_edges'] \
                    and np.sum(good[self.is_right]) > self.par['pca_min_edges'] \
                    if self.par['left_right_pca'] else np.sum(good) > self.par['pca_min_edges']

    # TODO: Consolidate the options in `add_user_traces` with this function to
    # enable more prediction options, not just the PCA.
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
            trace_add = np.zeros((self.nspec,_side.size), dtype='float')
            for s,p in zip([-1,1], [self.left_pca,self.right_pca]):
                indx = _side == s
                if not np.any(indx):
                    continue
                trace_add[:,indx] = p.predict(_edge_cen[indx])
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
        (:class:`~pypeit.par.pypeitpar.EdgeTracePar`) are
        `fit_min_spec_length`, `left_right_pca`, `pca_n`,
        `pca_var_percent`, `pca_function`, `pca_order`, `pca_sigrej`,
        `pca_maxrej`, and `pca_maxiter`.

        Args:
            use_center (:obj:`bool`, optional):
                Use the center measurements for the PCA decomposition
                instead of the functional fit to those data. This is
                only relevant if both are available. If no fits have
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
        use_trace = np.sum(np.logical_not(bpm), axis=0) > 0

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
        self.edge_fit = self.predict_traces(trace_ref, side=side)

        # TODO: Compare with the fit data. Remove traces where the mean
        # offset between the PCA prediction and the measured centroids
        # are larger than some threshold?

        # Log what was done
        self.log += [inspect.stack()[0][3]]

    def peak_refine(self, rebuild_pca=False, debug=False):
        """
        Refine the trace by isolating peaks and troughs in the Sobel-filtered
        image.

        This function *requires* that the PCA model exists; see
        :func:`build_pca` or :func:`pca_refine`. It is also primarily a wrapper
        for :func:`~pypeit.core.trace.peak_trace`. See the documentation of that
        function for the explanation of the algorithm.

        If the left and right traces have separate PCA decompositions, this
        function makes one call to :func:`~pypeit.core.trace.peak_trace` for
        each side.  Otherwise, a single call is made to
        :func:`~pypeit.core.trace.peak_trace` where both the peak and troughs in
        :attr:`sobelsig` are detected and traced.

        Optionally, the code will match and compare the traces found and fit by
        :func:`~pypeit.core.trace.peak_trace` to the original traces.  If the
        RMS difference between the matched traces is large, they can be removed
        (see ``trace_rms_tol`` in :class:`~pypeit.par.pypeitpar.EdgeTracePar`).

        Note that this effectively reinstantiates much of the object attributes,
        including :attr:`traceid`, :attr:`edge_cen`, :attr:`edge_err`,
        :attr:`edge_msk`, :attr:`edge_img`, :attr:`edge_fit`, and
        :attr:`fittype`.

        Used parameters from :attr:`par`
        (:class:`~pypeit.par.pypeitpar.EdgeTracePar`) are ``left_right_pca``,
        ``edge_thresh``, ``smash_range``, ``edge_detect_clip``,
        ``trace_median_frac``, ``trace_thresh``, ``trace_rms_tol``,
        ``fit_function``, ``fit_order``, ``fwhm_uniform``, ``niter_uniform``,
        ``fwhm_gaussian``, ``niter_gaussian``, ``fit_maxdev``, and
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
        bpm = self.traceimg.select_flag(flag='BPM') if self.tracebpm is None else self.tracebpm

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
            _sobelsig = trace.prepare_sobel_for_trace(self.sobelsig, bpm=bpm, boxcar=5, side=None)

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
            
        if self.par['trace_rms_tol'] is not None:
            # Get the PCA reference row.  The PCA *must* have been defined to get
            # this far (see pcatype test at the beginning of the function).
            reference_row = self.left_pca.reference_row if self.par['left_right_pca'] \
                                else self.pca.reference_row
            # Match the new trace positions to the input ones
            gpm = np.logical_not(self.bitmask.flagged(self.edge_msk[reference_row]))
            # TODO: Include a tolerance here?  Needs testing with more datasets.
            peak_indx = slitdesign_matching.match_positions_1D(
                            self.edge_fit[reference_row][gpm],
                            fit[reference_row])

            # Determine the RMS difference between the input and output traces.
            # This allows us to compare traces that had already been identified
            # to their new measurements resulting from peak_trace, and remove
            # them if they are too discrepant from their original form.  This is
            # largely meant to find and remove poorly constrained traces, where
            # the polynomial fit goes wonky.
            diff = fit - self.edge_fit.T[gpm][peak_indx].T
            rms = np.sqrt(np.mean((diff - np.mean(diff, axis=0)[None,:])**2, axis=0))

            # Report
            msgs.info('-'*30)
            msgs.info('Matched spatial locations and RMS difference along spectral direction')
            msgs.info(f' {"OLD":>8} {"NEW":>8} {"RMS":>8}')
            msgs.info(' '+'-'*8+' '+'-'*8+' '+'-'*8)
            for i in range(len(peak_indx)):
                if peak_indx[i] < 0:
                    continue
                msgs.info(f' {self.edge_fit[reference_row][gpm][peak_indx[i]]:8.1f}'
                          f' {fit[reference_row][i]:8.1f} {rms[i]:8.3f}')

            # Select traces below the RMS tolerance or that were newly
            # identified by peak_trace.  I.e., this will *not* catch newly
            # identified traces found by peak_trace that are also poorly
            # constrained!
            indx = (rms < self.par['trace_rms_tol']) | (peak_indx == -1)
            if not np.all(indx):
                msgs.info(f'Removing {indx.size - np.sum(indx)} trace(s) due to large RMS '
                          'difference with previous trace locations.')
                fit = fit[:,indx]
                cen = cen[:,indx]
                err = err[:,indx]
                msk = msk[:,indx]
                nleft -= np.sum(np.where(np.logical_not(indx))[0] < nleft)

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
        (:class:`~pypeit.par.pypeitpar.EdgeTracePar`) are
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
        trace_ref[np.logical_not(add_edge)] = trace_cen[reference_row,:]
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
        (:class:`~pypeit.par.pypeitpar.EdgeTracePar`). No limit is
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
        # Check vector size
        if trace_cen.shape[0] != self.nspec:
            msgs.error('Traces have incorrect length.')
        _buffer = self.par['det_buffer']
        if _buffer < 0:
            msgs.warn('Buffer must be greater than 0; ignoring.')
            _buffer = 0

        if self.par['max_nudge'] is not None:
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
        (:class:`~pypeit.par.pypeitpar.EdgeTracePar`) are
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

        Returns:
            :obj:`bool`: Returns the status of the syncing. True means
            success.
        """
        # Remove any fully masked traces. Keeps any inserted or
        # box-slit traces.
        self.clean_traces(rebuild_pca=rebuild_pca)

        # Make sure there are still traces left
        if self.is_empty:
            if not self.par['bound_detector']:
                return False
            msgs.warn('No traces left!  Left and right edges placed at detector boundaries.')
            self.bound_detector()

        # Make sure that the traces are sorted spatially
        self.spatial_sort()

        # If the traces are already synced, check them and log the
        # function as completed
        if self.is_synced \
                and self.check_synced(rebuild_pca=rebuild_pca and self.pcatype is not None):
            if self.log is not None:
                self.log += [inspect.stack()[0][3]]
            return True

        # Edges are currently not synced, so check the input
        if self.par['sync_predict'] not in ['pca', 'nearest', 'auto']:
            msgs.error('Unknown trace mode: {0}'.format(self.par['sync_predict']))
        if self.par['sync_predict'] == 'pca' and self.pcatype is None:
            msgs.error('The PCA decomposition does not exist.  Either run self.build_pca or use '
                       'a different trace_mode.')

        # Find the edges to add, what side they're on, and where to
        # insert them into the existing trace array
        side, add_edge, add_indx = self._get_insert_locations()
        if not np.any(add_edge):
            # No edges to add
            return True

        # Report
        msgs.info('-'*50)
        msgs.info('{0:^50}'.format('Synchronizing left and right traces'))
        msgs.info('-'*50)
        msgs.info('Found {0} left and {1} right trace(s) to add.'.format(
                    np.sum((side == -1) & add_edge), np.sum((side == 1) & add_edge)))

        # Allow the edges to be synced, even if a fit hasn't been done yet
        trace_cen = self.edge_cen if self.edge_fit is None else self.edge_fit

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
            trace_add = np.zeros((self.nspec, np.sum(add_edge)), dtype=float)
            trace_add[:,0] = trace_cen[:,0] + offset
            self.insert_traces(side[add_edge], trace_add, loc=add_indx[add_edge], mode='sync')
            return True

        # NOTE: The maximum number of iterations is hard-coded for now.  Testing
        # is needed to know if we need this to be any larger than 3.
        maxiter = 3
        i = 0
        while i < maxiter:
            msgs.info(f'Beginning syncing iteration : {i+1} (of at most {maxiter})')

            # Get the traces
            trace_cen = self.edge_cen if self.edge_fit is None else self.edge_fit

            # Find the edges to add, what side they're on, and where to insert
            # them into the existing trace array.  This is done again, in case
            # we're going through another iteration
            side, add_edge, add_indx = self._get_insert_locations()

            # Get the reference locations for the new edges
            trace_ref = self._get_reference_locations(trace_cen, add_edge)

            # Determine which sync_predict to use
            if self.par['sync_predict'] == 'pca' \
                    or (self.par['sync_predict'] == 'auto' and self.can_pca()):
                _sync_predict = 'pca'
            else:
                _sync_predict = 'nearest'

            # Predict the traces either using the PCA or using the nearest slit edge
            if _sync_predict == 'pca':
                trace_add = self.predict_traces(trace_ref[add_edge], side=side[add_edge])
            elif _sync_predict == 'nearest':
                # Index of trace nearest the ones to add
                # TODO: Force it to use the nearest edge of the same side;
                # i.e., when inserting a new right, force it to use the
                # nearest right instead of the nearest left?
                nearest = utils.nearest_unmasked(np.ma.MaskedArray(trace_ref, mask=add_edge))
                # Indices of the original traces
                indx = np.zeros(len(add_edge), dtype=int)
                indx[np.logical_not(add_edge)] = np.arange(self.ntrace)
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
                msgs.error('Catastrophic error in left-right synchronization.  Edge order is not '
                           'correctly sorted.')
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
                self.show(title='includes inserted traces before checking the sync', flag='any')

            # Check the full synchronized list and log completion of the
            # method
            if self.check_synced(rebuild_pca=rebuild_pca):
                break

            i += 1
            if i == maxiter:
                msgs.error('Fatal left-right trace de-synchronization error.')

        if self.log is not None:
            self.log += [inspect.stack()[0][3]]
        return True

    def add_user_traces(self, user_traces, method='straight'):
        """
        Add traces for user-defined slits.

        Args:
            user_slits (:obj:`list`):
                A list of lists with the coordinates for the new traces.  For
                each new slit, the list provides the spectral coordinate at
                which the slit edges are defined and the left and right spatial
                pixels that the traces should pass through.  I.e., ``[664, 323,
                348]`` mean construct a left edge that passes through pixel
                ``(664,323)`` (ordered spectral then spatial) and a right edge
                that passes through pixel ``(664,348)``.

            method (:obj:`str`, optional):
                The method used to construct the traces.  Options are:

                    - ``'straight'``: Simply insert traces that have a constant
                      spatial position as a function of spectral pixel.

                    - ``'nearest'``: Constrain the trace to follow the same form
                      as an existing trace that is nearest the provided new
                      trace coordinates.

                    - ``'pca'``: Use the PCA decomposition of the traces to
                      predict the added trace.  If the PCA does not currently
                      exist, the function will try to (re)build it.

        """
        #msgs.info("Adding new slits at x0, x1 (left, right)".format(x_spat0, x_spat1))
        # Number of added slits
        n_add = len(user_traces)
        # Add two traces for each slit, one left and one right
        side = np.tile([-1,1], (1,n_add)).ravel()
        # Reformat the user-defined input into an array of spectral and spatial
        # coordiates for the new traces
        new_trace_coo = np.array(user_traces, dtype=float)
        new_trace_coo = np.insert(new_trace_coo, 2, new_trace_coo[:,0], axis=1).reshape(-1,2)
        if method == 'straight':
            # Just repeat the spatial positions
            new_traces = np.tile(new_trace_coo[:,1], (self.nspec,1))
        elif method == 'nearest':
            if self.is_empty:
                msgs.error('No edge traces currently exist.  Cannot insert user slits with a '
                           'shape based on the nearest existing slit edges!  '
                           'Set add_predict = straight.')
            # Use the measured edges if the functional forms don't exist (yet)
            trace_cen = self.edge_cen if self.edge_fit is None else self.edge_fit
            # Find the trace nearest to the one to be inserted
            nearest = np.array([np.argmin(np.absolute(trace_cen[int(coo[0])] - coo[1]))
                                    for coo in new_trace_coo])
            # Construct the new traces by offseting the existing ones
            new_traces = trace_cen[:,nearest] + new_trace_coo[:,1] \
                            - trace_cen[new_trace_coo[:,0].astype(int), nearest]
        elif method == 'pca':
            if self.is_empty:
                msgs.error('No edge traces currently exist.  Cannot insert user slits with a '
                           'shape based on the PCA decomposition of the existing slit edges!  '
                           'Set add_predict = straight.')
            if self.pcatype is None:
                if not self.can_pca():
                    msgs.error('PCA does not exist and cannot be constructed!  Cannot insert user '
                               'slits with a shape based on the PCA decomposition of the existing '
                               'slit edges!  Use add_predict = straight or nearest.')
                self.build_pca()
            # Use the pca to predict the traces at the requested spatial positions
            trace_ref = new_trace_coo[:,1]
            new_traces = self.predict_traces(trace_ref, side=side)
            # NOTE: The users can pick spectral rows in the definition of the
            # new slits that are not the same as the reference row at which the
            # PCA is defined.  This iteration in the prediction tries to adjust
            # the reference spatial position for the prediction so that the
            # predicted trace runs through the user-specific spatial position at
            # the user-specified spectral position.
            trace_ref += new_trace_coo[:,1] \
                            - new_traces[new_trace_coo[:,0].astype(int),np.arange(n_add*2)]
            new_traces = self.predict_traces(trace_ref, side=side)
        else:
            msgs.error(f'Unknown method for adding user slit: {method}')

        # Insert
        self.insert_traces(side, new_traces, mode='user')
        # Sync
        self.check_synced(rebuild_pca=False)

    def insert_traces(self, side, trace_cen, loc=None, mode='user', resort=True, nudge=True):
        r"""
        Insert/append a set of edge traces.

        New traces to add are first nudged away from the detector
        edge (see :func:`nudge_traces`) according to parameters
        `max_nudge` and `det_buffer` from :attr:`par`
        (:class:`~pypeit.par.pypeitpar.EdgeTracePar`). They are then
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
                    - ``'order'``: Traces are the expected location of an
                      echelle order.

            resort (:obj:`bool`, optional):
                Resort the traces in the spatial dimension; see
                :func:`spatial_sort`.
            nudge (:obj:`bool`, optional):
                Allow the traces to be nudged away from the detector
                edge according to :attr:`par` and
                :func:`nudge_traces`.

        """
        # TODO: When inserting traces and echelle orders are already matched,
        # the length of the orderid vector is not longer valid.  For now, just
        # remove any existing array and warn the user they they'll need to
        # rematch the orders.
        if self.orderid is not None:
            msgs.warn('Inserting traces invalidates order matching.  Removing.')
            self.orderid = None

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

        msgs.info(f'Inserting {ntrace} new traces.')

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
        elif mode == 'order':
            mask = self.bitmask.turn_on(mask, 'ORDERINSERT')

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
        if self.maskdef_id is not None:
            self.maskdef_id = np.insert(self.maskdef_id, loc, -99)  # we don't know the slitmask id

        if resort:
            self.spatial_sort()

    def bound_detector(self):
        """
        Insert traces at both detector boundaries.

        Accounting for the requested detector buffer, a left and
        right trace are placed at the detector edges. Traces are
        masked as user-inserted.

        Only used parameter from :attr:`par`
        (:class:`~pypeit.par.pypeitpar.EdgeTracePar`) is
        `det_buffer`.
        """
        # find where sobelsig is not masked
        sobelsig_nomask = np.where(self.sobelsig[int(self.nspec/2.), :] != 0)[0]

        self.insert_traces(np.array([-1,1]),
                           np.array([np.full(self.nspec, sobelsig_nomask[0]+self.par['det_buffer'], dtype='float'),
                                     np.full(self.nspec, sobelsig_nomask[-1]-self.par['det_buffer'],
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
                used. See :func:`pypeit.bitmask.BitMask.flagged`.
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
        # TODO - break up this method (too long)
        Match slit info from the mask design data to the traced slits.

        Use of this method requires:
            - a PCA decomposition is available,
            - :attr:`spectrograph` has a viable `get_maskdef_slitedges` method
              to read edge trace locations predicted by the slitmask design. This data can be pulled
              from one of the files used to construct the trace image.

        The method uses a collection of scripts in pypeit.core.slitdesign_matching
        which are taken from DEEP2 IDL-based pipeline for DEIMOS data.

        NOTE this method was first written to deal with edge traces predicted by DEIMOS optical model, but
        subsequently adapted to be used with other spectrographs that don't use optical models.



        Args:
            debug (:obj:`bool`, optional):
                Run in debug mode.
        """

        # Remove any fully masked traces. Keeps any inserted or
        # box-slit traces.
        self.clean_traces(rebuild_pca=True)

        # Check that there are still traces to match!
        if self.is_empty:
            msgs.warn('No edges traced. Slitmask matching cannot be performed')
            return

        # `traceimg` must have knowledge of the flat frame that built it
        maskfiles = self.traceimg.files[0] if self.par['maskdesign_filename'] is None \
            else self.par['maskdesign_filename']
        self.maskfile = maskfiles[0] if isinstance(maskfiles, list) else maskfiles
        omodel_bspat, omodel_tspat, sortindx, self.slitmask = \
            self.spectrograph.get_maskdef_slitedges(
                ccdnum=self.traceimg.detector.det, 
                binning=self.traceimg.detector.binning, 
                filename=maskfiles, 
                trc_path = os.path.dirname(self.traceimg.files[0]),
                debug=debug)

        if omodel_bspat[omodel_bspat!=-1].size < 3:
            msgs.warn('Less than 3 slits are expected on this detector, slitmask matching cannot be performed')
            # update minimum_slit_gap and minimum_slit_length_sci par
            # this will allow to catch the boxslit, since in this case slitmask matching is not performed
            self.par = self.spectrograph.update_edgetracepar(self.par)
            return

        # reference row
        bpm = self.bitmask.flagged(self.edge_msk, self.bitmask.bad_flags)
        # TODO make reference row an attribute of EdgeTraceSet
        reference_row = trace.most_common_trace_row(bpm) if self.pcatype is None \
            else (self.left_pca.reference_row if self.par['left_right_pca'] else self.pca.reference_row)
        spat_bedge = self.edge_fit[reference_row, self.is_left]
        spat_tedge = self.edge_fit[reference_row, self.is_right]

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
            plt.show()
        msgs.info('SLIT_MATCH: RMS residuals for left and right edges: {}, {} pixels'.format(sigres_b, sigres_t))


        # We compute the predicted edge positions from the optical model after the x-correlation with the traced edges
        # bottom edges
        bot_edge_pred = omodel_bspat.copy()
        # Predictions that are outside the detector have values = -1.
        bot_edge_pred[omodel_bspat!=-1] = coeff_b[0] + coeff_b[1] * omodel_bspat[omodel_bspat!=-1] if not switched else\
                                          coeff_b[0] + coeff_b[1] * omodel_tspat[omodel_bspat!=-1]
        # top edges
        top_edge_pred = omodel_tspat.copy()
        # Predictions that are outside the detector have values = -1.
        top_edge_pred[omodel_tspat!=-1] = coeff_t[0] + coeff_t[1]*omodel_tspat[omodel_tspat!=-1] if not switched else\
                                          coeff_t[0] + coeff_t[1]*omodel_bspat[omodel_tspat!=-1]

        # Find if there are missing traces.
        # Need exactly one occurrence of each index in "need"
        buffer = 3*self.par['det_buffer'] + 1
        need = (top_edge_pred > buffer) & (bot_edge_pred < (self.traceimg.shape[1] - 1 - buffer)) & \
               ((omodel_bspat != -1) | (omodel_tspat != -1))

        # bottom edges
        needadd_b = need.copy()
        needadd_b[ind_b] = False
        needind_b = np.where(needadd_b)[0]  # edges we are missing

        # top edges
        needadd_t = need.copy()
        needadd_t[ind_t] = False
        needind_t = np.where(needadd_t)[0]  # edges we are missing

        if (needind_b.size > 0) | (needind_t.size > 0):
            msgs.warn('Missing edge traces: {} left and {} right'.format(needind_b.shape[0], needind_t.shape[0]))

        if debug:
            slitdesign_matching.plot_matches(self.edge_fit[:,self.is_left], ind_b, bot_edge_pred, reference_row,
                                             self.slitmask.slitindx, nspat=self.nspat, duplicates=dupl_b,
                                             missing=needind_b, edge='left')
            slitdesign_matching.plot_matches(self.edge_fit[:,self.is_right], ind_t, top_edge_pred, reference_row,
                                             self.slitmask.slitindx, nspat=self.nspat, duplicates=dupl_t,
                                             missing=needind_t, edge='right')

        # Put duplicate flags for left and right traces together in one array
        if np.any(dupl_b) or np.any(dupl_t):
            edges_dupl = np.zeros(self.ntrace, dtype=bool)
            if np.any(dupl_b):
                edges_dupl[self.is_left] = dupl_b
            if np.any(dupl_t):
                edges_dupl[self.is_right] = dupl_t

            # Remove duplicate match
            msgs.info('Removing duplicate matches: {} left and {} right'.format(ind_b[dupl_b].size, ind_t[dupl_t].size))
            self.remove_traces(edges_dupl, rebuild_pca=True)
            ind_b = ind_b[np.logical_not(dupl_b)]
            ind_t = ind_t[np.logical_not(dupl_t)]

        # RE-CHECK for missing traces after removing the duplicates
        # bottom edges
        needadd_b = need.copy()
        needadd_b[ind_b] = False
        needind_b = np.where(needadd_b)[0]  # edges we are missing

        # top edges
        needadd_t = need.copy()
        needadd_t[ind_t] = False
        needind_t = np.where(needadd_t)[0]  # edges we are missing

        #  The code below is to add traces that are predicted but not found.
        # Left edges
        if needind_b.size > 0:
            msgs.info('Adding {} left missing edge(s)'.format(needind_b.size))
            # Append the missing indices and re-sort all
            ind_b = np.append(ind_b, needind_b)
            sortind_b = np.argsort(utils.index_of_x_eq_y(self.slitmask.slitid[sortindx],
                                                         self.slitmask.slitid[ind_b], strict=True), kind='stable')
            ind_b = ind_b[sortind_b]
            for i in range(bot_edge_pred[needind_b].size):
                # check if the trace that will be added is off detector
                if bot_edge_pred[needind_b][i] < self.par['det_buffer']:
                    # add a trace parallel to the detector edge (i.e., no pca)
                    missing_left_traces = np.full(self.nspec, self.par['det_buffer'])
                    self.insert_traces(-1, missing_left_traces, mode='mask', nudge=False)
                else:
                    # We try to ensure that the left edge is inserted after the right edge of previous slit
                    # Find index of larger closest trace
                    indx_larger = np.where(self.edge_fit[reference_row, :] > bot_edge_pred[needind_b][i])[0]
                    if (indx_larger.size > 0) & (indx_larger.size < self.traceid.size):
                        indx = indx_larger[0]
                        # If close trace is "right" but it's too close to the new left trace,
                        # add the new trace after "right" one
                        if (self.traceid[indx] > 0) and \
                                ((self.edge_fit[reference_row, indx] - bot_edge_pred[needind_b][i]) < 5):
                            bot_edge_pred[needind_b[i]] = self.edge_fit[reference_row, indx] + 1
                    if self.pcatype is None:
                        # deal with the case when the pca is not available
                        nearest = np.argmin(np.absolute(self.edge_fit[reference_row, :] - bot_edge_pred[needind_b][i]))
                        missing_left_traces = self.edge_fit[:, nearest]-(self.edge_fit[reference_row, nearest] - bot_edge_pred[needind_b][i])
                    else:
                        missing_left_traces = self.predict_traces(bot_edge_pred[needind_b][i], side=-1)
                    self.insert_traces(-1, missing_left_traces, mode='mask', nudge=False)

        if needind_t.size > 0:
            # Right edges
            msgs.info('Adding {} right missing edge(s)'.format(needind_t.size))
            # Append the missing indices and re-sort all
            ind_t = np.append(ind_t, needind_t)
            sortind_t = np.argsort(utils.index_of_x_eq_y(self.slitmask.slitid[sortindx],
                                                         self.slitmask.slitid[ind_t], strict=True), kind='stable')
            ind_t = ind_t[sortind_t]
            for i in range(top_edge_pred[needind_t].size):
                # check if the trace that will be added is off detector
                if top_edge_pred[needind_t][i] > self.nspat - self.par['det_buffer']:
                    # add a trace parallel to the detector edge (i.e., no pca)
                    missing_right_traces = np.full(self.nspec, self.nspat - self.par['det_buffer'])
                    self.insert_traces(1, missing_right_traces, mode='mask', nudge=False)
                else:
                    # We try to ensure that the right edge is inserted before the left edge of following slit
                    # Find index of smaller closest trace
                    indx_smaller = np.where(self.edge_fit[reference_row, :] < top_edge_pred[needind_t][i])[0]
                    if (indx_smaller.size > 0) & (indx_smaller.size < self.traceid.size):
                        indx = indx_smaller[-1]
                        # If close trace is "left" but it's too close to the new left trace,
                        # add new trace before "left" one
                        if (self.traceid[indx] < 0) and \
                                ((top_edge_pred[needind_t][i] - self.edge_fit[reference_row, indx]) < 5):
                            top_edge_pred[needind_t[i]] = self.edge_fit[reference_row, indx] - 1
                    if self.pcatype is None:
                        # deal with the case when the pca is not available
                        nearest = np.argmin(np.absolute(self.edge_fit[reference_row, :] - top_edge_pred[needind_t[i]]))
                        missing_right_traces = self.edge_fit[:, nearest]-(self.edge_fit[reference_row, nearest] - top_edge_pred[needind_t[i]])
                    else:
                        missing_right_traces = self.predict_traces(top_edge_pred[needind_t][i], side=1)
                    self.insert_traces(1, missing_right_traces, mode='mask', nudge=False)

        if debug:
            slitdesign_matching.plot_matches(self.edge_fit[:, self.is_left], ind_b, bot_edge_pred, reference_row,
                                                 self.slitmask.slitindx, nspat=self.nspat, edge='left')
            slitdesign_matching.plot_matches(self.edge_fit[:, self.is_right], ind_t, top_edge_pred, reference_row,
                                                 self.slitmask.slitindx, nspat=self.nspat, edge='right')

        self.maskdef_id = np.zeros(self.ntrace, dtype=int)
        self.maskdef_id[self.is_left] = self.slitmask.slitid[ind_b]
        self.maskdef_id[self.is_right] = self.slitmask.slitid[ind_t]

        # flag as 'BOXSLIT' the edges that are predicted to be alignment boxes from optical model.
        align_slit = np.zeros(self.ntrace, dtype=bool)
        align_slit[self.is_left] = self.slitmask.alignment_slit[ind_b]
        align_slit[self.is_right] = self.slitmask.alignment_slit[ind_t]
        self.edge_msk[:, align_slit] = self.bitmask.turn_on(self.edge_msk[:, align_slit], 'BOXSLIT')

        # Propagate the coefficients, `coeff_b` and `coeff_t`, of the x-correlation and the
        # left and right spatial position of the slit edges from optical model (before x-correlation)
        # with the purpose to fill a table with the information on slitmask design matching. The table will
        # be filled out at the very end of the slit tracing process, in `get_slits()`.
        self.cc_params_b = coeff_b[0], coeff_b[1], sigres_b
        self.cc_params_t = coeff_t[0], coeff_t[1], sigres_t
        self.omodel_bspat = omodel_bspat
        self.omodel_tspat = omodel_tspat

        # Make sure the traces are synchronized
        if self.traceid[-1] < 0:
            self.edge_msk[:, -1] = self.bitmask.turn_on(self.edge_msk[:, -1], 'OFFDETECTOR')
        if self.traceid[0] > 0:
            self.edge_msk[:, 0] = self.bitmask.turn_on(self.edge_msk[:, 0], 'OFFDETECTOR')

        # sync
        self.sync()

        # remove traces with a mismatch in the maskdef_id (it's better to remove the traces
        # rather than propagating the mismatch)
        diff_maskdef_id = np.repeat(self.maskdef_id[self.is_left] - self.maskdef_id[self.is_right] != 0, 2)
        if np.any(diff_maskdef_id):
            self.edge_msk[:, diff_maskdef_id] = self.bitmask.turn_on(self.edge_msk[:, diff_maskdef_id], 'SYNCERROR')

        # remove short slits
        platescale = parse.parse_binning(self.traceimg.detector.binning)[1] * self.traceimg.detector['platescale']
        # minimum slit length
        min_slitlen = self.par['minimum_slit_length'] / platescale if \
            self.par['minimum_slit_length'] is not None else 20.
        slit_length = np.median(np.squeeze(np.diff(self.edge_fit.reshape(self.nspec, -1, 2), axis=-1)), axis=0)
        ind_short = np.repeat(slit_length < min_slitlen, 2)
        if np.any(ind_short):
            msgs.info('Rejecting {0} traces that are too short.'.format(np.sum(ind_short)))
            self.edge_msk[:, ind_short] = self.bitmask.turn_on(self.edge_msk[:, ind_short], 'SHORTSLIT')

        # force clean_traces for certain flags, because if a slit has one of these flags
        # but has also a bitmask.exclude_flags it will not be removed
        if not np.all(self.bitmask.flagged(self.edge_msk, self.bitmask.bad_flags)):
            self.clean_traces(force_flag=['SYNCERROR', 'OFFDETECTOR', 'SHORTSLIT'], rebuild_pca=True,
                              sync_mode='both', assume_synced=True)

        if self.is_synced:
            msgs.info('LEFT AND RIGHT EDGES SYNCHRONIZED AFTER MASK DESIGN MATCHING')
        else:
            msgs.warn('LEFT AND RIGHT EDGES *NOT* SYNCHRONIZED AFTER MASK DESIGN MATCHING')

    def _fill_design_table(self, maskdef_id, cc_params_b, cc_params_t, omodel_bspat, omodel_tspat, spat_id):
        """
        Fill :attr:`design` based on the results of the design
        registration.

        The :attr:`design` is an `astropy.table.Table`_ with 13 columns:
            - 'TRACEID': Trace ID Number
            - 'TRACESROW': Spectral row for provided left and right edges
            - 'TRACELPIX': Spatial pixel coordinate for left edge
            - 'TRACERPIX': Spatial pixel coordinate for right edge
            - 'MASKDEF_ID': Slit ID Number from slit-mask design
            - 'SLITLMASKDEF': Left edge of the slit in pixel from slitmask design before x-correlation
            - 'SLITRMASKDEF': Right edge of the slit in pixel from slitmask design before x-correlation
            - 'SLITRA': Right ascension of the slit center (deg)
            - 'SLITDEC': Declination of the slit center (deg)
            - 'SLITLEN': Slit length (arcsec)
            - 'SLITWID': Slit width (arcsec)
            - 'SLITPA': Slit position angle on sky (deg from N through E)
            - 'ALIGN': Slit used for alignment (1-yes; 0-no), not target observations.
        And three `.meta` info:
            - 'MASKFILE': name of file with the slitmask info
            - 'MASKOFFL', 'MASKOFFR': The coefficient 'offset' of the x-correlation between edges predicted by
                                      the slitmask design and the one traced on the image. One value per
                                      each edge side.
            - 'MASKSCLL', 'MASKSCLR': The coefficient 'scale' of the x-correlation between edges predicted by
                                      the slitmask design and the one traced on the image. One value per
                                      each edge side
            - 'MASKRMSL', 'MASKRMSR': The RMS of the x-correlation between edges predicted by the slitmask design
                                      and the one traced on the image. One value per each edge side

        Args:
            maskdef_id (`numpy.ndarray`_):
                Slit ID number from slit-mask design matched to traced slits.
            cc_params_b, cc_params_t (:obj:`tuple`):
                Three parameters of the cross-correlation (2 coefficients and RMS) between slit-mask design
                and traced edges for the left and right edges.
            omodel_bspat, omodel_tspat (`numpy.ndarray`_):
                Left and right spatial position of the slit edges from optical model
            spat_id (`numpy.ndarray`_):
                ID assigned by PypeIt to each slit. same as in `SlitTraceSet`.


        """
        # Check that slitmask is initiated
        if self.slitmask is None:
            msgs.error('Unable to read slitmask design info')
        # as reference row we use the midpoint in the spectral direction
        reference_row = self.edge_fit[:, 0].size // 2

        # matched index for the slit-mask design data.
        ind = utils.index_of_x_eq_y(self.slitmask.slitid, maskdef_id, strict=False)
        # if not all the element of self.slitmask.slitid[ind] are equal to maskdef_id, keep only the
        # elements that are equal (matched)
        matched = np.where(self.slitmask.slitid[ind] == maskdef_id)[0]
        ind = ind[matched]
        # Number of slits
        nslits = matched.size

        # Instantiate as an empty table
        self.design = EdgeTraceSet.empty_design_table(rows=nslits)
        # Save the fit parameters and the source file as table metadata
        self.design.meta['MASKFILE'] = self.maskfile
        self.design.meta['MASKOFFL'] = cc_params_b[0]
        self.design.meta['MASKOFFR'] = cc_params_t[0]
        self.design.meta['MASKSCLL'] = cc_params_b[1]
        self.design.meta['MASKSCLR'] = cc_params_t[1]
        self.design.meta['MASKRMSL'] = cc_params_b[2]
        self.design.meta['MASKRMSR'] = cc_params_t[2]
        # Fill the columns
        self.design['TRACEID'] = np.arange(nslits, dtype=self.design['TRACEID'].dtype)
        self.design['TRACESROW'] = np.full(nslits, reference_row,
                                           dtype=self.design['TRACESROW'].dtype)
        self.design['TRACELPIX'] = self.edge_fit[reference_row,self.traceid<0][matched].astype(
                                        dtype=self.design['TRACELPIX'].dtype)
        self.design['TRACERPIX'] = self.edge_fit[reference_row,self.traceid>0][matched].astype(
                                        dtype=self.design['TRACERPIX'].dtype)
        self.design['SPAT_ID'] = spat_id[matched].astype(dtype=self.design['SPAT_ID'].dtype)
        self.design['MASKDEF_ID'] = self.slitmask.slitid[ind].astype(
                                        dtype=self.design['MASKDEF_ID'].dtype)
        self.design['SLITLMASKDEF'] = omodel_bspat[ind].astype(dtype=self.design['SLITLMASKDEF'].dtype)
        self.design['SLITRMASKDEF'] = omodel_tspat[ind].astype(dtype=self.design['SLITRMASKDEF'].dtype)
        if self.slitmask.onsky is not None:
            for i,key in enumerate(['SLITRA', 'SLITDEC', 'SLITLEN', 'SLITWID', 'SLITPA']):
                self.design[key] = self.slitmask.onsky[ind,i].astype(
                                        dtype=self.design[key].dtype)
        self.design['ALIGN'] = self.slitmask.alignment_slit[ind].astype(
                                        dtype=self.design['ALIGN'].dtype)

    def _fill_objects_table(self, maskdef_id):
        """
        Fill :attr:`objects` based on the result of the design
        registration.

        The :attr:`objects` is an `astropy.table.Table`_ with 5 columns:
            - 'OBJID': Object ID Number
            - 'OBJRA': Right ascension of the object (deg)
            - 'OBJDEC': Declination of the object (deg)
            - 'OBJNAME': Object name assigned by the observer
            - 'OBJMAG': Object magnitude provided by the observer
            - 'OBJMAG_BAND': Band of the magnitude provided by the observer
            - 'MASKDEF_ID': Slit ID Number from slit-mask design
            - 'OBJ_TOPDIST': Projected distance (in arcsec) of the object from the left
               edge of the slit (in PypeIt orientation)
            - 'OBJ_BOTDIST': Projected distance (in arcsec) of the object from the right
               edge of the slit (in PypeIt orientation)
            - 'TRACEID': Row index that matches 'TRACEID' in the design table


        Args:
            maskdef_id (`numpy.ndarray`_):
                Slit ID number from slit-mask design matched to traced slits.
        """
        # Check that slitmask is initiated
        if self.slitmask is None:
            msgs.error('Unable to read slitmask design info')

        if self.slitmask.objects is None:
            # No object data available in slit mask design object
            self.objects = None
            return

        # The index in the objects table are found by mapping the slit
        # index of each object in the design file to the slit index
        # included in the registration
        obj_index = utils.index_of_x_eq_y(self.slitmask.objects[:,0], maskdef_id, strict=False)
        # if not all the element of self.slitmask.objects[obj_index,0] are equal to maskdef_id,
        # keep only the elements that are equal (matched)
        matched = np.where(self.slitmask.objects[obj_index,0].astype(int) == maskdef_id)[0]
        obj_index = obj_index[matched]

        # Number of objects
        nobj = len(obj_index)
        # Instantiate an empty table
        self.objects = EdgeTraceSet.empty_objects_table(rows=nobj)
        # Fill the columns
        for i,key in enumerate(['MASKDEF_ID', 'OBJID', 'OBJRA', 'OBJDEC', 'OBJNAME', 'OBJMAG', 'OBJMAG_BAND',
                                'OBJ_TOPDIST', 'OBJ_BOTDIST']):
            self.objects[key] = self.slitmask.objects[obj_index,i].astype(dtype=self.objects[key].dtype)

        # SLITINDX is the index of the slit in the `design` table, not
        # in the original slit-mask design data
        self.objects['TRACEID'] = utils.index_of_x_eq_y(self.objects['MASKDEF_ID'],
                                                         self.design['MASKDEF_ID'], strict=True)

    def order_refine(self, debug=False):
        """
        For echelle spectrographs, attempt to add any orders that are not
        present in the current set of edges.
        """
        if self.spectrograph.pypeline != 'Echelle':
            msgs.warn('Parameter add_missed_orders only valid for Echelle spectrographs.')
            return
        
        if not self.can_pca():
            msgs.error('Refining the orders currently requires a PCA decomposition of the '
                       'order edges.  Ensure that the calibrations.slitedges.auto_pca parameter '
                       'is True and that there are sufficient edges to create the PCA as set by '
                       'the calibrations.slitedges.pca_min_edges parameter.  If performing a '
                       'PCA of the left and right traces independently, this minimum number '
                       'must be available for both left and right traces.')

        # Update the PCA
        self.build_pca()

        reference_row = self.left_pca.reference_row if self.par['left_right_pca'] \
                            else self.pca.reference_row
        if self.spectrograph.ech_fixed_format:
            add_left, add_right = self.order_refine_fixed_format(reference_row, debug=debug)
            rmtraces = None
        else:
            # TODO: `bracket` is hard-coded!  Currently I expect we always want
            # to set bracket=True, but we should plan to revisit this and maybe
            # expose as a user parameter.
            bracket = True
            add_left, add_right, rmtraces \
                    = self.order_refine_free_format(reference_row, bracket=bracket, debug=debug)

        if add_left is None or add_right is None:
            msgs.info('No additional orders found to add')
            return
        
        if rmtraces is not None:
            self.remove_traces(rmtraces, rebuild_pca=True)

        # Get the predicted traces
        side = np.append(np.full(add_left.size, -1, dtype=int),
                         np.full(add_right.size, 1, dtype=int))
        missed_traces = self.predict_traces(np.append(add_left, add_right), side=side)

        # Insert them
        self.insert_traces(side, missed_traces, mode='order', nudge=False)

        # If fixed-format, rematch the orders
        if self.spectrograph.ech_fixed_format:
            self.match_order(reference_row=reference_row)

    def order_refine_fixed_format(self, reference_row, debug=False):
        """
        Refine the order locations for fixed-format Echelles.
        """
        if not self.spectrograph.ech_fixed_format:
            msgs.error('order_refine_fixed_format can only be used with fixed-format Echelles!')

        # TODO:
        #   - What happens if *more* edges are detected than there are archived
        #     order positions?
        #   - When an order is added, the edges are placed at the expected
        #     position plot the measured offset.  But the trace prediction
        #     requires the spatial positoion at the relevant reference row.
        #     This needs to be checked.

        # First match the expected orders
        spat_offset = self.match_order(reference_row=reference_row)

        available_orders = self.orderid[1::2]
        missed_orders = np.setdiff1d(self.spectrograph.orders, available_orders)
        if missed_orders.size == 0:
            # No missing orders, we're done
            return None, None

        # Find the indices of the missing orders
        missed_orders_indx = utils.index_of_x_eq_y(self.spectrograph.orders, missed_orders)

        # Get the spatial positions of the new left and right order edges
        add_right_edges = (self.spectrograph.order_spat_pos[missed_orders_indx]
                            + self.spectrograph.order_spat_width[missed_orders_indx]/2.
                            + spat_offset) * self.nspat

        add_left_edges = (self.spectrograph.order_spat_pos[missed_orders_indx]
                            - self.spectrograph.order_spat_width[missed_orders_indx]/2.
                            + spat_offset) * self.nspat
        
        return add_left_edges, add_right_edges

    # NOTE: combined_order_tol is effectively hard-coded; i.e., the current code
    # always uses the default when calling this function.
    def order_refine_free_format(self, reference_row, combined_order_tol=1.8, bracket=True,
                                 debug=False):
        """
        Refine the order locations for "free-format" Echelles.

        Traces must be synced before calling this function.

        The procedure is as follows:

            - The function selects the good traces and calculates the width of
              each order and the gap between each order and fits them with
              Legendre polynomials (using the polynomial orders set by the
              ``order_width_poly`` and ``order_gap_poly`` parameters); 5-sigma
              outliers are removed from the fit.

            - Based on this fit, the code adds missed orders, both interspersed
              with detected orders and extrapolated over the full spatial range
              of the detector/mosaic.  The spatial extent over which this
              prediction is performed is set by ``order_spat_range`` and can be
              limited by any resulting overlap in the prediction, as set by
              ``max_overlap``.

            - Any detected "orders" that are actually the adjoining of one or
              more orders are flagged for rejection.

        Args:
            reference_row (:obj:`int`):
                The index of the spectral pixel (row) in the set of left and
                right traces at which to predict the positions of the missed
                orders.  Nominally, this is the reference row used for the
                construction of the trace PCA.
            combined_order_tol (:obj:`float`, optional):
                For orders that are very nearly overlapping, the automated edge
                tracing can often miss the right and left edges of two adjacent
                orders.  This leads to the detected edges of two adjacent orders
                being combined into a single order.  This value sets the maximum
                ratio of the width of any given detected order to the polynomial
                fit to the order width as a function of spatial position on the
                detector.
            bracket (:obj:`bool`, optional):
                Bracket the added orders with one additional order on either side.
                This can be useful for dealing with predicted overlap.
            debug (:obj:`bool`, optional):
                Run in debug mode.

        Returns:
            :obj:`tuple`: Three `numpy.ndarray`_ objects that provide (1,2) the
            left and right edges of orders to be added to the set of edge traces
            and (3) a boolean array indicating which of the existing traces
            should be removed.
        """
        # Select the left/right traces
        # TODO: This is pulled from get_slits.  Maybe want a function for this.
        gpm = self.synced_selection(self.good_traces(), mode='neither')
        # Save the list of good left edges in case we need to remove any
        left_gpm = gpm & self.is_left
        left = self.edge_fit[:,left_gpm]
        right_gpm = gpm & self.is_right
        right = self.edge_fit[:,right_gpm]

        # Use the trace locations at the middle of the spectral shape of the
        # detector/mosaic
        left = left[reference_row]
        right = right[reference_row]

        # Get the order centers, widths, and gaps
        cen = (right + left)/2
        width = right - left
        gap = left[1:] - right[:-1]

        # Create the polynomial models.
        width_fit = fitting.robust_fit(cen, width, self.par['order_width_poly'],
                                       function='legendre', lower=self.par['order_fitrej'],
                                       upper=self.par['order_fitrej'], maxiter=5, sticky=True)
        # Connection of center to gap uses the gap spatially *after* the order.
        gap_fit = fitting.robust_fit(cen[:-1], gap, self.par['order_gap_poly'],
                                     function='legendre', lower=self.par['order_fitrej'],
                                     upper=self.par['order_fitrej'], maxiter=5, sticky=True)
        
        # Ideally, measured widths/gaps should be rejected for one of the
        # following reasons:
        #   - The width is too large because gaps were missed (i.e. multiple
        #     orders were combined)
        #   - The width is too small because order overlap was detected and
        #     removed.
        #   - The gap is too large because orders were missed

        # In the case when the user does not reject "outliers", we still reject
        # orders that we expected to be cases where multiple orders have been
        # combined
        bad_order = width / width_fit.eval(cen) > combined_order_tol
        if self.par['order_outlier'] is not None:
            # Exclude "outliers"
            resid = np.absolute(width_fit.yval - width_fit.eval(width_fit.xval))
            bad_order |= (resid/width_fit.calc_fit_rms() > self.par['order_outlier'])
            # TODO: The gaps for HIRES can have *very* large residuals.  Using
            # the gaps to identify outliers would remove many orders that
            # probably shouldn't be removed.
#            resid = np.absolute(gap_fit.yval - gap_fit.eval(gap_fit.xval))
#            bad_order[:-1] |= (resid/gap_fit.calc_fit_rms() > self.par['order_outlier'])

        # And sets flags used to remove them, in favor of replacing them with
        # the predicted locations of the individual orders.
        rmtraces = np.zeros(left_gpm.size, dtype=bool)
        rmtraces[np.where(left_gpm)[0][bad_order]] = True
        rmtraces = self.synced_selection(rmtraces, mode='both')
        
        # Interpolate any missing orders
        # TODO: Expose tolerances to the user?
        good_order = np.logical_not(bad_order)
        order_cen, order_missing \
                = trace.find_missing_orders(cen[good_order], width_fit, gap_fit)
        if np.sum(order_missing) > order_missing.size // 2:
            msgs.warn('Found more missing orders than detected orders.  Check the order '
                      'refinement QA file!  The code will continue, but you likely need to adjust '
                      'your edge-tracing parameters.')

        # Extrapolate orders; this includes one additional order to either side
        # of the spatial extent set by rng.
        rng = [0., float(self.nspat)] if self.par['order_spat_range'] is None \
                    else self.par['order_spat_range']
        lower_order_cen, upper_order_cen \
                    = trace.extrapolate_orders(cen[good_order], width_fit, gap_fit,
                                               rng[0], rng[1], bracket=bracket)

        # Combine the results
        order_cen = np.concatenate((lower_order_cen, order_cen, upper_order_cen))
        order_missing = np.concatenate((np.ones(lower_order_cen.size, dtype=bool),
                                        order_missing, np.ones(upper_order_cen.size, dtype=bool)))

        # If nothing is missing, return
        if not np.any(order_missing):
            return None, None, None

        # QA Plot
        ofile = None if debug else self.qa_path / 'PNGs' \
                    / f'{self.get_path().name.split(".")[0]}_orders_qa.png'
        # TODO: Making this directory should probably be done elsewhere
        if ofile is not None and not ofile.parent.is_dir():
            ofile.parent.mkdir(parents=True)
        self.order_refine_free_format_qa(cen, bad_order, width, gap, width_fit, gap_fit,
                                         order_cen, order_missing, bracket=bracket, ofile=ofile)

        # Return the coordinates for the left and right edges to add
        add_width = width_fit.eval(order_cen[order_missing])
        add_left = order_cen[order_missing] - add_width / 2
        add_right = order_cen[order_missing] + add_width / 2

        # Join the added edges with the existing ones
        _left = np.append(add_left, left[good_order])
        # Create a sorting vector
        srt = np.argsort(_left)
        # Create a vector that will reverse the sorting
        isrt = np.argsort(srt)
        # Join and sort the right edges
        _right = np.append(add_right, right[good_order])[srt]
        # Sort the left edges
        _left = _left[srt]

        # NOTE: Although I haven't tested this, I think this approach works best
        # under the assumption that the overlap *decreases* from small pixel
        # numbers to large pixel numbers.  This should be true if the pypeit
        # convention is maintained with blue orders toward small pixel values
        # and red orders at large pixel values.

        # Deal with overlapping orders among the ones to be added.  The edges
        # are adjusted equally on both sides to avoid changing the order center
        # and exclude the overlap regions from the reduction.
        if np.all(_left[1:] - _right[:-1] > 0):
            if bracket:
                add_left, add_right = self._handle_bracketing_orders(add_left, add_right)
            # There is no overlap, so just return the orders to add
            return add_left, add_right, rmtraces

        # Used to remove orders that have too much overlap
        nord = _left.size
        ok_overlap = np.ones(nord, dtype=bool)

        # Loop sequentially so that each pair is updated as the loop progresses
        for i in range(1, nord):
            # *Negative* of the gap; i.e., positives values means there's
            # overlap
            ngap = _right[i-1] - _left[i]
            if ngap > 0:
                if self.par['max_overlap'] is not None:
                    ok_overlap[i-1] = 2*ngap/(_right[i-1] - _left[i-1]) < self.par['max_overlap']
                    ok_overlap[i] = 2*ngap/(_right[i] - _left[i]) < self.par['max_overlap']
                # Adjust both order edges to avoid the overlap region but
                # keep the same center coordinate
                _left[i-1] += ngap
                _right[i-1] -= ngap
                _left[i] += ngap
                _right[i] -= ngap

        # For any *existing* traces that were adjusted because of the overlap,
        # this applies the adjustment to the `edge_fit` data.
        # NOTE: This only adjusts the "fit" locations (edge_fit), *not* the
        # measured centroid locations (edge_cen).  This should not cause
        # problems because, e.g., the `get_slits` function uses `edge_fit`.
        nadd = add_left.size
        left_indx = np.where(left_gpm)[0][good_order]
        offset = _left[isrt][nadd:] - left
        self.edge_fit[:,left_indx] += offset[None,:]
        right_indx = np.where(right_gpm)[0][good_order]
        offset = _right[isrt][nadd:] - right
        self.edge_fit[:,right_indx] += offset[None,:]

        # Get the adjusted traces to add.  Note this currently does *not* change
        # the original traces
        ok_overlap = ok_overlap[isrt][:nadd]
        add_left = _left[isrt][:nadd][ok_overlap]
        add_right = _right[isrt][:nadd][ok_overlap]

        if bracket:
            add_left, add_right = self._handle_bracketing_orders(add_left, add_right)
        return add_left, add_right, rmtraces

    @staticmethod
    def _handle_bracketing_orders(add_left, add_right):
        """
        Utility function to remove added orders that bracket the left and right
        edge of the detector, used to handle overlap.

        Args:
            add_left (`numpy.ndarray`_):
                List of left edges to add
            add_right (`numpy.ndarray`_):
                List of right edges to add

        Returns:
            :obj:`tuple`: The two `numpy.ndarray`_ objects after removing the
            bracketing orders.
        """
        nadd = add_left.size
        if nadd < 2:
            # TODO: The code should not get here!  If it does, we need to
            # figure out why and fix it.
            msgs.error('CODING ERROR: Order bracketing failed!')
        if nadd == 2:
            return None, None
        return add_left[1:-1], add_right[1:-1]
    
    def order_refine_free_format_qa(self, cen, bad_order, width, gap, width_fit, gap_fit,
                                    order_cen, order_missing, bracket=False, ofile=None):
        """
        Create the QA plot for order modeling.

        Args:
            cen (`numpy.ndarray`_):
                Spatial centers of the detected orders.
            bad_order (`numpy.ndarray`_):
                Boolean array selecting "orders" that have been flagged as
                outliers.
            width (`numpy.ndarray`_):
                Measured order spatial widths in pixels.
            gap (`numpy.ndarray`_):
                Measured order gaps in pixels.
            width_fit (:class:`~pypeit.core.fitting.PypeItFit`):
                Model of the order width as a function of the order center.
            gap_fit (:class:`~pypeit.core.fitting.PypeItFit`):
                Model of the order gap *after* each order as a function of the order
                center.
            order_cen (`numpy.ndarray`_):
                Spatial centers of all "individual" orders.
            order_missing (`numpy.ndarray`_):
                Boolean array selecting "individual" orders that were not traced
                by the automated tracing and flagged as missing.  See
                :func:`~pypeit.core.trace.find_missing_orders` and
                :func:`~pypeit.core.trace.extrapolate_orders`.
            bracket (:obj:`bool`, optional):
                Flag that missing orders have been bracketed by additional
                orders in an attempt to deal with overlap regions.
            ofile (:obj:`str`, `Path`_, optional):
                Path for the QA figure file.  If None, the plot is shown in a
                matplotlib window.
        """
        # Setup
        w_resid = width - width_fit.eval(cen)
        w_rms = width_fit.calc_fit_rms()
        med_wr = np.median(w_resid)
        mad_wr = np.median(np.absolute(w_resid - med_wr))

        w_out = bad_order & width_fit.gpm.astype(bool)
        w_rej = np.logical_not(bad_order) & np.logical_not(width_fit.gpm)
        w_outrej = bad_order & np.logical_not(width_fit.gpm)
        w_good = np.logical_not(w_out | w_rej | w_outrej)

        g_cen = cen[:-1]
        g_bad_order = bad_order[:-1]
        g_resid = gap - gap_fit.eval(g_cen)
        g_rms = gap_fit.calc_fit_rms()
        med_gr = np.median(g_resid)
        mad_gr = np.median(np.absolute(g_resid - med_gr))

        g_out = g_bad_order & gap_fit.gpm.astype(bool)
        g_rej = np.logical_not(g_bad_order) & np.logical_not(gap_fit.gpm)
        g_outrej = g_bad_order & np.logical_not(gap_fit.gpm)
        g_good = np.logical_not(g_out | g_rej | g_outrej)

        # Set the spatial limits based on the extent of the order centers and/or
        # the detector spatial extent
        sx = min(0, np.amin(order_cen))
        ex = max(self.nspat, np.amax(order_cen))
        buf = 1.1
        xlim = [(sx * (1 + buf) + ex * (1 - buf))/2, (sx * (1 - buf) + ex * (1 + buf))/2]

        # Set the residual plot limits based on the median and median absolute
        # deviation
        width_lim = np.array([med_wr - 20*mad_wr, med_wr + 20*mad_wr])
        gap_lim = np.array([med_gr - 20*mad_gr, med_gr + 20*mad_gr])

        # Sample the width and gap models over the full width of the plots
        mod_cen = np.linspace(*xlim, 100)
        mod_width = width_fit.eval(mod_cen)
        mod_gap = gap_fit.eval(mod_cen)

        # Create the plot
        w,h = plt.figaspect(1)
        fig = plt.figure(figsize=(1.5*w,1.5*h))

        # Plot the data and each fit
        ax = fig.add_axes([0.10, 0.35, 0.8, 0.6])
        ax.minorticks_on()
        ax.tick_params(which='both', direction='in', top=True, right=True)
        ax.grid(True, which='major', color='0.7', zorder=0, linestyle='-')
        ax.set_xlim(xlim)
        ax.xaxis.set_major_formatter(ticker.NullFormatter())

        # Set the plot title
        title = 'Order prediction model'
        if bracket:
            title += ' (bracketed)'
        ax.text(0.5, 1.02, title, ha='center', va='center', transform=ax.transAxes)

        # Plot the detector bounds
        ax.axvline(0, color='k', ls='--', lw=2)
        ax.axvline(self.nspat, color='k', ls='--', lw=2)

        # Models
        ax.plot(mod_cen, mod_width, color='C0', alpha=0.3, lw=3, zorder=3)
        ax.plot(mod_cen, mod_gap, color='C2', alpha=0.3, lw=3, zorder=3)

        # Measurements included in the fit
        ax.scatter(cen[w_good], width[w_good],
                   marker='.', color='C0', s=50, lw=0, label='fitted widths', zorder=4)
        if np.any(w_rej):
            # Rejected but not considered an outlier
            ax.scatter(cen[w_rej], width[w_rej],
                       marker='x', color='C1', s=50, lw=1, label='rej widths', zorder=4)
        if np.any(w_out):
            # Outlier but not rejected
            ax.scatter(cen[w_out], width[w_out],
                       marker='^', facecolor='none', edgecolor='C1', s=50, lw=1,
                       label='outlier widths', zorder=4)
        if np.any(w_outrej):
            # Both outlier and rejected
            ax.scatter(cen[w_outrej], width[w_outrej],
                       marker='^', facecolor='C1', s=50, lw=1, label='rej,outlier widths', zorder=4)
        # Orders to add
        ax.scatter(order_cen[order_missing], width_fit.eval(order_cen[order_missing]),
                   marker='s', facecolor='none', edgecolor='C0', s=80, lw=1,
                   label='missing widths', zorder=3)

        # Same as above but for gaps
        ax.scatter(g_cen[g_good], gap[g_good],
                   marker='.', color='C2', s=50, lw=0, label='fitted gaps', zorder=4)
        if np.any(g_rej):
            ax.scatter(g_cen[g_rej], gap[g_rej],
                       marker='x', color='C4', s=50, lw=1, label='rej gaps', zorder=4)
        if np.any(g_out):
            ax.scatter(g_cen[g_out], gap[g_out],
                       marker='^', facecolor='none', edgecolor='C4', s=50, lw=1,
                       label='outlier gaps', zorder=4)
        if np.any(g_outrej):
            ax.scatter(g_cen[g_outrej], gap[g_outrej],
                       marker='^', facecolor='C4', s=50, lw=1, label='rej,outlier gaps', zorder=4)
        ax.scatter(order_cen[order_missing], gap_fit.eval(order_cen[order_missing]),
                   marker='s', facecolor='none', edgecolor='C2', s=80, lw=1,
                   label='missing gaps', zorder=3)

        # Add the y label and legend
        ax.set_ylabel('Order Width/Gap [pix]')
        ax.legend()

        # Plot the width residuals
        ax = fig.add_axes([0.10, 0.25, 0.8, 0.1])
        ax.minorticks_on()
        ax.tick_params(which='both', direction='in', top=True, right=False)
        ax.set_xlim(xlim)
        ax.set_ylim(width_lim)
        ax.xaxis.set_major_formatter(ticker.NullFormatter())

        # Plot the detector bounds
        ax.axvline(0, color='k', ls='--', lw=2)
        ax.axvline(self.nspat, color='k', ls='--', lw=2)

        # Model is at 0 residual
        ax.axhline(0, color='C0', alpha=0.3, lw=3, zorder=2)
        # Measurements included in the fit
        ax.scatter(cen[w_good], w_resid[w_good], marker='.', color='C0', s=50, lw=0, zorder=4)
        # Rejected but not considered an outlier
        ax.scatter(cen[w_rej], w_resid[w_rej], marker='x', color='C1', s=50, lw=1, zorder=4)
        # Outlier but not rejected
        ax.scatter(cen[w_out], w_resid[w_out],
                   marker='^', facecolor='none', edgecolor='C1', s=50, lw=1, zorder=4)
        # Both outlier and rejected
        ax.scatter(cen[w_outrej], w_resid[w_outrej],
                   marker='^', facecolor='C1', s=50, lw=1, zorder=4)

        # Add the label
        ax.set_ylabel(r'$\Delta$Width')

        # Add a right axis that gives the residuals normalized by the rms; use
        # this to set the grid.
        axt = ax.twinx()
        axt.minorticks_on()
        axt.tick_params(which='both', direction='in')
        axt.grid(True, which='major', color='0.7', zorder=0, linestyle='-')
        axt.set_xlim(xlim)
        axt.set_ylim(width_lim / w_rms)
        axt.set_ylabel(r'$\Delta$/RMS')

        # Plot the gap residuals
        ax = fig.add_axes([0.10, 0.15, 0.8, 0.1])
        ax.minorticks_on()
        ax.tick_params(which='both', direction='in', top=True, right=False)
        ax.set_xlim(xlim)
        ax.set_ylim(gap_lim)

        # Plot the detector bounds
        ax.axvline(0, color='k', ls='--', lw=2)
        ax.axvline(self.nspat, color='k', ls='--', lw=2)

        # Model is at 0 residual
        ax.axhline(0, color='C2', alpha=0.3, lw=3, zorder=2)
        # Measurements included in the fit
        ax.scatter(g_cen[g_good], g_resid[g_good],
                   marker='.', color='C2', s=50, lw=0, zorder=4)
        # Rejected but not considered an outlier
        ax.scatter(g_cen[g_rej], g_resid[g_rej],
                   marker='x', color='C4', s=50, lw=1, zorder=4)
        # Outlier but not rejected
        ax.scatter(g_cen[g_out], g_resid[g_out],
                   marker='^', facecolor='none', edgecolor='C4', s=50, lw=1, zorder=4)
        # Both outlier and rejected
        ax.scatter(g_cen[g_outrej], g_resid[g_outrej],
                   marker='^', facecolor='C4', s=50, lw=1, zorder=4)
        
        # Add the axis labels
        ax.set_ylabel(r'$\Delta$Gap')
        ax.set_xlabel('Spatial pixel')

        # Add a right axis that gives the residuals normalized by the rms; use
        # this to set the grid.
        axt = ax.twinx()
        axt.minorticks_on()
        axt.tick_params(which='both', direction='in')
        axt.grid(True, which='major', color='0.7', zorder=0, linestyle='-')
        axt.set_xlim(xlim)
        axt.set_ylim(gap_lim / g_rms)
        axt.set_ylabel(r'$\Delta$/RMS')

        if ofile is None:
            plt.show()
        else:
            fig.canvas.print_figure(ofile, bbox_inches='tight')
            msgs.info(f'Missing order QA written to: {ofile}')
        fig.clear()
        plt.close(fig)

    def slit_spatial_center(self, normalized=True, spec=None, use_center=False, 
                            include_box=False):
        """
        Return the spatial coordinate of the center of each slit.

        The slit edges must be left-right synchronized.

        Args:
            normalized (:obj:`bool`, optional):
                Return coordinates normalized by the size of the
                detector.
            spec (:obj:`int`, optional):
                Spectral position (row) at which to return the spatial position.
                If ``None``, use the PCA reference row if a PCA exists or the
                central row (i.e., ``self.nspec//2``), otherwise.
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

        # Set the spectral position to use as a reference.
        _spec = spec
        if _spec is None:
            if self.pcatype is None:
                _spec = self.nspec//2
            else:
                _spec = self.left_pca.reference_row if self.par['left_right_pca'] \
                            else self.pca.reference_row

        # Synced, spatially sorted traces are always ordered in left,
        # right pairs
        trace_cen = self.edge_cen[:,gpm] if self.edge_fit is None or use_center \
                        else self.edge_fit[:,gpm]
        slit_c = np.ma.MaskedArray(np.zeros(self.nslits, dtype=float))
        slit_c[np.logical_not(good_slit)] = np.ma.masked
        slit_c[good_slit] = np.mean(trace_cen[_spec,:].reshape(-1,2), axis=1)
        return slit_c/self.nspat if normalized else slit_c

    def match_order(self, reference_row=None):
        """
        Match synchronized slits to the expected echelle orders.

        This function will fault if called for non-Echelle spectrographs!

        For Echelle spectrographs, this finds the best matching order for each
        left-right trace pair; the traces must have already been synchronized
        into left-right pairs. Currently, this is a very simple, non-optimized
        match:

            - The closest order from ``self.spectrograph.order_spat_pos`` is
              matched to each slit.
            - Any slit that cannot be matched to an order -- either because
              there are more "slits" than orders or because the separation is
              larger than the provided tolerance -- is flagged as
              ``ORDERMISMATCH``.
            - A warning is issued if the number of valid matches is not
              identical to the number of expected orders
              (``self.spectrograph.norders``). The warning includes the list of
              orders that were not identified.

        The match tolerance is set by the parameter ``order_match``.  An offset
        can be applied to improve the match via the parameter ``order_offset``;
        i.e., this should minimize the difference between the expected order
        positions and ``self.slit_spatial_center() + self.par['order_offset']``.
        Both ``order_match`` and ``order_offset`` are given in fractions of the
        detector size along the spatial axis.

        The result of this method is to instantiate :attr:`orderid`.

        Args:
            reference_row (:obj:`int`, optional):
                The spectral pixel (row) used to generate spatial positions of
                the orders to match against the expected positions.  If
                ``None``, use the PCA reference row if a PCA exists or the
                central row (i.e., ``self.nspec//2``), otherwise.

        Returns:
            :obj:`float`: The median offset in the relative of the detector size
            between the archived order positions and those measured via the edge
            tracing.

        Raises:
            PypeItError:
                Raised if the number of orders, the order number, or their
                spatial locations are not defined for an Echelle spectrograph.
                Also raised if the edge traces are not synced.
        """

        if self.spectrograph.norders is None:
            msgs.error('Coding error: norders not defined for {0}!'.format(
                        self.spectrograph.__class__.__name__))
        if self.spectrograph.orders is None:
            msgs.error('Coding error: orders not defined for {0}!'.format(
                        self.spectrograph.__class__.__name__))
        if self.spectrograph.order_spat_pos is None:
            msgs.error('Coding error: order_spat_pos not defined for {0}!'.format(
                       self.spectrograph.__class__.__name__))
        if not self.is_synced:
            msgs.error('EdgeTraceSet must be synced to match to orders.')

        offset = self.par['order_offset']
        if offset is None:
            offset = 0.0

        # Get the order centers in fractions of the detector width.  This
        # requires the slits to be synced (checked above)!  Masked elements in
        # slit_cen are for bad slits or syncing.
        slit_cen = self.slit_spatial_center(spec=reference_row) + offset
        good_sync = np.logical_not(np.ma.getmaskarray(slit_cen))
        # "slit_indx" matches the "slit" index to the order number.  I.e.,
        # "slit_indx" has one element per expected order, and the value at a
        # given position is the index of the paired traces for the relevant
        # order.
        slit_indx = slitdesign_matching.match_positions_1D(
                slit_cen.data[good_sync],           # (Good) Measured positions
                self.spectrograph.order_spat_pos,   # Expected positions
                tol=self.par['order_match'])        # Matching tolerance

        # Boolean array selecting found orders
        fnd = slit_indx > -1
        missed_orders = self.spectrograph.orders[np.logical_not(fnd)]
        if not np.all(fnd):
            msgs.warn(f'Did not find all orders!  Missing orders: {missed_orders}')

        # Flag paired edges that were not matched to a known order
        nomatch = np.setdiff1d(np.arange(np.sum(good_sync)), slit_indx[fnd])
        if nomatch.size > 0:
            msgs.warn(f'Flagging {nomatch.size} trace pairs as not being matched to an order.')
            # Create a vector that selects the appropriate traces.  This
            # *assumes* that the traces are left-right syncronized and the order
            # has not changed between the order of the traces in the relevant
            # array and how the centers of the synced traces are computed
            # (slit_spatial_center).
            flag = np.append(2*nomatch, 2*nomatch+1)
            self.edge_msk[:,flag] = self.bitmask.turn_on(self.edge_msk[:,flag], 'ORDERMISMATCH')

        # Minimum separation between the order and its matching slit;
        # keep the signed value for reporting, but use the absolute
        # value of the difference for vetting below.
        # NOTE: This includes indices for orders that were not found.  This is
        # largely for book-keeping purposes in the print statement below.
        sep = self.spectrograph.order_spat_pos - slit_cen.data[good_sync][slit_indx]
        med_offset = np.median(sep[fnd])

        # Report
        msgs.info(f'Median offset is {med_offset:.3f}.')
        msgs.info('After offsetting, order-matching separations are:')
        msgs.info(f' {"ORDER":>6} {"PAIR":>4} {"SEP":>6}')
        msgs.info(f' {"-"*6} {"-"*4} {"-"*6}')
        for i in range(self.spectrograph.norders):
            if fnd[i]:
                msgs.info(f' {self.spectrograph.orders[i]:>6} {i+1:>4} {sep[i]-med_offset:6.3f}')
            else:
                msgs.info(f' {self.spectrograph.orders[i]:>6} {"N/A":>4} {"MISSED":>6}')
        msgs.info(f' {"-"*6} {"-"*4} {"-"*6}')

        # Instantiate the order ID; 0 means the order is unassigned
        self.orderid = np.zeros(self.nslits*2, dtype=int)
        raw_indx = np.arange(self.nslits)[good_sync][slit_indx[fnd]]
        indx = (2*raw_indx[:,None] + np.tile(np.array([0,1]), (np.sum(fnd),1))).ravel()
        self.orderid[indx] = (np.array([-1,1])[None,:]*self.spectrograph.orders[fnd,None]).ravel()

        return med_offset

    def get_slits(self):
        """
        Use the data to instatiate the relevant
        :class:`~pypeit.slittrace.SlitTraceSet` object.

        The traced edges must have first been organized into slits;
        see :func:`sync`.

        For fixed-format echelle spectrographs, the method automatically calls
        :func:`match_order` and will only include the "slits" that have been
        correctly matched to known orders.
        
        The :class:`~pypeit.slittrace.SlitTraceSet` object will use calibration
        naming keys used by this parent :class:`EdgeTraceSet` object.

        The :class:`~pypeit.slittrace.SlitTraceSet` object will have
        an `astropy.table.Table`_ resulted from merging :attr:`design`
        and :attr:`objects` together.

        Returns:
            :class:`~pypeit.slittrace.SlitTraceSet`: Object holding
            the slit traces.
        """
        if not self.is_synced:
            msgs.error('Edges must be synced to construct SlitTraceSet object.')

        # For echelle spectrographs, match the left-right trace pairs
        # to echelle orders
        if self.spectrograph.pypeline == 'Echelle' and self.spectrograph.ech_fixed_format:
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
        if self.par['trim_spec'] is not None:
            trim_low, trim_high = self.par['trim_spec']
            specmin = np.asarray([trim_low]*nslit,dtype=np.float64)
            specmax = np.asarray([self.nspec-trim_high]*nslit,dtype=np.float64)
        elif self.spectrograph.spec_min_max is None or ech_order is None:
            specmin = np.asarray([-np.inf]*nslit)
            specmax = np.asarray([np.inf]*nslit)
        else:
            specmin, specmax = self.spectrograph.spec_min_max
            indx = utils.index_of_x_eq_y(self.spectrograph.orders, ech_order)
            specmin = specmin[indx]/binspec
            specmax = specmax[indx]/binspec

        # Define spat_id (in the same way is done in SlitTraceSet) to add it in the astropy table. It
        # is useful to have spat_id in that table.
        center = (left + right) / 2
        spat_id = np.round(center[int(np.round(0.5 * self.nspec)), :]).astype(int)

        if self.maskdef_id is not None:
            # Check for mismatched `maskdef_id` in the left and right edges
            mkd_id_mismatch = self.maskdef_id[self.is_left] != self.maskdef_id[self.is_right]
            if np.any(mkd_id_mismatch):
                msgs.warn("Mismatched `maskdefId` in left and right traces for {}/{} slits. ".format(
                          self.maskdef_id[self.is_left][mkd_id_mismatch].size, self.nslits) +
                          "Choosing the left edge `maskdefId` if it is not -99, otherwise choosing right one")
            _maskdef_id = self.maskdef_id[gpm & self.is_left]
            # choose the value of `maskdef_id` that is not -99.
            mkd_id_bad = _maskdef_id == -99
            # this may not work if the corresponding right edge is also -99. Assuming this is not the case
            _maskdef_id[mkd_id_bad] = self.maskdef_id[gpm & self.is_right][mkd_id_bad]
            if np.any(_maskdef_id == -99):
                msgs.warn("{} slits do not have `maskdefId` assigned.".format(_maskdef_id[_maskdef_id == -99].size) +
                          "They will not be included in the design table")

            # Store the matched slit-design and object information in a table.
            self._fill_design_table(_maskdef_id, self.cc_params_b, self.cc_params_b, self.omodel_bspat,
                                    self.omodel_tspat, spat_id)
            self._fill_objects_table(_maskdef_id)
            # TODO - instead of merging the two tables, just create a single one
            _merged_designtab = table.join(self.design, self.objects, keys=['TRACEID'])
            _merged_designtab.remove_column('MASKDEF_ID_2')
            _merged_designtab.rename_column('MASKDEF_ID_1', 'MASKDEF_ID')
            # One more item
            _posx_pa = float(self.slitmask.posx_pa)
        else:
            _maskdef_id = None
            _merged_designtab = None
            _posx_pa = None

        # Instantiate and return
        slits = slittrace.SlitTraceSet(left, right, self.spectrograph.pypeline,
                                       detname=self.traceimg.detector.name, nspat=self.nspat,
                                       PYP_SPEC=self.spectrograph.name, specmin=specmin,
                                       specmax=specmax, binspec=binspec, binspat=binspat,
                                       pad=self.par['pad'], mask_init=slit_msk,
                                       maskdef_id=_maskdef_id, maskdef_designtab=_merged_designtab,
                                       maskdef_posx_pa=_posx_pa, maskfile=self.maskfile,
                                       ech_order=ech_order)
        # Copy the CalibFrame internals (also identical to the self.traceimg internals)
        slits.copy_calib_internals(self)

        return slits


