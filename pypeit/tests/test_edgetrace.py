"""
Tests of the slit-edge tracing (EdgeTraceSet).
"""
import numpy as np
import pytest

from pypeit.core import flexure
from pypeit.core import trace
from pypeit.edgetrace import EdgeTraceSet
from pypeit.images.buildimage import TraceImage
from pypeit.images.detector_container import DetectorContainer
from pypeit.spectrographs.util import load_spectrograph
from pypeit.tests import tstutils


def _slicer_trace_image(nslit=6, width=60., gap=4., sig=0.8, nspec=300, level=2e4, shift=0.,
                        seed=3):
    """
    Synthetic trace image of closely packed, slightly curved slits (as for an
    image slicer).  The true edges are at the steepest flux gradient (the 50%
    points of error-function edges).  The image of the slits can be shifted
    spatially by ``shift`` pixels (e.g. by spatial flexure); the returned edges
    are not shifted.
    """
    spec = np.arange(nspec)
    curve = 1.5 * ((spec - nspec / 2) / nspec) ** 2
    lefts = np.array([30.37 + i * (width + gap) + curve for i in range(nslit)]).T
    rights = lefts + width
    nspat = int(rights.max()) + 30
    xx = np.arange(nspat)[None, :]
    img = np.zeros((nspec, nspat))
    for i in range(nslit):
        img += level * (trace._edge_model(xx - lefts[:, i:i + 1] - shift, 0., sig, 1)
                        + trace._edge_model(xx - rights[:, i:i + 1] - shift, 0., sig, -1) - 1)
    img = np.clip(img, 0, None) + 50.
    rng = np.random.default_rng(seed)
    return rng.normal(img, np.sqrt(img)), lefts, rights


def test_refine_edges():
    # With refine_edges (the default for SlicerIFU spectrographs), the slit
    # edges are at the steepest flux gradient; the traced edges are not
    img, lefts, rights = _slicer_trace_image()
    spectrograph = load_spectrograph('keck_kcwi')
    par = spectrograph.default_pypeit_par()['calibrations']['slitedges']
    assert par['refine_edges'], 'refine_edges should be True for SlicerIFU spectrographs'
    offsets = {}
    for refine in [False, True]:
        par['refine_edges'] = refine
        traceimg = TraceImage(img, ivar=1 / img)
        traceimg.detector = DetectorContainer(**tstutils.default_detector())
        edges = EdgeTraceSet(traceimg, spectrograph, par, auto=True)
        left, right = edges.edge_fit[:, edges.is_left], edges.edge_fit[:, edges.is_right]
        assert left.shape == lefts.shape and right.shape == rights.shape, \
            'Incorrect number of slits'
        offsets[refine] = (np.median(left - lefts, axis=0), np.median(right - rights, axis=0))
    assert np.all(offsets[False][0] < -0.3) and np.all(offsets[False][1] > 0.3), \
        'Test is not sensitive: the traced edges should be outside the steepest gradient'
    assert np.all(np.absolute(offsets[True][0]) < 0.05) \
           and np.all(np.absolute(offsets[True][1]) < 0.05), \
        'The refined edges should be at the steepest flux gradient'
    # If the edges cannot be refined (here, because the slits are not visible),
    # the traced edges are kept
    traced = edges.edge_fit.copy()
    edges.traceimg.image = np.random.default_rng(4).normal(50., 7., img.shape)
    edges.edge_refine()
    assert np.array_equal(edges.edge_fit, traced), 'Edges that cannot be refined should be kept'


def test_refine_edges_spat_flexure():
    # The edge-based spatial flexure is measured relative to the steepest flux
    # gradient, so it is only zero for the trace frame itself if the traced
    # edges are refined
    img, _, _ = _slicer_trace_image()
    spectrograph = load_spectrograph('keck_kcwi')
    par = spectrograph.default_pypeit_par()['calibrations']['slitedges']
    for refine in [False, True]:
        par['refine_edges'] = refine
        traceimg = TraceImage(img, ivar=1 / img)
        traceimg.detector = DetectorContainer(**tstutils.default_detector())
        slits = EdgeTraceSet(traceimg, spectrograph, par, auto=True).get_slits()
        spat_flexure = flexure.spat_flexure_shift(img, slits, method='edge', maxlag=10,
                                                  sigdetect=10.)
        if not refine:
            assert np.all(spat_flexure[:, 0] > 0.3) and np.all(spat_flexure[:, 1] < -0.3), \
                'Test is not sensitive: the traced edges should be outside the steepest gradient'
            continue
        assert np.all(np.absolute(spat_flexure) < 0.02), 'The trace frame should have no flexure'
        # A frame shifted by 1.6 pixels
        shifted, _, _ = _slicer_trace_image(shift=1.6, seed=5)
        spat_flexure = flexure.spat_flexure_shift(shifted, slits, method='edge', maxlag=10,
                                                  sigdetect=10.)
        assert np.allclose(spat_flexure, 1.6, atol=0.02), 'Incorrect spatial flexure'


def _longslit_image(edges, sig=1.0, level=2e4, shift=0., seed=1, nspec=400, nspat=400):
    """
    Synthetic image of slits whose edges can be beyond the detector (e.g. a
    long slit cut by the detector boundary).  ``edges`` is a list of (left,
    right) positions, with None for an edge beyond the detector.  The edges
    are slightly curved, and the image of the slits can be shifted spatially
    by ``shift`` pixels.
    """
    spec = np.arange(nspec)
    curve = 2.0 * ((spec - nspec / 2) / nspec) ** 2
    xx = np.arange(nspat)[None, :]
    img = np.zeros((nspec, nspat))
    for left, right in edges:
        slit = np.ones((nspec, nspat))
        if left is not None:
            slit *= trace._edge_model(xx - left - curve[:, None] - shift, 0., sig, 1)
        if right is not None:
            slit *= trace._edge_model(xx - right - curve[:, None] - shift, 0., sig, -1)
        img += slit
    img = level * img + 50.
    rng = np.random.default_rng(seed)
    return rng.normal(img, np.sqrt(img))


def _trace_slits(img, refine=False, bound=False):
    """Trace the slits with the Keck/LRIS (red) long-slit edge-tracing parameters."""
    spectrograph = load_spectrograph('keck_lris_red')
    par = spectrograph.default_pypeit_par()['calibrations']['slitedges']
    par['refine_edges'] = refine
    par['sync_predict'] = 'nearest'
    par['edge_thresh'] = 1000. if bound else 50.
    par['bound_detector'] = bound
    traceimg = TraceImage(img, ivar=1 / img)
    traceimg.detector = DetectorContainer(**tstutils.default_detector())
    return EdgeTraceSet(traceimg, spectrograph, par, auto=True).get_slits()


def test_traced_edges():
    # The edges inserted when only one edge of a slit is found, or when no
    # edges are found, are flagged as not traced
    for edges, bound, expected in [([(None, 200.37)], False, [[False, True]]),
                                   ([(180.62, None)], False, [[True, False]]),
                                   ([(100.3, 300.71)], False, [[True, True]]),
                                   ([(None, None)], True, [[False, False]])]:
        slits = _trace_slits(_longslit_image(edges, seed=2), bound=bound)
        assert np.array_equal(slits.traced_edges(), expected), \
            f'Incorrect traced-edge flags for slit edges {edges}'


@pytest.mark.parametrize('method', ['detector', 'slit', 'edge'])
def test_spat_flexure_single_edge(method):
    # A long slit with only its right edge on the detector.  The flexure is
    # measured from the refined traced edge, for all methods, and is assigned
    # to both edges of the slit
    edges = [(None, 200.37)]
    slits = _trace_slits(_longslit_image(edges, seed=2), refine=True)
    assert np.array_equal(slits.traced_edges(), [[False, True]]), 'Left edge should be inserted'
    # The last shift moves the inserted left edge off the detector
    for i, shift in enumerate([0.8, -6.6]):
        sci = _longslit_image(edges, sig=1.6, level=500., shift=shift, seed=10 + i)
        spat_flexure, gpm = flexure.spat_flexure_shift(sci, slits, method=method, maxlag=20,
                                                       sigdetect=5., return_gpm=True)
        assert np.allclose(spat_flexure, shift, atol=0.05), \
            f'Incorrect spatial flexure for a single slit edge and method={method}'
        assert np.all(gpm), 'The spatial flexure should be measured'
    assert slits.left_init.min() + spat_flexure[0, 0] < 0, \
        'Test is not sensitive: the inserted edge should be off the detector'
    slitmask = slits.slit_img(initial=True, spat_flexure=spat_flexure)
    assert np.all(slitmask[:, 0] == slits.spat_id[0]), \
        'The slit should extend to the detector edge'


@pytest.mark.parametrize('method', ['detector', 'edge'])
def test_spat_flexure_no_traced_edges(method):
    # No slit edges were found, so the slits bound the detector: the spatial
    # flexure cannot be measured.  The frame has a step in flux close to the
    # detector edge (e.g. vignetting), which must not be mistaken for a slit
    # edge.
    slits = _trace_slits(_longslit_image([(None, None)], seed=2), bound=True)
    assert np.array_equal(slits.traced_edges(), [[False, False]]), 'Edges should be inserted'
    sci = _longslit_image([(12.3, None)], level=500., seed=10)
    spat_flexure, gpm = flexure.spat_flexure_shift(sci, slits, method=method, maxlag=20,
                                                   sigdetect=5., return_gpm=True)
    assert np.all(spat_flexure == 0.) and not np.any(gpm), \
        'No spatial flexure should be applied without traced slit edges'


def test_spat_flexure_detector_balanced():
    # Three left and two right traced edges (one slit is cut by the detector
    # boundary), traced without refinement, so that the traced (and detected)
    # left and right edges are offset from the steepest gradient in opposite
    # directions by about a pixel.  The detector-level flexure averages the
    # left- and right-edge cross-correlation peaks, which cancels the offsets
    # whatever the number of edges on each side.
    edges = [(None, 70.4), (120.2, 220.8), (270.6, 360.1)]
    slits = _trace_slits(_longslit_image(edges, seed=2))
    assert np.array_equal(slits.traced_edges(), [[False, True], [True, True], [True, True]]), \
        'Incorrect traced-edge flags'
    for i, shift in enumerate([0.3, 2.4]):
        sci = _longslit_image(edges, level=500., shift=shift, seed=10 + i)
        spat_flexure = flexure.spat_flexure_shift(sci, slits, method='detector', maxlag=20,
                                                  sigdetect=5.)
        assert np.allclose(spat_flexure, shift, atol=0.1), 'Incorrect detector-level flexure'
