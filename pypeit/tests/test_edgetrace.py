"""
Tests of the slit-edge tracing (EdgeTraceSet).
"""
import numpy as np

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
