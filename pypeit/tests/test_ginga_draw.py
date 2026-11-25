"""
Tests for the in-process ginga drawing helpers.
"""
import numpy as np
import pytest

from ginga.AstroImage import AstroImage

from pypeit import PypeItError
from pypeit.display import ginga_draw


def _kinds(overlay):
    objs = overlay.get_objects()
    return [o.kind for o in objs], objs


def _edges(nspec=1000, nslits=3):
    spat = np.arange(nslits, dtype=float)*100 + 50
    left = np.tile(spat, (nspec,1))
    right = left + 40
    return left, right


def test_slit_overlay_synced():
    left, right = _edges()
    overlay = ginga_draw.slit_overlay(
        left, right, slit_ids=[10, 20, 30], maskdef_ids=[7, 8, 9], pstep=50
    )
    kinds, objs = _kinds(overlay)
    assert kinds.count('path') == 6, 'There should be one path per left and right edge'
    assert kinds.count('text') == 6, 'There should be a slit ID and a maskdef ID per slit'
    texts = [o.text for o in objs if o.kind == 'text']
    assert texts[:2] == ['S10', '7'], 'Slit IDs should be prefixed with S, maskdef IDs not'
    path = objs[0]
    assert len(path.points) == 20, 'Edges should be decimated by pstep'
    assert path.points[0] == (50., 0.), 'Without rotation, x is spatial and y spectral'
    assert path.color == 'green' and objs[3].color == 'magenta', \
        'Left edges should be green and right edges magenta'


def test_slit_overlay_rotate_and_no_maskdef():
    left, right = _edges(nslits=2)
    overlay = ginga_draw.slit_overlay(left, right, maskdef_ids=[1], rotate=True)
    kinds, objs = _kinds(overlay)
    assert kinds.count('text') == 2, \
        'maskdef IDs with the wrong length should be ignored, leaving one label per slit'
    assert objs[0].points[0] == (0., 50.), 'With rotation, x is spectral and y spatial'


def test_slit_overlay_unsynced():
    left, right = _edges(nslits=2)
    overlay = ginga_draw.slit_overlay(left, right[:,:1], synced=False)
    kinds, objs = _kinds(overlay)
    assert kinds.count('path') == 3, 'Each unsynced edge should be drawn'
    assert kinds.count('text') == 3, 'Each unsynced edge should be labeled'


def test_slit_overlay_errors():
    left, right = _edges()
    with pytest.raises(PypeItError):
        ginga_draw.slit_overlay(left, right[:500])
    with pytest.raises(PypeItError):
        ginga_draw.slit_overlay(left, right[:,:2])
    with pytest.raises(PypeItError):
        ginga_draw.slit_overlay(left, right, slit_ids=[1, 2])


def test_trace_overlay():
    trace = np.tile(np.array([100., 200., 300.]), (500,1))
    overlay = ginga_draw.trace_overlay(
        trace, trc_name=['a', 'b', 'c'], maskdef_extr=[False, True, False],
        manual_extr=[False, False, True]
    )
    kinds, objs = _kinds(overlay)
    assert kinds == ['path', 'text']*3, 'Each trace should have a path and a label'
    colors = [o.color for o in objs if o.kind == 'path']
    assert colors == ['orange', '#f0e442', '#33ccff'], \
        'Traces should be colored by how they were extracted'
    assert objs[1].text == 'a' and objs[1].y == 250., \
        'Labels should be placed at the middle of the trace'


def test_make_image():
    img = np.arange(12, dtype=float).reshape(3,4)
    image = ginga_draw.make_image(img, 'sci')
    assert type(image) is AstroImage, 'Without a wavelength image, a plain AstroImage is made'
    assert image.get('name') == 'sci', 'The image should carry the requested name'
    image = ginga_draw.make_image(img, 'sci', waveimg=img*10)
    assert isinstance(image, ginga_draw.SlitImage), 'A wavelength image should give a SlitImage'
    with pytest.raises(PypeItError):
        ginga_draw.make_image(img, 'sci', waveimg=img.T)


def test_tilt_overlays():
    from astropy.table import Table
    tbl = Table()
    tbl['good2dfit_lid'] = [np.array([0, 0, 0, 1, 1, 1])]
    tbl['good2dfit_spat'] = [np.arange(6, dtype=float)]
    tbl['good2dfit_tilt'] = [np.arange(6, dtype=float) + 100]
    tbl['goodpix_spat'] = [np.array([1., 2.])]
    tbl['goodpix_tilt'] = [np.array([3., 4.])]
    tbl['badpix_spat'] = [np.array([5.])]
    tbl['badpix_tilt'] = [np.array([6.])]
    tbl['bad2dfit_spat'] = [np.array([], dtype=float)]
    tbl['bad2dfit_tilt'] = [np.array([], dtype=float)]
    overlays = ginga_draw.tilt_overlays(tbl, 1000, pstep=1)
    assert set(overlays.keys()) == {'tilt_fits', 'tilt_masked', 'tilt_legend'}, \
        'Traced pixels are off by default, and overlays without data are omitted'
    assert len(overlays['tilt_fits'][1].get_objects()) == 2, 'There should be one path per line'
    legend = [o.text for o in overlays['tilt_legend'][1].get_objects()]
    assert legend == ['Masked pixel', 'Good tilt fit'], \
        'The legend should list only the overlays that are shown'
    overlays = ginga_draw.tilt_overlays(tbl, 1000, show_traces=True)
    assert len(overlays['tilt_pixels'][1].get_objects()) == 2, \
        'Each traced pixel should be marked when requested'
