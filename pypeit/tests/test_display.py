"""
Module to run tests on the ginga display tools
"""
import json
import socket

from IPython import embed

import numpy as np
import pytest

from astropy.io import fits

from pypeit.display import plugins_available
from pypeit.display import spec2d_data
from pypeit.display import viewer


def test_plugin():
    success, report = plugins_available(return_report=True)
    assert success, report


class MockViewer:
    """Minimal stand-in for a ginga viewer used to test the overlays."""
    def __init__(self, datarect, scale=1.0):
        self.datarect = datarect
        self.scale = scale

    def get_datarect(self):
        return self.datarect

    def get_scale_xy(self):
        return (self.scale, self.scale)

    def get_window_size(self):
        return (800, 600)


def _traces(ntrace, nspec=4096, spacing=10.):
    y = np.arange(nspec)
    return np.column_stack([20 + spacing*i + 2*np.sin(y/500) for i in range(ntrace)])


def test_trace_overlay_culls_and_decimates():
    from pypeit.display.overlays import TraceOverlay
    trc = _traces(20)
    overlay = TraceOverlay(trc, labels=[f'T{i}' for i in range(20)], max_vertices=50)
    # View a region including only some of the traces, and part of the rows
    objs = overlay.build(MockViewer((0., 1000., 75., 2000.), scale=2.))
    paths = [o for o in objs if o.kind == 'path']
    texts = [o for o in objs if o.kind == 'text']
    # Traces at 20, 30, ..., 80 (plus the margin) are in view
    assert 5 <= len(paths) < 20
    assert len(texts) == len(paths)
    # The number of vertices is limited
    assert all(len(p.points) <= 52 for p in paths)
    # Only the rows near the view are drawn
    assert all(np.min(np.asarray(p.points)[:, 1]) > 900 for p in paths)


def test_trace_overlay_level_of_detail():
    from pypeit.display.overlays import TraceOverlay
    trc = _traces(300, spacing=7.)
    vwr = MockViewer((0., 0., 2200., 4096.), scale=0.2)
    overlay = TraceOverlay(trc, labels=[f'F{i}' for i in range(300)], fontsize=12., owner=vwr)
    # Zoomed out: only the note is drawn
    objs = overlay.build(vwr)
    assert len(objs) == 1 and objs[0].kind == 'text' and '300' in objs[0].text
    # ... and nothing in other viewers (e.g., the Pan plugin)
    assert len(overlay.build(MockViewer((0., 0., 2200., 4096.), scale=0.2))) == 0
    # Zoomed in, but traces too close for labels
    vwr.datarect, vwr.scale = (500., 1000., 700., 1200.), 1.5
    objs = overlay.build(vwr)
    assert len(objs) > 0 and all(o.kind == 'path' for o in objs)
    # Zoomed in enough to show labels
    vwr.datarect, vwr.scale = (500., 1000., 600., 1100.), 4.
    objs = overlay.build(vwr)
    assert any(o.kind == 'text' for o in objs)


def test_trace_overlay_declutter():
    from pypeit.display.overlays import TraceOverlay
    # Traces 5 pixels apart; at a scale of 2, labels (fontsize 12) need
    # ~14.4 screen pixels, so only every other label is drawn
    trc = _traces(20, spacing=5.)
    vwr = MockViewer((0., 0., 200., 4096.), scale=2.)
    overlay = TraceOverlay(trc, labels=[f'T{i}' for i in range(20)], fontsize=12., owner=vwr)
    texts = [o for o in overlay.build(vwr) if o.kind == 'text']
    assert 5 <= len(texts) <= 10
    x = np.sort([t.x for t in texts])
    assert np.all(np.diff(x) * 2. >= 14.4)
    # No labels in other viewers
    other = MockViewer((0., 0., 200., 4096.), scale=2.)
    assert not any(o.kind == 'text' for o in overlay.build(other))


def test_trace_overlay_nearest():
    from pypeit.display.overlays import TraceOverlay
    trc = _traces(10)
    overlay = TraceOverlay(trc, labels=[f'T{i}' for i in range(10)])
    vwr = MockViewer((0., 0., 200., 4096.), scale=1.)
    assert overlay.nearest(vwr, trc[1000, 3] + 1., 1000.) == 3
    assert overlay.nearest(vwr, 1000., 1000.) is None
    assert overlay.nearest(vwr, trc[0, 3], -100.) is None
    # The hovered trace is highlighted and labeled
    overlay.hover = 3
    objs = overlay.build(vwr)
    assert any(o.kind == 'path' and o.color == 'cyan' for o in objs)
    assert any(o.kind == 'text' and o.text == 'T3' and o.color == 'cyan' for o in objs)


def test_marker_overlay():
    from pypeit.display.overlays import MarkerOverlay
    rng = np.random.default_rng(42)
    x, y = rng.uniform(0, 1000, size=(2, 20000))
    overlay = MarkerOverlay(x, y, max_markers=1000)
    # Only markers in view are selected, at most max_markers
    indx = overlay.select(MockViewer((0., 0., 500., 500.), scale=4.))
    assert 0 < indx.size <= 1000
    assert np.all(x[indx] <= 500) and np.all(y[indx] <= 500)
    # Zoomed out, at most one marker per screen patch (3 pixels) is selected
    indx = overlay.select(MockViewer((0., 0., 1000., 1000.), scale=0.1))
    assert indx.size <= (1000 * 0.1 / 3 + 1)**2
    cells = np.floor(np.column_stack([x[indx], y[indx]]) * 0.1 / 3).astype(int)
    assert len(np.unique(cells, axis=0)) == indx.size


class MockShell:
    """Minimal stand-in for the ginga shell used to test the command server."""
    class ds:
        toplevels = []

    def __init__(self):
        self.errors = []

    def show_error(self, msg):
        self.errors.append(msg)


def _request(pv, payload):
    """Send a request to a viewer's command handler and return the reply."""
    client, server = socket.socketpair()
    with client:
        client.sendall(payload)
        client.shutdown(socket.SHUT_WR)
        # NOTE: As in PypeItViewer._poll_server, the connection is closed
        # after it is handled.
        with server:
            pv._handle_connection(server)
        reply = client.recv(4096)
    return json.loads(reply.decode()) if len(reply) > 0 else None


def test_command_protocol():
    # NOTE: Bypass __init__, which would create the ginga GUI
    pv = viewer.PypeItViewer.__new__(viewer.PypeItViewer)
    pv.shell = MockShell()
    pv._commands = {}
    received = {}
    pv.register_command('test', lambda **kwargs: received.update(kwargs))

    msg = json.dumps(dict(cmd='test', kwargs=dict(file='spec2d.fits', det=1))) + '\n'
    assert _request(pv, msg.encode()) == dict(status='ok')
    assert received == dict(file='spec2d.fits', det=1)

    # Unknown commands and malformed requests are rejected
    msg = json.dumps(dict(cmd='not_a_command', kwargs={})) + '\n'
    assert _request(pv, msg.encode())['status'] == 'error'
    assert _request(pv, b'not json\n')['status'] == 'error'

    # Errors executing a command are reported by the viewer
    pv.register_command('fail', lambda **kwargs: 1/0)
    msg = json.dumps(dict(cmd='fail', kwargs={})) + '\n'
    assert _request(pv, msg.encode()) == dict(status='ok')
    assert len(pv.shell.errors) == 1

    # Connections that only check if the viewer is listening get no reply
    assert _request(pv, b'') is None


def test_send_command_no_viewer():
    # Find a free port
    with socket.socket() as s:
        s.bind(('127.0.0.1', 0))
        port = s.getsockname()[1]
    assert not viewer.is_listening(port=port)
    assert not viewer.send_command('spec2d', port=port, file='spec2d.fits')


def test_parse_detname(tmp_path):
    ofile = tmp_path / 'spec2d_test.fits'
    hdr = fits.Header()
    hdr['HIERARCH ALLSPEC2D_DETS'] = 'DET03,MSC01'
    fits.HDUList([fits.PrimaryHDU(header=hdr)]).writeto(ofile)
    assert spec2d_data.get_detnames(ofile) == ['DET03', 'MSC01']
    assert spec2d_data.parse_detname(ofile) == 'DET03'
    assert spec2d_data.parse_detname(ofile, det=2) == 'DET02'
    assert spec2d_data.parse_detname(ofile, det='2') == 'DET02'
    assert spec2d_data.parse_detname(ofile, det='MSC01') == 'MSC01'


def test_cuts():
    rng = np.random.default_rng(42)
    img = rng.normal(size=(200, 200))
    gpm = np.ones(img.shape, dtype=bool)
    lo, hi = spec2d_data._cuts(img, gpm)
    assert -1.5 < lo < -0.5 and 3.5 < hi < 4.5
    assert spec2d_data._cuts(img, np.zeros(img.shape, dtype=bool)) is None


def test_ginga_hidpi_plot_patch():
    # NOTE: Importing the overlays chooses ginga's widget toolkit, which must be
    # done before importing ginga's plot viewer.
    from pypeit.display import overlays
    from ginga.plot.PlotView import PlotViewBase
    viewer._patch_ginga_hidpi_plots()
    assert PlotViewBase._pypeit_hidpi_patch
    assert PlotViewBase.configure is PlotViewBase.set_window_size
    assert PlotViewBase.configure_window is PlotViewBase.set_window_size


def test_spec2d_make_channel():
    rng = np.random.default_rng(42)
    shape = (100, 50)
    sciimg = rng.normal(size=shape) + 10.
    gpm = np.ones(shape, dtype=bool)
    data = spec2d_data.Spec2DDisplay(file=None, detname='DET01', detnames=['DET01'],
                                     prefix='A_', removetrace=True, sciimg=sciimg,
                                     skymodel=np.full(shape, 10.), objmodel=None,
                                     ivarmodel=np.ones(shape), img_gpm=gpm, model_gpm=gpm,
                                     bpmmask=np.zeros(shape, dtype=np.int16))
    chan = data.make_channel(0)
    assert chan.name == 'A_sciimg-DET01' and chan.show_traces and chan.key == 0
    chan = data.make_channel(1)
    assert chan.name == 'A_skysub-DET01' and not chan.show_traces
    assert np.allclose(chan.image, sciimg - 10.)
    assert data.make_channel(2).cuts == (-5.0, 5.0)
    # No object model, so no residual image
    assert data.make_channel(3) is None
    chan = data.make_channel('mask')
    assert chan.name == 'A_MASK' and not chan.show_slits


def test_viewer_headless():
    # Smoke test of the in-process viewer.  Under pytest (and without a
    # display), show() must not wait for the window to be closed.
    import time
    from pypeit.display import display
    rng = np.random.default_rng(42)
    img = rng.normal(size=(200, 100))
    waveimg = np.tile(np.linspace(4000., 5000., 200)[:, None], (1, 100))
    left = np.column_stack([np.full(200, 10.), np.full(200, 50.)])
    right = left + 30.
    vwr, ch = display.show_image(img, chname='test', waveimg=waveimg, cuts=(-1., 1.), clear=True)
    display.show_slits(vwr, ch, left, right, slit_ids=np.array([1, 2]))
    display.show_trace(vwr, ch, left + 15., np.array(['obj1', 'obj2']))
    vwr2, ch2 = display.show_image(img, chname='test2', wcs_match=True)
    assert vwr2 is vwr
    t = time.perf_counter()
    display.show()
    assert time.perf_counter() - t < 30.
    assert {'test', 'test2'} <= set(vwr.channel_names())
    # Left edges, right edges, slit IDs, and traces
    assert len(ch.fitsimage.get_canvas().get_tags_by_tag_pfx(display.TAG_PREFIX)) == 4
    # The channels are linked
    ch.fitsimage.set_pan(20., 100.)
    assert np.allclose(ch2.fitsimage.get_pan()[:2], (20., 100.))


def test_is_onespec_file(tmp_path):
    from pypeit.onespec import OneSpec, is_onespec_file
    wave = np.linspace(4000., 5000., 100)
    spec = OneSpec(wave, wave, np.ones_like(wave), PYP_SPEC='shane_kast_blue',
                   ivar=np.ones_like(wave), mask=np.ones(wave.size, dtype=int), ext_mode='OPT',
                   fluxed=False)
    ofile = tmp_path / 'onespec.fits'
    spec.to_file(str(ofile), overwrite=True)
    assert is_onespec_file(ofile)
    # A file that is not a OneSpec file
    other = tmp_path / 'other.fits'
    fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(np.zeros((2, 2)))]).writeto(other)
    assert not is_onespec_file(other)


def test_spec1d_key_bindings_documented():
    # The bindings of the spec1d mode must be listed in the text shown by the
    # "Keys" button of Spec1dView (the docstring of Spec1DMode)
    import logging
    # NOTE: Importing the overlays chooses ginga's widget toolkit, which must be
    # done before importing ginga's plot modes.
    from pypeit.display import overlays
    from pypeit.display.pypeit_modes import Spec1DMode, key_bindings_text

    class _Viewer:
        def get_logger(self):
            return logging.getLogger('test')

    actions = Spec1DMode(_Viewer(), settings={}).actions
    text = key_bindings_text()
    for action, key, description in [('kp_help', 'h', '* h :'),
                                     ('kp_smooth_more', 's', '* s :'),
                                     ('kp_smooth_less', 'a', '* a :'),
                                     ('kp_no_smooth', 'A', '* A (shift + a) :'),
                                     ('kp_cut_x_range', 'singlequote', '* singlequote :'),
                                     ('kp_fitting', 'f', '* f :'),
                                     ('kp_clear_marks', 'c', '* c :')]:
        assert f'spec1d+{key}' in actions[action], action
        assert description in text, description
