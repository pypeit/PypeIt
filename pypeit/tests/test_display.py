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
    objs = overlay.build(MockViewer((0., 0., 500., 500.)))
    assert 0 < len(objs) <= 1000
    assert all(o.x <= 500 and o.y <= 500 for o in objs)


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
