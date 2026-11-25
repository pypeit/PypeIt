"""
Tests for the PypeIt ginga view plugins, without a running viewer.
"""
import argparse

import ginga.toolkit
# The ginga app selects its GUI toolkit before any ginga widgets are imported;
# do the same here (see test_qlview_spectrographs.py).
ginga.toolkit.use('qt')

import numpy as np
import pytest

from pypeit.display import remote_handlers
views = pytest.importorskip('pypeit.display.views')


class FakeWidget:
    def __init__(self, index=0, state=False):
        self.index = index
        self.state = state

    def get_index(self):
        return self.index

    def get_state(self):
        return self.state


def _plugin(cls, args, options=None):
    # Bypass the ginga plugin constructor, which needs a running viewer
    plugin = object.__new__(cls)
    plugin.state = remote_handlers.ViewState('cmd', args, '/a/file.fits')
    plugin.state.options = {} if options is None else options
    plugin.w = {}
    return plugin


def _spec2d_args(**kwargs):
    args = dict(det=None, channels=None, showmask=None, show_wave_rms=False,
                ignore_extract_mask=False)
    args.update(kwargs)
    return argparse.Namespace(**args)


def test_spec2d_options():
    plugin = _plugin(views.Spec2dView, _spec2d_args(),
                     options=dict(detname='DET01', detectors=['DET01', 'MSC02']))
    plugin.w['det'] = FakeWidget(index=1)
    plugin.w['channels'] = [FakeWidget(state=s) for s in (True, False, False, True)]
    plugin.w['showmask'] = FakeWidget(state=True)
    plugin.w['show_wave_rms'] = FakeWidget(state=True)
    plugin.w['ignore_extract_mask'] = FakeWidget(state=False)
    args = plugin.options_args(_spec2d_args())
    assert args.det == 'MSC02', 'The selected detector should be used'
    assert args.channels == '0,3', 'Only the checked channels should be shown'
    assert args.showmask == [], 'Checking Show mask should show the mask bits'
    assert args.show_wave_rms is True, 'The wavelength RMS option should be set'

    plugin.w['channels'] = [FakeWidget(state=True) for i in range(4)]
    plugin.w['showmask'] = FakeWidget(state=False)
    args = plugin.options_args(_spec2d_args(showmask=['BPM']))
    assert args.channels is None, 'All channels checked should give the default'
    assert args.showmask is None, 'Unchecking Show mask should remove the mask channel'

    plugin.w['showmask'] = FakeWidget(state=True)
    args = plugin.options_args(_spec2d_args(showmask=['BPM']))
    assert args.showmask == ['BPM'], 'The mask flags given on the command line should be kept'


def test_flats_options():
    plugin = _plugin(views.FlatsView, argparse.Namespace(type='all'))
    plugin.w['type'] = FakeWidget(index=2)
    assert plugin.options_args(argparse.Namespace(type='all')).type == 'illum', \
        'The selected flat type should be used'


def test_slit_rows():
    left = np.tile(np.array([10., 50.]), (100,1))
    right = left + 20.
    rows = remote_handlers.slit_rows(left, right, slit_ids=[20, 60], maskdef_ids=[5])
    assert [r['id'] for r in rows] == [20, 60], 'Rows should carry the slit IDs'
    assert rows[1]['x'] == 60. and rows[1]['width'] == 20. and rows[1]['y'] == 50., \
        'Rows should give the slit center and width at the middle of the spectrum'
    assert rows[0]['maskdef_id'] is None, 'maskdef IDs of the wrong length should be ignored'


def test_make_overlays():
    left = np.tile(np.array([10., 50.]), (100,1))
    overlays = remote_handlers.make_overlays(slits=dict(left=left, right=left + 20.))
    assert list(overlays.keys()) == ['slit_edges', 'slit_ids'], \
        'Without traces, only the slit overlays should be built'
    kinds = {key: {o.kind for o in ov.canvas.get_objects()} for key, ov in overlays.items()}
    assert kinds == {'slit_edges': {'path'}, 'slit_ids': {'text'}}, \
        'Edges and IDs should be separate overlays so they can be toggled separately'


class FakeShell:
    def show_help_text(self, name, text):
        self.name = name
        self.text = text


@pytest.mark.parametrize('cls', [views.Spec2dView, views.EdgesView, views.FlatsView,
                                 views.CalibView, views.TiltsView, views.FileView])
def test_help_is_unique_rst(cls):
    plugin = object.__new__(cls)
    plugin.fv = FakeShell()
    plugin.help()
    assert plugin.fv.name == cls.__name__, 'The help tab should be named after the plugin'
    assert plugin.fv.text.startswith(f'{cls.__name__}\n{"=" * len(cls.__name__)}\n'), \
        'The help should start with the plugin name as its title'
    assert cls.__doc__.split('\n')[1].strip() in plugin.fv.text, \
        'The help should include the docstring of the plugin class'
    assert 'Control panel' in plugin.fv.text and '.. include' not in plugin.fv.text, \
        'The help should end with the shared panel description, without Sphinx directives'
