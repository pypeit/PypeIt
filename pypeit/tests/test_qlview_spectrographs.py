"""
Characterization tests for the quicklook-viewer hooks on the spectrograph classes.

Results are compared by column *display name* so that the internal attribute
names used to key each column are free to change.
"""
import logging
import subprocess
import sys

from astropy.io import fits
import pytest

from pypeit.display.qlview import calib_utils, spectrograph_support
from pypeit.spectrographs.util import load_spectrograph, spectrograph_classes

# ---------------------------------------------------------------------------
# Adapters onto the API under test
# ---------------------------------------------------------------------------

_LOGGER = logging.getLogger(__name__)


def _columns(name, mode):
    return [d for d, _ in spectrograph_support.build_columns(load_spectrograph(name), mode)]


def _view(name, mode, path):
    spec = load_spectrograph(name)
    info = spectrograph_support.get_header_info(spec, str(path), mode, _LOGGER)
    return {d: str(info[a]) for d, a in spectrograph_support.build_columns(spec, mode)
            if a in info}


def _header_name(name):
    return load_spectrograph(name).header_name


def _det_prefix(name):
    return spectrograph_support.det_label(load_spectrograph(name), '01')[:-2]


# ---------------------------------------------------------------------------
# Fixtures and expected values
# ---------------------------------------------------------------------------

_HEADER_KEYS = [
    'OBJECT', 'TARGNAME', 'FRAMENO', 'FRAMENUM', 'KOAIMTYP', 'IMTYPE', 'MASKNAME', 'SLMSKNAM',
    'SLITNAME', 'OBSMODE', 'EXPTIME', 'ELAPTIME', 'TRUITIME', 'TTIME', 'GRATENAM', 'DWFILNAM',
    'DECKNAME', 'XDISPERS', 'GRISNAME', 'GRANAME', 'DICHNAME', 'PATTERN', 'FRAMEID', 'SCIFILT1',
    'SCIFILT2', 'FIL1NAME', 'FILTER', 'FILTER1', 'FILTER2', 'MGTNAME', 'SLIT',
]

_CONFIG_KEYS = ['decker', 'decker_secondary', 'dispname', 'filter1', 'filter2', 'dichroic',
                'slitwid']

NAMES = ['keck_deimos', 'keck_hires', 'keck_lris_blue', 'keck_lris_red_mark4', 'keck_mosfire',
         'keck_nirspec_high']

EXPECTED = {
    'keck_deimos': dict(
        header='DEIMOS', prefix='MSC',
        raw_cols=['Type', 'Frame No', 'Name', 'Object', 'Img Type', 'Mask Name', 'Grating',
                  'Filter', 'Exp Time', 'Last Changed'],
        red_cols=['Type', 'Name', 'Mask/Slit', 'Grating', 'Blocking Filter', 'Last Changed'],
        raw_full={'Frame No': 'v_FRAMENO', 'Object': 'v_TARGNAME', 'Img Type': 'v_KOAIMTYP',
                  'Mask Name': 'v_SLMSKNAM', 'Grating': 'v_GRATENAM', 'Filter': 'v_DWFILNAM',
                  'Exp Time': 'v_ELAPTIME'},
        red_full={'Mask/Slit': 'v_MASKNAME', 'Grating': 'v_GRATENAM',
                  'Blocking Filter': 'v_SLITNAME'},
        red_dir={'Mask/Slit': 'c_decker', 'Grating': 'c_dispname',
                 'Blocking Filter': 'c_filter1'},
    ),
    'keck_hires': dict(
        header='HIRES', prefix='MSC',
        raw_cols=['Type', 'Frame No', 'Name', 'Object', 'Img Type', 'Decker', 'XDisp',
                  'Exp Time', 'Last Changed'],
        red_cols=['Type', 'Name', 'Decker', 'XDisp', 'Filter', 'Last Changed'],
        raw_full={'Frame No': 'v_FRAMENO', 'Object': 'v_TARGNAME', 'Img Type': 'v_KOAIMTYP',
                  'Decker': 'v_DECKNAME', 'XDisp': 'v_XDISPERS', 'Exp Time': 'v_ELAPTIME'},
        red_full={'Decker': 'v_DECKNAME', 'XDisp': 'v_XDISPERS', 'Filter': 'v_FIL1NAME'},
        red_dir={'Decker': 'c_decker', 'XDisp': 'c_dispname', 'Filter': 'c_filter1'},
    ),
    'keck_lris_blue': dict(
        # NOTE: The original viewer used 'MSC' here, but LRIS Blue is displayed
        # with mosaic=False, so 'DET' is correct.
        header='LRISBLUE', prefix='DET',
        raw_cols=['Type', 'Frame No', 'Name', 'Object', 'Img Type', 'Slit/Mask', 'Grism',
                  'Dichroic', 'Exp Time', 'Last Changed'],
        red_cols=['Type', 'Name', 'Slit/Mask', 'Grating/Grism', 'Dichroic', 'Last Changed'],
        raw_full={'Frame No': 'v_FRAMENO', 'Object': 'v_TARGNAME', 'Img Type': 'v_KOAIMTYP',
                  'Slit/Mask': 'v_SLITNAME', 'Grism': 'v_GRISNAME', 'Dichroic': 'v_DICHNAME',
                  'Exp Time': 'v_EXPTIME'},
        red_full={'Slit/Mask': 'v_SLITNAME', 'Grating/Grism': 'v_GRISNAME',
                  'Dichroic': 'v_DICHNAME'},
        red_dir={'Slit/Mask': 'c_decker', 'Grating/Grism': 'c_dispname',
                 'Dichroic': 'c_dichroic'},
    ),
    'keck_lris_red_mark4': dict(
        header='LRIS', prefix='DET',
        raw_cols=['Type', 'Frame No', 'Name', 'Object', 'Img Type', 'Slit/Mask', 'Grating',
                  'Dichroic', 'Exp Time', 'Last Changed'],
        red_cols=['Type', 'Name', 'Slit/Mask', 'Grating/Grism', 'Dichroic', 'Last Changed'],
        raw_full={'Frame No': 'v_FRAMENO', 'Object': 'v_TARGNAME', 'Img Type': 'v_KOAIMTYP',
                  'Slit/Mask': 'v_SLITNAME', 'Grating': 'v_GRANAME', 'Dichroic': 'v_DICHNAME',
                  'Exp Time': 'v_EXPTIME'},
        red_full={'Slit/Mask': 'v_SLITNAME', 'Grating/Grism': 'v_GRANAME',
                  'Dichroic': 'v_DICHNAME'},
        red_dir={'Slit/Mask': 'c_decker', 'Grating/Grism': 'c_dispname',
                 'Dichroic': 'c_dichroic'},
    ),
    'keck_mosfire': dict(
        header='MOSFIRE', prefix='DET',
        raw_cols=['Type', 'Frame No', 'Name', 'Dither Pos', 'Object', 'Img Type', 'Mask Name',
                  'Obs Mode', 'Exp Time', 'Last Changed'],
        red_cols=['Type', 'Name', 'CSU Mask', 'Filter', 'Dispname', 'Slit Width',
                  'Last Changed'],
        raw_full={'Frame No': 'v_FRAMENO', 'Dither Pos': 'v_FRAMEID', 'Object': 'v_OBJECT',
                  'Img Type': 'v_KOAIMTYP', 'Mask Name': 'v_MASKNAME', 'Obs Mode': 'v_OBSMODE',
                  'Exp Time': 'v_EXPTIME'},
        red_full={'CSU Mask': 'v_MASKNAME', 'Filter': 'v_FILTER', 'Dispname': 'v_FILTER2',
                  'Slit Width': 'v_MGTNAME'},
        red_dir={'CSU Mask': 'c_decker_secondary', 'Filter': 'c_filter1',
                 'Dispname': 'c_dispname', 'Slit Width': 'c_slitwid'},
    ),
    'keck_nirspec_high': dict(
        header='NIRSPEC', prefix='DET',
        raw_cols=['Type', 'Frame No', 'Name', 'Object', 'Img Type', 'Filter 1', 'Filter 2',
                  'Slit', 'Exp Time', 'Last Changed'],
        red_cols=['Type', 'Name', 'Slit', 'Filter 1', 'Filter 2', 'Last Changed'],
        raw_full={'Frame No': 'v_FRAMENUM', 'Object': 'v_TARGNAME', 'Img Type': 'v_IMTYPE',
                  'Filter 1': 'v_SCIFILT1', 'Filter 2': 'v_SCIFILT2', 'Slit': 'v_SLITNAME',
                  'Exp Time': 'v_TRUITIME'},
        red_full={'Slit': 'v_SLITNAME', 'Filter 1': 'v_SCIFILT1', 'Filter 2': 'v_SCIFILT2'},
        red_dir={'Slit': 'c_decker', 'Filter 1': 'c_filter1', 'Filter 2': 'c_filter2'},
    ),
}


@pytest.fixture
def fits_files(tmp_path):
    full = fits.Header()
    for k in _HEADER_KEYS:
        full[k] = f'v_{k}'
    full_path = tmp_path / 'full.fits'
    empty_path = tmp_path / 'empty.fits'
    fits.PrimaryHDU(header=full).writeto(full_path)
    fits.PrimaryHDU().writeto(empty_path)
    return full_path, empty_path


def _write_calib_dir(root, name):
    """Write a minimal PypeIt setup directory, ``{name}_A/``, with a .pypeit file."""
    setup_dir = root / f'{name}_A'
    setup_dir.mkdir()
    cfg = '\n'.join(f'    {k}: c_{k}' for k in _CONFIG_KEYS)
    (setup_dir / f'{name}_A.pypeit').write_text(
        '[rdx]\n'
        f'    spectrograph = {name}\n'
        '\n'
        'setup read\n'
        'Setup A:\n'
        '  --:\n'
        f'{cfg}\n'
        'setup end\n'
        '\n'
        'data read\n'
        ' path ' + str(root) + '\n'
        '|    filename |   frametype |\n'
        '| raw_01.fits |     science |\n'
        'data end\n'
    )
    return setup_dir


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('name', NAMES)
def test_columns(name):
    assert _columns(name, 'raw') == EXPECTED[name]['raw_cols']
    assert _columns(name, 'reduced') == EXPECTED[name]['red_cols']


@pytest.mark.parametrize('name', NAMES)
def test_header_name_and_det_prefix(name):
    assert _header_name(name) == EXPECTED[name]['header']
    assert _det_prefix(name) == EXPECTED[name]['prefix']


@pytest.mark.parametrize('name', NAMES)
def test_raw_info(name, fits_files):
    full, empty = fits_files
    assert _view(name, 'raw', full) == EXPECTED[name]['raw_full']
    expected_empty = {k: 'N/A' for k in EXPECTED[name]['raw_full']}
    assert _view(name, 'raw', empty) == expected_empty


def test_mosfire_stare_dither(tmp_path):
    hdr = fits.Header({'PATTERN': 'Stare', 'FRAMEID': 'A'})
    path = tmp_path / 'stare.fits'
    fits.PrimaryHDU(header=hdr).writeto(path)
    assert _view('keck_mosfire', 'raw', path)['Dither Pos'] == 'N/A'


@pytest.mark.parametrize('name', NAMES)
def test_reduced_info_fits(name, fits_files):
    full, empty = fits_files
    assert _view(name, 'reduced', full) == EXPECTED[name]['red_full']
    expected_empty = {k: 'N/A' for k in EXPECTED[name]['red_full']}
    assert _view(name, 'reduced', empty) == expected_empty


@pytest.mark.parametrize('name', NAMES)
def test_reduced_info_calib_dir(name, tmp_path):
    setup_dir = _write_calib_dir(tmp_path, name)
    assert _view(name, 'reduced', setup_dir) == EXPECTED[name]['red_dir']


def test_reduced_info_symlinked_calib_dir(tmp_path):
    target = _write_calib_dir(tmp_path, 'keck_deimos')
    link = tmp_path / 'links' / 'keck_deimos_B'
    link.parent.mkdir()
    link.symlink_to(target)
    assert _view('keck_deimos', 'reduced', link) == EXPECTED['keck_deimos']['red_dir']


def test_setup_config_layouts(tmp_path):
    # Configuration keys directly under the setup name, with a block-style
    # per-detector sub-block that must not leak into the configuration.
    setup_dir = tmp_path / 'keck_deimos_A'
    setup_dir.mkdir()
    (setup_dir / 'keck_deimos_A.pypeit').write_text(
        'setup read\n'
        'Setup A:\n'
        '  decker: c_decker\n'
        '  dispname: c_dispname\n'
        "  '01':\n"
        '    binning: 2,2\n'
        'setup end\n'
    )
    config = calib_utils.read_pypeit_setup_config(str(setup_dir), 'keck_deimos', _LOGGER)
    assert config == {'decker': 'c_decker', 'dispname': 'c_dispname'}


def test_reduced_info_unrelated_dir(tmp_path):
    other = tmp_path / 'Calibrations'
    other.mkdir()
    assert _view('keck_deimos', 'reduced', other) \
            == {'Mask/Slit': 'N/A', 'Grating': 'N/A', 'Blocking Filter': 'N/A'}


# ---------------------------------------------------------------------------
# Viewer support
# ---------------------------------------------------------------------------

def test_supported_spectrographs():
    assert spectrograph_support.supported_spectrographs() == ('keck_deimos', 'keck_mosfire')
    for name in spectrograph_support.supported_spectrographs():
        spec = load_spectrograph(name)
        assert spec.header_name is not None
        assert spec.qlview_label is not None
        for mode in ['raw', 'reduced']:
            cols = spectrograph_support.build_columns(spec, mode)
            assert all(isinstance(d, str) and isinstance(a, str) for d, a in cols)
            assert len({a for _, a in cols}) == len(cols)


@pytest.mark.parametrize('name', sorted(spectrograph_classes().keys()))
def test_qlview_raw_info_empty_header(name):
    # Every spectrograph must tolerate missing header keywords
    spec = load_spectrograph(name)
    info = spec.qlview_raw_info(fits.Header())
    assert all(info[key] == 'N/A' for _, key in spec.qlview_raw_columns())


def test_load_qlview_spectrograph():
    assert load_spectrograph('keck_mosfire').name \
            == spectrograph_support.load_qlview_spectrograph('keck_mosfire').name
    for name in ['keck_lris_blue', 'DEIMOS', 'not_a_spectrograph']:
        with pytest.raises(ValueError):
            spectrograph_support.load_qlview_spectrograph(name)


def test_match_header_name():
    assert spectrograph_support.match_header_name(' deimos ') == 'keck_deimos'
    assert spectrograph_support.match_header_name('MOSFIRE') == 'keck_mosfire'
    # LRIS is not supported by the viewer
    assert spectrograph_support.match_header_name('LRIS') is None
    assert spectrograph_support.label('keck_deimos') == 'DEIMOS'


def test_no_viewer_import_from_spectrographs():
    # The spectrograph classes must not depend on the viewer or on
    # pypeit.scripts.ql, which itself imports the spectrographs.
    code = ('import sys; import pypeit.spectrographs.spectrograph; '
            'from pypeit.spectrographs.util import load_spectrograph; '
            '[load_spectrograph(n) for n in ["keck_deimos", "keck_mosfire"]]; '
            'bad = [m for m in sys.modules '
            'if m.startswith("pypeit.display") or m == "pypeit.scripts.ql"]; '
            'print(bad); assert not bad')
    subprocess.run([sys.executable, '-c', code], check=True)


# ---------------------------------------------------------------------------
# HTTP server
# ---------------------------------------------------------------------------

@pytest.fixture
def http_client():
    pytest.importorskip('flask')
    from pypeit.display.qlview.servers import HTTPserver
    HTTPserver.app.config['TESTING'] = True
    return HTTPserver.app.test_client()


def test_http_header_info(http_client, fits_files):
    full, _ = fits_files
    resp = http_client.get('/api/header_info',
                           query_string={'path': str(full), 'instrument': 'keck_deimos'})
    assert resp.status_code == 200
    assert resp.get_json()['MASKNAME'] == 'v_SLMSKNAM'


@pytest.mark.parametrize('instrument', ['DEIMOS', 'keck_lris_blue', ''])
def test_http_header_info_bad_instrument(http_client, fits_files, instrument):
    full, _ = fits_files
    resp = http_client.get('/api/header_info',
                           query_string={'path': str(full), 'instrument': instrument})
    assert resp.status_code == 400


# ---------------------------------------------------------------------------
# Calibration version checks
# ---------------------------------------------------------------------------

def _write_calib_header(path, cls_name, dmodver, verspyp='9.9.9'):
    """Write a header-only stand-in for a calibration file."""
    ext = fits.ImageHDU(name='SLITS')
    ext.header['DMODCLS'] = cls_name
    ext.header['DMODVER'] = dmodver
    fits.HDUList([fits.PrimaryHDU(header=fits.Header({'VERSPYP': verspyp})), ext]).writeto(path)


def test_check_calib_versions(tmp_path):
    from pypeit.slittrace import SlitTraceSet
    from pypeit.wavecalib import WaveCalib
    _write_calib_header(tmp_path / 'Slits_A_0_DET01.fits', 'SlitTraceSet', '0.0.1', '1.0.0')
    _write_calib_header(tmp_path / 'WaveCalib_A_0_DET01.fits', 'WaveCalib', WaveCalib.version)
    (tmp_path / 'Tilts_A_0_DET01.fits').write_text('not a fits file')

    mismatches = calib_utils.check_calib_versions(str(tmp_path))
    assert mismatches == [dict(file='Slits_A_0_DET01.fits', cls='SlitTraceSet',
                               file_version='0.0.1', code_version=SlitTraceSet.version,
                               pypeit_version='1.0.0')]
    # Only the requested file types are checked
    assert calib_utils.check_calib_versions(str(tmp_path), prefixes=['WaveCalib']) == []


def test_check_calib_versions_current(tmp_path):
    from pypeit.slittrace import SlitTraceSet
    _write_calib_header(tmp_path / 'Slits_A_0_DET01.fits', 'SlitTraceSet', SlitTraceSet.version)
    assert calib_utils.check_calib_versions(str(tmp_path)) == []
    assert calib_utils.check_calib_versions(str(tmp_path / 'missing')) == []


def test_warn_outdated_calibrations(monkeypatch):
    qlview = pytest.importorskip('pypeit.display.qlview.qlview')
    shown = []
    monkeypatch.setattr(qlview.QtGui.QMessageBox, 'warning',
                        lambda parent, title, msg: shown.append((title, msg)))

    class _Plugin:
        logger = _LOGGER

    mismatches = [dict(file=f'Slits_A_0_MSC0{i}.fits.gz', cls='SlitTraceSet',
                       file_version='1.1.5', code_version='1.1.6', pypeit_version='1.18.2')
                  for i in range(1, 11)]
    qlview.QLView._warn_outdated_calibrations(_Plugin(), '/cals', mismatches=mismatches)
    title, msg = shown[-1]
    assert title == 'Calibrations Out of Date'
    assert 'Slits_A_0_MSC01.fits.gz: written by PypeIt 1.18.2' in msg
    assert 'SlitTraceSet 1.1.5; this PypeIt expects 1.1.6' in msg
    assert '... and 2 more' in msg

    qlview.QLView._warn_outdated_calibrations(_Plugin(), '/cals',
                                              error=ValueError('version mismatch'))
    assert 'version mismatch' in shown[-1][1]


# ---------------------------------------------------------------------------
# Calibration-directory selection
# ---------------------------------------------------------------------------

@pytest.fixture
def calib_tree(tmp_path):
    """reduced/keck_deimos_A/{keck_deimos_A.pypeit, Calibrations/Slits_A_0_MSC01.fits}"""
    setup_dir = tmp_path / 'reduced' / 'keck_deimos_A'
    calib_dir = setup_dir / 'Calibrations'
    calib_dir.mkdir(parents=True)
    (setup_dir / 'keck_deimos_A.pypeit').write_text('')
    (calib_dir / 'Slits_A_0_MSC01.fits').write_text('')
    (tmp_path / 'reduced' / 'Calibrations_old').mkdir()
    return setup_dir, calib_dir


def test_resolve_calib_dirs(calib_tree):
    setup_dir, calib_dir = calib_tree
    expected = (setup_dir, calib_dir)
    for path in [setup_dir, f'{setup_dir}/', f'{setup_dir}/*', calib_dir, f'{calib_dir}/*',
                 calib_dir / 'Slits_A_0_MSC01.fits', setup_dir / 'keck_deimos_A.pypeit']:
        assert calib_utils.resolve_calib_dirs(str(path)) == expected, path


def test_resolve_calib_dirs_invalid(calib_tree):
    setup_dir, _ = calib_tree
    # The directory holding the calibration sets is not itself a set
    for path in [None, '', str(setup_dir.parent), f'{setup_dir.parent}/*',
                 str(setup_dir.parent / 'Calibrations_old'), str(setup_dir / 'missing')]:
        assert calib_utils.resolve_calib_dirs(path) is None, path


def test_update_calib_buttons(calib_tree):
    qlview = pytest.importorskip('pypeit.display.qlview.qlview')
    from pypeit.display.qlview.state import QLViewState
    setup_dir, calib_dir = calib_tree

    class _Widget:
        def __init__(self): self.enabled, self.text = None, ''
        def set_enabled(self, val): self.enabled = val
        def set_text(self, val): self.text = val
        def get_text(self): return self.text

    class _Plugin:
        _current_calib_dirs = qlview.QLView._current_calib_dirs
        _update_calib_buttons = qlview.QLView._update_calib_buttons
        _set_cal_status = qlview.QLView._set_cal_status
        def __init__(self):
            self.state = QLViewState()
            self.spectrograph = load_spectrograph('keck_deimos')
            self.reduced_btn, self.show_wavelengths_btn = _Widget(), _Widget()
            self.cal_status_label, self.reduced_text_entry = _Widget(), _Widget()
            self._showing_calib_hint = False

    plugin = _Plugin()
    # The parent of the calibration sets: disabled, with a hint
    plugin.state.reduced_filepath = f'{setup_dir.parent}/*'
    plugin._update_calib_buttons()
    assert not plugin.reduced_btn.enabled and not plugin.show_wavelengths_btn.enabled
    assert 'keck_deimos_A' in plugin.cal_status_label.text

    # Calibrations/ itself: enabled, and the hint is cleared
    plugin.state.reduced_filepath = str(calib_dir)
    plugin._update_calib_buttons()
    assert plugin.reduced_btn.enabled and plugin.show_wavelengths_btn.enabled
    assert plugin.cal_status_label.text == ''
    assert plugin._current_calib_dirs() == (setup_dir, calib_dir)

    # Other status messages are not cleared
    plugin.cal_status_label.text = 'Calibrations matched'
    plugin.state.reduced_filepath = str(setup_dir)
    plugin._update_calib_buttons()
    assert plugin.cal_status_label.text == 'Calibrations matched'

    # A status message that replaces the hint is not cleared by a later
    # valid selection
    plugin.state.reduced_filepath = f'{setup_dir.parent}/*'
    plugin._update_calib_buttons()
    plugin._set_cal_status('No calibrations found')
    plugin.state.reduced_filepath = str(setup_dir)
    plugin._update_calib_buttons()
    assert plugin.cal_status_label.text == 'No calibrations found'

    # A typed path is used when nothing has been selected
    plugin.state.reduced_filepath = None
    plugin.reduced_text_entry.text = str(setup_dir)
    assert plugin._current_calib_dirs() == (setup_dir, calib_dir)


def test_suggest_calibrations_search_root(calib_tree, monkeypatch):
    # Matching searches for setup directories *inside* the search root, so
    # selecting a setup directory must not make it the root
    qlview = pytest.importorskip('pypeit.display.qlview.qlview')
    from pypeit.display.qlview.state import QLViewState
    setup_dir, calib_dir = calib_tree
    roots = []
    monkeypatch.setattr(qlview, 'recommend_calibrations',
                        lambda name, raw, root, log: roots.append(root) or [])

    class _Thread:
        def __init__(self, target, **kwargs): self.target = target
        def start(self): self.target()
    monkeypatch.setattr(qlview.threading, 'Thread', _Thread)

    class _Widget:
        def __init__(self): self.text = ''
        def set_text(self, val): self.text = val
        def get_text(self): return self.text

    class _Plugin:
        _current_calib_dirs = qlview.QLView._current_calib_dirs
        _get_tree_base_dir = qlview.QLView._get_tree_base_dir
        _set_cal_status = qlview.QLView._set_cal_status
        _on_cal_not_found = qlview.QLView._on_cal_not_found
        logger = _LOGGER
        def __init__(self, selection):
            self.state = QLViewState()
            self.state.reduced_filepath = selection
            self.spectrograph = load_spectrograph('keck_mosfire')
            self.reduced_text_entry, self.cal_status_label = _Widget(), _Widget()
            self._showing_calib_hint = False
            self.fv = type('FV', (), {'gui_do': staticmethod(lambda fn, *a: fn(*a))})()

    for selection in [f'{setup_dir.parent}/*', str(setup_dir), f'{setup_dir}/*', str(calib_dir),
                      str(calib_dir / 'Slits_A_0_MSC01.fits')]:
        roots.clear()
        qlview.QLView._suggest_calibrations(_Plugin(selection), 'raw.fits')
        assert roots == [str(setup_dir.parent)], selection
