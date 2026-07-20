"""
Tests for the PypeIt Dashboard (Stage 0 — walking skeleton).

They are split into:

* Headless tests of the Qt-free model layer (``read_header_info``,
  ``HeaderInfo``), which run with plain ``pytest`` and need no display.
* Structural widget tests of the views, which use ``pytest-qt`` with the
  offscreen Qt platform (mirroring the ``setup_gui`` test approach).  They
  assert the *structure* of the window (the tab bar, the header fields),
  not its pixel-level appearance.

All are CI-safe: they require no RAW_DATA and use a minimal bundled
``.pypeit`` fixture.
"""

import os
# Ensure Qt can initialize without a display before any Qt import.
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

from pathlib import Path

import pytest

# qtpy/pyqt6 are dependencies, but skip cleanly if the GUI stack is absent.
pytest.importorskip('qtpy')

from pypeit.dashboard.model import read_header_info, HeaderInfo, \
    PYPELINE_DISPLAY


def data_path(filename):
    """
    Return the path to a file in the test ``files`` directory.

    Args:
        filename (:obj:`str`): The file basename.

    Returns:
        str: The absolute path to the test file.
    """
    return str(Path(__file__).parent / 'files' / filename)


# -----------------------------------------------------------------------
# Headless model-layer tests (no Qt, no display)
# -----------------------------------------------------------------------

def test_read_header_info_fields():
    """
    The header metadata is parsed correctly from a minimal ``.pypeit``
    file (R6 fields).
    """
    pyfile = data_path('dashboard_shane_kast_blue.pypeit')
    info = read_header_info(pyfile)
    assert isinstance(info, HeaderInfo)
    assert info.pypeit_file == 'dashboard_shane_kast_blue.pypeit'
    assert info.spectrograph == 'shane_kast_blue'
    assert info.setup == 'A'
    # shane_kast_blue is a long-slit (MultiSlit) spectrograph.
    assert info.path == 'MultiSlit'
    # redux_dir defaults to the directory containing the .pypeit file.
    assert info.redux_dir == str(Path(pyfile).parent.resolve())


def test_read_header_info_redux_path_override(tmp_path):
    """
    An explicit ``redux_path`` overrides the default directory.
    """
    pyfile = data_path('dashboard_shane_kast_blue.pypeit')
    info = read_header_info(pyfile, redux_path=str(tmp_path))
    assert info.redux_dir == str(tmp_path.resolve())


def test_read_header_info_missing_file():
    """
    A named-but-absent ``.pypeit`` file raises a clear error (R11 edge
    case).
    """
    with pytest.raises(FileNotFoundError):
        read_header_info('/no/such/file.pypeit')


def test_pypeline_display_maps_ifu():
    """
    The internal 'SlicerIFU' pypeline is displayed as 'IFU' (S0-Q7).
    """
    assert PYPELINE_DISPLAY['SlicerIFU'] == 'IFU'
    assert PYPELINE_DISPLAY['MultiSlit'] == 'MultiSlit'
    assert PYPELINE_DISPLAY['Echelle'] == 'Echelle'


# -----------------------------------------------------------------------
# Structural widget tests (pytest-qt, offscreen)
# -----------------------------------------------------------------------

def _example_header_info():
    """
    Build a representative :class:`HeaderInfo` for the widget tests.

    Returns:
        :class:`HeaderInfo`: Example header metadata.
    """
    return read_header_info(data_path('dashboard_shane_kast_blue.pypeit'))


def _example_model():
    """
    Build a `DashboardModel` from the minimal ``.pypeit`` fixture (no state
    file → "not started"), for the Stage 0 window/header tests.

    Returns:
        :class:`~pypeit.dashboard.model.DashboardModel`: The example model.
    """
    return dash_model.DashboardModel(
        data_path('dashboard_shane_kast_blue.pypeit'), derive=False)


def test_main_window_tab_bar(qtbot):
    """
    The main window has the Status | Calibrations | Science tab bar
    (R3/C1/R15), in that order.
    """
    from pypeit.dashboard.view.main_window import (DashboardMainWindow,
                                                   TAB_LABELS)
    window = DashboardMainWindow(_example_model())
    qtbot.addWidget(window)
    tabs = [window.tab_widget.tabText(i)
            for i in range(window.tab_widget.count())]
    assert tabs == list(TAB_LABELS)
    assert tabs == ['Status', 'Calibrations', 'Science']


def test_header_banner_fields(qtbot):
    """
    The header banner exposes all five R6 fields, populated from the
    metadata.
    """
    from pypeit.dashboard.view.main_window import DashboardMainWindow
    window = DashboardMainWindow(_example_model())
    qtbot.addWidget(window)
    labels = window.header.value_labels
    assert set(labels.keys()) == {'pypeit_file', 'spectrograph', 'setup',
                                  'path', 'redux_dir'}
    assert labels['spectrograph'].text() == 'shane_kast_blue'
    assert labels['setup'].text() == 'A'
    assert labels['path'].text() == 'MultiSlit'


def test_main_window_title_and_render(qtbot, tmp_path):
    """
    The window has the expected title and renders offscreen without error
    (a basic smoke test of the walking skeleton).
    """
    from pypeit.dashboard.view.main_window import DashboardMainWindow
    window = DashboardMainWindow(_example_model())
    qtbot.addWidget(window)
    window.resize(1650, 900)
    assert window.windowTitle() == 'PypeIt Dashboard'
    # grab() exercises the full paint path; saving proves it produced pixels.
    out = tmp_path / 'render.png'
    assert window.grab().save(str(out))
    assert out.is_file() and out.stat().st_size > 0


def test_main_window_failed_run_marks_build_channel(qtbot):
    """
    **Bug report 000 regression.** A (Re)Build run that exits non-zero must
    leave a "failed" message on the Build channel — and the still-active live
    monitor (the ``.log`` mtime stays "recent" for seconds after a crash) must
    not overwrite it back to "Monitoring…"/"Idle".  A new run starting clears
    the sticky failure.
    """
    from pypeit.dashboard.view.main_window import DashboardMainWindow
    window = DashboardMainWindow(_example_model())
    qtbot.addWidget(window)
    # A failed run (non-zero exit) flags the Build channel as failed.
    window._on_run_finished(1)
    assert window._last_run_failed
    assert 'fail' in window.activity.build_message().lower()
    # A subsequent live-monitor tick must NOT clobber the failure message.
    window._on_state_changed()
    assert 'fail' in window.activity.build_message().lower()
    # Nor does the lock releasing reset it to "Idle".
    window._on_lock_changed(False)
    assert 'fail' in window.activity.build_message().lower()
    # A new run starting clears the sticky failure and resumes monitoring.
    window._on_lock_changed(True)
    assert not window._last_run_failed
    assert 'Monitoring' in window.activity.build_message()
    # A successful run does not flag a failure.
    window._on_run_finished(0)
    assert not window._last_run_failed


# -----------------------------------------------------------------------
# Stage 1: headless state-data-layer tests (DashboardModel, no Qt/RAW_DATA)
# -----------------------------------------------------------------------

from pypeit.dashboard import model as dash_model
from pypeit.dashboard import palette


def _make_redux(tmp_path, case):
    """
    Stage a reduction directory: the minimal ``.pypeit`` fixture as
    ``shane_kast_blue_A.pypeit`` plus, when ``case`` is given, the matching
    state fixture as ``shane_kast_blue_A_state.json`` (the name the model
    derives).  Returns the ``.pypeit`` path.

    Args:
        tmp_path (`pathlib.Path`_): Pytest temporary directory.
        case (:obj:`str`, optional): Fixture case name (e.g. ``healthy``),
            or ``None`` to stage no state file.

    Returns:
        str: Path to the staged ``.pypeit`` file.
    """
    pf = tmp_path / 'shane_kast_blue_A.pypeit'
    pf.write_text(
        Path(data_path('dashboard_shane_kast_blue.pypeit')).read_text())
    if case is not None:
        sj = tmp_path / 'shane_kast_blue_A_state.json'
        sj.write_text(
            Path(data_path(f'dashboard_state_{case}.json')).read_text())
    return str(pf)


def test_model_load_healthy(tmp_path):
    """
    A healthy state loads from the state file and normalizes correctly.
    """
    m = dash_model.DashboardModel(_make_redux(tmp_path, 'healthy'),
                                  derive=False)
    assert m.load_status == dash_model.LOAD_STATE_FILE
    assert m.is_started()
    table = m.status_table()
    assert table.colnames == dash_model.STATUS_COLUMNS
    # A successful, required, in-pipeline step.
    wv = table[table['step'] == 'wv_calib'][0]
    assert wv['status'] == 'success'
    assert bool(wv['required']) is True    # normalized to a real bool
    assert bool(wv['in_pipeline']) is True
    assert wv['output_file'] == 'WaveCalib_A_0_DET01.fits'
    # 'align' is not in the MultiSlit pipeline -> in_pipeline False, absent.
    align = table[table['step'] == 'align'][0]
    assert not bool(align['in_pipeline'])
    assert align['status'] == 'absent'
    assert align['required'] is None
    # 'bpm' is never tracked in the state, so it is not a row.
    assert 'bpm' not in set(table['step'])


def test_model_step_order_path_aware(tmp_path):
    """
    The step order matches the MultiSlit pipeline and omits ``bpm`` by
    default.
    """
    m = dash_model.DashboardModel(_make_redux(tmp_path, 'healthy'),
                                  derive=False)
    order = m.step_order()
    assert order == ['bias', 'dark', 'slits', 'arc', 'tiltimg', 'wv_calib',
                     'tilts', 'scattlight', 'flats']
    assert 'bpm' not in order
    assert 'bpm' in m.step_order(include_bpm=True)


def test_model_calib_det_pairs(tmp_path):
    """
    The (calib_id, det) enumeration returns the single kast pair, det raw.
    """
    m = dash_model.DashboardModel(_make_redux(tmp_path, 'healthy'),
                                  derive=False)
    assert m.calib_det_pairs() == [(0, 1)]


def test_model_failed_step(tmp_path):
    """
    A failed step is normalized as ``status == 'fail'``.
    """
    m = dash_model.DashboardModel(_make_redux(tmp_path, 'failed'),
                                  derive=False)
    wv = m.status_table()
    row = wv[wv['step'] == 'wv_calib'][0]
    assert row['status'] == 'fail'


def test_model_not_started_with_empty_state(tmp_path):
    """
    A state file with no calibration entries reads as "not started"
    (loaded, but empty).
    """
    m = dash_model.DashboardModel(_make_redux(tmp_path, 'not_started'),
                                  derive=False)
    assert m.load_status == dash_model.LOAD_STATE_FILE
    assert not m.is_started()
    assert len(m.status_table()) == 0
    assert m.calib_det_pairs() == []


def test_model_not_started_no_state_file(tmp_path):
    """
    No state file and deriving disabled -> ``not_started`` (R11), no crash.
    """
    m = dash_model.DashboardModel(_make_redux(tmp_path, None), derive=False)
    assert m.load_status == dash_model.LOAD_NOT_STARTED
    assert not m.is_started()
    assert len(m.status_table()) == 0


def test_model_malformed_state(tmp_path):
    """
    A malformed state file is reported as ``malformed`` without raising
    (R11).
    """
    m = dash_model.DashboardModel(_make_redux(tmp_path, 'malformed'),
                                  derive=False)
    assert m.load_status == dash_model.LOAD_MALFORMED
    assert m.run_state is None
    assert len(m.status_table()) == 0


def test_model_file_not_found(tmp_path):
    """
    A missing ``.pypeit`` file is reported, not raised (R11).
    """
    m = dash_model.DashboardModel(str(tmp_path / 'nope.pypeit'),
                                  derive=False)
    assert m.load_status == dash_model.LOAD_FILE_NOT_FOUND
    assert m.header_info is None
    assert len(m.status_table()) == 0


# --- palette (Dashboard-wide status palette) ---------------------------

def test_palette_categories():
    """
    ``classify`` maps (required, status, in_pipeline) to the design's
    palette categories.
    """
    assert palette.classify(True, 'success', True) == palette.SUCCESS
    assert palette.classify(True, 'complete', True) == palette.SUCCESS
    assert palette.classify(True, 'running', True) == palette.RUNNING
    assert palette.classify(True, 'fail', True) == palette.FAIL
    assert palette.classify(True, 'undone', True) == palette.REQUIRED_UNDONE
    assert palette.classify(False, 'undone', True) == palette.OPTIONAL
    # Not in the spectrograph's pipeline -> dimmed, regardless of status.
    assert palette.classify(True, 'success', False) == palette.NOT_USED


def test_palette_colors_and_glyphs_match_design():
    """
    The light palette hex + glyph match the design doc table exactly.
    """
    expected = {
        (True, 'success', True): ('#2E7D32', '✓'),
        (True, 'running', True): ('#EF6C00', '⏳'),
        (True, 'fail', True): ('#C62828', '✗'),
        (True, 'undone', True): ('#FFFFFF', '○'),
        (False, 'undone', True): ('#9E9E9E', '–'),
        (True, 'success', False): ('#BDBDBD', '–'),
    }
    for (req, status, in_pipe), (hexcol, glyph) in expected.items():
        style = palette.step_style(req, status, in_pipe, theme='light')
        assert style.color == hexcol
        assert style.glyph == glyph
    # Dark theme returns the dark hex set (different fills, same glyph).
    dark = palette.step_style(True, 'success', True, theme='dark')
    assert dark.color == palette.DARK_COLORS[palette.SUCCESS]
    assert dark.glyph == '✓'


def test_palette_worst_category():
    """
    ``worst_category`` picks the most severe category; optional/not-used
    never worsen an otherwise-successful cell (R17).
    """
    assert palette.worst_category(
        [palette.SUCCESS, palette.FAIL]) == palette.FAIL
    assert palette.worst_category(
        [palette.SUCCESS, palette.RUNNING]) == palette.RUNNING
    assert palette.worst_category(
        [palette.SUCCESS, palette.REQUIRED_UNDONE]) \
        == palette.REQUIRED_UNDONE
    # optional / not_used do not worsen a successful cell.
    assert palette.worst_category(
        [palette.SUCCESS, palette.OPTIONAL, palette.NOT_USED]) \
        == palette.SUCCESS
    assert palette.worst_category([]) == palette.NOT_USED


# -----------------------------------------------------------------------
# Stage 2: Status view (pytest-qt, offscreen) — structural tests
# -----------------------------------------------------------------------

def _status_view(case, tmp_path, qtbot):
    """
    Build a `StatusView` over a staged fixture reduction.

    Args:
        case (:obj:`str`): Fixture case (e.g. ``healthy``), or ``None``.
        tmp_path (`pathlib.Path`_): Pytest temp dir.
        qtbot: pytest-qt fixture.

    Returns:
        tuple: ``(StatusView, DashboardModel)``.
    """
    from pypeit.dashboard.view.status_view import StatusView
    model = dash_model.DashboardModel(_make_redux(tmp_path, case),
                                      derive=False)
    view = StatusView(model)
    qtbot.addWidget(view)
    return view, model


def test_status_view_healthy_table(tmp_path, qtbot):
    """
    The healthy Status view shows the scoped table (Step | Required |
    Status | Output) with one row per pipeline step.
    """
    view, model = _status_view('healthy', tmp_path, qtbot)
    table = view._table
    assert table is not None
    headers = [table.horizontalHeaderItem(c).text()
               for c in range(table.columnCount())]
    assert headers == ['Step', 'Required', 'Status', 'Output']
    assert table.rowCount() == len(model.step_order())
    # The wv_calib row reads as a success (✓).
    steps = [table.item(r, 0).text() for r in range(table.rowCount())]
    wv = steps.index('wv_calib')
    assert '✓' in table.item(wv, 2).text()
    assert table.item(wv, 3).text() == 'WaveCalib_A_0_DET01.fits'


def test_status_view_scope_dropdowns(tmp_path, qtbot):
    """
    The scope drop-downs are populated from the model's (calib, det) pairs.
    """
    view, model = _status_view('healthy', tmp_path, qtbot)
    assert view._group_combo.count() == 1
    assert view._group_combo.itemData(0) == 0
    assert view._det_combo.count() == 1
    # det is formatted readably by the view (DET01 for a single detector).
    assert view._det_combo.itemText(0) == 'DET01'


def test_status_view_failed_row(tmp_path, qtbot):
    """
    A failed step renders with the fail glyph (✗).
    """
    view, _ = _status_view('failed', tmp_path, qtbot)
    table = view._table
    steps = [table.item(r, 0).text() for r in range(table.rowCount())]
    wv = steps.index('wv_calib')
    assert '✗' in table.item(wv, 2).text()


def test_status_view_not_started_message(tmp_path, qtbot):
    """
    A not-started (empty) state renders an edge message and no table (R11).
    """
    view, _ = _status_view('not_started', tmp_path, qtbot)
    assert view._table is None


def test_status_view_malformed_message(tmp_path, qtbot):
    """
    A malformed state renders an edge message and no table (R11).
    """
    view, model = _status_view('malformed', tmp_path, qtbot)
    assert model.load_status == dash_model.LOAD_MALFORMED
    assert view._table is None


# -----------------------------------------------------------------------
# Stage 3: model detail accessors, palette skip, inspect builders
# -----------------------------------------------------------------------

def _multidet_model(tmp_path):
    """
    A `DashboardModel` over the multi-slit/multi-detector fixture.
    """
    return dash_model.DashboardModel(_make_redux(tmp_path, 'multidet'),
                                     derive=False)


def test_model_step_metrics(tmp_path):
    """
    step_metrics returns entry-level metrics per step (Stage 3).
    """
    m = _multidet_model(tmp_path)
    assert m.step_metrics('bias', 0, 1)['mean'] == 1001.0
    assert m.step_metrics('slits', 0, 1)['nslits'] == 3
    fm = m.step_metrics('flats', 0, 1)
    assert fm['corrections'] == ['pixelflat', 'spat_illum']
    assert fm['pixelflat_source'] == 'raw'


def test_model_slit_table(tmp_path):
    """
    slit_table returns per-slit rows; flats carries per-correction metrics
    and a `skip` slit (Stage 3).
    """
    m = _multidet_model(tmp_path)
    wv = m.slit_table('wv_calib', 0, 1)
    assert len(wv) == 3
    assert 'rms' in wv[0] and wv[0]['status'] == 'success'
    flats = m.slit_table('flats', 0, 1)
    assert [r['status'] for r in flats] == ['success', 'success', 'skip']
    assert set(flats[0]['corrections']) == {'pixelflat', 'spat_illum'}
    # Steps without per-slit data return an empty list.
    assert m.slit_table('bias', 0, 1) == []


def test_model_output_and_calib_paths(tmp_path):
    """
    output_path joins Calibrations/<output_file>; calib_file_path builds a
    file by naming convention (Edges for slits).
    """
    m = _multidet_model(tmp_path)
    out = m.output_path('wv_calib', 0, 1)
    assert out.name == 'WaveCalib_A_0_DET01.fits'
    assert out.parent.name == 'Calibrations'
    edges = m.calib_file_path('Edges', 0, 1, ext='fits.gz')
    assert edges.name == 'Edges_A_0_DET01.fits.gz'


def test_palette_slit_style_skip():
    """
    The per-slit palette maps `skip` to a distinct category with its glyph
    (S3-Q15).
    """
    assert palette.slit_style('skip').category == palette.SKIP
    assert palette.slit_style('skip').glyph == '⊘'
    assert palette.slit_style('success').category == palette.SUCCESS
    assert palette.slit_style('fail').category == palette.FAIL
    assert palette.SKIP in palette.LIGHT_COLORS
    assert palette.SKIP in palette.DARK_COLORS


def test_inspect_output_commands(tmp_path):
    """
    inspect.output_command builds the right argv per step, incl. the
    non-uniform cases (S3-Q12/Q13).
    """
    from pypeit.dashboard import inspect as dash_inspect
    m = _multidet_model(tmp_path)
    # Processed calibration image → opened directly in ginga (Stage 4 R2 #1).
    bias = dash_inspect.output_command(m, 'bias', 0, 1)
    assert bias[0] == 'ginga'
    assert bias[1].endswith('Bias_A_0_DET01.fits')
    arc = dash_inspect.output_command(m, 'arc', 0, 1)
    assert arc[0] == 'ginga' and arc[1].endswith('Arc_A_0_DET01.fits')
    # slits → pypeit_chk_edges on the Edges file, not the Slits output.
    slits = dash_inspect.output_command(m, 'slits', 0, 1)
    assert slits[0] == 'pypeit_chk_edges'
    assert slits[1].endswith('Edges_A_0_DET01.fits.gz')
    # wv_calib has no standalone viewer (Round-2 #4) → command is None.
    assert dash_inspect.output_command(m, 'wv_calib', 0, 1) is None
    # An input frame → pypeit_view_fits <spec> <file> --proc (+ --det).
    vin = dash_inspect.view_input_command(m, 'b3.fits.gz', det=1)
    assert vin == ['pypeit_view_fits', 'shane_kast_blue', 'b3.fits.gz',
                   '--proc', '--det', '1']


# -----------------------------------------------------------------------
# Stage 3: Calibrations view (pytest-qt, offscreen) — structural tests
# -----------------------------------------------------------------------

def _calib_view(case, tmp_path, qtbot):
    """
    Build a `CalibrationsView` over a staged fixture reduction.
    """
    from pypeit.dashboard.view.calibrations_view import CalibrationsView
    model = dash_model.DashboardModel(_make_redux(tmp_path, case),
                                      derive=False)
    view = CalibrationsView(model)
    qtbot.addWidget(view)
    return view, model


def test_calib_view_button_row_path_aware(tmp_path, qtbot):
    """
    The step-button row contains the spectrograph's steps in order, with
    `bpm` omitted (C3/C12/C13).
    """
    view, model = _calib_view('multidet', tmp_path, qtbot)
    assert list(view._step_buttons.keys()) == model.step_order()
    assert 'bpm' not in view._step_buttons


def test_calib_view_detector_selector_multidet(tmp_path, qtbot):
    """
    The detector selector lists both detectors (C2, multi-detector).
    """
    view, _ = _calib_view('multidet', tmp_path, qtbot)
    assert view._group_combo.count() == 1
    assert view._det_combo.count() == 2


def test_calib_view_per_slit_table(tmp_path, qtbot):
    """
    Selecting `wv_calib` shows a per-slit table (one row per slit), and
    `flats` shows per-correction columns incl. a skipped slit (C11, S3-Q14).
    """
    from qtpy.QtWidgets import QTableWidget
    view, _ = _calib_view('multidet', tmp_path, qtbot)
    view._select_step('wv_calib')
    tables = view.findChildren(QTableWidget)
    assert tables and tables[-1].rowCount() == 3   # 3 slits
    # Flats: per-correction columns (Slit, Status, then mean/rms per corr).
    view._select_step('flats')
    tbl = view.findChildren(QTableWidget)[-1]
    headers = [tbl.horizontalHeaderItem(c).text()
               for c in range(tbl.columnCount())]
    assert headers[:2] == ['Slit', 'Status']
    assert any('pixelflat' in h for h in headers)
    # The third slit is skipped.
    statuses = [tbl.item(r, 1).text() for r in range(tbl.rowCount())]
    assert any('skipped' in s for s in statuses)


def test_calib_view_not_started_message(tmp_path, qtbot):
    """
    A not-started state shows a message and builds no step buttons (R11).
    """
    view, _ = _calib_view('not_started', tmp_path, qtbot)
    assert view._step_buttons == {}


# -----------------------------------------------------------------------
# Stage 3b: echelle (keck_nires) coverage
# -----------------------------------------------------------------------

def _make_echelle_redux(tmp_path):
    """
    Stage an echelle reduction: the minimal ``keck_nires`` ``.pypeit`` fixture
    plus the synthesized echelle state, named so the model finds the state.

    Args:
        tmp_path (`pathlib.Path`_): Pytest temporary directory.

    Returns:
        str: Path to the staged ``.pypeit`` file.
    """
    pf = tmp_path / 'keck_nires_A.pypeit'
    pf.write_text(Path(data_path('dashboard_keck_nires.pypeit')).read_text())
    (tmp_path / 'keck_nires_A_state.json').write_text(
        Path(data_path('dashboard_state_echelle.json')).read_text())
    return str(pf)


def test_echelle_pipeline_and_orders(tmp_path):
    """
    The echelle reduction loads as the Echelle pipeline and exposes per-order
    rows.
    """
    m = dash_model.DashboardModel(_make_echelle_redux(tmp_path), derive=False)
    assert m.header_info.path == 'Echelle'
    assert m.calib_det_pairs() == [(0, 1)]
    assert len(m.slit_table('wv_calib', 0, 1)) == 5


def test_echelle_slits_edges_all_naming(tmp_path):
    """
    For shared (``*_all_*``) echelle calibrations, the slits viewer resolves
    the Edges file from the recorded Slits output (Stage 3b fix).
    """
    from pypeit.dashboard import inspect as dash_inspect
    m = dash_model.DashboardModel(_make_echelle_redux(tmp_path), derive=False)
    target = dash_inspect.output_target(m, 'slits', 0, 1)
    assert target.name == 'Edges_A_all_DET01.fits.gz'
    cmd = dash_inspect.output_command(m, 'slits', 0, 1)
    assert cmd[0] == 'pypeit_chk_edges'
    assert cmd[1].endswith('Edges_A_all_DET01.fits.gz')


def test_echelle_per_order_label(tmp_path, qtbot):
    """
    The per-order table is labeled "Order" (not "Slit") for the Echelle
    pipeline (S3b-Q1).
    """
    from pypeit.dashboard.view.calibrations_view import CalibrationsView
    m = dash_model.DashboardModel(_make_echelle_redux(tmp_path), derive=False)
    view = CalibrationsView(m)
    qtbot.addWidget(view)
    view._select_step('wv_calib')
    from qtpy.QtWidgets import QTableWidget
    table = view.findChildren(QTableWidget)[-1]
    assert table.horizontalHeaderItem(0).text() == 'Order'
    assert table.rowCount() == 5


# -----------------------------------------------------------------------
# Stage 4: Execution, locking & (Re)Build (C10, C15, X1–X3)
# -----------------------------------------------------------------------

def test_inspect_run_command(tmp_path):
    """
    run_command builds the pypeit_run_to_calibstep argv with --calib_group,
    --det, and --redux_path; mosaics become a parenthesized tuple (C10,
    S4-Q3).
    """
    from pypeit.dashboard import inspect as dash_inspect
    m = _multidet_model(tmp_path)
    argv = dash_inspect.run_command(m, 'wv_calib', 0, 1)
    assert argv[0] == 'pypeit_run_to_calibstep'
    assert argv[2] == 'wv_calib'
    assert argv[argv.index('--calib_group') + 1] == '0'
    assert argv[argv.index('--det') + 1] == '1'
    assert '--redux_path' in argv
    # Detector-argument formatting for eval_detectors: int vs mosaic tuple.
    assert dash_inspect._det_run_arg(3) == '3'
    assert dash_inspect._det_run_arg((1, 2)) == '(1,2)'


def test_model_step_output_files(tmp_path):
    """
    step_output_files names the existing Calibrations/ outputs a (re)build
    overwrites — both Slits_* and Edges_* for slits, the single output
    otherwise; empty when nothing is on disk (S4-Q4).
    """
    m = _multidet_model(tmp_path)
    # Nothing on disk yet → a fresh build, nothing to overwrite.
    assert m.step_output_files('slits', 0, 1) == []
    m.calib_dir.mkdir(parents=True, exist_ok=True)
    slits = m.output_path('slits', 0, 1)
    edges = slits.with_name(slits.name.replace('Slits', 'Edges', 1))
    slits.write_text('x')
    edges.write_text('x')
    names = {p.name for p in m.step_output_files('slits', 0, 1)}
    assert names == {slits.name, edges.name}
    # A non-slits step returns just its single output when present.
    wv = m.output_path('wv_calib', 0, 1)
    wv.write_text('x')
    assert [p.name for p in m.step_output_files('wv_calib', 0, 1)] \
        == [wv.name]


def test_runlock_state_machine(tmp_path, qtbot):
    """
    RunLock locks for a Dashboard-launched run and for an externally-written
    .log (recent mtime), emitting lockChanged on transitions (X1).
    """
    import os
    from pypeit.dashboard.runlock import RunLock

    # The pure recency test.
    assert RunLock._is_recent(100.0, 105.0, 10.0)
    assert not RunLock._is_recent(100.0, 200.0, 10.0)

    # Dashboard-launched run: locks + emits on each transition.
    lock = RunLock()
    seen = []
    lock.lockChanged.connect(lambda v: seen.append(v))
    assert not lock.is_locked()
    lock.set_dashboard_running(True)
    lock.set_dashboard_running(True)        # no-op (already running)
    lock.set_dashboard_running(False)
    assert not lock.is_locked()
    assert seen == [True, False]

    # External run: a freshly written .log locks; an old one does not.
    log_file = tmp_path / 'run.log'
    ext = RunLock(log_path=log_file)
    ext.poll()
    assert not ext.is_locked()              # no log yet
    log_file.write_text('running')
    ext.poll()
    assert ext.is_locked()                  # recent mtime
    os.utime(log_file, (1.0, 1.0))          # far in the past
    ext.poll()
    assert not ext.is_locked()


def test_calib_view_rebuild_button(tmp_path, qtbot):
    """
    The detail panel shows a "(Re)Build" control (C10) that the run lock
    enables/disables (X1).
    """
    view, _ = _calib_view('multidet', tmp_path, qtbot)
    view._select_step('wv_calib')
    button = view._rebuild_button
    assert button is not None and '(Re)Build' in button.text()
    assert button.isEnabled()               # idle
    view.set_locked(True)
    assert not button.isEnabled()           # locked
    view.set_locked(False)
    assert button.isEnabled()


def test_calib_view_rebuild_disabled_when_locked(tmp_path, qtbot):
    """
    With the run lock already engaged, the (Re)Build control builds disabled
    (X1).
    """
    from pypeit.dashboard.view.calibrations_view import CalibrationsView
    from pypeit.dashboard.runlock import RunLock
    model = dash_model.DashboardModel(_make_redux(tmp_path, 'multidet'),
                                      derive=False)
    lock = RunLock()
    lock.set_dashboard_running(True)
    view = CalibrationsView(model, run_lock=lock)
    qtbot.addWidget(view)
    view._select_step('wv_calib')
    assert not view._rebuild_button.isEnabled()


def test_calib_view_rebuild_clobbers_and_launches(tmp_path, qtbot,
                                                   monkeypatch):
    """
    On confirmation, (Re)Build removes the selected step's existing output(s)
    in code and launches pypeit_run_to_calibstep (C10, X3, S4-Q3).
    """
    view, model = _calib_view('multidet', tmp_path, qtbot)
    # Stage an existing wv_calib output to be clobbered.
    model.calib_dir.mkdir(parents=True, exist_ok=True)
    wv = model.output_path('wv_calib', 0, 1)
    wv.write_text('x')

    captured = {}

    class _FakeLauncher:
        def run(self, argv, description=None, on_finished=None):
            captured['argv'] = argv
            captured['description'] = description

    view._launcher = _FakeLauncher()
    # Auto-confirm the clobber dialog (no modal UI in the test).
    monkeypatch.setattr(view, '_confirm_rebuild', lambda *a, **k: True)
    view._select_step('wv_calib')
    view._on_rebuild('wv_calib', 0, 1)
    # Moved aside (crash-safe clobber): the original is gone, a .dashboard_bak
    # holds it until the run finishes (no shell rm).
    assert not wv.exists()
    assert (wv.parent / (wv.name + '.dashboard_bak')).exists()
    assert captured['argv'][0] == 'pypeit_run_to_calibstep'
    assert 'wv_calib' in captured['description']


def test_calib_view_rebuild_restores_on_failure(tmp_path, qtbot):
    """
    A failed (re)build restores the moved-aside output so the calibration is
    not lost (Round-1 #2).
    """
    view, model = _calib_view('multidet', tmp_path, qtbot)
    model.calib_dir.mkdir(parents=True, exist_ok=True)
    wv = model.output_path('wv_calib', 0, 1)
    bak = wv.parent / (wv.name + '.dashboard_bak')
    bak.write_text('orig')                  # the moved-aside original
    view._after_rebuild(1, 'wv_calib', [(bak, wv)])   # non-zero → restore
    assert wv.exists() and not bak.exists()
    assert wv.read_text() == 'orig'


def test_calib_view_rebuild_drops_backup_on_success(tmp_path, qtbot):
    """
    A successful (re)build drops the backup and keeps the freshly built output
    (Round-1 #2).
    """
    view, model = _calib_view('multidet', tmp_path, qtbot)
    model.calib_dir.mkdir(parents=True, exist_ok=True)
    wv = model.output_path('wv_calib', 0, 1)
    bak = wv.parent / (wv.name + '.dashboard_bak')
    bak.write_text('old')
    wv.write_text('new')                    # the rebuild produced a fresh one
    view._after_rebuild(0, 'wv_calib', [(bak, wv)])   # success → drop backup
    assert wv.exists() and not bak.exists()
    assert wv.read_text() == 'new'


def test_calib_view_rebuild_cancel_keeps_output(tmp_path, qtbot, monkeypatch):
    """
    Cancelling the clobber confirmation leaves the output untouched and
    launches nothing (X2).
    """
    view, model = _calib_view('multidet', tmp_path, qtbot)
    model.calib_dir.mkdir(parents=True, exist_ok=True)
    wv = model.output_path('wv_calib', 0, 1)
    wv.write_text('x')

    launched = {'called': False}

    class _FakeLauncher:
        def run(self, *a, **k):
            launched['called'] = True

    view._launcher = _FakeLauncher()
    monkeypatch.setattr(view, '_confirm_rebuild', lambda *a, **k: False)
    view._select_step('wv_calib')
    view._on_rebuild('wv_calib', 0, 1)
    assert wv.exists()                      # not clobbered
    assert launched['called'] is False


def test_calib_view_output_filename_shown(tmp_path, qtbot):
    """
    The detail panel shows the output filename, even for wv_calib (no viewer)
    (Round-2 #3).
    """
    from qtpy.QtWidgets import QLabel
    view, _ = _calib_view('multidet', tmp_path, qtbot)
    view._select_step('wv_calib')
    labels = [w.text() for w in view.findChildren(QLabel)]
    assert any(t.startswith('Output:') and 'WaveCalib' in t for t in labels)


def test_calib_view_rebuild_locked_visual_cue(tmp_path, qtbot):
    """
    While locked, the (Re)Build button turns orange, shows a "run in progress"
    label, and is disabled (the lock's visual clue); idle restores the blue
    "(Re)Build <step>" enabled button (Round-3 #1/#2).
    """
    view, _ = _calib_view('multidet', tmp_path, qtbot)
    view._select_step('tilts')
    button = view._rebuild_button
    assert button.isEnabled() and button.text() == '(Re)Build tilts'
    view.set_locked(True)
    assert not button.isEnabled()
    assert 'progress' in button.text().lower()
    assert palette.LIGHT_COLORS[palette.RUNNING] in button.styleSheet()
    view.set_locked(False)
    assert button.isEnabled() and button.text() == '(Re)Build tilts'


def test_calib_view_rebuild_keeps_step_selected(tmp_path, qtbot):
    """
    On completion the rebuilt step stays selected, not reset to the first step
    (Round-2 #2).
    """
    view, _ = _calib_view('multidet', tmp_path, qtbot)
    view._select_step('tilts')
    view._after_rebuild(0, 'tilts', [])
    assert view._selected_step == 'tilts'


# -----------------------------------------------------------------------
# Stage 5: Monitoring (live updates) — R14
# -----------------------------------------------------------------------

def test_activity_bar_two_channels(qtbot):
    """
    The ActivityBar has independent Build and Inspection channels (S5-Q9).
    """
    from pypeit.dashboard.view.activity import ActivityBar
    bar = ActivityBar()
    qtbot.addWidget(bar)
    bar.set_build('building', busy=True)
    bar.set_inspection('viewing a frame')
    assert bar.build_message() == 'building'
    assert bar.inspection_message() == 'viewing a frame'
    # Channels are independent: setting one leaves the other.
    bar.set_build('monitoring')
    assert bar.inspection_message() == 'viewing a frame'
    bar.idle()
    assert bar.build_message() == 'Idle' and bar.inspection_message() == '—'


def test_runlock_state_changed_signal(tmp_path, qtbot):
    """
    RunLock emits `stateChanged` when the state file's mtime advances **while a
    run is active**, and not when idle or unchanged (R14, S5-Q3/Q8).
    """
    import os
    from pypeit.dashboard.runlock import RunLock
    state = tmp_path / 'x_state.json'
    state.write_text('{}')
    lock = RunLock(state_path=state)        # no log → lock only via dashboard
    fired = []
    lock.stateChanged.connect(lambda: fired.append(True))

    # Idle: a state change does NOT fire (monitoring only while active).
    os.utime(state, (1000.0, 1000.0))
    lock.poll()
    assert fired == []

    # Active (a Dashboard-launched run) + the state mtime advances → fires once.
    lock.set_dashboard_running(True)
    os.utime(state, (2000.0, 2000.0))
    lock.poll()
    assert fired == [True]
    # No further change → no new fire.
    lock.poll()
    assert fired == [True]


def test_main_window_live_refresh(tmp_path, qtbot):
    """
    While a run is active, a state-file change live-refreshes both views
    (new model) and preserves the selected step (R14).
    """
    import os
    from pypeit.dashboard.view.main_window import DashboardMainWindow
    model = dash_model.DashboardModel(_make_redux(tmp_path, 'multidet'),
                                      derive=False)
    window = DashboardMainWindow(model)
    qtbot.addWidget(window)
    window.tab_widget.setCurrentIndex(1)
    window.calibrations_view._select_step('wv_calib')

    # Make the run "active" (a recent .log) and bump the state-file mtime.
    window.model.log_path.write_text('running')
    far_future = 4102444800.0               # year 2100, > the last-seen mtime
    os.utime(window.model.state_path, (far_future, far_future))

    before = window.model
    window.run_lock.poll()                  # lock engages + stateChanged fires
    # The view live-refreshed to a fresh model, keeping the selected step.
    assert window.model is not before
    assert window.calibrations_view._selected_step == 'wv_calib'


# -----------------------------------------------------------------------
# Stage 6: Science frames (R15/R18)
# -----------------------------------------------------------------------

def _science_model(tmp_path):
    """
    A `DashboardModel` over the synthesized science fixture.
    """
    return dash_model.DashboardModel(_make_redux(tmp_path, 'science'),
                                     derive=False)


def test_model_science_accessors(tmp_path):
    """
    The model exposes the science table + per-slit/per-object detail (Stage 6).
    """
    m = _science_model(tmp_path)
    assert m.has_science()
    sci = m.science_table()
    assert sci.colnames == dash_model.SCIENCE_COLUMNS
    assert len(sci) == 2
    assert set(sci['objtype']) == {'science', 'standard'}
    # Per-slit + per-object detail for the science frame.
    slits = m.science_slit_table('b27-J1217p3905_KASTb', 1)
    assert slits and slits[0]['status'] == 'success' and slits[0]['nobj'] == 1
    objs = m.science_object_table('b27-J1217p3905_KASTb', 1)
    assert objs and objs[0]['snr_find'] == 142.5 and objs[0]['extracted']


def test_inspect_science_commands(tmp_path):
    """
    inspect builds the spec2d/spec1d viewers (per-object --obj) and the
    pypeit_reduce_by_step (Re)Build argv (Stage 6, S6-Q10).
    """
    from pypeit.dashboard import inspect as dash_inspect
    m = _science_model(tmp_path)
    # 2D viewer.
    s2d = dash_inspect.spec2d_command('spec2d_b27.fits', det=1)
    assert s2d[0] == 'pypeit_show_2dspec' and s2d[-1] == '1'
    # Object-name reconstruction + per-object 1D viewer.
    name = dash_inspect.science_object_name(424.0, 175, 'DET01')
    assert name == 'SPAT0424-SLIT0175-DET01'
    s1d = dash_inspect.spec1d_command('spec1d_b27.fits', obj_name=name)
    assert s1d[0] == 'pypeit_show_1dspec'
    assert s1d[s1d.index('--obj') + 1] == name
    # (Re)Build via reduce_by_step uses the recorded raw frame.
    run = dash_inspect.science_run_command(m, 'b27-J1217p3905_KASTb', 1,
                                           'findobj')
    assert run[0] == 'pypeit_reduce_by_step'
    assert 'b27.fits.gz' in run and run[3] == 'findobj'
    assert run[run.index('--det') + 1] == '1'


def _science_view(tmp_path, qtbot, run_lock=None):
    """
    Build a `ScienceView` over the science fixture.
    """
    from pypeit.dashboard.view.science_view import ScienceView
    model = _science_model(tmp_path)
    view = ScienceView(model, run_lock=run_lock)
    qtbot.addWidget(view)
    return view, model


def test_science_view_table_and_detail(tmp_path, qtbot):
    """
    The Science view lists per-frame rows and, on selection, shows per-slit +
    per-object detail (R15/R18).
    """
    from qtpy.QtWidgets import QTableWidget
    view, _ = _science_view(tmp_path, qtbot)
    assert len(view._rows) == 2                      # two frames
    headers = [view._table.horizontalHeaderItem(c).text()
               for c in range(view._table.columnCount())]
    assert headers[:4] == ['Frame', 'Calib', 'Detector', 'Type']
    assert 'extract' in headers
    # First row auto-selected → detail panel has per-slit + per-object tables.
    view._select_pair(view._rows[0])
    tables = view.findChildren(QTableWidget)
    assert len(tables) >= 3       # the frame table + per-slit + per-object


def test_science_view_rebuild_buttons(tmp_path, qtbot):
    """
    The Science detail offers (Re)Build controls, prerequisite-gated and
    lock-aware (S6-Q6/Q12).
    """
    from pypeit.dashboard.runlock import RunLock
    lock = RunLock()
    view, _ = _science_view(tmp_path, qtbot, run_lock=lock)
    view._select_pair(view._rows[0])
    # All three steps available (prereqs succeeded) and enabled when idle.
    assert set(view._rebuild_buttons) == {'process', 'findobj', 'extract'}
    btn, avail = view._rebuild_buttons['extract']
    assert avail and btn.isEnabled() and btn.text() == '(Re)Build extract'
    # Locking turns them orange + disabled.
    view.set_locked(True)
    btn, _ = view._rebuild_buttons['extract']
    assert not btn.isEnabled() and 'progress' in btn.text().lower()


def test_main_window_has_science_tab(tmp_path, qtbot):
    """
    The window now has a third Science tab (S6-Q1).
    """
    from pypeit.dashboard.view.main_window import (DashboardMainWindow,
                                                   TAB_LABELS)
    assert TAB_LABELS == ('Status', 'Calibrations', 'Science')
    model = _science_model(tmp_path)
    window = DashboardMainWindow(model)
    qtbot.addWidget(window)
    assert window.tab_widget.count() == 3
    assert window.tab_widget.tabText(2) == 'Science'


# -- Stage 6 Round-1 Modifications ---------------------------------------

def test_palette_selection_style_neutral():
    """
    The neutral selection style names a soft blue-grey and targets the
    selected-item rule (Round-1 #2).
    """
    light = palette.selection_style('light')
    dark = palette.selection_style('dark')
    assert 'item:selected' in light
    assert palette.SELECTION_COLORS['light'].lower() in light.lower()
    assert palette.SELECTION_COLORS['dark'].lower() in dark.lower()


def test_science_view_neutral_selection(tmp_path, qtbot):
    """
    The Science per-frame table carries the neutral selection stylesheet, so a
    selected frame never inherits the theme's (red) highlight (Round-1 #2).
    """
    view, _ = _science_view(tmp_path, qtbot)
    ss = view._table.styleSheet()
    assert 'item:selected' in ss
    assert (palette.SELECTION_COLORS['light'].lower() in ss.lower()
            or palette.SELECTION_COLORS['dark'].lower() in ss.lower())


def test_status_view_science_navigator(tmp_path, qtbot):
    """
    The Status view shows a clickable science-navigator cell per (frame, det);
    clicking one emits scienceFrameActivated (Round-1 #1).
    """
    from pypeit.dashboard.view.status_view import StatusView, ScienceNavCell
    model = _science_model(tmp_path)
    view = StatusView(model)
    qtbot.addWidget(view)
    cells = view.findChildren(ScienceNavCell)
    assert len(cells) == 2                        # one per (frame, det)
    captured = []
    view.scienceFrameActivated.connect(
        lambda f, d: captured.append((f, d)))
    cells[0].clicked.emit(cells[0]._frame, cells[0]._det)
    assert captured and captured[0][0] in {'b27-J1217p3905_KASTb',
                                           'b24-Feige66_KASTb'}


def test_status_navigator_switches_to_science_tab(tmp_path, qtbot):
    """
    Activating a science-navigator cell switches to the Science tab and selects
    that frame (Round-1 #1).
    """
    from pypeit.dashboard.view.main_window import DashboardMainWindow
    from pypeit.dashboard.view.status_view import ScienceNavCell
    model = _science_model(tmp_path)
    window = DashboardMainWindow(model)
    qtbot.addWidget(window)
    window.tab_widget.setCurrentIndex(0)
    cells = window.status_view.findChildren(ScienceNavCell)
    assert cells
    cells[0].clicked.emit(cells[0]._frame, cells[0]._det)
    assert window.tab_widget.currentWidget() is window.science_view


def _stage_science_qa(model):
    """
    Create a few science QA PNGs for the b27 frame on DET01.  Returns the
    frame basename.
    """
    detname = model.det_name(1)
    qa = model.redux_dir / 'QA' / 'PNGs'
    qa.mkdir(parents=True, exist_ok=True)
    frame = 'b27-J1217p3905_KASTb'
    names = [f'{frame}_{detname}_S0175_obj_prof.png',
             f'{frame}_{detname}_S0175_obj_trace.png',
             f'{frame}_global_{detname}_global_{detname}_S0175'
             f'_spec_flex_corr.png',
             f'b24-Feige66_KASTb_{detname}_S0175_obj_prof.png']  # other frame
    for n in names:
        (qa / n).write_bytes(b'')
    return frame


def test_model_science_qa_files(tmp_path):
    """
    The model globs the per-(frame, det) QA PNGs (all of them) and maps the
    per-object obj_prof/obj_trace (Round-1 #3, S6-Q15).
    """
    m = _science_model(tmp_path)
    frame = _stage_science_qa(m)
    files = m.science_qa_files(frame, 1)
    assert len(files) == 3                        # b27's 3, not the b24 one
    assert all(frame in f.name for f in files)
    assert any('spec_flex' in f.name for f in files)   # flexure included
    objqa = m.science_object_qa_files(frame, 1, 175)
    assert 'obj_prof' in objqa['obj_prof'].name
    assert 'obj_trace' in objqa['obj_trace'].name


def test_science_view_qa_list_and_object_qa(tmp_path, qtbot):
    """
    The Science detail panel exposes a QA list (all PNGs) and per-object
    obj_prof/obj_trace columns (Round-1 #3).
    """
    from qtpy.QtWidgets import QListWidget
    from pypeit.dashboard.view.science_view import ScienceView
    model = _science_model(tmp_path)
    frame = _stage_science_qa(model)
    view = ScienceView(model)
    qtbot.addWidget(view)
    view._select_pair((frame, 1))
    lists = view.findChildren(QListWidget)
    assert lists and lists[0].count() == 3        # all 3 QA PNGs listed
    assert 'obj_prof' in view._obj_qa_cols
    assert 'obj_trace' in view._obj_qa_cols


# -- Stage 6 Round-3 (planned science frames + calib-gated (Re)Build) -----

def test_model_calibrations_ready(tmp_path):
    """
    calibrations_ready() is True only when a group/det's required calibs are
    all successful (Round-3 #2).
    """
    m = _science_model(tmp_path)                  # healthy calibs
    assert m.calibrations_ready(0, 1) is True
    # An unknown group has no calib entries → not ready.
    assert m.calibrations_ready(99, 1) is False
    # A failed-calibration reduction is not ready.
    mf = dash_model.DashboardModel(_make_redux(tmp_path, 'failed'),
                                   derive=False)
    pairs = mf.calib_det_pairs()
    if pairs:
        g, d = pairs[0]
        assert mf.calibrations_ready(g, d) is False


def test_science_rebuild_gated_on_calibrations(tmp_path, qtbot):
    """
    The science (Re)Build is disabled when the frame's calibrations are not
    built (Round-3 #2).
    """
    from pypeit.dashboard.view.science_view import ScienceView
    model = _science_model(tmp_path)
    view = ScienceView(model)
    qtbot.addWidget(view)
    view._select_pair(view._rows[0])
    # Healthy calibs → process (Re)Build available.
    assert view._rebuild_buttons['process'][1] is True
    # Monkeypatch calibrations to "not ready" and rebuild the detail: the
    # (Re)Build must now be unavailable.
    model.calibrations_ready = lambda *a, **k: False
    view._build_detail(*view._rows[0])
    assert view._rebuild_buttons['process'][1] is False


def test_seed_planned_science():
    """
    seed_planned_science() creates an undone entry per (frame, det) from the
    metadata, with raw files, for science and standard frames (Round-3 #2).
    """
    from pypeit.state.run_state import RunPypeItState
    from pypeit.state import science_status

    class _FakeFitstbl:
        """Minimal stand-in for PypeItMetaData for the seeder."""
        _names = ['b27-sci', 'b24-std']
        _files = ['b27.fits.gz', 'b24.fits.gz']

        def find_frames(self, ftype):
            if ftype == 'science':
                return [True, False]
            if ftype == 'standard':
                return [False, True]
            return [False, False]

        def construct_basename(self, i):
            return self._names[i]

        def find_frame_calib_groups(self, i):
            return [0]

        def keys(self):
            return ['filename', 'comb_id']

        def __len__(self):
            return 2

        def __getitem__(self, key):
            return {'filename': self._files, 'comb_id': [-1, -1]}[key]

    rs = RunPypeItState(pypeit_file='x.pypeit', current_step='init',
                        current_det=-1, current_calibID=-1)
    science_status.seed_planned_science(rs, _FakeFitstbl(), {0: [1]})
    tbl = rs.get_science_status()
    assert len(tbl) == 2
    assert set(tbl['objtype']) == {'science', 'standard'}
    # Planned frames are all 'undone' with the raw frame recorded.
    assert set(tbl['process']) == {'undone'}
    sci = rs.science_entry('b27-sci', 1)
    assert sci.raw_files == ['b27.fits.gz'] and sci.calib_id == 0


def test_planned_science_persists_on_state_load(tmp_path):
    """
    Round-4: loading a state file that has calibrations but **no** science
    (as a calibration (re)build leaves it) still shows the planned science
    frames, re-seeded from the cached planned-frame list.
    """
    pf = _make_redux(tmp_path, 'healthy')         # calibs, no science entries
    # Warm the planned-frame cache (as the launch/derive would have), then load
    # the calib-only state file the way the post-(re)build refresh does.
    dash_model._PLANNED_SCIENCE_CACHE[str(Path(pf))] = [
        {'frame': 'b188-sci', 'objtype': 'science', 'calib_id': 0,
         'comb_id': None, 'raw_files': ['b188.fits.gz']}]
    try:
        m = dash_model.DashboardModel(pf, derive=False)
        assert m.load_status == dash_model.LOAD_STATE_FILE
        sci = m.science_table()
        assert len(sci) > 0 and 'b188-sci' in set(sci['frame'])
        # The planned frame is undone and carries its raw file (so it can be
        # (re)built once the calibrations are ready).
        row = sci[sci['frame'] == 'b188-sci'][0]
        assert row['process'] == 'undone'
        entry = m.science_frame_entry('b188-sci', 1)
        assert entry.raw_files == ['b188.fits.gz']
    finally:
        dash_model._PLANNED_SCIENCE_CACHE.pop(str(Path(pf)), None)


def test_merge_from_disk_preserves_other_portion(tmp_path):
    """
    Round-5: a step-runner's fresh state merges the existing on-disk
    calibrations + science, so writing it does not blank out the other
    portion (which disabled all (Re)Build buttons).
    """
    from pypeit.state.run_state import RunPypeItState
    # A state file with BOTH calibrations and science (the science fixture).
    sj = tmp_path / 'x_state.json'
    sj.write_text(Path(data_path('dashboard_state_science.json')).read_text())
    # A fresh run_state, as pypeit_reduce_by_step has: no science yet.
    fresh = RunPypeItState(pypeit_file=str(tmp_path / 'x.pypeit'),
                           current_step='init', current_det=-1,
                           current_calibID=-1)
    assert fresh.get_science_status() is None
    fresh.merge_from_disk()
    # Science carried forward from disk (so a calib build keeps it; Round-4),
    # and calibrations carried forward (so a science build keeps them; Round-5).
    sci = fresh.get_science_status()
    assert sci is not None and len(sci) == 2
    assert fresh.get_status() is not None and len(fresh.get_status()) > 0


def test_inspect_run_pypeit_command(tmp_path):
    """
    inspect builds the full-reduction `run_pypeit -o` argv (Round-6).
    """
    from pypeit.dashboard import inspect as dash_inspect
    m = _science_model(tmp_path)
    argv = dash_inspect.run_pypeit_command(m)
    assert argv[0] == 'run_pypeit' and '-o' in argv
    assert argv[argv.index('--redux_path') + 1] == str(m.redux_dir)


def test_science_view_run_pypeit_button(tmp_path, qtbot):
    """
    The Science view offers a view-level "Run PypeIt" button, lock-aware
    (Round-6).
    """
    from pypeit.dashboard.runlock import RunLock
    lock = RunLock()
    view, _ = _science_view(tmp_path, qtbot, run_lock=lock)
    assert view._run_pypeit_button is not None
    assert view._run_pypeit_available is True
    assert view._run_pypeit_button.isEnabled()
    assert view._run_pypeit_button.text() == 'Run PypeIt'
    # Locking turns it orange + disabled, like the (Re)Build controls.
    view.set_locked(True)
    assert not view._run_pypeit_button.isEnabled()
    assert 'progress' in view._run_pypeit_button.text().lower()
