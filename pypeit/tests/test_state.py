"""
Unit tests for the PypeIt reduction-state class and its I/O.

These exercise :class:`pypeit.state.run_state.RunPypeItState` and the
``pypeit.state`` helpers in isolation.  They are CI-safe: they require no
RAW_DATA and write only to pytest's ``tmp_path``.

Coverage includes regression tests for the safety fixes made to the state
machinery:

* **C1** — the ``status`` Literal must reject the invalid string
  ``'failed'`` and a ``'fail'`` status must survive a write/load round trip.
* **C3** — a mosaic detector (a tuple, stored by pydantic as a list) must
  match an existing entry so repeated updates do not spawn duplicates.
* **C5** — ``safe_write`` / ``safe_update_calib`` must log-and-continue
  (never raise) on failure.
* **S1** — ``get_status`` / ``print_status`` must list steps in the
  pipeline processing order (``slits`` before ``arc``).
"""

from pathlib import Path

import pytest

from pydantic import ValidationError

from pypeit.state import run_state
from pypeit.state import science_status
from pypeit import PypeItError
from pypeit.calibrations import Calibrations, MultiSlitCalibrations


def make_state(tmp_path, name='test.pypeit'):
    """
    Build a minimal, empty :class:`~pypeit.state.run_state.RunPypeItState` whose
    output file lives under ``tmp_path``.

    Args:
        tmp_path (:obj:`pathlib.Path`):
            pytest temporary directory fixture.
        name (:obj:`str`, optional):
            Basename of the ``.pypeit`` file to associate with the state.

    Returns:
        :class:`~pypeit.state.run_state.RunPypeItState`: A fresh state object with
        its ``path`` (output JSON) set under ``tmp_path``.
    """
    pyfile = str(tmp_path / name)
    s = run_state.RunPypeItState(pypeit_file=pyfile, current_step='init',
                             current_det=-1, current_calibID=-1)
    # Write the JSON beside the .pypeit file
    s.path = str(tmp_path / name.replace('.pypeit', '_state.json'))
    return s


# -----------------------------------------------------------------------
# Construction and the outfile property
# -----------------------------------------------------------------------

def test_construction_defaults():
    """
    A freshly built state has empty per-step lists and the expected
    required fields.
    """
    s = run_state.RunPypeItState(pypeit_file='foo.pypeit', current_step='init',
                             current_det=-1, current_calibID=-1)
    assert s.previous_step == 'none'
    for step in run_state.calib_classes:
        assert getattr(s, step) == []


def test_outfile_property():
    """
    ``outfile`` derives from the ``.pypeit`` name unless ``path`` is set.
    """
    s = run_state.RunPypeItState(pypeit_file='/a/b/run.pypeit',
                             current_step='init', current_det=-1,
                             current_calibID=-1)
    assert s.outfile == '/a/b/run_state.json'
    s.path = '/elsewhere/state.json'
    assert s.outfile == '/elsewhere/state.json'


# -----------------------------------------------------------------------
# same_det helper (C3)
# -----------------------------------------------------------------------

def test_same_det_scalars_and_mosaics():
    """
    ``same_det`` compares single detectors and mosaics regardless of
    list/tuple storage.
    """
    assert run_state.same_det(1, 1)
    assert not run_state.same_det(1, 2)
    # A mosaic stored as a list must match the same mosaic as a tuple
    assert run_state.same_det([1, 5], (1, 5))
    assert run_state.same_det((1, 5), [1, 5])
    assert not run_state.same_det((1, 5), (2, 6))
    # A scalar is never equal to a mosaic
    assert not run_state.same_det(1, (1, 5))


# -----------------------------------------------------------------------
# update_calib
# -----------------------------------------------------------------------

def test_update_calib_create_and_update():
    """
    The first ``update_calib`` creates the entry; subsequent calls update
    the same entry in place.
    """
    s = run_state.RunPypeItState(pypeit_file='x.pypeit', current_step='init',
                             current_det=-1, current_calibID=-1)
    s.update_calib('arc', 0, 1, 'status', 'running')
    assert len(s.arc) == 1
    s.update_calib('arc', 0, 1, 'status', 'success')
    s.update_calib('arc', 0, 1, 'output_file', '/path/Arc.fits')
    # Still a single entry, now carrying both updates
    assert len(s.arc) == 1
    assert s.arc[0].status == 'success'
    assert s.arc[0].output_file == '/path/Arc.fits'


def test_update_calib_input_files_stored_as_list():
    """
    ``input_files`` defaults to ``None``, so a passed list is stored
    wholesale (this mirrors the real usage, where ``self.raw_files`` — a
    list — is handed in by ``base_state``).
    """
    s = run_state.RunPypeItState(pypeit_file='x.pypeit', current_step='init',
                             current_det=-1, current_calibID=-1)
    s.update_calib('arc', 0, 1, 'input_files', ['/a.fits', '/b.fits'])
    assert s.arc[0].input_files == ['/a.fits', '/b.fits']


def test_update_calib_input_files_repeat_no_nesting(tmp_path):
    """
    Regression: calling ``update_calib`` with ``input_files`` more than
    once must *replace* the list with a flat list, not nest the new list
    inside the old one (the latter corrupts the field and breaks
    ``load()`` validation).
    """
    s = make_state(tmp_path)
    # First call stores the list; a second call (e.g. a re-run of the
    # step) must not produce a nested list
    s.update_calib('arc', 0, 1, 'input_files', ['/a.fits'])
    s.update_calib('arc', 0, 1, 'input_files', ['/a.fits', '/b.fits'])
    assert s.arc[0].input_files == ['/a.fits', '/b.fits']
    # Every element must be a string (no nested list), so it round-trips
    assert all(isinstance(f, str) for f in s.arc[0].input_files)
    s.write()
    loaded = s.load()
    assert loaded.arc[0].input_files == ['/a.fits', '/b.fits']


def test_update_calib_list_field_appends():
    """
    Assigning a scalar to a field that already holds a list appends rather
    than overwrites.  ``FlatsState.corrections`` defaults to ``[]``, so
    successive scalar updates accumulate.
    """
    s = run_state.RunPypeItState(pypeit_file='x.pypeit', current_step='init',
                             current_det=-1, current_calibID=-1)
    s.update_calib('flats', 0, 1, 'corrections', 'pixelflat')
    s.update_calib('flats', 0, 1, 'corrections', 'spat_illum')
    assert s.flats[0].corrections == ['pixelflat', 'spat_illum']


def test_update_calib_tracks_step_transitions():
    """
    ``current_step`` / ``previous_step`` are updated on a step change.
    """
    s = run_state.RunPypeItState(pypeit_file='x.pypeit', current_step='init',
                             current_det=-1, current_calibID=-1)
    s.update_calib('arc', 0, 1, 'status', 'running')
    assert s.current_step == 'arc'
    assert s.previous_step == 'init'
    s.update_calib('tilts', 0, 1, 'status', 'running')
    assert s.current_step == 'tilts'
    assert s.previous_step == 'arc'


def test_update_calib_unknown_step_is_noop():
    """
    An unrecognized step is silently ignored (no entry, no crash).
    """
    s = run_state.RunPypeItState(pypeit_file='x.pypeit', current_step='init',
                             current_det=-1, current_calibID=-1)
    # 'bpm' is not tracked in calib_classes
    s.update_calib('bpm', 0, 1, 'status', 'running')
    assert all(getattr(s, step) == [] for step in run_state.calib_classes)


def test_update_calib_per_slit():
    """
    A per-slit update creates the slit sub-entry and sets its metric.
    """
    s = run_state.RunPypeItState(pypeit_file='x.pypeit', current_step='init',
                             current_det=-1, current_calibID=-1)
    s.update_calib('wv_calib', 0, 1, 'status', 'success', slit=196)
    s.update_calib('wv_calib', 0, 1, 'rms', 0.06, slit=196)
    entry = s.wv_calib[0]
    assert 196 in entry.slits
    assert entry.slits[196].status == 'success'
    assert entry.slits[196].rms == pytest.approx(0.06)


def test_update_calib_mosaic_det_single_entry():
    """
    **C3 regression.** Repeated updates with a mosaic ``det`` (tuple) must
    update a single entry, not create duplicates, and all values must land
    on that one entry.
    """
    s = run_state.RunPypeItState(pypeit_file='x.pypeit', current_step='init',
                             current_det=-1, current_calibID=-1)
    det = (1, 5)
    s.update_calib('bias', 0, det, 'status', 'running')
    s.update_calib('bias', 0, det, 'status', 'success')
    s.update_calib('bias', 0, det, 'mean', 1.0)
    assert len(s.bias) == 1
    assert s.bias[0].status == 'success'
    assert s.bias[0].mean == pytest.approx(1.0)
    # Distinct mosaics remain distinct entries
    s.update_calib('bias', 0, (2, 6), 'status', 'running')
    assert len(s.bias) == 2


# -----------------------------------------------------------------------
# Flats state (corrections, provenance, grouped files, per-slit metrics)
# -----------------------------------------------------------------------

def test_flats_slit_model_and_skip_status():
    """
    ``FlatsSlit`` accepts the flats-specific ``'skip'`` status and holds a
    dict of per-correction ``FlatCorrectionMetric``; ``'flats'`` is wired
    into ``slit_classes``.
    """
    m = run_state.FlatCorrectionMetric(mean=1.0, rms=0.02)
    fs = run_state.FlatsSlit(status='skip', corrections={'pixelflat': m})
    assert fs.status == 'skip'
    assert fs.corrections['pixelflat'].rms == pytest.approx(0.02)
    assert run_state.slit_classes['flats'] is run_state.FlatsSlit
    # An invalid slit status is rejected
    with pytest.raises(ValidationError):
        run_state.FlatsSlit(status='running')


def test_flats_per_slit_metrics_via_update_calib():
    """
    The per-slit ``corrections`` dict is set through the existing ``slit=``
    path of ``update_calib`` (no special machinery), and an entry-level
    ``corrections`` list records which corrections are present.
    """
    s = run_state.RunPypeItState(pypeit_file='x.pypeit', current_step='init',
                             current_det=-1, current_calibID=-1)
    # Entry-level: which corrections exist (list set wholesale)
    s.update_calib('flats', 0, 1, 'corrections', ['pixelflat', 'spat_illum'])
    # Per-slit status + per-correction metrics
    s.update_calib('flats', 0, 1, 'status', 'success', slit=196)
    metrics = {'pixelflat': run_state.FlatCorrectionMetric(mean=1.0, rms=0.026),
               'spat_illum': run_state.FlatCorrectionMetric(mean=1.0, rms=0.014)}
    s.update_calib('flats', 0, 1, 'corrections', metrics, slit=196)
    entry = s.flats[0]
    assert entry.corrections == ['pixelflat', 'spat_illum']
    assert entry.slits[196].status == 'success'
    spat = entry.slits[196].corrections['spat_illum']
    assert spat.rms == pytest.approx(0.014)


def test_flats_state_full_roundtrip(tmp_path):
    """
    A fully populated ``FlatsState`` (grouped input files, provenance,
    corrections, per-slit metrics) survives a write/load round trip.
    """
    s = make_state(tmp_path)
    s.update_calib('flats', 0, 1, 'corrections', ['pixelflat', 'spat_illum'])
    s.update_calib('flats', 0, 1, 'pixelflat_source', 'raw')
    s.update_calib('flats', 0, 1, 'pixelflat_files', ['/a.fits', '/b.fits'])
    s.update_calib('flats', 0, 1, 'illumflat_files', ['/a.fits'])
    s.update_calib('flats', 0, 1, 'input_files', ['/a.fits', '/b.fits'])
    s.update_calib('flats', 0, 1, 'status', 'success', slit=196)
    metric = run_state.FlatCorrectionMetric(mean=1.0, rms=0.03)
    s.update_calib('flats', 0, 1, 'corrections', {'pixelflat': metric},
                   slit=196)
    s.write()
    loaded = s.load()
    fe = loaded.flats[0]
    assert fe.corrections == ['pixelflat', 'spat_illum']
    assert fe.pixelflat_source == 'raw'
    assert fe.pixelflat_files == ['/a.fits', '/b.fits']
    assert fe.illumflat_files == ['/a.fits']
    assert fe.slits[196].status == 'success'
    assert fe.slits[196].corrections['pixelflat'].mean == pytest.approx(1.0)
    assert fe.slits[196].corrections['pixelflat'].rms == pytest.approx(0.03)


# -----------------------------------------------------------------------
# Status Literal (C1)
# -----------------------------------------------------------------------

def test_status_literal_rejects_failed():
    """
    **C1 regression.** ``'failed'`` is not a valid status; ``'fail'`` is.
    Validation (i.e. load) must reject the former.
    """
    with pytest.raises(ValidationError):
        run_state.BiasCalibState(calib_id=0, det=1, status='failed')
    # The corrected value validates fine
    ok = run_state.BiasCalibState(calib_id=0, det=1, status='fail')
    assert ok.status == 'fail'


# -----------------------------------------------------------------------
# write / load round trip
# -----------------------------------------------------------------------

def test_write_load_roundtrip(tmp_path):
    """
    A state written to disk reloads to an equivalent object.
    """
    s = make_state(tmp_path)
    s.update_calib('arc', 0, 1, 'status', 'success')
    s.update_calib('arc', 0, 1, 'output_file', '/path/Arc.fits')
    s.update_calib('wv_calib', 0, 1, 'status', 'success', slit=196)
    s.update_calib('wv_calib', 0, 1, 'rms', 0.06, slit=196)
    s.write()
    assert Path(s.outfile).is_file()

    loaded = s.load()
    assert loaded.arc[0].status == 'success'
    assert loaded.arc[0].output_file == '/path/Arc.fits'
    assert loaded.wv_calib[0].slits[196].rms == pytest.approx(0.06)
    # Full model equality is the strongest check
    assert loaded.model_dump() == s.model_dump()


def test_load_fail_status_roundtrip(tmp_path):
    """
    **C1 regression (end-to-end).** A ``'fail'`` status survives a
    write/load round trip without raising.
    """
    s = make_state(tmp_path)
    s.update_calib('bias', 0, 1, 'status', 'fail')
    s.write()
    loaded = s.load()
    assert loaded.bias[0].status == 'fail'


def test_load_missing_file_returns_self(tmp_path):
    """
    ``load`` is a no-op (returns the same object) when no file exists.
    """
    s = make_state(tmp_path)
    assert not Path(s.outfile).exists()
    assert s.load() is s


def test_load_mosaic_det_roundtrip(tmp_path):
    """
    A mosaic detector round-trips and re-matches after load (the stored
    list re-matches the tuple via ``same_det``).
    """
    s = make_state(tmp_path)
    s.update_calib('bias', 0, (1, 5), 'status', 'success')
    s.write()
    loaded = s.load()
    # The reloaded entry still matches the mosaic and updates in place
    loaded.update_calib('bias', 0, (1, 5), 'mean', 2.0)
    assert len(loaded.bias) == 1
    assert loaded.bias[0].mean == pytest.approx(2.0)


# -----------------------------------------------------------------------
# safe_write / safe_update_calib (C5)
# -----------------------------------------------------------------------

def test_safe_write_success_and_failure(tmp_path):
    """
    **C5.** ``safe_write`` returns True on success and False (without
    raising) when the path is unwritable.
    """
    s = make_state(tmp_path)
    assert s.safe_write() is True
    assert Path(s.outfile).is_file()
    # Point at a non-existent directory to force a failure
    s.path = str(tmp_path / 'no_such_dir' / 'state.json')
    assert s.safe_write() is False


def test_safe_update_calib_success_and_failure(tmp_path):
    """
    **C5.** ``safe_update_calib`` returns True on success and False on a
    bad update, never raising.
    """
    s = make_state(tmp_path)
    assert s.safe_update_calib('arc', 0, 1, 'status', 'success') is True
    assert s.arc[0].status == 'success'
    # An invalid key triggers the guarded failure path
    assert s.safe_update_calib('arc', 0, 1, 'no_such_field', 1) is False


# -----------------------------------------------------------------------
# get_status / print_status (S1 ordering)
# -----------------------------------------------------------------------

def test_get_status_empty_is_none():
    """
    ``get_status`` returns None when there are no entries.
    """
    s = run_state.RunPypeItState(pypeit_file='x.pypeit', current_step='init',
                             current_det=-1, current_calibID=-1)
    assert s.get_status() is None


def test_get_status_columns_and_values():
    """
    ``get_status`` returns the expected columns and reports required /
    status / output basename per step.
    """
    s = run_state.RunPypeItState(pypeit_file='x.pypeit', current_step='init',
                             current_det=-1, current_calibID=-1)
    s.update_calib('arc', 0, 1, 'required', True)
    s.update_calib('arc', 0, 1, 'status', 'success')
    s.update_calib('arc', 0, 1, 'output_file', '/p/Arc_A_0_DET01.fits')
    df = s.get_status()
    assert df.colnames == ['calibration_group', 'detector', 'steps',
                                'required', 'status', 'output_file']
    arc_row = df[df['steps'] == 'arc'][0]
    assert arc_row['required'] == 'True'
    assert arc_row['status'] == 'success'
    assert arc_row['output_file'] == 'Arc_A_0_DET01.fits'
    # A step with no entry shows placeholders
    bias_row = df[df['steps'] == 'bias'][0]
    assert bias_row['status'] == '--'
    assert bias_row['output_file'] == '--'


def test_get_status_ordering_matches_processing_order():
    """
    **S1 regression.** The step order in ``get_status`` matches the
    pipeline processing order (``MultiSlitCalibrations.default_steps()``,
    with ``bpm`` dropped and ``align`` appended) and is independent of the
    order entries were added.
    """
    s = run_state.RunPypeItState(pypeit_file='x.pypeit', current_step='init',
                             current_det=-1, current_calibID=-1)
    # Deliberately add flats before slits to prove ordering is not
    # insertion order
    s.update_calib('flats', 0, 1, 'status', 'success')
    s.update_calib('slits', 0, 1, 'status', 'success')
    s.update_calib('arc', 0, 1, 'status', 'success')
    steps = list(s.get_status()['steps'])
    # slits must come before arc (the bug fixed in S1)
    assert steps.index('slits') < steps.index('arc')
    # And the full order matches calib_classes == default_steps order
    assert steps == list(run_state.calib_classes.keys())
    expected = [st for st in MultiSlitCalibrations.default_steps()
                if st != 'bpm'] + ['align']
    assert steps == expected


def test_print_status_runs(capsys):
    """
    ``print_status`` prints a table (smoke test; no entries -> friendly
    message).
    """
    s = run_state.RunPypeItState(pypeit_file='x.pypeit', current_step='init',
                             current_det=-1, current_calibID=-1)
    s.print_status()
    out = capsys.readouterr().out
    assert 'No calibration state entries found.' in out
    # With an entry, the step name appears
    s.update_calib('arc', 0, 1, 'status', 'success')
    s.print_status()
    out = capsys.readouterr().out
    assert 'arc' in out


# -----------------------------------------------------------------------
# Science-frame state
# -----------------------------------------------------------------------

def test_add_or_get_science_creates_once():
    """
    ``add_or_get_science`` creates an entry on first call and returns the
    same one thereafter (matched by (frame, det), mosaic-aware).
    """
    s = run_state.RunPypeItState(pypeit_file='x.pypeit', current_step='init',
                             current_det=-1, current_calibID=-1)
    e1 = s.add_or_get_science('b27', 1, calib_id=0, objtype='science')
    e2 = s.add_or_get_science('b27', 1)
    assert e1 is e2
    assert len(s.science) == 1
    # A mosaic tuple stored as a list still matches
    s.add_or_get_science('b50', (1, 5), objtype='standard')
    again = s.add_or_get_science('b50', (1, 5))
    assert len(s.science) == 2
    assert run_state.same_det(again.det, (1, 5))


def test_safe_update_science_sets_step_and_fields():
    """
    ``safe_update_science`` sets a step status and top-level fields, and
    returns False (no raise) on a bad field.
    """
    s = run_state.RunPypeItState(pypeit_file='x.pypeit', current_step='init',
                             current_det=-1, current_calibID=-1)
    e = s.safe_update_science('b27', 1, step='findobj', status='success',
                              nobj=3, calib_id=0)
    assert e.findobj.status == 'success'
    assert e.nobj == 3
    # Unknown field -> guarded failure, returns None, does not raise
    assert s.safe_update_science('b27', 1, badfield=1) is None


def test_science_step_status_literal():
    """
    ``ScienceStep`` accepts the four science statuses and ``ScienceSlit``
    accepts ``'skip'``; invalid values are rejected.
    """
    assert run_state.ScienceStep(status='running').status == 'running'
    assert run_state.ScienceSlit(status='skip').status == 'skip'
    with pytest.raises(ValidationError):
        run_state.ScienceStep(status='skip')          # not valid for a step


def test_science_full_roundtrip(tmp_path):
    """
    A populated science entry (steps, products, per-slit, per-object)
    survives a write/load round trip.
    """
    s = make_state(tmp_path)
    e = s.add_or_get_science('b27', 1, calib_id=0, objtype='science',
                             comb_id=0)
    for step in run_state.science_steps:
        getattr(e, step).status = 'success'
    e.nobj = 1
    e.spec1d_file = 'Science/spec1d_b27.fits'
    e.spec2d_file = 'Science/spec2d_b27.fits'
    e.slits[175] = run_state.ScienceSlit(status='success', nobj=1)
    e.objects.append(run_state.ScienceObj(objid=1, slitid=175, spat_pixpos=175.5,
                                      fwhm=5.1, snr_find=142.5, s2n=9.3,
                                      sign=1, extracted=True))
    s.write()
    loaded = s.load()
    le = loaded.science[0]
    assert le.objtype == 'science'
    assert le.extract.status == 'success'
    assert le.objects[0].snr_find == pytest.approx(142.5)
    assert le.objects[0].s2n == pytest.approx(9.3)
    assert le.objects[0].extracted is True
    assert le.slits[175].nobj == 1


def test_get_science_status_table():
    """
    ``get_science_status`` returns one row per (frame, det) with the four
    macro-step columns; None when empty.
    """
    s = run_state.RunPypeItState(pypeit_file='x.pypeit', current_step='init',
                             current_det=-1, current_calibID=-1)
    assert s.get_science_status() is None
    s.safe_update_science('b27', 1, step='process', status='success')
    df = s.get_science_status()
    assert df.colnames == ['frame', 'detector', 'calib', 'objtype',
                                'process', 'findobj', 'skysub', 'extract',
                                'nobj', 'spec2d', 'spec1d']
    row = df[0]
    assert row['process'] == 'success'
    assert row['findobj'] == 'undone'


# -----------------------------------------------------------------------
# run_the_steps failure marking (bug report 000)
# -----------------------------------------------------------------------

class _StepRunner:
    """
    Minimal stand-in that drives :meth:`Calibrations.run_the_steps` with
    controllable per-step behavior, so the failure-marking path is unit-
    testable without a full :class:`~pypeit.calibrations.Calibrations`
    instance (which needs raw data, a spectrograph, metadata, ...).

    A ``get_<step>`` for ``raising_step`` raises a :class:`PypeItError`
    (mimicking an unrecoverable error mid-build, as in bug report 000); all
    other ``get_<step>`` / ``<step>_state`` methods are no-ops.
    """

    # Borrow the real method under test.
    run_the_steps = Calibrations.run_the_steps

    def __init__(self, state_obj, steps, raising_step):
        """

        Args:
            state_obj (:class:`~pypeit.state.run_state.RunPypeItState`): The state to
                update.
            steps (:obj:`list`): Calibration step names to iterate.
            raising_step (:obj:`str`): The step whose ``get_`` method raises.
        """
        self.state = state_obj
        self.det = 1
        self.calib_ID = 0
        self.steps = steps
        self.success = False
        self.failed_step = None
        self._raising_step = raising_step

    def __getattr__(self, name):
        """
        Supply ``get_<step>`` / ``<step>_state`` methods on demand.

        Args:
            name (:obj:`str`): The attribute being looked up.

        Returns:
            callable: The synthesized step method.
        """
        if name.startswith('get_'):
            step = name[len('get_'):]

            def _run(force=None, step=step):
                if step == self._raising_step:
                    raise PypeItError(f'boom in {step}')
            return _run
        if name.endswith('_state'):
            return lambda: None
        raise AttributeError(name)


def test_run_the_steps_marks_fail_on_exception(tmp_path):
    """
    **Bug report 000 regression.** When a calibration step raises mid-build,
    :meth:`Calibrations.run_the_steps` must mark that step ``'fail'`` (and
    persist it) before re-raising — never leave it stuck at ``'running'``, or
    the Dashboard/``pypeit_status`` would show it as in-progress forever.
    """
    s = make_state(tmp_path)
    runner = _StepRunner(s, ['bias', 'flats'], raising_step='flats')
    # The run still aborts (the error propagates with a non-zero exit).
    with pytest.raises(PypeItError):
        runner.run_the_steps()
    # The failing step is recorded as 'fail', not left at 'running'.
    assert s.flats[0].status == 'fail'
    # Localized: the preceding step is not marked failed.
    assert s.bias[0].status != 'fail'
    # And the failure was persisted to disk, so a reload sees it too.
    loaded = s.load()
    assert loaded.flats[0].status == 'fail'


# -----------------------------------------------------------------------
# science_status pure helpers (CI-safe; no RAW_DATA / products needed)
# -----------------------------------------------------------------------

def test_science_status_safe_float():
    """
    ``_safe_float`` casts numbers/strings and returns None for None or
    uncastable input, never raising.
    """
    assert science_status._safe_float(3) == pytest.approx(3.0)
    assert science_status._safe_float('2.5') == pytest.approx(2.5)
    assert science_status._safe_float(None) is None
    assert science_status._safe_float('not-a-number') is None
    assert science_status._safe_float([1, 2]) is None


def test_science_status_det_from_name():
    """
    ``_det_from_name`` parses ``DET##`` to a 1-indexed int and returns None
    for mosaics / unparseable names.
    """
    assert science_status._det_from_name('DET01') == 1
    assert science_status._det_from_name('DET07') == 7
    assert science_status._det_from_name('MSC01') is None   # mosaic
    assert science_status._det_from_name('DETxx') is None   # unparseable
    assert science_status._det_from_name(None) is None


def test_science_status_basename_from_product():
    """
    ``_basename_from_product`` extracts ``<basename>`` from
    ``<prefix>_<basename>.fits`` and returns None on a mismatch.
    """
    assert science_status._basename_from_product(
        '/redux/Science/spec2d_b27-sci.fits', 'spec2d') == 'b27-sci'
    assert science_status._basename_from_product(
        'spec1d_J1234+5678_DEIMOS.fits', 'spec1d') == 'J1234+5678_DEIMOS'
    # Wrong prefix / not a .fits -> None
    assert science_status._basename_from_product(
        'spec2d_b27.fits', 'spec1d') is None
    assert science_status._basename_from_product(
        'spec2d_b27.txt', 'spec2d') is None


def test_science_status_exposure_ids():
    """
    ``_exposure_ids`` returns ``(comb_id, bkg_id)`` with negative (unset)
    values mapped to None, and ``(None, None)`` if a column is missing.
    """
    tbl = {'comb_id': [5, 6], 'bkg_id': [-1, 9]}
    assert science_status._exposure_ids(tbl, [0]) == (5, None)
    assert science_status._exposure_ids(tbl, [1]) == (6, 9)
    # Missing column -> graceful (None, None), never raises.
    assert science_status._exposure_ids({'comb_id': [5]}, [0]) == (None, None)


def test_derive_science_from_disk_empty_dir_is_noop(tmp_path):
    """
    ``derive_science_from_disk`` on a reduction directory with no
    ``Science/``/``Intermediate/`` products leaves the state's science empty
    and returns it (the graceful no-product path).
    """
    s = make_state(tmp_path)
    out = science_status.derive_science_from_disk(s, str(tmp_path))
    assert out is s
    assert s.get_science_status() is None
