"""
Test the logic of the pypeit_preprocess_jwst script that decides which files
to process and how the CRDS environment is set, without running the JWST
pipeline.
"""
import json
import os
import sys
import types
from pathlib import Path

import pytest

from pypeit import PypeItError
from pypeit.scripts.preprocess_jwst import PreprocessJWST, JWST_CRDS_SERVER, CRDS_SERVERLESS

ROOT = 'jw01967012001_03102_00001_nrs1'
CRDS_VARIABLES = ['CRDS_PATH', 'CRDS_CONTEXT', 'CRDS_SERVER_URL', 'CRDS_READONLY_CACHE']


class PipelineRecord:
    # What the stand-in JWST pipelines were called with: the files processed and
    # the step parameters used by each pipeline, and the CRDS environment
    # variables at each call
    def __init__(self):
        self.files = {'Detector1Pipeline': [], 'Spec2Pipeline': []}
        self.steps = {'Detector1Pipeline': [], 'Spec2Pipeline': []}
        self.env = []


class FakePipeline:
    # Stand-in for a JWST pipeline: it only records its calls in `record` (set
    # by the `pipeline` fixture), so that neither the jwst package nor data are
    # needed
    name = None
    record = None

    @classmethod
    def call(cls, input_file, **kwargs):
        cls.record.files[cls.name].append(Path(input_file).name)
        cls.record.steps[cls.name].append(kwargs.get('steps'))
        cls.record.env.append({k: os.environ.get(k) for k in CRDS_VARIABLES})


class FakeDetector1Pipeline(FakePipeline):
    name = 'Detector1Pipeline'


class FakeSpec2Pipeline(FakePipeline):
    name = 'Spec2Pipeline'


def resolve_crds_context(offline):
    # Stand-in for PreprocessJWST._resolve_crds_context, which needs CRDS (and,
    # if not offline, the CRDS server)
    return os.environ.get('CRDS_CONTEXT', 'jwst_9999.pmap')


def crds_not_configured(offline):
    # Stand-in for PreprocessJWST._resolve_crds_context when CRDS cannot be
    # configured
    raise RuntimeError('CRDS could not be configured')


@pytest.fixture
def pipeline(monkeypatch, tmp_path):
    # Replace the JWST pipeline with the stand-ins, with a new record for each test
    record = PipelineRecord()
    monkeypatch.setattr(FakePipeline, 'record', record)
    module = types.ModuleType('jwst.pipeline')
    module.Detector1Pipeline = FakeDetector1Pipeline
    module.Spec2Pipeline = FakeSpec2Pipeline
    monkeypatch.setitem(sys.modules, 'jwst', types.ModuleType('jwst'))
    monkeypatch.setitem(sys.modules, 'jwst.pipeline', module)
    monkeypatch.setattr(PreprocessJWST, '_resolve_crds_context',
                        staticmethod(resolve_crds_context))
    # Start from a known CRDS environment, restored after the test
    monkeypatch.setenv('CRDS_PATH', str(tmp_path / 'crds'))
    for key in CRDS_VARIABLES[1:]:
        monkeypatch.delenv(key, raising=False)
    monkeypatch.chdir(tmp_path)
    return record


def _run(tmp_path, *options):
    raw_dir = tmp_path / 'raw'
    raw_dir.mkdir(exist_ok=True)
    (raw_dir / f'{ROOT}_uncal.fits').touch()
    args = PreprocessJWST.parse_args([str(raw_dir), '-o', str(tmp_path / 'out'),
                                      '--log_file', 'None', *options])
    PreprocessJWST.main(args)


def _make_outputs(tmp_path, *suffixes):
    out = tmp_path / 'out'
    out.mkdir(exist_ok=True)
    for suffix in suffixes:
        (out / f'{ROOT}_{suffix}.fits').touch()


def test_skip_completed_stages(pipeline, tmp_path):
    # Both stages already completed: nothing is run
    _make_outputs(tmp_path, 'rate', 'assign_wcs', 'cal')
    _run(tmp_path)
    assert pipeline.files == {'Detector1Pipeline': [], 'Spec2Pipeline': []}


def test_rerun_incomplete_stage2(pipeline, tmp_path):
    # Stage 2 started (assign_wcs exists) but did not produce the cal file: it
    # must be run again
    _make_outputs(tmp_path, 'rate', 'assign_wcs')
    _run(tmp_path)
    assert pipeline.files == {'Detector1Pipeline': [], 'Spec2Pipeline': [f'{ROOT}_rate.fits']}


def test_overwrite(pipeline, tmp_path):
    # The overwrite flags force both stages to run again
    _make_outputs(tmp_path, 'rate', 'assign_wcs', 'cal')
    _run(tmp_path, '--overwrite_stage1', '--overwrite_stage2')
    assert pipeline.files == {'Detector1Pipeline': [f'{ROOT}_uncal.fits'],
                              'Spec2Pipeline': [f'{ROOT}_rate.fits']}


def test_crds_environment(pipeline, tmp_path, monkeypatch):
    # The stand-in pipeline does not write a cal file, so Stage 2 runs every time
    _make_outputs(tmp_path, 'rate')
    before = {k: os.environ.get(k) for k in CRDS_VARIABLES}

    # Defaults: the CRDS server is set, the context is the CRDS default
    _run(tmp_path)
    assert pipeline.env[-1] == dict(before, CRDS_SERVER_URL=JWST_CRDS_SERVER)

    # Options: all set for the run
    _run(tmp_path, '--crds_path', str(tmp_path / 'other'), '--crds_context', 'jwst_1413.pmap')
    assert pipeline.env[-1] == dict(before, CRDS_PATH=str(tmp_path / 'other'),
                                    CRDS_CONTEXT='jwst_1413.pmap',
                                    CRDS_SERVER_URL=JWST_CRDS_SERVER)

    # Offline: no server, read-only cache
    _run(tmp_path, '--crds_offline')
    assert pipeline.env[-1] == dict(before, CRDS_SERVER_URL=CRDS_SERVERLESS,
                                    CRDS_READONLY_CACHE='1')

    # The environment is restored after each run
    assert {k: os.environ.get(k) for k in CRDS_VARIABLES} == before

    # Variables already set are used, and restored after a run that overrides them
    monkeypatch.setenv('CRDS_SERVER_URL', 'https://mirror.example')
    monkeypatch.setenv('CRDS_CONTEXT', 'jwst-latest')
    _run(tmp_path)
    assert pipeline.env[-1]['CRDS_SERVER_URL'] == 'https://mirror.example'
    assert pipeline.env[-1]['CRDS_CONTEXT'] == 'jwst-latest'
    _run(tmp_path, '--crds_context', 'jwst_1413.pmap', '--crds_offline')
    assert pipeline.env[-1]['CRDS_CONTEXT'] == 'jwst_1413.pmap'
    assert pipeline.env[-1]['CRDS_SERVER_URL'] == CRDS_SERVERLESS
    assert os.environ['CRDS_SERVER_URL'] == 'https://mirror.example'
    assert os.environ['CRDS_CONTEXT'] == 'jwst-latest'
    assert 'CRDS_READONLY_CACHE' not in os.environ


def test_crds_log(pipeline, tmp_path, monkeypatch, caplog):
    # The log reports the CRDS server and context used, and where they come from
    _make_outputs(tmp_path, 'rate')
    for options, expected in [
            ((), [f'CRDS server: {JWST_CRDS_SERVER} (default)',
                  'CRDS context: jwst_9999.pmap (from the CRDS default']),
            (('--crds_context', 'jwst_1413.pmap', '--crds_offline'),
             ['CRDS server: none; only the local CRDS cache is used',
              'CRDS context: jwst_1413.pmap (from --crds_context)'])]:
        caplog.clear()
        with caplog.at_level('INFO'):
            _run(tmp_path, *options)
        for text in expected:
            assert text in caplog.text, caplog.text


def test_crds_not_configured(pipeline, tmp_path, monkeypatch):
    # Errors explain what to check, and nothing is run
    monkeypatch.setattr(PreprocessJWST, '_resolve_crds_context',
                        staticmethod(crds_not_configured))
    with pytest.raises(PypeItError, match='Check that the CRDS server'):
        _run(tmp_path)
    with pytest.raises(PypeItError, match='With --crds_offline, the local CRDS cache'):
        _run(tmp_path, '--crds_offline')
    assert pipeline.files == {'Detector1Pipeline': [], 'Spec2Pipeline': []}


def test_crds_offline_empty_cache(tmp_path, monkeypatch):
    # With the real CRDS: offline with an empty cache fails with an explanation
    pytest.importorskip('crds')
    monkeypatch.setenv('CRDS_PATH', str(tmp_path))
    monkeypatch.setenv('CRDS_SERVER_URL', CRDS_SERVERLESS)
    monkeypatch.setenv('CRDS_READONLY_CACHE', '1')
    with pytest.raises(RuntimeError, match='CRDS could not load its configuration'):
        PreprocessJWST._resolve_crds_context(True)


def test_crds_environment_restored_on_error(pipeline, tmp_path):
    # The CRDS environment is restored even if the run fails
    before = {k: os.environ.get(k) for k in CRDS_VARIABLES}
    with pytest.raises(Exception):
        PreprocessJWST.main(PreprocessJWST.parse_args(
            [str(tmp_path / 'missing'), '-o', str(tmp_path / 'out'), '--log_file', 'None',
             '--crds_context', 'jwst_1413.pmap', '--crds_offline']))
    assert {k: os.environ.get(k) for k in CRDS_VARIABLES} == before


def _write_config(tmp_path, params):
    config = tmp_path / 'spec2_user.json'
    config.write_text(json.dumps(params))
    return config


def test_pypeit_products_added(pipeline, tmp_path):
    # A user configuration without the parameters that save the products
    # needed by PypeIt: they are added
    _make_outputs(tmp_path, 'rate')
    config = _write_config(tmp_path, {'extract_1d': {'skip': False}})
    _run(tmp_path, '--spec2_config', str(config))
    steps = pipeline.steps['Spec2Pipeline'][-1]
    assert steps['extract_1d'] == {'skip': False}
    assert steps['assign_wcs']['save_results'] is True
    assert steps['flat_field']['save_interpolated_flat'] is True


def test_pypeit_products_forced(pipeline, tmp_path, caplog):
    # A user configuration that turns off the products needed by PypeIt: they
    # are turned back on, with a warning, and the saved configuration reflects
    # what was run
    _make_outputs(tmp_path, 'rate')
    config = _write_config(tmp_path, {'assign_wcs': {'save_results': False},
                                      'flat_field': {'save_interpolated_flat': False}})
    with caplog.at_level('WARNING'):
        _run(tmp_path, '--spec2_config', str(config))
    steps = pipeline.steps['Spec2Pipeline'][-1]
    assert steps['assign_wcs']['save_results'] is True
    assert steps['flat_field']['save_interpolated_flat'] is True
    assert 'setting assign_wcs.save_results = True' in caplog.text
    assert 'setting flat_field.save_interpolated_flat = True' in caplog.text
    saved = json.loads((tmp_path / 'out' / 'spec2_config.json').read_text())
    assert saved['assign_wcs']['save_results'] is True
    assert saved['flat_field']['save_interpolated_flat'] is True


def test_pypeit_products_step_skipped(pipeline, tmp_path, caplog):
    # Skipping a step that writes a product needed by PypeIt gives a warning
    _make_outputs(tmp_path, 'rate')
    config = _write_config(tmp_path, {'flat_field': {'skip': True}})
    with caplog.at_level('WARNING'):
        _run(tmp_path, '--spec2_config', str(config))
    assert 'The flat_field step is skipped' in caplog.text
