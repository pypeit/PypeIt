"""
Test the logic of the pypeit_preprocess_jwst script that decides which files
to process, without running the JWST pipeline.
"""
import sys
import types

import pytest

from pypeit.scripts.preprocess_jwst import PreprocessJWST

ROOT = 'jw01967012001_03102_00001_nrs1'


@pytest.fixture
def pipeline_calls(monkeypatch, tmp_path):
    # Replace the JWST pipeline with stand-ins that only record the files they
    # are called on, so that neither the jwst package nor data are needed
    calls = {'Detector1Pipeline': [], 'Spec2Pipeline': []}

    def _pipeline(name):
        def call(input_file, **kwargs):
            calls[name].append(str(input_file))
        return type(name, (), {'call': staticmethod(call)})

    module = types.ModuleType('jwst.pipeline')
    module.Detector1Pipeline = _pipeline('Detector1Pipeline')
    module.Spec2Pipeline = _pipeline('Spec2Pipeline')
    monkeypatch.setitem(sys.modules, 'jwst', types.ModuleType('jwst'))
    monkeypatch.setitem(sys.modules, 'jwst.pipeline', module)
    # The script sets CRDS_PATH; make sure it is restored
    monkeypatch.setenv('CRDS_PATH', str(tmp_path / 'crds'))
    monkeypatch.chdir(tmp_path)
    return calls


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


def test_skip_completed_stages(pipeline_calls, tmp_path):
    # Both stages already completed: nothing is run
    _make_outputs(tmp_path, 'rate', 'assign_wcs', 'cal')
    _run(tmp_path)
    assert pipeline_calls == {'Detector1Pipeline': [], 'Spec2Pipeline': []}


def test_rerun_incomplete_stage2(pipeline_calls, tmp_path):
    # Stage 2 started (assign_wcs exists) but did not produce the cal file: it
    # must be run again
    _make_outputs(tmp_path, 'rate', 'assign_wcs')
    _run(tmp_path)
    assert pipeline_calls['Detector1Pipeline'] == []
    assert [s.split('/')[-1] for s in pipeline_calls['Spec2Pipeline']] == [f'{ROOT}_rate.fits']


def test_overwrite(pipeline_calls, tmp_path):
    # The overwrite flags force both stages to run again
    _make_outputs(tmp_path, 'rate', 'assign_wcs', 'cal')
    _run(tmp_path, '--overwrite_stage1', '--overwrite_stage2')
    assert [s.split('/')[-1] for s in pipeline_calls['Detector1Pipeline']] == [f'{ROOT}_uncal.fits']
    assert [s.split('/')[-1] for s in pipeline_calls['Spec2Pipeline']] == [f'{ROOT}_rate.fits']
