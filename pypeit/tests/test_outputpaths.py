"""
Module to test pypeit.pkg.outputpaths.PypeItOutputPaths
"""
import logging

import pytest

from pypeit import PypeItPathError
from pypeit import outputfiles
from pypeit.pkg.outputpaths import PypeItOutputPaths
from pypeit.par.pypeitpar import PypeItPar


def test_default_resolution(tmp_path):
    # An unconfigured instance still resolves every property to the
    # expected Path, but nothing is created on disk and nothing is ready.
    op = PypeItOutputPaths(redux_path=tmp_path)

    assert op.redux == tmp_path
    assert op.science == tmp_path / 'Science'
    assert op.qa == tmp_path / 'QA'
    assert op.qa_pngs == tmp_path / 'QA' / 'PNGs'
    assert op.calibrations == tmp_path / 'Calibrations'
    assert op.collate == tmp_path
    assert op.coadd_science == tmp_path / 'Science_coadd'
    assert op.coadd_qa == tmp_path / 'QA_coadd'
    assert op.coadd_qa_pngs == tmp_path / 'QA_coadd' / 'PNGs'

    assert not op.configured
    assert not any(rec.ready for rec in op._paths.values())
    assert not tmp_path.joinpath('Science').exists()


def test_explicit_dirs(tmp_path):
    redux_path = tmp_path / 'redux'
    collate_path = tmp_path / 'collate_out'
    op = PypeItOutputPaths(redux_path=redux_path, scidir='sci', qadir='qa',
                            calib_dir='calib', coadd_suffix='_2d',
                            collate_outdir=collate_path)

    assert op.science == redux_path / 'sci'
    assert op.qa == redux_path / 'qa'
    assert op.calibrations == redux_path / 'calib'
    assert op.coadd_science == redux_path / 'sci_2d'
    assert op.coadd_qa == redux_path / 'qa_2d'
    assert op.collate == collate_path

    # collate is independent of redux, not aliased to it
    assert op._paths['collate'].parent != op._paths['redux'].parent \
        or op._paths['collate'].name != op._paths['redux'].name


def test_configure_creates_on_access(tmp_path, caplog):
    op = PypeItOutputPaths()
    op.configure(redux_path=tmp_path)

    assert not op._paths['science'].ready

    with caplog.at_level(logging.INFO, logger='pypeit'):
        science_dir = op.science
    assert science_dir == tmp_path / 'Science'
    assert science_dir.is_dir()
    assert op._paths['science'].ready
    n_msgs = sum('science' in rec.message for rec in caplog.records)
    assert n_msgs == 1

    # A second access does not re-create or re-log.
    caplog.clear()
    with caplog.at_level(logging.INFO, logger='pypeit'):
        _ = op.science
    assert not any('science' in rec.message for rec in caplog.records)


def test_nested_pngs_creates_parent(tmp_path):
    op = PypeItOutputPaths()
    op.configure(redux_path=tmp_path)

    # Access qa_pngs directly, without ever touching qa first.
    qa_pngs_dir = op.qa_pngs
    assert qa_pngs_dir == tmp_path / 'QA' / 'PNGs'
    assert qa_pngs_dir.is_dir()
    assert (tmp_path / 'QA').is_dir()
    assert op._paths['qa_pngs'].ready

    # qa's own record is independent -- its parent/name were set directly
    # by configure(), not derived from qa_pngs, so it is not marked ready
    # just because the directory it points to happens to already exist.
    assert not op._paths['qa'].ready


def test_configure_twice_raises(tmp_path):
    op = PypeItOutputPaths()
    op.configure(redux_path=tmp_path)
    assert op.configured

    with pytest.raises(PypeItPathError):
        op.configure(redux_path=tmp_path / 'other')

    # Also raises even when the requested values are identical -- no
    # idempotent no-op.
    with pytest.raises(PypeItPathError):
        op.configure(redux_path=tmp_path)


def test_dryrun_configure_repeatable(tmp_path):
    op = PypeItOutputPaths(dryrun=True)
    op.configure(redux_path=tmp_path / 'first')
    assert op.science == tmp_path / 'first' / 'Science'

    # A second configure() call is allowed because dryrun is True, and
    # fully resets state to reflect the new redux_path.
    op.configure(redux_path=tmp_path / 'second')
    assert op.science == tmp_path / 'second' / 'Science'

    # Nothing is ever created on disk or marked ready in dry-run mode.
    assert not op._paths['science'].ready
    assert not (tmp_path / 'second' / 'Science').exists()


def test_configure_logs(caplog):
    op = PypeItOutputPaths()
    with caplog.at_level(logging.INFO, logger='pypeit'):
        op.configure()
    n_msgs = sum('Output paths configured' in rec.message for rec in caplog.records)
    assert n_msgs == 1


def test_configure_from_par(tmp_path):
    par = PypeItPar()
    par['rdx']['redux_path'] = str(tmp_path)

    op = PypeItOutputPaths()
    op.configure(par)

    assert op.science == outputfiles.science_path(par)
    assert op.qa == tmp_path / par['rdx']['qadir']
    assert op.calibrations == tmp_path / par['calibrations']['calib_dir']


def test_derive_not_implemented():
    op = PypeItOutputPaths()
    with pytest.raises(NotImplementedError):
        op.derive()
