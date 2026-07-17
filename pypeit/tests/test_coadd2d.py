"""
Module to test coadd2d.py
"""
import pytest

from pypeit import coadd2d
from pypeit.par.pypeitpar import PypeItPar
from pypeit.pkg.outputpaths import PypeItOutputPaths


@pytest.fixture
def output_paths(tmp_path, monkeypatch):
    """
    Provide a dedicated, real (non-dry-run) `PypeItOutputPaths` instance,
    monkeypatched into `coadd2d`, so tests don't depend on (or interfere
    with) the package-level singleton's configuration state.
    """
    op = PypeItOutputPaths()
    op.configure(redux_path=str(tmp_path))
    monkeypatch.setattr(coadd2d, 'outputPaths', op)
    return op


def test_output_paths_no_par_mutation(output_paths):
    # par is no longer read or mutated by output_paths -- the qadir value
    # (and everything else) must come back bit-for-bit unchanged.
    par = PypeItPar()
    qadir = par['rdx']['qadir']
    coadd2d.CoAdd2D.output_paths(['spec2d_file.fits'], par)
    assert par['rdx']['qadir'] == qadir, 'output_paths must not mutate par'


def test_output_paths_never_configures(output_paths, monkeypatch):
    # output_paths must only read already-resolved outputPaths properties --
    # it must never call configure() itself (that's the calling script's
    # job, exactly once, per guideline 3).
    def _raise_if_called(*args, **kwargs):
        raise AssertionError('output_paths must not call configure()')
    monkeypatch.setattr(output_paths, 'configure', _raise_if_called)

    par = PypeItPar()
    coadd2d.CoAdd2D.output_paths(['spec2d_file.fits'], par)


def test_output_paths_matches_coadd_properties(output_paths):
    par = PypeItPar()
    sci_dir, qa_dir = coadd2d.CoAdd2D.output_paths(['spec2d_file.fits'], par)
    assert sci_dir == str(output_paths.coadd_science)
    assert qa_dir == str(output_paths.coadd_qa_pngs)
