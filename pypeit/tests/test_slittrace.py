"""
Module to run tests on SlitTraceSet
"""
from pathlib import Path

from IPython import embed

import numpy as np

from astropy.io import fits
import pytest

from pypeit import PypeItError
from pypeit.slittrace import SlitTraceSet, SlitTraceBitMask
from pypeit.par.pypeitpar import EdgeTracePar
from pypeit.tests.tstutils import data_output_path


def test_bits():
    # Make sure bits are correct
    bm = SlitTraceBitMask()
    assert bm.bits['SHORTSLIT'] == 0, 'Bits changed'
    assert bm.bits['BOXSLIT'] == 1, 'Bits changed'
    assert bm.bits['USERIGNORE'] == 2, 'Bits changed'
    assert bm.bits['BADWVCALIB'] == 3, 'Bits changed'
    assert bm.bits['BADTILTCALIB'] == 4, 'Bits changed'
    assert bm.bits['BADALIGNCALIB'] == 5, 'Bits changed'
    assert bm.bits['SKIPFLATCALIB'] == 6, 'Bits changed'
    assert bm.bits['BADFLATCALIB'] == 7, 'Bits changed'
    assert bm.bits['BADREDUCE'] == 8, 'Bits changed'
    assert bm.bits['BADSKYSUB'] == 9, 'Bits changed'
    assert bm.bits['BADEXTRACT'] == 10, 'Bits changed'


def test_init():
    slits = SlitTraceSet(left_init=np.full((1000,3), 2, dtype=float),
                         right_init=np.full((1000,3), 8, dtype=float),
                         pypeline='MultiSlit',
                         nspat=10, PYP_SPEC='dummy')

    left, right, _ = slits.select_edges()
    center = (left+right)/2
    assert np.all(center == 5), 'Bad center'


@pytest.mark.parametrize('nslits', [1, 2, 3])
def test_traced_edges(nslits, tmp_path):
    left = np.tile(2. + 10. * np.arange(nslits), (100, 1))
    right = left + 6.
    slits = SlitTraceSet(left, right, 'MultiSlit', nspat=10 * nslits, PYP_SPEC='dummy')
    assert slits.traced_edges().shape == (nslits, 2), 'Incorrect shape'
    assert np.all(slits.traced_edges()), 'All edges should be traced by default'
    left_traced = np.arange(nslits) > 0
    right_traced = np.arange(nslits) < nslits - 1
    slits = SlitTraceSet(left, right, 'MultiSlit', nspat=10 * nslits, PYP_SPEC='dummy',
                         left_traced=left_traced, right_traced=right_traced)
    assert np.array_equal(slits.traced_edges(), np.column_stack((left_traced, right_traced))), \
        'Incorrect traced-edge flags'
    # I/O
    slits.set_paths(tmp_path, 'A', '1', 'DET01')
    slits.to_file()
    _slits = SlitTraceSet.from_file(slits.get_path())
    assert np.array_equal(_slits.traced_edges(), slits.traced_edges()), \
        'Bad traced-edge flags read'
    # Flags for one side only, or for the wrong number of slits
    with pytest.raises(PypeItError):
        SlitTraceSet(left, right, 'MultiSlit', nspat=10 * nslits, left_traced=left_traced)
    with pytest.raises(PypeItError):
        SlitTraceSet(left, right, 'MultiSlit', nspat=10 * nslits,
                     left_traced=np.ones(nslits + 1, dtype=bool),
                     right_traced=np.ones(nslits + 1, dtype=bool))


def test_io():

    slits = SlitTraceSet(np.full((1000,3), 2, dtype=float), np.full((1000,3), 8, dtype=float),
                         'MultiSlit', nspat=10, PYP_SPEC='dummy')
    slits.set_paths(data_output_path(''), 'A', '1', 'DET01')
    ofile = Path(slits.get_path()).absolute()

    # Try to save it
    slits.to_file(overwrite=True)
    assert ofile.exists(), 'File not written'

    # Instantiate an empty SlitTraceSet with the same file, and
    # indicate it should be reused
    _slits = SlitTraceSet.from_file(ofile)
    assert np.array_equal(_slits.left_init, np.full((1000,3), 2, dtype=float)), 'Bad left read'
    # And that it's the same as the existing one
    assert np.array_equal(_slits.left_init, slits.left_init), 'Bad left read'

    ofile.unlink()


def test_pad_par():
    # The slitedges `pad` parameter must be stored as a header-safe list (not a
    # numpy array), since it is written to a FITS header via ParSet.to_header.

    # A scalar is expanded to a symmetric 2-element list
    par = EdgeTracePar(pad=5)
    assert isinstance(par['pad'], list), 'pad should be stored as a list'
    assert par['pad'] == [5.0, 5.0], 'Scalar pad not expanded correctly'

    # A 2-element input (e.g. `pad = 8,3` in a .pypeit file) is preserved in order
    par = EdgeTracePar(pad=[8, 3])
    assert isinstance(par['pad'], list), 'pad should be stored as a list'
    assert par['pad'] == [8.0, 3.0], 'Asymmetric pad not preserved'

    # A single-element list is expanded to a symmetric 2-element list
    par = EdgeTracePar(pad=[4])
    assert par['pad'] == [4.0, 4.0], 'Single-element pad not expanded correctly'

    # Writing to a FITS header must not raise (arrays are illegal header values)
    hdr = fits.Header()
    par.to_header(hdr)

    # ... and the value round-trips back to the same list
    par_rt = EdgeTracePar.from_header(hdr)
    assert par_rt['pad'] == [4.0, 4.0], 'pad did not round-trip through the header'


def test_pad():
    # SlitTraceSet requires pad as an ndarray (as the consumers now supply via
    # np.asarray) and normalizes it to a 2-element array.

    # A single-element array is expanded to a symmetric 2-element array
    slits = SlitTraceSet(np.full((1000, 3), 2, dtype=float), np.full((1000, 3), 8, dtype=float),
                         'MultiSlit', nspat=10, PYP_SPEC='dummy',
                         pad=np.asarray([3.0], dtype=float))
    assert isinstance(slits.pad, np.ndarray), 'pad should be an ndarray on SlitTraceSet'
    assert np.array_equal(slits.pad, np.array([3.0, 3.0])), 'Single-element pad not expanded correctly'

    # A 2-element list (as delivered by EdgeTracePar) is preserved in order
    slits = SlitTraceSet(np.full((1000, 3), 2, dtype=float), np.full((1000, 3), 8, dtype=float),
                         'MultiSlit', nspat=10, PYP_SPEC='dummy',
                         pad=np.asarray([8.0, 3.0], dtype=float))
    assert np.array_equal(slits.pad, np.array([8.0, 3.0])), 'Asymmetric pad not preserved'

    # pad survives a write/read round-trip to disk
    slits.set_paths(data_output_path(''), 'A', '1', 'DET01')
    ofile = Path(slits.get_path()).absolute()
    slits.to_file(overwrite=True)
    assert ofile.exists(), 'File not written'
    _slits = SlitTraceSet.from_file(ofile)
    assert np.array_equal(_slits.pad, np.array([8.0, 3.0])), 'pad did not round-trip through file'

    ofile.unlink()


def test_slit_img_split_overlap():
    # Three slits: the first two are 3.5 pixels apart, the last two overlap
    nspec, nspat = 20, 60
    left = np.tile([5.5, 21.0, 33.5], (nspec, 1))
    right = np.tile([17.5, 35.5, 50.5], (nspec, 1))
    slits = SlitTraceSet(left, right, 'MultiSlit', nspat=nspat, PYP_SPEC='dummy')
    spat_id = slits.spat_id
    # Without padding, the pixels where the last two slits overlap are
    # assigned to the last slit, unless the overlap is split at its midpoint
    # (x = 34.5)
    nopad = slits.slit_img(pad=0)
    assert np.all(nopad[:, 18:21] == -1), 'Gap pixels should not be assigned'
    assert np.all(nopad[:, 34:36] == spat_id[2]), 'Overlap pixels should be overwritten'
    split = slits.slit_img(pad=0, split_overlap=True)
    assert np.all(split[:, 33:36] == spat_id[[1, 1, 2]]), 'Overlap not split at its midpoint'
    # With padding, the padding of the second slit also covers the right edge
    # of the first slit, unless the overlap is split: then, the pixels on the
    # first slit stay on it, and the gap between the slits (17.5 to 21.0) is
    # split at its midpoint (x = 19.25)
    pad = 5
    padded = slits.slit_img(pad=pad)
    assert np.all(padded[:, 17:23] == spat_id[1]), 'Padded pixels should be overwritten'
    split = slits.slit_img(pad=pad, split_overlap=True)
    assert np.all(split[:, 15:23] == spat_id[[0, 0, 0, 0, 0, 1, 1, 1]]), \
        'Padded overlap not split at the midpoint of the gap'
    assert np.all(split[:, 33:37] == spat_id[[1, 1, 2, 2]]), \
        'Padded overlap not split at the midpoint of the overlapping slits'
    # Pixels in only one (padded) slit are not affected
    nslit = np.sum([slits.slit_img(pad=pad, slitidx=i) != -1 for i in range(slits.nslits)],
                   axis=0)
    assert np.array_equal(split[nslit < 2], padded[nslit < 2]), 'Pixels in a single slit changed'


def test_io_single():
    # NOTE this is just a test string. The file itself is not actually read, so it is not required to run the test.
    file = '/home/xavier/Projects/PypeIt-development-suite/RAW_DATA/keck_deimos/830G_M_8500/DE.20100913.57006.fits.gz'
    slits = SlitTraceSet(np.full((1000, 1), 2, dtype=float), np.full((1000, 1), 8, dtype=float),
                         'MultiSlit',
                         nspat=10, PYP_SPEC='dummy',
                         maskfile=file)
    slits.set_paths(data_output_path(''), 'A', '1', 'DET01')
    ofile = Path(slits.get_path()).absolute()

    # Try to save it
    slits.to_file()

    # And read it back in
    _slits = SlitTraceSet.from_file(ofile)

    assert np.array_equal(_slits.left_init, np.full((1000, 1), 2, dtype=float)), 'Bad left read'
    # And that it's the same as the existing one
    assert np.array_equal(_slits.left_init, slits.left_init), 'Bad left read'

    ofile.unlink()




def _maskinfo_objects(slits, positions, spat_flexure):
    """Detected objects (target and serendipitous source) in the first slit."""
    from pypeit import specobj, specobjs
    sobjs = specobjs.SpecObjs()
    for ii, pos in enumerate(positions):
        sobj = specobj.SpecObj('MultiSlit', 'DET01', SLITID=slits.spat_id[0])
        sobj.OBJID = ii + 1
        sobj.SPAT_PIXPOS = pos + (0. if spat_flexure is None else np.mean(spat_flexure[0]))
        sobj.smash_peakflux = 100.
        sobj.FWHM = 3.
        sobjs.add_sobj(sobj)
    return sobjs


def test_assign_maskinfo_spat_flexure():
    # The RA and Dec assigned to a serendipitous source must not depend on the
    # spatial flexure of the frame
    from astropy.table import Table
    nspec, platescale = 100, 0.1
    left = np.tile([50., 150.], (nspec, 1))
    designtab = Table({'MASKDEF_ID': [1001, 1002], 'SLITRA': [150.0, 150.0],
                       'SLITDEC': [10.0, 10.001], 'SLITPA': [30.0, 30.0],
                       'OBJRA': [150.0001, 150.0002], 'OBJDEC': [10.0001, 10.0011],
                       'OBJNAME': ['targ1', 'targ2'], 'OBJMAG': [20., 21.],
                       'OBJMAG_BAND': ['r', 'r']})
    designtab.meta['MASKRMSL'] = 0.
    slits = SlitTraceSet(left_init=left, right_init=left + 60., pypeline='MultiSlit',
                         detname='DET01', nspat=250, PYP_SPEC='dummy',
                         maskdef_id=np.array([1001, 1002]), maskdef_designtab=designtab,
                         maskdef_objpos=np.array([20., 30.]), maskdef_offset=0.5)
    slits.maskdef_slitcen = slits.center.copy()
    # Target (at its expected position) and a serendipitous source, relative
    # to the slit in the frame
    positions = left[0, 0] + np.array([20.5, 45.])
    radec = []
    for spat_flexure in [None, np.array([[2.5, 2.5], [1., 1.]])]:
        sobjs = _maskinfo_objects(slits, positions, spat_flexure)
        sobjs = slits.assign_maskinfo(sobjs, platescale, spat_flexure)
        assert list(sobjs.MASKDEF_OBJNAME) == ['targ1', 'SERENDIP'], 'Incorrect object matching'
        radec.append((sobjs[1].RA, sobjs[1].DEC))
    assert np.allclose(radec[0], radec[1], rtol=0, atol=1e-9), \
        'The RA/Dec of the serendipitous source depends on the spatial flexure'
