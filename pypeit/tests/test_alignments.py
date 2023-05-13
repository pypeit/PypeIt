"""
Module to run tests on Alignment frames
"""
from pathlib import Path

from IPython import embed

import numpy as np

from astropy.io import fits

from pypeit import alignframe
from pypeit.tests.tstutils import data_path


def test_alignments():
    nspat, nspec, nslits, nalign = 100, 1000, 10, 5
    tmp = np.ones((nspec, nspat)) * 10.
    traces = np.zeros((nspec, nalign, nslits))
    instant_dict = dict(alignframe=tmp,
                        nspec=nspec,
                        nalign=nalign,
                        nslits=nslits,
                        traces=traces,
                        PYP_SPEC='keck_kcwi',
                        spat_id=np.arange(nslits))

    alignments = alignframe.Alignments(**instant_dict)
    alignments.set_paths(data_path(''), 'A', '1', 'DET01')
    ofile = Path(alignments.get_path()).resolve()

    # I/O
    alignments.to_file(overwrite=True)
    assert ofile.exists(), 'File not written'
    with fits.open(ofile) as hdu:
        assert hdu[1].name == 'ALIGN', 'Extension name changed'
    _alignments = alignframe.Alignments.from_file(str(ofile))

    # Test
    for key in instant_dict.keys():
        if isinstance(instant_dict[key], np.ndarray):
            assert np.array_equal(alignments[key], _alignments[key])
        else:
            assert alignments[key] == _alignments[key]

    # Clean-up
    ofile.unlink()

