"""
Module to run tests on Alignment frames
"""
from pathlib import Path

from IPython import embed

import numpy as np

from astropy.io import fits

from pypeit import alignframe
from pypeit.tests.tstutils import data_output_path


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
    alignments.set_paths(data_output_path(''), 'A', '1', 'DET01')
    ofile = Path(alignments.get_path()).absolute()

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



def test_alignments_spat_flexure(tmp_path):
    # The spatial flexure of the alignment frame is written and read
    nspat, nspec, nslits, nalign = 100, 200, 3, 5
    spat_flexure = np.array([[0.5, 0.7], [1.2, 1.3], [-0.4, 0.1]])
    alignments = alignframe.Alignments(alignframe=np.ones((nspec, nspat)), nspec=nspec,
                                       nalign=nalign, nslits=nslits,
                                       traces=np.zeros((nspec, nalign, nslits)),
                                       PYP_SPEC='keck_kcwi', spat_id=np.arange(nslits),
                                       spat_flexure=spat_flexure)
    alignments.set_paths(tmp_path, 'A', '1', 'DET01')
    alignments.to_file(overwrite=True)
    _alignments = alignframe.Alignments.from_file(alignments.get_path())
    assert np.allclose(_alignments.spat_flexure, spat_flexure), 'Spatial flexure changed on I/O'
    # The spatial flexure is optional
    alignments.spat_flexure = None
    alignments.to_file(overwrite=True)
    assert alignframe.Alignments.from_file(alignments.get_path()).spat_flexure is None
