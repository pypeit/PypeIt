"""
Module to run tests on Alignment frames
"""
import os
import pytest
import numpy as np

from pypeit import alignframe


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


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

    # I/O
    outfile = data_path('tst_alignments.fits')
    alignments.to_file(outfile, overwrite=True)
    _alignments = alignframe.Alignments.from_file(outfile)
    # Test
    for key in instant_dict.keys():
        if isinstance(instant_dict[key], np.ndarray):
            assert np.array_equal(alignments[key], _alignments[key])
        else:
            assert alignments[key] == _alignments[key]
    # Now delete
    os.remove(outfile)
