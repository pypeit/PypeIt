"""
Module to run tests on FlatField class
"""
import os

import pytest

from IPython import embed

import numpy as np


from pypeit import flatfield
from pypeit import bspline

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_flatimages():
    tmp = np.ones((1000, 100)) * 10.
    x = np.random.rand(500)
    # Create bspline
    spat_bspline1 = bspline.bspline(x, bkspace=0.01*(np.max(x)-np.min(x)))
    spat_bspline2 = bspline.bspline(x, bkspace=0.01*(np.max(x)-np.min(x)))
    instant_dict = dict(pixelflat_raw=tmp,
                        pixelflat_norm=np.ones_like(tmp),
                        pixelflat_model=None,
                        pixelflat_spat_bsplines=np.asarray([spat_bspline1, spat_bspline2]),
                        pixelflat_spec_illum=None,
                        illumflat_raw=tmp,
                        illumflat_spat_bsplines=np.asarray([spat_bspline1, spat_bspline2]),
                        spat_id=np.asarray([100, 200]),
                        PYP_SPEC="specname")

    flatImages = flatfield.FlatImages(**instant_dict)
    assert flatImages.pixelflat_model is None
    assert flatImages.pixelflat_spec_illum is None
    assert flatImages.pixelflat_spat_bsplines is not None

    # I/O
    outfile = data_path('tst_flatimages.fits')
    flatImages.to_master_file(outfile)
    _flatImages = flatfield.FlatImages.from_file(outfile)

    # Test
    for key in instant_dict.keys():
        if key == 'pixelflat_spat_bsplines':
            np.array_equal(flatImages[key][0].breakpoints,
                           _flatImages[key][0].breakpoints)
            continue
        if key == 'illumflat_spat_bsplines':
            np.array_equal(flatImages[key][0].breakpoints,
                           _flatImages[key][0].breakpoints)
            continue
        if isinstance(instant_dict[key], np.ndarray):
            assert np.array_equal(flatImages[key],_flatImages[key])
        else:
            assert flatImages[key] == _flatImages[key]

    os.remove(outfile)
