"""
Module to run tests on PypeItImage class
"""
from pathlib import Path
import os

from IPython import embed

import numpy as np

from astropy.io import fits

from pypeit.images import pypeitimage
from pypeit.images import imagebitmask
from pypeit.tests.tstutils import data_path


def test_full():
    pypeitImage = pypeitimage.PypeItImage(np.ones((1000, 1000)))
    pypeitImage.reinit_mask()
    # Full datamodel
#    full_datamodel = pypeitImage.full_datamodel()
#    assert 'gain' in full_datamodel.keys()
    assert 'detector' in pypeitImage.keys()

    # I/O
    outfile = data_path('tst_pypeitimage.fits')
    pypeitImage.to_file(outfile, overwrite=True)
    _pypeitImage = pypeitimage.PypeItImage.from_file(outfile)

    # Cleanup
    os.remove(outfile)

    # Test
    assert isinstance(_pypeitImage.image, np.ndarray)
    assert _pypeitImage.ivar is None


def test_bitmask():
    bm = imagebitmask.ImageBitMask()

    # Check that a few of the bit numbers are correct (i.e., that the order
    # didn't get messed up)
    assert bm.bits['BPM'] == 0, 'BPM bit number changed'
    assert bm.bits['CR'] == 1, 'CR bit number changed'
    assert bm.bits['OFFSLITS'] == 4, 'OFFSLITS bit number changed'
    assert bm.bits['EXTRACT'] == 8, 'EXTRACT bit number changed'


def test_bitmaskarray():
    shape = (10,10)
    mask = imagebitmask.ImageBitMaskArray(shape)

    # Check the shape
    assert mask.shape == shape, 'Shape is incorrect'

    # Check the convenience flag access
    assert isinstance(mask.bpm, np.ndarray), 'Bad bpm property'
    # Check all are initiated to False
    assert np.array_equal(np.zeros(shape, dtype=bool), mask.bpm), \
            'Bad instantiation; should all be false'

    # Random number generator
    rng = np.random.default_rng()

    # Select random elements in the 2D array
    bpm_indx = np.unravel_index(rng.integers(low=0, high=np.prod(shape), size=20), shape)
    # Flag them as 'BPM'
    mask.turn_on(bpm_indx, 'BPM')
    # Check they were flagged correctly
    bpm_mask = np.zeros(shape, dtype=bool)
    bpm_mask[bpm_indx] = True
    assert np.array_equal(bpm_mask, mask.flagged(flag='BPM')), 'Bad BPM flagging'
    # Check the convenience functionality
    assert np.array_equal(mask.bpm, mask.flagged(flag='BPM')), 'Bad convenience access to BPM'

    # Add some cosmic-ray hits
    cr_indx = np.unravel_index(rng.integers(low=0, high=np.prod(shape), size=20), shape)
    mask.turn_on(cr_indx, 'CR')
    # Check they were flagged correctly
    cr_mask = np.zeros(shape, dtype=bool)
    cr_mask[cr_indx] = True
    assert np.array_equal(cr_mask, mask.cr), 'Bad CR flagging'

    # Make sure the bpm mask didn't change
    assert np.array_equal(bpm_mask, mask.bpm), 'BPM should not change'

    # Check combined masking
    assert np.array_equal(bpm_mask | cr_mask, mask.flagged(flag=['BPM', 'CR'])), \
            'Combined masking is wrong'


def test_bitmaskarray_io():
    path = Path(data_path('test.fits')).resolve()
    if path.exists():
        path.unlink()

    # Create a new mask, flag some bits, and write it
    shape = (10,10)
    mask = imagebitmask.ImageBitMaskArray(shape)
    rng = np.random.default_rng()
    bpm_indx = np.unravel_index(rng.integers(low=0, high=np.prod(shape), size=20), shape)
    mask.turn_on(bpm_indx, 'BPM')
    mask.to_file(str(path))

    # Open it directly
    with fits.open(path) as hdu:
        # Check the number of extensions
        assert len(hdu) == 2, 'Should be two extensions'
        # Check the type
        assert hdu[1].data.dtype.type == mask.bitmask.minimum_dtype(), 'Type mismatch'
        # Check the data
        assert np.array_equal(hdu[1].data, mask.mask), 'Writing changed the data'

    # Read using the DataContainer methods
    _mask = imagebitmask.ImageBitMaskArray.from_file(str(path))
    # Check the data
    assert np.array_equal(_mask.mask, mask.mask), 'Read in data is wrong'

    # Clean-up
    path.unlink()

