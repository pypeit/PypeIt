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
    shape = (100,100)
    pypeitImage = pypeitimage.PypeItImage(np.ones(shape), ivar=np.ones(shape))
    # Full datamodel
    assert 'detector' in pypeitImage.keys(), 'Detector somehow missing!'

    # I/O
    outfile = data_path('tst_pypeitimage.fits')
    pypeitImage.to_file(outfile, overwrite=True)
    _pypeitImage = pypeitimage.PypeItImage.from_file(outfile)

    # Cleanup
    os.remove(outfile)

    # Test
    assert np.array_equal(_pypeitImage.image, np.ones(shape)), 'image array changed'
    assert np.array_equal(_pypeitImage.ivar, np.ones(shape)), 'ivar array changed'
    assert isinstance(_pypeitImage.fullmask, imagebitmask.ImageBitMaskArray), 'mask type changed'
    assert np.array_equal(_pypeitImage.fullmask.mask, pypeitImage.fullmask.mask), \
                'mask array changed'

def test_sub():
    shape = (10,10)
    # Create two images
    img1 = pypeitimage.PypeItImage(np.ones(shape), ivar=np.ones(shape))
    img2 = pypeitimage.PypeItImage(np.ones(shape), ivar=np.ones(shape))

    # Random number generator (set the seed so that the performance is
    # deterministic)
    rng = np.random.default_rng(99)

    # Select random elements in the 2D array to flag as bad pixels
    indx = np.unravel_index(rng.integers(low=0, high=np.prod(shape), size=20), shape)
    img1.update_mask('BPM', indx=indx)
    # Select a different set for the 2nd image
    indx = np.unravel_index(rng.integers(low=0, high=np.prod(shape), size=20), shape)
    img2.update_mask('BPM', indx=indx)

    # Select random elements in the 2D array to flag as cosmic rays
    indx = np.unravel_index(rng.integers(low=0, high=np.prod(shape), size=20), shape)
    img1.update_mask('CR', indx=indx)
    # Select a different set for the 2nd image
    indx = np.unravel_index(rng.integers(low=0, high=np.prod(shape), size=20), shape)
    img2.update_mask('CR', indx=indx)

    diff = img1 - img2

    assert np.array_equal(diff.image, np.zeros(shape)), 'Bad subtraction'
    assert np.array_equal(diff.ivar, np.full(shape, 0.5)), 'Bad error propagation'
    assert np.array_equal(diff.fullmask.bpm, img1.fullmask.bpm | img2.fullmask.bpm), \
                'Bad BPM propagation'
    assert np.array_equal(diff.fullmask.cr, img1.fullmask.cr | img2.fullmask.cr), \
                'Bad CR propagation'


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

    # Random number generator (set the seed so that the performance is
    # deterministic)
    rng = np.random.default_rng(99)

    # Select random elements in the 2D array
    bpm_indx = np.unravel_index(rng.integers(low=0, high=np.prod(shape), size=20), shape)
    # Flag them as 'BPM'
    mask.turn_on('BPM', select=bpm_indx)
    # Check they were flagged correctly
    bpm_mask = np.zeros(shape, dtype=bool)
    bpm_mask[bpm_indx] = True
    assert np.array_equal(bpm_mask, mask.flagged(flag='BPM')), 'Bad BPM flagging'
    # Check the convenience functionality
    assert np.array_equal(mask.bpm, mask.flagged(flag='BPM')), 'Bad convenience access to BPM'

    # Add some cosmic-ray hits
    cr_indx = np.unravel_index(rng.integers(low=0, high=np.prod(shape), size=20), shape)
    mask.turn_on('CR', select=cr_indx)
    # Check they were flagged correctly
    cr_mask = np.zeros(shape, dtype=bool)
    cr_mask[cr_indx] = True
    assert np.array_equal(cr_mask, mask.cr), 'Bad CR flagging'

    # Make sure the bpm mask didn't change
    assert np.array_equal(bpm_mask, mask.bpm), 'BPM should not change'

    # Check combined masking
    assert np.array_equal(bpm_mask | cr_mask, mask.flagged(flag=['BPM', 'CR'])), \
            'Combined masking is wrong'

    # Flag everything as saturated
    mask.turn_on('SATURATION')
    assert np.all(mask.saturation), 'All should have been flagged as saturated.'


def test_bitmaskarray_io():
    path = Path(data_path('test.fits')).resolve()
    if path.exists():
        path.unlink()

    # Create a new mask, flag some bits, and write it
    shape = (10,10)
    mask = imagebitmask.ImageBitMaskArray(shape)
    # Random number generator (set the seed so that the performance is
    # deterministic)
    rng = np.random.default_rng(99)
    bpm_indx = np.unravel_index(rng.integers(low=0, high=np.prod(shape), size=20), shape)
    mask.turn_on('BPM', select=bpm_indx)
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

