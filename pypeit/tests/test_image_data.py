"""
Tests for the image preparation functions of the ginga display scripts.
"""
from astropy.io import fits
import numpy as np
import pytest

from pypeit import PypeItError
from pypeit.display import image_data


def test_pixflat_images(tmp_path):
    ofile = tmp_path / 'pixelflat_test.fits'
    hdus = [fits.PrimaryHDU()]
    for det in (1, 2):
        hdus.append(fits.ImageHDU(data=np.full((4,3), float(det)),
                                  name=f'DET0{det}-PIXELFLAT_NORM'))
    fits.HDUList(hdus).writeto(ofile)
    images = image_data.pixflat_images(ofile)
    assert [name for _, name in images] == ['DET01-PIXELFLAT_NORM', 'DET02-PIXELFLAT_NORM'], \
        'All detectors should be read by default'
    images = image_data.pixflat_images(ofile, det=[2, 5])
    assert len(images) == 1 and images[0][0][0,0] == 2., \
        'Only the requested detectors in the file should be read'
    with pytest.raises(PypeItError, match='not found'):
        image_data.pixflat_images(ofile, det=[7])
