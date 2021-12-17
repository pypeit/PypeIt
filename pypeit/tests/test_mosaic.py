"""
Module to run tests on pypeit.core.mosaic methods.
"""
import os
from IPython import embed

import numpy as np
from astropy.io import fits

from pypeit.spectrographs.util import load_spectrograph
from pypeit.images.rawimage import RawImage
from pypeit.tests.tstutils import dev_suite_required, data_path
from pypeit.core import mosaic
from pypeit.images.mosaic import Mosaic


def test_transform():
    shape = (256,512)
    tform = mosaic.build_image_mosaic_transform(shape, (10.,10.), 0., (1,1))
    expected = np.array([[1., 0., 10.],
                         [0., 1., 10.],
                         [0., 0., 1.]])
    assert np.array_equal(tform, expected), 'Transform is incorrect.'

    tform = mosaic.build_image_mosaic_transform(shape, (0.,0.), 45., (1,1))
    assert np.allclose(np.absolute(tform[:2,:2].ravel()), np.repeat([np.sqrt(0.5)], 4)), \
            'Rotation transform is incorrect.'

    shape = (128, 512)
    tform = mosaic.build_image_mosaic_transform(shape, (0.,0.), 45., (2,1))

    assert np.allclose(np.diag(tform)[:2], np.repeat([np.sqrt(0.5)], 2)) and \
            np.allclose(np.diag(np.fliplr(tform[:2,:2])), np.sqrt(0.5)*np.array([0.5,-2])), \
            'Bad scale transform'

@dev_suite_required
def test_io():
    file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'gemini_gmos', 'GN_HAM_R400_885',
                        'N20190205S0035.fits')

    # Load the spectrograph
    spec = load_spectrograph('gemini_gmos_north_ham')
    msc = spec.get_mosaic_par((1,2,3), hdu=fits.open(file))

    test_file = data_path('tmp_mosaic.fits')
    msc.to_file(test_file, overwrite=True)

    _msc = Mosaic.from_file(test_file)

    # Check a few attributes are equal
    assert np.array_equal(msc.tform, _msc.tform), 'Bad transform read'
    assert _msc.shape == msc.shape, 'Bad shape read'
    assert _msc.name == msc.name, 'Bad name setup'

    assert len(_msc.detectors) == len(msc.detectors), 'Bad number of detectors'

    assert _msc.detectors[0]['dataext'] == msc.detectors[0]['dataext'], 'Bad read dataext'
    assert np.array_equal(_msc.detectors[1]['gain'], msc.detectors[1]['gain']), 'Bad read gain'
    assert _msc.detectors[2]['binning'] == msc.detectors[2]['binning'], 'Bad read binning'
    assert np.array_equal(_msc.detectors[1]['datasec'], msc.detectors[1]['datasec']), \
            'Bad read datasec'

    os.remove(test_file)


@dev_suite_required
def test_gemini_gmos():
    """
    The DRAGONS output was constructed as follows:

    .. code-block:: python

        from astropy.io import fits

        import astrodata
        from geminidr.gmos.lookups import geometry_conf
        from geminidr.gemini.lookups.keyword_comments import keyword_comments
        from gempy.gemini import gemini_tools
        from gempy.library import transform

        # Read in the data
        ad = astrodata.open(file)

        # There must be a better way to do this, but this converts the type
        # to float
        ad = ad * 1.

        # Trim it
        gemini_tools.trim_to_data_section(adinput=ad, keyword_comments=keyword_comments)

        # Add the transformation information
        transform.add_mosaic_wcs(ad, geometry_conf)

        # And perform the transformations to create the mosaic image.
        mosaic_dragons = transform.resample_from_wcs(ad, 'mosaic', order=0)[0].data

        fits.HDUList([fits.PrimaryHDU(data=mosaic_dragons.astype(np.uint16))]
                     ).writeto(dragons_file, overwrite=True)

    It's worth noting that this reproduces the DRAGONS result exactly
    *specifically for this case*, but that's not true for some of the other
    setups.  E.g., for some reason, the red chip in the Gemini South data sets
    in the dev-suite is slightly off; the mosaic images are within ~1e-10 when
    the order is >0, but there's different behavior in the nearest-grid-point
    result leading to significant differences.  We're using exactly the same
    transformations, so I don't understand why this should happen.
    """
    dragons_file = data_path('GN_HAM_R400_885_N20190205S0035_dragons_mosaic.fits')
    file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'gemini_gmos', 'GN_HAM_R400_885',
                        'N20190205S0035.fits')

    # Load the spectrograph
    spec = load_spectrograph('gemini_gmos_north_ham')
    msc = (1,2,3)
    msc_par = spec.get_mosaic_par(msc, hdu=fits.open(file))

    # Load the images and trim and orient them
    imgs = [None]*spec.ndet
    for i, det in enumerate(msc):
        imgs[i] = RawImage(file, spec, det)
        imgs[i].trim()
        imgs[i].orient()
        imgs[i] = imgs[i].image[0]

    mosaic_pypeit, mosaic_ivar, mosaic_npix, _ = mosaic.build_image_mosaic(imgs, msc_par.tform)

    _mosaic_pypeit = np.fliplr(mosaic_pypeit.T).astype(np.uint16)
    mosaic_dragons = fits.open(dragons_file)[0].data

    assert np.array_equal(mosaic_dragons, _mosaic_pypeit[1:,:-2]), 'Bad mosaic'


