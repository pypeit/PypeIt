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
            np.allclose(np.diag(np.fliplr(tform[:2,:2])), np.sqrt(0.5)*np.array([-0.5,2])), \
            'Bad scale transform'

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
    """
    dragons_file = data_path('GN_HAM_R400_885_N20190205S0035_dragons_mosaic.fits')
    file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'gemini_gmos', 'GN_HAM_R400_885',
                        'N20190205S0035.fits')

    # Load the spectrograph
    spec = load_spectrograph('gemini_gmos_north_ham')
    msc = (1,2,3)
    msc_par = spec.get_mosaic_par(msc, hdu=fits.open(file))

#    expected_shape, shift, rot = spec.get_mosaic_par(msc)
    # Load the images and trim and orient them
    imgs = [None]*spec.ndet
#    tforms = [None]*spec.ndet
    for i, det in enumerate(msc):
        imgs[i] = RawImage(file, spec, det)
        imgs[i].trim()
        imgs[i].orient()
        imgs[i] = imgs[i].image
#        binning = tuple(int(b) for b in imgs[det].detector.binning.split(','))
#        shape = tuple(n // b for n, b in zip(expected_shape, binning))
#        tforms[det] = mosaic.build_image_mosaic_transform(shape, shift[det], rot[det], binning)

#    input_imgs = [img.image for img in imgs]
#    mosaic_pypeit, mosaic_ivar, mosaic_npix = mosaic.build_image_mosaic(input_imgs, tforms)

    mosaic_pypeit, mosaic_ivar, mosaic_npix = mosaic.build_image_mosaic(imgs, msc_par.tform)

    _mosaic_pypeit = np.fliplr(mosaic_pypeit.T).astype(np.uint16)
    mosaic_dragons = fits.open(dragons_file)[0].data
    assert np.array_equal(mosaic_dragons, _mosaic_pypeit[1:-1,1:-1]), 'Bad mosaic'


test_gemini_gmos()


