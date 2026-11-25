"""
Tests for the image combination in :class:`~pypeit.images.combineimage.CombineImage`.
"""
import numpy as np

from pypeit.core import procimg
from pypeit.images import detector_container
from pypeit.images import pypeitimage
from pypeit.images.combineimage import CombineImage
from pypeit.par import pypeitpar
from pypeit.tests import tstutils


def _scaled_image(rng, shape, rn2, scale):
    """
    Build a processed image whose variance arrays follow the conventions of
    :func:`~pypeit.images.rawimage.RawImage.process`: ``rn2img`` in detector
    units and ``base_var`` including the count scaling.
    """
    rn2img = np.full(shape, rn2, dtype=float)
    img_scale = np.full(shape, scale, dtype=float)
    base_var = procimg.base_variance(rn2img, count_scale=img_scale)
    image = scale * rng.normal(loc=100., scale=1., size=shape)
    img = pypeitimage.PypeItImage(
        image, ivar=np.ones(shape, dtype=float), rn2img=rn2img, base_var=base_var,
        img_scale=img_scale, exptime=10., noise_floor=0., shot_noise=False,
        amp_img=np.ones(shape, dtype=int),
        detector=detector_container.DetectorContainer(**tstutils.default_detector())
    )
    img.filename = 'test.fits'
    return img


def test_combine_base_var_scaling():
    rng = np.random.default_rng(seed=42)
    shape = (20, 10)
    rn2 = 9.
    scale = 2.
    nimg = 3
    images = [_scaled_image(rng, shape, rn2, scale) for _ in range(nimg)]

    par = pypeitpar.ProcessImagesPar(
        combine='mean', clip=False, shot_noise=False, noise_floor=0.
    )
    comb = CombineImage(images, par).run()

    # The mean of N frames each with base variance s^2 V_rn has base variance
    # s^2 V_rn / N.
    expected = scale**2 * rn2 / nimg
    assert np.allclose(comb.base_var, expected), \
        'Combined base variance should be s^2 V_rn / N, without double scaling.'
    assert np.allclose(comb.rn2img, expected), \
        'Combined readnoise variance should be s^2 V_rn / N.'
    assert np.allclose(comb.img_scale, scale / nimg), \
        'Combined image scaling should be the mean scaling divided by N.'
    assert np.allclose(comb.ivar, 1 / expected), \
        'Without shot noise, the combined inverse variance should be 1 / base_var.'
