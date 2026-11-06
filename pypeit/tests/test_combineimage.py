"""
Tests for the CombineImage class
"""
import numpy as np

from pypeit.images import pypeitimage, combineimage
from pypeit.par.pypeitpar import ProcessImagesPar
from pypeit.spectrographs.util import load_spectrograph


def _processed_image(spat_flexure, rng, detector):
    """Build a minimal processed image with the provided spatial flexure."""
    img = 100. + rng.normal(0, 1, (20, 30))
    return pypeitimage.PypeItImage(img, ivar=np.ones_like(img), rn2img=np.ones_like(img),
                                   base_var=np.ones_like(img), img_scale=np.ones_like(img),
                                   amp_img=np.ones(img.shape, dtype=int),
                                   bpm=np.zeros(img.shape, dtype=bool), detector=detector,
                                   spat_flexure=spat_flexure, exptime=10., units='e-',
                                   noise_floor=0.0, shot_noise=False)


def test_combine_spat_flexure():
    rng = np.random.default_rng(1)
    detector = load_spectrograph('shane_kast_blue').get_detector_par(1)
    par = ProcessImagesPar(combine='mean')
    flex = np.array([[1.0, 2.0], [3.0, 4.0]])

    # The combined spatial flexure is the mean of the input spatial flexures
    images = [_processed_image(f, rng, detector) for f in [flex, 3*flex]]
    comb = combineimage.CombineImage(images, par).run()
    assert np.allclose(comb.spat_flexure, 2*flex), 'Incorrect combined spatial flexure'

    # Images without a spatial flexure (e.g., if a manual spatial flexure was
    # only provided for some of the frames) are assumed to have no flexure
    images = [_processed_image(f, rng, detector) for f in [flex, None, flex]]
    comb = combineimage.CombineImage(images, par).run()
    assert comb.spat_flexure is not None, 'Partial spatial flexure values should not be discarded'
    assert np.allclose(comb.spat_flexure, 2*flex/3), 'Incorrect combined spatial flexure'

    # No spatial flexure for any of the images
    images = [_processed_image(None, rng, detector) for _ in range(2)]
    comb = combineimage.CombineImage(images, par).run()
    assert comb.spat_flexure is None, 'Combined spatial flexure should be None'
