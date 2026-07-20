"""
Module to run tests on FlatField class
"""
from pathlib import Path

from IPython import embed

import numpy as np

from astropy.io import fits

from pypeit import flatfield
from pypeit.core import bspline
from pypeit.spectrographs.util import load_spectrograph
from pypeit.tests.tstutils import data_output_path


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
    flatImages.set_paths(data_output_path(''), 'A', '1', 'DET01')

    # I/O
    ofile = Path(flatImages.get_path()).absolute()
    flatImages.to_file(overwrite=True)
    assert ofile.exists(), 'File not written'

    _flatImages = flatfield.FlatImages.from_file(ofile)

    # Test
    for key in instant_dict.keys():
        if key == 'pixelflat_spat_bsplines':
            assert np.array_equal(flatImages[key][0].breakpoints,
                                  _flatImages[key][0].breakpoints), 'pixelflat breakpoints changed'
            continue
        if key == 'illumflat_spat_bsplines':
            assert np.array_equal(flatImages[key][0].breakpoints,
                                  _flatImages[key][0].breakpoints), 'illumflat breakpoints changed'
            continue
        if isinstance(instant_dict[key], np.ndarray):
            assert np.array_equal(flatImages[key],_flatImages[key])
        else:
            assert flatImages[key] == _flatImages[key]

    ofile.unlink()


def test_fit_det_response():
    spec = load_spectrograph('keck_kcwi')
    # Generate a good pixel mask
    frsize = 4100
    gpm = np.ones((frsize,frsize), dtype=bool)
    # Generate a fake image
    sinemodel = lambda xx, yy, amp, scl, phase, wavelength, angle: 1 + (amp + xx * scl) * np.sin(
                2 * np.pi * (xx * np.cos(angle*np.pi / 180.0) + yy * np.sin(angle*np.pi / 180.0)) / wavelength + phase)
    x = np.arange(frsize)
    y = np.arange(frsize)
    xx, yy = np.meshgrid(x, y, indexing='ij')
    amp, scale, wavelength, phase, angle = 0.02, 0.0, 1.41*frsize/31.5, 0.0, -45.34
    img = sinemodel(xx, yy, amp, scale, phase, wavelength, angle)
    model = spec.fit_2d_det_response(img, gpm)
    assert np.allclose(img, model, atol=0.001), 'structure fitting failed.'


def test_load_pixflat_preserves_calib_key(tmp_path):
    """
    ``load_pixflat`` merges a user/archived pixel flat into the existing
    ``FlatImages``.  The merge rebuilds the object from datamodel fields
    only, so the calibration identifiers must be carried over explicitly;
    otherwise ``calib_key`` is None and ``get_path()`` raises (this broke
    the flats state recording for keck_lris_blue, which uses an archived
    pixel flat).
    """
    spec = load_spectrograph('keck_lris_blue')
    cal_dir = str(tmp_path)
    # Minimal archived pixel-flat file: one extension per detector, named
    # '<detname>-PIXELFLAT_NORM', as load_pixflat expects.
    pixflat_file = str(tmp_path / 'mypixflat.fits')
    fits.HDUList([
        fits.PrimaryHDU(),
        fits.ImageHDU(np.full((10, 10), 0.5, dtype=float),
                      name='DET01-PIXELFLAT_NORM')]).writeto(pixflat_file,
                                                             overwrite=True)
    # An existing FlatImages with calibration identifiers set
    flats = flatfield.FlatImages(pixelflat_norm=np.ones((10, 10)))
    flats.calib_key = 'B_1_DET01'
    flats.calib_dir = cal_dir
    merged = flatfield.load_pixflat(pixflat_file, spec, 1, flats,
                                    calib_dir=cal_dir, chk_version=False)
    # The merged result keeps the user pixel flat...
    assert np.allclose(merged.pixelflat_norm, 0.5)
    # ...and, crucially, the calibration identifiers so get_path() works.
    assert merged.calib_key == 'B_1_DET01'
    assert merged.calib_dir == cal_dir
    assert Path(merged.get_path()).name == 'Flat_B_1_DET01.fits'


