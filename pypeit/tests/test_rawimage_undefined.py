"""
Tests of the handling of FITS BLANK values and undefined (non-finite) raw pixel
values.
"""
from functools import partial

import numpy as np
import pytest

from astropy.io import fits

from pypeit import dataPaths, io
from pypeit.images import combineimage, pypeitimage, rawimage
from pypeit.spectrographs.util import load_spectrograph

# Raw (untrimmed, unoriented) pixel coordinates in the Kast test frame (350 x
# 2112).  Amplifier 1 data occupy columns 0-1023, and its overscan occupies
# columns 2049-2079.
DATA_PIX = (100, 500)
OSCAN_PIX = (120, 2060)


def _kast_par():
    """Return the arc processing parameters for Kast, without bias subtraction."""
    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.default_pypeit_par()['calibrations']['arcframe']['process']
    par['use_biasimage'] = False
    return spectrograph, par


def _get_rawimage_with_nan(get_rawimage, pixels, raw_file, det):
    """Call a spectrograph's ``get_rawimage`` and set the selected pixels to NaN."""
    detector, raw, hdu, exptime, datasec, oscansec = get_rawimage(raw_file, det)
    raw = raw.copy()
    for pix in pixels:
        raw[pix] = np.nan
    return detector, raw, hdu, exptime, datasec, oscansec


def _process(spectrograph, par, monkeypatch, pixels):
    """Process the Kast test frame after setting the selected raw pixels to NaN."""
    with monkeypatch.context() as m:
        m.setattr(spectrograph, 'get_rawimage',
                  partial(_get_rawimage_with_nan, spectrograph.get_rawimage, pixels))
        raw = rawimage.RawImage(str(dataPaths.tests.get_file_path('b24.fits.gz')),
                                spectrograph, 1)
    img = raw.process(par)
    raw.hdu.close()
    return raw, img


def test_undefined_raw_pixels(monkeypatch):
    spectrograph, par = _kast_par()
    ref_raw, ref = _process(spectrograph, par, monkeypatch, [])
    raw, img = _process(spectrograph, par, monkeypatch, [DATA_PIX, OSCAN_PIX])

    assert not np.any(ref_raw.raw_undef), 'Reference frame should not have undefined pixels'
    assert np.all(np.isfinite(img.image)), 'Undefined raw pixels propagated to processed image'

    # Only the data pixel survives trimming; it must be flagged as BPM and
    # IS_NAN
    is_nan = img.select_flag(flag='IS_NAN')
    assert np.sum(is_nan) == 1, 'Undefined data pixel not flagged'
    assert np.array_equal(is_nan, raw.raw_undef[0]), 'IS_NAN does not match the processed mask'
    assert img.select_flag(flag='BPM')[is_nan].all(), 'Undefined pixel should also be BPM'
    assert not np.any(ref.select_flag(flag='IS_NAN')), 'Reference frame should not be flagged'

    # The undefined overscan pixel should have a negligible effect on the
    # overscan subtraction (one of 31 pixels in its row is ignored).
    gpm = np.logical_not(is_nan)
    assert np.allclose(img.image[gpm], ref.image[gpm], atol=1.), \
            'Undefined overscan pixel significantly altered the overscan subtraction'


@pytest.fixture
def kast_blank_file(tmp_path):
    """
    Copy of the Kast frame (pseudo-unsigned 16-bit) with a FITS BLANK card.
    Astropy >= 8.0.3 converts BLANK pixels to NaN unless ``ignore_blank=True``.
    """
    ofile = tmp_path / 'b24_blank.fits'
    with fits.open(dataPaths.tests.get_file_path('b24.fits.gz'),
                   do_not_scale_image_data=True) as hdu:
        hdu[0].header['BLANK'] = -32768
        # Stored value of -32768 is 0 after applying BZERO = 32768
        hdu[0].data[DATA_PIX] = -32768
        hdu[0].data[OSCAN_PIX] = -32768
        hdu.writeto(ofile)
    return ofile


def test_fits_open_ignores_blank(kast_blank_file):
    with io.fits_open(kast_blank_file) as hdu:
        data = hdu[0].data
        assert not np.any(np.isnan(data)), 'BLANK pixels should not be converted to NaN'
        assert data[DATA_PIX] == 0 and data[OSCAN_PIX] == 0, \
                'BLANK pixels should be read as their stored (scaled) value'


def test_blank_card(kast_blank_file):
    # BLANK pixels are read as 0 counts (independent of the astropy version),
    # so they are neither undefined nor flagged.
    spectrograph, par = _kast_par()
    raw = rawimage.RawImage(str(kast_blank_file), spectrograph, 1)
    img = raw.process(par)
    raw.hdu.close()

    assert not np.any(raw.raw_undef), 'BLANK pixels should not be undefined'
    assert np.all(np.isfinite(img.image)), 'BLANK pixels propagated to processed image'
    assert not np.any(img.select_flag(flag='IS_NAN')), 'BLANK pixels should not be flagged'


def test_build_mask_keeps_is_nan():
    img = pypeitimage.PypeItImage(np.ones((5, 5), dtype=float),
                                  ivar=np.ones((5, 5), dtype=float))
    indx = np.zeros((5, 5), dtype=bool)
    indx[2, 3] = True
    img.update_mask('IS_NAN', indx=indx)
    img.update_mask('SATURATION', indx=indx)
    img.build_mask()
    assert np.array_equal(img.select_flag(flag='IS_NAN'), indx), 'IS_NAN not preserved'
    assert not np.any(img.select_flag(flag='SATURATION')), 'Other flags should be erased'


@pytest.mark.parametrize('method', ['mean', 'median'])
def test_combine_undefined(monkeypatch, method):
    spectrograph, par = _kast_par()
    par['combine'] = method
    # The data pixel is undefined in both frames; a second pixel is undefined
    # in only one frame.
    other_pix = (200, 700)
    raw1, img1 = _process(spectrograph, par, monkeypatch, [DATA_PIX, other_pix])
    raw2, img2 = _process(spectrograph, par, monkeypatch, [DATA_PIX])

    comb = combineimage.CombineImage([img1, img2], par).run()
    assert np.all(np.isfinite(comb.image)), 'Undefined pixels propagated to combined image'
    is_nan = comb.select_flag(flag='IS_NAN')
    assert np.array_equal(is_nan, raw2.raw_undef[0]), \
            'Only pixels undefined in all frames should be flagged as IS_NAN'
    assert comb.select_flag(flag='BPM')[is_nan].all(), 'Undefined pixel should also be BPM'
    only_one = raw1.raw_undef[0] & np.logical_not(raw2.raw_undef[0])
    assert np.sum(only_one) == 1, 'Expected one pixel undefined in a single frame'
    assert not comb.select_flag(flag='IS_NAN')[only_one].any() \
            and not comb.select_flag(flag='BPM')[only_one].any(), \
            'Pixel undefined in only one frame should be filled by the other frame'
