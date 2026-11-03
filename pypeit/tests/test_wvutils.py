
from IPython import embed
import numpy as np

from pypeit.core.wavecal import waveio
from pypeit.core.wavecal import wvutils

def test_arc_lines_from_spec():

    sigdetect = 10.
    nonlinear_counts = 8e4
    fwhm = 3.
    reid_arxiv = 'keck_deimos_830G.fits'
    det = 3
    wave_start = 6522.7

    wave, flux, *_ = waveio.load_template(reid_arxiv, det)
    s = np.argmin(np.absolute(wave - wave_start))
    assert s == 4186, 'Template changed'

    wave = wave[s:s+4096]
    flux = flux[s:s+4096]

    tcent, ecent, cut_tcent, icut, spec_cont_sub = wvutils.arc_lines_from_spec(
        flux[:2048], sigdetect=sigdetect, nonlinear_counts=nonlinear_counts, fwhm=fwhm,
    )
    assert len(tcent) == 0, 'should initially find no good lines'

    tcent, ecent, cut_tcent, icut, spec_cont_sub = wvutils.arc_lines_from_spec(
        flux[:2048], sigdetect=sigdetect, nonlinear_counts=nonlinear_counts, fwhm=fwhm,
        good_frac=0.5, fwhm_incr=1.1, max_good_iter=5,
    )
    assert len(tcent) == 16, 'iteratively increasing the FWHM should yield many lines'

    tcent, ecent, cut_tcent, icut, spec_cont_sub = wvutils.arc_lines_from_spec(
        flux[2048:], sigdetect=sigdetect, nonlinear_counts=nonlinear_counts, fwhm=fwhm,
    )
    assert len(tcent) == 0, 'should initially find no good lines'

    tcent, ecent, cut_tcent, icut, spec_cont_sub = wvutils.arc_lines_from_spec(
        flux[2048:], sigdetect=sigdetect, nonlinear_counts=nonlinear_counts, fwhm=fwhm,
        good_frac=0.5, fwhm_incr=1.1, max_good_iter=5,
    )
    assert len(tcent) == 30, 'iteratively increasing the FWHM should yield many lines'


def _synthetic_arc_pair(shift, stretch, nspec=1024, fwhm=3.0, nlines=25, seed=42):
    """
    Build a reference arc and an "observed" arc related by
    :func:`~pypeit.core.wavecal.wvutils.shift_and_stretch`; i.e., reference
    pixel j appears at observed pixel shift + stretch*j.
    """
    rng = np.random.default_rng(seed=seed)
    x = np.arange(nspec, dtype=float)
    cen = rng.uniform(-0.1*nspec, 1.1*nspec, nlines)
    amp = rng.lognormal(mean=0., sigma=1., size=nlines)
    sig = fwhm / 2.355

    def _arc(pix):
        return np.sum(amp[:,None] * np.exp(-0.5*((pix[None,:] - cen[:,None])/sig)**2), axis=0)

    ref = _arc(x) + 10.
    obs = _arc((x - shift) / stretch) + 10. + rng.normal(loc=0., scale=0.02, size=nspec)
    return obs, ref


def _endpoint_errors(result, shift, stretch, nspec):
    """Pixel offsets of the recovered transform at both ends of the spectrum."""
    _, shift_out, stretch_out, _, _, _, _ = result
    return (
        abs(shift_out - shift),
        abs(shift_out + stretch_out * (nspec - 1) - shift - stretch * (nspec - 1))
    )


def test_xcorr_shift_stretch_linear():
    nspec, fwhm = 1024, 3.0
    shift, stretch = 37.3, 1.012
    obs, ref = _synthetic_arc_pair(shift, stretch, nspec=nspec, fwhm=fwhm)
    result = wvutils.xcorr_shift_stretch(obs, ref, fwhm=fwhm, stretch_func='linear')
    assert result[0] == 1, 'Shift/stretch search should succeed'
    assert result[3] == 0.0, 'Linear stretch should not have a quadratic term'
    assert result[4] > 0.9, 'Cross-correlation should be high'
    assert max(_endpoint_errors(result, shift, stretch, nspec)) < 0.2, \
            'Did not recover the input shift and stretch'


def test_xcorr_shift_stretch_quadratic():
    nspec, fwhm = 1024, 3.0
    shift, stretch = -52.6, 0.981
    obs, ref = _synthetic_arc_pair(shift, stretch, nspec=nspec, fwhm=fwhm, seed=7)
    result = wvutils.xcorr_shift_stretch(obs, ref, fwhm=fwhm, stretch_func='quadratic')
    assert result[0] == 1, 'Shift/stretch search should succeed'
    assert result[4] > 0.9, 'Cross-correlation should be high'
    assert max(_endpoint_errors(result, shift, stretch, nspec)) < 0.5, \
            'Did not recover the input shift and stretch'


def test_xcorr_shift_stretch_deterministic():
    obs, ref = _synthetic_arc_pair(20.4, 1.03, seed=3)
    result1 = wvutils.xcorr_shift_stretch(obs, ref, fwhm=3.0, stretch_func='linear')
    result2 = wvutils.xcorr_shift_stretch(obs, ref, fwhm=3.0, stretch_func='linear')
    assert result1 == result2, 'Shift/stretch search should be deterministic'


def test_xcorr_shift_stretch_no_lines():
    flat = np.ones(512, dtype=float)
    result = wvutils.xcorr_shift_stretch(flat, flat, fwhm=3.0, stretch_func='linear')
    assert result[0] == 0, 'Should fail if no lines are detected'
    assert all(r is None for r in result[1:]), 'Failure should return None values'
