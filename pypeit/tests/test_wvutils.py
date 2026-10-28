
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
