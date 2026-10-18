"""
Module to run tests on WaveCalib and WaveFit
"""
from pathlib import Path
from types import SimpleNamespace
from IPython import embed

import pytest

import numpy as np
from astropy.table import Table

from pypeit.core.wavecal import wv_fitting
from pypeit.core.wavecal import autoid
from pypeit.core.wavecal import wvutils
from pypeit.core import arc
from pypeit.core import fitting
from pypeit import wavecalib
from pypeit import slittrace
from pypeit.tests.tstutils import data_output_path
from pypeit.core.wavecal import waveio
from pypeit.par import pypeitpar


@pytest.mark.remote_data
def test_load_template():
    # Test full read
    ret = waveio.load_template('gemini_gmos_r831_ham.fits', 1)
    assert len(ret) == 7, 'Number of returned objects changed'

    # Only catch the spectrum, and test its contents
    wave, flux = waveio.load_template('gemini_gmos_r831_ham.fits', 1)[:2]
    assert wave[0] > 6200, 'Bad wavelength read'
    assert np.amax(flux) > 1e5, 'Bad flux read'


def test_wavefit_hduprefix():
    spat_id = 175
    prefix = wv_fitting.WaveFit.hduext_prefix_from_spatid(spat_id)
    _spat_id = wv_fitting.WaveFit.parse_spatid_from_hduext(prefix)
    assert spat_id == _spat_id, 'Bad parse'


def test_wavefit():
    "Fuss with the WaveFit DataContainer"
    out_file = Path(data_output_path('test_wavefit.fits')).absolute()
    if out_file.exists():
        out_file.unlink()
    pypeitFit = fitting.PypeItFit(fitc=np.arange(5).astype(float))
    #
    ions = np.asarray(['CdI', 'HgI'])
    ion_bits = np.zeros(len(ions), dtype=wv_fitting.WaveFit.bitmask.minimum_dtype())
    for kk,ion in enumerate(ions):
        ion_bits[kk] = wv_fitting.WaveFit.bitmask.turn_on(ion_bits[kk], ion)
    # Do it
    waveFit = wv_fitting.WaveFit(99, pypeitfit=pypeitFit, pixel_fit=np.arange(10).astype(float),
                                 wave_fit=np.linspace(1.,10.,10), sigrej=3.,
                                 ion_bits=ion_bits)

    # Write
    waveFit.to_file(out_file)
    # Read
    waveFit2 = wv_fitting.WaveFit.from_file(out_file)
    assert np.array_equal(waveFit.pypeitfit.fitc, waveFit2.pypeitfit.fitc)
    # Write again
    waveFit2.to_file(out_file, overwrite=True)
    # And one more read
    waveFit2b = wv_fitting.WaveFit.from_file(out_file)
    assert np.array_equal(waveFit.pypeitfit.fitc, waveFit2b.pypeitfit.fitc)
    # Finish
    out_file.unlink()

    # No fit
    waveFit3 = wv_fitting.WaveFit(99, pypeitfit=None, pixel_fit=np.arange(10).astype(float),
                                 wave_fit=np.linspace(1.,10.,10), sigrej=3.,
                                 ion_bits=ion_bits)
    waveFit3.to_file(out_file)
    waveFit4 = wv_fitting.WaveFit.from_file(out_file)
    assert waveFit4.pypeitfit is None

    # Clean up
    out_file.unlink()


def test_wavecalib():
    "Fuss with the WaveCalib DataContainer"
    # Pieces
    pypeitFit = fitting.PypeItFit(fitc=np.arange(5).astype(float), xval=np.linspace(1,100., 100))
    # 2D fit
    pypeitFit2 = fitting.PypeItFit(fitc=np.linspace((1,2),(10,20),10),
                                   xval=np.linspace(1,100., 100),
                                   x2=np.linspace(1, 100., 100))
    waveFit = wv_fitting.WaveFit(232, pypeitfit=pypeitFit, pixel_fit=np.arange(10).astype(float),
                                 wave_fit=np.linspace(1.,10.,10))

    waveCalib = wavecalib.WaveCalib(wv_fits=np.asarray([waveFit]),
                                    nslits=1, spat_ids=np.asarray([232]),
                                    wv_fit2d=np.array([pypeitFit2]),
                                    fwhm_map=np.array([pypeitFit2]))
    waveCalib.set_paths(data_output_path(''), 'A', '1', 'DET01')

    ofile = Path(waveCalib.get_path()).absolute()

    # Write
    waveCalib.to_file(overwrite=True)
    assert ofile.exists(), 'File not written'

    # Read
    waveCalib2 = wavecalib.WaveCalib.from_file(ofile)

    # Test
    assert np.array_equal(waveCalib.spat_ids, waveCalib2.spat_ids), 'Bad spat_ids'
    assert np.array_equal(waveCalib.wv_fits[0].pypeitfit.fitc,
                          waveCalib2.wv_fits[0].pypeitfit.fitc), 'Bad fitc'
    assert np.array_equal(waveCalib.wv_fit2d[0].xval, waveCalib2.wv_fit2d[0].xval)
    assert np.array_equal(waveCalib.fwhm_map[0].xval, waveCalib2.fwhm_map[0].xval)

    # Write again!
    waveCalib2.to_file(overwrite=True)

    # Finish
    ofile.unlink()

    # With None (failed wave) and no FWHM mapping
    spat_ids = np.asarray([232, 949])
    waveCalib3 = wavecalib.WaveCalib(wv_fits=np.asarray([waveFit, wv_fitting.WaveFit(949)]),
                                     nslits=2, spat_ids=spat_ids,
                                     wv_fit2d=np.array([pypeitFit2, pypeitFit2]))
    waveCalib3.set_paths(data_output_path(''), 'A', '1', 'DET01')
    waveCalib3.to_file(overwrite=True)
    waveCalib4 = wavecalib.WaveCalib.from_file(ofile)

    # Check masking
    slits = slittrace.SlitTraceSet(left_init=np.full((1000,2), 2, dtype=float),
                         right_init=np.full((1000,2), 8, dtype=float),
                         pypeline='MultiSlit', spat_id=spat_ids,
                         nspat=2, PYP_SPEC='dummy')
    slits.mask_wvcalib(waveCalib3)
    assert slits.bitmask.flagged(slits.mask[1], flag='BADWVCALIB')

    # Finish
    ofile.unlink()


def test_wvcalib_no2d():
    """ Test WaveCalib without a 2D fit (not echelle) """
    # Pieces
    pypeitFit = fitting.PypeItFit(fitc=np.arange(5).astype(float), xval=np.linspace(1,100., 100))
    waveFit = wv_fitting.WaveFit(232, pypeitfit=pypeitFit, pixel_fit=np.arange(10).astype(float),
                                 wave_fit=np.linspace(1.,10.,10))

    # With None (failed wave)
    spat_ids = np.asarray([232, 949])
    waveCalib = wavecalib.WaveCalib(wv_fits=np.asarray([waveFit, wv_fitting.WaveFit(949)]),
                                    nslits=2, spat_ids=spat_ids)
    waveCalib.set_paths(data_output_path(''), 'A', '1', 'DET01')
    ofile = Path(waveCalib.get_path()).absolute()

    waveCalib.to_file(overwrite=True)
    assert ofile.exists(), 'File not written'

    waveCalib2 = wavecalib.WaveCalib.from_file(ofile)
    assert len(waveCalib2.wv_fits) == 2, 'Should be two wavelength fits'
    assert waveCalib2.wv_fits[0].spat_id == 232, 'Wrong spat_id'
    assert np.array_equal(waveCalib2.wv_fits[0].wave_fit, waveCalib.wv_fits[0].wave_fit), \
            'Bad wave_fit read'

    # Finish
    ofile.unlink()



def test_order_for_dof():
    # Enough lines for the requested order
    assert wv_fitting.order_for_dof(4, 10) == 4
    # Order reduced to leave 2 degrees of freedom
    assert wv_fitting.order_for_dof(4, 6) == 3
    assert wv_fitting.order_for_dof(4, 5) == 2
    assert wv_fitting.order_for_dof(4, 4) == 1
    # Not enough lines for a linear fit with 2 degrees of freedom
    assert wv_fitting.order_for_dof(4, 3) is None
    # Without the degrees of freedom requirement
    assert wv_fitting.order_for_dof(4, 5, min_dof=0) == 4
    assert wv_fitting.order_for_dof(4, 2, min_dof=0) == 1


def _synthetic_line_list():
    """A simple line list for the synthetic arc spectra used below."""
    wave = np.linspace(5050., 5950., 19)
    return Table(dict(ion=np.full(wave.size, 'ArI'), wave=wave, NIST=np.ones(wave.size, dtype=int),
                      Source=np.full(wave.size, 'NIST')))


def test_iterative_fitting_dof():
    # Cubic wavelength solution, and the pixel positions of the lines
    nspec = 1000
    llist = _synthetic_line_list()
    pix = np.arange(nspec, dtype=float)
    wave = 5000. + 1.0*pix + 1e-4*pix**2 - 5e-8*pix**3
    tcent = np.interp(llist['wave'], wave, pix)
    spec = np.zeros(nspec)
    disp = 1.0

    # Plenty of lines: the requested order is used
    fit = wv_fitting.iterative_fitting(spec, tcent, np.arange(tcent.size), llist['wave'].data, llist,
                                       disp, n_first=2, n_final=4, input_only=True)
    assert fit.pypeitfit.order[0] == 4, 'Order should not be reduced'

    # Only 5 lines: the order is reduced to keep 2 degrees of freedom, so the
    # fit does not pass exactly through all lines
    ifit = np.array([2, 6, 9, 12, 16])
    fit = wv_fitting.iterative_fitting(spec, tcent, ifit, llist['wave'].data[ifit], llist, disp,
                                       n_first=2, n_final=4, input_only=True)
    assert fit.pypeitfit.order[0] == 2, 'Order should be reduced to nlines-3'
    assert fit.rms > 0, 'Fit should not pass exactly through the lines'

    # Only 3 lines: the fit fails
    ifit = np.array([2, 9, 16])
    fit = wv_fitting.iterative_fitting(spec, tcent, ifit, llist['wave'].data[ifit], llist, disp,
                                       n_first=2, n_final=4, input_only=True)
    assert fit is None, 'Fit should fail with too few lines'

    # Same lines, without the degrees of freedom requirement: the order is
    # reduced to nlines-1, as used when re-fitting echelle orders
    fit = wv_fitting.iterative_fitting(spec, tcent, ifit, llist['wave'].data[ifit], llist, disp,
                                       n_first=2, n_final=4, input_only=True, min_dof=0)
    assert fit.pypeitfit.order[0] == 2, 'Order should be reduced to nlines-1'


def test_fit_slit_min_dof():
    # Three identified lines of a cubic wavelength solution
    nspec = 1000
    llist = _synthetic_line_list()
    pix = np.arange(nspec, dtype=float)
    wave = 5000. + 1.0*pix + 1e-4*pix**2 - 5e-8*pix**3
    tcent = np.interp(llist['wave'], wave, pix)
    ids = np.zeros(tcent.size)
    use = np.array([2, 9, 16])
    ids[use] = llist['wave'][use]
    mask = np.zeros(tcent.size, dtype=bool)
    mask[use] = True
    patt_dict = dict(mask=mask, IDs=list(ids), bdisp=1.0)
    spec = np.zeros(nspec)
    # The default requires 2 degrees of freedom, so the fit fails.  NOTE: The
    # small match_toler prevents more lines from being added to the fit.
    fit = wv_fitting.fit_slit(spec, patt_dict, tcent, llist, n_final=4, match_toler=0.01)
    assert fit is None, 'Fit should fail with too few lines'
    # Without the requirement, the fit succeeds with a reduced order
    fit = wv_fitting.fit_slit(spec, patt_dict, tcent, llist, n_final=4, match_toler=0.01, min_dof=0)
    assert fit is not None and fit.pypeitfit.order[0] == 2


def _synthetic_arcs(rng, nspec=1000, shifts=(0., 25., -30., 40.), noise=(0.5, 2., 2., 2.),
                    min_amp=0.2):
    """
    Synthetic arc spectra of several slits that are shifted with respect to
    each other, with different noise levels.
    """
    line_pix = np.linspace(60., 1940., 40)
    amp = 100.*rng.uniform(min_amp, 1.0, line_pix.size)
    sig = 1.5
    pix = np.arange(nspec, dtype=float)
    spec = np.zeros((nspec, len(shifts)))
    for i, (shift, _noise) in enumerate(zip(shifts, noise)):
        # Pixel in slit i where each line falls: shift and a small stretch
        lp = (line_pix - 500.)*1.002 + shift
        for p, a in zip(lp, amp):
            spec[:, i] += a*np.exp(-0.5*((pix - p)/sig)**2)
        spec[:, i] += 10. + rng.normal(0., _noise, nspec)
    return spec, line_pix, amp


def test_select_reference_slit():
    rng = np.random.default_rng(42)
    # Same wavelength coverage, but slit 0 has the lowest noise, so the most
    # detected lines
    spec, _, _ = _synthetic_arcs(rng, shifts=(0., 0., 0., 0.), noise=(0.5, 5., 5., 5.),
                                 min_amp=0.02)
    assert autoid.select_reference_slit(spec, np.arange(spec.shape[1]), fwhm=3.5) == 0
    # Exclude slit 0
    assert autoid.select_reference_slit(spec, np.arange(1, spec.shape[1]), fwhm=3.5) != 0


def test_stack_arc_spectra():
    rng = np.random.default_rng(42)
    spec, line_pix, amp = _synthetic_arcs(rng, noise=(2., 2., 2., 2.))
    nspec, nslit = spec.shape
    # Add a slit with no arc lines; it should not be included in the stack
    spec = np.column_stack([spec, 10. + rng.normal(0., 2., nspec)])
    ref_idx = 1
    stack, offset, used, ref_pix = autoid.stack_arc_spectra(spec, np.arange(spec.shape[1]), ref_idx,
                                                            fwhm=3.5, cc_thresh=0.7)
    assert np.array_equal(used, [True, True, True, True, False]), 'Wrong slits included'
    # The pixel mapping onto the reference slit recovers the input shifts
    # (relative to the reference slit; the stretch is the same for all slits)
    assert np.all(np.isnan(ref_pix[:, 4])), 'Mapping should be undefined for the excluded slit'
    for slit, shift in zip(range(4), (0., 25., -30., 40.)):
        assert np.allclose(ref_pix[:, slit] - np.arange(nspec), 25. - shift, atol=0.2), \
            f'Bad mapping for slit {slit}'
    # The stack covers all slits, so it extends beyond the reference slit
    assert offset < 0 and stack.size + offset > nspec
    # The lines are at the same pixels as in the reference slit
    _, _, cen, _, _, good, _, _ = arc.detect_lines(stack, sigdetect=10., fwhm=3.5)
    ref_lp = (line_pix - 500.)*1.002 + 25.
    ref_lp = ref_lp[(ref_lp > 20) & (ref_lp < nspec - 20)]
    stack_lp = cen[good] + offset
    assert np.all([np.min(np.abs(stack_lp - p)) < 0.2 for p in ref_lp]), 'Lines not aligned'
    # The noise in the stack is lower than in the reference slit.  Compare
    # with the noiseless spectrum of the reference slit, in the region covered
    # by all slits.
    pix = np.arange(nspec, dtype=float)
    model = np.sum([a*np.exp(-0.5*((pix - p)/1.5)**2)
                    for p, a in zip((line_pix - 500.)*1.002 + 25., amp)], axis=0)
    _, _, _, _, ref_cs = wvutils.arc_lines_from_spec(spec[:, ref_idx], fwhm=3.5)
    _, _, _, _, stack_cs = wvutils.arc_lines_from_spec(stack, fwhm=3.5)
    common = slice(100, 900)
    stack_resid = stack_cs[-offset:-offset+nspec][common] - model[common]
    ref_resid = ref_cs[common] - model[common]
    assert np.std(stack_resid) < 0.7*np.std(ref_resid), 'Stack noise is too high'


def _synthetic_arc_from_lines(wave, lines, amps, noise, rng, sig_pix=1.0):
    """Synthetic arc spectrum with Gaussian lines at the given wavelengths."""
    disp = np.gradient(wave)
    pix = np.arange(wave.size, dtype=float)
    lpix = np.interp(lines, wave, pix, left=np.nan, right=np.nan)
    spec = np.full(wave.size, 10.)
    for p, a in zip(lpix, amps):
        if np.isfinite(p):
            spec += a*np.exp(-0.5*((pix - p)/sig_pix)**2)
    return spec + rng.normal(0., noise, wave.size)


def test_full_template_stack():
    rng = np.random.default_rng(1)
    llist, _, _ = waveio.load_line_lists(['ArI'], include_unknown=False)
    lines = np.asarray(llist['wave'].data)
    amps = rng.uniform(0.01, 1.0, lines.size)**2

    # High S/N template covering a wide wavelength range
    tpix = np.arange(3000, dtype=float)
    twave = 6800. + 0.5*tpix + 2e-5*tpix**2
    tspec = _synthetic_arc_from_lines(twave, lines, 1e4*amps, 1., rng)
    tfile = Path(data_output_path('test_full_template_stack.fits')).absolute()
    tbl = Table(dict(wave=twave, flux=tspec))
    tbl.meta['BINSPEC'] = 1
    tbl.write(tfile, overwrite=True)

    # Low S/N arc spectra of 4 slits with the same dispersion but different
    # wavelength coverage
    nspec = 2000
    pix = np.arange(nspec, dtype=float)
    shifts = [0., 10., -15., 20.]
    wave = np.column_stack([7000. + s + 0.5*pix + 2e-5*pix**2 for s in shifts])
    spec = np.column_stack([_synthetic_arc_from_lines(wave[:, i], lines, 300*amps, 3., rng)
                            for i in range(len(shifts))])

    par = pypeitpar.WavelengthSolutionPar(method='full_template', lamps=['ArI'], reid_arxiv=str(tfile))
    fwhms = np.full(len(shifts), 2.35)
    wvcalib, stack_fit, used = autoid.full_template_stack(spec, ['ArI'], par, np.arange(len(shifts)), 1, 1,
                                                          0, fwhms)
    tfile.unlink()
    assert stack_fit is not None, 'Stack should be calibrated'
    assert np.all(used), 'All slits should be included in the stack'
    for i in range(len(shifts)):
        dev = np.abs(wvcalib[str(i)].wave_soln - wave[:, i]) / np.gradient(wave[:, i])
        assert np.max(dev) < 0.3, f'Bad wavelength solution for slit {i}'


def test_mask_bad_rms():
    # Three slits: a good solution, a solution with a poor RMS, and a slit
    # that is already masked
    pypeitFit = fitting.PypeItFit(fitc=np.arange(5).astype(float), xval=np.linspace(1, 100., 100))
    spat_ids = np.array([100, 200, 300])
    wv_fits = np.asarray([wv_fitting.WaveFit(s, pypeitfit=pypeitFit, pixel_fit=np.arange(10).astype(float),
                                             wave_fit=np.linspace(1., 10., 10), rms=rms)
                          for s, rms in zip(spat_ids, [0.1, 2.0, 5.0])])
    # Build a minimal BuildWaveCalib object
    bwc = wavecalib.BuildWaveCalib.__new__(wavecalib.BuildWaveCalib)
    bwc.wv_calib = wavecalib.WaveCalib(wv_fits=wv_fits, nslits=3, spat_ids=spat_ids)
    bwc.par = pypeitpar.WavelengthSolutionPar(rms_thresh_frac_fwhm=0.15)
    bwc.measured_fwhms = np.array([4.0, 4.0, 4.0])  # RMS threshold = 0.6 pixels
    bwc.wvc_bpm = np.array([False, False, True])
    bwc.slits = SimpleNamespace(slitord_id=spat_ids)

    bad_rms = bwc.mask_bad_rms()
    assert np.array_equal(bad_rms, [False, True, False]), 'Wrong slits rejected'
    assert np.array_equal(bwc.wvc_bpm, [False, True, True]), 'Wrong slits masked'
    # The good solution is kept, the poor one is replaced by an empty WaveFit,
    # and the masked slit is not changed
    assert bwc.wv_calib.wv_fits[0].pypeitfit is not None
    assert bwc.wv_calib.wv_fits[1].pypeitfit is None and bwc.wv_calib.wv_fits[1].spat_id == 200
    assert bwc.wv_calib.wv_fits[2].pypeitfit is not None
