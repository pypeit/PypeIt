"""
Tests for Magellan/LDSS3 file handling, detector parameters and frame typing.
"""
import numpy as np
import pytest

from astropy.io import fits
from astropy.table import Table

from pypeit import PypeItError
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core.wavecal import waveio


SPEC = load_spectrograph('magellan_ldss3')

# Small synthetic amplifier: 16 data + 4 overscan columns, 8 spectral rows
NDATA, NOSCAN, NSPEC = 16, 4, 8


def write_amp(path, amp, speed='Fast', gain='Low'):
    """Write a synthetic LDSS3 amplifier file and return the data written."""
    rng = np.random.default_rng(amp)
    data = rng.integers(0, 1000, size=(NSPEC, NDATA + NOSCAN)).astype(np.uint16)
    hdr = fits.Header()
    hdr['OPAMP'] = amp
    hdr['DATASEC'] = f'[1:{NDATA},1:{NSPEC}]'
    hdr['BIASSEC'] = f'[{NDATA+1}:{NDATA+NOSCAN},1:{NSPEC}]'
    hdr['SPEED'] = speed
    hdr['GAIN'] = gain
    hdr['EGAIN'] = 1.6 if amp == 1 else 1.5
    hdr['ENOISE'] = 4.7
    hdr['BINNING'] = '1x1'
    hdr['EXPTIME'] = 10.
    hdr['INSTRUME'] = 'LDSS3-C'
    fits.PrimaryHDU(data=data, header=hdr).writeto(path)
    return data


def make_dir(root, amps_per_exposure):
    """Write synthetic files; ``amps_per_exposure`` maps exposure number to amplifiers."""
    root.mkdir()
    for num, amps in amps_per_exposure.items():
        for amp in amps:
            write_amp(root / f'ccd{num:04d}c{amp}.fits', amp)
    return root

def test_amp_files_gz(tmp_path):
    root = tmp_path / 'raw'
    root.mkdir()
    for amp in (1, 2):
        write_amp(root / f'ccd0001c{amp}.fits.gz', amp)
    # An uncompressed copy of amplifier 1 is not paired with the gzipped amplifier 2
    write_amp(root / 'ccd0001c1.fits', 1)
    assert [f.name for f in SPEC.amp_files(root / 'ccd0001c1.fits.gz')] \
            == ['ccd0001c1.fits.gz', 'ccd0001c2.fits.gz']
    assert SPEC.amp_files(root / 'ccd0001c1.fits')[1] is None

def test_parse_amp_file():
    assert SPEC.parse_amp_file('ccd0042c1.fits') == ('ccd0042', 1, '.fits')
    assert SPEC.parse_amp_file('/a/b/CCD0042c2.fits.gz') == ('CCD0042', 2, '.fits.gz')
    assert SPEC.parse_amp_file('ccd0042.fits') is None
    assert SPEC.parse_amp_file('merged_c1.fits') is None


def test_find_raw_files_multiple_dirs(tmp_path):
    make_dir(tmp_path / 'night1', {1: [1, 2]})
    make_dir(tmp_path / 'night2', {1: [1, 2]})
    with pytest.raises(PypeItError):
        SPEC.find_raw_files([tmp_path / 'night1', tmp_path / 'night2'])

def test_find_raw_files_pairs(tmp_path):
    root = make_dir(tmp_path / 'raw', {1: [1, 2], 2: [1, 2], 3: [1]})
    files = SPEC.find_raw_files(root)
    # Exposure 3 has no c2 file and is skipped; only the c1 files are returned
    assert [f.name for f in files] == ['ccd0001c1.fits', 'ccd0002c1.fits']


@pytest.mark.parametrize('amp', [1, 2])
def test_find_raw_files_single_amp(tmp_path, amp):
    root = make_dir(tmp_path / 'raw', {1: [amp], 2: [amp]})
    files = SPEC.find_raw_files(root)
    assert [f.name for f in files] == [f'ccd0001c{amp}.fits', f'ccd0002c{amp}.fits']

@pytest.mark.parametrize('arxiv', ['magellan_ldss3_vph_all_HeINeIArI.fits',
                                   'magellan_ldss3_vph_blue_HeINeIArI.fits',
                                   'magellan_ldss3_vph_red_HeINeIArI.fits',
                                   'magellan_ldss3_vph_red_sky.fits'])
def test_load_templates(arxiv):
    wave, flux = waveio.load_template(arxiv, 1)[:2]
    assert wave.size == flux.size > 0

def test_load_oh_line_list():
    assert len(waveio.load_line_list('OH_LDSS3_vac', use_ion=True)) == 300

def test_amp_files_opamp_mismatch(tmp_path):
    root = tmp_path / 'raw'
    root.mkdir()
    # File named as amplifier 1 but OPAMP says amplifier 2
    write_amp(root / 'ccd0001c1.fits', 2)
    assert SPEC.amp_files(root / 'ccd0001c1.fits') == [None, None]


def test_get_rawimage_two_amps(tmp_path):
    root = tmp_path / 'raw'
    root.mkdir()
    d1 = write_amp(root / 'ccd0001c1.fits', 1)
    d2 = write_amp(root / 'ccd0001c2.fits', 2)
    det, img, hdu, exptime, datasec, oscansec = SPEC.get_rawimage(root / 'ccd0001c1.fits', 1)
    hdu.close()

    width = NDATA + NOSCAN
    assert img.shape == (NSPEC, 2 * width)
    assert exptime == 10.
    # [ overscan_1 | data_1 | flip(data_2) | flip(overscan_2) ]
    assert np.array_equal(img[:, :NOSCAN], d1[:, NDATA:])
    assert np.array_equal(img[:, NOSCAN:width], d1[:, :NDATA])
    assert np.array_equal(img[:, width:width+NDATA], d2[:, :NDATA][:, ::-1])
    assert np.array_equal(img[:, width+NDATA:], d2[:, NDATA:][:, ::-1])
    assert np.array_equal(np.unique(datasec[0]), [0, 1, 2])
    assert np.all(datasec[:, NOSCAN:width] == 1) and np.all(datasec[:, width:width+NDATA] == 2)
    assert np.all(oscansec[:, :NOSCAN] == 1) and np.all(oscansec[:, width+NDATA:] == 2)

    assert det.numamplifiers == 2
    assert np.allclose(det.gain, [1.65, 1.44])
    assert np.allclose(det.ronoise, [6.4, 6.7])
    assert det.saturation == 65535.
    assert det.nonlinear == 0.95


def test_get_rawimage_amp2_only(tmp_path):
    root = tmp_path / 'raw'
    root.mkdir()
    d2 = write_amp(root / 'ccd0001c2.fits', 2)
    det, img, hdu, _, datasec, oscansec = SPEC.get_rawimage(root / 'ccd0001c2.fits', 1)
    hdu.close()

    width = NDATA + NOSCAN
    # Amplifier 2 stays in the high-spatial half, flipped, but is labelled 1
    assert img.shape == (NSPEC, 2 * width)
    assert np.all(img[:, :width] == 0)
    assert np.array_equal(img[:, width:width+NDATA], d2[:, :NDATA][:, ::-1])
    assert np.all(datasec[:, width:width+NDATA] == 1) and np.all(datasec[:, :width] == 0)
    assert np.all(oscansec[:, width+NDATA:] == 1)

    assert det.numamplifiers == 1
    assert np.allclose(det.gain, [1.44])
    assert np.allclose(det.ronoise, [6.7])


def test_get_detector_par_modes(tmp_path):
    root = tmp_path / 'raw'
    root.mkdir()
    write_amp(root / 'ccd0001c1.fits', 1, speed='Slow')
    write_amp(root / 'ccd0001c2.fits', 2, speed='Slow')
    write_amp(root / 'ccd0002c1.fits', 1, gain='High')
    write_amp(root / 'ccd0002c2.fits', 2, gain='High')

    hdrs = [fits.getheader(root / f'ccd0001c{a}.fits') for a in (1, 2)]
    det = SPEC.get_detector_par(1, amp_headers=hdrs)
    assert np.allclose(det.gain, [0.15, 0.13])
    assert np.allclose(det.ronoise, [4.4, 5.3])

    # The readout mode is also read from a single hdu
    with fits.open(root / 'ccd0001c1.fits') as hdu:
        det = SPEC.get_detector_par(1, hdu=hdu)
    assert np.allclose(det.gain, [0.15])

    # Uncharacterized modes fall back on the header cards
    hdrs = [fits.getheader(root / f'ccd0002c{a}.fits') for a in (1, 2)]
    det = SPEC.get_detector_par(1, amp_headers=hdrs)
    assert np.allclose(det.gain, [1.6, 1.5])
    assert np.allclose(det.ronoise, [4.7, 4.7])

    with pytest.raises(PypeItError):
        SPEC.get_detector_par(2)


def test_check_frame_type():
    tbl = Table({'exptime': [1., 70., 1200., 30., 0., 600., 10.],
                 'idname': ['Object', 'Flat', 'Object', 'Object', 'Bias', 'Object', 'Object'],
                 'target': ['HeNeAr', 'ES1_34_QlQh', 'science ATM3a2_1', 'LTT6248', 'bias',
                            'acq field', 'arc'],
                 'dispname': ['VPH-All', 'VPH-All', 'VPH-All', 'VPH-All', 'Open', 'VPH-All',
                              'Open']})
    par = SPEC.default_pypeit_par()
    cal = par['calibrations']
    assert np.array_equal(SPEC.check_frame_type('arc', tbl, cal['arcframe']['exprng']),
                          [1, 0, 0, 0, 0, 0, 0])
    assert np.array_equal(SPEC.check_frame_type('trace', tbl, cal['traceframe']['exprng']),
                          [0, 1, 0, 0, 0, 0, 0])
    assert np.array_equal(SPEC.check_frame_type('science', tbl, par['scienceframe']['exprng']),
                          [0, 0, 1, 0, 0, 0, 0])
    assert np.array_equal(SPEC.check_frame_type('standard', tbl, cal['standardframe']['exprng']),
                          [0, 0, 0, 1, 0, 0, 0])
    assert np.array_equal(SPEC.check_frame_type('bias', tbl, cal['biasframe']['exprng']),
                          [0, 0, 0, 0, 1, 0, 0])
