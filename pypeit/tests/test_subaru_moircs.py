"""
Unit tests for the Subaru/MOIRCS spectrograph class, using synthetic
headers and files only.
"""
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from pypeit import PypeItError
from pypeit.spectrographs.util import load_spectrograph

SPEC = load_spectrograph('subaru_moircs')


def _write_frame(path, det_id, exp_id):
    hdr = fits.Header()
    hdr['DET-ID'] = det_id
    hdr['EXP-ID'] = exp_id
    fits.writeto(path, np.zeros((4, 4), dtype=np.float32), hdr)
    return path


def test_companion_match(tmp_path):
    f1 = _write_frame(tmp_path / 'MCSP00237323.fits', 1, 'MCSE00237323')
    f2 = _write_frame(tmp_path / 'MCSP00237324.fits', 2, 'MCSE00237323')
    assert SPEC.companion_file(f1) == f2


def test_companion_missing(tmp_path):
    f1 = _write_frame(tmp_path / 'MCSP00237323.fits', 1, 'MCSE00237323')
    with pytest.raises(PypeItError, match='Missing the chip-2 file'):
        SPEC.companion_file(f1)


def test_companion_expid_mismatch(tmp_path):
    f1 = _write_frame(tmp_path / 'MCSP00237323.fits', 1, 'MCSE00237323')
    _write_frame(tmp_path / 'MCSP00237324.fits', 2, 'MCSE00237325')
    with pytest.raises(PypeItError, match='not the chip-2 companion'):
        SPEC.companion_file(f1)


def test_companion_not_chip1(tmp_path):
    f2 = _write_frame(tmp_path / 'MCSP00237324.fits', 2, 'MCSE00237323')
    with pytest.raises(PypeItError, match='not a MOIRCS chip-1 file'):
        SPEC.companion_file(f2)


def test_frame_typing():
    datatyp = ['DOMEFLAT', 'DOMEFLAT', 'DOMEFLAT', 'OBJECT']
    obj = ['DOMEFLAT', 'DOMEFLAT_OFF', 'MASKIMAGE', 'COSMOS2']
    hdrs = [fits.Header({'DATA-TYP': d, 'OBJECT': o})
            for d, o in zip(datatyp, obj)]
    fitstbl = Table()
    fitstbl['idname'] = [SPEC.compound_meta([h], 'idname') for h in hdrs]
    fitstbl['lampstat01'] = [SPEC.compound_meta([h], 'lampstat01')
                             for h in hdrs]
    fitstbl['exptime'] = [5., 5., 5., 180.]

    assert list(fitstbl['idname']) \
        == ['DOMEFLAT', 'DOMEFLAT_OFF', 'MASKIMAGE', 'OBJECT']
    assert list(fitstbl['lampstat01']) == ['on', 'off', 'off', 'off']
    assert np.array_equal(SPEC.check_frame_type('pixelflat', fitstbl),
                          [True, False, False, False])
    assert np.array_equal(SPEC.check_frame_type('trace', fitstbl),
                          [True, False, False, False])
    assert np.array_equal(SPEC.check_frame_type('lampoffflats', fitstbl),
                          [False, True, False, False])
    assert np.array_equal(SPEC.check_frame_type('science', fitstbl),
                          [False, False, False, True])
    assert np.array_equal(SPEC.check_frame_type('arc', fitstbl),
                          [False, False, False, True])
    # The mask image is not typed at all
    for ftype in ['pixelflat', 'lampoffflats', 'science', 'arc', 'bias',
                  'dark', 'standard']:
        assert not SPEC.check_frame_type(ftype, fitstbl)[2]


def test_dither_parsing():
    hdr_a = fits.Header({'K_DITPAT': 'LINE2', 'K_DITCNT': 1.,
                         'K_DITWID': 3.})
    hdr_b = fits.Header({'K_DITPAT': 'LINE2', 'K_DITCNT': 2.,
                         'K_DITWID': 3.})
    hdr_n = fits.Header({'K_DITPAT': 'NONE', 'K_DITCNT': 0.,
                         'K_DITWID': 0.})
    assert SPEC.compound_meta([hdr_a], 'dithpat') == 'LINE2'
    assert SPEC.compound_meta([hdr_a], 'dithpos') == 'A'
    assert SPEC.compound_meta([hdr_b], 'dithpos') == 'B'
    assert SPEC.compound_meta([hdr_a], 'dithoff') == 1.5
    assert SPEC.compound_meta([hdr_b], 'dithoff') == -1.5
    assert SPEC.compound_meta([hdr_n], 'dithpat') == 'none'
    assert SPEC.compound_meta([hdr_n], 'dithpos') == 'none'
    assert SPEC.compound_meta([hdr_n], 'dithoff') == 0.


def test_comb_group_abba():
    tbl = Table()
    tbl['frametype'] = ['arc,science,tilt'] * 4 + ['pixelflat']
    tbl['setup'] = ['A'] * 5
    tbl['dithpat'] = ['LINE2'] * 4 + ['none']
    tbl['dithpos'] = ['A', 'B', 'B', 'A', 'none']
    tbl['mjd'] = [0., 1., 2., 3., -1.]
    tbl['comb_id'] = [1, 2, 3, 4, -1]
    tbl['bkg_id'] = [-1] * 5
    tbl = SPEC.get_comb_group(tbl)
    # Each frame is paired with the closest frame at the other position
    assert list(tbl['bkg_id']) == [2, 1, 4, 3, -1]
    assert list(tbl['comb_id']) == [1, 2, 3, 4, -1]


def test_binning():
    hdr = fits.Header({'BIN-FCT1': 2, 'BIN-FCT2': 1})
    # BIN-FCT1 is along the dispersion axis: binning is 'spec,spat'
    assert SPEC.compound_meta([hdr], 'binning') == '2,1'
    assert SPEC.compound_meta([fits.Header()], 'binning') == '1,1'


def test_valid_detector():
    assert SPEC.valid_configuration_values() == {'detector': ['1']}
    assert SPEC.ndet == 2
