"""
Module to run tests on calibration frames
"""
from pathlib import Path

from IPython import embed

import pytest

from pypeit.pypmsgs import PypeItError
from pypeit.calibframe import CalibFrame
from pypeit import io
from pypeit.tests.tstutils import data_path


class NoTypeCalibFrame(CalibFrame):
    version = '1.0.0'


class MissingPYPSPECCalibFrame(CalibFrame):
    version = '1.0.0'
    calib_type = 'Junk'
    datamodel = {'Junk': dict(otype=str)}


class MinimalCalibFrame(CalibFrame):
    # These are the minimal things that need to be defined to actually
    # instantiate a derive class of CalibFrame
    version = '1.0.0'
    calib_type = 'Minimal'


def test_implementation_faults():
    # CalibFrame cannot be instantiated by itself because a version of the
    # datamodel does not exist.
    with pytest.raises(ValueError):
        calib = CalibFrame()
    # Any calibframe needs to define its type
    with pytest.raises(PypeItError):
        calib = NoTypeCalibFrame()
    # The elements of the base-class datamodel *must* exist in all its derived
    # classes.
    with pytest.raises(PypeItError):
        calib = MissingPYPSPECCalibFrame()


def test_init():
    calib = MinimalCalibFrame()
    odir = Path(data_path('')).resolve()
    calib.set_paths(odir, 'A', '1', 'DET01')
    ofile = Path(calib.get_path()).name
    assert ofile == 'Minimal_A_1_DET01.fits', 'Wrong file name'
    opath = Path(calib.get_path()).resolve() # Now with the full path
    assert opath.parent == odir, 'Wrong parent directory'
    assert opath.name == ofile, 'Wrong file name'

    calib.set_paths(odir, 'A', ['1','2'], 'DET01')
    ofile = Path(calib.get_path()).name
    assert ofile == 'Minimal_A_1-2_DET01.fits', 'Wrong file name'


def test_io():
    calib = MinimalCalibFrame()
    odir = Path(data_path('')).resolve()
    calib.set_paths(odir, 'A', '1', 'DET01')
    calib.PYP_SPEC = 'this is a test'
    opath = Path(calib.get_path()).resolve()
    calib.to_file(overwrite=True)

    with io.fits_open(str(opath)) as hdu:
        assert hdu[1].header['CALIBTYP'] == calib.calib_type, 'Calibration type incorrect'
        assert len(hdu) == 2, 'Should be two HDUs'
        assert hdu[1].name == '', 'HDU 1 name should be empty'
        assert all([h.data is None for h in hdu]), 'File should have no data'

    _calib = MinimalCalibFrame.from_file(str(opath))
    assert _calib.calib_key == calib.calib_key, 'Error parsing calib key from file name'
    assert _calib.calib_dir == calib.calib_dir, 'Error parsing calib directory from file name'
    assert _calib.PYP_SPEC == calib.PYP_SPEC, 'Datamodel read error'

    opath.unlink()

    calib.set_paths(odir, 'A', [1,2], 'DET01')
    calib.to_file(overwrite=True)
    opath = Path(calib.get_path()).resolve()
    with io.fits_open(str(opath)) as hdu:
        assert hdu[1].header['CALIBID'] == '1,2', 'Calibration ID incorrect'
    _calib = MinimalCalibFrame.from_file(str(opath))

    assert _calib.calib_id == calib.calib_id, 'Calibration ID incorrect'
    assert _calib.calib_id == ['1', '2'], 'Calibration ID incorrect'

    opath.unlink()


def test_construct_calib_key():
    key = CalibFrame.construct_calib_key('A', '1', 'DET01')
    assert key == 'A_1_DET01', 'Key changed'
    key = CalibFrame.construct_calib_key('A', ['1','2'], 'DET01')
    assert key == 'A_1-2_DET01', 'Key changed'
    key = CalibFrame.construct_calib_key('A', 'all', 'DET01')
    assert key == 'A_all_DET01', 'Key changed'


def test_ingest_calib_id():
    assert CalibFrame.ingest_calib_id('all') == ['all'], 'Bad ingest'
    assert CalibFrame.ingest_calib_id(['all', 1]) == ['all'], 'Bad ingest'
    assert CalibFrame.ingest_calib_id('1,2') == ['1', '2'], 'Bad ingest'
    assert CalibFrame.ingest_calib_id(['1,2']) == ['1', '2'], 'Bad ingest'
    assert CalibFrame.ingest_calib_id([1, 2]) == ['1', '2'], 'Bad ingest'
    assert CalibFrame.ingest_calib_id([2, 1]) == ['1', '2'], 'Bad ingest'
    assert CalibFrame.ingest_calib_id([2, 1, 2]) == ['1', '2'], 'Bad ingest'
    assert CalibFrame.ingest_calib_id(['1', 2]) == ['1', '2'], 'Bad ingest'
    assert CalibFrame.ingest_calib_id(['1', '2']) == CalibFrame.ingest_calib_id([2, 1]), \
            'Bad ingest'


def test_parse_key_dir():
    calib = MinimalCalibFrame()
    odir = Path(data_path('')).resolve()
    calib.set_paths(odir, 'A', '1', 'DET01')
    calib.PYP_SPEC = 'this is a test'
    opath = Path(calib.get_path()).resolve()
    calib.to_file(overwrite=True)

    key, _odir = CalibFrame.parse_key_dir(str(opath), from_filename=True)
    assert key == calib.calib_key, 'Key parsed incorrectly'
    assert _odir == calib.calib_dir, 'Key parsed incorrectly'
    key, _dir = CalibFrame.parse_key_dir(str(opath))
    assert key == calib.calib_key, 'Key parsed incorrectly'
    assert _odir == calib.calib_dir, 'Key parsed incorrectly'
    with io.fits_open(str(opath)) as hdu:
        key, _dir = CalibFrame.parse_key_dir(hdu[1].header)
        assert key == calib.calib_key, 'Key parsed incorrectly'
        assert _odir == calib.calib_dir, 'Key parsed incorrectly'

    opath.unlink()


def test_hdr():
    calib = MinimalCalibFrame()
    odir = Path(data_path('')).resolve()
    calib.set_paths(odir, 'A', '1', 'DET01')

    hdr = calib._base_header()
    assert 'CALIBTYP' in hdr, 'Missing keyword'
    assert hdr['CALIBTYP'] == calib.calib_type
    assert 'CALIBDIR' in hdr, 'Missing keyword'
    assert hdr['CALIBDIR'] == calib.calib_dir
    assert 'CALIBKEY' in hdr, 'Missing keyword'
    assert hdr['CALIBKEY'] == calib.calib_key
    assert 'CALIBID' in hdr, 'Missing keyword'
    assert hdr['CALIBID'] == ','.join(calib.calib_id)


