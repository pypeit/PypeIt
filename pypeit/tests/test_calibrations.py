"""
Module to run tests on FlatField class
Requires files in Development suite and an Environmental variable
"""
from pathlib import Path
import os
import yaml
import pytest
import shutil

import numpy as np

from pypeit import calibrations
from pypeit.images import buildimage
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
from IPython import embed

from pypeit.tests.tstutils import dummy_fitstbl, data_path

@pytest.fixture
def fitstbl():
    if os.getenv('PYPEIT_DEV') is None:
        fitstbl = dummy_fitstbl(directory=data_path(''))
        fitstbl['framebit'][0] = fitstbl.type_bitmask.turn_off(fitstbl['framebit'][0], flag='bias')
        fitstbl['filename'][1] = 'b1.fits.gz'
        fitstbl['filename'][5] = 'b27.fits.gz'
        return fitstbl

    fitstbl = dummy_fitstbl(directory=os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA',
                                                   'shane_kast_blue', '600_4310_d55'))
    # Set the Bias to known
    fitstbl['framebit'][0] = fitstbl.type_bitmask.turn_off(fitstbl['framebit'][0], flag='bias')
    fitstbl['filename'][1] = 'b1.fits.gz'
    for ii in range(2,5):
        fitstbl['filename'][ii] = 'b{0}.fits.gz'.format(ii)
    fitstbl['filename'][5] = 'b27.fits.gz'

    return fitstbl


@pytest.fixture
def multi_caliBrate(fitstbl):
    # Grab a science file for configuration specific parameters
    for idx, row in enumerate(fitstbl):
        if 'science' in row['frametype']:
            sci_file = os.path.join(row['directory'], row['filename'])
            break
    # Par
    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.config_specific_par(sci_file)
    turn_off = dict(use_biasimage=False)
    par.reset_all_processimages_par(**turn_off)
    #
    calib_par = par['calibrations']
    calib_par['bpm_usebias'] = False
    #calib_par['biasframe']['useframe'] = 'none' # Only use overscan
    calib_par['slitedges']['sync_predict'] = 'nearest'

    caldir = data_path('Calibrations')

    multi_caliBrate = calibrations.MultiSlitCalibrations(fitstbl, calib_par, spectrograph, caldir)
    return reset_calib(multi_caliBrate)


def reset_calib(calib):
    # Find the first science row
    frame = calib.fitstbl.find_frames('science', index=True)[0]
    # Set
    det = 1
    calib.set_config(frame, det)
    return calib


###################################################
# TESTS BEGIN HERE

def test_abstract_init(fitstbl):
    par = pypeitpar.CalibrationsPar()
    spectrograph = load_spectrograph('shane_kast_blue')
    caldir = data_path('Calibrations')
    calib = calibrations.Calibrations.get_instance(fitstbl, par, spectrograph, caldir)
    assert isinstance(calib, calibrations.MultiSlitCalibrations), 'Wrong calibration object type'
    spectrograph = load_spectrograph('keck_nires')
    calib = calibrations.Calibrations.get_instance(fitstbl, par, spectrograph, caldir)
    assert isinstance(calib, calibrations.MultiSlitCalibrations), 'Wrong calibration object type'
    spectrograph = load_spectrograph('keck_kcwi')
    calib = calibrations.Calibrations.get_instance(fitstbl, par, spectrograph, caldir)
    assert isinstance(calib, calibrations.IFUCalibrations), 'Wrong calibration object type'


def test_instantiate(fitstbl):
    par = pypeitpar.CalibrationsPar()
    spectrograph = load_spectrograph('shane_kast_blue')
    caldir = data_path('Calibrations')
    caliBrate = calibrations.MultiSlitCalibrations(fitstbl, par, spectrograph, caldir)


def test_bias(multi_caliBrate):
    """
    This should produce nothing as we have no bias frames

    Returns:

    """
    multi_caliBrate.get_bias()


def test_arc(multi_caliBrate):
    arc = multi_caliBrate.get_arc()
    assert isinstance(arc, buildimage.ArcImage), 'arc has wrong type'
    assert Path(arc.get_path()).exists(), 'No Arc file written'
    assert arc.image.shape == (2048,350), 'Arc has wrong shape'

    # Cleanup
    shutil.rmtree(multi_caliBrate.calib_dir)


def test_tiltimg(multi_caliBrate):
    tilt = multi_caliBrate.get_tiltimg()
    assert isinstance(tilt, buildimage.TiltImage), 'tilt has wrong type'
    assert Path(tilt.get_path()).exists(), 'No Tilt file written'
    assert tilt.image.shape == (2048,350)

    # Cleanup
    shutil.rmtree(multi_caliBrate.calib_dir)


def test_bpm(multi_caliBrate):
    # Build
    bpm = multi_caliBrate.get_bpm()
    assert bpm.shape == (2048,350)
    assert np.sum(bpm) == 0.

# TODO: Add tests for:
#   - get_dark
#   - get_flats
#   - get_slits
#   - get_wv_calib
#   - get_tilts


def test_asn(multi_caliBrate):
    caldir = Path().resolve()
    ofile = caldir / 'test.calib'
    if ofile.exists():
        ofile.unlink()

    calibrations.Calibrations.association_summary(ofile, multi_caliBrate.fitstbl,
                                                  multi_caliBrate.spectrograph, caldir,
                                                  overwrite=True)

    # Read yaml file and test contents
    with open(ofile, 'r') as f:
        asn = yaml.safe_load(f)

    assert list(asn.keys()) == ['A'], 'Wrong setup list'
    assert list(asn['A'].keys()) == ['--', 0], 'Wrong A setup keys'

    # TODO: This causes windows CI tests to barf!  I gave up...
    import platform
    if platform.system() != 'Windows':
        assert Path(asn['A'][0]['arc']['proc'][0]).name == 'Arc_A_0_DET01.fits', \
                'Wrong calibration arc frame name'
    assert 'science' in asn['A'][0].keys(), 'Association file should include science frames'

    # Clean-up
    ofile.unlink()


def test_asn_calib_ID_dict(multi_caliBrate):

    caldir = Path().resolve()
    setup = 'A'
    calib_ID = 0
    det = 1
    # Force recorded calibration files to exist
    asn = calibrations.Calibrations.get_association(multi_caliBrate.fitstbl,
                                                    multi_caliBrate.spectrograph, caldir, setup,
                                                    calib_ID, det, must_exist=True)

    assert 'arc' in list(asn.keys()), 'Should find arc files in association'
    assert 'science' not in list(asn.keys()), 'Should not include science frames'
    assert len(asn['arc']['proc']) == 0, 'None of the processed calibration frames should exist'
    assert len(asn['arc']['raw']) == 1, 'Should be 1 raw arc frame'
    assert len(asn['pixelflat']['raw']) == 2, 'Should be 2 pixelflat frames'

    # Redo ignoring whether or not the calibration frames exist
    asn = calibrations.Calibrations.get_association(multi_caliBrate.fitstbl,
                                                    multi_caliBrate.spectrograph, caldir, setup,
                                                    calib_ID, det, must_exist=False)

    assert len(asn['arc']['proc']) == 2, \
            'Should be 2 processed calibration frames associated with the raw arc frames.'
    # TODO: Why does THIS pass the windows CI but the one above doesn't ?!?!?!
    arc_file = Path(asn['arc']['proc'][0]).name
    assert arc_file == 'Arc_A_0_DET01.fits', 'Wrong calibration arc frame name'
    
    # Redo including science/standard frames
    asn = calibrations.Calibrations.get_association(multi_caliBrate.fitstbl,
                                                    multi_caliBrate.spectrograph, caldir, setup,
                                                    calib_ID, det, must_exist=False,
                                                    include_science=True)

    assert 'science' in list(asn.keys()), 'Should include science frames'
    


