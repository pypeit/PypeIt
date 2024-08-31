from pathlib import Path
import yaml
import pytest
import shutil

from IPython import embed

import numpy as np

from pypeit import dataPaths
from pypeit import calibrations
from pypeit import pypeitsetup
from pypeit.images import buildimage
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph

from pypeit.tests.tstutils import data_output_path

@pytest.fixture
def fitstbl():

    # Get the files
    file_names = [
        'b1.fits.gz',    # arc
        'b11.fits.gz',   # trace
        'b21.fits.gz',   # bias
        'b24.fits.gz',   # standard
        'b27.fits.gz'    # science
    ]
    files = [dataPaths.tests.get_file_path(f, to_pkg='symlink') for f in file_names]

    setupc = pypeitsetup.PypeItSetup(files, spectrograph_name='shane_kast_blue')
    setupc.build_fitstbl(files)
    setupc.fitstbl.finalize_usr_build(None, 'A')
    return setupc.fitstbl

@pytest.fixture
def multi_caliBrate(fitstbl):
    # Grab a science file for configuration specific parameters
    indx = fitstbl.find_frames('science', index=True)[0]
    sci_file = fitstbl.frame_paths(indx)
    # Par
    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.config_specific_par(sci_file)
    par.reset_all_processimages_par(use_biasimage=False)
    #
    calib_par = par['calibrations']
    calib_par['bpm_usebias'] = False
    calib_par['slitedges']['sync_predict'] = 'nearest'

    multi_caliBrate = calibrations.MultiSlitCalibrations(fitstbl, calib_par, spectrograph,
                                                         data_output_path('Calibrations'))
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
    caldir = data_output_path('Calibrations')
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
    caldir = data_output_path('Calibrations')
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
    caldir = Path().absolute()
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

    caldir = Path().absolute()
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
    assert len(asn['pixelflat']['raw']) == 1, 'Should be 1 pixelflat frames'

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
    


