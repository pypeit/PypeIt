# Module to run tests on WaveCalib class
#   Requires files in Development suite and an Environmental variable
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# TEST_UNICODE_LITERALS

import os

import pytest
import glob
import numpy as np

from astropy.table import Table

from pypit import wavecalib
from pypit import traceslits
from pypit import arcimage

# These tests are not run on Travis
if os.getenv('PYPIT_DEV') is None:
    skip_test=True
else:
    skip_test=False

def chk_for_files(root):
    files = glob.glob(root+'*')
    if len(files) == 0:
        return False
    else:
        return True

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_user_redo():
    if skip_test:
        assert True
        return
    # Check for files
    wvcalib_file = os.path.join(os.getenv('PYPIT_DEV'), 'Cooked', 'WaveCalib',
                                'MasterWaveCalib_ShaneKastBlue_A.json')
    assert chk_for_files(wvcalib_file)
    # Instantiate
    waveCalib = wavecalib.WaveCalib(None, spectrograph='shane_kast_blue')
    waveCalib.load_wv_calib(wvcalib_file)
    # Do it
    waveCalib.arcparam['min_ampl'] = 1000.
    new_wv_calib = waveCalib.calibrate_spec(0)
    # Test
    assert new_wv_calib['rms'] < 0.035

def test_step_by_step():
    if skip_test:
        assert True
        return

    root_path = os.path.join(os.getenv('PYPIT_DEV'), 'Cooked', 'MF')
    setup = 'A_01_aa'

    # Load up the Masters
    AImg = arcimage.ArcImage('shane_kast_blue', setup=setup, root_path=root_path, mode='reuse')
    msarc, header, _ = AImg.load_master_frame()
    TSlits = traceslits.TraceSlits.from_master_files(os.path.join(AImg.directory_path,
                                                                  'MasterTrace_A_01_aa'))
    TSlits._make_pixel_arrays()
    fitstbl = Table.read(os.path.join(AImg.directory_path, 'shane_kast_blue_setup_A.fits'))

    # Instantiate
    waveCalib = wavecalib.WaveCalib(msarc, spectrograph='shane_kast_blue', setup=setup,
                                    root_path=root_path, mode='reuse', fitstbl=fitstbl,
                                    sci_ID=1, det=1)
    # Extract arcs
    arccen, maskslits = waveCalib._extract_arcs(TSlits.lcen, TSlits.rcen, TSlits.pixlocn)
    assert arccen.shape == (2048,1)
    # Arcparam
    arcparam = waveCalib._load_arcparam()
    assert isinstance(arcparam, dict)
    # wv_calib
    wv_calib = waveCalib._build_wv_calib('arclines', skip_QA=True)
    assert isinstance(wv_calib, dict)
    # Master
    # This doesn't save steps nor arcparam which *is* done in the master() call
    waveCalib.save_master(wv_calib, outfile=data_path('tmp.json'))


def test_one_shot():
    if skip_test:
        assert True
        return

    root_path = os.path.join(os.getenv('PYPIT_DEV'), 'Cooked', 'MF')
    setup = 'A_01_aa'

    # Load up the Masters
    AImg = arcimage.ArcImage('shane_kast_blue', setup=setup, root_path=root_path, mode='reuse')
    msarc, header, _ = AImg.load_master_frame()
    TSlits = traceslits.TraceSlits.from_master_files(os.path.join(AImg.directory_path,
                                                                  'MasterTrace_A_01_aa'))
    TSlits._make_pixel_arrays()
    fitstbl = Table.read(os.path.join(AImg.directory_path, 'shane_kast_blue_setup_A.fits'))
    # Do it
    waveCalib = wavecalib.WaveCalib(msarc, spectrograph='shane_kast_blue', setup=setup,
                                    root_path=root_path, fitstbl=fitstbl, sci_ID=1, det=1)
    wv_calib2, _ = waveCalib.run(TSlits.lcen, TSlits.rcen, TSlits.pixlocn, skip_QA=True)
    #
    assert 'arcparam' in wv_calib2.keys()


