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

from pypeit import wavecalib
from pypeit import traceslits
from pypeit import arcimage

from pypeit.core.wavecal import autoid
import json
import h5py

# These tests are not run on Travis
if os.getenv('PYPEIT_DEV') is None:
    skip_test = True
else:
    skip_test = False


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files/wavecalib')
    return os.path.join(data_dir, filename)

master_dir = data_path('MF_shane_kast_blue') if os.getenv('PYPEIT_DEV') is None \
    else os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'MF_shane_kast_blue')

def chk_for_files(root):
    files = glob.glob(root+'*')
    if len(files) == 0:
        return False
    else:
        return True


def test_user_redo():
    if skip_test:
        assert True
        return
    # Check for files
    wvcalib_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'WaveCalib',
                                'MasterWaveCalib_ShaneKastBlue_A.json')
    assert chk_for_files(wvcalib_file)
    # Instantiate
    waveCalib = wavecalib.WaveCalib(None, spectrograph='shane_kast_blue')
    waveCalib.load_wv_calib(wvcalib_file)
    # Do it
    waveCalib.arcparam['min_ampl'] = 1000.
    new_wv_calib = waveCalib.calibrate_spec(0)
    # Test
    assert new_wv_calib['0']['rms'] < 0.035


def test_step_by_step():
    if skip_test:
        assert True
        return

    root_path = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'MF')
    setup = 'A_01_aa'

    # Load up the Masters
    AImg = arcimage.ArcImage('shane_kast_blue', setup=setup, master_dir=master_dir, mode='reuse')
    msarc, header, _ = AImg.load_master_frame()
    TSlits = traceslits.TraceSlits.from_master_files(os.path.join(master_dir,
                                                                  'MasterTrace_A_01_aa'))
    TSlits._make_pixel_arrays()
    fitstbl = Table.read(os.path.join(master_dir, 'shane_kast_blue_setup_A.fits'))

    # Instantiate
    waveCalib = wavecalib.WaveCalib(msarc, spectrograph='shane_kast_blue', setup=setup,
                                    master_dir=master_dir, mode='reuse', fitstbl=fitstbl,
                                    sci_ID=1, det=1)
    # Extract arcs
    arccen, maskslits = waveCalib._extract_arcs(TSlits.lcen, TSlits.rcen, TSlits.pixlocn)
    assert arccen.shape == (2048, 1)
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

    root_path = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'MF')
    setup = 'A_01_aa'

    # Load up the Masters
    AImg = arcimage.ArcImage('shane_kast_blue', setup=setup, master_dir=master_dir, mode='reuse')
    msarc, header, _ = AImg.load_master_frame()
    TSlits = traceslits.TraceSlits.from_master_files(os.path.join(master_dir,
                                                                  'MasterTrace_A_01_aa'))
    TSlits._make_pixel_arrays()
    fitstbl = Table.read(os.path.join(master_dir, 'shane_kast_blue_setup_A.fits'))
    # Do it
    waveCalib = wavecalib.WaveCalib(msarc, spectrograph='shane_kast_blue', setup=setup,
                                    master_dir=master_dir, fitstbl=fitstbl, sci_ID=1, det=1)
    wv_calib2, _ = waveCalib.run(TSlits.lcen, TSlits.rcen, TSlits.pixlocn, skip_QA=True)
    #
    assert 'arcparam' in wv_calib2.keys()


def test_wavecalib_general():

    # LRISb 600/4000 with the longslit
    names = ['LRISb_600_4000_longslit']
    spec_files = ['lrisb_600_4000_PYPIT.json']
    all_lines = [['CdI', 'HgI', 'ZnI']]
    all_wvcen = [4400.]
    all_disp = [1.26]
    fidxs = [0]
    scores = [dict(rms=0.1, nxfit=13, nmatch=10)]

    # LRISb 400/3400 with the longslit
    names += ['LRISb_400_3400_longslit']
    spec_files += ['lrisb_400_3400_PYPIT.json']
    all_lines += [['NeI', 'ArI', 'CdI', 'KrI', 'XeI', 'ZnI', 'HgI']]
    all_wvcen += [4400.]
    all_disp += [1.26]
    fidxs += [0]
    scores += [dict(rms=0.1, nxfit=13, nmatch=10)]

    '''
    # LRISb off-center
    # NeI lines are seen beyond the dicrhoic but are ignored
    names += ['LRISb_600_4000_red']
    src_files += ['LRISb_600_LRX.hdf5']  # Create with low_redux.py if needed
    all_wvcen += [5000.]
    all_disp += [1.26]
    all_lines += [['CdI', 'ZnI', 'HgI']]  #,'NeI']]
    fidxs += [18]
    scores += [dict(rms=0.08, nxfit=10, nmatch=10)]
    '''

    # LRISr 600/7500 longslit
    names += ['LRISr_600_7500_longslit']
    spec_files += ['lrisr_600_7500_PYPIT.json']
    all_wvcen += [7000.]
    all_disp += [1.6]
    all_lines += [['ArI', 'HgI', 'KrI', 'NeI', 'XeI']]
    fidxs += [-1]
    scores += [dict(rms=0.08, nxfit=30, nmatch=50)]

    '''
    # LRISr 900/XX00 longslit -- blue
    names += ['LRISr_900_XX00_longslit']
    src_files += ['lrisr_900_XX00_PYPIT.json']
    all_wvcen += [5800.]
    all_disp += [1.08]
    all_lines += [['ArI', 'HgI', 'KrI', 'NeI', 'XeI', 'CdI', 'ZnI']]
    fidxs += [-1]
    scores += [dict(rms=0.08, nxfit=10, nmatch=10)]
    '''

    # LRISr 400/8500 longslit -- red
    names += ['LRISr_400_8500_longslit']
    spec_files += ['lrisr_400_8500_PYPIT.json']
    all_wvcen += [8000.]
    all_disp += [2.382]
    all_lines += [['ArI', 'HgI', 'KrI', 'NeI', 'XeI']]
    fidxs += [-1]
    scores += [dict(rms=0.1, nxfit=40, nmatch=40)]

    # Kastb 600 grism
    names += ['KASTb_600_standard']
    spec_files += ['kastb_600_PYPIT.json']
    all_lines += [['CdI', 'HeI', 'HgI']]
    all_wvcen += [4400.]
    all_disp += [1.02]
    fidxs += [0]
    scores += [dict(rms=0.1, nxfit=13, nmatch=8)]

    # Kastr 600/7500 grating
    names += ['KASTr_600_7500_standard']
    spec_files += ['kastr_600_7500_PYPIT.json']
    all_lines += [['ArI', 'NeI', 'HgI']]
    all_wvcen += [6800.]
    all_disp += [2.345]
    fidxs += [0]
    scores += [dict(rms=0.1, nxfit=10, nmatch=8)]

    # Keck DEIMOS
    names += ['keck_deimos_830g_l']
    spec_files += ['keck_deimos_830g_l_PYPIT.json']
    all_lines += [['ArI', 'NeI', 'KrI', 'XeI']]
    all_wvcen += [7450.]
    all_disp += [0.467]
    fidxs += [0]
    scores += [dict(rms=0.1, nxfit=20, nmatch=20)]


    # Favored parameters (should match those in the defaults)
    min_ampl = 1000.

    # Run it
    for name, spec_file, lines, wvcen, disp, score, fidx in zip(
            names, spec_files, all_lines, all_wvcen, all_disp, scores, fidxs):
        # Load spectrum
        exten = spec_file.split('.')[-1]
        if exten == 'json':
            with open(data_path(spec_file), 'r') as f:
                pypit_fit = json.load(f)
            try:
                # Old format
                spec = np.array(pypit_fit['spec'])
            except KeyError:
                spec = np.array(pypit_fit['0']['spec'])
        elif exten == 'hdf5':
            hdf = h5py.File(data_path(spec_file), 'r')
            spec = hdf['arcs/{:d}/spec'.format(fidx)].value

        arcfitter = autoid.General(spec.reshape((spec.size, 1)), lines, min_ampl=min_ampl, rms_threshold=score['rms'])
        final_fit = arcfitter._all_final_fit

        # Score
        grade = True
        slit = '0'
        if final_fit[slit]['rms'] > score['rms']:
            grade = False
            print("Solution for {:s} failed RMS!!".format(name))
        if len(final_fit[slit]['xfit']) < score['nxfit']:
            grade = False
            print("Solution for {:s} failed N xfit!!".format(name))
#        if patt_dict[slit]['nmatch'] < score['nmatch']:
#            grade = False
#            print("Solution for {:s} failed N match!!".format(name))
        assert grade

