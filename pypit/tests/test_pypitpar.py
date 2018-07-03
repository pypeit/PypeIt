# Module to run tests on PypitPar classes
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import pytest

from pypit.par import pypitpar

def test_run():
    pypitpar.RunPar()

def test_overscan():
    pypitpar.OverscanPar()

def test_flatfield():
    pypitpar.FlatFieldPar()

def test_flexure():
    pypitpar.FlexurePar()

def test_wavelengthcalibration():
    pypitpar.WavelengthCalibrationPar()

def test_fluxcalibration():
    pypitpar.FluxCalibrationPar()

def test_skysubtraction():
    pypitpar.SkySubtractionPar()

def test_spectrographs():
    s = pypitpar.ReducePar.valid_spectrographs()
    assert 'KECK_LRISb' in s, 'Expected to find KECK_LRISb as a spectrograph!'

def test_reduce():
    pypitpar.ReducePar()

def test_combine():
    pypitpar.CombineFramesPar()

def test_framegroup_types():
    t = pypitpar.FrameGroupPar.valid_frame_types()
    assert 'bias' in t, 'Expected to find \'bias\' in list of valid frame types'

def test_framegroup():
    pypitpar.FrameGroupPar()

def test_wavelengthsolution():
    pypitpar.WavelengthSolutionPar()

def test_pca():
    pypitpar.PCAPar()

def test_traceslits():
    pypitpar.TraceSlitsPar()

def test_tracetilts():
    pypitpar.TraceTiltsPar()

def test_traceobjects():
    pypitpar.TraceObjectsPar()

def test_manualextraction():
    pypitpar.ManualExtractionPar()

def test_extractobjects():
    pypitpar.ExtractObjectsPar()

def test_detector():
    pypitpar.DetectorPar()

def test_instrument():
    pypitpar.InstrumentPar()

def test_framefits():
    pypitpar.FrameFitsPar()

def test_frameid():
    pypitpar.FrameIDPar()

def test_pypit():
    pypitpar.PypitPar()

def test_fromcfgfile():
    pypitpar.PypitPar.from_cfg_file()

def test_writecfg():
    default_file = 'default.cfg'
    if os.path.isfile(default_file):
        os.remove(default_file)
    pypitpar.PypitPar.from_cfg_file().to_config(default_file)
    assert os.path.isfile(default_file), 'No file written!'
    os.remove(default_file)

def test_readcfg():
    default_file = 'default.cfg'
    if not os.path.isfile(default_file):
        pypitpar.PypitPar.from_cfg_file().to_config(default_file)
    pypitpar.PypitPar.from_cfg_file('default.cfg', expand_spectrograph=False)
    os.remove(default_file)

def test_mergecfg():
    # Create a file with the defaults
    user_file = 'user_adjust.cfg'
    if os.path.isfile(user_file):
        os.remove(user_file)
    p = pypitpar.PypitPar.from_cfg_file()

    # Make some modifications
    p['rdx']['spectrograph'] = 'KECK_LRISb'
    p['rdx']['pipeline'] = 'ARMS'
    p['biasgroup']['useframe'] = 'overscan'

    # Write the modified config
    p.to_config(user_file)
    assert os.path.isfile(user_file), 'User file was not written!'

    # Read it in as a config to merge with the defaults
    p = pypitpar.PypitPar.from_cfg_file(merge_with=user_file)

    # Check the values are correctly read in
    assert p['rdx']['spectrograph'] == 'KECK_LRISb', 'Test spectrograph is incorrect!'
    assert p['rdx']['pipeline'] == 'ARMS', 'Test pipeline is incorrect!'
    assert p['biasgroup']['useframe'] == 'overscan', 'Test biasgroup:useframe is incorrect!'

    # Clean-up
    os.remove(user_file)

