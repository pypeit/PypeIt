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
#    specdir = os.path.abspath(__file__)
#    for i in range(3):
#        specdir = os.path.split(specdir)[0]
#    usr_example = os.path.join(specdir, 'doc', 'nb', 'keck_lris_blue_long_400_3400_d560.cfg')
    usr_example = 'usr_merge_test.cfg'
    p = pypitpar.PypitPar.from_cfg_file(merge_with=usr_example)
    assert p['rdx']['spectrograph'] == 'KECK_LRISb', 'Test spectrograph is incorrect!'
    assert p['rdx']['pipeline'] == 'ARMS', 'Test pipeline is incorrect!'
    assert p['biasgroup']['useframe'] == 'overscan', 'Test biasgroup:useframe is incorrect!'

