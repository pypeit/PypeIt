"""
Module to run tests on WaveTilts class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import numpy as np


from pypeit.tests.tstutils import dev_suite_required, load_kast_blue_masters, cooked_required
from pypeit import wavetilts
from pypeit.core import tracewave, pixels
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


@pytest.fixture
@cooked_required
def master_dir():
    return os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Shane_Kast_blue')

@cooked_required
def test_instantiate_from_master(master_dir):
    master_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Shane_Kast_blue', 'MasterTilts_A_1_01.fits')
    waveTilts = wavetilts.WaveTilts.from_master_file(master_file)
    assert isinstance(waveTilts.tilts_dict, dict)

@cooked_required
def test_step_by_step(master_dir):
    # Masters
    spectrograph = load_spectrograph('shane_kast_blue')
    msarc, edges = load_kast_blue_masters(aimg=True, edges=True)
    # Instantiate
    master_key = 'A_1_01'
    parset = spectrograph.default_pypeit_par()
    par = parset['calibrations']['tilts']
    wavepar = parset['calibrations']['wavelengths']
    waveTilts = wavetilts.WaveTilts(msarc, edges.convert_to_tslits_dict(), spectrograph, par,
                                    wavepar, det=1, master_key=master_key, master_dir=master_dir,
                                    reuse_masters=True)
    # Extract arcs
    arccen, arccen_bpm, maskslits = waveTilts.extract_arcs()#waveTilts.slitcen, waveTilts.slitmask, waveTilts.inmask)
    assert arccen.shape == (2048,1)
    # Tilts in the slit
    slit = 0
    waveTilts.slitmask = pixels.tslits2mask(waveTilts.tslits_dict)
    thismask = waveTilts.slitmask == slit
    waveTilts.lines_spec, waveTilts.lines_spat = waveTilts.find_lines(arccen[:, slit], waveTilts.slitcen[:, slit], slit)

    trcdict = waveTilts.trace_tilts(waveTilts.msarc.image, waveTilts.lines_spec,
                                    waveTilts.lines_spat, thismask, slit)
    assert isinstance(trcdict, dict)
    # 2D Fit
    spat_order = waveTilts._parse_param(waveTilts.par, 'spat_order', slit)
    spec_order = waveTilts._parse_param(waveTilts.par, 'spec_order', slit)
    coeffs = waveTilts.fit_tilts(trcdict, thismask, waveTilts.slitcen[:, slit], spat_order, spec_order,slit, doqa=False)
    tilts = tracewave.fit2tilts(waveTilts.slitmask_science.shape, coeffs, waveTilts.par['func2d'])
    assert np.max(tilts) < 1.01


@cooked_required
def test_run(master_dir):
    # Masters
    spectrograph = load_spectrograph('shane_kast_blue')
    msarc, edges = load_kast_blue_masters(aimg=True, edges=True)
    # Instantiate
    master_key = 'A_1_01'
    spectrograph.detector[0]['saturation'] = 60000.
    spectrograph.detector[0]['nonlinear'] = 0.9
    par = pypeitpar.WaveTiltsPar()
    wavepar = pypeitpar.WavelengthSolutionPar()
    waveTilts = wavetilts.WaveTilts(msarc, edges.convert_to_tslits_dict(), spectrograph, par,
                                    wavepar, det=1, master_key=master_key, master_dir=master_dir,
                                    reuse_masters=True)
    # Run
    tilts_dict, mask = waveTilts.run(doqa=False)
    assert isinstance(tilts_dict['tilts'], np.ndarray)

