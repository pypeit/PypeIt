# Module to run tests on WaveTilts class
#   Requires files in Development suite and an Environmental variable
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# TEST_UNICODE_LITERALS

import os

import pytest
import numpy as np


from pypeit.tests.tstutils import dev_suite_required, load_kast_blue_masters
from pypeit import wavetilts
from pypeit.core import tracewave, pixels
from pypeit.par import pypeitpar


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


@pytest.fixture
@dev_suite_required
def master_dir():
    # Any test that uses this directory also requires the DevSuite!
#    return data_path('MF_shane_kast_blue') if os.getenv('PYPEIT_DEV') is None \
#            else os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'MF_shane_kast_blue')
    return os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'MF_shane_kast_blue')


@dev_suite_required
def test_step_by_step(master_dir):
    # Masters
    spectrograph, msarc, tslits_dict = load_kast_blue_masters(get_spectrograph=True, aimg=True,tslits=True)
    # Instantiate
    master_key = 'A_1_01'
    parset = spectrograph.default_pypeit_par()
    par = parset['calibrations']['tilts']
    waveTilts = wavetilts.WaveTilts(msarc, spectrograph=spectrograph, par=par, det=1, master_key=master_key,
                                    master_dir=master_dir,reuse_masters=True, tslits_dict=tslits_dict)
    # Extract arcs
    arccen, maskslits = waveTilts.extract_arcs(waveTilts.slitcen, waveTilts.slitmask, msarc, waveTilts.inmask)
    assert arccen.shape == (2048,1)
    # Tilts in the slit
    slit = 0
    waveTilts.slitmask = pixels.tslits2mask(waveTilts.tslits_dict)
    thismask = waveTilts.slitmask == slit
    waveTilts.lines_spec, waveTilts.lines_spat = waveTilts.find_lines(arccen[:, slit], waveTilts.slitcen[:, slit], slit)

    trcdict = waveTilts.trace_tilts(waveTilts.msarc, waveTilts.lines_spec, waveTilts.lines_spat, thismask, slit)
    assert isinstance(trcdict, dict)
    # 2D Fit
    spat_order = waveTilts._parse_param(waveTilts.par, 'spat_order', slit)
    spec_order = waveTilts._parse_param(waveTilts.par, 'spec_order', slit)
    coeffs = waveTilts.fit_tilts(trcdict, thismask, waveTilts.slitcen[:, slit], spat_order, spec_order,slit, doqa=False)
    tilts = tracewave.fit2tilts(waveTilts.slitmask_science.shape, coeffs, waveTilts.par['func2d'])
    assert np.max(tilts) < 1.01


@dev_suite_required
def test_run(master_dir):
    # Masters
    spectrograph, msarc, tslits_dict = load_kast_blue_masters(get_spectrograph=True, aimg=True,
                                                         tslits=True)
    # Instantiate
    master_key = 'A_1_01'
    spectrograph.detector[0]['saturation'] = 60000.
    spectrograph.detector[0]['nonlinear'] = 0.9
    par = pypeitpar.WaveTiltsPar()
    waveTilts = wavetilts.WaveTilts(msarc, spectrograph=spectrograph, par=par, det=1, master_key=master_key,
                                    master_dir=master_dir, reuse_masters=True, tslits_dict=tslits_dict)
    # Run
    tilts_dict, mask = waveTilts.run(doqa=False)
    assert isinstance(tilts_dict['tilts'], np.ndarray)

