# Module to run tests on WaveTilts class
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

from pypit.tests import tstutils
from pypit import wavetilts
from pypit.par import pypitpar

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


def test_step_by_step():
    if skip_test:
        assert True
        return
    # Masters
    spectrograph, msarc, TSlits \
            = tstutils.load_kast_blue_masters(get_spectrograph=True, aimg=True, tslits=True)
    # Instantiate
    setup = 'A_01_aa'
    spectrograph.detector[0]['saturation'] = 60000.
    spectrograph.detector[0]['nonlinear'] = 0.9
    root_path = data_path('MF') if os.getenv('PYPIT_DEV') is None \
                    else os.path.join(os.getenv('PYPIT_DEV'), 'Cooked', 'MF')
    par = pypitpar.TraceTiltsPar()
    waveTilts = wavetilts.WaveTilts(msarc, spectrograph=spectrograph, par=par, det=1, setup=setup,
                                    root_path=root_path, mode='reuse', pixlocn=TSlits.pixlocn,
                                    tslits_dict=TSlits.tslits_dict)
    # Extract arcs
    arccen, maskslits = waveTilts._extract_arcs()
    assert arccen.shape == (2048,1)
    # Tilts in the slit
    slit = 0
    waveTilts.par['tracethresh'] = 500.   # Lowers the thershold amplitude of the arc lines used
    trcdict = waveTilts._trace_tilts(slit)
    assert isinstance(trcdict, dict)
    # Analyze
    waveTilts.par['order'] = 3
    badrows = waveTilts._analyze_lines(slit)
    assert badrows == 0
    # 2D Fit
    waveTilts.par['yorder'] = 4
    waveTilts.par['func2D'] = 'legendre'
    tilts = waveTilts._fit_tilts(slit, doqa=False)
    assert np.max(tilts) < 1.0001


def test_run():
    if skip_test:
        assert True
        return

    # Masters
    spectrograph, msarc, TSlits \
            = tstutils.load_kast_blue_masters(get_spectrograph=True, aimg=True, tslits=True)
    # Instantiate
    setup = 'A_01_aa'
    spectrograph.detector[0]['saturation'] = 60000.
    spectrograph.detector[0]['nonlinear'] = 0.9
    root_path = data_path('MF') if os.getenv('PYPIT_DEV') is None \
                    else os.path.join(os.getenv('PYPIT_DEV'), 'Cooked', 'MF')
    par = pypitpar.TraceTiltsPar()
    waveTilts = wavetilts.WaveTilts(msarc, spectrograph=spectrograph, par=par, det=1, setup=setup,
                                    root_path=root_path, mode='reuse', pixlocn=TSlits.pixlocn,
                                    tslits_dict=TSlits.tslits_dict)
    # Run
    tilts, mask = waveTilts.run(doqa=False)
    assert isinstance(tilts, np.ndarray)

