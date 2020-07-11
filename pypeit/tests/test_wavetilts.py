"""
Module to run tests on WaveTilts and BuildWaveTilts classes
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import numpy as np


from pypeit.tests.tstutils import dev_suite_required, load_kast_blue_masters, cooked_required
from pypeit import wavetilts
from pypeit import slittrace
from pypeit.core import tracewave, pixels
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


@pytest.fixture
@cooked_required
def master_dir():
    return os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'shane_kast_blue')

instant_dict = dict(coeffs=np.ones((6,4,1)),
                    nslit=1,
                    spat_order=np.array([3]),
                    spec_order=np.array([5]),
                    spat_id=np.array([150]),
                    func2d='legendre2d')

# Test WaveTilts
def test_wavetilts():
    #
    wvtilts = wavetilts.WaveTilts(**instant_dict)
    # I/O
    outfile = data_path('tst_wavetilts.fits')
    wvtilts.to_file(outfile, overwrite=True)
    _wvtilts = wavetilts.WaveTilts.from_file(outfile)

    # Test
    for key in instant_dict.keys():
        if isinstance(instant_dict[key], np.ndarray):
            assert np.array_equal(wvtilts[key],_wvtilts[key])
        else:
            assert wvtilts[key] == _wvtilts[key]
    # Write again
    wvtilts.to_file(outfile, overwrite=True)
    os.remove(outfile)


@cooked_required
def test_instantiate_from_master(master_dir):
    master_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'shane_kast_blue',
                               'MasterTilts_A_1_01.fits')
    slit_master_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'shane_kast_blue',
                               'MasterSlits_A_1_01.fits.gz')
    slits = slittrace.SlitTraceSet.from_file(slit_master_file)
    waveTilts = wavetilts.WaveTilts.from_file(master_file)
    tilts = waveTilts.fit2tiltimg(slits.slit_img())
    assert isinstance(tilts, np.ndarray)


# Test rebuild tilts with a flexure offset
@cooked_required
def test_flexure(master_dir):
    flexure = 1.
    master_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'shane_kast_blue', 'MasterTilts_A_1_01.fits')
    #master_file = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT', 'shane_kast_blue', '600_4310_d55',
    #                           'shane_kast_blue_A','Masters',
    #                           'MasterTilts_A_1_01.fits')
    waveTilts = wavetilts.WaveTilts.from_file(master_file)
    # Need slitmask
    slit_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'shane_kast_blue', 'MasterSlits_A_1_01.fits.gz')
    #slit_file =os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT', 'shane_kast_blue', '600_4310_d55',
    #                                                  'shane_kast_blue_A','Masters',  'MasterSlits_A_1_01.fits.gz')
    slits = slittrace.SlitTraceSet.from_file(slit_file)
    slitmask = slits.slit_img(flexure=flexure)
    # Do it
    new_tilts = waveTilts.fit2tiltimg(slitmask, flexure=flexure)
    # Test?

# Test BuildWaveTilts
#  To maintain this test means copying in the code from run over-and-again.  No thanks..
#@cooked_required
#def test_step_by_step(master_dir):
#    # Masters
#    spectrograph = load_spectrograph('shane_kast_blue')
#    mstilt, edges = load_kast_blue_masters(mstilt=True, edges=True)
#    # Instantiate
#    parset = spectrograph.default_pypeit_par()
#    par = parset['calibrations']['tilts']
#    wavepar = parset['calibrations']['wavelengths']
#    buildwaveTilts = wavetilts.BuildWaveTilts(mstilt, edges.get_slits(), spectrograph, par, wavepar, det=1)
#    # Extract arcs
#    arccen, arccen_bpm = buildwaveTilts.extract_arcs()
#    assert arccen.shape == (2048,1)
#    # Tilts in the slit
#    slit_idx = 0
#    buildwaveTilts.slitmask = buildwaveTilts.slits.slit_img()
#    thismask = buildwaveTilts.slitmask == buildwaveTilts.slits.spat_id[slit_idx]
#    buildwaveTilts.lines_spec, buildwaveTilts.lines_spat \
#            = buildwaveTilts.find_lines(arccen[:, slit_idx], buildwaveTilts.slitcen[:, slit_idx], slit)
#
#    trcdict = buildwaveTilts.trace_tilts(buildwaveTilts.mstilt.image, buildwaveTilts.lines_spec,
#                                    buildwaveTilts.lines_spat, thismask, slit)
#    assert isinstance(trcdict, dict)
#    # 2D Fit
#    spat_order = buildwaveTilts._parse_param(buildwaveTilts.par, 'spat_order', slit)
#    spec_order = buildwaveTilts._parse_param(buildwaveTilts.par, 'spec_order', slit)
#    coeffs = buildwaveTilts.fit_tilts(trcdict, thismask, buildwaveTilts.slitcen[:, slit], spat_order,
#                                 spec_order,slit, doqa=False)
#    tilts = tracewave.fit2tilts(buildwaveTilts.slitmask_science.shape, coeffs, buildwaveTilts.par['func2d'])
#    assert np.max(tilts) < 1.01


@cooked_required
def test_run(master_dir):
    # Masters
    spectrograph = load_spectrograph('shane_kast_blue')
    mstilt, edges = load_kast_blue_masters(mstilt=True, edges=True)
    # Instantiate
    #spectrograph.detector[0]['saturation'] = 60000.
    #spectrograph.detector[0]['nonlinear'] = 0.9
    par = pypeitpar.WaveTiltsPar()
    wavepar = pypeitpar.WavelengthSolutionPar()
    slits = edges.get_slits()
    buildwaveTilts = wavetilts.BuildWaveTilts(mstilt, slits, spectrograph, par, wavepar, det=1)
    # Run
    waveTilts = buildwaveTilts.run(doqa=False)
    assert isinstance(waveTilts.fit2tiltimg(slits.slit_img()), np.ndarray)


#if __name__ == '__main__':
#    test_instantiate_from_master(master_dir())
