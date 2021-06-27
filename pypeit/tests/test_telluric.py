import os
from pkg_resources import resource_filename

from IPython import embed

import numpy as np

from pypeit.core import telluric
from pypeit.tests.tstutils import telluric_required, tell_test_grid

@telluric_required
def test_telluric_init():

    wave = np.linspace(8250, 25200, 4400)
    flux = np.ones(wave.size, dtype=float)
    ivar = np.ones(wave.size, dtype=float)
    gpm = np.ones(wave.size, dtype=bool)

    obj_params = {'pca_file': os.path.join(resource_filename('pypeit', 'data'), 'telluric',
                                           'models', 'qso_pca_1200_3100.fits'),
                  'npca': 8, 'z_qso': 7.52, 'delta_zqso': 0.1, 'lbound_norm': 0.1,
                  'ubound_norm': 3.0, 'tell_norm_thresh': 0.9,
                  'output_meta_keys': ('pca_file', 'npca', 'z_qso', 'delta_zqso', 'lbound_norm',
                                       'ubound_norm', 'tell_norm_thresh'),
                  'debug_init': False}


    t = telluric.Telluric(wave, flux, ivar, gpm, tell_test_grid, obj_params,
                          telluric.init_qso_model, telluric.eval_qso_model)

    assert t.model['WAVE'].shape == (1,4400), 'bad wavelength array shape'
    assert t.model['TELL_THETA'].shape == (1,7), 'bad or changed telluric parameter shape'
    assert t.model['OBJ_THETA'].shape == (1,9), 'bad or changed object parameter shape'
    assert t.exptime is None, \
                'exptime was not part of obj_params and output_meta_keys, so it should ' \
                'not have been defined.'
    assert t.z_qso == 7.52, 'QSO redshift incorrect'
    

@telluric_required
def test_telluric_io():

    # Remove any existing file from previous runs that were interrupted
    test_file = 'test_telluric.fits'
    if os.path.isfile(test_file):
        os.remove(test_file) 

    wave = np.linspace(8250, 25200, 4400)
    flux = np.ones(wave.size, dtype=float)
    ivar = np.ones(wave.size, dtype=float)
    gpm = np.ones(wave.size, dtype=bool)

    obj_params = {'pca_file': os.path.join(resource_filename('pypeit', 'data'), 'telluric',
                                           'models', 'qso_pca_1200_3100.fits'),
                  'npca': 8, 'z_qso': 7.52, 'delta_zqso': 0.1, 'lbound_norm': 0.1,
                  'ubound_norm': 3.0, 'tell_norm_thresh': 0.9,
                  'output_meta_keys': ('pca_file', 'npca', 'z_qso', 'delta_zqso', 'lbound_norm',
                                       'ubound_norm', 'tell_norm_thresh'),
                  'debug_init': False}

    t = telluric.Telluric(wave, flux, ivar, gpm, tell_test_grid, obj_params,
                          telluric.init_qso_model, telluric.eval_qso_model)
    t.to_file(test_file)
    _t = telluric.Telluric.from_file(test_file)

    assert t.telgrid == _t.telgrid, 'telgrid file name changed'
    assert t.pca_file == _t.pca_file, 'pca_file changed'
    assert t.z_qso == _t.z_qso, 'QSO redshift changed'
    assert t.model.keys() == _t.model.keys(), 'model table columns changed'
    assert np.all(t.model == _t.model), 'model data changed'

    # Clean-up
    os.remove(test_file) 

# TODO: Additional tests?

