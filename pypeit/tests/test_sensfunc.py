import os

from IPython import embed

import numpy as np

from astropy.io import fits
from astropy import table

from pypeit import sensfunc
from pypeit.spectrographs.util import load_spectrograph
from pypeit.tests.tstutils import cooked_required

@cooked_required
def sensfunc_io_ir_test():
    sensfunc_io_tester('IR')

@cooked_required
def sensfunc_io_uvis_test():
    sensfunc_io_tester('UVIS')

def sensfunc_io_tester(algorithm):

    # Remove any existing file from previous runs that were interrupted
    test_file = 'test_sens.fits'
    if os.path.isfile(test_file):
        os.remove(test_file) 

    # Random spec1d file of a standard from cooked
    spec1dfile = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science',
                              'spec1d_b24-Feige66_KASTb_20150520T041246.960.fits')

    # Determine the spectrograph
    header = fits.getheader(spec1dfile)
    spectrograph = load_spectrograph(header['PYP_SPEC'])
    par = spectrograph.default_pypeit_par()
    par['sensfunc']['algorithm'] = algorithm

    # Instantiate the relevant class for the requested algorithm
    sensobj = sensfunc.SensFunc.get_instance(spec1dfile, test_file, par['sensfunc'])

    assert sensobj.std_dict['name'] == 'FEIGE66', 'incorrect standard star found'

    # Assign junk values just to test that I/O works
    sensobj.sens = sensobj.empty_sensfunc_table(*sensobj.wave_cnts.T.shape)
    sensobj.wave = sensobj.wave_cnts.copy()
    sensobj.zeropoint = sensobj.counts.copy()
    sensobj.throughput = np.ones(sensobj.counts.shape, dtype=float)

    assert sensobj.sens['SENS_WAVE'].shape == sensobj.wave_cnts.T.shape, 'shape is wrong'

    # Write it
    sensobj.to_file(test_file, overwrite=True)

    # Check the extensions
    hdu = fits.open(test_file)
    ext = [h.name for h in hdu]
    assert ext == ['PRIMARY', 'SENS', 'WAVE', 'ZEROPOINT', 'THROUGHPUT'], \
                'incorrect extensions written'
    del hdu

    # Read it back in and check that the data are the same
    _sensobj = sensfunc.SensFunc.from_file(test_file)

    assert np.array_equal(_sensobj.wave, sensobj.wave), 'I/O error'
    assert np.array_equal(_sensobj.zeropoint, sensobj.zeropoint), 'I/O error'
    assert isinstance(_sensobj.sens, table.Table), 'sens table has wrong type'
    assert _sensobj.algorithm == sensobj.algorithm, \
        'algorithm mismatch after writing to and reading from disk'

    os.remove(test_file) 
