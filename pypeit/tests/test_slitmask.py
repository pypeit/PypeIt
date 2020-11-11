"""
Module to run tests on PypeItPar classes
"""
import os
import numpy

import pytest

from astropy.io import fits

from pypeit import slittrace
from pypeit import specobjs
from pypeit.spectrographs.keck_deimos import KeckDEIMOSSpectrograph
from pypeit.tests.tstutils import dev_suite_required

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_assign_maskinfo():

    # Spectrograph
    spec = KeckDEIMOSSpectrograph()
    det = 3
    maskdef_file = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_deimos',
                                '830G_M_8500', 'DE.20100913.22358.fits.gz')
    hdul = fits.open(maskdef_file)
    det_par = spec.get_detector_par(hdul, det=det)

    # Load
    slitmask = spec.get_slitmask(maskdef_file)
    slits_file = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT', 'keck_deimos',
                              '830G_M_8500', 'Masters', 'MasterSlits_A_1_03.fits.gz')
    slits = slittrace.SlitTraceSet.from_file(slits_file)

    specobjs_file = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT', 'keck_deimos',
                                 '830G_M_8500', 'Science',
                                 'spec1d_DE.20100913.22358-CFHQS1_DEIMOS_2010Sep13T061231.334.fits')
    sobjs = specobjs.SpecObjs.from_fitsfile(specobjs_file)
    # Init at null
    for sobj in sobjs:
        sobj.MASKOBJ_NAME = None

    # Run me
    slitmask.assign_maskinfo(slits, sobjs, det_par['platescale'])

    # Test
    assert sobjs[0].MASKOBJ_NAME == 'ero89'

    # Write sobjs
    sobjs.write_to_fits({}, data_path('tst_sobjs.fits'))


@dev_suite_required
def test_deimosslitmask():
    f = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_deimos', '830G_M_8500',
                     'DE.20100913.22358.fits.gz')
    spec = KeckDEIMOSSpectrograph()
    spec.get_slitmask(f)
    assert spec.slitmask.nslits == 106, 'Incorrect number of slits read!'

