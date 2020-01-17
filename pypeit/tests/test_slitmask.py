"""
Module to run tests on PypeItPar classes
"""
import os
import numpy

import pytest

from pypeit.spectrographs.keck_deimos import KeckDEIMOSSpectrograph
from pypeit.tests.tstutils import dev_suite_required

@dev_suite_required
def test_deimosslitmask():
    f = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_deimos', '830G_M_8500',
                     'DE.20100913.22358.fits.gz')
    spec = KeckDEIMOSSpectrograph()
    spec.get_slitmask(f)
    assert spec.slitmask.nslits == 106, 'Incorrect number of slits read!'

