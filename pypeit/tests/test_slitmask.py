# Module to run tests on PypeItPar classes
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import numpy

import pytest

from pypeit.spectrographs.keck_deimos import KeckDEIMOSSpectrograph

# These tests are not run on Travis
if os.getenv('PYPEIT_DEV') is None:
    skip_test=True
else:
    skip_test=False

def test_deimosslitmask():
    if skip_test:
        return
    f = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Keck_DEIMOS', '830G_M_8500',
                     'DE.20100913.22358.fits.gz')
    spec = KeckDEIMOSSpectrograph()
    spec.get_slitmask(f)
    assert spec.slitmask.nslits == 106, 'Incorrect number of slits read!'

