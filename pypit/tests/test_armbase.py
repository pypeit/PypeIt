# Module to run tests on armbase module
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# TEST_UNICODE_LITERALS

import numpy as np
import pytest

from pypit import arparse
from pypit import arsciexp
from pypit import armbase


def test_update_masters():
    pass
'''    
    # Dummy self
    arparse.dummy_settings(spectrograph='shane_kast_blue', set_idx=True)
    slf1 = arsciexp.dummy_self()
    slf1._idx_arcs = np.array([0,1])
    slf2 = arsciexp.dummy_self()
    sciexp = [slf1, slf2]
    #  Not actually filling anything
    armbase.UpdateMasters(sciexp, 0, 1, 'arc')

'''
