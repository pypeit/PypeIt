# Module to run tests on ararclines
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import pytest

from pypit import pyputils
msgs = pyputils.get_dummy_logger()
from pypit import ararclines as alines
from pypit import arutils
#from pypit import arwave as arwv

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_load_linelist():
    """ Touches nearly all the methods in ararclines
    Returns
    -------

    """
    # Init
    arutils.dummy_settings()
    # Load
    alist = alines.load_arcline_list(None,None,['CuI','ArI','NeI'],None,
                                 modify_parse_dict=dict(NeI={'min_wave': 3000.},ArI={'min_intensity': 399.}))
    # Min NeI
    NeI = alist['Ion'] == 'NeI'
    np.testing.assert_allclose(np.min(alist['wave'][NeI]), 3455.1837999999998)

