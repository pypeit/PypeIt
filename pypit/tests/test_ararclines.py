# Module to run tests on ararclines


import numpy as np
import pytest

from pypit import pyputils
msgs = pyputils.get_dummy_logger()
from pypit import ararclines as alines
#from pypit import arwave as arwv

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_load_linelist():
    """ Touches nearly all the methods in ararclines
    Returns
    -------

    """
    # Load
    alist = alines.load_arcline_list(None,None,['CuI','ArI','NeI'],None,
                                 modify_parse_dict=dict(NeI={'min_wave': 3000.},ArI={'min_intensity': 399.}))
    # Min NeI
    np.testing.assert_allclose(np.min(alist['wave'][alist['Ion'] == 'NeI']),
                               3455.1837999999998)

