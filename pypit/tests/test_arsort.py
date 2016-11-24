# Module to run tests on arsort

import pytest


from pypit import pyputils
msgs = pyputils.get_dummy_logger()
from pypit import arsort
from pypit import arutils
from pypit import arparse as settings

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_match_science():
    """ Test match science routine
    """
    # Load settings
    arutils.dummy_settings(spectrograph='kast_blue')
    # Generate fitsdict
    fitsdict = pyputils.get_dummy_fitsdict()
    pytest.set_trace()
    #

