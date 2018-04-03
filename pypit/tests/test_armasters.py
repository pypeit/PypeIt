# Module to run tests on armasters
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pytest

from pypit import pyputils
msgs = pyputils.get_dummy_logger()#develop=True)
from pypit import armasters

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


#@pytest.fixture
#def fitsdict():
#    return arutils.dummy_fitsdict()

def test_master_name():
    """ Test master name method
    """
    types = ['bias', 'badpix', 'trace', 'normpixelflat', 'arc', 'wave', 'wv_calib', 'tilts']
    suff = ['Bias', 'BadPix', 'Trace', 'FlatField', 'Arc', 'Wave', 'WaveCalib', 'Tilts']
    for isuff,itype in zip(suff,types):
        if itype == 'wv_calib':
            exten = 'json'
        else:
            exten = 'fits'
        assert armasters.master_name(itype, '01', mdir='MasterFrames') == 'MasterFrames/Master{:s}_01.{:s}'.format(isuff,exten)

