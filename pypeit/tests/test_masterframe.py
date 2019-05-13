"""
Module to run tests on armasters
"""
import os
import numpy as np
import pytest

from pypeit import msgs
from pypeit import masterframe
from pypiet import io

#@pytest.fixture
#def fitsdict():
#    return arutils.dummy_fitsdict()

#def test_master_name():
#    """ Test master name method
#    """
#    types = ['bias', 'badpix', 'trace', 'pixelflat', 'arc', 'wave', 'wv_calib', 'tilts']
#    suff = ['Bias', 'BadPix', 'Trace', 'PixelFlat', 'Arc', 'Wave', 'WaveCalib', 'Tilts']
#    for isuff,itype in zip(suff,types):
#        if itype == 'wv_calib':
#            exten = '.json'
##        elif itype == 'trace':
##            exten = ''
#        else:
#            exten = '.fits'
#        assert masterframe.master_name(itype, '01', mdir='MasterFrames') == 'MasterFrames/Master{:s}_01{:s}'.format(isuff,exten)

def data_root():
    return os.path.join(os.path.dirname(__file__), 'files')

# Could break this into a bunch of separate tests...
def test_master_io():
    """
    Test the save and load methods by saving a fake master and reading it back in.
    """
    # Init
    master_dir = data_root()
    mf = masterframe.MasterFrame('Test', master_dir=master_dir, master_key='A_01_1',
                                 reuse_masters=True)
    # In case of previous test failure
    if os.path.isfile(mf.file_path):
        os.remove(mf.file_path)
    # Check the filename
    assert mf.file_name == 'MasterTest_A_01_1.fits', 'Incorrect master frame file name'
    # No file so load should return None
    assert mf.load('JUNK') is None, 'Load for non-existent frame should return None'
    # Should also return None for the header
    assert mf.load('JUNK', return_header=True) == (None, None), \
            'Load for non-existent frame should return None'
    # Save a fake file
    data = np.arange(10)
    extnames = 'JUNK'
    raw_files = ['fake1.fits', 'fake2.fits']
    steps = ['foo', 'bar']
    mf.save(data, extnames, overwrite=False, raw_files=raw_files, steps=steps)
    assert os.path.isfile(mf.file_path), 'No file written'
    # Try to load it
    _data = mf.load('JUNK')
    assert np.array_equal(data, _data), 'Data written/read incorrectly'
    # Try to load it with the header
    _data, _hdr = mf.load('JUNK', return_header=True)
    assert np.array_equal(data, _data), 'Data written/read incorrectly'
    assert _hdr['FRAMETYP'] == 'Test', 'Incorrect master type'
    assert _hdr['STEPS'].split(',') == steps, 'Steps written incorrectly'
    assert io.parse_hdr_key_group(_hdr, prefix='F') == raw_files, 'Did not correctly read files'
    # Loading if reuse_masters is false should yield None
    mf.reuse_masters = False
    assert mf.load('JUNK') is None, 'Load when reuse_masters=False should return None'
    # Clean up
    os.remove(mf.file_path)

