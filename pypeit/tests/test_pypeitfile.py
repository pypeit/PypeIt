""" Tests on PypeItFile """

import pytest

from pypeit import pypeitfile
from pypeit.tests.tstutils import data_path

import configobj

#def test_read_pypeit_file():
    # Read the PypeIt file
    #cfg, data, frametype, usrdata, setups, _ \
pypeItFile = pypeitfile.PypeItFile.from_file(
                data_path('example_pypeit_file.pypeit'))

'''
tmp = {'rdx': {'spectrograph': 'keck_hires'}}
confObj = configobj.ConfigObj(tmp)

tmp2 = ['[rdx]', 'spectrograph = shane_kast_blue']
confObj2 = configobj.ConfigObj(tmp2)
'''

pytest.set_trace()