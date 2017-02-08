# Module to run tests on ararclines


import os
import numpy as np
import pytest

from linetools import utils as ltu

from pypit import pyputils
msgs = pyputils.get_dummy_logger()
from pypit import arparse as settings
from pypit import arvcorr as py_arvcorr
from pypit import arutils as arut


#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)

hdr_mjd = 57783.269661
hdr_exptime = 1200.
RA = '07:06:23.45'
DEC = '+30:20:50.5'
coord = ltu.radec_to_coord((RA,DEC))
hdr_ra = coord.ra.value
hdr_dec = coord.dec.value
hdr_equ = 2000.
hdr_lon = 155.47833            # Longitude of the telescope (NOTE: West should correspond to positive longitudes)
hdr_lat = 19.82833              # Latitude of the telescope
hdr_alt = 4160.0               # Elevation of the telescope (in m)

def test_jd_to_date():
    hdr_jd = hdr_mjd + 2400000.5
    year, month, day, ut = py_arvcorr.jd_to_date(hdr_jd)
    assert year == 2017
    assert month == 1
    assert day == 30
    assert np.isclose(ut, 6.471863999962807)

def test_vhelio():
    """ Test the full helio-centric calculation
    """
    hdr_jd = hdr_mjd + 2400000.5
    vhel = py_arvcorr.vhelio(hdr_jd, hdr_exptime, hdr_ra, hdr_dec, hdr_equ, hdr_lat, hdr_lon, hdr_alt)
    # IDL
    # vhel = x_keckhelio(106.59770833333332, 30.34736111111111, 2000., jd=2457783.769661)
    #    vrotate = -0.25490532
    pytest.set_trace()

'''
def test_vhelio():
    """ Test the full helio-centric calculation
    Returns
    -------

    """
    # Initialize some settings
    arut.dummy_settings()
    # Load Dummy self
    slf = arut.dummy_self()
    settings.argflag['run']['spectrograph'] = 'kast_blue'
    settings.spect['arc'] = {}
    settings.spect['arc']['index'] = [[0]]
    fitsdict = arut.dummy_fitsdict()
    # Run
    arcparm = pyarc.setup_param(slf, 0, 1, fitsdict)
    for key in ['llist','disp','wvmnx']:
        assert key in arcparm
'''
