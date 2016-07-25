# Module to run tests on arsave


import numpy as np
import os
import pytest

#from astropy.io import fits
from astropy import units as u

from pypit import pyputils
msgs = pyputils.get_dummy_logger()
from pypit import arspecobj as pyp_sobj
from pypit import arsave as arsv

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)

def mk_specobj():
    # specobj
    npix = 100
    specobj = pyp_sobj.SpecObjExp((100,100), 'Kast', 0, 0, (0.4,0.6), 0.5, 0.5, objtype='science')
    specobj.boxcar = dict(wave=np.arange(npix)*u.AA, counts=np.ones(npix)*5)
    specobj.optimal = dict(wave=np.arange(npix)*u.AA, counts=np.ones(npix)*4.4)
    specobj.trace = np.arange(npix) / npix
    # Return
    return specobj

def test_save1d():
    """ save1d to FITS and HDF5
    """
    from pypit import arutils as arut
    # Dummy self
    slf = arut.dummy_self()
    slf._specobjs = []
    slf._specobjs.append([mk_specobj()])
    # Write to FITS
    arsv.save_1d_spectra_fits(slf)

"""
def test_save1d():
    """ save1d to FITS and HDF5
    """
    from pypit import arutils as arut
    # Dummy self
    slf = arut.dummy_self()
    # specobj
    npix = 100
    specobj = pyp_sobj.SpecObjExp((100,100), 'Kast', 0, 0, (0.4,0.6), 0.5, 0.5, objtype='science')
    specobj.boxcar = dict(wave=np.arange(npix)*u.AA, counts=np.ones(npix)*5)
    specobj.optimal = dict(wave=np.arange(npix)*u.AA, counts=np.ones(npix)*4.4)
    specobj.trace = np.arange(npix) / npix
    slf._specobjs = []
    slf._specobjs.append([specobj])
    # Write to HDF5
    arsv.save_1d_spectra_hdf5(slf)

"""
