# Module to run tests on arsave
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import numpy as np
import pytest

from astropy import units
from astropy.io import fits

from pypeit import specobjs
from pypeit import metadata
from pypeit.core import save

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def mk_specobj(flux=5, objid=500):
    # specobj
    npix = 100
    specobj = specobjs.SpecObj((100,100), 0, (0.4,0.6), objtype='science',
                               spat_pixpos=300)
    specobj.boxcar = dict(wave=np.arange(npix)*units.AA, counts=np.ones(npix)*flux)
    specobj.optimal = dict(wave=np.arange(npix)*units.AA, counts=np.ones(npix)*flux-0.5)
    specobj.objid = objid
    specobj.trace_spat = np.arange(npix) / npix
    specobj.fwhmfit = np.arange(npix) / npix
    # Return
    return specobj


def test_save2d_fits():
    #settings.dummy_settings()
    #fitsdict = arutils.dummy_fitsdict(nfile=1, spectrograph='none', directory=data_path(''))
    fitstbl = metadata.dummy_fitstbl(directory=data_path(''))
    # Kludge
    fitstbl.table.remove_column('filename')
    fitstbl['filename'] = 'b1.fits.gz'
    # Settings
    #settings.argflag['run']['directory']['science'] = data_path('')
    setup = 'A_01_aa'
    spectrograph = 'shane_kast_blue'
    # Fill with dummy images
    dum = np.ones((100,100))
    sci_dict = {}
    sci_dict[0] = {}
    sci_dict[0]['sciframe'] = dum
    sci_dict[0]['finalvar'] = dum * 2
    sci_dict[0]['finalsky'] = dum + 0.1
    basename = 'test'
    scidx = 5
    save.save_2d_images(sci_dict, fitstbl, scidx, 0, setup, data_path('MF')+'_'+spectrograph,
                        data_path(''), basename)
    # Read and test
    head0 = fits.getheader(data_path('spec2d_test.fits'))
    assert head0['PYPCNFIG'] == 'A'
    assert head0['PYPCALIB'] == 'aa'
    assert 'PYPEIT' in head0['PIPELINE']


def test_save1d_fits():
    """ save1d to FITS and HDF5
    """
    # Init
    fitstbl = metadata.dummy_fitstbl(spectrograph='shane_kast_blue', directory=data_path(''))
    sobj = mk_specobj()
    specObjs = specobjs.SpecObjs([sobj])
    # Write to FITS
    save.save_1d_spectra_fits(specObjs, fitstbl[5], data_path('tst.fits'))


# NEEDS REFACTORING
#def test_save1d_hdf5():
#    """ save1d to FITS and HDF5
#    """
#    # Dummy self
#    fitstbl = arsort.dummy_fitstbl(spectrograph='shane_kast_blue', directory=data_path(''))
#    slf = arsciexp.dummy_self(fitstbl=fitstbl)
#    # specobj
#    slf._specobjs = []
#    slf._specobjs.append([])
#    slf._specobjs[0].append([mk_specobj(objid=455), mk_specobj(flux=3., objid=555)])
#    # Write to HDF5
#    arsave.save_1d_spectra_hdf5(slf, fitstbl)

