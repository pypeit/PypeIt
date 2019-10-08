"""
Module to run tests on arsave
"""
import os

import numpy as np
import pytest

from astropy import units
from astropy.io import fits

from pypeit.core import save

from pypeit.tests.tstutils import dummy_fitstbl
from pypeit.spectrographs import util

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def mk_specobj(flux=5, objid=500):
    # specobj
    npix = 100
    specobj = specobj.SpecObj((100,100), 0, (0.4,0.6), objtype='science',
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
    fitstbl = dummy_fitstbl(directory=data_path(''))
    # Kludge
    fitstbl.table.remove_column('filename')
    fitstbl['filename'] = 'b1.fits.gz'
    # Settings
    #settings.argflag['run']['directory']['science'] = data_path('')
    spectrograph = util.load_spectrograph('shane_kast_blue')
    # Fill with dummy images
    dum = np.ones((100,100))
    sci_dict = {}
    sci_dict[0] = {}
    sci_dict[0]['sciframe'] = dum
    sci_dict[0]['finalvar'] = dum * 2
    sci_dict[0]['finalsky'] = dum + 0.1

    sci_dict['meta'] = {}
    sci_dict['meta']['vel_corr'] = 0.
    sci_dict['meta']['ir_redux'] = False

    basename = 'test'
    scidx = 5
    path = fitstbl['directory'][scidx]
    ifile = fitstbl['filename'][scidx]
    rawfile = os.path.join(path, ifile)
    master_dir = data_path('MF')+'_'+'shane_kast_blue'
    outfile = data_path('') + 'spec2d_{:s}.fits'.format(basename)
    # Create a dummy master_key_dict
    master_key_dict = dict(frame='', bpm='bpmkey',bias='',arc='',trace='',flat='')
    raw_hdr = fits.open(rawfile)[0].header
    save.save_2d_images(sci_dict, raw_hdr, spectrograph, master_key_dict, master_dir, outfile)
    # Read and test
    head0 = fits.getheader(data_path('spec2d_test.fits'))
    assert head0['PYPMFDIR'] == master_dir
    assert head0['BPMMKEY'] == 'bpm'            # See save_2d_images; removes last 3 characters
    assert 'PYPEIT' in head0['PIPELINE']


