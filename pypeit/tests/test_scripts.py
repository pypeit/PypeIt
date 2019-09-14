"""
Module to run tests on scripts
"""
import os
import sys
import glob
import shutil

import pytest

import matplotlib
matplotlib.use('agg')  # For Travis

from pypeit import msgs
from pypeit.scripts import setup, show_1dspec, coadd_1dspec, chk_edges, view_fits, chk_flats
from pypeit.tests.tstutils import dev_suite_required, cooked_required
from pypeit import ginga

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


'''
def test_view_fits():
    """ Only test the list option
    """
    spec_file = data_path('spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits')
    pargs = view_fits.parser([spec_file, '--list'])
'''

@cooked_required
def test_show_1dspec():
    spec_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science',
                             'spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits')
    # Just list
    pargs = show_1dspec.parser([spec_file, '--list'])
    show_1dspec.main(pargs, unit_test=True)


@cooked_required
def test_chk_edges():
    mstrace_root = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Trace',
                                'MasterTrace_KeckLRISr_400_8500_det1.fits')
    # Ginga needs to be open in RC mode
    ginga.connect_to_ginga(raise_err=True, allow_new=True)
    #
    pargs = chk_edges.parser([mstrace_root])
    chk_edges.main(pargs)

@cooked_required
def test_chk_flat():
    mstrace_root = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Shane_Kast_blue',
                                'MasterFlat_A_1_01.fits')
    # Ginga needs to be open in RC mode
    ginga.connect_to_ginga(raise_err=True, allow_new=True)
    #
    pargs = chk_flats.parser([mstrace_root])
    chk_flats.main(pargs)

def test_coadd():
    coadd_file = data_path('coadd_UGC3672A_red.yaml')
    args = coadd_1dspec.parser([coadd_file])
    # Main
    gparam, ex_value, flux_value, iobj, outfile, files, local_kwargs \
            = coadd_1dspec.main(args, unit_test=True, path=data_path('./'))
    # Test
    assert len(gparam) == 2
    assert isinstance(gparam, dict)
    assert ex_value == 'opt'
    assert flux_value is True
    assert iobj == 'O210-S1467-D02-I0012'
    assert outfile == 'UGC3672A_r.fits'
    assert len(files) == 4
    assert isinstance(local_kwargs, dict)
    assert 'otol' in list(local_kwargs.keys())
    assert 'scale_method' in list(gparam.keys())



def test_coadd2():
    """ Test using a list of object names
    """
    coadd_file = data_path('coadd_UGC3672A_red_objlist.yaml')
    args = coadd_1dspec.parser([coadd_file])
    # Main
    gparam, ex_value, flux_value, iobj, outfile, files, obj_kwargs \
            = coadd_1dspec.main(args, unit_test=True, path=data_path('./'))
    # Test
    assert len(iobj) == len(files)
    # Crash it
    coadd_file = data_path('coadd_UGC3672A_red_badlist.yaml')
    args = coadd_1dspec.parser([coadd_file])
    with pytest.raises(IOError):
        gparam, ex_value, flux_value, iobj, outfile, files, _ \
                = coadd_1dspec.main(args, unit_test=True, path=data_path('./'))

