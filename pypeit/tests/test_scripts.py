"""
Module to run tests on scripts
"""
import os
import sys
import glob
import shutil

from configobj import ConfigObj

import pytest

import matplotlib
matplotlib.use('agg')  # For Travis

from pypeit.par.util import parse_pypeit_file, make_pypeit_file
from pypeit import msgs
from pypeit.pypeitsetup import PypeItSetup
from pypeit.scripts import setup, run_pypeit, show_1dspec, coadd_1dspec, chk_edges, view_fits
from pypeit.tests.tstutils import dev_suite_required, cooked_required
from pypeit import ginga

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


#def test_arcid_plot():
#    json_file = data_path('LRISb_600_WaveCalib_01.json')
#    pargs = arcid_plot.parser([json_file, 'LRISb', 'tmp.pdf'])
#    # Run
#    arcid_plot.main(pargs)


@dev_suite_required
def test_run_pypeit():
    # Get the directories
    rawdir = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA/Shane_Kast_blue/600_4310_d55/')
    assert os.path.isdir(rawdir), 'Incorrect raw directory'

    # Just get a few files
    testrawdir = os.path.join(rawdir, 'TEST')
    if os.path.isdir(testrawdir):
        shutil.rmtree(testrawdir)
    os.makedirs(testrawdir)
    files = [ 'b21.fits.gz', 'b22.fits.gz', 'b23.fits.gz', 'b27.fits.gz', 'b1.fits.gz',
              'b11.fits.gz', 'b12.fits.gz', 'b13.fits.gz' ]
    for f in files:
        shutil.copy(os.path.join(rawdir, f), os.path.join(testrawdir, f))

    outdir = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT_TEST')

    # For previously failed tests
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)

    # Run the setup
    args = setup.parser(['-r', testrawdir, '-s', 'shane_kast_blue', '-c all', '-o',
                                        '--output_path', outdir])
    setup.main(args)

    # Change to the configuration directory and set the pypeit file
    configdir = os.path.join(outdir, 'shane_kast_blue_A')
    pyp_file = os.path.join(configdir, 'shane_kast_blue_A.pypeit')
    assert os.path.isfile(pyp_file), 'PypeIt file not written.'

    # This data does not have enough traces to do the PCA, so the
    # prediction for the trace has to use the 'nearest' mode, instead of
    # the 'pca' mode.  So we need to modify the parameters used and
    # re-write the pypeit file
    ps = PypeItSetup.from_pypeit_file(pyp_file)
    ps.run()
    # Get rid of the existing calibrations keywords because they just
    # set the number
    cfg = ConfigObj(ps.user_cfg)
    del cfg['calibrations']
    cfg['calibrations'] = {}
    cfg['calibrations']['slitedges'] = {}
    cfg['calibrations']['slitedges']['sync_predict'] = 'nearest'
    ps.fitstbl.write_pypeit(os.path.split(pyp_file)[0], cfg_lines=cfg.write())

    # Perform the original reductions
    args = run_pypeit.parser([pyp_file, '-o'])
    run_pypeit.main(args)

    # Now try to reuse the old masters
    args = run_pypeit.parser([pyp_file, '-o', '-m'])
    run_pypeit.main(args)

    # Now try not overwriting and using the old masters
    args = run_pypeit.parser([pyp_file, '-m'])
    run_pypeit.main(args)

    # Clean-up
    shutil.rmtree(outdir)
    shutil.rmtree(testrawdir)


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


def test_view_fits():
    """ Only test the list option
    """
    spec_file = data_path('spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits')
    pargs = view_fits.parser([spec_file, '--list'])


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


