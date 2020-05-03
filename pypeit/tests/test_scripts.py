"""
Module to run tests on scripts
"""
import os
import shutil

import matplotlib
matplotlib.use('agg')  # For Travis

from astropy.io import fits

from pypeit.scripts import setup, show_1dspec, coadd_1dspec, chk_edges, view_fits, chk_flats
from pypeit.scripts import trace_edges, run_pypeit, ql_mos, show_2dspec
from pypeit.tests.tstutils import dev_suite_required, cooked_required
from pypeit import edgetrace
from pypeit import ginga


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@dev_suite_required
def test_quicklook():
    # The following needs the LRISb calibration files to be
    # found in a CALIBS/ folder in the DEV Suite.  You can get
    # that folder here:
    # https://drive.google.com/drive/folders/1NSg5Rmr8hD_1-fOchQc3WXjt59D6f9od?usp=sharing
    calib_dir = os.path.join(os.environ['PYPEIT_DEV'], 'CALIBS')
    if not os.path.isdir(calib_dir):
        raise IOError("You need to get the CALIBS folder as described above!!")

    # Define the output directories (HARDCODED!!)
    outdir = os.path.join(os.getcwd(), 'keck_lris_blue_A')
    # Remove them if they already exist
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)

    # Raw path
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_lris_blue',
                         'long_600_4000_d560')
    ql_mos.main(ql_mos.parser(['keck_lris_blue', droot, 'b150910_2033.fits.gz',
                               'b150910_2051.fits.gz', 'b150910_2070.fits.gz', '--det=2',
                               '--user_pixflat={0}'.format(
                                   os.path.join(calib_dir,
                                        'PYPEIT_LRISb_pixflat_B600_2x2_17sep2009.fits.gz'))]))

@dev_suite_required
def test_trace_edges():
    # Define the output directories (HARDCODED!!)
    setupdir = os.path.join(os.getcwd(), 'setup_files')
    outdir = os.path.join(os.getcwd(), 'shane_kast_blue_A')
    masterdir = os.path.join(os.getcwd(), 'shane_kast_blue_A', 'Masters')
    # Remove them if they already exist
    if os.path.isdir(setupdir):
        shutil.rmtree(setupdir)
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)

    # Perform the setup
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/shane_kast_blue/600_4310_d55')
    droot += '/'
    setup.main(setup.parser(['-r', droot, '-s', 'shane_kast_blue', '-c', 'all']))

    # Generate the Masters folder
    os.mkdir(masterdir)

    # Define the pypeit file (HARDCODED!!)
    pypeit_file = os.path.join(outdir, 'shane_kast_blue_A.pypeit')

    # Run the tracing
    trace_edges.main(trace_edges.parser(['-f', pypeit_file]))

    # Define the edges master file (HARDCODED!!)
    trace_file = os.path.join(outdir, 'Masters', 'MasterEdges_A_1_01.fits.gz')

    # Check that the correct number of traces were found
    edges = edgetrace.EdgeTraceSet.from_file(trace_file)
    assert edges.ntrace == 2, 'Did not find the expected number of traces.'

    # Clean up
    shutil.rmtree(setupdir)
    shutil.rmtree(outdir)


@cooked_required
def test_show_1dspec():
    spec_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science',
                             'spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits')
    # Just list
    pargs = show_1dspec.parser([spec_file, '--list'])
    show_1dspec.main(pargs)

@cooked_required
def test_show_2dspec():
    droot = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked')
    spec2d_file = os.path.join(droot, 'Science',
                             'spec2d_b27-J1217p3905_KASTb_2015May20T045733.560.fits')
    # Ginga needs to be open in RC mode
    ginga.connect_to_ginga(raise_err=True, allow_new=True)
    # Save
    cdir = os.getcwd()
    os.chdir(droot)
    # List
    pargs = show_2dspec.parser([spec2d_file, '--list'])
    show_2dspec.main(pargs)
    # Show
    pargs = show_2dspec.parser([spec2d_file])
    show_2dspec.main(pargs)
    # Go back
    os.chdir(cdir)

@cooked_required
def test_chk_edges():
    mstrace_root = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Trace',
                                'MasterEdges_KeckLRISr_400_8500_det1.fits.gz')
    # Ginga needs to be open in RC mode
    ginga.connect_to_ginga(raise_err=True, allow_new=True)
    #
    pargs = chk_edges.parser([mstrace_root])
    chk_edges.main(pargs)


def test_view_fits():
    """ Only test the list option
    """
    spec_file = data_path('spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits')
    pargs = view_fits.parser([spec_file, '--list', 'shane_kast_blue'])


@cooked_required
def test_chk_flat():
    mstrace_root = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'shane_kast_blue',
                                'MasterFlat_A_1_01.fits')
    # Ginga needs to be open in RC mode
    ginga.connect_to_ginga(raise_err=True, allow_new=True)
    #
    pargs = chk_flats.parser([mstrace_root])
    chk_flats.main(pargs)


def test_coadd1d_1():
    """
    Test basic coadd using shane_kast_blue
    """
    # NOTE: flux_value is False
    parfile = 'coadd1d.par'
    if os.path.isfile(parfile):
        os.remove(parfile)
    coadd_ofile = data_path('J1217p3905_coadd.fits')
    if os.path.isfile(coadd_ofile):
        os.remove(coadd_ofile)

    coadd_ifile = data_path('shane_kast_blue.coadd1d')
    coadd_1dspec.main(coadd_1dspec.parser([coadd_ifile, '--test_spec_path', data_path('')]))

    hdu = fits.open(coadd_ofile)
    assert hdu[1].header['EXT_MODE'] == 'OPT'
    assert hdu[1].header['FLUXED'] is False

    # Clean up
    os.remove(parfile)
    os.remove(coadd_ofile)


def test_coadd1d_2():
    """
    Test combining Echelle
    """
    # NOTE: flux_value is False
    parfile = 'coadd1d.par'
    if os.path.isfile(parfile):
        os.remove(parfile)
    coadd_ofile = data_path('pisco_coadd.fits')
    if os.path.isfile(coadd_ofile):
        os.remove(coadd_ofile)

    coadd_ifile = data_path('gemini_gnirs_32_sb_sxd.coadd1d')
    coadd_1dspec.main(coadd_1dspec.parser([coadd_ifile, '--test_spec_path', data_path('')]))

    hdu = fits.open(coadd_ofile)
    assert hdu[1].header['EXT_MODE'] == 'OPT'
    assert hdu[1].header['FLUXED'] is False

    # Clean up
    os.remove(parfile)
    os.remove(coadd_ofile)

# TODO: Include tests for coadd2d, sensfunc, flux_calib
