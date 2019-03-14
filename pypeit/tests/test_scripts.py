"""
Module to run tests on scripts
"""
import os
import sys
import glob

import pytest

import matplotlib
matplotlib.use('agg')  # For Travis

from pypeit import msgs
from pypeit import scripts
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

def test_run_pypeit():

                    command_line = ['pypeit_setup', '-r', rawdir, '-s', instr.lower(), '-c all', '-o',
                                    '--output_path', wdir ]

                wdir = os.path.join(wdir, instr.lower()+'_A')
                print('Finished running pypeit on {0} --- '.format(outfile_root),
                      file=sys.stderr, end='')

                # Run pypeit on _A.pypeit
                pyp_file = glob.glob(os.path.join(wdir, '*_A.pypeit'))
                if len(pyp_file) != 1:
                    print('Could not find expected pypeit file: {0}'.format(pyp_file))
                    print("\x1B[" + "1;31m" + "FAILED" + "\x1B[" + "0m", file=sys.stderr)
                    passed = False
                else:
                    pyp_file = os.path.split(pyp_file[0])[1]
                    if retval == 0:
                        print("\x1B[" + "1;32m" + "PASSED" + "\x1B[" + "0m", file=sys.stderr)
                        npass += 1
                subprocess.call(['tail', '-2', logfile])
                print("\n", file=sys.stderr)

                ntest += 1
                logfile = os.path.join(wdir, outfile_root+'.test.log')
                print('Running pypeit on {:s} --- '.format(pyp_file), file=sys.stderr, end='')
                with open(logfile, 'w') as f:
                    print('Directory: {0}'.format(wdir))
                    command_line = [ 'run_pypeit', pyp_file, '-o' ]
                    print('Command line: {0}'.format(' '.join(command_line)))
                    retval = subprocess.call(command_line, stderr=f, cwd=wdir)
                    npass = report_test(retval, npass)
                    subprocess.call(['tail', '-2', logfile])
                    print("\n", file=sys.stderr)


@cooked_required
def test_show_1dspec():
    spec_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science',
                                'spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits')
    # Just list
    pargs = scripts.show_1dspec.parser([spec_file, '--list'])
    scripts.show_1dspec.main(pargs, unit_test=True)


@cooked_required
def test_chk_edges():
    mstrace_root = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Trace',
                                'MasterTrace_KeckLRISr_400_8500_det1.fits')
    # Ginga needs to be open in RC mode
    ginga.connect_to_ginga(raise_err=True, allow_new=True)
    #
    pargs = scripts.chk_edges.parser([mstrace_root])
    scripts.chk_edges.main(pargs)


def test_view_fits():
    """ Only test the list option
    """
    spec_file = data_path('spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits')
    pargs = scripts.view_fits.parser([spec_file, '--list'])


def test_coadd():
    coadd_file = data_path('coadd_UGC3672A_red.yaml')
    args = scripts.coadd_1dspec.parser([coadd_file])
    # Main
    gparam, ex_value, flux_value, iobj, outfile, files, local_kwargs \
            = scripts.coadd_1dspec.main(args, unit_test=True, path=data_path('./'))
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
    args = scripts.coadd_1dspec.parser([coadd_file])
    # Main
    gparam, ex_value, flux_value, iobj, outfile, files, obj_kwargs \
            = scripts.coadd_1dspec.main(args, unit_test=True, path=data_path('./'))
    # Test
    assert len(iobj) == len(files)
    # Crash it
    coadd_file = data_path('coadd_UGC3672A_red_badlist.yaml')
    args = scripts.coadd_1dspec.parser([coadd_file])
    with pytest.raises(IOError):
        gparam, ex_value, flux_value, iobj, outfile, files, _ \
                = scripts.coadd_1dspec.main(args, unit_test=True, path=data_path('./'))

