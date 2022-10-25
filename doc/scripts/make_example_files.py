"""
Dynamically build example files included in the documentation.
"""

from pathlib import Path
import sys
import os
import time
import shutil 
import glob

from pkg_resources import resource_filename

from IPython import embed

from pypeit.tests.tstutils import data_path
from pypeit.scripts import setup
from pypeit import pypeitsetup

#-----------------------------------------------------------------------------

def make_example_kast_pypeit_file():

    oroot = os.path.join(os.path.split(os.path.abspath(resource_filename('pypeit', '')))[0],
                                       'doc', 'include')

    droot = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA/shane_kast_blue/600_4310_d55')
    droot += '/'
    pargs = setup.Setup.parse_args(['-r', droot, '-s', 'shane_kast_blue', '-c', 'all', '-d', oroot])
    setup.Setup.main(pargs)

    shutil.rmtree(os.path.join(oroot, 'setup_files'))

    ofile = os.path.join(oroot, 'shane_kast_blue_A.pypeit.rst')
    with open(ofile, 'w') as f:
        with open(os.path.join(oroot, 'shane_kast_blue_A', 'shane_kast_blue_A.pypeit'), 'r') as p:
            lines = p.readlines()
        f.write('.. code-block:: console\n')
        f.write('\n')
        for l in lines:
            f.write('    '+l)
        f.write('\n\n')

    shutil.rmtree(os.path.join(oroot, 'shane_kast_blue_A'))


def make_example_deimos_pypeit_file():

    oroot = Path(resource_filename('pypeit', '')).resolve().parent / 'doc' / 'include'
    if (oroot / 'setup_files').exists():
        shutil.rmtree(oroot / 'setup_files')

    droot = Path(os.getenv('PYPEIT_DEV')).resolve() / 'RAW_DATA' / 'keck_deimos' / '1200G_M_7750'
    
    pargs = setup.Setup.parse_args(['-r', str(droot), '-s', 'keck_deimos', '-c', 'all',
                                    '-d', str(oroot)])
    setup.Setup.main(pargs)

    shutil.rmtree(oroot / 'setup_files')

    ofile = oroot / 'keck_deimos_A.pypeit.rst'
    with open(ofile, 'w') as f:
        with open(oroot / 'keck_deimos_A' / 'keck_deimos_A.pypeit', 'r') as p:
            lines = p.readlines()
        f.write('.. code-block:: console\n')
        f.write('\n')
        for l in lines:
            f.write('    '+l)
        f.write('\n\n')

    shutil.rmtree(oroot / 'keck_deimos_A')


def make_example_gnirs_pypeit_files():

    oroot = Path(resource_filename('pypeit', '')).resolve().parent / 'doc' / 'include'
    if (oroot / 'setup_files').exists():
        shutil.rmtree(oroot / 'setup_files')

    # Create the default pypeit file
    droot = Path(os.getenv('PYPEIT_DEV')).resolve() / 'RAW_DATA' / 'gemini_gnirs' / '32_SB_SXD'
    
    pargs = setup.Setup.parse_args(['-r', str(droot), '-s', 'gemini_gnirs', '-b', '-c', 'A',
                                    '-d', str(oroot)])
    setup.Setup.main(pargs)

    shutil.rmtree(oroot / 'setup_files')

    ofile = oroot / 'gemini_gnirs_A.pypeit.rst'
    with open(ofile, 'w') as f:
        with open(oroot / 'gemini_gnirs_A' / 'gemini_gnirs_A.pypeit', 'r') as p:
            lines = p.readlines()
        f.write('.. code-block:: console\n')
        f.write('\n')
        for l in lines:
            f.write('    '+l)
        f.write('\n\n')

    shutil.rmtree(oroot / 'gemini_gnirs_A')

    # Copy over the one that is actually used by the dev-suite
    dev = Path(os.getenv('PYPEIT_DEV')).resolve() \
                / 'pypeit_files' / 'gemini_gnirs_32_sb_sxd.pypeit'

    ofile = oroot / 'gemini_gnirs_A_corrected.pypeit.rst'
    with open(ofile, 'w') as f:
        with open(dev, 'r') as p:
            lines = p.readlines()
        f.write('.. code-block:: console\n')
        f.write('\n')
        for l in lines:
            f.write('    '+l)
        f.write('\n\n')


def make_example_sorted_file():

    root = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_deimos')
    files = glob.glob(os.path.join(root, '830G_L_8100', '*fits*'))
    files += glob.glob(os.path.join(root, '830G_L_8400', '*fits*'))

    ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_deimos')
    ps.run(setup_only=True)

    sfile = os.path.abspath('keck_deimos.sorted')
    oroot = os.path.join(os.path.split(os.path.abspath(resource_filename('pypeit', '')))[0],
                                       'doc', 'include')
    ofile = os.path.join(oroot, 'keck_deimos.sorted.rst')
    with open(ofile, 'w') as f:
        with open(sfile, 'r') as p:
            lines = p.readlines()
        f.write('.. code-block:: console\n')
        f.write('\n')
        for l in lines:
            f.write('    '+l)
        f.write('\n\n')

    os.remove(sfile)

def make_meta_examples():

    ofile = Path(resource_filename('pypeit', '')).resolve().parent \
                / 'doc' / 'include' / 'deimos_meta_key_map.rst'
    otmp = ofile.parent / 'tmp_meta'
    if otmp.exists():
        otmp.unlink()

    stdout = sys.__stdout__
    with open(otmp, 'w') as sys.stdout:
        from pypeit.spectrographs.util import load_spectrograph
        spec = load_spectrograph('keck_deimos')
        spec.meta_key_map()
    sys.stdout = stdout

    with open(ofile, 'w') as f:
        with open(otmp, 'r') as p:
            lines = p.readlines()
        f.write('.. code-block:: console\n')
        f.write('\n')
        for l in lines:
            f.write('    '+l)
        f.write('\n\n')

    if otmp.exists():
        otmp.unlink()

if __name__ == '__main__':
    t = time.perf_counter()
    print('Making shane_kast_blue_A.pypeit.rst')
    make_example_kast_pypeit_file()
    print('Making keck_deimos_A.pypeit.rst')
    make_example_deimos_pypeit_file()
    print('Making gemini_gnirs_A.pypeit.rst')
    make_example_gnirs_pypeit_files()
    print('Making keck_deimos.sorted.rst')
    make_example_sorted_file()
    print('Make meta examples')
    make_meta_examples()
    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))


