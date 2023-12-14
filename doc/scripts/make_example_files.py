"""
Dynamically build example files included in the documentation.
"""

from importlib import resources
import os
import pathlib
import shutil 
import sys
import time

from pypeit.scripts import setup
from pypeit import pypeitsetup

from IPython import embed

PYP_ROOT = resources.files('pypeit').parent
DEV_ROOT = pathlib.Path(os.getenv('PYPEIT_DEV')).resolve()

#-----------------------------------------------------------------------------

def make_example_kast_pypeit_file(version, date):

    oroot = PYP_ROOT / 'doc' / 'include'
    droot = DEV_ROOT / 'RAW_DATA' / 'shane_kast_blue' / '600_4310_d55'
    
    pargs = setup.Setup.parse_args(['-r', str(droot), '-s', 'shane_kast_blue', '-c', 'all',
                                    '-d', str(oroot),
                                    '--version_override', version, 
                                    '--date_override', date])
    setup.Setup.main(pargs)

    ofile = oroot / 'shane_kast_blue_A.pypeit.rst'
    with open(ofile, 'w') as f:
        with open(oroot / 'shane_kast_blue_A' / 'shane_kast_blue_A.pypeit', 'r') as p:
            lines = p.readlines()
        f.write('.. code-block:: console\n')
        f.write('\n')
        for l in lines:
            f.write('    '+l)
        f.write('\n\n')

    shutil.rmtree(oroot / 'shane_kast_blue_A')


def make_example_deimos_pypeit_file(version, date):

    oroot = PYP_ROOT / 'doc' / 'include'
    droot = DEV_ROOT / 'RAW_DATA' / 'keck_deimos' / '1200G_M_7750'
    
    pargs = setup.Setup.parse_args(['-r', str(droot), '-s', 'keck_deimos', '-c', 'all',
                                    '-d', str(oroot),
                                    '--version_override', version, 
                                    '--date_override', date])
    setup.Setup.main(pargs)

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


def make_example_gnirs_pypeit_files(version, date):

    oroot = PYP_ROOT / 'doc' / 'include'

    # Create the default pypeit file
    droot = DEV_ROOT / 'RAW_DATA' / 'gemini_gnirs_echelle' / '32_SB_SXD'
    
    pargs = setup.Setup.parse_args(['-r', str(droot), '-s', 'gemini_gnirs_echelle', '-b',
                                    '-c', 'A', '-d', str(oroot),
                                    '--version_override', version, 
                                    '--date_override', date])
    setup.Setup.main(pargs)

    ofile = oroot / 'gemini_gnirs_echelle_A.pypeit.rst'
    with open(ofile, 'w') as f:
        with open(oroot / 'gemini_gnirs_echelle_A' / 'gemini_gnirs_echelle_A.pypeit', 'r') as p:
            lines = p.readlines()
        f.write('.. code-block:: console\n')
        f.write('\n')
        for l in lines:
            f.write('    '+l)
        f.write('\n\n')

    shutil.rmtree(oroot / 'gemini_gnirs_echelle_A')

    # Copy over the one that is actually used by the dev-suite
    dev = DEV_ROOT / 'pypeit_files' / 'gemini_gnirs_echelle_32_sb_sxd.pypeit'

    ofile = oroot / 'gemini_gnirs_echelle_A_corrected.pypeit.rst'
    with open(ofile, 'w') as f:
        with open(dev, 'r') as p:
            lines = p.readlines()
        f.write('.. code-block:: console\n')
        f.write('\n')
        for l in lines:
            f.write('    '+l)
        f.write('\n\n')


def make_example_nires_pypeit_files(version, date):

    oroot = PYP_ROOT / 'doc' / 'include'

    # Create the default pypeit file
    droot = DEV_ROOT / 'RAW_DATA' / 'keck_nires' / 'ABBA_wstandard'
    
    pargs = setup.Setup.parse_args(['-r', str(droot), '-s', 'keck_nires', '-b', '-c', 'A',
                                    '-d', str(oroot),
                                    '--version_override', version, 
                                    '--date_override', date])
    setup.Setup.main(pargs)

    ofile = oroot / 'keck_nires_A.pypeit.rst'
    with open(ofile, 'w') as f:
        with open(oroot / 'keck_nires_A' / 'keck_nires_A.pypeit', 'r') as p:
            lines = p.readlines()
        f.write('.. code-block:: console\n')
        f.write('\n')
        for l in lines:
            f.write('    '+l)
        f.write('\n\n')

    shutil.rmtree(oroot / 'keck_nires_A')

    # Copy over the one that is actually used by the dev-suite
    dev = DEV_ROOT / 'pypeit_files' / 'keck_nires_abba_wstandard.pypeit'

    ofile = oroot / 'keck_nires_A_corrected.pypeit.rst'
    with open(ofile, 'w') as f:
        with open(dev, 'r') as p:
            lines = p.readlines()
        f.write('.. code-block:: console\n')
        f.write('\n')
        for l in lines:
            f.write('    '+l)
        f.write('\n\n')


def make_example_sorted_file():

    root = DEV_ROOT / 'RAW_DATA' / 'keck_deimos'
    files = list(root.joinpath('830G_L_8100').glob('*fits*'))
    files += list(root.joinpath('830G_L_8400').glob('*fits*'))

    ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_deimos')
    ps.run(setup_only=True)

    # Write the sorted file,
    sorted_file = pathlib.Path().resolve() / ps.pypeit_file.replace('.pypeit', '.sorted')
    ps.fitstbl.write_sorted(sorted_file)

    oroot = PYP_ROOT / 'doc' / 'include'
    ofile = oroot / 'keck_deimos.sorted.rst'
    with open(ofile, 'w') as f:
        with open(sorted_file, 'r') as p:
            lines = p.readlines()
        f.write('.. code-block:: console\n')
        f.write('\n')
        for l in lines:
            f.write('    '+l)
        f.write('\n\n')

    os.remove(sorted_file)

def make_meta_examples():
    oroot = PYP_ROOT / 'doc' / 'include'
    ofile = oroot / 'deimos_meta_key_map.rst'
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
    make_example_kast_pypeit_file('1.12.2', '2023-04-05T22:42:29.971')
    print('Making keck_deimos_A.pypeit.rst')
    make_example_deimos_pypeit_file('1.12.2', '2023-04-05T22:42:29.971')
    print('Making gemini_gnirs files')
    make_example_gnirs_pypeit_files('1.12.2', '2023-04-05T22:42:29.971')
    print('Making keck_nires files')
    make_example_nires_pypeit_files('1.12.2', '2023-04-05T22:42:29.971')
    print('Making keck_deimos.sorted.rst')
    make_example_sorted_file()
    print('Make meta examples')
    make_meta_examples()
    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))


