"""
Dynamically build example files included in the documentation.
"""

from pathlib import Path
import sys
import os
import time
import shutil 
import glob

from pypeit import dataPaths
from pypeit.scripts import setup
from pypeit import pypeitsetup
from pypeit.inputfiles import PypeItFile
from pypeit.cache import git_most_recent_tag

from IPython import embed

from pypeit.tests.tstutils import data_path
from pypeit.scripts import setup
from pypeit import pypeitsetup

#-----------------------------------------------------------------------------

def obscure_path(paths):
    """
    Removes the local root path to "PypeIt-development-suite/" and replaces it
    with "/path/to".
    """
    return [str(pathlib.Path('/path/to') / pathlib.Path(p).absolute().relative_to(DEV_ROOT.parent))
                for p in paths]


def verbatim_to_rst(inp, out):
    with open(out, 'w') as f:
        with open(inp, 'r') as p:
            lines = p.readlines()
        f.write('.. code-block:: console\n')
        f.write('\n')
        for l in lines:
            f.write('    '+l)
        f.write('\n\n')


def extinction_files():
    inp = dataPaths.extinction.get_file_path('extinction_curves.txt')
    ofile = PYP_ROOT / 'doc' / 'include' / f'{inp.name}.rst'
    verbatim_to_rst(inp, ofile)


def make_example_kast_pypeit_file(version, date):

    # Set the paths for the rst file, ...
    oroot = PYP_ROOT / 'doc' / 'include'
    # ... with the data, and ...
    droot = DEV_ROOT / 'RAW_DATA' / 'shane_kast_blue' / '600_4310_d55'
    # ... for the setup files.
    sroot = oroot / 'shane_kast_blue_A'
    if sroot.is_dir():
        # Remove pre-existing directory
        shutil.rmtree(sroot)

    # Run setup 
    pargs = setup.Setup.parse_args(['-r', str(droot), '-s', 'shane_kast_blue', '-c', 'all',
                                    '-d', str(oroot),
                                    '--version_override', version, 
                                    '--date_override', date])
    setup.Setup.main(pargs)

    # Read the PypeItFile
    pypeit_file = sroot / 'shane_kast_blue_A.pypeit'
    p = PypeItFile.from_file(pypeit_file)
    # Abstract the directory so that it isn't system/developer dependent
    p.file_paths = obscure_path(p.file_paths)
    # Overwrite the file
    p.write(pypeit_file, version_override=version, date_override=date)

    ofile = oroot / 'shane_kast_blue_A.pypeit.rst'
    verbatim_to_rst(pypeit_file, ofile)

    shutil.rmtree(sroot)


def make_example_deimos_pypeit_file(version, date):

    oroot = PYP_ROOT / 'doc' / 'include'
    droot = DEV_ROOT / 'RAW_DATA' / 'keck_deimos' / '1200G_M_7750'
    sroot = oroot / 'keck_deimos_A'
    if sroot.is_dir():
        # Remove pre-existing directory
        shutil.rmtree(sroot)

    pargs = setup.Setup.parse_args(['-r', str(droot), '-s', 'keck_deimos', '-c', 'all',
                                    '-d', str(oroot),
                                    '--version_override', version, 
                                    '--date_override', date])
    setup.Setup.main(pargs)

    # Read the PypeItFile
    pypeit_file = sroot / 'keck_deimos_A.pypeit'
    p = PypeItFile.from_file(pypeit_file)
    # Abstract the directory so that it isn't system/developer dependent
    p.file_paths = obscure_path(p.file_paths)
    # Overwrite the file
    p.write(pypeit_file, version_override=version, date_override=date)

    ofile = oroot / 'keck_deimos_A.pypeit.rst'
    verbatim_to_rst(pypeit_file, ofile)

    shutil.rmtree(sroot)


def make_example_gnirs_pypeit_files(version, date):

    oroot = PYP_ROOT / 'doc' / 'include'
    droot = DEV_ROOT / 'RAW_DATA' / 'gemini_gnirs_echelle' / '32_SB_SXD'
    sroot = oroot / 'gemini_gnirs_echelle_A'
    if sroot.is_dir():
        # Remove pre-existing directory
        shutil.rmtree(sroot)
    
    pargs = setup.Setup.parse_args(['-r', str(droot), '-s', 'gemini_gnirs_echelle', '-b',
                                    '-c', 'A', '-d', str(oroot),
                                    '--version_override', version, 
                                    '--date_override', date])
    setup.Setup.main(pargs)

    # Read the PypeItFile
    pypeit_file = sroot / 'gemini_gnirs_echelle_A.pypeit'
    p = PypeItFile.from_file(pypeit_file)
    # Abstract the directory so that it isn't system/developer dependent
    p.file_paths = obscure_path(p.file_paths)
    # Overwrite the file
    p.write(pypeit_file, version_override=version, date_override=date)

    ofile = oroot / 'gemini_gnirs_echelle_A.pypeit.rst'
    verbatim_to_rst(pypeit_file, ofile)

    shutil.rmtree(sroot)

    # Copy over the one that is actually used by the dev-suite
    dev = Path(os.getenv('PYPEIT_DEV')).resolve() \
                / 'pypeit_files' / 'gemini_gnirs_echelle_32_sb_sxd.pypeit'

    ofile = oroot / 'gemini_gnirs_echelle_A_corrected.pypeit.rst'
    verbatim_to_rst(dev, ofile)


def make_example_nires_pypeit_files(version, date):

    oroot = PYP_ROOT / 'doc' / 'include'
    droot = DEV_ROOT / 'RAW_DATA' / 'keck_nires' / 'ABBA_wstandard'
    sroot = oroot / 'keck_nires_A'
    if sroot.is_dir():
        # Remove pre-existing directory
        shutil.rmtree(sroot)

    pargs = setup.Setup.parse_args(['-r', str(droot), '-s', 'keck_nires', '-b', '-c', 'A',
                                    '-d', str(oroot),
                                    '--version_override', version, 
                                    '--date_override', date])
    setup.Setup.main(pargs)

    # Read the PypeItFile
    pypeit_file = sroot / 'keck_nires_A.pypeit'
    p = PypeItFile.from_file(pypeit_file)
    # Abstract the directory so that it isn't system/developer dependent
    p.file_paths = obscure_path(p.file_paths)
    # Overwrite the file
    p.write(pypeit_file, version_override=version, date_override=date)

    ofile = oroot / 'keck_nires_A.pypeit.rst'
    verbatim_to_rst(pypeit_file, ofile)

    shutil.rmtree(sroot)

    # Copy over the one that is actually used by the dev-suite
    dev = Path(os.getenv('PYPEIT_DEV')).resolve() \
                / 'pypeit_files' / 'keck_nires_abba_wstandard.pypeit'

    ofile = oroot / 'keck_nires_A_corrected.pypeit.rst'
    verbatim_to_rst(dev, ofile)


def make_example_sorted_file():

    root = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_deimos')
    files = glob.glob(os.path.join(root, '830G_L_8100', '*fits*'))
    files += glob.glob(os.path.join(root, '830G_L_8400', '*fits*'))

    ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_deimos')
    ps.run(setup_only=True)

    # Write the sorted file,
    sorted_file = pathlib.Path().absolute() / ps.pypeit_file.replace('.pypeit', '.sorted')
    ps.fitstbl.write_sorted(sorted_file)

    oroot = Path(resource_filename('pypeit', '')).resolve().parent / 'doc' / 'include'
    ofile = oroot / 'keck_deimos.sorted.rst'
    verbatim_to_rst(sorted_file, ofile)

    os.remove(sorted_file)

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

    verbatim_to_rst(otmp, ofile)

    if otmp.exists():
        otmp.unlink()

if __name__ == '__main__':
    t = time.perf_counter()
    tag, date = git_most_recent_tag()
    print('Making shane_kast_blue_A.pypeit.rst')
    make_example_kast_pypeit_file(tag, date)
    print('Making keck_deimos_A.pypeit.rst')
    make_example_deimos_pypeit_file(tag, date)
    print('Making gemini_gnirs files')
    make_example_gnirs_pypeit_files(tag, date)
    print('Making keck_nires files')
    make_example_nires_pypeit_files(tag, date)
    print('Making keck_deimos.sorted.rst')
    make_example_sorted_file()
    print('Make meta examples')
    make_meta_examples()
    print('Make extinction file')
    extinction_files()
    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))


