"""
Dynamically build example files included in the documentation.
"""

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

def make_example_pypeit_file():

    oroot = os.path.join(os.path.split(os.path.abspath(resource_filename('pypeit', '')))[0],
                                       'doc', 'include')

    droot = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA/shane_kast_blue/600_4310_d55')
    droot += '/'
    pargs = setup.parse_args(['-r', droot, '-s', 'shane_kast_blue', '-c', 'all', '-d', oroot])
    setup.main(pargs)

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


if __name__ == '__main__':
    t = time.perf_counter()
    print('Making shane_kast_blue_A.pypeit.rst')
    make_example_pypeit_file()
    print('Making keck_deimos.sorted.rst')
    make_example_sorted_file()
    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))


