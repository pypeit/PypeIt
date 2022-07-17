import os
import shutil

from IPython import embed

import pytest

import numpy as np

#from pypeit.par.util import parse_pypeit_file
from pypeit.tests.tstutils import data_path
from pypeit.metadata import PypeItMetaData
from pypeit.spectrographs.util import load_spectrograph
from pypeit.scripts.setup import Setup
from pypeit.inputfiles import PypeItFile

def test_read_combid():

    # ------------------------------------------------------------------
    # In case of failed tests
    setup_dir = data_path('setup_files')
    if os.path.isdir(setup_dir):
        shutil.rmtree(setup_dir)
    config_dir = data_path('shane_kast_blue_A')
    if os.path.isdir(config_dir):
        shutil.rmtree(config_dir)
    # ------------------------------------------------------------------

    # Generate the pypeit file with the comb_id
    droot = data_path('b')
    pargs = Setup.parse_args(['-r', droot, '-s', 'shane_kast_blue', '-c=all', '-b',
                             '--extension=fits.gz', '--output_path={:s}'.format(data_path(''))])
    Setup.main(pargs)
    shutil.rmtree(setup_dir)

    pypeit_file = os.path.join(config_dir, 'shane_kast_blue_A.pypeit')
    pypeItFile = PypeItFile.from_file(pypeit_file)

    # Get the spectrograph
    spectrograph = None
    for l in pypeItFile.cfg_lines:
        if 'spectrograph' in l:
            spectrograph = load_spectrograph(l.split(' ')[-1])
            break
    assert spectrograph is not None, 'Did not appropriately read spectrograph'

    # Set the metadata
    pmd = PypeItMetaData(spectrograph, spectrograph.default_pypeit_par(), 
                         files=pypeItFile.filenames,
                         usrdata=pypeItFile.data, strict=False)

    indx = pmd['filename'] == 'b27.fits.gz'
    assert pmd['comb_id'][indx] == [1], 'Incorrect combination group ID'
    assert pmd['comb_id'][np.where(~indx)[0]][0] == -1, 'Incorrect combination group ID'

    shutil.rmtree(config_dir)
