from pathlib import Path
import shutil
import os

from IPython import embed

import pytest

import numpy as np

#from pypeit.par.util import parse_pypeit_file
from pypeit.tests.tstutils import data_path, make_fake_fits_files
from pypeit.metadata import PypeItMetaData
from pypeit.spectrographs.util import load_spectrograph
from pypeit.scripts.setup import Setup
from pypeit.inputfiles import PypeItFile
from astropy.table import Table


def test_read_combid():

    # ------------------------------------------------------------------
    # In case of failed tests
    config_dir = Path(data_path('shane_kast_blue_A')).resolve()
    if config_dir.exists():
        shutil.rmtree(config_dir)
    # ------------------------------------------------------------------

    # Generate the pypeit file with the comb_id
    droot = data_path('b')
    pargs = Setup.parse_args(['-r', droot, '-s', 'shane_kast_blue', '-c', 'all', '-b',
                             '--extension', 'fits.gz', '--output_path', f'{config_dir.parent}'])
    Setup.main(pargs)

    pypeit_file = config_dir / 'shane_kast_blue_A.pypeit'
    pypeItFile = PypeItFile.from_file(str(pypeit_file))

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

    b27_indx = pmd['filename'] == 'b27.fits.gz'
    b24_indx = pmd['filename'] == 'b24.fits.gz'
    assert pmd['comb_id'][b27_indx] > 0, 'Science file should have a combination group ID'
    assert pmd['comb_id'][b24_indx] > 0, 'Standard file should have a combination group ID'
    assert pmd['comb_id'][b27_indx] != pmd['comb_id'][b24_indx], 'Science and standard should not have same combination group ID'
    no_combid_indx = np.logical_not(b27_indx | b24_indx)
    assert pmd['comb_id'][np.where(no_combid_indx)[0]][0] == -1, 'Incorrect combination group ID'

    shutil.rmtree(config_dir)


def test_nirspec_lamps():
    # Load the spectrograph
    spectrograph = load_spectrograph("keck_nirspec_low")
    # Setup a fake table with information about files
    fitstbl = Table(names=('fakename', 'lampstat01', 'lampstat02', 'lampstat03', 'lampstat04', 'lampstat05', 'lampstat06'), dtype=('S', 'd', 'd', 'd', 'd', 'd', 'd'))
    fitstbl.add_row(('off_01', 0, 0, 0, 0, 0, 0))
    fitstbl.add_row(('off_02', 0, 0, 0, 0, 0, 0))
    fitstbl.add_row(('arcs_01', 0, 0, 0, 0, 1, 0))
    fitstbl.add_row(('arcs_02', 0, 0, 1, 0, 0, 0))
    fitstbl.add_row(('arcs_03', 1, 1, 1, 1, 1, 0))
    fitstbl.add_row(('dome_01', 0, 0, 0, 0, 0, 1))
    fitstbl.add_row(('dome_02', 0, 0, 0, 0, 0, 1))
    fitstbl.add_row(('dome_03', 0, 0, 0, 0, 0, 1))
    # Check off
    tst = spectrograph.lamps(fitstbl, 'off')
    assert np.array_equal(tst, np.array([True, True, False, False, False,  False,  False,  False]))
    # Check arcs
    tst = spectrograph.lamps(fitstbl, 'arcs')
    assert np.array_equal(tst, np.array([False, False, True, True, True, False, False, False]))
    # Check dome
    tst = spectrograph.lamps(fitstbl, 'dome')
    assert np.array_equal(tst, np.array([False, False, False, False, False,  True,  True,  True]))


def test_setup_iter():

    gen = PypeItMetaData.configuration_generator()
    assert next(gen) == 'A', 'First setup identifier changed'

    end = False
    while not end:
        try:
            setup = next(gen)
        except StopIteration:
            end = True

    assert setup == 'ZZ', 'Last setup identifier changed'
    assert len(list(PypeItMetaData.configuration_generator())) \
                == PypeItMetaData.maximum_number_of_configurations(), \
                'Number of configuration identifiers changed'


def test_multiple_setups():
    filelist = make_fake_fits_files()
    spectrograph = load_spectrograph("shane_kast_blue")
    # Set the metadata
    fitstbl = PypeItMetaData(spectrograph, spectrograph.default_pypeit_par(), files=filelist, strict=True)
    fitstbl.get_frame_types()
    cfgs = fitstbl.unique_configurations()
    fitstbl.set_configurations(configs=cfgs)
    # Now do some checks
    for ff in range(len(fitstbl)):
        if fitstbl[ff]['frametype'] == 'bias':
            assert(len(fitstbl[ff]['setup'].split(",")) == 2)  # Two configurations for the bias frames
        else:
            assert (len(fitstbl[ff]['setup'].split(",")) == 1)  # One configuration for everything else
    # Remove the created files
    for fil in filelist:
        os.remove(fil)
