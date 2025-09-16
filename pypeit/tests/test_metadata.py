from pathlib import Path
import shutil
import os

from IPython import embed

import pytest

import numpy as np

from pypeit.tests import tstutils
from pypeit.metadata import PypeItMetaData
from pypeit.spectrographs.util import load_spectrograph
from pypeit.scripts.setup import Setup
from pypeit.inputfiles import PypeItFile
from astropy.table import Table


def test_read_combid():

    # ------------------------------------------------------------------
    # In case of failed tests
    config_dir = Path(tstutils.data_output_path('shane_kast_blue_A')).absolute()
    if config_dir.exists():
        shutil.rmtree(config_dir)
    # ------------------------------------------------------------------

    tstutils.install_shane_kast_blue_raw_data()

    # Generate the pypeit file with the comb_id
    droot = tstutils.data_output_path('b')
    pargs = Setup.parse_args(['-r', droot, '-s', 'shane_kast_blue', '-c', 'all', '-b',
                              '--output_path', f'{config_dir.parent}'])
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
    no_combid_indx = np.logical_not(b27_indx | b24_indx)
    b27_indx = np.where(b27_indx)[0][0]
    b24_indx = np.where(b24_indx)[0][0]

    assert pmd['comb_id'][b27_indx] > 0, 'Science file should have a combination group ID'
    assert pmd['comb_id'][b24_indx] > 0, 'Standard file should have a combination group ID'
    assert pmd['comb_id'][b27_indx] != pmd['comb_id'][b24_indx], 'Science and standard should not have same combination group ID'
    assert pmd['comb_id'][np.where(no_combid_indx)[0][0]] == -1, 'Incorrect combination group ID'

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

def test_soar_goodman_metadata():
    # Load the spectrograph
    spectrograph = load_spectrograph("soar_goodman_red")
    # Setup a fake table with information about files
    names = ('ra','dec','target','decker','binning','exptime','mjd','airmass','dispname','mode','dispangle','idname','lampstat01','lampstat02','lampstat03','lampstat04','lampstat05','lampstat06','lampstat07','lampstat08','directory','filename','manual','shift')
    dtype = (np.float64, np.float64, str, str, str, np.float64, np.float64, np.float64, str, str, np.float64, str, str, str, str, str, str, str, str, str, str, str, str, np.int32)
    fitstbl = Table(names=names, dtype=dtype)
    # Add four dummy rows for a standard, science, dome flat, and arc lamp frame
    # Standard
    fitstbl.add_row((298.6849,0.2755,'hr7596','1.0_LONG_SLIT','2,2',0.4,60930.9714,1.32,'400_SYZY','400_M1',5.7992,'SPECTRUM','FALSE','FALSE','FALSE','FALSE','FALSE','FALSE','FALSE','FALSE','/home','hr7596.fits.fz','',0))
    # Science
    fitstbl.add_row((5.0208,8.5102,'supernova','1.0_LONG_SLIT','2,2',400.0,60931.1597,1.46,'400_SYZY','400_M1',5.7999,'SPECTRUM','FALSE','FALSE','FALSE','FALSE','FALSE','FALSE','FALSE','FALSE','/home','supernova.fits.fz','',0))
    # Arc lamp
    fitstbl.add_row((5.0212,8.5116,'supernova_comp','1.0_LONG_SLIT','2,2',0.5,60931.1544,1.48,'400_SYZY','400_M1',5.7995,'ARC','TRUE','TRUE','FALSE','FALSE','FALSE','FALSE','FALSE','FALSE','/home','supernova.fits.fz','',0))
    # Flat
    fitstbl.add_row((148.3446,-17.6491,'dflat','1.0_LONG_SLIT','2,2',7.0,60930.7888,1.66,'400_SYZY','400_M1',5.7997,'LAMPFLAT','FALSE','FALSE','FALSE','FALSE','FALSE','FALSE','FALSE','TRUE','/home','dflat.fits.fz','',0))

    # Check standard
    tst = spectrograph.check_frame_type('standard',fitstbl)
    assert np.array_equal(tst, np.array([True, False, False, False]))
    # Check science - note that the standard will also pass this test
    tst = spectrograph.check_frame_type('science',fitstbl)
    assert np.array_equal(tst, np.array([False, True, False, False]))
    # Check arc
    tst = spectrograph.check_frame_type('arc',fitstbl)
    assert np.array_equal(tst, np.array([False, False, True, False]))
    # Check flat
    tst = spectrograph.check_frame_type('pixelflat',fitstbl)
    assert np.array_equal(tst, np.array([False, False, False, True]))


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
    filelist = tstutils.make_fake_fits_files()
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
