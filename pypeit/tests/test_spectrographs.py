"""
Module to test spectrograph read functions
"""
import copy
import os
import pathlib

import pytest

from pypeit import dataPaths
from pypeit.pypmsgs import PypeItError
from pypeit import spectrographs
from pypeit.spectrographs.util import load_spectrograph
from pypeit import pypeitsetup
from pypeit.tests import tstutils
from pypeit.tests.tstutils import data_output_path

from IPython import embed


def test_shanekastblue():
    s = spectrographs.shane_kast.ShaneKastBlueSpectrograph()
    example_file = dataPaths.tests.get_file_path('b1.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Shane Kast blue read.'
    det=1
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (350, 2112)
    assert bpm.shape == (2048,350)


def test_select_detectors_pypeit_file():
    # Generate a PypeIt file
    tstutils.install_shane_kast_blue_raw_data()
    pypeItFile = tstutils.make_shane_kast_blue_pypeitfile()
    pypeit_file = data_output_path('test.pypeit')
    pypeItFile.write(pypeit_file)

    # Perform the setup
    setup = pypeitsetup.PypeItSetup.from_pypeit_file(pypeit_file)
    par, spectrograph, fitstbl = setup.run()

    assert spectrograph.select_detectors(subset=par['rdx']['detnum']) == [1], \
            'Incorrect detectors selected.'

    # Clean-up
    os.remove(pypeit_file)


def test_select_detectors_mosaic():

    spec = load_spectrograph('gemini_gmos_north_ham')

    # Invalid detector
    with pytest.raises(PypeItError):
        spec.select_detectors(subset=4)
    # Invalid mosaic
    with pytest.raises(PypeItError):
        spec.select_detectors(subset=(2,3))

    # Valid
    assert spec.select_detectors() == [1,2,3], 'Bad detector selection'
    # Valid
    assert spec.select_detectors(subset=[3,(1,2,3)]) == [3,(1,2,3)], 'Bad detector selection'


def test_list_detectors_deimos():
    deimos = load_spectrograph('keck_deimos')
    dets = deimos.list_detectors()
    assert dets.ndim == 2, 'DEIMOS has a 2D array of detectors'
    assert dets.size == 8, 'DEIMOS has 8 detectors'
    mosaics = deimos.list_detectors(mosaic=True)
    assert mosaics.ndim == 1, 'Mosaics are listed as 1D arrays'
    assert mosaics.size == 4, 'DEIMOS has 4 predefined mosaics'


def test_list_detectors_mosfire():
    mosfire = load_spectrograph('keck_mosfire')
    dets = mosfire.list_detectors()
    assert dets.ndim == 1, 'MOSFIRE has a 1D array of detectors'
    assert dets.size == 1, 'MOSFIRE has 1 detector'
    with pytest.raises(PypeItError):
        mosaics = mosfire.list_detectors(mosaic=True)


def test_list_detectors_mods():
    mods = load_spectrograph('lbt_mods1r')
    dets = mods.list_detectors()
    assert dets.ndim == 1, 'MODS1R has a 1D array of detectors'
    assert dets.size == 1, 'MODS1R has 1 detector'
    with pytest.raises(PypeItError):
        mosaics = mods.list_detectors(mosaic=True)


def test_list_detectors_hires():
    hires = load_spectrograph('keck_hires')
    dets = hires.list_detectors()
    assert dets.ndim == 1, 'HIRES has a 1D array of detectors'
    assert dets.size == 3, 'HIRES has 3 detectors'
    mosaics = hires.list_detectors(mosaic=True)
    assert mosaics.ndim == 1, 'Mosaics are listed as 1D arrays'
    assert mosaics.size == 1, 'HIRES has 1 predefined mosaic'


def test_configs():

    spec = load_spectrograph('keck_deimos')
    cfg1 = dict(amp='"SINGLE:B"',
                binning='1,1',
                decker='LongMirr',
                dispangle=8099.98291016,
                dispname='830G',
                filter1='OG550')
    cfg2 = dict(amp='"SINGLE:B"',
                binning='1,1',
                decker='LongMirr',
                dispangle=8399.93554688,
                dispname='830G',
                filter1='OG550')

    assert spec.same_configuration([cfg1,cfg1]), 'Configurations should be the same'
    assert not spec.same_configuration([cfg1,cfg2]), 'Configurations should be different'

    cfg3 = copy.deepcopy(cfg1)
    cfg3['dispangle'] *= (1.+spec.meta['dispangle']['rtol']/2)

    assert spec.same_configuration([cfg1,cfg3]), \
        'Configurations should be the same within tolerance'

    cfg3 = copy.deepcopy(cfg1)
    cfg3['dispangle'] *= (1.+2*spec.meta['dispangle']['rtol'])

    assert not spec.same_configuration([cfg1,cfg3]), \
        'Configurations should not be the same within tolerance'


def test_atmext():
    
    spec = load_spectrograph('keck_deimos')
    atmext = spec.get_atmospheric_extinction('closest')
    assert atmext.file == 'mkoextinct.dat', 'Found wrong extinction file'

    spec = load_spectrograph('shane_kast_blue')
    atmext = spec.get_atmospheric_extinction('closest')
    assert atmext.file == 'mthamextinct.dat', 'Found wrong extinction file'

    # Override the file
    atmext = spec.get_atmospheric_extinction('mkoextinct.dat')
    assert atmext.file == 'mkoextinct.dat', 'Used wrong extinction file'

def test_load_spectrograph():

    # Basic test
    spec = load_spectrograph('shane_kast_blue')
    assert isinstance(spec, spectrographs.spectrograph.Spectrograph), 'Not a Spectrograph class'

    # Load using existing class
    spec2 = load_spectrograph(spec)
    assert isinstance(spec2, spectrographs.spectrograph.Spectrograph), 'Not a Spectrograph class'

    # Load using a single processed data file
    raw_file = dataPaths.tests.get_file_path('spec1d_b28.fits', to_pkg='symlink')
    spec3 = load_spectrograph(raw_file)
    assert isinstance(spec3, spectrographs.spectrograph.Spectrograph), 'Not a Spectrograph class'
    assert spec3.name == 'shane_kast_blue'
    assert spec3.allowed_extensions == ['.fits', '.fits.gz'], 'Found wrong extensions'

    # None in --> None out
    spec4 = load_spectrograph(None)
    assert spec4 is None

    # Test the allowed extensions for an oddball spectrograph
    spec5 = load_spectrograph('soar_goodman_red')
    assert spec5.allowed_extensions == [".fz"], 'Found wrong extensions'

    # Call as it from a post-processing script
    spec6 = load_spectrograph('soar_goodman_red', pypeit_fits=True)
    assert spec6.allowed_extensions == [".fits"], 'Postproc scripts only allow .fits'

    # Call using instance and from a post-processing script
    spec7 = load_spectrograph(spec5, pypeit_fits=True)
    assert spec7.allowed_extensions == [".fits"], 'Postproc scripts only allow .fits'

    # Call using a single processed data file, and from a post-processing script
    spec8 = load_spectrograph(raw_file, pypeit_fits=True)
    assert spec8.allowed_extensions == [".fits"], 'Postproc scripts only allow .fits'


@pytest.fixture
def fitstbl():

    # Get the files
    file_names = [
        'b1.fits.gz',    # arc
        'b11.fits.gz',   # trace
        'b21.fits.gz',   # bias
        'b24.fits.gz',   # standard
        'b27.fits.gz'    # science
    ]
    files = [dataPaths.tests.get_file_path(f, to_pkg='symlink') for f in file_names]

    setupc = pypeitsetup.PypeItSetup(files, spectrograph_name='shane_kast_blue')
    setupc.build_fitstbl(files)
    setupc.fitstbl.finalize_usr_build(None, 'A')
    return setupc.fitstbl

def test_config_specific_par(fitstbl):
    # Grab a science file for configuration specific parameters
    indx = fitstbl.find_frames('science', index=True)[0]
    sci_file = fitstbl.frame_paths(indx)
 
    # Load the parameters based on the (raw) science file
    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.config_specific_par(sci_file)

    # Load the parameters based on the fitstbl object
    _ = spectrograph.config_specific_par(fitstbl.get_row_for_filename(sci_file))
    
    # Check the value of configuration-dependent `reid_arxiv`
    assert par['calibrations']['wavelengths']['reid_arxiv'] == 'shane_kast_blue_600.fits'

    # Change the ``dispname`` value in the fitstbl, and make sure the par changed
    ft2 = fitstbl.get_row_for_filename(sci_file)
    ft2['dispname'] = '452/3306'
    par = spectrograph.config_specific_par(ft2)
    assert par['calibrations']['wavelengths']['reid_arxiv'] == 'shane_kast_blue_452.fits'
