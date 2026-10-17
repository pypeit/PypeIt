"""
Module to test spectrograph read functions
"""
import copy
import os
import pathlib

import numpy as np
import pytest
import astropy.table

from pypeit import dataPaths
from pypeit import PypeItError
from pypeit import spectrographs
from pypeit.spectrographs.util import load_spectrograph
from pypeit import pypeitsetup
from pypeit.tests import tstutils
from pypeit.tests.tstutils import data_output_path

from IPython import embed


@pytest.mark.remote_data
def test_shanekastblue():
    s = spectrographs.shane_kast.ShaneKastBlueSpectrograph()
    example_file = dataPaths.tests.get_file_path('b1.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Shane Kast blue read.'
    det=1
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (350, 2112)
    assert bpm.shape == (2048,350)


@pytest.mark.remote_data
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
    assert spec.select_detectors(subset=[3,(1,2,3)]) == [3,(1,2,3)], 'Bad mix detector/mosaic selection'

    # String input that is *not* slitspatnum
    spec = load_spectrograph('keck_deimos') 
    assert spec.select_detectors(subset='3') == [3]
    assert spec.select_detectors(subset="3,(1,5)") == [3,(1,5)], 'Bad string of mix detector/mosaic selection'
    assert spec.select_detectors(subset="[3,(1,5)]") == [3,(1,5)], 'Bad string of mix detector/mosaic selection'

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

@pytest.mark.remote_data
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
    assert spec5.allowed_extensions == [".fits.fz"], 'Found wrong extensions'

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

@pytest.mark.remote_data
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

def test_apf_levy_final_config_frametypes():
    """
    Test the final_config_frametypes method for APF Levy spectrograph.
    
    The method should change 'pixelflat,trace' frames to 'pixelflat' when:
    - There are 'trace' frames (narrowflat) in the table
    - AND the setup decker is '3.0'
    """
    # Load the spectrograph
    spec = load_spectrograph('apf_levy')

    # Test case 1: decker is '3.0' with both narrowflat and wideflat frames
    # Should change wideflat to pixelflat
    table1 = astropy.table.Table()
    table1['frametype'] = ['pixelflat,trace', 'trace', 'science', 'arc']
    table1['filename'] = ['file1.fits', 'file2.fits', 'file3.fits', 'file4.fits']
    setup1 = {'decker': '3.0', 'binning': '1,1'}

    spec.final_config_frametypes(setup1, table1)

    # Check that wideflat frames were changed to pixelflat
    assert table1['frametype'][0] == 'pixelflat', \
        "Wideflat frame should be changed to pixelflat when decker is 3.0 and narrowflat exists"

    # Test case 2: decker is '8.0'
    # Should NOT change wideflat frames even if narrowflat exists
    table2 = astropy.table.Table()
    table2['frametype'] = ['pixelflat,trace', 'trace', 'science']
    table2['filename'] = ['file1.fits', 'file2.fits', 'file3.fits']
    setup2 = {'decker': '8.0', 'binning': '1,1'}

    spec.final_config_frametypes(setup2, table2)

    # Check that wideflat frames were NOT changed
    assert table2['frametype'][0] == 'pixelflat,trace', \
        "Wideflat frame should NOT be changed when decker is not 3.0"
    assert table2['frametype'][1] == 'trace', \
        "Narrowflat frame should remain unchanged"


def test_ifu_get_datacube_bins_wavelength_axis_first():
    """`get_datacube_bins()` must return its three bin-edge arrays in
    `(spec_bins, ybins, xbins)` order -- i.e. wavelength first -- for every
    `SlicerIFU` spectrograph that implements it.

    `pypeit.core.datacube.subpixellate()` always treats `bins[0]` as the
    wavelength axis (see `outshape = (bins[0].size-1, ...)`), regardless of
    what `get_datacube_bins()` actually returns; the per-pixel voxel
    coordinates it bins against are independently, always ordered
    `(wave_pix, dec_pix, ra_pix)`. `keck_kcwi.py`'s `get_datacube_bins()` was
    updated to this order when `datacube.py`'s cube-axis convention was
    reworked, but `gemini_gnirs.py`'s `GNIRSIFUSpectrograph` and
    `gtc_osiris.py`'s `GTCMAATSpectrograph` were not (a TODO to that effect
    was left in `keck_kcwi.py` but never actioned).
    `num_wave` and `slitlength` are chosen distinct from every spectrograph's
    (fixed) number of IFU slices, so a bin-order swap is guaranteed to be
    caught by a bin-count mismatch rather than accidentally matching.
    """
    num_wave = 1234
    slitlength = 50
    minmax = np.array([[-25.0, 25.0]])

    for spec_name, num_slices in [
        ('keck_kcwi', 24), ('keck_kcrm', 24), ('gemini_gnirs_ifu', 21), ('gtc_maat', 23)
    ]:
        spectrograph = load_spectrograph(spec_name)
        spec_bins, ybins, xbins = spectrograph.get_datacube_bins(slitlength, minmax, num_wave)
        assert spec_bins.size - 1 == num_wave, \
            f"{spec_name}: get_datacube_bins()[0] should be the {num_wave} wavelength bins, " \
            f"not the {spec_bins.size - 1} slice bins"
        assert ybins.size - 1 == slitlength, \
            f"{spec_name}: get_datacube_bins()[1] should be the {slitlength} along-slit bins"
        assert xbins.size - 1 == num_slices, \
            f"{spec_name}: get_datacube_bins()[2] should be the {num_slices} slice bins, " \
            f"not the {xbins.size - 1} wavelength bins"


def test_ldt_rimas(tmp_path, monkeypatch):
    """Exercise RIMAS frame typing, mode validation, and dither-sidecar behavior."""
    rimas = load_spectrograph('ldt_rimas_vph')

    # Flat exposure ranges across the supported configurations.
    accepted = [
        ('Vph300', "1.2'' long", 'HK', 5.72),
        ('Vph300', "1.2'' long", 'YJ', 60.06),
        ('Vph300', "0.6''", 'HK', 30.0),
        ('Vph300', "0.6''", 'YJ', 60.06),
        ('Vph300', "1.0''", 'HK', 15.0),
        ('Vph300', "1.0''", 'YJ', 60.06),
        ('Vph300', "2.0''", 'HK', 9.0),
        ('Vph300', "2.0''", 'YJ', 30.0),
        ('Vph30', "1.2'' long", 'HK', 5.72),
        ('Vph30', "1.2'' long", 'YJ', 28.6),
        ('Vph30', "0.6''", 'HK', 3.0),
        ('Vph30', "0.6''", 'YJ', 9.0),
        ('Vph30', "1.0''", 'YJ', 9.0),
        ('Vph30', "2.0''", 'YJ', 6.0),
        ('Grism', "0.6''", 'HK', 60.06),
        ('Grism', "0.6''", 'YJ', 120.12),
        ('Grism', "1.0''", 'HK', 30.0),
        ('Grism', "1.0''", 'YJ', 120.12),
        ('Grism', "2.0''", 'HK', 15.0),
        ('Grism', "2.0''", 'YJ', 120.12),
    ]
    rejected = [
        ('Vph30', "1.0''", 'HK', 16.0),
        ('Vph30', "2.0''", 'HK', 11.0),
        ('Vph300', "1.2'' long", 'HK', 60.06),
    ]
    fitstbl = astropy.table.Table(
        rows=accepted + rejected, names=('dispname', 'decker', 'arm', 'exptime')
    )
    fitstbl['idname'] = ['DOME_FLAT'] * len(fitstbl)
    fitstbl['filter1'] = ['open'] * len(fitstbl)
    expected = np.array([True] * len(accepted) + [False] * len(rejected))
    assert np.array_equal(rimas.check_frame_type('pixelflat', fitstbl), expected)
    fitstbl['idname'] = ['DOME_BACKGROUND'] * len(fitstbl)
    assert np.array_equal(rimas.check_frame_type('lampoffflats', fitstbl), expected)

    # Darks are retained only when their arm and raw geometry match the mode.
    fitstbl = astropy.table.Table(
        rows=[
            ('vph30.fits', 'Vph30', "1.2'' long", 'Vph30', 'HK', '4096,4230'),
            ('vph300.fits', 'Vph300', "1.2'' long", 'Vph300', 'YJ', '500,2000'),
            ('grism.fits', 'Grism', "1.0''", 'Grism', 'HK', '2000,2000'),
            ('dark_vph30.fits', 'Grism', "1.0''", 'blank', 'HK', '4096,4230'),
            ('dark_grism.fits', 'Vph30', "1.2'' long", 'blank', 'HK', '2000,2000'),
            ('dark_wrong_arm.fits', 'Grism', "1.0''", 'blank', 'YJ', '4096,4230'),
        ],
        names=('filename', 'dispname', 'decker', 'filter1', 'arm', 'rawshape'),
    )
    vph = rimas.validate_fitstbl(fitstbl)
    assert set(vph['filename']) == {'vph30.fits', 'vph300.fits', 'dark_vph30.fits'}
    grism = load_spectrograph('ldt_rimas_grism').validate_fitstbl(fitstbl)
    assert set(grism['filename']) == {'grism.fits', 'dark_grism.fits'}

    sidecar_name = 'rimas_dither_corrections.ecsv'
    sidecar_columns = ('filename', 'dithpat', 'dithnum', 'dithpos', 'dithseq', 'source')

    def write_sidecar(directory, rows):
        sidecar = directory / sidecar_name
        astropy.table.Table(rows=rows, names=sidecar_columns).write(
            sidecar, format='ascii.ecsv', overwrite=True
        )
        return sidecar

    # Prefer the setup directory, then fall back to a raw-data directory.
    setup_dir, raw_dir = tmp_path / 'setup', tmp_path / 'raw'
    setup_dir.mkdir()
    raw_dir.mkdir()
    raw_sidecar = write_sidecar(
        raw_dir, [('science.fits', 'ABBA', 2, 'B', 20, 'raw log')]
    )
    setup_sidecar = write_sidecar(
        setup_dir, [('science.fits', 'ABBA', 1, 'A', 10, 'setup log')]
    )
    discovery_table = astropy.table.Table(
        rows=[('science.fits', str(raw_dir))], names=('filename', 'directory')
    )
    monkeypatch.chdir(setup_dir)
    assert rimas._load_dither_corrections(discovery_table) == {
        'science.fits': ('ABBA', 'A', 10)
    }
    setup_sidecar.unlink()
    assert raw_sidecar.is_file()
    assert rimas._load_dither_corrections(discovery_table) == {
        'science.fits': ('ABBA', 'B', 20)
    }

    # Reject files that do not have the schema produced by the observer tool.
    monkeypatch.chdir(tmp_path)
    astropy.table.Table(
        rows=[('science.fits', 1, 10)], names=('filename', 'dithnum', 'dithseq')
    ).write(tmp_path / sidecar_name, format='ascii.ecsv')
    with pytest.raises(PypeItError, match='must contain exactly these columns'):
        rimas._load_dither_corrections(
            astropy.table.Table(rows=[('science.fits',)], names=('filename',))
        )

    def science_table(nframe, pattern, positions, mjd=None):
        filenames = [f'science_{i}.fits' for i in range(nframe)]
        table = astropy.table.Table()
        table['filename'] = filenames
        table['directory'] = [str(raw_dir)] * nframe
        table['frametype'] = ['science'] * nframe
        table['setup'] = ['A'] * nframe
        table['target'] = ['target'] * nframe
        table['dithpat'] = np.array([pattern] * nframe, dtype='U8')
        table['dithpos'] = np.array(positions, dtype='U8')
        table['mjd'] = np.arange(nframe, dtype=float) if mjd is None else mjd
        table['comb_id'] = np.arange(nframe, dtype=int)
        table['bkg_id'] = np.full(nframe, -1, dtype=int)
        return filenames, table

    # Correct two ABBA sequences while keeping unlisted header-based frames separate.
    filenames, fitstbl = science_table(
        6, 'ABBA', ['B', 'A', 'A', 'B', 'A', 'B'], mjd=[0, 10, 1, 2, 20, 21]
    )
    write_sidecar(
        tmp_path,
        [
            (filenames[0], 'ABBA', 1, 'A', 10, 'observing log'),
            (filenames[1], 'ABBA', 2, 'B', 10, 'observing log'),
            (filenames[2], 'ABBA', 2, 'B', 20, 'observing log'),
            (filenames[3], 'ABBA', 4, 'A', 20, 'observing log'),
        ],
    )
    result = rimas.get_comb_group(fitstbl)
    assert list(result['dithpos']) == ['A', 'B', 'B', 'A', 'A', 'B']
    assert list(result['comb_id']) == [0, 1, 2, 3, 4, 5]
    assert list(result['bkg_id']) == [1, 0, 3, 2, 5, 4]

    # A complete corrected ABBA sequence is reduced as AA-BB and BB-AA.
    filenames, fitstbl = science_table(4, 'ABBA', ['A', 'B', 'B', 'A'])
    write_sidecar(
        tmp_path,
        [
            (filenames[0], 'ABBA', 1, 'A', 1, 'observing log'),
            (filenames[1], 'ABBA', 2, 'B', 1, 'observing log'),
            (filenames[2], 'ABBA', 3, 'B', 1, 'observing log'),
            (filenames[3], 'ABBA', 4, 'A', 1, 'observing log'),
        ],
    )
    result = rimas.get_comb_group(fitstbl)
    assert list(result['comb_id']) == [0, 1, 1, 0]
    assert list(result['bkg_id']) == [1, 0, 0, 1]

    # The sidecar can replace an incorrect header pattern with On/Off metadata.
    filenames, fitstbl = science_table(2, 'NONE', ['Unknown'] * 2)
    write_sidecar(
        tmp_path,
        [
            (filenames[0], 'ONOFF', 1, 'On', 1, 'observing log'),
            (filenames[1], 'ONOFF', 2, 'Off', 1, 'observing log'),
        ],
    )
    result = rimas.get_comb_group(fitstbl)
    assert list(result['dithpat']) == ['OnOff', 'OnOff']
    assert list(result['dithpos']) == ['On', 'Off']
    assert list(result['comb_id']) == [0, 1]
    assert list(result['bkg_id']) == [1, 0]
