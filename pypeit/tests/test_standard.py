
from IPython import embed
import pytest

import numpy as np

from astropy import units
from astropy import coordinates

from pypeit.core import standard
from pypeit import PypeItError

def test_mab_to_cgs():
    wave = 5000.
    mAB = 14.
    cgs = standard.mAB_to_cgs(wave, mAB) * 1e17
    assert isinstance(cgs, (float, np.floating)), 'Should return scalar'
    assert np.absolute(cgs - 1094.41) < 0.01, 'Bad value'

    wave = [5000., 5500.]
    _cgs = standard.mAB_to_cgs(wave, mAB) * 1e17
    assert isinstance(_cgs, np.ndarray), 'Should return array'
    assert np.absolute(cgs - _cgs[0]) < 1e-10, 'Bad value'

    wave = 5000.
    mAB = [14., 14.5]
    _cgs = standard.mAB_to_cgs(wave, mAB) * 1e17
    assert isinstance(_cgs, np.ndarray), 'Should return array'
    assert np.absolute(cgs - _cgs[0]) < 1e-10, 'Bad value'

    wave = [5000., 5500.]
    mAB = [14., 14.5]
    _cgs = standard.mAB_to_cgs(wave, mAB) * 1e17
    assert isinstance(_cgs, np.ndarray), 'Should return array'
    assert np.absolute(cgs - _cgs[0]) < 1e-10, 'Bad value'

    # Broadcast error
    wave = [5000., 5500., 6000.]
    mAB = [14., 14.5]
    with pytest.raises(ValueError):
        _cgs = standard.mAB_to_cgs(wave, mAB) * 1e17


def test_archive_entry():
    archives = {
        'blackbody': 'BB113337-110529',
        'calspec': 'HD106252',
        'esofil': 'HR4468',
        'ing': 'SP0946+139',
        'noao': 'EG81',
        'xshooter': 'GD153',
    }
    for archive, name in archives.items():
        row = standard.archive_entry(archive, name)
        assert row['Name'] == name, 'Name mismatch'

    # Test bogus archive
    with pytest.raises(PypeItError):
        row = standard.archive_entry('bogus', 'junk')

    # Test inability to find a match
    with pytest.raises(PypeItError):
        row = standard.archive_entry('xshooter', 'junk')


def test_nearest_standard():
    # Test known values
    archives = {
        'blackbody': 'BB113337-110529',
        'calspec': 'HD106252',
        'esofil': 'HR4468',
        'ing': 'SP0946+139',
        'noao': 'EG81',
        'xshooter': 'GD153',
    }
    for archive, nearest in archives.items():
        sep, row = standard.nearest_archive_entry(archive, 180.0, 0.0)
        assert row['Name'] == nearest, f'Bad match.  Expected {nearest}, found {row["Name"]}.'
        assert sep > 1 * units.deg, 'None of the "nearest" matches are close.'

    # Test coordinate types
    sep, row = standard.nearest_archive_entry('blackbody', 180.0, 0.0)

    _sep, _row = standard.nearest_archive_entry('blackbody', '12:00:00.0', '00:00:00.0')
    assert np.absolute(sep.value - _sep.value) < 1e-10, 'Separations should be identical'
    assert row['Name'] == _row['Name'], 'Should match the same object'

    _sep, _row = standard.nearest_archive_entry('blackbody', '180.', '0.', unit=units.deg)
    assert np.absolute(sep.value - _sep.value) < 1e-10, 'Separations should be identical'
    assert row['Name'] == _row['Name'], 'Should match the same object'

    # Can't handle arrays
    with pytest.raises(PypeItError):
        sep, row = standard.nearest_archive_entry('blackbody', [120., 180.], [0., 0.])

    # Test bogus archive
    with pytest.raises(PypeItError):
        sep, row = standard.nearest_archive_entry('bogus', 180.0, 0.0)


def test_archive_classes():
    classes = standard.archived_flux_classes()
    assert len(classes) == 5, 'Number of classes changed'


def test_archived_standards():

    archive_classes = standard.archived_flux_classes()
    archive_names = {
        'calspec': {
            'name': 'HD106252',
            'len' : 6815,
            'wave0': 1711.28,
            'med_flux': 13174.91,
            'tol': 700.,
        },
        'esofil': {
            'name': 'HR4468',
            'len' : 426,
            'wave0': 3300.0,
            'med_flux': 2555500.0,
            'tol': 700.,
        },
        'ing': {
            'name': 'SP0946+139',
            'len' : 68,
            'wave0': 3080.0,
            'med_flux': 162979.32,
            'tol': 2150.,
        },
        'noao': {
            'name': 'EG81',
            'len' : 96,
            'wave0': 3200.0,
            'med_flux': 1685.17,
            'tol': 950.,
        },
        'xshooter': {
            'name': 'GD153',
            'len' : 174678,
            'wave0': 3000.92,
            'med_flux': 52.23,
            'tol': 1600.,
        },
    }
    for archive, meta in archive_names.items():
        sep, _name, file = archive_classes[archive].nearest_standard(180., 0.)
        assert _name == meta['name'], 'Found a different object.'
        file_path = archive_classes[archive].path.get_file_path(file)
        assert file_path.is_file(), 'File is not available.'
        spec = archive_classes[archive](file)

        _sep, row = standard.nearest_archive_entry(archive, 180., 0.)
        assert sep.value - _sep.value < 1e-10, 'Should match to the same object'
        assert row['Name'] == _name, 'Should match to the same object'

        # Should not match within tolerance
        with pytest.raises(PypeItError):
            spec = archive_classes[archive].from_coordinates(180., 0.)

        # Increase the tolerance
        spec = archive_classes[archive].from_coordinates(180., 0., tol=meta['tol'])
        assert np.absolute(spec.wave[0] - meta['wave0']) < 0.01, 'Bad wave'
        assert np.absolute(np.median(spec.flux) - meta['med_flux']) < 0.01, 'Bad flux'
        assert spec.meta['source'] == archive, 'Bad source'

        # Test name match
        spec = archive_classes[archive].from_name(meta['name'])
        assert spec.meta['Name'] == meta['name'], 'Name mismatch'
        assert np.absolute(spec.wave[0] - meta['wave0']) < 0.01, 'Bad wave'
        assert np.absolute(np.median(spec.flux) - meta['med_flux']) < 0.01, 'Bad flux'
        assert spec.meta['source'] == archive, 'Bad source'

        # Improve the coordinates
        coo = coordinates.SkyCoord(
            row['RA_2000'], row['DEC_2000'], unit=(units.hourangle, units.deg)
        )
        spec = archive_classes[archive].from_coordinates(coo.ra.value, coo.dec.value)
        assert np.absolute(spec.wave[0] - meta['wave0']) < 0.01, 'Bad wave'
        assert np.absolute(np.median(spec.flux) - meta['med_flux']) < 0.01, 'Bad flux'
        assert spec.meta['source'] == archive, 'Bad source'

        # Should work for just the file name, pulling from the archive
        spec = archive_classes[archive](file)
        # Or provided the full path
        _spec = archive_classes[archive](file_path)
        assert np.array_equal(spec.wave, _spec.wave), 'Both reads should give the same result'
        assert np.array_equal(spec.flux, _spec.flux), 'Both reads should give the same result'

        # Test the values
        assert len(spec.wave) == meta['len']
        assert np.absolute(spec.wave[0] - meta['wave0']) < 0.01, 'Bad wave'
        assert np.absolute(np.median(spec.flux) - meta['med_flux']) < 0.01, 'Bad flux'


def test_blackbody():

    assert standard.BlackbodyStandard.model_type == 'blackbody', 'Name changed'

    # Test spectrum
    bb = standard.BlackbodyStandard(4.0, 1e4)
    assert len(bb.wave) == 250880, 'Default length of spectrum changed'
    assert np.absolute(bb.wave[0] - 912.) < 1e-10, 'Initial wavelength changed'
    assert np.absolute((bb.flux[0] - 1.063)/bb.flux[0]) < 1e-4, 'Bad flux calculation'
    assert bb.meta['source'] == 'blackbody', 'Bad source'

    _bb = standard.BlackbodyStandard(4.0, 1e4, wave=bb.wave[:10])
    assert len(_bb.wave) == 10, 'Length wrong'
    assert np.array_equal(_bb.wave, bb.wave[:10]), 'Bad wavelength array'
    assert np.array_equal(_bb.flux, bb.flux[:10]), 'Bad fluxes'

    sep, name, a, teff = standard.BlackbodyStandard.nearest_blackbody_coeffs(180., 0.)
    assert name == 'BB113337-110529', 'Bad nearest'
    assert np.absolute(a - 3.5) < 0.01, 'Bad a for nearest'

    # Not within the defined tolerance
    with pytest.raises(PypeItError):
        bb = standard.BlackbodyStandard.from_coordinates(180., 0.)

    # Increase the tolerance
    bb = standard.BlackbodyStandard.from_coordinates(180., 0., tol=800.)
    assert np.absolute((bb.flux[0] - 3.59)/bb.flux[0]) < 0.01, 'Bad fluxes'
    assert bb.meta['source'] == 'blackbody', 'Bad source'

    # Use the name
    bb = standard.BlackbodyStandard.from_name('BB113337-110529')
    assert np.absolute((bb.flux[0] - 3.59)/bb.flux[0]) < 0.01, 'Bad fluxes'
    assert bb.meta['source'] == 'blackbody', 'Bad source'

    # Get the nearest set of coordinates and to make sure the search finds
    # something within the tolerance
    sep, row = standard.nearest_archive_entry('blackbody', 180.0, 0.0)
    coo = coordinates.SkyCoord(row['RA_2000'], row['DEC_2000'], unit=(units.hourangle, units.deg))
    bb = standard.BlackbodyStandard.from_coordinates(coo.ra.value, coo.dec.value)
    assert len(bb.wave) == 250880, 'Default length of spectrum changed'
    assert np.absolute((bb.flux[0] - 3.59)/bb.flux[0]) < 0.01, 'Bad fluxes'
    _bb = standard.BlackbodyStandard.from_coordinates(
        coo.ra.value, coo.dec.value, wave=bb.wave[:10]
    )
    assert len(_bb.wave) == 10, 'Length wrong'
    assert np.array_equal(_bb.wave, bb.wave[:10]), 'Bad wavelength array'
    assert np.array_equal(_bb.flux, bb.flux[:10]), 'Bad fluxes'
    _bb = standard.BlackbodyStandard.from_name(bb.meta['Name'], wave=bb.wave[:10])
    assert len(_bb.wave) == 10, 'Length wrong'
    assert np.array_equal(_bb.wave, bb.wave[:10]), 'Bad wavelength array'
    assert np.array_equal(_bb.flux, bb.flux[:10]), 'Bad fluxes'


def test_kurucz_models():
    spec = standard.KuruczModelStandard(14., 'B0')

    assert spec.size == 1221, 'Spectrum has the wrong length'
    assert np.absolute(spec.wave[0] - 90.9) < 0.01, 'Bad wave'
    assert np.absolute(np.median(spec.flux) - 90.1) < 0.01, 'Bad flux'
    assert spec.meta['source'] == 'Kurucz', 'Bad model source'

    # Unavailable spectral type
    with pytest.raises(PypeItError):
        spec = standard.KuruczModelStandard(14., 'K0')


def test_vega_model():
    spec = standard.VegaStandard(14.0)
    assert spec.size == 21617, 'Spectrum has the wrong length'
    assert np.absolute(spec.wave[0] - 900.092) < 0.01, 'Bad wave'
    assert np.absolute(np.median(spec.flux) - 2.11) < 0.01, 'Bad flux'
    assert spec.meta['source'] == 'Vega', 'Bad model source'


def test_phoenix_model():
    spec = standard.PhoenixStandard(14.0)
    assert spec.size == 1554127, 'Spectrum has the wrong length'
    assert np.absolute(spec.wave[0] - 2000.1) < 0.01, 'Bad wave'
    assert np.absolute(np.median(spec.flux) - 201.3) < 0.01, 'Bad flux'
    assert spec.meta['source'] == 'PHOENIX', 'Bad model source'


def test_pseudo_model():
    spec = standard.PseudoStandard()
    assert spec.size == 48000, 'Spectrum has the wrong length'
    assert np.absolute(spec.wave[0] - 2000.0) < 0.01, 'Bad wave'
    assert np.absolute(np.median(spec.flux) - 1.0) < 0.01, 'Bad flux'
    assert spec.meta['source'] == 'pseudo', 'Bad model source'

    _spec = standard.PseudoStandard(wave=spec.wave[:10])
    assert _spec.size == 10, 'Spectrum has the wrong length'
    assert np.absolute(_spec.wave[0] - 2000.0) < 0.01, 'Bad wave'
    assert np.absolute(np.median(_spec.flux) - 1.0) < 0.01, 'Bad flux'


def test_archive_sets():
    archive = standard.get_archive_sets()
    assert len(archive) == 5, 'Default list of archives changed'
    assert archive[0] == 'xshooter', 'Order changed'

    archive = standard.get_archive_sets(archives=['esofil', 'calspec'])
    assert len(archive) == 2, 'Incorrect number of archives'
    assert archive[1] == 'calspec', 'Order incorrect'

    # Allow a single archive
    archive = standard.get_archive_sets(archives=['calspec'])
    assert len(archive) == 1, 'Incorrect number of archives'
    assert archive[0] == 'calspec', 'Order incorrect'

    # Ignore invalid archives
    archive = standard.get_archive_sets(archives=['esofil', 'junk'])
    assert len(archive) == 1, 'Incorrect number of archives'
    assert archive[0] == 'esofil', 'Order incorrect'

    # Fault if none are valid
    with pytest.raises(PypeItError):
        archive = standard.get_archive_sets(archives=['junk1', 'junk2'])


def test_get_archive_standard():

    # No star within the tolerance
    with pytest.raises(PypeItError):
        spec = standard.get_archive_standard(180., 0.)

    # Force it to use the blackbody "archive"
    ra, dec = '12:45:35.622', '+42:38:24.675'
    spec = standard.get_archive_standard(ra, dec, archives='blackbody')
    assert isinstance(spec, standard.BlackbodyStandard), 'Has the wrong type'
    assert spec.meta['Name'] == 'BB124535+423824', 'Found the wrong object'

    # Make sure it finds the same if the archive is not specified
    spec = standard.get_archive_standard(ra, dec)
    assert isinstance(spec, standard.BlackbodyStandard), 'Has the wrong type'
    assert spec.meta['Name'] == 'BB124535+423824', 'Found the wrong object'

    # Make sure it faults if blackbody is excluded
    with pytest.raises(PypeItError):
        spec = standard.get_archive_standard(ra, dec, archives=standard.get_archive_sets())

    # Find an XShooter star
    ra, dec = '12:57:02.34', '+22:01:52.7'
    spec = standard.get_archive_standard(ra, dec)
    assert isinstance(spec, standard.XShooterFluxStandard), 'Has the wrong type'
    assert spec.meta['Name'] == 'GD153', 'Found the wrong object'

    # Force it to find the calspec version of GD153
    spec = standard.get_archive_standard(ra, dec, archives='calspec')
    assert isinstance(spec, standard.CalSpecFluxStandard), 'Has the wrong type'
    assert spec.meta['Name'] == 'GD153', 'Found the wrong object'

    # Find an calspec star
    ra, dec = '12:53:15.053', '-18:31:20.01'
    spec = standard.get_archive_standard(ra, dec)
    assert isinstance(spec, standard.CalSpecFluxStandard), 'Has the wrong type'
    assert spec.meta['Name'] == 'HD111980', 'Found the wrong object'

    # Find an esofil star
    ra, dec = '12:06:47.25', '11:40:12.7'
    spec = standard.get_archive_standard(ra, dec)
    assert isinstance(spec, standard.ESOFilFluxStandard), 'Has the wrong type'
    assert spec.meta['Name'] == 'Feige56', 'Found the wrong object'

    # Find an NOAO star
    ra, dec = '11:37:05.06', '+29:47:58.2'
    spec = standard.get_archive_standard(ra, dec)
    assert isinstance(spec, standard.NOAOFluxStandard), 'Has the wrong type'
    assert spec.meta['Name'] == 'GD140', 'Found the wrong object'

    # Find an ING star
    ra, dec = '04:40:39.32', '+08:40:45.3'
    spec = standard.get_archive_standard(ra, dec)
    assert isinstance(spec, standard.INGFluxStandard), 'Has the wrong type'
    assert spec.meta['Name'] == 'HZ15', 'Found the wrong object'


def test_get_model_standard():

    # Get the TSpecTool model
    spec = standard.get_model_standard('A0', 14.)
    assert isinstance(spec, standard.VegaStandard), 'Wrong type'
    assert spec.meta['V_mag'] == 14., 'Wrong metadata'

    # Get the PHOENIX model
    spec = standard.get_model_standard('PHOENIX', 14.)
    assert isinstance(spec, standard.PhoenixStandard), 'Wrong type'
    assert spec.meta['V_mag'] == 14., 'Wrong metadata'

    # Get the flat continuum
    spec = standard.get_model_standard('NONE', 14.)
    assert isinstance(spec, standard.PseudoStandard), 'Wrong type'

    # Get a Kurucz model
    spec = standard.get_model_standard('G0', 14.)
    assert isinstance(spec, standard.KuruczModelStandard), 'Wrong type'
    assert spec.meta['V_mag'] == 14., 'Wrong metadata'
    assert spec.meta['Sp'] == 'G0', 'Wrong spectral type'

    # Make sure it faults on an unknown spectral type
    with pytest.raises(PypeItError):
        spec = standard.get_model_standard('K0', 14.)


def test_get_standard_spectrum():

    # Get the TSpecTool model
    spec = standard.get_standard_spectrum(spectral_type='A0', V_mag=14.)
    assert isinstance(spec, standard.VegaStandard), 'Wrong type'
    assert spec.meta['V_mag'] == 14., 'Wrong metadata'

    # Get an XShooter star
    ra, dec = '12:57:02.34', '+22:01:52.7'
    spec = standard.get_standard_spectrum(ra=ra, dec=dec)
    assert isinstance(spec, standard.XShooterFluxStandard), 'Has the wrong type'
    assert spec.meta['Name'] == 'GD153', 'Found the wrong object'

    # Make sure it faults if nothing is provided
    with pytest.raises(PypeItError):
        spec = standard.get_standard_spectrum()

