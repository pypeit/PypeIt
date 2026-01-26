"""
Implements classes for flux standard spectra used for flux calibration.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from astropy import constants
from astropy import coordinates
from astropy import table
from astropy import units
from astropy.io import fits
from astropy.io import ascii
from IPython import embed
import numpy as np

from pypeit import msgs
from pypeit import dataPaths
from pypeit.pypmsgs import PypeItError
from pypeit.core import spectrum
from pypeit.core.meta import convert_radec
from pypeit.core.wave import airtovac
from pypeit.utils import all_subclasses


def mAB_to_cgs(wave, mAB):
    r"""
    Convert AB magnitudes to :math:`F_\lambda` in the cgs units :math:`{\rm
    erg/cm}^2{\rm/s}/\AA`.

    Parameters
    ----------
    wave: scalar-like, array-like
        Vacuum Wavelength in Angstrom.
    mAB : scalar-like, array-like
        AB magnitudes.  If array-like, must be possible to broadcast to match
        ``wave``.

    Returns
    -------
    float, `numpy.ndarray`_
        Flux density in cgs units.  Returned as scalar or array depending on
        input.
    """
    _mAB = mAB if isinstance(mAB, (float, np.floating, int, np.integer)) else np.asarray(mAB)
    _wave = wave if isinstance(wave, (float, np.floating, int, np.integer)) else np.asarray(wave)
    return 10**((-48.6-_mAB)/2.5) * 3e18 / _wave**2


def archive_entry(archive, name):
    """
    Find the row of data in the specified archive associated with the named
    source.

    Parameters
    ----------
    archive : str
        Name of the archive to search; see :func:`get_archive_sets`.  This
        function also works for `archive='blackbody'`.
    name : str
        The name of the archive source.  Must be an exact match.

    Returns
    -------
    `astropy.table.Row`_
        Single table row with the data from the archive.
    """
    # Set the path (creates a new PypeItDataPath object)
    stds_path = dataPaths.standards / archive
    # Get the file
    star_file = stds_path.get_file_path(f'{archive}_info.txt')
    if not star_file.is_file():
        msgs.error(f'File does not exist!: {star_file}')

    star_tbl = table.Table.read(star_file, comment='#', format='ascii')
    idx = np.where(star_tbl['Name'] == name)[0]
    if len(idx) != 1:
        msgs.error(f'{name} is not a named source in {star_file}.')
    return star_tbl[idx[0]]


def nearest_archive_entry(archive, ra, dec, unit=None):
    """
    Find the row of data in the specified archive with coordinates nearest
    to the provided coordinates.

    Parameters
    ----------
    archive : str
        Name of the archive to search; see :func:`get_archive_sets`.  This
        function also works for `archive='blackbody'`.
    ra, dec : float, str
        On-sky coordinates.  If ``units`` are None, the coordinates are assumed
        to be in degree if provided as floats, and they are assumed to be in
        hours and degrees if provided as (e.g., sexagesimal) strings.
    unit : str, tuple, optional
        Units for the on-sky coordinates.  See ``ra`` and ``dec`` for the
        default behavior.

    Returns
    -------
    float
        Separation between the nearest standard star and the provided coordinates.
    `astropy.table.Row`_
        Single table row with the data from the archive.
    """
    # Instatiate the coordinates
    if unit is None:
        _unit = units.deg if isinstance(ra, (float, np.floating, int, np.integer)) \
                    else (units.hourangle, units.deg)
    else:
        _unit = unit
    obj_coord = coordinates.SkyCoord([ra], [dec], unit=_unit)
    if obj_coord.size > 1:
        msgs.error('Matching to archive can only be done one object at a time.')

    # Set the path (creates a new PypeItDataPath object)
    stds_path = dataPaths.standards / archive
    # Get the file
    star_file = stds_path.get_file_path(f"{archive}_info.txt")
    if not star_file.is_file():
        msgs.error(f"File does not exist!: {star_file}")

    star_tbl = table.Table.read(star_file, comment='#', format='ascii')
    star_coords = coordinates.SkyCoord(star_tbl['RA_2000'], star_tbl['DEC_2000'],
                                       unit=(units.hourangle, units.deg))
    # This returns 
    idx, d2d = coordinates.match_coordinates_sky(obj_coord, star_coords, nthneighbor=1)[:2]
    return d2d[0], star_tbl[idx[0]]


class ArchivedFluxStandard(spectrum.Spectrum):
    """
    Abstract class used to provide common methods for all archive standards.
    """
    
    archive = None
    """
    Archive identifier
    """

    path = None
    """
    Root with data files
    """

    @classmethod
    def nearest_standard(cls, ra, dec, unit=None):
        """
        Find the standard star with an archived flux-calibration spectrum
        nearest the provided set of coordinates.

        Parameters
        ----------
        ra, dec : float, str
            On-sky coordinates.  If ``units`` are None, the coordinates are assumed
            to be in degree if provided as floats, and they are assumed to be in
            hours and degrees if provided as (e.g., sexagesimal) strings.
        unit : str, tuple, optional
            Units for the on-sky coordinates.  See ``ra`` and ``dec`` for the
            default behavior.

        Returns
        -------
        float
            Separation between the nearest standard star and the provided coordinates.
        str
            Name of the object in the archive.
        str
            Name of the file with the archived spectrum.
        """
        sep, row = nearest_archive_entry(cls.archive, ra, dec, unit=unit)
        return sep, row['Name'], row['File']

    @classmethod
    def found_match(cls, ra, dec, tol=20., unit=None):
        """
        Check if there is a match to the provided coordinates within the
        tolerance.

        Parameters
        ----------
        ra, dec : float, str
            On-sky coordinates.  If ``units`` are None, the coordinates are assumed
            to be in degree if provided as floats, and they are assumed to be in
            hours and degrees if provided as (e.g., sexagesimal) strings.
        tol : float, optional
            Tolerance for coordinate matching in arcmin
        unit : str, tuple, optional
            Units for the on-sky coordinates.  See ``ra`` and ``dec`` for the
            default behavior.

        Returns
        -------
        bool
            Flag that an appropriate match was found.
        """
        sep, row = nearest_archive_entry(cls.archive, ra, dec, unit=unit)
        return sep < tol * units.arcmin
    
    @classmethod
    def from_coordinates(cls, ra, dec, tol=20., unit=None):
        """
        Instantiate the class using the spectrum for the object closest to the
        provided set of coordinates.

        Parameters
        ----------
        ra, dec : float, str
            On-sky coordinates.  If ``units`` are None, the coordinates are assumed
            to be in degree if provided as floats, and they are assumed to be in
            hours and degrees if provided as (e.g., sexagesimal) strings.
        tol : float, optional
            Tolerance for coordinate matching in arcmin
        unit : str, tuple, optional
            Units for the on-sky coordinates.  See ``ra`` and ``dec`` for the
            default behavior.
        """
        sep, row = nearest_archive_entry(cls.archive, ra, dec, unit=unit)
        if sep > tol * units.arcmin:
            msgs.error(f'Closest object ({row["Name"]}) in archive "{cls.archive}" is '
                       f'separated by {sep.to("arcmin").value:.2f} arcmin, which is beyond the '
                       f'required tolerance ({tol} arcmin).')
        else:
            msgs.info(f'Using object {row["Name"]} in archive "{cls.archive}" '
                      f'(sep = {sep.to("arcmin").value:.2f} arcmin)')
        return cls(row['File'], meta=cls._init_meta(row=row))
    
    @classmethod
    def from_name(cls, name):
        """
        Instantiate from the name of an archive source.

        Parameters
        ----------
        name : str
            Name of the source in the archive data table.  Must be an exact
            match.
        """
        row = archive_entry(cls.archive, name)
        return cls(row['File'], meta=cls._init_meta(row=row))
    
    @classmethod
    def _init_meta(cls, row=None):
        """
        Instantiate the metadata.
        """
        # Add all of the tabulated data to the metadata for this spectrum
        meta = {} if row is None else dict(row)
        # Also add the "source"
        meta['source'] = cls.archive

        # If the coordinates are in the row, add entries that convert the
        # coordinates to degrees
        if 'RA_2000' in meta and 'DEC_2000' in meta:
            meta['ra_deg'], meta['dec_deg'] = convert_radec(meta['RA_2000'], meta['DEC_2000'])

        return meta


class CalSpecFluxStandard(ArchivedFluxStandard):
    """
    Container class for a "calspec" standard star spectrum.
    """
    archive = 'calspec'
    path = dataPaths.standards / archive

    def __init__(self, file, meta=None):
        self.file = self.path.get_file_path(file)
        with fits.open(self.file) as hdu:
            wave = hdu[1].data['WAVELENGTH']
            flux = hdu[1].data['FLUX'] * 1e17
        super().__init__(wave, flux, meta=meta)


class ESOFilFluxStandard(ArchivedFluxStandard):
    """
    Container class for an "esofil" standard star spectrum.
    """
    archive = 'esofil'
    path = dataPaths.standards / archive

    def __init__(self, file, meta=None):
        self.file = self.path.get_file_path(file)
        std_spec = table.Table.read(self.file, format='ascii')
        wave = std_spec['col1']
        flux = std_spec['col2'] * 10    # Convert from 1e-16 to 1e-17 erg/s/cm^2/Angstrom

        # At this low resolution, best to throw out entries affected by A and B-band absorption
        gpm = np.logical_not((wave > 7551.) & (wave < 7749.))
        super().__init__(wave[gpm], flux[gpm], meta=meta)


class INGFluxStandard(ArchivedFluxStandard):
    """
    Container class for an "ing" standard star spectrum.
    """
    archive = 'ing'
    path = dataPaths.standards / archive

    def __init__(self, file, meta=None):
        self.file = self.path.get_file_path(file)
        std_spec = table.Table.read(self.file, format='ascii')
        wave = std_spec['col1']
        flux = mAB_to_cgs(std_spec['col1'], std_spec['col2']) * 1e17

        # At this low resolution, best to throw out entries affected by A and B-band absorption
        gpm = np.logical_not((wave > 7551.) & (wave < 7749.))
        super().__init__(wave[gpm], flux[gpm], meta=meta)


class NOAOFluxStandard(ArchivedFluxStandard):
    """
    Container class for an "noao" standard star spectrum.
    """
    archive = 'noao'
    path = dataPaths.standards / archive

    def __init__(self, file, meta=None):

        self.file = self.path.get_file_path(file)
        std_spec = table.Table.read(self.file, format='ascii')
        wave = std_spec['col1']
        flux = mAB_to_cgs(std_spec['col1'], std_spec['col2']) * 1e17

        # At this low resolution, best to throw out entries affected by A and B-band absorption
        gpm = np.logical_not((wave > 7551.) & (wave < 7749.))
        super().__init__(wave[gpm], flux[gpm], meta=meta)


class XShooterFluxStandard(ArchivedFluxStandard):
    """
    Container class for an "xshooter" standard star spectrum.
    """
    archive = 'xshooter'
    path = dataPaths.standards / archive

    def __init__(self, file, meta=None):
        self.file = self.path.get_file_path(file)
        std_spec = table.Table.read(self.file, format='ascii')
        # XShooter standard files use air wavelengths, convert them to vacuum
        wave = airtovac(std_spec['col1'] * units.AA).value
        flux = std_spec['col2'] * 1e17
        super().__init__(wave, flux, meta=meta)


def archived_flux_classes():
    """
    Construct a dictionary with the set of classes that subclass from
    :class:`~pypeit.core.standard.ArchivedFluxStandard`.

    Returns
    -------
    dict
        Dictionary with keys that identify the archive name and class.
    """
    # Recursively collect all subclasses
    c = np.array(list(all_subclasses(ArchivedFluxStandard)))
    # Select classes with a defined archive name
    c = c[[cls.archive is not None for cls in c]]
    # Construct a dictionary with the archive and class
    srt = np.argsort(np.array([cls.archive for cls in c]))
    return {cls.archive:cls for cls in c[srt]}


class ModelFluxStandard(spectrum.Spectrum):
    """
    Base class for "model-based" flux standard spectra.
    """
    model_type = None
    required_metadata = ['Name', 'File', 'ra_deg', 'dec_deg']

    @classmethod
    def _init_meta(cls, row=None):
        """
        Instantiate the metadata.
        """
        # Add all of the tabulated data to the metadata for this spectrum
        meta = {} if row is None else dict(row)
        # Also add the "source"
        meta['source'] = cls.model_type

        # If the coordinates are in the row, add entries that convert the
        # coordinates to degrees
        if 'RA_2000' in meta and 'DEC_2000' in meta:
            meta['ra_deg'], meta['dec_deg'] = convert_radec(meta['RA_2000'], meta['DEC_2000'])

        # Add in required meta
        for key in cls.required_metadata:
            if key not in meta.keys():
                meta[key] = None
        return meta


class BlackbodyStandard(ModelFluxStandard):
    """
    Generate a blackbody spectrum based on the normalisation and effective
    temperature.  See Suzuki & Fukugita, 2018, AJ, 156, 219:
    https://ui.adsabs.harvard.edu/abs/2018AJ....156..219S/abstract

    Parameters
    ----------
    a : float
        flux normalisation factor (dimensionless)
    teff : float
        Effective temperature of the blackbody in Kelvin
    wave : array-like, optional
        Vacuum wavelength in angstroms at which to calculate the blackbody flux.
        If None, the default wavelength range is set to 912 - 26000 Angstrom at
        a step of 0.1 Angstrom.
    meta : dict, optional
        The metadata to keep with the spectrum.
    """

    model_type = 'blackbody'
    """
    Identifier for the type of model spectrum.
    """

    def __init__(self, a, teff, wave=None, meta=None):
        # TODO: Simplify the unit stuff here!
        if wave is None:
            resln = 0.1  # Resolution to generate the blackbody spectrum
            _wave = np.arange(912.0, 26000.0, resln) * units.AA
        else:
            _wave = np.asarray(wave) * units.AA
        _teff = teff * units.K
        # Calculate the function
        flam = (a * 2 * constants.h * constants.c**2 / _wave**5) / (
            np.exp((constants.h * constants.c / 
                   (_wave * constants.k_B * _teff)).to(units.m/units.m).value
            ) - 1.0
        )
        # Convert to 1e-17 erg/s/cm^2/Angstrom, and apply the "BB_SCALE_FACTOR" (1e-23)
        flam = flam.to(units.erg / units.s / units.cm ** 2 / units.AA).value * 1e-6
        super().__init__(_wave.value, flam, meta=self._init_meta(row=meta))

    @classmethod
    def nearest_blackbody_coeffs(cls, ra, dec, unit=None):
        """
        Find the entry in the blackbody reference table nearest to the provided
        coordinates and return the angular separation and the coefficients
        needed to construct the model.

        Parameters
        ----------
        ra, dec : float, str
            On-sky coordinates.  If ``units`` are None, the coordinates are assumed
            to be in degree if provided as floats, and they are assumed to be in
            hours and degrees if provided as (e.g., sexagesimal) strings.
        unit : str, tuple, optional
            Units for the on-sky coordinates.  See ``ra`` and ``dec`` for the
            default behavior.

        Returns
        -------
        sep : Quantity
            Angular separation between the provided coordinates and the nearest
            entry in the blackbody reference table.
        name : str
            Name of the blackbody reference object.
        a, teff : float
            Parameters used to construct the blackbody spectrum.
        """
        sep, row = nearest_archive_entry(cls.model_type, ra, dec, unit=unit)
        return sep, row['Name'], row['a_x10m23'], row['T_K']

    # TODO: Consolidate this with ArchivedFluxStandard?
    @classmethod
    def found_match(cls, ra, dec, tol=20., unit=None):
        """
        Check if there is a match to the provided coordinates within the
        tolerance.

        Parameters
        ----------
        ra, dec : float, str
            On-sky coordinates.  If ``units`` are None, the coordinates are assumed
            to be in degree if provided as floats, and they are assumed to be in
            hours and degrees if provided as (e.g., sexagesimal) strings.
        tol : float, optional
            Tolerance for coordinate matching in arcmin
        unit : str, tuple, optional
            Units for the on-sky coordinates.  See ``ra`` and ``dec`` for the
            default behavior.

        Returns
        -------
        bool
            Flag that an appropriate match was found.
        """
        sep, row = nearest_archive_entry(cls.model_type, ra, dec, unit=unit)
        return sep < tol * units.arcmin

    @classmethod
    def from_coordinates(cls, ra, dec, tol=20., unit=None, wave=None):
        """
        Instantiate the class using coefficients for the object closest to the
        provided set of coordinates.

        Parameters
        ----------
        ra, dec : float, str
            On-sky coordinates.  If ``units`` are None, the coordinates are assumed
            to be in degree if provided as floats, and they are assumed to be in
            hours and degrees if provided as (e.g., sexagesimal) strings.
        tol : float, optional
            Tolerance for coordinate matching in arcmin
        unit : str, tuple, optional
            Units for the on-sky coordinates.  See ``ra`` and ``dec`` for the
            default behavior.
        wave : array-like, optional
            Vacuum wavelength in angstroms at which to calculate the blackbody flux.
            If None, the default wavelength range is set to 912 - 26000 Angstrom at
            a step of 0.1 Angstrom.
        """
        sep, row = nearest_archive_entry(cls.model_type, ra, dec, unit=unit)
        if sep > tol * units.arcmin:
            msgs.error(f'Closest object ({row["Name"]}) is separated by {sep.to("arcmin").value} '
                       f'arcmin, which is beyond the required tolerance ({tol} arcmin).')
        return cls(row['a_x10m23'], row['T_K'], wave=wave, meta=cls._init_meta(row=row))

    @classmethod
    def from_name(cls, name, wave=None):
        """
        Instantiate from the name of an archive source.

        Parameters
        ----------
        name : str
            Name of the source in the archive data table.  Must be an exact
            match.
        """
        row = archive_entry(cls.model_type, name)
        return cls(row['a_x10m23'], row['T_K'], wave=wave, meta=cls._init_meta(row=row))


class KuruczModelStandard(ModelFluxStandard):
    """
    The Kurucz stellar model for a given apparent magnitude and spectral type.

    The spectrum is instantiated as follows:

        - Get the temperature, logg, and bolometric luminosity from the
          Schmidt-Kaler (1982) table for the provided spectral type.

        - Find the nearest neighbor in the Kurucz stellar atmosphere ATLAS.

        - Convert the units for the wavelength and flux.

    .. warning::

        Spectra can currently only be generated for spectral types provided in
        the Schmidt-Kaler (1982) table.  See
        `here <https://github.com/pypeit/PypeIt/blob/release/pypeit/data/standards/kurucz93/schmidt-kaler_table.txt>`__
        for available spectral types.

    Parameters
    ----------
    V_mag : float
        Apparent magnitude of the star in the V band (Vega system).
    spectral_type : str
        Stellar spectral type
    """

    model_type = 'Kurucz'
    """
    Identifier for the type of model spectrum.
    """

    def __init__(self, V_mag, spectral_type):

        # Load Schmidt-Kaler (1982) table
        sk82_file = dataPaths.standards.get_file_path('kurucz93/schmidt-kaler_table.txt')
        sk82_tab = ascii.read(
            sk82_file,
            names=('Sp', 'logTeff', 'Teff', '(B-V)_0', 'M_V', 'B.C.', 'M_bol', 'L/L_sol')
        )

        # Match input type.
        # TODO: currently this only works on select stellar types. Add ability to
        # interpolate across types.
        indx = np.where(spectral_type == sk82_tab['Sp'])[0]
        if len(indx) != 1:
            msgs.error(
                f'Provided spectral type {spectral_type} not available in Schmidt-Kaler (1982) '
                'table.  See the KuruczModelStandard API.'
            )
        indx = indx[0]

        # log(g) of the Sun
        logg_sol = np.log10(6.67259e-8) + np.log10(1.989e33) - 2.0 * np.log10(6.96e10)

        # Relation between radius, temp, and bolometric luminosity
        logR = 0.2 * (42.26 - sk82_tab['M_bol'][indx] - 10.0 * sk82_tab['logTeff'][indx])

        # Mass-bolometric luminosity relation from schimdt-kaler p28 valid for M_bol < 7.5
        logM = 0.46 - 0.10 * sk82_tab['M_bol'][indx]
        logg = logM - 2.0 * logR + logg_sol
        M_V = sk82_tab['M_V'][indx]
        Teff = sk82_tab['Teff'][indx]

        # Distance modulus
        logd = 0.2 * (V_mag - M_V) + 1.0
        D = constants.pc.cgs.value * 10. ** logd
        R = constants.R_sun.cgs.value * 10. ** logR

        # Factor converts the Kurucz surface flux densities to flux observed on Earth
        flux_factor = (R / D) ** 2

        # Grab closest T in Kurucz SEDs
        T1 = 3000. + np.arange(28) * 250
        T2 = 10000. + np.arange(6) * 500
        T3 = 13000. + np.arange(22) * 1000
        T4 = 35000. + np.arange(7) * 2500
        Tk = np.concatenate([T1, T2, T3, T4])
        indT = np.argmin(np.abs(Tk - Teff))

        # Grab closest g in Kurucz SEDs
        loggk = np.arange(11) * 0.5
        indg = np.argmin(np.abs(loggk - logg))
        gdict = { 0: 'g00', 1: 'g05', 2: 'g10', 3: 'g15',  4: 'g20', 5: 'g25',
                  6: 'g30', 7: 'g35', 8: 'g40', 9: 'g45', 10: 'g50'}

        # Grab Kurucz filename
        meta = dict(sk82_tab[indx])
        meta['V_mag'] = V_mag
        std_file = dataPaths.standards.get_file_path(f'kurucz93/kp00/kp00_{int(Tk[indT])}.fits.gz')
        with fits.open(std_file) as hdu:
            return super().__init__(
                hdu[1].data['WAVELENGTH'],
                hdu[1].data[gdict[indg]] * flux_factor * 1e17,
                meta=self._init_meta(row=meta)
            )


class VegaStandard(ModelFluxStandard):
    """
    Provides a Vega spectrum from TSpecTool.

    Parameters
    ----------
    V_mag : float
        The V-band magnitude for the star.
    """

    model_type = 'Vega'
    """
    Identifier for the type of model spectrum.
    """

    def __init__(self, V_mag):
        file = dataPaths.standards.get_file_path('vega_tspectool_vacuum.dat')
        data = table.Table.read(file, comment='#', format='ascii')
        return super().__init__(
            data['col1'],
            data['col2'] * 10**(0.4*(0.03-V_mag)) * 1e17,
            meta=self._init_meta(row={'V_mag': V_mag})
        )


class PhoenixStandard(ModelFluxStandard):
    """
    Provides the PHOENIX spectrum.

    Parameters
    ----------
    V_mag : float
        The V-band magnitude for the star.
    """

    model_type = 'PHOENIX'
    """
    Identifier for the type of model spectrum.
    """

    def __init__(self, V_mag):
        file = dataPaths.standards.get_file_path('PHOENIX_10000K_4p0.fits')
        data = table.Table.read(file, format='fits')
        return super().__init__(
            data['Wavelength'],
            data['Flux'] * 10**(0.4*(0.03-V_mag)) * 1e6,
            meta=self._init_meta(row={'V_mag': V_mag})
        )
    

class PseudoStandard(ModelFluxStandard):
    """
    Provides a unity continuum spectrum.

    Parameters
    ----------
    wave : array-like, optional
        The wavelength array to use.  If None, wavelengths range from 0.2-5
        micron in steps of 1 Angstrom.
    """

    model_type = 'pseudo'
    """
    Identifier for the type of model spectrum.
    """

    def __init__(self, wave=None):
        _wave = np.arange(2000,50000,1.0) if wave is None else np.asarray(wave)
        return super().__init__(_wave, np.ones(_wave.shape, dtype=float), meta=self._init_meta())


def get_archive_sets(archives=['xshooter', 'calspec', 'esofil', 'noao', 'ing']):
    """
    Helper function to setup the prioritized list of archive sets to search
    through when matching a set of coordinates to a file containing the flux
    standard data.

    Parameters
    ----------
    archives : array-like, str, optional
        Name of the archives to search, in a prioritized order.  If None, all
        archives are searched.

    Returns
    -------
    `numpy.ndarray`_
        The list of standard sets to search

    Raises
    ------
    PypeItError
        Raised if none of the provided sets are recognized.
    """
    archive_classes = archived_flux_classes()

    _archives = np.asarray([archives] if isinstance(archives, str) else archives)
    good = np.ones(len(_archives), dtype=bool)
    for i, s in enumerate(_archives):
        if s not in archive_classes.keys():
            msgs.warn(f'{s} is not a recognized archive of standard spectra.  Ignoring.')
            good[i] = False
    if not any(good):
        msgs.error('None of the provided standard spectra archives are valid.  Try using '
                   'the default list.')
    return _archives[good]


def get_archive_standard(ra, dec, tol=20., unit=None, archives='default', check=False):
    """
    Attempt to find and return an archive flux calibration spectrum that is
    closest to the provided coordinates.

    The archives searched *always* start with the empirical archives (see
    :func:`~pypeit.core.standard.get_archive_standard`).  If no archive match is
    found, the function attempts (unless explicitly excluded via the
    ``archives`` argument) to find a suitable set of blackbody parameters; see
    :class:`~pypeit.core.standard.BlackbodyStandard`.
    
    Parameters
    ----------
    ra, dec : float, str
        On-sky coordinates.  If ``units`` are None, the coordinates are assumed
        to be in degree if provided as floats, and they are assumed to be in
        hours and degrees if provided as (e.g., sexagesimal) strings.
    tol : float, optional
        Tolerance for coordinate matching in arcmin
    unit : str, tuple, optional
        Units for the on-sky coordinates.  See ``ra`` and ``dec`` for the
        default behavior.
    archives : array-like, str, optional
        Name of the archives to search, in a prioritized order.  If
        ``'default'``, all archives are searched.  To only search for suitable
        blackbody parameters, set ``archives='blackbody'``.
    check : :obj:`bool`, optional
        Only check if a standard matches the provided coordinates, as opposed to
        also reading the spectral data.

    Returns
    -------
    spectrum.Spectrum, bool
        If ``check`` is True, a flag is returned indicating if a standard
        matches the provided coordinates.  Otherwise, the standard spectrum is
        returned.
    """
    archive_classes = archived_flux_classes()
    # NOTE: The if statement below has to use `'default' in archives` to allow
    # for archives to be provided as numpy array; i.e., the statement `archives
    # == default` yields a boolean array.  The statement as is works for
    # strings, lists, and arrays.  If `archives` is provided as a list that
    # includes "default", it will effectively ignore anything that isn't in the
    # default list of archives...
    try:
        _archives = (
            get_archive_sets() if 'default' in archives else get_archive_sets(archives=archives)
        )
    except PypeItError:
        # It is possible for this to fail if the only archive set is 'blackbody'
        _archives = np.asarray([])

    for key in _archives:
        if check:
            if archive_classes[key].found_match(ra, dec, tol=tol, unit=unit):
                return True
        else:
            try:
                return archive_classes[key].from_coordinates(ra, dec, tol=tol, unit=unit)
            except PypeItError:
                # Ignore PypeItErrors, assuming they're because there was no object
                # within tol
                continue

    # Try to find a nearby blackbody
    # NOTE: The use of "in" here follows the same reason as above.  And note
    # this is using `archives` (i.e., what was provided to the function) not
    # `_archives` (i.e., how the function parses the input).
    if 'default' in archives or 'blackbody' in archives:
        if check:
            return BlackbodyStandard.found_match(ra, dec, tol=tol, unit=unit)

        _archives = np.append(_archives, ['blackbody'])
        try:
            return BlackbodyStandard.from_coordinates(ra, dec, tol=tol, unit=unit)
        except PypeItError:
            # Ignore PypeItErrors, assuming they're because there was no object
            # within tol
            pass
    
    if check:
        return False

    # Unable to find a standard within the provided tolerance.  Find the closest
    # one, report it, and fault.
    res = np.asarray(list([nearest_archive_entry(key, ra, dec, unit=unit) for key in _archives]))
    indx = np.argmin(res[:,0])
    sep, row = res[indx]
    msgs.error(f'Unable to find a standard star within {tol:.1f} arcmin of RA={ra}, DEC={dec} in '
               f'the following archives: {_archives}.  The nearest object is {row["Name"]} in '
               f'{_archives[indx]} at RA={row["RA_2000"]}, DEC={row["DEC_2000"]}, separated by '
               f'{sep.to("arcmin").value:.1f} arcmin.')


def get_model_standard(spectral_type, V_mag):
    """
    Return a model flux standard based on the spectral type and V-band magnitude.

    The models returned are as follows:

    If ``spectral type`` is:

        - "A0": the TSpecTool spectrum of Vega is returned (see
          :class:`VegaStandard`);

        - "PHOENIX": a PHOENIX model of a Teff = 10k, log(g) = 4. star is
          returned (see :class:`PhoenixStandard`);

        - "NONE": A constant unity spectrum is returned (see
          :class:`PseudoStandard`);

        - otherwise, a Kurucz model is returned (see
          :class:`KuruczModelStandard`).

    Parameters
    ----------
    spectral_type : str
        The spectral type of the star or the signifier of the spectrum to use.
        See above.
    V_mag : float
        The V-band magnitude for the star.

    Returns
    -------
    spectrum.Spectrum
        The standard spectrum.
    """
    match spectral_type:
        case 'A0':
            return VegaStandard(V_mag)
        case 'PHOENIX':
            return PhoenixStandard(V_mag)
        case 'NONE':
            return PseudoStandard()
        case _:
            return KuruczModelStandard(V_mag, spectral_type)


def get_standard_spectrum(spectral_type=None, V_mag=None, ra=None, dec=None, tol=20., unit=None,
                          archives='default'):
    """
    Return a standard spectrum.

    Must provide either the spectral type and the magnitude (for use with
    :func:`~pypeit.core.standard.get_model_standard`) or a set of on-sky
    coordinates (for use with
    :func:`~pypeit.core.standard.get_archive_standard`).  If all four are
    provided, the spectral type and magnitude take precedence.

    Parameters
    ----------
    spectral_type : str, optional
        The spectral type of the star or the signifier of the spectrum to use.
        See :func:`~pypeit.core.standard.get_model_standard`.
    V_mag : float, optional
        The V-band magnitude for the star.
    ra, dec : float, str
        On-sky coordinates.  If ``units`` are None, the coordinates are assumed
        to be in degree if provided as floats, and they are assumed to be in
        hours and degrees if provided as (e.g., sexagesimal) strings.
    tol : float, optional
        Tolerance for coordinate matching in arcmin
    unit : str, tuple, optional
        Units for the on-sky coordinates.  See ``ra`` and ``dec`` for the
        default behavior.
    archives : array-like, str, optional
        Name of the archives to search, in a prioritized order.  If
        ``'default'``, all archives are searched.

    Returns
    -------
    spectrum.Spectrum
        Standard star spectrum
    """
    if spectral_type is not None and V_mag is not None:
        return get_model_standard(spectral_type, V_mag)
    if ra is None or dec is None:
        msgs.error('Insufficient data provided to determine the appropriate standard spectrum.  '
                   'Provide either the coordinates of the standard or a stellar type and '
                   'magnitude.')
    return get_archive_standard(ra, dec, tol=tol, unit=unit, archives=archives)

