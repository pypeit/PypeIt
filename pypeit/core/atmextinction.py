""" Atmospheric extinction class

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from IPython import embed

import numpy as np

from scipy import interpolate

from astropy import units
from astropy import coordinates
from astropy import table

from pypeit import log
from pypeit import dataPaths
from pypeit import PypeItError
from pypeit import utils


class AtmosphericExtinction:
    """
    Class used to load, generate, and apply atmospheric extinction curves.

    Parameters
    ----------
    wave : array-like
        Vacuum wavelengths in angstroms
    mag_ext : array-like
        Magnitudes of extinction at 1 airmass.  Must be the same length as
        ``wave``.
    assume_sorted : :obj:`bool`, optional
        Assume the wavelength vector is sorted
    file : str, optional
        *Used for informational purposes only.*  If the data were read from a
        file, this provides the name of the source file.

    Raises
    ------
    PypeItError :
        Raised if the length of ``wave`` and ``mag_ext`` is not identical, or if
        they are multidimensional.
    """
    def __init__(self, wave, mag_ext, assume_sorted=True, file=None):

        if len(wave) != len(mag_ext):
            raise PypeItError('Wavelength and extinction vectors must have the same length.')

        self.wave = np.asarray(wave, dtype=float)
        self.mag_ext = np.asarray(mag_ext, dtype=float)

        if self.wave.ndim != 1 or self.mag_ext.ndim != 1:
            raise PypeItError('Atmospheric extinction must be 1D.')

        if not assume_sorted:
            srt = np.argsort(self.wave)
            self.wave = self.wave[srt]
            self.mag_ext = self.mag_ext[srt]

        self.file = file

    @property
    def size(self):
        """Return the length of the extinction curve."""
        return self.wave.size
    
    @staticmethod
    def closest_extinction_file(longitude, latitude, toler=5.):
        """
        Find the extinction reference file provided by PypeIt that is closest to
        the provided coordinates.

        Parameters
        ----------
        longitude : float
            Geocentric longitude in degrees.
        latitude : float
            Geocentric latitude in degrees.
        toler : float, optional
            Tolerance for matching detector to site in degrees
        """
        # Observation coordinates
        obs_coord = coordinates.SkyCoord(longitude, latitude, frame='gcrs', unit=units.deg)
        # Read list
        extinct_summ = dataPaths.extinction.get_file_path('extinction_curves.txt')
        extinct_files = table.Table.read(extinct_summ, comment='#', format='ascii')
        # Coords
        ext_coord = coordinates.SkyCoord(extinct_files['Lon'], extinct_files['Lat'], frame='gcrs',
                                         unit=units.deg)
        # Match
        idx, d2d, _ = coordinates.match_coordinates_sky(obs_coord, ext_coord, nthneighbor=1)
        if d2d < toler * units.deg:
            return extinct_files[int(idx)]['File']

        # Crash with a helpful error message
        raise PypeItError(
            f'No atmospheric extinction file was found within {toler} degrees of observation at '
            f'lon = {longitude:.1f} lat = {latitude:.1f}.'
        )

    @classmethod
    def from_coordinates(cls, longitude, latitude, toler=5.):
        """
        Instantiate the class using the most appropriate extinction curve given
        the available curves provided by PypeIt.

        Parameters
        ----------
        longitude : float
            Geocentric longitude in degrees.
        latitude : float
            Geocentric latitude in degrees.
        toler : float, optional
            Tolerance for matching detector to site in degrees
        """
        try:
            extinct_file = cls.closest_extinction_file(longitude, latitude, toler=toler)
        except PypeItError as e:
            raise PypeItError(
                f'{e}  You may select a specific extinction file (e.g., KPNO) for use by adding '
                'an ``extinct_file`` to your pypeit_sensfunc or pypeit_fluxcalib input file.  '
                'See instructions at'
                'https://pypeit.readthedocs.io/en/latest/fluxing.html#extinction-correction.'
            )

        log.info(f'Using {extinct_file} for extinction corrections.')
        return cls.from_file(extinct_file)

    @classmethod
    def from_file(cls, extinct_file):
        """
        Load an extinction curve from a file.

        Parameters
        ----------
        extinct_file : :obj:`str`
            Name of a local file or a file distributed by PypeIt.
        """
        _file = dataPaths.extinction.get_file_path(extinct_file)
        data = table.Table.read(_file, comment='#', format='ascii', names=('iwave', 'mag_ext'))
        return cls(data['iwave'], data['mag_ext'], file=extinct_file)

    def correction_factor(self, wave, airmass=1.):
        """
        Return the multiplicative correction factor to apply to fluxes to remove
        atmospheric extinction.

        .. warning::

            Spectral regions outside of the bounds of the atmospheric extinction
            curve are set to the nearest value.

        Parameters
        ----------
        wave : array-like
            Vacuum wavelengths of the *observed* spectrum.
        airmass : float, optional
            Airmass of the observation

        Returns
        -------
        `numpy.ndarray`_
            The correction factor at each wavelength.  Shape matches ``wave``.
        """
        # Warn if extrapolation is necessary
        if np.amin(wave) < np.amin(self.wave) or np.amax(wave) > np.amax(self.wave):
            log.warning(
                'Spectral regions outside of the bounds of the atmospheric extinction curve are '
                'set to the nearest value.'
            )
        # Setup the interpolator
        _mag_ext = interpolate.interp1d(
            self.wave, self.mag_ext, bounds_error=False,
            fill_value=(self.mag_ext[0],self.mag_ext[-1])
        )
        return 10**(0.4 * _mag_ext(wave) * airmass)

    @staticmethod 
    def correct(flux, factor, ivar=None):
        """
        Correct a spectrum for atmospheric extinction.

        Parameters
        ----------
        flux : array-like
            Flux values.
        factor : array-like
            The correction factor to apply: correct flux = flux * factor.  Shape
            must match ``flux``; see :func:`correction_factor`.
        ivar : array-like, optional
            Inverse variance in the flux.  Shape must match ``flux``.  If None,
            uncertainties are not returned.

        Returns
        -------
        corrected_flux : `numpy.ndarray`_
            The corrected fluxes
        corrected_ivar : `numpy.ndarray`_
            Inverse variances of the corrected flux.  If ``ivar`` is None, this
            is not returned.
        """
        # NOTE: These `asarray` statements do not always copy the input
        _flux = np.asarray(flux)
        _factor = np.asarray(factor)
        if _flux.size != _factor.size:
            raise PypeItError('Flux and correction factor arrays must have the same size.')

        if ivar is None:
            return _flux * _factor

        _ivar = np.asarray(ivar)
        if _ivar.size != _flux.size:
            raise PypeItError('Inverse variance and flux arrays must have the same size.')

        return _flux * _factor, _ivar * utils.inverse(_factor**2)

