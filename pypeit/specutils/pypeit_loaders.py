"""
.. include:: ../include/links.rst

Data parsers for use with `specutils`_.

.. include:: ../include/specutils_usage.rst

Version History
---------------

- 2022-02-04: Initial Version (tbowers)
- 2022-09-16: Correct an import error and add module docstring (tbowers)
- 2023-03-09: Moved into the main pypeit repo and refactored (KBW)
- 2023-06-23: Added sensible error message for incorrect spec1d loading (tbowers)

"""

from IPython import embed

import numpy as np

import astropy.io.fits
import astropy.nddata
import astropy.units

try:
    from specutils.io.parsing_utils import read_fileobj_or_hdulist
    from specutils.io.registers import data_loader
    from specutils import Spectrum1D, SpectrumList
except ModuleNotFoundError:
    raise ModuleNotFoundError('Unable to import specutils.  Install pypeit with the specutils '
                              'option to use the pypeit.specutils module.')

from pypeit import __version__
from pypeit.pypmsgs import PypeItError
from pypeit import msgs
from pypeit import specobjs
from pypeit import onespec
from pypeit import utils


def _enforce_monotonic_wavelengths(wave, flux, ivar, strict=True):
    """
    Force the spectrum to have a monotonically increasing wavelength vector.

    Parameters
    ----------
    wave : `numpy.ndarray`_
        Spectrum wavelengths
    flux : `numpy.ndarray`_
        Spectrum flux
    ivar : `numpy.ndarray`_
        Spectrum inverse variance. Can be None or the standard deviation.  The
        only operation on this and the ``flux`` vector is to downselect the
        monotonically increasing values.
    strict : bool, optional
        Check that the wavelength vector is monotonically increasing.  If not,
        raise an error (as would be done by the `specutils.SpectrumList`_ class).
        If False, wavelengths that are *not* monotonically increasing are masked
        in the construction of the returned `specutils.SpectrumList`_ object.

    Returns
    -------
    wave : `numpy.ndarray`_
        Edited wavelength vector.  This may be an unchanged reference to the
        original vector.
    flux : `numpy.ndarray`_
        Edited flux vector.  This may be an unchanged reference to the original
        vector.
    ivar : `numpy.ndarray`_
        Edited inverse variance vector.  This may be an unchanged reference to
        the original vector.
    """
    indx = np.diff(wave) > 0
    if np.all(indx):
        # Wavelengths are monotonic, so we're done
        return wave, flux, ivar

    if strict:
        # Wavelengths are not monotonic, but the user expects them to be, so
        # fault.
        msgs.error('Wavelengths are not monotonically increasing!  Circumvent this fault by '
                   'setting strict=False, but BEWARE that this is likely the result of an '
                   'error in the data reduction!')

    # Wavelengths are not monotonic, but the user wants to keep going.
    msgs.warn('Wavelengths are not monotonically increasing!  Strict was set to False, so '
              'measurements after a negative step in wavelength are removed from the constructed '
              'spectrum.  BEWARE that this is likely the result of an error in the data '
              'reduction!')

    # NOTE: This is the brute force approach.  If this becomes something that we
    # want to be more acceptable, we should consider instead fitting a low-order
    # polynomial to the pixel vs. wavelength function and rejecting strong
    # outliers.
    pix = np.arange(wave.size)
    indx = np.append([True], indx)
    while not np.all(indx):
        pix = pix[indx]
        wave = wave[indx]
        indx = np.append([True], np.diff(wave) > 0)

    return wave, flux[pix], None if ivar is None else ivar[pix]


# Identifier Functions =======================================================#
def _identify_pypeit(*args, **kwargs):
    """
    Check if a file is a PypeIt output file, in the most general sense.

    This currently only checks if ``VERSPYP`` is in the primary header.
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdu:
        # Check for header keywords that should be unique to pypeit
        return 'VERSPYP' in hdu[0].header


def identify_pypeit_spec1d(origin, *args, **kwargs):
    """
    Check if a file is a PypeIt spec1d file.

    In addition to checking if it is a PypeIt file (see
    :func:`_identify_pypeit`), this also checks that the datamodel classes are
    correct.
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdu:
        return _identify_pypeit(*args, **kwargs) \
                    and hdu[0].header.get('DMODCLS') == 'SpecObjs' \
                    and hdu[1].header.get('DMODCLS') == 'SpecObj'


def identify_pypeit_onespec(origin, *args, **kwargs):
    """
    Check if a file is a PypeIt file that follows the
    :class:`pypeit.onespec.OneSpec` datamodel.

    In addition to checking if it is a PypeIt file (see
    :func:`_identify_pypeit`), this also checks that the datamodel classes are
    correct.
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdu:
        return _identify_pypeit(*args, **kwargs) \
                    and hdu[1].header.get('DMODCLS') == 'OneSpec'


# Loader Functions ===========================================================#
@data_loader('PypeIt spec1d',
             identifier=identify_pypeit_spec1d,
             extensions=["fits"],
             priority=10,
             dtype=SpectrumList,
             autogenerate_spectrumlist=False)
def pypeit_spec1d_loader(filename, extract=None, fluxed=True, strict=True, chk_version=True,
                         **kwargs):
    """
    Load spectra from a PypeIt spec1d file into a SpectrumList.

    Parameters
    ----------
    filename : str
        The path to the FITS file
    extract : str, optional
        The extraction used to produce the spectrum.  Must be either None,
        ``'BOX'`` (for a boxcar extraction), or ``'OPT'`` for optimal
        extraction.  If None, the optimal extraction will be returned, if it
        exists, otherwise the boxcar extraction will be returned.
    fluxed : bool, optional
        If True, return the flux-calibrated spectrum, if it exists.  If the flux
        calibration hasn't been performed or ``fluxed=False``, the spectrum is
        returned in counts.
    strict : bool, optional
        Check that the wavelength vector is monotonically increasing.  If not,
        raise an error (as would be done by the `specutils.SpectrumList`_ class).
        If False, wavelengths that are *not* monotonically increasing are masked
        in the construction of the returned `specutils.SpectrumList`_ object.
    chk_version : :obj:`bool`, optional
        When reading in existing files written by PypeIt, perform strict version
        checking to ensure a valid file.  If False, the code will try to keep
        going, but this may lead to faults and quiet failures.  User beware!
    kwargs : dict, optional
        **Ignored!**  Used to catch spurious arguments passed to the base class
        that are ignored by this function.

    Returns
    -------
    spec : `specutils.SpectrumList`_
        Contains all spectra in the PypeIt spec1d file
    """
    # Try to load the file and ignoring any version mismatch
    try:
        sobjs = specobjs.SpecObjs.from_fitsfile(filename, chk_version=chk_version)
    except PypeItError:
        file_pypeit_version = astropy.io.fits.getval(filename, 'VERSPYP', 'PRIMARY')
        msgs.error(f'Unable to ingest {filename.name} using pypeit.specobjs module from your version '
                   f'of PypeIt ({__version__}).  The version used to write the file is '
                   f'{file_pypeit_version}.  If these are different, you may need to re-reduce '
                   'your data using your current PypeIt version or install the matching version '
                   'of PypeIt (e.g., pip install pypeit==1.11.0).')

    spec = []
    for sobj in sobjs:
        # Check that the file has the requested data
        _ext, _cal = sobj.best_ext_match(extract=extract, fluxed=fluxed)
        _wave, _flux, _ivar, _gpm = sobj.get_box_ext(fluxed=_cal) if _ext == 'BOX' \
                                        else sobj.get_opt_ext(fluxed=_cal)
        if not np.all(_gpm):
            msgs.warn(f'Ignoring {np.sum(np.logical_not(_gpm))} masked pixels.')
        if not np.any(_gpm):
            msgs.warn(f'Spectrum {sobj.NAME} is fully masked and will be ignored!')
            continue
        _wave, _flux, _ivar = _enforce_monotonic_wavelengths(_wave[_gpm], _flux[_gpm], _ivar[_gpm],
                                                             strict=strict)
        _sigma = np.sqrt(utils.inverse(_ivar))
        flux_unit = astropy.units.Unit("1e-17 erg/(s cm^2 Angstrom)" if _cal else "electron")
        spec += \
            [Spectrum1D(flux=astropy.units.Quantity(_flux * flux_unit),
                        uncertainty=astropy.nddata.StdDevUncertainty(_sigma * flux_unit),
                        meta={'name': sobj.NAME, 'extract': _ext, 'fluxed': _cal},
                        spectral_axis=astropy.units.Quantity(_wave * astropy.units.angstrom),
                        velocity_convention="doppler_optical",
                        bin_specification="centers")]
    return SpectrumList(spec)


@data_loader('PypeIt onespec',
             identifier=identify_pypeit_onespec,
             extensions=["fits"],
             priority=10,
             dtype=Spectrum1D)
def pypeit_onespec_loader(filename, grid=False, strict=True, chk_version=True, **kwargs):
    """
    Load a spectrum from a PypeIt OneSpec file into a Spectrum1D object.

    Parameters
    ----------
    filename : str
        The path to the FITS file
    grid : bool, optional
        Use the uniform grid wavelengths, instead of the contribution-weighted
        center.
    strict : bool, optional
        Check that the wavelength vector is monotonically increasing.  If not,
        raise an error (as would be done by the `specutils.Spectrum1D`_ class).
        If False, wavelengths that are *not* monotonically increasing are masked
        in the construction of the returned `specutils.Spectrum1D`_ object.
    chk_version : :obj:`bool`, optional
        When reading in existing files written by PypeIt, perform strict version
        checking to ensure a valid file.  If False, the code will try to keep
        going, but this may lead to faults and quiet failures.  User beware!
    kwargs : dict, optional
        **Ignored!**  Used to catch spurious arguments passed to the base class
        that are ignored by this function.

    Returns
    -------
    spec : `specutils.Spectrum1D`_
        Spectrum in the PypeIt OneSpec file
    """
    # Try to load the file and ignoring any version mismatch
    try:
        spec = onespec.OneSpec.from_file(filename, chk_version=chk_version)
    except PypeItError:
        file_pypeit_version = astropy.io.fits.getval(filename, 'VERSPYP', 'PRIMARY')
        msgs.error(f'Unable to ingest {filename.name} using pypeit.specobjs module from your version '
                   f'of PypeIt ({__version__}).  The version used to write the file is '
                   f'{file_pypeit_version}.  If these are different, you may need to re-reduce '
                   'your data using your current PypeIt version or install the matching version '
                   'of PypeIt (e.g., pip install pypeit==1.11.0).')

    flux_unit = astropy.units.Unit("1e-17 erg/(s cm^2 Angstrom)" if spec.fluxed else "ct/s")
    wave = spec.wave_grid_mid if grid else spec.wave
    wave, flux, sigma = _enforce_monotonic_wavelengths(wave, spec.flux, spec.sigma, strict=strict)

    # If the input filename is actually a string, assign it as the spectrum
    # name.  Otherwise, try assuming it's a _io.FileIO object, and if that
    # doesn't work assign an empty string as the name.
    if isinstance(filename, str):
        name = filename
    else:
        try:
            name = filename.name    # Needed for _io.FileIO objects
        except AttributeError:
            name = ''

    # TODO We should be dealing with masking here.
    return Spectrum1D(flux=astropy.units.Quantity(flux * flux_unit),
                      uncertainty=None if spec.sigma is None 
                                  else astropy.nddata.StdDevUncertainty(sigma * flux_unit),
                      meta={'name': name, 'extract': spec.ext_mode, 'fluxed': spec.fluxed,
                            'grid': grid},
                      spectral_axis=astropy.units.Quantity(wave * astropy.units.angstrom),
                      velocity_convention="doppler_optical",
                      bin_specification="centers")

# Warning Function ===========================================================#
@data_loader('PypeIt spec1d nolist',
             identifier=identify_pypeit_spec1d,
             extensions=["fits"],
             priority=10,
             dtype=Spectrum1D,
             autogenerate_spectrumlist=False)
def pypeit_spec1d_loader_nolist(filename, extract=None, fluxed=True, **kwargs):
    """
    Sensible error message if a user tries to load spectra from a PypeIt spec1d
    file into a Spectrum1D.

    This is not allowed because spec1d files may contain mutliple spectra.  This
    function accepts all arguments as the SpectrumList version, but only outputs
    a PypeIt Error with a sensible message.

    This avoids receiving unhelpful error messages such as::

        OSError: Could not identify column containing the wavelength, frequency or energy

    Parameters
    ----------
    filename : str
        The path to the FITS file
    extract : str, optional
        The extraction used to produce the spectrum.  Must be either None,
        ``'BOX'`` (for a boxcar extraction), or ``'OPT'`` for optimal
        extraction.  If None, the optimal extraction will be returned, if it
        exists, otherwise the boxcar extraction will be returned.
    fluxed : bool, optional
        If True, return the flux-calibrated spectrum, if it exists.  If the flux
        calibration hasn't been performed or ``fluxed=False``, the spectrum is
        returned in counts.
    """
    msgs.error(f'The spec1d file {filename.name} cannot be ingested into a Spectrum1D object.'
               f'{msgs.newline()}Please use the SpectrumList object for spec1d files.')
