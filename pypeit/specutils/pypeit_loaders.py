"""
.. include:: ../include/links.rst

Data parsers for use with `specutils`_.

.. include:: ../include/specutils_usage.rst

Version History
---------------

- 2022-02-04: Initial Version (tbowers)
- 2022-09-16: Correct an import error and add module docstring (tbowers)
- 2023-03-09: Moved into the main pypeit repo and refactored (KBW)

"""

from IPython import embed

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
             dtype=SpectrumList)
def pypeit_spec1d_loader(filename, extract=None, fluxed=True, **kwargs):
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

    Returns
    -------
    spec : SpectrumList
        Contains all spectra in the PypeIt spec1d file
    """
    # Try to load the file and ignoring any version mismatch
    try:
        sobjs = specobjs.SpecObjs.from_fitsfile(filename, chk_version=False)
    except PypeItError:
        file_pypeit_version = astropy.io.fits.getval(filename, 'VERSPYP', 'PRIMARY')
        msgs.error(f'Unable to ingest {filename} using pypeit.specobjs module from your version '
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

        flux_unit = astropy.units.Unit("1e-17 erg/(s cm^2 Angstrom)" if _cal else "electron")
        spec += \
            [Spectrum1D(flux=astropy.units.Quantity(_flux * flux_unit),
                        uncertainty=astropy.nddata.InverseVariance(_ivar / flux_unit**2),
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
def pypeit_onespec_loader(filename, grid=False, **kwargs):
    """
    Load a spectrum from a PypeIt OneSpec file into a Spectrum1D object.

    Parameters
    ----------
    filename : str
        The path to the FITS file
    grid : bool, optional
        Use the uniform grid wavelengths, instead of the contribution-weighted
        center.

    Returns
    -------
    spec : Spectrum1D
        Spectrum in the PypeIt OneSpec file
    """
    # Try to load the file and ignoring any version mismatch
    try:
        spec = onespec.OneSpec.from_file(filename)
    except PypeItError:
        file_pypeit_version = astropy.io.fits.getval(filename, 'VERSPYP', 'PRIMARY')
        msgs.error(f'Unable to ingest {filename} using pypeit.specobjs module from your version '
                   f'of PypeIt ({__version__}).  The version used to write the file is '
                   f'{file_pypeit_version}.  If these are different, you may need to re-reduce '
                   'your data using your current PypeIt version or install the matching version '
                   'of PypeIt (e.g., pip install pypeit==1.11.0).')

    flux_unit = astropy.units.Unit("1e-17 erg/(s cm^2 Angstrom)" if spec.fluxed else "ct/s")
    wave = spec.wave_grid_mid if grid else spec.wave
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

    return Spectrum1D(flux=astropy.units.Quantity(spec.flux * flux_unit),
                      uncertainty=None if spec.ivar is None 
                        else astropy.nddata.InverseVariance(spec.ivar / flux_unit**2),
                      meta={'name': name, 'extract': spec.ext_mode, 'fluxed': spec.fluxed,
                            'grid': grid},
                      spectral_axis=astropy.units.Quantity(wave * astropy.units.angstrom),
                      velocity_convention="doppler_optical",
                      bin_specification="centers")



