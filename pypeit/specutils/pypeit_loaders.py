# -*- coding: utf-8 -*-
#
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 04-Feb-2022
#
#  @author: tbowers

""" Specutils loader for PypeIt spec1d data

This is a custom loader for importing PypeIt spec1d data files into specutils
for analysis.  The result is either a Spectrum1D object (for one extracted
object in a given spec1d file) or a SpectrumList (containing all extracted
objects in the spec1d file).

Until this loader is incorportated into Specutils proper, it may be used by
copying it into the user's Specutils cache (nominally ~/.specutils/).

Version History:
    2022-02-04: Initial Version
    2022-09-16: Correct an import error and add module docstring
"""

import astropy.io.fits
import astropy.nddata
import astropy.table
import astropy.units

from specutils.io.parsing_utils import read_fileobj_or_hdulist
from specutils.io.registers import data_loader
from specutils import Spectrum1D, SpectrumList

from pypeit import specobjs



# Identifier Functions =======================================================#
# Identify a PypeIt file in the most general sense
def _identify_pypeit(*args, **kwargs):
    """
    Check if a file is a PypeIt output file.
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdu:
        # Check for header keywords that should be unique to pypeit
        print(f'is pypeit: {"VERSPYP" in hdu[0].header and "PYPELINE" in hdu[0].header}')
        return 'VERSPYP' in hdu[0].header and 'PYPELINE' in hdu[0].header


# Define the base PypeIt identifier based on FITS header information.
def _identify_pypeit_spec1d(*args, **kwargs):
    """
    Test if the file is a PypeIt spec1d file.
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdu:
        return _identify_pypeit(*args, **kwargs) \
                    and hdu[0].header.get('DMODCLS') == 'SpecObjs' \
                    and hdu[1].header.get('DMODCLS') == 'SpecObj'


def identify_pypeit_multislit(origin, *args, **kwargs):
    """
    Check whether the given file is a PypeIt spec1d from the MultiSlit pipeline
    with a single extracted spectrum.
    """
    is_pypeit = _identify_pypeit_spec1d(*args, **kwargs)
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return (
            is_pypeit
            and hdulist[0].header.get("PYPELINE") == "MultiSlit"
            and hdulist[0].header.get("NSPEC") == 1
        )


def identify_pypeit_multislit_list(origin, *args, **kwargs):
    """
    Check whether the given file is a PypeIt spec1d from the MultiSlit pipeline
    with multiple extracted spectra.
    """
    is_pypeit = _identify_pypeit_spec1d(*args, **kwargs)
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return (
            is_pypeit
            and hdulist[0].header.get("PYPELINE") == "MultiSlit"
            and hdulist[0].header.get("NSPEC") > 1
        )


def identify_pypeit_echelle(origin, *args, **kwargs):
    """
    Check whether the given file is a PypeIt spec1d from the Echelle pipeline.
    NOTE: Loader functionality not yet implemented.
    """
    is_pypeit = _identify_pypeit_spec1d(*args, **kwargs)
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return is_pypeit and hdulist[0].header.get("PYPELINE") == "Echelle"


def identify_pypeit_ifu(origin, *args, **kwargs):
    """
    Check whether the given file is a PypeIt spec1d from the IFU pipeline.
    NOTE: Loader functionality not yet implemented.
    """
    is_pypeit = _identify_pypeit_spec1d(*args, **kwargs)
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return is_pypeit and hdulist[0].header.get("PYPELINE") == "IFU"


# Loader Functions ===========================================================#
@data_loader(
    "PypeIt MultiSlit",
    identifier=identify_pypeit_multislit,
    extensions=["fits"],
    priority=10,
    dtype=Spectrum1D,
)
def pypeit_multislit_single_loader(file_name, **kwargs):
    """
    Loader for PypeIt MultiSlit single-spectrum files.

    Parameters
    ----------
    filename : str
        The path to the FITS file

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file.
    """
    spectrum_list = _pypeit_multislit_loader(file_name, **kwargs)
    if len(spectrum_list) == 1:
        return spectrum_list[0]
    if len(spectrum_list) > 1:
        raise RuntimeError(
            f"Input data has {len(spectrum_list)} spectra. "
            "Use SpectrumList.read() instead."
        )
    raise RuntimeError(f"Input data has {len(spectrum_list)} spectra.")


@data_loader(
    "PypeIt MultiSlit list",
    identifier=identify_pypeit_multislit_list,
    extensions=["fits"],
    priority=10,
    dtype=SpectrumList,
)
def pypeit_multislit_multi_loader(filename, **kwargs):
    """
    Loader for PypeIt MultiSlit multiple-spectra files.

    Parameters
    ----------
    filename : str
        The path to the FITS file

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file.
    """
    return _pypeit_multislit_loader(filename, **kwargs)


@data_loader(
    "PypeIt Echelle",
    identifier=identify_pypeit_echelle,
    extensions=["fits"],
    priority=10,
)
def pypeit_echelle_loader(file_name, **kwargs):
    """
    Loader for PypeIt Echelle spectra files.
    """
    raise NotImplementedError("PypeIt Echelle loading not yet implemented.")


@data_loader(
    "PypeIt IFU", identifier=identify_pypeit_ifu, extensions=["fits"], priority=10
)
def pypeit_ifu_loader(file_name, **kwargs):
    """
    Loader for PypeIt IFU spectra files.
    """
    raise NotImplementedError("PypeIt IFU loading not yet implemented.")


# Internal loading functions =================================================#
def _pypeit_multislit_loader(file_name, **kwargs):
    """
    Do the heavy lifting for MultiSlit spectra, to be returned as a
    SpectrumList().
    """
    with astropy.io.fits.open(file_name, **kwargs) as hdulist:

        # Loop through the HDUs looking for spectra
        spectra = []
        for hdu in hdulist:

            # Skip non-spectral HDUs
            # All PypeIt spectra HDUs have EXTNAME starting with 'SPAT'
            if "EXTNAME" not in hdu.header or hdu.header["EXTNAME"][:4] != "SPAT":
                continue

            # Read in this BinTable
            spec_obj = astropy.table.Table.read(hdu)

            # Combine the primary header with the header for this BinTable
            meta = {"header": hdulist[0].header + hdu.header}

            # Check for Optimal extraction... else use Boxcar
            extract_type = "OPT" if "OPT_COUNTS" in spec_obj.colnames else "BOX"
            # Check for fluxed spectrum, else use counts
            if f"{extract_type}_FLAM" in spec_obj.colnames:
                flux_type = "FLAM"
                flux_unit = astropy.units.Unit("1e-17 erg/(s cm^2 Angstrom)")
            else:
                flux_type = "COUNTS"
                flux_unit = astropy.units.Unit("electron")

            data = astropy.units.Quantity(
                spec_obj[f"{extract_type}_{flux_type}"] * flux_unit
            )
            uncert = astropy.nddata.InverseVariance(
                spec_obj[f"{extract_type}_{flux_type}_IVAR"] / flux_unit**2
            )

            # Wavelength
            wavl = astropy.units.Quantity(
                spec_obj[f"{extract_type}_WAVE"] * astropy.units.Unit("angstrom")
            )

            # Package the spectrum as a Spectrum1d() object
            spec = Spectrum1D(
                flux=data,
                uncertainty=uncert,
                meta=meta,
                spectral_axis=wavl,
                velocity_convention="doppler_optical",
                bin_specification="centers",
            )
            spectra.append(spec)

        # Package and return
        return SpectrumList(spectra)

