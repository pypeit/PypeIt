"""
Quicklook viewer configuration for the Keck/LRIS blue and red arms.
"""

from __future__ import annotations

import os
from typing import Dict

import numpy as np
from astropy.io import fits

from pypeit.images import buildimage
from pypeit.spectrographs.util import load_spectrograph

from .base import Instrument


class LRISBlue(Instrument):
    """Keck LRIS Blue channel — multi-slit, 2-detector mosaic.
    
    Untested!
    """

    instrume_value = "LRISBLUE"
    detector_prefix = "MSC"

    def __init__(self, logger) -> None:
        """Initialise the LRIS Blue instrument with Keck-LRIS–Blue–specific column definitions.

        Defines raw columns for grism (``GRISNAME``) and dichroic
        (``DICHNAME``), and sets LRIS Blue–specific reduced columns with
        slit/mask, grating/grism, and dichroic fields.

        Parameters
        ----------
        logger : logging.Logger
            Ginga application logger.
        """
        super().__init__(logger)
        self.pypeit_name = "keck_lris_blue"
        self.columns["raw"] = [
            ("Type", "icon"),
            ("Frame No", "FRAMENO"),
            ("Name", "name"),
            ("Object", "OBJECT"),
            ("Img Type", "IMTYPE"),
            ("Slit/Mask", "MASKNAME"),
            ("Grism", "GRISNAME"),
            ("Dichroic", "DICHNAME"),
            ("Exp Time", "EXPTIME"),
            ("Last Changed", "st_mtime_str"),
        ]
        self.columns["reduced"] = [
            ("Type", "icon"),
            ("Name", "name"),
            ("Slit/Mask", "SLITNAME"),
            ("Grating/Grism", "DISPNAME"),
            ("Dichroic", "DICHNAME"),
            ("Last Changed", "st_mtime_str"),
        ]

    def get_raw_info(self, path: str) -> Dict[str, object]:
        """Read display metadata from an LRIS Blue raw FITS file.

        Overrides :meth:`Instrument.get_raw_info` to populate LRIS Blue
        column keys: slit/mask (``SLITNAME``), grism (``GRISNAME``), and
        dichroic (``DICHNAME``).

        Parameters
        ----------
        path : str
            Absolute path to an LRIS Blue raw FITS file.

        Returns
        -------
        dict
            Keys: ``OBJECT`` (``TARGNAME``), ``IMTYPE`` (``KOAIMTYP``),
            ``MASKNAME`` (``SLITNAME``), ``GRISNAME``, ``DICHNAME``,
            plus all keys from :meth:`Instrument._read_header_fields`.
        """
        with fits.open(path) as hdul:
            hdr = hdul[0].header
            info = self._read_header_fields(hdr)
            info["OBJECT"] = hdr.get("TARGNAME", "N/A")
            info["IMTYPE"] = hdr.get("KOAIMTYP", "N/A")
            info["MASKNAME"] = hdr.get("SLITNAME", "N/A")
            info["GRISNAME"] = hdr.get("GRISNAME", "N/A")
            info["DICHNAME"] = hdr.get("DICHNAME", "N/A")
            return info

    def get_display_image(self, raw_path: str) -> np.ndarray:
        """Build an overscan-subtracted, oriented display image for LRIS Blue.

        Delegates to PypeIt's ``buildimage_fromlist`` with ``biasframe``
        processing parameters applied to the 2-detector mosaic.

        Parameters
        ----------
        raw_path : str
            Absolute path to an LRIS Blue raw FITS file.

        Returns
        -------
        numpy.ndarray
            Processed 2-D image array.
        """
        spec = load_spectrograph("keck_lris_blue")
        par = spec.default_pypeit_par()["calibrations"]["biasframe"]
        img = buildimage.buildimage_fromlist(spec, 1, par, [raw_path], mosaic=False)
        return img.image

    def get_reduced_info(self, path: str) -> Dict[str, object]:
        """Read display metadata from an LRIS Blue calibration directory or reduced file.

        Overrides :meth:`Instrument.get_reduced_info`.  For calibration
        directories the ``.pypeit`` setup file is parsed.

        Parameters
        ----------
        path : str
            Absolute path to a calibration directory (e.g. ``keck_lris_blue_A/``)
            or a reduced FITS file.

        Returns
        -------
        dict
            Keys: ``SLITNAME`` (PypeIt ``decker``), ``DISPNAME`` (grism,
            PypeIt ``dispname``), ``DICHNAME`` (dichroic, PypeIt
            ``dichroic``).
        """
        if os.path.isdir(path):
            cfg = self._read_pypeit_setup_config(path)
            return {
                "SLITNAME": cfg.get("decker", "N/A"),
                "DISPNAME": cfg.get("dispname", "N/A"),
                "DICHNAME": cfg.get("dichroic", "N/A"),
            }
        try:
            with fits.open(path) as hdul:
                hdr = hdul[0].header
                return {
                    "SLITNAME": hdr.get("SLITNAME", "N/A"),
                    "DISPNAME": hdr.get("GRISNAME", "N/A"),
                    "DICHNAME": hdr.get("DICHNAME", "N/A"),
                }
        except Exception:
            return {}


class LRISRed(Instrument):
    """Keck LRIS Red channel (Mark4 detector) — multi-slit, single detector.
    
    Untested!"""

    instrume_value = "LRIS"
    detector_prefix = "DET" #

    def __init__(self, logger) -> None:
        """Initialise the LRIS Red instrument with Keck-LRIS–Red–specific column definitions.

        Defines raw columns for grating (``GRANAME``) and dichroic
        (``DICHNAME``), and sets LRIS Red–specific reduced columns with
        slit/mask, grating/grism, and dichroic fields.

        Parameters
        ----------
        logger : logging.Logger
            Ginga application logger.
        """
        super().__init__(logger)
        self.pypeit_name = "keck_lris_red_mark4"
        self.columns["raw"] = [
            ("Type", "icon"),
            ("Frame No", "FRAMENO"),
            ("Name", "name"),
            ("Object", "OBJECT"),
            ("Img Type", "IMTYPE"),
            ("Slit/Mask", "MASKNAME"),
            ("Grating", "GRANAME"),
            ("Dichroic", "DICHNAME"),
            ("Exp Time", "EXPTIME"),
            ("Last Changed", "st_mtime_str"),
        ]
        self.columns["reduced"] = [
            ("Type", "icon"),
            ("Name", "name"),
            ("Slit/Mask", "SLITNAME"),
            ("Grating/Grism", "DISPNAME"),
            ("Dichroic", "DICHNAME"),
            ("Last Changed", "st_mtime_str"),
        ]

    def get_raw_info(self, path: str) -> Dict[str, object]:
        """Read display metadata from an LRIS Red raw FITS file.

        Overrides :meth:`Instrument.get_raw_info` to populate LRIS Red column
        keys: slit/mask (``SLITNAME``), grating (``GRANAME``), and dichroic
        (``DICHNAME``).

        Parameters
        ----------
        path : str
            Absolute path to an LRIS Red raw FITS file.

        Returns
        -------
        dict
            Keys: ``OBJECT`` (``TARGNAME``), ``IMTYPE`` (``KOAIMTYP``),
            ``MASKNAME`` (``SLITNAME``), ``GRANAME``, ``DICHNAME``,
            plus all keys from :meth:`Instrument._read_header_fields`.
        """
        with fits.open(path) as hdul:
            hdr = hdul[0].header
            info = self._read_header_fields(hdr)
            info["OBJECT"] = hdr.get("TARGNAME", "N/A")
            info["IMTYPE"] = hdr.get("KOAIMTYP", "N/A")
            info["MASKNAME"] = hdr.get("SLITNAME", "N/A")
            info["GRANAME"] = hdr.get("GRANAME", "N/A")
            info["DICHNAME"] = hdr.get("DICHNAME", "N/A")
            return info

    def get_display_image(self, raw_path: str) -> np.ndarray:
        """Build an overscan-subtracted, oriented display image for LRIS Red.

        Delegates to PypeIt's ``buildimage_fromlist`` with ``biasframe``
        processing parameters for the Mark4 single-detector configuration.

        Parameters
        ----------
        raw_path : str
            Absolute path to an LRIS Red raw FITS file.

        Returns
        -------
        numpy.ndarray
            Processed 2-D image array.
        """
        spec = load_spectrograph("keck_lris_red_mark4")
        par = spec.default_pypeit_par()["calibrations"]["biasframe"]
        img = buildimage.buildimage_fromlist(spec, 1, par, [raw_path], mosaic=False)
        return img.image

    def get_reduced_info(self, path: str) -> Dict[str, object]:
        """Read display metadata from an LRIS Red calibration directory or reduced file.

        Overrides :meth:`Instrument.get_reduced_info`.  For calibration
        directories the ``.pypeit`` setup file is parsed.

        Parameters
        ----------
        path : str
            Absolute path to a calibration directory (e.g.
            ``keck_lris_red_mark4_A/``) or a reduced FITS file.

        Returns
        -------
        dict
            Keys: ``SLITNAME`` (PypeIt ``decker``), ``DISPNAME`` (grating,
            PypeIt ``dispname``), ``DICHNAME`` (dichroic, PypeIt
            ``dichroic``).
        """
        if os.path.isdir(path):
            cfg = self._read_pypeit_setup_config(path)
            return {
                "SLITNAME": cfg.get("decker", "N/A"),
                "DISPNAME": cfg.get("dispname", "N/A"),
                "DICHNAME": cfg.get("dichroic", "N/A"),
            }
        try:
            with fits.open(path) as hdul:
                hdr = hdul[0].header
                return {
                    "SLITNAME": hdr.get("SLITNAME", "N/A"),
                    "DISPNAME": hdr.get("GRANAME", "N/A"),
                    "DICHNAME": hdr.get("DICHNAME", "N/A"),
                }
        except Exception:
            return {}
