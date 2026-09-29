"""
Quicklook viewer configuration for Keck/HIRES.
"""

from __future__ import annotations

import os
from typing import Dict

import numpy as np
from astropy.io import fits

from .base import Instrument


class HIRES(Instrument):
    """Keck HIRES — UV/optical echelle, 3-detector mosaic.

    PypeIt marks HIRES as ``supported = False`` so reductions are not
    expected, but the file browser and calibration-directory viewer work
    normally for header inspection.
    """

    instrume_value = "HIRES"
    detector_prefix = "MSC"

    def __init__(self, logger) -> None:
        """Initialise the HIRES instrument with Keck-HIRES–specific column definitions.

        Defines raw columns for decker (``DECKNAME``) and cross-disperser
        (``XDISPERS``), and sets instrument-specific reduced columns.
        Because PypeIt marks HIRES as unsupported, ``get_display_image``
        falls back to a raw pixel read of extension 1.

        Parameters
        ----------
        logger : logging.Logger
            Ginga application logger.
        """
        super().__init__(logger)
        self.pypeit_name = "keck_hires"
        raw_cols = [
            ("Type", "icon"),
            ("Frame No", "FRAMENO"),
            ("Name", "name"),
            ("Object", "OBJECT"),
            ("Img Type", "IMTYPE"),
            ("Decker", "DECKNAME"),
            ("XDisp", "XDISPERS"),
            ("Exp Time", "EXPTIME"),
            ("Last Changed", "st_mtime_str"),
        ]
        self.columns["raw"] = raw_cols
        self.columns["reduced"] = [
            ("Type", "icon"),
            ("Name", "name"),
            ("Decker", "DECKNAME"),
            ("XDisp", "XDISPERS"),
            ("Filter", "FILTER1"),
            ("Last Changed", "st_mtime_str"),
        ]

    def get_raw_info(self, path: str) -> Dict[str, object]:
        """Read display metadata from a HIRES raw FITS file.

        Overrides :meth:`Instrument.get_raw_info` to populate HIRES-specific
        column keys: decker (``DECKNAME``), cross-disperser (``XDISPERS``),
        and elapsed time (``ELAPTIME``).

        Parameters
        ----------
        path : str
            Absolute path to a HIRES raw FITS file.

        Returns
        -------
        dict
            Keys: ``OBJECT`` (``TARGNAME`` → ``OBJECT``), ``IMTYPE``
            (``KOAIMTYP``), ``DECKNAME``, ``XDISPERS``, ``EXPTIME``
            (``ELAPTIME``).
        """
        with fits.open(path) as hdul:
            hdr = hdul[0].header
            info = self._read_header_fields(hdr)
            info["OBJECT"] = hdr.get("TARGNAME", hdr.get("OBJECT", "N/A"))
            info["IMTYPE"] = hdr.get("KOAIMTYP", "N/A")
            info["DECKNAME"] = hdr.get("DECKNAME", "N/A")
            info["XDISPERS"] = hdr.get("XDISPERS", "N/A")
            info["EXPTIME"] = hdr.get("ELAPTIME", "N/A")
            return info

    def get_display_image(self, raw_path: str) -> np.ndarray:
        """Display the first science extension of a HIRES raw file.

        HIRES has three detectors but PypeIt does not currently support it,
        so we fall back to a simple raw-pixel read of extension 1.
        """
        with fits.open(raw_path) as hdul:
            # Extension 0 is the primary (empty); science data start at 1
            data = hdul[1].data
        if data is None:
            return np.zeros((100, 100), dtype=float)
        return data.astype(float)

    def get_reduced_info(self, path: str) -> Dict[str, object]:
        """Read display metadata from a HIRES calibration directory or reduced file.

        Overrides :meth:`Instrument.get_reduced_info`.  For calibration
        directories the ``.pypeit`` setup file is parsed.

        Parameters
        ----------
        path : str
            Absolute path to a calibration directory (e.g.
            ``keck_hires_A/``) or a reduced FITS file.

        Returns
        -------
        dict
            Keys: ``DECKNAME`` (PypeIt ``decker``), ``XDISPERS`` (PypeIt
            ``dispname``), ``FILTER1`` (cross-disperser filter, PypeIt
            ``filter1`` / ``FIL1NAME``).
        """
        if os.path.isdir(path):
            cfg = self._read_pypeit_setup_config(path)
            return {
                "DECKNAME": cfg.get("decker", "N/A"),
                "XDISPERS": cfg.get("dispname", "N/A"),
                "FILTER1": cfg.get("filter1", "N/A"),
            }
        try:
            with fits.open(path) as hdul:
                hdr = hdul[0].header
                return {
                    "DECKNAME": hdr.get("DECKNAME", "N/A"),
                    "XDISPERS": hdr.get("XDISPERS", "N/A"),
                    "FILTER1": hdr.get("FIL1NAME", "N/A"),
                }
        except Exception:
            return {}
