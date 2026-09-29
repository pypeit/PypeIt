"""
Quicklook viewer configuration for Keck/MOSFIRE.
"""

from __future__ import annotations

import os
from typing import Dict

import numpy as np
from astropy.io import fits

from pypeit.images import buildimage
from pypeit.spectrographs.util import load_spectrograph

from .base import Instrument


class MOSFIRE(Instrument):
    instrume_value = "MOSFIRE"
    detector_prefix = "DET"

    def __init__(self, logger) -> None:
        """Initialise the MOSFIRE instrument with Keck-MOSFIRE–specific column definitions.

        Defines raw columns including a ``DITHER_POS`` column decoded from the
        ``PATTERN``/``FRAMEID`` FITS headers, and sets MOSFIRE-specific reduced
        columns with CSU mask, filter, and dispname fields.

        Parameters
        ----------
        logger : logging.Logger
            Ginga application logger.
        """
        super().__init__(logger)
        self.pypeit_name = "keck_mosfire"
        self.columns["raw"] = [
            ("Type", "icon"),
            ("Frame No", "FRAMENO"),
            ("Name", "name"),
            ("Dither Pos", "DITHER_POS"),
            ("Object", "OBJECT"),
            ("Img Type", "IMTYPE"),
            ("Mask Name", "MASKNAME"),
            ("Obs Mode", "OBSMODE"),
            ("Exp Time", "EXPTIME"),
            ("Last Changed", "st_mtime_str"),
        ]
        self.columns["reduced"] = [
            ("Type", "icon"),
            ("Name", "name"),
            ("CSU Mask", "MASKNAME"),
            ("Filter", "FILTER1"),
            ("Dispname", "FILTER2"),
            ("Slit Width", "SLITWIDTH"),
            ("Last Changed", "st_mtime_str"),
        ]

    def get_raw_info(self, path: str) -> Dict[str, object]:
        """Read display metadata from a MOSFIRE raw FITS file.

        Overrides :meth:`Instrument.get_raw_info` to add the ``DITHER_POS``
        field: ``"N/A"`` when the pattern is ``"Stare"``, otherwise the
        ``FRAMEID`` value (e.g. ``"A"``, ``"B"``).

        Parameters
        ----------
        path : str
            Absolute path to a MOSFIRE raw FITS file.

        Returns
        -------
        dict
            All keys from :meth:`Instrument._read_header_fields` plus
            ``DITHER_POS``.
        """
        with fits.open(path) as hdul:
            hdr = hdul[0].header
            info = self._read_header_fields(hdr)
            pattern = hdr.get("PATTERN", "")
            info["DITHER_POS"] = "N/A" if pattern == "Stare" else hdr.get("FRAMEID", "N/A")
            return info

    def get_display_image(self, raw_path: str) -> np.ndarray:
        """Build an overscan-subtracted, oriented display image.

        Uses PypeIt's standard ``buildimage_fromlist`` pipeline with
        ``biasframe`` processing parameters (overscan subtraction, trimming,
        orientation — no dark or flat calibration).
        """
        spec = load_spectrograph("keck_mosfire")
        par = spec.default_pypeit_par()['calibrations']['biasframe']
        img = buildimage.buildimage_fromlist(spec, 1, par, [raw_path], mosaic=False)
        return img.image

    def get_reduced_info(self, path: str) -> Dict[str, object]:
        """Read display metadata from a MOSFIRE calibration directory or reduced file.

        Overrides :meth:`Instrument.get_reduced_info`.  For calibration
        directories the ``.pypeit`` setup file is parsed via
        :meth:`_read_pypeit_setup_config`; for FITS files the primary header
        is used.

        Parameters
        ----------
        path : str
            Absolute path to a calibration directory (e.g. ``keck_mosfire_A/``)
            or a reduced FITS file.

        Returns
        -------
        dict
            Keys: ``MASKNAME`` (CSU mask, PypeIt ``decker_secondary``),
            ``FILTER1`` (bandpass filter, PypeIt ``filter1``),
            ``FILTER2`` (dispname/order-blocking), ``SLITWIDTH``
            (slit width, PypeIt ``slitwid``).
        """
        if os.path.isdir(path):
            # MOSFIRE configuration_keys: decker_secondary (CSU mask name),
            # slitlength, slitwid, dispname, filter1.
            cfg = self._read_pypeit_setup_config(path)
            return {
                # decker_secondary is the CSU mask name (e.g. "GS37", "LONGSLIT_46x0.7")
                "MASKNAME": cfg.get("decker_secondary", "N/A"),
                # filter1 is the primary bandpass filter (e.g. "K", "H", "J")
                "FILTER1": cfg.get("filter1", "N/A"),
                # dispname is the grating/order-blocking setting; often same as filter1
                "FILTER2": cfg.get("dispname", "N/A"),
                "SLITWIDTH": cfg.get("slitwid", "N/A"),
            }
        try:
            with fits.open(path) as hdul:
                hdr = hdul[0].header
                filter1 = hdr.get("FILTER", hdr.get("FILTER1", "N/A"))
                filter2 = hdr.get("FILTER2", "N/A")
                return {
                    "MASKNAME": hdr.get("MASKNAME", "N/A"),
                    "FILTER1": filter1,
                    "FILTER2": filter2,
                    "SLITWIDTH": hdr.get("MGTNAME", hdr.get("SLIT", hdr.get("SLITWIDTH", "N/A"))),
                }
        except Exception:
            return {}
