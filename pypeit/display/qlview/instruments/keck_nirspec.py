"""
Quicklook viewer configuration for Keck/NIRSPEC.
"""

from __future__ import annotations

import os
from typing import Dict

import numpy as np
from astropy.io import fits

from pypeit.images import buildimage
from pypeit.spectrographs.util import load_spectrograph

from .base import Instrument


class NIRSPEC(Instrument):
    """Keck NIRSPEC (post-2018 upgrade) — near-IR echelle, single detector.
    
    Untested!"""

    instrume_value = "NIRSPEC"
    detector_prefix = "DET"

    def __init__(self, logger) -> None:
        """Initialise the NIRSPEC instrument with Keck-NIRSPEC–specific column definitions.

        Defines raw columns for dual science filters (``SCIFILT1``,
        ``SCIFILT2``) and slit name (``SLITNAME``), and sets instrument-
        specific reduced columns.

        Parameters
        ----------
        logger : logging.Logger
            Ginga application logger.
        """
        super().__init__(logger)
        self.pypeit_name = "keck_nirspec_high"
        raw_cols = [
            ("Type", "icon"),
            ("Frame No", "FRAMENO"),
            ("Name", "name"),
            ("Object", "OBJECT"),
            ("Img Type", "IMTYPE"),
            ("Filter 1", "FILTER1"),
            ("Filter 2", "FILTER2"),
            ("Slit", "MASKNAME"),
            ("Exp Time", "EXPTIME"),
            ("Last Changed", "st_mtime_str"),
        ]
        self.columns["raw"] = raw_cols
        self.columns["reduced"] = [
            ("Type", "icon"),
            ("Name", "name"),
            ("Slit", "SLITNAME"),
            ("Filter 1", "FILTER1"),
            ("Filter 2", "FILTER2"),
            ("Last Changed", "st_mtime_str"),
        ]

    def get_raw_info(self, path: str) -> Dict[str, object]:
        """Read display metadata from a NIRSPEC raw FITS file.

        Overrides :meth:`Instrument.get_raw_info` to handle NIRSPEC-specific
        header keywords: ``FRAMENUM`` (not ``FRAMENO``), ``TRUITIME`` (ramp
        integration time, not ``ELAPTIME``), ``IMTYPE`` (not ``KOAIMTYP``),
        and dual science filters ``SCIFILT1``/``SCIFILT2``.

        Parameters
        ----------
        path : str
            Absolute path to a NIRSPEC raw FITS file.

        Returns
        -------
        dict
            Keys: ``FRAMENO`` (``FRAMENUM``), ``OBJECT`` (``TARGNAME``),
            ``IMTYPE``, ``EXPTIME`` (``TRUITIME``), ``MASKNAME``
            (``SLITNAME``), ``FILTER1`` (``SCIFILT1``), ``FILTER2``
            (``SCIFILT2``).
        """
        with fits.open(path) as hdul:
            hdr = hdul[0].header
            info = self._read_header_fields(hdr)
            # NIRSPEC uses FRAMENUM (not FRAMENO), TRUITIME (not ELAPTIME),
            # IMTYPE (not KOAIMTYP), and TARGNAME (not OBJECT)
            info["FRAMENO"] = hdr.get("FRAMENUM", "N/A")
            info["OBJECT"] = hdr.get("TARGNAME", "N/A")
            info["IMTYPE"] = hdr.get("IMTYPE", "N/A")
            info["EXPTIME"] = hdr.get("TRUITIME", "N/A")
            info["MASKNAME"] = hdr.get("SLITNAME", "N/A")
            info["FILTER1"] = hdr.get("SCIFILT1", "N/A")
            info["FILTER2"] = hdr.get("SCIFILT2", "N/A")
            return info

    def get_display_image(self, raw_path: str) -> np.ndarray:
        """Build an overscan-subtracted, oriented display image for NIRSPEC.

        Delegates to PypeIt's ``buildimage_fromlist`` with ``biasframe``
        processing parameters for the single-detector configuration.

        Parameters
        ----------
        raw_path : str
            Absolute path to a NIRSPEC raw FITS file.

        Returns
        -------
        numpy.ndarray
            Processed 2-D image array.
        """
        spec = load_spectrograph("keck_nirspec_high")
        par = spec.default_pypeit_par()["calibrations"]["biasframe"]
        img = buildimage.buildimage_fromlist(spec, 1, par, [raw_path], mosaic=False)
        return img.image

    def get_reduced_info(self, path: str) -> Dict[str, object]:
        """Read display metadata from a NIRSPEC calibration directory or reduced file.

        Overrides :meth:`Instrument.get_reduced_info`.  For calibration
        directories the ``.pypeit`` setup file is parsed.

        Parameters
        ----------
        path : str
            Absolute path to a calibration directory (e.g.
            ``keck_nirspec_high_A/``) or a reduced FITS file.

        Returns
        -------
        dict
            Keys: ``SLITNAME`` (PypeIt ``decker``), ``FILTER1`` (PypeIt
            ``filter1``), ``FILTER2`` (PypeIt ``filter2``).
        """
        if os.path.isdir(path):
            cfg = self._read_pypeit_setup_config(path)
            return {
                "SLITNAME": cfg.get("decker", "N/A"),
                "FILTER1": cfg.get("filter1", "N/A"),
                "FILTER2": cfg.get("filter2", "N/A"),
            }
        try:
            with fits.open(path) as hdul:
                hdr = hdul[0].header
                return {
                    "SLITNAME": hdr.get("SLITNAME", "N/A"),
                    "FILTER1": hdr.get("SCIFILT1", "N/A"),
                    "FILTER2": hdr.get("SCIFILT2", "N/A"),
                }
        except Exception:
            return {}
