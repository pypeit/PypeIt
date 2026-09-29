"""
Quicklook viewer configuration for Keck/DEIMOS.
"""

from __future__ import annotations

import os
from typing import Dict

import numpy as np
from astropy.io import fits

from pypeit.core.mosaic import build_image_mosaic
from pypeit.io import fits_open
from pypeit.spectrographs.keck_deimos import KeckDEIMOSSpectrograph, deimos_read_1chip

from .base import Instrument


class DEIMOS(Instrument):
    instrume_value = "DEIMOS"

    def __init__(self, logger) -> None:
        """Initialise the DEIMOS instrument with Keck-DEIMOS–specific column definitions.

        Overrides the base raw column list to use the DEIMOS FITS vocabulary
        (``SLMSKNAM``, ``GRATENAM``, ``DWFILNAM``, ``ELAPTIME``) and sets
        instrument-specific reduced columns.

        Parameters
        ----------
        logger : logging.Logger
            Ginga application logger.
        """
        super().__init__(logger)
        self.pypeit_name = "keck_deimos"
        # DEIMOS-specific raw columns: SLMSKNAM for mask, GRATENAM for grating,
        # DWFILNAM for blocking filter, TARGNAME for object, ELAPTIME for exp time.
        self.columns["raw"] = [
            ("Type", "icon"),
            ("Frame No", "FRAMENO"),
            ("Name", "name"),
            ("Object", "OBJECT"),
            ("Img Type", "IMTYPE"),
            ("Mask Name", "MASKNAME"),
            ("Grating", "GRATING"),
            ("Filter", "FILTER1"),
            ("Exp Time", "EXPTIME"),
            ("Last Changed", "st_mtime_str"),
        ]
        # Reduced columns: read from pypeit config keys (decker, dispname, filter1)
        self.columns["reduced"] = [
            ("Type", "icon"),
            ("Name", "name"),
            ("Mask/Slit", "MASKNAME"),
            ("Grating", "FILTER"),
            ("Blocking Filter", "SLITWIDTH"),
            ("Last Changed", "st_mtime_str"),
        ]

    def get_raw_info(self, path: str) -> Dict[str, object]:
        """Read display metadata from a DEIMOS raw FITS file.

        Overrides :meth:`Instrument.get_raw_info` to use the DEIMOS-specific
        FITS keyword names, which differ from the KOA defaults assumed by
        :meth:`Instrument._read_header_fields`.

        Parameters
        ----------
        path : str
            Absolute path to a DEIMOS raw FITS file.

        Returns
        -------
        dict
            Keys: ``OBJECT`` (``TARGNAME``), ``FRAMENO``, ``IMTYPE``
            (``KOAIMTYP``), ``MASKNAME`` (``SLMSKNAM``), ``GRATING``
            (``GRATENAM``), ``FILTER1`` (``DWFILNAM``), ``EXPTIME``
            (``ELAPTIME`` → ``EXPTIME``).
        """
        with fits.open(path) as hdul:
            hdr = hdul[0].header
            info = self._read_header_fields(hdr)
            # DEIMOS uses TARGNAME (not OBJECT), SLMSKNAM (not MASKNAME),
            # GRATENAM for grating, DWFILNAM for blocking filter,
            # ELAPTIME for exposure time, and KOAIMTYP for image type.
            info["OBJECT"] = hdr.get("TARGNAME", hdr.get("OBJECT", "N/A"))
            info["MASKNAME"] = hdr.get("SLMSKNAM", "N/A")
            info["GRATING"] = hdr.get("GRATENAM", "N/A")
            info["FILTER1"] = hdr.get("DWFILNAM", "N/A")
            info["IMTYPE"] = hdr.get("KOAIMTYP", "N/A")
            info["EXPTIME"] = hdr.get("ELAPTIME", hdr.get("EXPTIME", "N/A"))
            return info

    def get_display_image(self, raw_path: str) -> np.ndarray:
        """Build an overscan-subtracted display image using the PypeIt mosaic pipeline.

        Reads the 8 chips, subtracts the per-row median overscan bias,
        and assembles each of the four detector pairs (MSC01–MSC04) using
        :func:`pypeit.core.mosaic.build_image_mosaic` with the same affine
        transforms used during PypeIt reductions.  The four mosaics are then
        concatenated along the spatial axis to form the full display image.
        Padding is then added to each in order to make them the same visual size
        for concatenating and then rendering.

        The resulting image is in the same coordinate system as the
        :class:`~pypeit.slittrace.SlitTraceSet` objects produced by the
        reduction pipeline, so slit-trace overlays will be correctly registered.

        Parameters
        ----------
        raw_path : str
            Absolute path to the raw DEIMOS FITS file.

        Returns
        -------
        numpy.ndarray
            2-D float array with shape ``(nspec, 4*nspat_mosaic)`` where
            ``nspec`` and ``nspat_mosaic`` are determined by
            :func:`~pypeit.core.mosaic.prepare_mosaic`.
        """
        spectrograph = KeckDEIMOSSpectrograph()

        with fits_open(raw_path) as hdu:
            mosaic_images = []
            for mosaic_tuple in spectrograph.allowed_mosaics:
                det_blue, det_red = mosaic_tuple

                # Read trimmed, orientation-corrected data in (nspec, nspat) order.
                data_blue, oscan_blue = deimos_read_1chip(hdu, det_blue)
                data_red, oscan_red = deimos_read_1chip(hdu, det_red)

                # Per-spectral-row median overscan subtraction.
                data_blue = data_blue.astype(float)
                data_blue -= np.median(oscan_blue.astype(float), axis=1)[:, np.newaxis]
                data_red = data_red.astype(float)
                data_red -= np.median(oscan_red.astype(float), axis=1)[:, np.newaxis]

                # Build the mosaic using the same transforms as the reduction pipeline.
                msc = spectrograph.get_mosaic_par(mosaic_tuple, hdu=hdu)
                mosaic_img, _, _, _ = build_image_mosaic(
                    [data_blue, data_red], list(msc.tform)
                )
                mosaic_images.append(mosaic_img)

        # The four mosaics may differ slightly in nspec (axis 0) because each
        # MSC has a different rotation angle and prepare_mosaic computes the
        # bounding box independently.  Pad shorter mosaics with zeros so all
        # have the same nspec before concatenating along the spatial axis.
        max_nspec = max(img.shape[0] for img in mosaic_images)
        padded = []
        for img in mosaic_images:
            deficit = max_nspec - img.shape[0]
            if deficit > 0:
                img = np.pad(img, ((0, deficit), (0, 0)))
            padded.append(img)
        return np.concatenate(padded, axis=1)

    def get_display_image_simple(self, raw_path: str) -> np.ndarray:
        """Reference implementation: direct chip concatenation without mosaic transforms.

        This method preserves the original ``get_display_image`` implementation
        that was used before the mosaic-aware version.  It assembles all 8 chips
        into a 2×4 grid by simple concatenation (after overscan subtraction and
        trimming) without applying the per-detector affine transforms stored in
        :class:`~pypeit.spectrographs.keck_deimos.DEIMOSMosaicLookUp`.

        The result is *not* in the same coordinate frame as the PypeIt
        :class:`~pypeit.slittrace.SlitTraceSet`, so slit-trace overlays will be
        misregistered by up to ~30 pixels.  This method is kept as a reference
        to aid future development and debugging.

        Parameters
        ----------
        raw_path : str
            Absolute path to the raw DEIMOS FITS file.

        Returns
        -------
        numpy.ndarray
            2-D float array assembled as a simple 2×4 chip grid.
        """
        with fits_open(raw_path) as hdu:
            hdr0 = hdu[0].header
            binning = hdr0["BINNING"].split(",")
            precol = int(hdr0["PRECOL"]) // int(binning[0])
            postpix = int(hdr0["POSTPIX"]) // int(binning[0])

            chips = []
            for i in range(1, 9):
                data = hdu[i].data.astype(float)
                height, width = data.shape
                bias = np.median(data[:, width - postpix:], axis=1)
                data -= bias[:, np.newaxis]
                chips.append(data[:, precol: width - postpix])

        # Detectors 1–4: concatenate left-to-right, then flip the row upward
        r0 = np.flipud(np.concatenate(chips[:4], axis=1))
        # Detectors 5–8: flip each chip left-to-right, then concatenate
        r1 = np.concatenate([np.fliplr(c) for c in chips[4:]], axis=1)
        return np.concatenate((r1, r0), axis=0)

    def get_reduced_info(self, path: str) -> Dict[str, object]:
        """Read display metadata from a DEIMOS calibration directory or reduced file.

        Overrides :meth:`Instrument.get_reduced_info`.  For calibration
        *directories* the ``.pypeit`` setup file is the authoritative source
        and is read via :meth:`_read_pypeit_setup_config`.  For individual
        FITS files the primary header is used as a fallback.

        Parameters
        ----------
        path : str
            Absolute path to a calibration directory (e.g. ``keck_deimos_A/``)
            or a reduced FITS file.

        Returns
        -------
        dict
            Keys: ``MASKNAME`` (PypeIt ``decker`` / ``SLMSKNAM``),
            ``FILTER`` (grating name, PypeIt ``dispname`` / ``GRATENAM``),
            ``SLITWIDTH`` (blocking filter, PypeIt ``filter1``).
        """
        if os.path.isdir(path):
            # Prefer metadata from the pypeit file in this calibration directory.
            # DEIMOS configuration_keys: dispname (grating), decker (slit/mask),
            # binning, dispangle, amp, filter1.
            cfg = self._read_pypeit_setup_config(path)
            return {
                # decker is the slit-mask name for DEIMOS (e.g. "GS62" or "0.75arcsec")
                "MASKNAME": cfg.get("decker", "N/A"),
                # dispname is the grating name (e.g. "600ZD", "830G")
                "FILTER": cfg.get("dispname", "N/A"),
                # filter1 is the blocking filter (e.g. "GG455")
                "SLITWIDTH": cfg.get("filter1", "N/A"),
            }
        try:
            with fits.open(path) as hdul:
                hdr = hdul[0].header
                return {
                    "MASKNAME": hdr.get("MASKNAME", "N/A"),
                    "FILTER": hdr.get("GRATENAM", hdr.get("FILTER", "N/A")),
                    "SLITWIDTH": hdr.get("SLITNAME", hdr.get("SLITWIDTH", "N/A")),
                }
        except Exception:
            return {}
