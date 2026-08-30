"""
Classes used by the quicklook viewer to configure the frontend.

These classes serve two purposes:

1. If a class is in this file, the quicklook viewer will display it as an opion
in the instrument dropdown.

2. They provide configuration information needed to display instrument info,
for example column headers, rendering the raw image, etc.

Instrument classes must be added to the InstrumentRegistry class to be seen!
"""

from __future__ import annotations

import os
from typing import Dict, List

import numpy as np
from astropy.io import fits



class Instrument:
    """Base class for instrument-specific behavior."""

    pypeit_name: str = ""
    # TODO: This is keck specific behavior, should change this soon
    instrume_value: str = ""  # Expected value of the INSTRUME FITS keyword
    # TODO: should this default to DET?
    detector_prefix: str = "MSC"  # Prefix for --slitspatnum (MSC for mosaics, DET for single detectors)

    def __init__(self, logger) -> None:
        """Initialise the instrument with a logger and default column definitions.

        Parameters
        ----------
        logger : logging.Logger
            Ginga application logger used for debug/warning messages throughout
            the class hierarchy.

        Notes
        -----
        ``self.columns`` is a dict with keys ``"raw"`` and ``"reduced"``.
        Each value is a list of ``(display_name, attr_name)`` tuples that match
        the format expected by ``Ginga.gw.Widgets.TreeView.setup_table()``.
        Subclasses must populate these lists in their own ``__init__`` to suit
        the instrument's FITS header vocabulary.
        """
        self.logger = logger
        # Per-view column definitions: keys are "raw" and "reduced".
        # Each value is a list of (display_name, attr_name) tuples matching
        # the format expected by Ginga's TreeView.setup_table().
        # Concrete subclasses are responsible for defining both lists.
        self.columns: Dict[str, List] = {
            "raw": [],
            "reduced": [],
        }

    def get_display_image(self, raw_path: str) -> np.ndarray:
        """Return a 2-D float array suitable for display in the Ginga viewer.

        Generally, this should just call buildimage.buildimage_fromlist, but 
        custom implementations are allowed. Note that the coordinates of this
        output image must align with those used by PypeIt (i.e. within a spec2D
        object) in order for each rendered element to stay aligned.

        Echelle spectrographs currently must rotate their images such that the
        spectra axis runs vertically (which is typically 90 deg from how they
        are normally viewed) in order to keep the code general, this is a place
        for future improvement.

        Parameters
        ----------
        raw_path : str
            Absolute path to a raw FITS file for this instrument.

        Returns
        -------
        numpy.ndarray
            2-D array with shape ``(nrows, ncols)`` in display orientation
            (spatial axis along columns, spectral axis along rows). (Is this
            true? Echelle spectrographs are sideways, must we enforce this?)

        Raises
        ------
        NotImplementedError
            Subclasses must override this method.
        """
        raise NotImplementedError

    def get_raw_info(self, path: str) -> Dict[str, object]:
        """Read per-file display metadata from a raw FITS file.

        The returned dict is merged into the ``Bunch`` that populates a row in
        the raw-data tree view.  Keys must match the ``attr_name`` entries in
        ``self.columns["raw"]``.

        Parameters
        ----------
        path : str
            Absolute path to a raw FITS file.

        Returns
        -------
        dict
            Mapping of column attribute name → display value (typically a
            string or number).

        Raises
        ------
        NotImplementedError
            Subclasses must override this method.
        """
        raise NotImplementedError

    def get_reduced_info(self, path: str) -> Dict[str, object]:
        """Read metadata from a reduced FITS file or calibration directory.

        Called by :class:`~.backends.LocalFileBrowserBackend` when populating
        the reduced-data tree view.  For calibration *directories* the
        preferred source is the ``.pypeit`` setup file (via
        :meth:`_read_pypeit_setup_config`); for individual FITS files the
        primary FITS header is used.

        Parameters
        ----------
        path : str
            Absolute path to a reduced FITS file **or** a calibration
            directory (e.g. ``keck_mosfire_A/``).

        Returns
        -------
        dict
            Mapping of column attribute name → display value.  Keys must
            match the ``attr_name`` entries in ``self.columns["reduced"]``.
            Returns an empty dict when no metadata can be extracted.
        """
        return {}

    def _read_pypeit_setup_config(self, dirpath: str) -> Dict[str, object]:
        """Parse the first ``.pypeit`` file found in *dirpath* and return
        the flat instrument-configuration dict from its setup block.

        The setup block YAML looks like::

            Setup A:
              --:
                dispname: 600ZD
                decker: 0.75arcsec
                ...
              '01': {binning: 1,1, ...}

        The ``'--'`` entry holds the per-instrument configuration keys.
        Returns an empty dict when no pypeit file is found or parsing fails.

        This is largely taken from core PypeIt functionality, but I couldn't
        find exactly what I needed as something reasonably callable (could have
        missed it!) so there's some duplicated code here. This is potentially
        a source of errors if the pypeit file format is changed, since this
        won't pick up any parser updates.
        """
        import glob as _glob
        import re
        import yaml

        # Only attempt to parse directories whose name matches the PypeIt
        # calibration-set naming convention: {spectrograph}_{letter}
        # (e.g. keck_mosfire_A, keck_deimos_B).  This skips .., Calibrations/,
        # Science/, and any other unrelated directories without I/O overhead.
        dir_name = os.path.basename(os.path.normpath(dirpath))
        if not re.match(rf'^{re.escape(self.pypeit_name)}_[A-Za-z]$', dir_name):
            return {}

        # Search the directory directly, then fall back one level deeper so the
        # method works whether dirpath is a setup dir (keck_mosfire_A/) or its
        # parent (the reductions root).
        pypeit_files = sorted(_glob.glob(os.path.join(dirpath, "*.pypeit")))
        if not pypeit_files:
            pypeit_files = sorted(_glob.glob(os.path.join(dirpath, "*", "*.pypeit")))
        if not pypeit_files:
            self.logger.info(f"No .pypeit files found in or under: {dirpath}")
            return {}
        try:
            with open(pypeit_files[0]) as fh:
                content = fh.read()
            # Extract the text between "setup read" and "setup end"
            match = re.search(r'setup read\n(.*?)setup end', content,
                              re.DOTALL | re.IGNORECASE)
            if not match:
                self.logger.debug("No setup block found")
                return {}
            setup_text = match.group(1)
            parsed = yaml.safe_load(setup_text)
            if not isinstance(parsed, dict):
                self.logger.debug("Could not parse setup block")
                return {}
            # parsed: {'Setup A': {'--': {config}, '01': {...}, ...}}
            # OR for some instruments: {'Setup A': {config_key: val, ...}}
            first_setup = next(iter(parsed.values()))
            if not isinstance(first_setup, dict):
                self.logger.debug(f"Unexpected setup type {type(first_setup)} in {dirpath}")
                return {}

            self.logger.debug(f"setup keys in {dirpath}: {list(first_setup.keys())}")

            # Prefer the '--' sub-block (non-detector-specific config).  When that
            # is absent or empty, fall back to the top-level dict itself — some
            # instruments write config keys directly under the Setup name.
            inner = first_setup.get("--") or first_setup.get("-")
            if not inner or not isinstance(inner, dict):
                # Filter out detector-index keys ('01', '02', …) which are ints or
                # short digit strings; keep the instrument-config key/value pairs.
                inner = {
                    k: v for k, v in first_setup.items()
                    if not (isinstance(k, int) or (isinstance(k, str) and k.isdigit()))
                }

            result = {k: "N/A" if v is None else str(v) for k, v in inner.items()
                      if not isinstance(v, dict)}  # skip nested sub-blocks
            self.logger.debug(f"pypeit config for {dirpath}: {result}")
            return result
        except Exception as exc:
            self.logger.warning(f"Could not parse pypeit file in {dirpath}: {exc}")
            return {}

    def recommend_calibrations(self, raw_path: str, cal_root: str) -> List[str]:
        """Return candidate calibration directories ranked by compatibility.

        Uses PypeIt's own configuration-matching logic (``PypeItSetup`` +
        ``match_to_calibs``) so that instrument-specific ``configuration_keys``
        and tolerances are respected automatically.

        This is used to attempt to select calibrations for the user.

        Parameters
        ----------
        raw_path : str
            Path to the selected raw FITS file.
        cal_root : str
            Root directory to search for calibration directories.

        Returns
        -------
        list of str
            Calibration directory paths, best match first.  Returns an empty
            list when no match is found or when an error occurs.
        """
        if not self.pypeit_name:
            return []

        try:
            from pypeit.pypeitsetup import PypeItSetup
            from pypeit.scripts.ql import match_to_calibs
        except Exception as e:
            self.logger.warning(f"Could not import PypeIt calibration API: {e}")
            return []

        try:
            ps = PypeItSetup.from_rawfiles([raw_path], self.pypeit_name)
            ps.run(setup_only=True)
        except Exception as e:
            self.logger.warning(f"PypeItSetup failed for {raw_path}: {e}")
            return []

        try:
            matched = match_to_calibs(ps, cal_root)
        except Exception as e:
            self.logger.warning(f"match_to_calibs failed: {e}")
            return []

        results = []
        for setup_match in matched.values():
            if setup_match is None:
                continue
            calib_dir = setup_match.get("calib_dir")
            if calib_dir is not None:
                results.append(str(calib_dir))
        return results

    # TODO: should this be implemented at all, or just send it to each instrument
    # class?
    @staticmethod
    def _read_header_fields(header) -> Dict[str, object]:
        """Extract the common set of display fields from a raw FITS primary header.

        This is a convenience helper called by subclass ``get_raw_info``
        implementations.  It populates the keys shared by all instruments;
        subclasses should override individual entries afterward to handle
        instrument-specific keyword names.

        This super implementation is likely overkill and should just be
        implemented from scratch in each instrument class.

        Parameters
        ----------
        header : astropy.io.fits.Header
            Primary HDU header of a raw FITS file.

        Returns
        -------
        dict
            Keys: ``OBJECT``, ``FRAMENO``, ``IMTYPE``, ``MASKNAME``,
            ``OBSMODE``, ``EXPTIME``.  ``EXPTIME`` falls back through
            ``TTIME``, ``ITIME``, ``ETIME``, and ``ELAPTIME`` in that order.

        Examples
        --------
        Typical usage inside a subclass ``get_raw_info``::

            info = self._read_header_fields(hdr)
            info["MASKNAME"] = hdr.get("SLMSKNAM", "N/A")   # DEIMOS override
            return info
        """
        header_dict = {
            "OBJECT": header.get("OBJECT", "N/A"),
            "FRAMENO": header.get("FRAMENO", "N/A"),
            "IMTYPE": header.get("KOAIMTYP", "N/A"),
            "MASKNAME": header.get("MASKNAME", "N/A"),
            "OBSMODE": header.get("OBSMODE", "N/A"),
            "EXPTIME": header.get("EXPTIME", None),
        }
        if header_dict["EXPTIME"] is None:
            header_dict["EXPTIME"] = header.get("TTIME", None)
        if header_dict["EXPTIME"] is None:
            header_dict["EXPTIME"] = header.get("ITIME", None)
        if header_dict["EXPTIME"] is None:
            header_dict["EXPTIME"] = header.get("ETIME", None)
        if header_dict["EXPTIME"] is None:
            header_dict["EXPTIME"] = header.get("ELAPTIME", "N/A")
        return header_dict


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
        from pypeit.io import fits_open
        from pypeit.spectrographs.keck_deimos import (
            KeckDEIMOSSpectrograph,
            deimos_read_1chip,
        )
        from pypeit.core.mosaic import build_image_mosaic

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
        from pypeit.io import fits_open

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
        from pypeit.spectrographs.util import load_spectrograph
        from pypeit.images import buildimage

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


class NIRES(Instrument):
    """Keck NIRES — near-IR, fixed-format echelle, single detector.
    
    Untested!"""

    instrume_value = "NIRES"
    detector_prefix = "DET"

    def __init__(self, logger) -> None:
        """Initialise the NIRES instrument with Keck-NIRES–specific column definitions.

        Defines a custom raw column list that omits mask/obsmode columns
        (NIRES is fixed-format) and adds a ``DITHER_POS`` column decoded from
        the ``DPATNAME``/``DPATIPOS`` header pair.

        
        Parameters
        ----------
        logger : logging.Logger
            Ginga application logger.
        """
        super().__init__(logger)
        self.pypeit_name = "keck_nires"
        self.columns["raw"] = [
            ("Type", "icon"),
            ("Frame No", "FRAMENO"),
            ("Name", "name"),
            ("Dither Pos", "DITHER_POS"),
            ("Object", "OBJECT"),
            ("Obs Type", "IMTYPE"),
            ("Exp Time", "EXPTIME"),
            ("Last Changed", "st_mtime_str"),
        ]
        # Fixed-format echelle — no grating/filter variation to show
        self.columns["reduced"] = [
            ("Type", "icon"),
            ("Name", "name"),
            ("Mask Name", "MASKNAME"),
            ("Filter", "FILTER"),
            ("Slit Width", "SLITWIDTH"),
            ("Last Changed", "st_mtime_str"),
        ]

    def get_raw_info(self, path: str) -> Dict[str, object]:
        """Read display metadata from a NIRES raw FITS file.

        Overrides :meth:`Instrument.get_raw_info` to handle NIRES-specific
        header keywords: ``FRAMENUM`` (not ``FRAMENO``), ``ITIME`` (not
        ``ELAPTIME``), ``OBSTYPE`` (not ``KOAIMTYP``), and the
        ``DPATNAME``/``DPATIPOS`` dither-position pair.

        Parameters
        ----------
        path : str
            Absolute path to a NIRES raw FITS file.

        Returns
        -------
        dict
            Keys: ``FRAMENO`` (``FRAMENUM``), ``OBJECT`` (``TARGNAME``),
            ``IMTYPE`` (``OBSTYPE``), ``EXPTIME`` (``ITIME``),
            ``DITHER_POS`` (character decoded from ``DPATNAME[DPATIPOS-1]``).
        """
        with fits.open(path) as hdul:
            hdr = hdul[0].header
            info = self._read_header_fields(hdr)
            # NIRES uses FRAMENUM (not FRAMENO), ITIME (not ELAPTIME),
            # and OBSTYPE (not KOAIMTYP)
            info["FRAMENO"] = hdr.get("FRAMENUM", "N/A")
            info["IMTYPE"] = hdr.get("OBSTYPE", "N/A")
            info["EXPTIME"] = hdr.get("ITIME", "N/A")
            info["OBJECT"] = hdr.get("TARGNAME", "N/A")
            # Dither position: decode DPATIPOS (1-based index) via DPATNAME
            dpat = hdr.get("DPATNAME", "")
            dpos = hdr.get("DPATIPOS")
            if dpos is not None and dpat and 1 <= dpos <= len(dpat):
                info["DITHER_POS"] = dpat[dpos - 1]
            else:
                info["DITHER_POS"] = "N/A"
            return info

    def get_display_image(self, raw_path: str) -> np.ndarray:
        """Build an overscan-subtracted, oriented display image for NIRES.

        Delegates to PypeIt's ``buildimage_fromlist`` with ``biasframe``
        processing parameters (overscan subtraction, trimming, orientation).

        Parameters
        ----------
        raw_path : str
            Absolute path to a NIRES raw FITS file.

        Returns
        -------
        numpy.ndarray
            Processed 2-D image array.
        """
        from pypeit.spectrographs.util import load_spectrograph
        from pypeit.images import buildimage

        spec = load_spectrograph("keck_nires")
        par = spec.default_pypeit_par()["calibrations"]["biasframe"]
        img = buildimage.buildimage_fromlist(spec, 1, par, [raw_path], mosaic=False)
        return img.image

    def get_reduced_info(self, path: str) -> Dict[str, object]:
        """Read display metadata from a NIRES calibration directory or reduced file.

        Overrides :meth:`Instrument.get_reduced_info`.  NIRES is a
        fixed-format echelle so the reduced columns only carry slit/decker
        information parsed from the ``.pypeit`` file.

        Parameters
        ----------
        path : str
            Absolute path to a calibration directory (e.g. ``keck_nires_A/``)
            or a reduced FITS file.

        Returns
        -------
        dict
            Keys: ``MASKNAME`` (decker/slit), ``FILTER`` (dispname),
            ``SLITWIDTH`` (slit width).
        """
        if os.path.isdir(path):
            cfg = self._read_pypeit_setup_config(path)
            return {
                "MASKNAME": cfg.get("decker", "N/A"),
                "FILTER": cfg.get("dispname", "N/A"),
                "SLITWIDTH": cfg.get("slitwid", "N/A"),
            }
        try:
            with fits.open(path) as hdul:
                hdr = hdul[0].header
                return {
                    "MASKNAME": hdr.get("SLITNAME", "N/A"),
                    "FILTER": hdr.get("INSTR", "N/A"),
                    "SLITWIDTH": "N/A",
                }
        except Exception:
            return {}


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
        from pypeit.spectrographs.util import load_spectrograph
        from pypeit.images import buildimage

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
        from pypeit.spectrographs.util import load_spectrograph
        from pypeit.images import buildimage

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
        from pypeit.spectrographs.util import load_spectrograph
        from pypeit.images import buildimage

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


class InstrumentRegistry:
    def __init__(self, logger) -> None:
        """Initialise the registry and register all built-in instrument classes.

        Parameters
        ----------
        logger : logging.Logger
            Ginga application logger, forwarded to each
            :class:`Instrument` instance created via :meth:`create`.
        """
        self.logger = logger
        self._registry = {
            "DEIMOS": DEIMOS,
            # "HIRES": HIRES,
            # "LRIS Blue": LRISBlue,
            # "LRIS Red": LRISRed,
            "MOSFIRE": MOSFIRE,
            # "NIRES": NIRES,
            # "NIRSPEC": NIRSPEC,
        }

    def create(self, name: str) -> Instrument:
        """Instantiate and return the :class:`Instrument` for *name*.

        Parameters
        ----------
        name : str
            Display name as it appears in the instrument combo box
            (e.g. ``"DEIMOS"``, ``"LRIS Blue"``).

        Returns
        -------
        Instrument
            A freshly constructed instrument instance.

        Notes
        -----
        Falls back to :class:`DEIMOS` and logs an error when *name* is not
        found in the registry, so callers always receive a usable object.
        """
        cls = self._registry.get(name)
        if cls is None:
            self.logger.error(f"Instrument not recognized: {name}")
            cls = DEIMOS
        return cls(self.logger)

    def names(self) -> List[str]:
        """Return the list of registered instrument display names.

        Returns
        -------
        list of str
            Names in insertion order, matching the order shown in the
            instrument combo box.
        """
        return list(self._registry.keys())

    def instrume_values(self) -> List[tuple]:
        """Return ``(display_name, instrume_value)`` pairs for all registered instruments.

        Reads ``instrume_value`` directly from each class object rather than
        constructing an instance, making this safe to call when you only need
        the FITS keyword value for matching purposes.

        Returns
        -------
        list of (str, str)
            ``(display_name, instrume_value)`` tuples in registration order.
        """
        return [(name, cls.instrume_value) for name, cls in self._registry.items()]
