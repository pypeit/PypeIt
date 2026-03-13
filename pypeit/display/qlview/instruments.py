from __future__ import annotations

import os
from typing import Dict, List

import numpy as np
from astropy.io import fits


_BASE_RAW_COLUMNS = [
    ("Type", "icon"),
    ("Frame No", "FRAMENO"),
    ("Name", "name"),
    ("Object", "OBJECT"),
    ("Img Type", "IMTYPE"),
    ("Mask Name", "MASKNAME"),
    ("Obs Mode", "OBSMODE"),
    ("Exp Time", "EXPTIME"),
    ("Last Changed", "st_mtime_str"),
]

_BASE_REDUCED_COLUMNS = [
    ("Type", "icon"),
    ("Name", "name"),
    ("Mask Name", "MASKNAME"),
    ("Filter", "FILTER"),
    ("Slit Width", "SLITWIDTH"),
    ("Last Changed", "st_mtime_str"),
]

_MOSFIRE_REDUCED_COLUMNS = [
    ("Type", "icon"),
    ("Name", "name"),
    ("CSU Mask", "MASKNAME"),
    ("Filter", "FILTER1"),
    ("Dispname", "FILTER2"),
    ("Slit Width", "SLITWIDTH"),
    ("Last Changed", "st_mtime_str"),
]


class Instrument:
    """Base class for instrument-specific behavior."""

    pypeit_name: str = ""
    instrume_value: str = ""  # Expected value of the INSTRUME FITS keyword
    detector_prefix: str = "MSC"  # Prefix for --slitspatnum (MSC for mosaics, DET for single detectors)

    def __init__(self, logger) -> None:
        self.logger = logger
        # Per-view column definitions: keys are "raw" and "reduced".
        # Each value is a list of (display_name, attr_name) tuples matching
        # the format expected by Ginga's TreeView.setup_table().
        self.columns: Dict[str, List] = {
            "raw": list(_BASE_RAW_COLUMNS),
            "reduced": list(_BASE_REDUCED_COLUMNS),
        }

    def get_display_image(self, raw_path: str) -> np.ndarray:
        raise NotImplementedError

    def get_raw_info(self, path: str) -> Dict[str, object]:
        raise NotImplementedError

    def get_reduced_info(self, path: str) -> Dict[str, object]:
        """Read metadata from a reduced FITS file or calibration directory.

        Subclasses should override this to map instrument-specific keys.
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

    @staticmethod
    def _read_header_fields(header) -> Dict[str, object]:
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
        super().__init__(logger)
        self.pypeit_name = "keck_deimos"
        # Raw columns same as base; override reduced with DEIMOS-specific labels.
        # MASKNAME → decker (slit/mask), FILTER → dispname (grating),
        # SLITWIDTH → filter1 (blocking filter — most useful distinguishing column).
        self.columns["reduced"] = [
            ("Type", "icon"),
            ("Name", "name"),
            ("Mask/Slit", "MASKNAME"),
            ("Grating", "FILTER"),
            ("Blocking Filter", "SLITWIDTH"),
            ("Last Changed", "st_mtime_str"),
        ]

    def get_display_image(self, raw_path: str) -> np.ndarray:
        """Build an overscan-subtracted display image of all 8 DEIMOS detectors.

        Opens the file once via ``pypeit.io.fits_open``, reads all 8 chip
        extensions, subtracts the per-row median overscan bias, trims the
        pre- and post-scan columns, then assembles the chips into a 2×4 grid
        (detectors 1–4 on the bottom row, 5–8 on the top row) matching the
        DEIMOS focal-plane layout.
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

    def get_raw_info(self, path: str) -> Dict[str, object]:
        with fits.open(path) as hdul:
            return self._read_header_fields(hdul[0].header)

    def get_reduced_info(self, path: str) -> Dict[str, object]:
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
        super().__init__(logger)
        self.pypeit_name = "keck_mosfire"
        raw_cols = list(_BASE_RAW_COLUMNS)
        # Insert "Dither Pos" after "Name" (index 2) and before "Object"
        raw_cols.insert(3, ("Dither Pos", "DITHER_POS"))
        self.columns["raw"] = raw_cols
        self.columns["reduced"] = list(_MOSFIRE_REDUCED_COLUMNS)

    def get_raw_info(self, path: str) -> Dict[str, object]:
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


class InstrumentRegistry:
    def __init__(self, logger) -> None:
        self.logger = logger
        self._registry = {
            "DEIMOS": DEIMOS,
            "MOSFIRE": MOSFIRE,
        }

    def create(self, name: str) -> Instrument:
        cls = self._registry.get(name)
        if cls is None:
            self.logger.error(f"Instrument not recognized: {name}")
            cls = DEIMOS
        return cls(self.logger)

    def names(self) -> List[str]:
        return list(self._registry.keys())
