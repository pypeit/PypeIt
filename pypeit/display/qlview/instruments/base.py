"""
Base class defining the interface between an instrument and the quicklook viewer.
"""

from __future__ import annotations

import glob
import os
import re
from typing import Dict, List

import numpy as np
import yaml

from pypeit.pypeitsetup import PypeItSetup
from pypeit.scripts.ql import match_to_calibs


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
        pypeit_files = sorted(glob.glob(os.path.join(dirpath, "*.pypeit")))
        if not pypeit_files:
            pypeit_files = sorted(glob.glob(os.path.join(dirpath, "*", "*.pypeit")))
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
        # Instruments without a PypeIt spectrograph name can't be matched
        if not self.pypeit_name:
            return []

        # Run PypeIt's setup step on the single raw file.  This reads its
        # headers and assigns it to a setup (A, B, ...) using the
        # spectrograph's configuration_keys, exactly as run_pypeit would.
        try:
            ps = PypeItSetup.from_rawfiles([raw_path], self.pypeit_name)
            ps.run(setup_only=True)
        except Exception as e:
            self.logger.warning(f"PypeItSetup failed for {raw_path}: {e}", exc_info=True)
            return []

        # Compare that setup against every reduced calibration set found under
        # cal_root.  This is the same matching pypeit_ql uses to reuse existing
        # calibrations.
        try:
            matched = match_to_calibs(ps, cal_root)
        except Exception as e:
            self.logger.warning(f"match_to_calibs failed for {cal_root}: {e}", exc_info=True)
            return []

        # ``matched`` maps each setup in ``ps`` to either None (no compatible
        # calibrations) or a dict holding the matching ``calib_dir``.
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
