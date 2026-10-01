"""
Calibration-directory utilities for the quicklook viewer.

These functions are generic (not instrument-specific).  They live in the
viewer instead of :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
because :func:`recommend_calibrations` depends on :mod:`pypeit.scripts.ql`,
which itself imports the spectrograph classes.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Dict, List

import yaml

from pypeit.pypeitsetup import PypeItSetup
from pypeit.scripts.ql import match_to_calibs


def read_pypeit_setup_config(dirpath: str, spec_name: str, logger) -> Dict[str, str]:
    """Return the instrument configuration from the setup block of the
    ``.pypeit`` file in a calibration directory.

    Only the setup block is parsed, so this is fast and tolerant of problems
    elsewhere in the file.  The block looks like either::

        Setup A:
          --:
            dispname: 600ZD
            decker: 0.75arcsec
          '01': {binning: 1,1, ...}

    or, with the configuration keys directly under the setup name::

        Setup A:
          dispname: 600ZD
          decker: 0.75arcsec

    Parameters
    ----------
    dirpath : str
        Calibration directory.  Only directories named after the PypeIt
        calibration-set convention, ``{spec_name}_{letter}`` (e.g.
        ``keck_mosfire_A``), are parsed; anything else returns an empty dict.
        The ``.pypeit`` file is searched for in *dirpath* and then one level
        deeper.
    spec_name : str
        PypeIt spectrograph name (e.g. ``keck_deimos``).
    logger : logging.Logger
        Viewer logger.

    Returns
    -------
    dict
        Mapping of configuration key (e.g. ``decker``) to its value as a
        string (``"N/A"`` for empty values).  Empty when no ``.pypeit`` file
        is found or parsing fails.
    """
    # Skipping unrelated directories (.., Calibrations/, Science/, ...) avoids
    # unnecessary I/O while browsing.  Normalize without resolving symlinks,
    # so that a symlinked setup directory keeps its name.
    _dirpath = Path(os.path.normpath(dirpath))
    if not re.match(rf'^{re.escape(spec_name)}_[A-Za-z]$', _dirpath.name):
        return {}

    pypeit_files = sorted(_dirpath.glob('*.pypeit')) or sorted(_dirpath.glob('*/*.pypeit'))
    if not pypeit_files:
        logger.info(f"No .pypeit files found in or under: {dirpath}")
        return {}

    try:
        match = re.search(r'setup read\n(.*?)setup end', pypeit_files[0].read_text(),
                          re.DOTALL | re.IGNORECASE)
        if not match:
            logger.debug(f"No setup block found in {pypeit_files[0]}")
            return {}
        parsed = yaml.safe_load(match.group(1))
        if not isinstance(parsed, dict):
            logger.debug(f"Could not parse setup block in {pypeit_files[0]}")
            return {}
        # parsed: {'Setup A': {...}}
        first_setup = next(iter(parsed.values()))
        if not isinstance(first_setup, dict):
            logger.debug(f"Unexpected setup type {type(first_setup)} in {dirpath}")
            return {}

        # Prefer the '--' sub-block (non-detector-specific config); otherwise
        # use the top-level dict, without the detector-index keys ('01', ...).
        config = first_setup.get("--") or first_setup.get("-")
        if not isinstance(config, dict) or not config:
            config = {k: v for k, v in first_setup.items()
                      if not (isinstance(k, int) or (isinstance(k, str) and k.isdigit()))}

        # Nested (e.g. per-detector) blocks are skipped
        config = {k: "N/A" if v is None else str(v) for k, v in config.items()
                  if not isinstance(v, dict)}
    except Exception as exc:
        logger.warning(f"Could not parse pypeit file in {dirpath}: {exc}")
        return {}
    logger.debug(f"pypeit config for {dirpath}: {config}")
    return config


def recommend_calibrations(spec_name: str, raw_path: str, cal_root: str, logger) -> List[str]:
    """Return calibration directories compatible with a raw file.

    Uses PypeIt's own configuration matching (``PypeItSetup`` +
    ``match_to_calibs``), so the spectrograph's ``configuration_keys`` and
    tolerances are respected.

    Parameters
    ----------
    spec_name : str
        PypeIt spectrograph name (e.g. ``keck_deimos``).
    raw_path : str
        Path to the selected raw FITS file.
    cal_root : str
        Root directory to search for calibration directories.
    logger : logging.Logger
        Viewer logger.

    Returns
    -------
    list of str
        Calibration directory paths, best match first.  Empty when no match
        is found or when an error occurs.
    """
    # Run PypeIt's setup step on the single raw file, assigning it to a setup
    # (A, B, ...) exactly as run_pypeit would.
    try:
        ps = PypeItSetup.from_rawfiles([raw_path], spec_name)
        ps.run(setup_only=True)
    except Exception as e:
        logger.warning(f"PypeItSetup failed for {raw_path}: {e}", exc_info=True)
        return []

    # Compare that setup against every reduced calibration set under cal_root;
    # this is the same matching pypeit_ql uses to reuse existing calibrations.
    try:
        matched = match_to_calibs(ps, cal_root)
    except Exception as e:
        logger.warning(f"match_to_calibs failed for {cal_root}: {e}", exc_info=True)
        return []

    # ``matched`` maps each setup to either None (no compatible calibrations)
    # or a dict holding the matching ``calib_dir``.
    return [str(m['calib_dir']) for m in matched.values()
            if m is not None and m.get('calib_dir') is not None]
