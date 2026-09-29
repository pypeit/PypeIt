"""
Calibration-directory utilities for the quicklook viewer.

These functions are generic (not instrument-specific).  They live in the
viewer instead of :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
because :func:`recommend_calibrations` depends on :mod:`pypeit.scripts.ql`,
which itself imports the spectrograph classes.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List

from pypeit.inputfiles import PypeItFile
from pypeit.pypeitsetup import PypeItSetup
from pypeit.scripts.ql import match_to_calibs


def read_pypeit_setup_config(dirpath: str, spec_name: str, logger) -> Dict[str, str]:
    """Return the instrument configuration from the setup block of the
    ``.pypeit`` file in a calibration directory.

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
        string.  Empty when no ``.pypeit`` file is found or parsing fails.
    """
    # Skipping unrelated directories (.., Calibrations/, Science/, ...) avoids
    # unnecessary I/O while browsing.
    _dirpath = Path(dirpath)
    if not re.match(rf'^{re.escape(spec_name)}_[A-Za-z]$', _dirpath.resolve().name):
        return {}

    pypeit_files = sorted(_dirpath.glob('*.pypeit')) or sorted(_dirpath.glob('*/*.pypeit'))
    if not pypeit_files:
        logger.info(f"No .pypeit files found in or under: {dirpath}")
        return {}

    try:
        setup = PypeItFile.from_file(str(pypeit_files[0]), vet=False).setup
    except Exception as exc:
        logger.warning(f"Could not parse pypeit file in {dirpath}: {exc}")
        return {}
    if not setup:
        logger.debug(f"No setup block found in {pypeit_files[0]}")
        return {}

    # PypeItFile flattens the setup block; the setup name (``Setup A``) and
    # any sub-block headings (e.g. ``--``) are left as keys with None values.
    # Nested (e.g. per-detector) blocks are skipped.
    config = {k: str(v) for k, v in setup.items() if v is not None and not isinstance(v, dict)}
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
