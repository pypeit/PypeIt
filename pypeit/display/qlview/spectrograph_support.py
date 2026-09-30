"""
Interface between the quicklook viewer and the PypeIt spectrograph classes.

Instrument-specific behavior is provided by the ``qlview_*`` hooks of
:class:`~pypeit.spectrographs.spectrograph.Spectrograph`.  Spectrographs with
``qlview_supported = True`` are offered in the viewer's instrument selector.
Everywhere in the viewer (including the remote HTTP API), instruments are
identified by their PypeIt spectrograph name (e.g. ``keck_deimos``); the
``qlview_label`` is only used for display.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from astropy.io import fits

from pypeit.images.detector_container import DetectorContainer
from pypeit.images.mosaic import Mosaic
from pypeit.spectrographs.spectrograph import Spectrograph
from pypeit.spectrographs.util import load_spectrograph, spectrograph_classes

from .calib_utils import read_pypeit_setup_config

DEFAULT_SPECTROGRAPH = 'keck_deimos'
"""Spectrograph selected when the viewer starts."""


@lru_cache(maxsize=None)
def supported_spectrographs() -> Tuple[str, ...]:
    """Return the names of the spectrographs offered in the viewer, sorted by label.

    The result is cached; the set of spectrograph classes is fixed at import.

    Returns
    -------
    tuple of str
        PypeIt spectrograph names with ``qlview_supported = True``.
    """
    classes = [c for c in spectrograph_classes().values() if c.qlview_supported]
    return tuple(c.name for c in sorted(classes, key=lambda c: c.qlview_label or c.name))


def load_qlview_spectrograph(name: str) -> Spectrograph:
    """Instantiate a spectrograph supported by the viewer.

    Parameters
    ----------
    name : str
        PypeIt spectrograph name (e.g. ``keck_deimos``).

    Returns
    -------
    Spectrograph
        The spectrograph instance.

    Raises
    ------
    ValueError
        If *name* is not a spectrograph supported by the viewer.
    """
    if name not in supported_spectrographs():
        raise ValueError(f"'{name}' is not a spectrograph supported by the quicklook viewer; "
                         f"options are: {', '.join(supported_spectrographs())}")
    return load_spectrograph(name)


def label(spec: Spectrograph | str) -> str:
    """Return the name shown in the instrument selector.

    Parameters
    ----------
    spec : Spectrograph or str
        Spectrograph instance or PypeIt spectrograph name.

    Returns
    -------
    str
        The spectrograph's ``qlview_label``, or its name if it has none.
    """
    if isinstance(spec, str):
        spec = spectrograph_classes()[spec]
    return spec.qlview_label or spec.name


def match_header_name(instrume: str) -> Optional[str]:
    """Find the supported spectrograph whose ``header_name`` matches ``INSTRUME``.

    Parameters
    ----------
    instrume : str
        Value of the ``INSTRUME`` header keyword.

    Returns
    -------
    str or None
        Name of the first matching supported spectrograph (in the order of
        :func:`supported_spectrographs`), or None if there is no match.
    """
    instrume = instrume.strip().upper()
    classes = spectrograph_classes()
    for name in supported_spectrographs():
        header_name = classes[name].header_name
        if header_name is not None and header_name.upper() == instrume:
            return name
    return None


def build_columns(spec: Spectrograph, mode: str) -> List[Tuple[str, str]]:
    """Return the file-browser column definitions for *spec*.

    The spectrograph's instrument columns are combined with the viewer's
    file-system columns: the type icon first, then the file name (after the
    ``FRAMENO`` column, if the spectrograph has one), and the modification
    time last.

    Parameters
    ----------
    spec : Spectrograph
        Active spectrograph.
    mode : str
        ``"raw"`` or ``"reduced"``.

    Returns
    -------
    list of (str, str)
        ``(display_name, attr_name)`` tuples, as expected by Ginga's
        ``TreeView.setup_table()``.
    """
    cols = spec.qlview_raw_columns() if mode == 'raw' else spec.qlview_reduced_columns()
    cols = list(cols)
    name_at = 1 if cols and cols[0][1] == 'FRAMENO' else 0
    return [('Type', 'icon')] + cols[:name_at] + [('Name', 'name')] + cols[name_at:] \
            + [('Last Changed', 'st_mtime_str')]


def get_header_info(spec: Spectrograph, path: str, mode: str, logger) -> Dict[str, object]:
    """Read the file-browser column values for a raw file, reduced file, or
    calibration directory.

    Parameters
    ----------
    spec : Spectrograph
        Active spectrograph.
    path : str
        Path to a FITS file or, in ``"reduced"`` mode, a calibration
        directory (e.g. ``keck_deimos_A/``).
    mode : str
        ``"raw"`` or ``"reduced"``.
    logger : logging.Logger
        Viewer logger.

    Returns
    -------
    dict
        Mapping of column attribute name to display value.
    """
    if mode == 'raw':
        with fits.open(path) as hdul:
            return spec.qlview_raw_info(hdul[0].header)

    if Path(path).is_dir():
        # Calibration directories are described by the setup block of their
        # .pypeit file, whose keys are the reduced-column attribute names.
        cfg = read_pypeit_setup_config(path, spec.name, logger)
        return {key: cfg.get(key, 'N/A') for _, key in spec.qlview_reduced_columns()}
    try:
        with fits.open(path) as hdul:
            return spec.qlview_reduced_info(hdul[0].header)
    except Exception:
        return {}


def det_label(spec: Spectrograph, det_id: str) -> str:
    """Return the detector/mosaic name for a slit-trace index.

    Parameters
    ----------
    spec : Spectrograph
        Active spectrograph.
    det_id : str
        Zero-padded detector or mosaic index taken from a ``Slits`` file name
        (e.g. ``"01"``).

    Returns
    -------
    str
        E.g. ``"MSC01"`` for spectrographs with detector mosaics, otherwise
        ``"DET01"``.
    """
    prefix = Mosaic.name_prefix if spec.allowed_mosaics else DetectorContainer.name_prefix
    return f'{prefix}{det_id}'
