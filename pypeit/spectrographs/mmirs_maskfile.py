"""
Parser for MMT/MMIRS xfitmask ``.msk`` mask-design files.

.. include:: ../include/links.rst
"""
import re
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.coordinates import Angle
from astropy import units

from pypeit import PypeItError


def _sexed_ra_hours(s):
    """
    Convert a sexagesimal RA string in *hours* to degrees.

    Parameters
    ----------
    s : :obj:`str`
        Sexagesimal right ascension in hours (e.g. ``'17:22:20.2490'``).

    Returns
    -------
    :obj:`float`
        Right ascension in degrees.
    """
    return Angle(s, unit=units.hourangle).to('deg').value


def _sexed_dec_deg(s):
    """
    Convert a sexagesimal Dec string in *degrees* to degrees.

    Parameters
    ----------
    s : :obj:`str`
        Sexagesimal declination in degrees (e.g. ``'65:56:13.040'``).

    Returns
    -------
    :obj:`float`
        Declination in degrees.
    """
    return Angle(s, unit=units.deg).value


def read_mmirs_maskfile(path):
    """
    Parse an MMIRS ``.msk`` mask-design file.

    The ``.msk`` file (xfitmask output) is plain text with mixed tab/space
    separation and three sections: a mask-level key/value header, a
    ``GuideStars`` block (ignored), and a slit table.  The slit-table columns
    are ``slit ra dec x y target object type height width offset theta bbox
    polygon`` (field indices 0..13), where ``bbox`` (field 12) and ``polygon``
    (field 13) are space-separated numbers within a single tab-delimited field.

    Parameters
    ----------
    path : :obj:`str` or `Path`_
        Path to the ``.msk`` file.

    Returns
    -------
    header : :obj:`dict`
        Mask-level metadata: ``label``, ``ra_deg``, ``dec_deg``, ``pa``,
        ``scale``, ``arc2mm``, ``corners`` (list of 4 floats, mm), and, when
        present in the file, ``grism`` and ``filter``.
    slits : `astropy.table.Table`_
        One row per slit with columns ``slit`` (int), ``ra_deg`` (float),
        ``dec_deg`` (float), ``x_mm`` (float), ``y_mm`` (float), ``target``
        (int), ``object`` (str), ``type`` (str), ``height_mm`` (float),
        ``width_mm`` (float), ``offset_mm`` (float, target displacement from
        the slit center along the slit long axis), ``theta_deg`` (float, slit
        tilt), ``bbox`` (float, shape ``(N, 4)``), ``polygon`` (float, shape
        ``(N, 8)``).

    Raises
    ------
    PypeItError
        If the file has no slit-table header row, or no slit rows can be
        parsed.
    """
    lines = Path(path).read_text().splitlines()

    header = {}
    slit_hdr_idx = None
    for i, line in enumerate(lines):
        fields = re.split(r'\t+', line.rstrip())
        if not fields or fields[0] == '':
            continue
        key = fields[0].strip()
        if key == 'label' and len(fields) > 1:
            header['label'] = fields[1].strip()
        elif key == 'ra' and len(fields) > 1:
            header['ra_deg'] = _sexed_ra_hours(fields[1].strip())
        elif key == 'dec' and len(fields) > 1:
            header['dec_deg'] = _sexed_dec_deg(fields[1].strip())
        elif key in ('pa', 'scale', 'arc2mm') and len(fields) > 1:
            header[key] = float(fields[1].strip())
        elif key in ('grism', 'filter') and len(fields) > 1:
            header[key] = fields[1].strip()
        elif key == 'corners' and len(fields) >= 5:
            header['corners'] = [float(v) for v in fields[1:5]]
        elif key == 'slit' and slit_hdr_idx is None:
            slit_hdr_idx = i

    if slit_hdr_idx is None:
        raise PypeItError(f'No slit table found in mask file {path}')

    rows = []
    for line in lines[slit_hdr_idx + 1:]:
        fields = re.split(r'\t+', line.rstrip())
        if len(fields) < 14 or not re.match(r'^\s*\d+\s*$', fields[0]):
            continue
        try:
            rows.append(dict(
                slit=int(fields[0]),
                ra_deg=_sexed_ra_hours(fields[1].strip()),
                dec_deg=_sexed_dec_deg(fields[2].strip()),
                x_mm=float(fields[3]),
                y_mm=float(fields[4]),
                target=int(fields[5]),
                object=fields[6].strip(),
                type=fields[7].strip(),
                height_mm=float(fields[8]),
                width_mm=float(fields[9]),
                offset_mm=float(fields[10]),
                theta_deg=float(fields[11]),
                bbox=[float(v) for v in fields[12].split()],
                polygon=[float(v) for v in fields[13].split()],
            ))
        except (ValueError, IndexError):
            continue

    if len(rows) == 0:
        raise PypeItError(f'No slit rows parsed from mask file {path}')

    slits = Table()
    slits['slit'] = np.array([r['slit'] for r in rows], dtype=int)
    slits['ra_deg'] = np.array([r['ra_deg'] for r in rows], dtype=float)
    slits['dec_deg'] = np.array([r['dec_deg'] for r in rows], dtype=float)
    slits['x_mm'] = np.array([r['x_mm'] for r in rows], dtype=float)
    slits['y_mm'] = np.array([r['y_mm'] for r in rows], dtype=float)
    slits['target'] = np.array([r['target'] for r in rows], dtype=int)
    slits['object'] = np.array([r['object'] for r in rows], dtype=object)
    slits['type'] = np.array([r['type'] for r in rows], dtype=object)
    slits['height_mm'] = np.array([r['height_mm'] for r in rows], dtype=float)
    slits['width_mm'] = np.array([r['width_mm'] for r in rows], dtype=float)
    slits['offset_mm'] = np.array([r['offset_mm'] for r in rows], dtype=float)
    slits['theta_deg'] = np.array([r['theta_deg'] for r in rows], dtype=float)
    slits['bbox'] = np.array([r['bbox'] for r in rows], dtype=float)
    slits['polygon'] = np.array([r['polygon'] for r in rows], dtype=float)
    return header, slits
