"""
Utilities for loading and assembling wavelength images from PypeIt calibration
files for display in the quicklook plugin.

The primary entry points are:

* :func:`find_triplets` — discover matching WaveCalib/Tilts/Slits file triplets
  in a Calibrations directory.
* :func:`build_waveimg_mosaic` — load or compute a per-MSC wavelength image for
  each triplet and concatenate them into a single array whose spatial axis
  matches the concatenated display image produced by the instrument class.
"""

from __future__ import annotations

import glob
import logging
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# Type alias for a calibration triplet: (suffix, wv_file, tilt_file, slit_file)
Triplet = Tuple[str, str, str, str]


def _calib_suffix(filename: str) -> str:
    """Return the ``_{setup}_{calib_id}_{det}`` suffix from a calibration filename.

    Parameters
    ----------
    filename : str
        Basename or full path to a calibration file whose name begins with
        ``WaveCalib``, ``Tilts``, or ``Slits``.

    Returns
    -------
    str
        The suffix string (e.g. ``"_A_01_MSC01"``), or an empty string if the
        pattern is not found.
    """
    m = re.match(r'(?:WaveCalib|Tilts|Slits)(_[^.]+)', os.path.basename(filename))
    return m.group(1) if m else ""


def find_triplets(cal_path: str) -> List[Triplet]:
    """Find all matched WaveCalib/Tilts/Slits file triplets in *cal_path*.

    Parameters
    ----------
    cal_path : str
        Path to a ``Calibrations`` directory produced by ``run_pypeit``.

    Returns
    -------
    list of tuple
        Each element is ``(suffix, wv_file, tilt_file, slit_file)`` sorted by
        suffix so multi-MSC instruments are assembled in detector order.  Returns
        an empty list if no complete triplets are found.
    """
    wv_files   = sorted(glob.glob(os.path.join(cal_path, "WaveCalib_*.fits*")))
    tilt_files = sorted(glob.glob(os.path.join(cal_path, "Tilts_*.fits*")))
    slit_files = sorted(glob.glob(os.path.join(cal_path, "Slits_*.fits*")))

    tilt_by_suffix = {_calib_suffix(f): f for f in tilt_files}
    slit_by_suffix = {_calib_suffix(f): f for f in slit_files}

    triplets = sorted(
        [
            (sfx, wv_file, tilt_by_suffix[sfx], slit_by_suffix[sfx])
            for wv_file in wv_files
            if (sfx := _calib_suffix(wv_file)) in tilt_by_suffix
            and sfx in slit_by_suffix
        ],
        key=lambda t: t[0],
    )
    return triplets


def _build_one_mosaic(
    sfx: str,
    wv_file: str,
    tilt_file: str,
    slit_file: str,
    cal_path: str,
    log: Optional[logging.Logger] = None,
) -> Tuple[str, np.ndarray, list]:
    """Load or compute the wavelength image for a single mosaic.

    Checks for a pre-computed ``WaveImage{suffix}.fits`` file alongside the
    other calibration files.  If found it is loaded directly, which is much
    faster than rebuilding from WaveCalib + Tilts + Slits.  Falls back to
    computing on the fly when the file is absent.

    Parameters
    ----------
    sfx : str
        Calibration key suffix shared by all files in this triplet.
    wv_file : str
        Path to the ``WaveCalib`` FITS file.
    tilt_file : str
        Path to the ``Tilts`` FITS file.
    slit_file : str
        Path to the ``Slits`` FITS file.
    cal_path : str
        Calibrations directory; used to locate the pre-computed WaveImage file.
    log : logging.Logger, optional
        Logger instance.  Falls back to the module-level logger when omitted.

    Returns
    -------
    tuple
        ``(sfx, waveimg, rms_list)`` where *waveimg* is a float32 ndarray and
        *rms_list* is a list of ``(spat_id, rms_or_None)`` pairs.
    """
    from pypeit.wavecalib import WaveCalib
    from pypeit.wavetilts import WaveTilts
    from pypeit.slittrace import SlitTraceSet
    from astropy.io import fits as astropy_fits

    log = log or logger

    precomputed = os.path.join(cal_path, f"WaveImage{sfx}.fits")
    wvcalib = WaveCalib.from_file(wv_file)

    if os.path.exists(precomputed):
        log.info(f"Loading pre-computed wavelength image: {precomputed}")
        with astropy_fits.open(precomputed) as hdul:
            waveimg = hdul[0].data.astype(np.float32)
    else:
        log.info(f"No pre-computed wavelength image for {sfx}; building on the fly.")
        wavetilts = WaveTilts.from_file(tilt_file)
        slits = SlitTraceSet.from_file(slit_file)
        slitmask = slits.slit_img()
        tilts = wavetilts.fit2tiltimg(slitmask, flexure=wavetilts.spat_flexure)
        waveimg = wvcalib.build_waveimg(tilts, slits).astype(np.float32)

    rms_list = []
    if wvcalib.spat_ids is not None and wvcalib.wv_fits is not None:
        for spat_id, wvfit in zip(wvcalib.spat_ids, wvcalib.wv_fits):
            rms = None if (wvfit is None or wvfit.rms is None) else wvfit.rms
            rms_list.append((int(spat_id), rms))

    return sfx, waveimg, rms_list


def build_waveimg_mosaic(
    cal_path: str,
    matched_triplets: List[Triplet],
    log: Optional[logging.Logger] = None,
) -> Tuple[np.ndarray, list]:
    """Build the full-mosaic wavelength image from a list of calibration triplets.

    Each triplet is processed concurrently (one thread per MSC), then the
    resulting per-MSC arrays are padded to a common nspec and concatenated along
    the spatial axis, matching the display image layout produced by the
    instrument's ``get_display_image`` method.

    Parameters
    ----------
    cal_path : str
        Calibrations directory; forwarded to :func:`_build_one_mosaic`.
    matched_triplets : list of tuple
        Sorted list of ``(suffix, wv_file, tilt_file, slit_file)`` as returned
        by :func:`find_triplets`.
    log : logging.Logger, optional
        Logger instance.

    Returns
    -------
    waveimg : numpy.ndarray
        Float32 wavelength image in the same coordinate frame as the display
        image.
    rms_data : list of tuple
        Flat list of ``(spat_id, rms_or_None)`` pairs from all mosaics.
    """
    log = log or logger

    results: dict = {}
    with ThreadPoolExecutor(max_workers=len(matched_triplets)) as executor:
        futures = {
            executor.submit(_build_one_mosaic, sfx, wv, tilt, slit, cal_path, log): sfx
            for sfx, wv, tilt, slit in matched_triplets
        }
        for future in as_completed(futures):
            sfx, waveimg, rms_list = future.result()
            results[sfx] = (waveimg, rms_list)

    # Reassemble in sorted order to match the display image spatial layout.
    mosaic_waveimgs = []
    rms_data = []
    for sfx, _, _, _ in matched_triplets:
        waveimg, rms_list = results[sfx]
        mosaic_waveimgs.append(waveimg)
        rms_data.extend(rms_list)

    # Per-mosaic waveimgs may differ slightly in nspec due to rotation —
    # pad to a common height before concatenating along the spatial axis,
    # mirroring the display image assembly.
    if len(mosaic_waveimgs) > 1:
        max_nspec = max(img.shape[0] for img in mosaic_waveimgs)
        padded = []
        for img in mosaic_waveimgs:
            deficit = max_nspec - img.shape[0]
            if deficit > 0:
                img = np.pad(img, ((0, deficit), (0, 0)))
            padded.append(img)
        waveimg = np.concatenate(padded, axis=1)
    else:
        waveimg = mosaic_waveimgs[0]

    return waveimg, rms_data
