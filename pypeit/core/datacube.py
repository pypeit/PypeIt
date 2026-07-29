"""
Module containing routines used by 3D datacubes.

.. include:: ../include/links.rst
"""

import os
# import line_profiler
from pathlib import Path
from typing import TYPE_CHECKING

from astropy import wcs, units
from astropy.coordinates import AltAz, SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats, SigmaClip
from fast_histogram import histogramdd
from IPython import embed
import scipy.optimize as opt
from scipy import signal, ndimage
from scipy.interpolate import interp1d, griddata
from scipy.spatial import QhullError
import numpy as np

# NOTE: photutils is an optional dependency
try:
    from photutils.detection import DAOStarFinder
except ModuleNotFoundError:
    DAOStarFinder = None

from pypeit import log
from pypeit import PypeItError
from pypeit import utils
from pypeit import slittrace
from pypeit import spec2dobj
from pypeit import specobj
from pypeit import specobjs
from pypeit.core import coadd, sampling
from pypeit.display import display
from pypeit.images.imagebitmask import ImageBitMaskArray
from pypeit.onespec import OneSpec
from pypeit.spectrographs.util import load_spectrograph

if TYPE_CHECKING:
    from pypeit.spectrographs.spectrograph import Spectrograph


def resample_spec_to_grid(wave, flux, ivar, wave_grid, min_good=2, min_frac=0.5):
    """
    Resample one spectrum's flux and inverse variance onto a wavelength grid.

    The flux and error are resampled with :class:`~pypeit.sampling.Resample`,
    PypeIt's flux-conserving spectral resampler, rather than a plain linear
    interpolation.  The error (``1/sqrt(ivar)``) is propagated by ``Resample``
    in quadrature and converted back to inverse variance, and pixels whose
    native coverage falls short of ``min_frac`` are masked.  Samples that fall
    outside the spectrum's native wavelength coverage are left at zero.

    Parameters
    ----------
    wave, flux, ivar : `numpy.ndarray`_
        Native wavelength, flux, and inverse-variance arrays for a single
        spectrum (e.g. one fiber).  ``ivar`` may be ``None``, in which case
        the returned inverse variance is all zeros.
    wave_grid : `numpy.ndarray`_
        Monotonically increasing target wavelength grid.
    min_good : :obj:`int`, optional
        Minimum number of finite samples with positive wavelength required
        for the spectrum to contribute; below this, all-zero arrays are
        returned.
    min_frac : :obj:`float`, optional
        Minimum covering fraction (from ``Resample.outf``) for an output
        pixel to be flagged as covered and retain its resampled values.

    Returns
    -------
    flux_grid : `numpy.ndarray`_
        Flux resampled onto ``wave_grid`` (zero outside native coverage).
    ivar_grid : `numpy.ndarray`_
        Inverse variance resampled onto ``wave_grid`` (zero outside coverage
        or where the native inverse variance was non-positive).
    covered : `numpy.ndarray`_
        Boolean mask, ``True`` where ``wave_grid`` is covered by the
        spectrum's native wavelength range.
    """
    wave_grid = np.asarray(wave_grid, dtype=float)
    flux_grid = np.zeros(wave_grid.shape, dtype=float)
    ivar_grid = np.zeros(wave_grid.shape, dtype=float)
    covered = np.zeros(wave_grid.shape, dtype=bool)

    wave = np.asarray(wave, dtype=float)
    flux = np.asarray(flux, dtype=float)
    good = (wave > 0) & np.isfinite(flux)
    if np.count_nonzero(good) < min_good:
        return flux_grid, ivar_grid, covered

    srt = np.argsort(wave[good])
    w_s = wave[good][srt]
    f_s = flux[good][srt]
    if ivar is not None:
        iv_s = np.asarray(ivar, dtype=float)[good][srt]
        err = np.sqrt(utils.inverse(iv_s))
        # Mask native pixels with no inverse variance so they carry no weight.
        mask = np.logical_not(iv_s > 0)
    else:
        err = None
        mask = None

    # newLog=False: the wavelength grid is linear, so output pixel borders are
    # arithmetic (not geometric) midpoints of wave_grid.
    r = sampling.Resample(f_s, e=err, mask=mask, x=w_s, newx=wave_grid,
                          inLog=False, newLog=False, conserve=False, ext_value=0.0)

    covered = r.outf > min_frac
    flux_grid = np.where(covered, r.outy, 0.0)
    if err is not None:
        ivar_grid = np.where(covered, utils.inverse(r.oute)**2, 0.0)
    return flux_grid, ivar_grid, covered


def build_cube_from_spec1d(spec1d_file, spectrograph, targetx, targety,
                           boxcar=False, spatial_scale=0.27, method='linear',
                           output=None):
    """Build a single fiber-IFU datacube from one spec1d file.

    Reads the already-extracted 1D fiber spectra from a PypeIt spec1d file
    and builds a datacube.  Uses OPT (optimal) extraction by default, or
    BOX (boxcar) extraction if ``boxcar`` is True.  Sky is already
    subtracted by the pipeline, so no additional sky subtraction is
    performed here.

    The spectrograph must implement :func:`get_fiber_metadata`,
    :func:`get_science_fiber_layout_indices`, :func:`load_sky_layout`, and
    :func:`ifu_sky_wcs` to map fibers to sky positions and build the WCS.

    Parameters
    ----------
    spec1d_file : :obj:`str`
        Path to the spec1d FITS file.
    spectrograph : :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
        Spectrograph instance.
    targetx : `numpy.ndarray`_
        Fiber x positions on sky (arcsec).
    targety : `numpy.ndarray`_
        Fiber y positions on sky (arcsec).
    boxcar : :obj:`bool`, optional
        Use BOX (boxcar) extraction columns instead of OPT (optimal).
    spatial_scale : :obj:`float`, optional
        Output spatial pixel scale in arcsec.
    method : :obj:`str`, optional
        Spatial interpolation method (``'nearest'``, ``'linear'``, or
        ``'cubic'``).
    output : :obj:`str`, optional
        Output FITS filename.  If None, a name is auto-generated from the
        input file.
    """
    sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file)
    if sobjs.nobj == 0:
        log.warning(f"No objects in {Path(spec1d_file).name}, skipping")
        return

    # Choose extraction type
    prefix = 'BOX' if boxcar else 'OPT'
    log.info(f"  Using {prefix} extraction from spec1d")

    # ------------------------------------------------------------------
    # Step 1: Organize fiber spectra by detector
    # ------------------------------------------------------------------
    det_fiber_data = {}

    # Process each detector present in the spec1d file.
    for det_name in sorted(np.unique(sobjs.DET)):
        det_sobjs = sobjs[sobjs.DET == det_name]
        if len(det_sobjs) == 0:
            log.warning(f"  No objects for {det_name}, skipping")
            continue

        nfibers = len(det_sobjs)
        log.info(f"  {det_name}: {nfibers} fibers from spec1d")

        # Read extracted spectra
        wave_key = f'{prefix}_WAVE'
        flux_key = f'{prefix}_COUNTS'
        ivar_key = f'{prefix}_COUNTS_IVAR'

        # Determine spectral length from first object
        nspec = getattr(det_sobjs[0], wave_key).shape[0]

        fiber_flux = np.zeros((nfibers, nspec))
        fiber_ivar = np.zeros((nfibers, nspec))
        fiber_wave = np.zeros((nfibers, nspec))
        spat_ids = np.zeros(nfibers, dtype=int)
        slit_centers = np.zeros(nfibers, dtype=float)

        for i, sobj in enumerate(det_sobjs):
            fiber_wave[i] = getattr(sobj, wave_key)
            fiber_flux[i] = getattr(sobj, flux_key)
            fiber_ivar[i] = getattr(sobj, ivar_key)
            spat_ids[i] = sobj.SLITID
            slit_centers[i] = sobj.SPAT_PIXPOS

        # Get fiber metadata (pass float centers for sub-pixel accuracy)
        fiber_meta = spectrograph.get_fiber_metadata(
            int(det_name.replace('DET', '')), spat_ids,
            slit_centers=slit_centers)
        if fiber_meta is None:
            raise PypeItError(
                f"{spectrograph.name} does not implement get_fiber_metadata(); "
                "it is required to map fibers to sky positions.")

        det_fiber_data[det_name] = {
            'flux': fiber_flux,
            'ivar': fiber_ivar,
            'wave': fiber_wave,
            'fiber_meta': fiber_meta,
        }

    if len(det_fiber_data) == 0:
        log.warning(f"No detector data from "
                    f"{Path(spec1d_file).name}, skipping")
        return

    # ------------------------------------------------------------------
    # Build the datacube
    # ------------------------------------------------------------------
    with fits.open(spec1d_file) as hdu:
        raw_hdr = hdu[0].header

    build_cube_common(det_fiber_data, spectrograph, targetx, targety,
                      raw_hdr, spec1d_file, spatial_scale=spatial_scale,
                      method=method, output=output)


def build_cube_common(det_fiber_data, spectrograph, targetx, targety,
                      raw_hdr, input_file, spatial_scale=0.27,
                      method='linear', output=None):
    """Shared steps for building a fiber-IFU datacube from fiber spectra.

    The input fiber spectra are assumed to be already sky-subtracted
    (as produced by the PypeIt pipeline in the spec1d files).

    Parameters
    ----------
    det_fiber_data : :obj:`dict`
        Per-detector fiber data. Keys are detector names (e.g. ``'DET01'``),
        values are dicts with keys ``'flux'``, ``'ivar'``, ``'wave'``,
        ``'fiber_meta'``.
    spectrograph : :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
        Spectrograph instance.
    targetx : `numpy.ndarray`_
        Fiber x positions on sky (arcsec).
    targety : `numpy.ndarray`_
        Fiber y positions on sky (arcsec).
    raw_hdr : `astropy.io.fits.Header`_
        Primary header from the input file.
    input_file : :obj:`str`
        Path to the input file (for generating the output filename).
    spatial_scale : :obj:`float`, optional
        Output spatial pixel scale in arcsec.
    method : :obj:`str`, optional
        Spatial interpolation method (``'nearest'``, ``'linear'``, or
        ``'cubic'``).
    output : :obj:`str`, optional
        Output FITS filename.  If None, a name is auto-generated from the
        input file.
    """
    # ------------------------------------------------------------------
    # Identify sky vs. science fibers (flux is already sky-subtracted)
    # ------------------------------------------------------------------
    for det_name, data in det_fiber_data.items():
        fiber_meta = data['fiber_meta']
        sky_mask = fiber_meta['fiber_type'] == 'SKY'
        n_sky = np.sum(sky_mask)
        n_sci = np.sum(~sky_mask)
        log.info(f"  {det_name}: {n_sky} sky fibers, {n_sci} science fibers")

        data['sky_mask'] = sky_mask
        data['sci_mask'] = ~sky_mask

    # ------------------------------------------------------------------
    # Step 4: Wavelength linearization
    # ------------------------------------------------------------------
    # Find global wavelength range across all detectors
    all_waves = []
    for data in det_fiber_data.values():
        wave = data['wave']
        valid = wave > 0
        if np.any(valid):
            all_waves.extend([np.min(wave[valid]), np.max(wave[valid])])

    if len(all_waves) == 0:
        log.error("No valid wavelength data found. "
                  "Check that input files contain extracted fiber spectra.")
        return
    wave_min = min(all_waves[::2])
    wave_max = max(all_waves[1::2])

    # Use median dispersion for wavelength step
    dispersions = []
    for data in det_fiber_data.values():
        wave = data['wave']
        for i in range(wave.shape[0]):
            valid = wave[i] > 0
            if np.sum(valid) > 10:
                dw = np.diff(wave[i, valid])
                dw = dw[dw > 0]
                if len(dw) > 0:
                    dispersions.append(np.median(dw))
                break  # One fiber is enough per detector

    dwv = np.median(dispersions)
    n_wave = int(np.ceil((wave_max - wave_min) / dwv)) + 1
    wave_grid = np.linspace(wave_min, wave_min + (n_wave - 1) * dwv, n_wave)
    log.info(f"Wavelength grid: {wave_min:.1f} to {wave_grid[-1]:.1f} A, "
              f"dw={dwv:.3f} A, {n_wave} pixels")

    # Resample each fiber onto the common wavelength grid using the
    # flux-conserving resampler (shared with the 1D fiber extractor).
    for data in det_fiber_data.values():
        nfibers = data['flux'].shape[0]
        flux_resamp = np.zeros((nfibers, n_wave))
        ivar_resamp = np.zeros((nfibers, n_wave))

        for i in range(nfibers):
            flux_resamp[i], ivar_resamp[i], _ = resample_spec_to_grid(
                data['wave'][i], data['flux'][i], data['ivar'][i], wave_grid,
                min_good=10)

        data['flux_resamp'] = flux_resamp
        data['ivar_resamp'] = ivar_resamp

    # ------------------------------------------------------------------
    # Step 5: Combine both detectors
    # ------------------------------------------------------------------
    # Map science fibers to layout file positions
    sci_flux_list = []
    sci_ivar_list = []
    layout_idx_list = []

    for det_name in sorted(det_fiber_data.keys()):
        data = det_fiber_data[det_name]
        det_num = int(det_name.replace('DET', ''))

        sci_mask = data['sci_mask']
        fiber_meta = data['fiber_meta']
        layout_indices = spectrograph.get_science_fiber_layout_indices(
            det_num, fiber_meta['fiber_id'], fiber_meta['fiber_type'])

        sci_flux = data['flux_resamp'][sci_mask]
        sci_ivar = data['ivar_resamp'][sci_mask]
        sci_layout = layout_indices[sci_mask]

        # Remove any fibers with invalid layout indices
        valid = sci_layout >= 0
        sci_flux_list.append(sci_flux[valid])
        sci_ivar_list.append(sci_ivar[valid])
        layout_idx_list.append(sci_layout[valid])

    combined_flux = np.vstack(sci_flux_list)
    combined_ivar = np.vstack(sci_ivar_list)
    combined_layout = np.concatenate(layout_idx_list)

    n_sci_fibers = combined_flux.shape[0]
    log.info(f"Combined {n_sci_fibers} science fibers from "
              f"{len(det_fiber_data)} detector(s)")

    # Trim wavelength range to where a reasonable fraction of fibers
    # have valid data (avoids degenerate interpolation at edges)
    n_valid = np.sum((combined_flux != 0) | (combined_ivar > 0), axis=0)
    min_fibers = max(10, int(0.10 * n_sci_fibers))
    good_wave = n_valid >= min_fibers
    if not np.all(good_wave):
        first = np.argmax(good_wave)
        last = n_wave - 1 - np.argmax(good_wave[::-1])
        log.info(f"Trimming wavelength range: slices {first}-{last} of "
                  f"{n_wave} (>={min_fibers} fibers required)")
        wave_grid = wave_grid[first:last + 1]
        combined_flux = combined_flux[:, first:last + 1]
        combined_ivar = combined_ivar[:, first:last + 1]
        n_wave = len(wave_grid)

    # Load fiber sky positions
    fiber_x = targetx[combined_layout]
    fiber_y = targety[combined_layout]

    # ------------------------------------------------------------------
    # Step 6: Build datacube via spatial interpolation
    # ------------------------------------------------------------------
    # Fiber positions are in arcsec (from the spectrograph sky layout); the
    # output grid and the WCS share the same arcsec spatial scale.
    scl = spatial_scale

    # Compute grid dimensions from the fiber positions
    x_min, x_max = np.min(fiber_x), np.max(fiber_x)
    y_min, y_max = np.min(fiber_y), np.max(fiber_y)

    # Add small padding
    pad = scl
    x_min -= pad
    x_max += pad
    y_min -= pad
    y_max += pad

    nx = int(np.ceil((x_max - x_min) / scl)) + 1
    ny = int(np.ceil((y_max - y_min) / scl)) + 1

    log.info(f"Output cube dimensions: {nx} x {ny} x {n_wave}")

    # Build regular grid
    x_grid = np.linspace(x_min, x_min + (nx - 1) * scl, nx)
    y_grid = np.linspace(y_min, y_min + (ny - 1) * scl, ny)
    grid_x, grid_y = np.meshgrid(x_grid, y_grid, indexing='ij')

    # Interpolate at each wavelength
    points = np.column_stack([fiber_x, fiber_y])
    cube = np.zeros((nx, ny, n_wave), dtype=np.float32)
    var_cube = np.zeros((nx, ny, n_wave), dtype=np.float32)

    log.info(f"Interpolating {n_wave} wavelength slices using "
              f"method='{method}'...")
    for k in range(n_wave):
        if k % 500 == 0:
            log.info(f"  Wavelength slice {k}/{n_wave}")

        flux_slice = combined_flux[:, k]
        ivar_slice = combined_ivar[:, k]

        # Only interpolate fibers with valid data
        good = (flux_slice != 0) | (ivar_slice > 0)
        if np.sum(good) < 4:
            continue

        try:
            cube[:, :, k] = griddata(
                points[good], flux_slice[good], (grid_x, grid_y),
                method=method, fill_value=0.0)
        except QhullError:
            # Fall back to nearest-neighbor when valid points are
            # degenerate (e.g. collinear at spectral edges)
            cube[:, :, k] = griddata(
                points[good], flux_slice[good], (grid_x, grid_y),
                method='nearest', fill_value=0.0)

        # Interpolate variance
        var_slice = np.where(ivar_slice > 0, 1.0 / ivar_slice, 0.0)
        if np.any(var_slice[good] > 0):
            try:
                var_cube[:, :, k] = griddata(
                    points[good], var_slice[good], (grid_x, grid_y),
                    method=method, fill_value=0.0)
            except QhullError:
                var_cube[:, :, k] = griddata(
                    points[good], var_slice[good], (grid_x, grid_y),
                    method='nearest', fill_value=0.0)

    # ------------------------------------------------------------------
    # Step 7: Build WCS and write output
    # ------------------------------------------------------------------
    # Pointing and celestial CD matrix; the spectrograph owns the IFU sky
    # convention (see Spectrograph.ifu_sky_wcs) and the IFU mode metadata.
    coord, cd = spectrograph.ifu_sky_wcs(raw_hdr, scl)
    (cd11, cd12), (cd21, cd22) = cd
    ifu_meta = spectrograph.get_ifu_datacube_meta(raw_hdr)

    w = wcs.WCS(naxis=3)
    w.wcs.equinox = raw_hdr.get('EQUINOX', 2000.0)
    w.wcs.name = ifu_meta['name']
    w.wcs.radesys = 'ICRS'
    w.wcs.cname = ['RA', 'DEC', 'Wavelength']
    w.wcs.cunit = [units.degree, units.degree, units.Angstrom]
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN', 'WAVE']
    w.wcs.crval = [coord.ra.degree, coord.dec.degree, wave_grid[0]]
    w.wcs.crpix = [nx / 2.0, ny / 2.0, 1.0]
    w.wcs.cd = np.array([[cd11, cd12, 0.0],
                         [cd21, cd22, 0.0],
                         [0.0, 0.0, dwv]])
    w.wcs.lonpole = 180.0
    w.wcs.latpole = 0.0

    # Build output FITS.  Instrument identity comes from the spectrograph.
    hdr = w.to_header()
    hdr['INSTRUME'] = (spectrograph.camera, 'Instrument')
    hdr['TELESCOP'] = (spectrograph.telescope['name'], 'Telescope')
    hdr['IFUMODE'] = (ifu_meta['mode'], 'IFU mode')
    hdr['NFIBERS'] = (n_sci_fibers, 'Number of science fibers')
    hdr['SPATSCL'] = (scl, 'Spatial pixel scale [arcsec]')
    hdr['WAVEMIN'] = (wave_grid[0], 'Minimum wavelength [Angstrom]')
    hdr['WAVEMAX'] = (wave_grid[-1], 'Maximum wavelength [Angstrom]')
    hdr['WAVESTP'] = (dwv, 'Wavelength step [Angstrom]')
    hdr['INTERP'] = (method, 'Spatial interpolation method')

    # Copy useful keywords from raw header
    for key in ['OBJECT', 'EXPTIME', 'DATE-OBS', 'DISPERSE', 'FILTER']:
        if key in raw_hdr:
            hdr[key] = raw_hdr[key]

    # Output filename
    if output is not None:
        outfile = output
    else:
        base = Path(input_file).stem
        if 'spec1d_' in base:
            base = base.replace('spec1d_', 'cube_')
        outfile = base + '.fits'

    # Transpose from numpy (nx, ny, n_wave) to FITS order (n_wave, ny, nx)
    # so that NAXIS1=nx(RA), NAXIS2=ny(DEC), NAXIS3=n_wave(WAVE)
    cube = np.transpose(cube, (2, 1, 0))
    var_cube = np.transpose(var_cube, (2, 1, 0))

    primary = fits.PrimaryHDU(header=fits.Header())
    primary.header['AUTHOR'] = 'PypeIt'
    flux_hdu = fits.ImageHDU(data=cube, header=hdr, name='FLUX')
    var_hdu = fits.ImageHDU(data=var_cube, header=hdr, name='VAR')

    hdulist = fits.HDUList([primary, flux_hdu, var_hdu])
    hdulist.writeto(outfile, overwrite=True)
    log.info(f"Wrote datacube to {outfile}")
    log.info(f"Cube shape: {cube.shape}")


def project_to_sky(x_arcsec, y_arcsec, raw_hdr, spectrograph):
    """Project instrument-frame fiber offsets to sky coordinates.

    Uses the spectrograph's shared TAN/POSANG WCS convention
    (:meth:`~pypeit.spectrographs.spectrograph.Spectrograph.ifu_sky_wcs`) so
    that the extracted hex view aligns with the matching datacube.

    Parameters
    ----------
    x_arcsec, y_arcsec : `numpy.ndarray`_
        1D arrays of fiber offsets in instrument-frame arcseconds.
    raw_hdr : `astropy.io.fits.Header`_
        Primary header from the spec1d file, providing ``RA``, ``DEC``,
        and (optionally) ``POSANG``.
    spectrograph : :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
        Provides the ``ifu_sky_wcs`` reference-coordinate/CD-matrix helper.

    Returns
    -------
    ra, dec : `numpy.ndarray`_
        Fiber RA/Dec in degrees, same shape as the inputs.
    """
    # The inputs are already in arcsec and are fed directly to pixel_to_world
    # below as the WCS pixel coordinates, so the WCS scale is fixed at 1 arcsec
    # per unit step: a fiber offset of N arcsec maps to N arcsec on sky.  This
    # is a units identity, not the output sampling (cf. build_cube_common, which
    # passes the real spatial_scale).
    coord, cd = spectrograph.ifu_sky_wcs(raw_hdr, 1.0)

    w = wcs.WCS(naxis=2)
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.crval = [coord.ra.degree, coord.dec.degree]
    # crpix is FITS 1-indexed; pixel_to_world is 0-indexed.  Setting crpix=1
    # makes pixel (0, 0) == reference pixel == crval.
    w.wcs.crpix = [1.0, 1.0]
    w.wcs.cd = cd
    w.wcs.lonpole = 180.0
    w.wcs.latpole = 0.0

    sky = w.pixel_to_world(np.asarray(x_arcsec, dtype=float),
                           np.asarray(y_arcsec, dtype=float))
    return np.atleast_1d(sky.ra.degree), np.atleast_1d(sky.dec.degree)


def resample_and_combine(waves, fluxes, ivars):
    """Resample selected fibers onto a common grid and combine.

    Parameters
    ----------
    waves, fluxes, ivars : :obj:`list` of `numpy.ndarray`_
        Native per-fiber wavelength, flux, and inverse-variance arrays
        (parallel lists).  Wavelengths must be positive; the function
        sorts each fiber internally so callers need not pre-sort.

    Returns
    -------
    wave_out : `numpy.ndarray`_
        Common output wavelength grid (Angstrom).
    flux_out : `numpy.ndarray`_
        Summed flux across fibers, NaN-safe with count rescaling.
    ivar_out : `numpy.ndarray`_
        Inverse variance of the combined flux.

    Raises
    ------
    PypeItError
        If the selected fibers do not share any overlapping wavelength range.
    """
    valid = [(w, f, iv) for w, f, iv in zip(waves, fluxes, ivars)
             if np.any(w > 0)]
    if not valid:
        raise PypeItError("No valid wavelength data in selected fibers")

    # Common range = intersection of fibers' native ranges.
    wave_min = max(np.min(w[w > 0]) for w, _, _ in valid)
    wave_max = min(np.max(w[w > 0]) for w, _, _ in valid)
    if wave_max <= wave_min:
        raise PypeItError("Selected fibers have no overlapping wavelength "
                          f"range (min={wave_min:.2f}, max={wave_max:.2f})")

    # Median pixel width across selected fibers.
    diffs = []
    for w, _, _ in valid:
        good = w > 0
        if good.sum() > 1:
            d = np.diff(w[good])
            d = d[d > 0]
            if d.size:
                diffs.append(np.median(d))
    if not diffs:
        raise PypeItError("Could not determine wavelength dispersion from "
                          "selected fibers")
    dwv = np.median(diffs)
    n_wave = int(round((wave_max - wave_min) / dwv)) + 1
    wave_out = np.linspace(wave_min, wave_max, n_wave)

    n_fib = len(valid)
    flux_resamp = np.zeros((n_fib, n_wave))
    ivar_resamp = np.zeros((n_fib, n_wave))
    have_flux = np.zeros((n_fib, n_wave), dtype=bool)

    for i, (w, f, iv) in enumerate(valid):
        flux_resamp[i], ivar_resamp[i], have_flux[i] = \
            resample_spec_to_grid(w, f, iv, wave_out)

    # Coadd: sum flux (rescaling for partial wavelength coverage) and sum the
    # per-fiber variances (utils.inverse is zero where ivar <= 0, so masked
    # pixels drop out of the sum automatically).
    n_good = have_flux.sum(axis=0)
    rescale = n_fib / np.maximum(n_good, 1)
    flux_out = np.sum(np.where(have_flux, flux_resamp, 0.0), axis=0) * rescale

    var_out = np.sum(utils.inverse(ivar_resamp), axis=0)
    ivar_out = utils.inverse(var_out)

    return wave_out, flux_out, ivar_out


def load_fibers(sobjs, spectrograph, targetx, targety, prefix='OPT'):
    """Build per-fiber records from a SpecObjs object.

    Parameters
    ----------
    sobjs : :class:`~pypeit.specobjs.SpecObjs`-like
        Iterable of :class:`~pypeit.specobj.SpecObj`.  Must expose ``DET``,
        ``SLITID``, ``SPAT_PIXPOS``, and the relevant ``{prefix}_*`` arrays.
    spectrograph : :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
        Provides ``get_fiber_metadata`` and
        ``get_science_fiber_layout_indices``.
    targetx, targety : `numpy.ndarray`_
        Layout-file fiber positions in instrument-frame arcseconds, indexed
        by science layout index.
    prefix : {'OPT', 'BOX'}, optional
        Extraction column prefix (optimal or boxcar).

    Returns
    -------
    :obj:`list` of :obj:`dict`
        One record per surviving science fiber, with keys ``wave``, ``flux``,
        ``ivar``, ``x``, ``y``, ``fiber_id``, ``fiber_type``, ``det``.

    Raises
    ------
    PypeItError
        If no science fibers survive (e.g. spec1d contains only sky fibers).
    """
    wave_key = f'{prefix}_WAVE'
    flux_key = f'{prefix}_COUNTS'
    ivar_key = f'{prefix}_COUNTS_IVAR'

    fibers = []
    for det_name in sorted(np.unique(sobjs.DET)):
        det_mask = sobjs.DET == det_name
        det_sobjs = sobjs[det_mask]
        if len(det_sobjs) == 0:
            continue
        det_num = int(det_name.replace('DET', ''))
        spat_ids = np.array([s.SLITID for s in det_sobjs], dtype=int)
        slit_centers = np.array([s.SPAT_PIXPOS for s in det_sobjs],
                                dtype=float)
        meta = spectrograph.get_fiber_metadata(det_num, spat_ids,
                                               slit_centers=slit_centers)
        if meta is None:
            raise PypeItError(
                f"{spectrograph.name} does not implement get_fiber_metadata(); "
                "it is required to map fibers to sky positions.")
        fiber_ids = meta['fiber_id']
        fiber_types = np.asarray(meta['fiber_type'])
        layout = spectrograph.get_science_fiber_layout_indices(
            det_num, fiber_ids, fiber_types)
        for i, sobj in enumerate(det_sobjs):
            if fiber_types[i] == 'SKY' or layout[i] < 0:
                continue
            fibers.append({
                'wave': np.asarray(getattr(sobj, wave_key), dtype=float),
                'flux': np.asarray(getattr(sobj, flux_key), dtype=float),
                'ivar': np.asarray(getattr(sobj, ivar_key), dtype=float),
                'x': float(targetx[layout[i]]),
                'y': float(targety[layout[i]]),
                'fiber_id': int(fiber_ids[i]),
                'fiber_type': str(fiber_types[i]),
                'det': det_name,
            })

    if not fibers:
        raise PypeItError("No science fibers found in spec1d input")
    return fibers


def write_onespec(wave, flux, ivar, raw_hdr, pyp_spec, outfile):
    """Write a combined fiber spectrum as a :class:`~pypeit.onespec.OneSpec`.

    Parameters
    ----------
    wave, flux, ivar : `numpy.ndarray`_
        1D arrays of equal length (output wavelength grid, summed flux,
        inverse variance).
    raw_hdr : `astropy.io.fits.Header`_
        Primary header of the source spec1d, copied through to the output.
    pyp_spec : :obj:`str`
        Spectrograph short name, written as ``PYP_SPEC`` so downstream
        PypeIt scripts can re-load the file.
    outfile : :obj:`str`
        Output FITS path.
    """
    # `wave_grid_mid` is the (optional) uniformly-spaced grid used by 1D
    # coaddition; this writer is not a coadd, so leave it as None.
    one = OneSpec(wave=wave, wave_grid_mid=None, flux=flux, ivar=ivar,
                  PYP_SPEC=pyp_spec)
    # Do NOT set `one.head0` before `to_file` -- that triggers
    # `spectrograph.subheader_for_spec` which expects extra metadata keys
    # (e.g. target tables) not present in arbitrary input headers.  Pass
    # the raw header via `primary_hdr=` to copy keys through.
    one.to_file(outfile, primary_hdr=raw_hdr, overwrite=True)


def gaussian2D(tup, intflux, xo, yo, sigma_x, sigma_y, theta, offset):
    """
    Fit a 2D Gaussian function to an image.

    Args:
        tup (:obj:`tuple`):
            A two element tuple containing the x and y coordinates of each pixel
            in the image
        intflux (float):
            The Integrated flux of the 2D Gaussian
        xo (float):
            The centre of the Gaussian along the x-coordinate when z=0 (units of pixels)
        yo (float):
            The centre of the Gaussian along the y-coordinate when z=0 (units of pixels)
        sigma_x (float):
            The standard deviation in the x-direction (units of pixels)
        sigma_y (float):
            The standard deviation in the y-direction (units of pixels)
        theta (float):
            The orientation angle of the 2D Gaussian
        offset (float):
            Constant offset

    Returns:
        `numpy.ndarray`_: The 2D Gaussian evaluated at the coordinate (x, y)
    """
    # Extract the (x, y, z) coordinates of each pixel from the tuple
    (x, y) = tup
    # Ensure these are floating point
    xo = float(xo)
    yo = float(yo)
    # Account for a rotated 2D Gaussian
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    # Normalise so that the integrated flux is a parameter, instead of the amplitude
    norm = 1/(2*np.pi*np.sqrt(a*c-b*b))
    gtwod = offset + norm*intflux*np.exp(-(a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    return gtwod.ravel()


def fitGaussian2D(image, ivar=None, gpm=None, init_obj_position=None, 
                  fwhm=3.0, nsigma=5.0, mask_edge=0, median_filter=False, 
                  norm=False, pixelscale=None, verbose=False):
    """
    Fit a 2D Gaussian to an input image. It is recommended that the input image
    is scaled to a maximum value that is ~1, so that all fit parameters are of
    the same order of magnitude. Set norm=True if you do not care about the
    amplitude or integrated flux. Otherwise, make sure you scale the image by
    a known value prior to passing it into this function.

    Image coordinates are quoted as (x, y), matching the coordinates read from
    Ginga or DS9. In numpy terms, if the image has shape (ny, nx), a position
    (x, y) refers to image[y, x].

    Parameters
    ----------
    image : `numpy.ndarray`_
        A 2D input image
    ivar : `numpy.ndarray`_
        The inverse variance of the image. Optional. If not passed, the standard deviation computed
        from the image will be used to compute the inverse variance. Default is None.
    gpm : `numpy.ndarray`_, optional
        A good pixel mask. Pixels that are True are good. Default is None.
    init_obj_position : tuple, optional
        The initial guess for the object position in the image with format
        (x, y), matching Ginga or DS9 image coordinates. In numpy terms, if
        the image has shape (ny, nx), a position (x, y) refers to image[y, x].
        If set, the 2D Gaussian fit will be performed with the position
        constrainted to be within plus or minus fwhm/3 in x and y. If not set,
        the position will be determined by running DAOStarFinder on the image.
        Default is None.
    fwhm : float, optional
        The FWHM of the image in pixels. This is used to estimate the initial
        guess for the Gaussian fit, the fit bounds, and the median filter kernel
        width if median filtering is used. Default is 3.0 pixels.
    nsigma : float, optional
        The number of sigma to use when determining the threshold for the
        DAOStarFinder object, i.e. the threshold is nsigma*std_dev.  Default is 5.0.
    mask_edge : int, optional
        The number of pixels to mask at the edges of the image. Default is 4 pixels.
    median_filter : bool, optional
        If True, the object finding will be performed on a median filtered
        image with a kernel size of fwhm, instead of the image itself. Default is False.
    norm : bool, optional
        If True, the input image will be normalised to the maximum value
        of the input image.
    pixelscale : float, optional
        The plate scale of the image in arcseconds per pixel. This is only used to print the
        FWHM of the Gaussian to the screen in arcseconds.  Default is None, in which 
        case the FWHM will be printed in pixels.
    verbose : bool, optional
        If True, the DAOStarfinder properties of the brightest source will be printed to the screen. 

    Returns
    -------
    popt : `numpy.ndarray`_
       The optimum parameters of the Gaussian in the following order:
       integrated flux, x center, y center, sigma_x, sigma_y, theta, offset.
       See :func:`~pypeit.core.datacube.gaussian2D` for a more detailed
       description of the model.
    pcov : `numpy.ndarray`_
        Corresponding covariance matrix
    model : `numpy.ndarray`_
        The 2D Gaussian model evaluated at the input image pixel locations
    _init_obj_position : tuple
        If the init_obj_position input parameter is None, this will be the
        initial guess for the object position in (x, y) image coordinates,
        determined by running DAOStarFinder on the image. Otherwise it will be
        the input value.
    flux_opt : float
        The optimally extracted object flux of the brightest source in the image
    sigma_opt : float
        The optimally extracted one sigma error of the object flux of the brightest source in the image.       
    """
    _gpm = np.ones_like(image, dtype=bool) if gpm is None else gpm
    fwhm2sigma = 1.0 / (2 * np.sqrt(2 * np.log(2)))
    sigma = fwhm*fwhm2sigma
    # Normalise if requested
    wlscl = np.max(image) if norm else 1.0

    ## Find the objects
    ximg = np.tile(np.arange(image.shape[1]), (image.shape[0], 1))
    yimg = np.tile(np.arange(image.shape[0]), (image.shape[1], 1)).T
    edgemask = (ximg < mask_edge) | (ximg >= image.shape[1] - mask_edge) | \
                (yimg < mask_edge) | (yimg >= image.shape[0] - mask_edge)
    totmask = edgemask | np.logical_not(_gpm)

    if ivar is None:
        mean, median, std = sigma_clipped_stats(image[np.logical_not(totmask)], sigma=3.0)
        if std > 0:
            _ivar = np.full_like(image, 1.0/std**2)
        else:
            log.warning('Could not measure standard deviation from image.  Assuming 1.')
            _ivar = np.ones_like(image)
    else:
        _ivar = ivar

    if init_obj_position is None: 
        if DAOStarFinder is None:
            raise PypeItError(
                'Requires optional photutils (>=3.0.0) dependency to proceed.  Try to reinstall pypeit '
                'including the datacube dependencies; e.g., pip install "pypeit[datacube]".'
            )
        if median_filter:
            int_kernel = np.clip(round(fwhm), 3, None)
            if int_kernel % 2 == 0:
                int_kernel += 1 if fwhm > int_kernel else -1
            objfind_image = signal.medfilt2d(image, kernel_size=int_kernel)
            mean_objfind, median_objfind, std_objfind = sigma_clipped_stats(
                objfind_image[np.logical_not(totmask)], sigma=3.0)
            ivar_objfind = np.full_like(image, 1.0/std_objfind**2)
        else:
            objfind_image = image
            ivar_objfind = _ivar
            mean_objfind, median_objfind, std_objfind = sigma_clipped_stats(
                objfind_image[np.logical_not(totmask)], sigma=3.0)        

        # Create a border mask to exclude junk at the edges
        daofind = DAOStarFinder(
            fwhm=fwhm, threshold=nsigma, sharpness_range=(0.2, 2.0),
            exclude_border=False, n_brightest=1)
        # switched exclude_border to False since we use the edgemask now
        sources = daofind((objfind_image - median_objfind)*np.sqrt(ivar_objfind), mask=totmask)
        if sources is None:
            display.show_image((objfind_image*np.logical_not(totmask)*np.sqrt(ivar_objfind)),
                            chname='S/N objfind_image', cuts=(-2.0, 5.0))
            raise PypeItError(
                "No sources found in the image. Try lowering the significance threshold, "
                f"nsigma = {nsigma:.1f} or adjust the DAOStarFinder parameters."
            )
        if verbose:
            log.info('DAOStarFinder brightest source properties')
            for col in sources.colnames:
                if col not in ('id', 'npix'):
                    sources[col].info.format = '%.2f'  # for consistent table output
            sources.pprint(max_width=76)

        _init_obj_position = sources['x_centroid'][0], sources['y_centroid'][0]
    else:
        _init_obj_position = init_obj_position
        
    initial_guess = (1, _init_obj_position[0], _init_obj_position[1], fwhm*fwhm2sigma, fwhm*fwhm2sigma, 0, 0)
    bounds = ([0,      _init_obj_position[0]-fwhm/3.0, _init_obj_position[1]-fwhm/3.0, fwhm/6.0, fwhm/6.0, -np.pi, -np.inf],
              [np.inf, _init_obj_position[0]+fwhm/3.0, _init_obj_position[1]+fwhm/3.0, fwhm    , fwhm    , np.pi , np.inf])

    # Perform the fit
    # TODO :: May want to generate the image on a finer pixel scale first
    # TODO JFH: The 2D Gaussian fitting should be using the noise and the gpm. This should be
    # implemented with scipy.optimize and a loss function instead of curve_fit
    # Setup the coordinates
    x = np.linspace(0, image.shape[1] - 1, image.shape[1])
    y = np.linspace(0, image.shape[0] - 1, image.shape[0])
    xx, yy = np.meshgrid(x, y, indexing='xy')
    popt, pcov = opt.curve_fit(gaussian2D, (xx, yy), image.ravel() / wlscl,
                               bounds=bounds, p0=initial_guess)
    _, xobj, yobj, sigma_x_gauss, sigma_y_gauss, theta_gauss, _ = popt
    log.info("Gaussian fit gives:")
    log.info("--------------------------------")
    if pixelscale is not None: 
        log.info(f"FWHM_x: {sigma_x_gauss*pixelscale/fwhm2sigma:.2f} arcsec")
        log.info(f"FWHM_y: {sigma_y_gauss*pixelscale/fwhm2sigma:.2f} arcsec")
    else: 
        log.info(f"FWHM_x: {sigma_x_gauss/fwhm2sigma:.2f} pixels")
        log.info(f"FWHM_y: {sigma_y_gauss/fwhm2sigma:.2f} pixels")
    log.info(f"Theta: {np.degrees(theta_gauss):.2f} degrees")
    log.info("--------------------------------")    
    
    # Generate a best fit model
    model = gaussian2D((xx, yy), *popt).reshape(image.shape) * wlscl

    # Optimally extract the object flux using the Gaussian fit
    # Create an apodization window about the object position
    radius = np.sqrt((xx - xobj)**2 + (yy - yobj)**2)
    apodization_window = (radius <= 5.0*sigma)
    # Compute the Gaussian model without an offset to use as the optimal extraction profile
    popt_no_offset = popt.copy()
    popt_no_offset[-1] = 0.0
    gauss_profile = np.clip(gaussian2D((xx, yy), *popt_no_offset).reshape(image.shape) * wlscl, 0.0, None)
    optkern = gauss_profile*apodization_window
    optkern /= np.sum(optkern)
    # Optimally extract at the object position
    mean_image, median_image, std_image = sigma_clipped_stats(image[np.logical_not(totmask)], sigma=3.0)
    image_skysub = image - median_image
    ivar_denom = np.sum(_gpm*optkern)
    var_denom = utils.inverse(ivar_denom)
    ivar_num = np.sum(_gpm*_ivar*optkern**2)
    ivar_opt = ivar_num * var_denom
    flux_opt = np.sum(_gpm*_ivar*image_skysub*optkern) * utils.inverse(ivar_num)
    tot_weight = np.sum(_gpm*_ivar*optkern)
    sigma_opt = np.sqrt(utils.inverse(ivar_opt))
    # Print out a report for the S/N of the optimally extracted object
    log.info(f"Optimal extraction of the brightest object gives")
    log.info(f"     -----------------------------")
    log.info(f"     | (x, y)  = {xobj:6.2f}, {yobj:6.2f}  |")
    log.info(f"     |   Flux  = {flux_opt:7.3f}         |")
    log.info(f"     |   Sigma = {sigma_opt:7.3f}         |")
    log.info(f"     |   S/N   = {flux_opt/sigma_opt:6.2f}          |")
    log.info(f"     -----------------------------")
 
    return popt, pcov, model, _init_obj_position, flux_opt, sigma_opt


def dar_fitfunc(radec, coord_ra, coord_dec, datfit, wave, obstime, location, pressure,
                temperature, rel_humidity):
    """
    Generates a fitting function to calculate the offset due to differential
    atmospheric refraction

    Args:
        radec (tuple):
            A tuple containing two floats representing the shift in ra and dec
            due to DAR.
        coord_ra (float):
            RA in degrees
        coord_dec (float):
            Dec in degrees
        datfit (`numpy.ndarray`_):
            The RA and DEC that the model needs to match
        wave (float):
            Wavelength to calculate the DAR
        location (`astropy.coordinates.EarthLocation`_):
            observatory location
        pressure (float):
            Outside pressure at `location`
        temperature (float):
            Outside ambient air temperature at `location`
        rel_humidity (float):
            Outside relative humidity at `location`. This should be between 0 to 1.

    Returns:
        float: chi-squared difference between datfit and model
    """
    (diff_ra, diff_dec) = radec
    # Generate the coordinate with atmospheric conditions
    coord_atmo = SkyCoord(coord_ra + diff_ra, coord_dec + diff_dec, unit=(units.deg, units.deg))
    coord_altaz = coord_atmo.transform_to(AltAz(obstime=obstime, location=location, obswl=wave,
                                          pressure=pressure, temperature=temperature,
                                          relative_humidity=rel_humidity))
    # Return chi-squared value
    return np.sum((np.array([coord_altaz.alt.value, coord_altaz.az.value])-datfit)**2)


def correct_grating_shift(wave_eval, wave_curr, spl_curr, wave_ref, spl_ref, order=2):
    """
    Using spline representations of the blaze profile, calculate the grating
    correction that should be applied to the current spectrum (suffix ``curr``)
    relative to the reference spectrum (suffix ``ref``). The grating correction
    is then evaluated at the wavelength array given by ``wave_eval``.

    Args:
        wave_eval (`numpy.ndarray`_):
            Wavelength array to evaluate the grating correction
        wave_curr (`numpy.ndarray`_):
            Wavelength array used to construct spl_curr
        spl_curr (`scipy.interpolate.interp1d`_):
            Spline representation of the current blaze function (based on the illumflat).
        wave_ref (`numpy.ndarray`_):
            Wavelength array used to construct spl_ref
        spl_ref (`scipy.interpolate.interp1d`_):
            Spline representation of the reference blaze function (based on the illumflat).
        order (int):
            Polynomial order used to fit the grating correction.

    Returns:
        `numpy.ndarray`_: The grating correction to apply
    """
    log.info("Calculating the grating correction")
    # Calculate the grating correction
    grat_corr_tmp = spl_curr(wave_eval) / spl_ref(wave_eval)
    # Determine the useful overlapping wavelength range
    minw, maxw = max(np.min(wave_curr), np.min(wave_ref)), max(np.min(wave_curr), np.max(wave_ref))
    # Perform a low-order polynomial fit to the grating correction (should be close to linear)
    wave_corr = (wave_eval - minw) / (maxw - minw)  # Scale wavelengths to be of order 0-1
    wblz = np.where((wave_corr > 0.1) & (wave_corr < 0.9))  # Remove the pixels that are within 10% of the edges
    coeff_gratcorr = np.polyfit(wave_corr[wblz], grat_corr_tmp[wblz], order)
    grat_corr = np.polyval(coeff_gratcorr, wave_corr)
    # Return the estimates grating correction
    return grat_corr


def extract_point_source(
    wave, flxcube, ivarcube, bpmcube, wcscube, exptime, min_frac_use=0.05, whitelight_range=None,
    fluxed=False, subpixel=20, boxcar_radius=None, fwhm=1.5, skysub_resid=True, snr_thresh=5.0,
    manual_position=None, opt_prof_method='fit_gauss', spectrograph='keck_kcrm', show_qa=False
):
    """
    Extract a spectrum of a standard star from a datacube

    Parameters
    ----------
    wave : `numpy.ndarray`_
        Wavelength array for the datacube
    flxcube : `numpy.ndarray`_
        Datacube of the flux
    ivarcube : `numpy.ndarray`_
        Datacube of the inverse variance
    bpmcube : `numpy.ndarray`_
        Datacube of the bad pixel mask
    wcscube : `astropy.wcs.WCS`_
        WCS of the datacube
    exptime : float
        Exposure time listed in the header of the datacube
    min_frac_use : float, optional
        Minimum accepted value for the sum of the normalized object profile
        across the spatial direction.  For each spectral pixel, if the majority
        of the object profile has been masked, i.e., the sum of the normalized
        object profile across the spatial direction is less than `min_frac_use`,
        the optimal extraction will also be masked. The default value is 0.05.
    whitelight_range (None, list, optional):
        A two element list that specifies the minimum and maximum
        wavelengths (in Angstroms) to use when constructing the white light
        image (format is: [min_wave, max_wave]). If None, the cube will be
        collapsed over the full wavelength range. If a list is provided an
        either element of the list is None, then the minimum/maximum
        wavelength range of that element will be set by the minimum/maximum
        wavelength of all_wave.
    fluxed : bool, optional
        Is the datacube fluxed?
    subpixel : int, optional
        Number of pixels to subpixelate spectrum when creating mask
    boxcar_radius : float, optional
        Radius of the circular boxcar (in arcseconds) to use for the extraction.
        If None, the radius is set to 4 times the sigma of the 2D Gaussian
        fit to the whitelight image.
    fwhm : float, optional
        FWHM of the PSF in arcseconds. Use to determine the degree of smoothing
        of the whitelight image, the kernel size for the initial object finding,
        and the bounds of the parameters for the 2D Gaussian fit.  Note that if
        the opt_prof_method is set to 'user_gauss', this parameter will be also
        be used as the FWHM of the the 2D (symmetric) Gaussian spatial profile
        for optimal extraction.  Default is 1.5 arcseconds.
    no_sksyub : bool, optional
        If True, the residual sky will not be subtracted from the datacube or
        the whitelight image. Default is False.
    snr_thresh : float, optional
        The signal-to-noise ratio threshold to use when determining the initial
        object position in the whitelight image with DAOStarFinder (this is the
        nsigma parameter in fitGaussian2D). Default is 5.0 
    manual_position : tuple, optional
        Manual position of the object in user-facing cube spatial coordinates,
        i.e. (x, y). Default is None, which means that the position will be
        determined from the whitelight image. This follows the Ginga/DS9 image
        viewer convention: if the image has shape (ny, nx), a position (x, y)
        refers to image[y, x].
    opt_prof_method : str, optional

        The method to be used to determine the object spatial profile for
        optimal extraction.  Options are ``'fit_gauss'``, ``'user_gauss'``, or
        ``'whitelight'``. The default is ``'fit_gauss'``.  Behavior is as
        follows:

            - ``'fit_gauss'``:  Use the 2D Gaussian (possibly asymmetric)
              Gaussian fit to the whitelight image which was used to determine
              the object position.  This creates a model using
              func:`pypeit.core.datacube.fitGaussian2D` but the offset is set to
              zero.

            - ``'user_gauss'``: Use a 2D symmetric Gaussian profile. The FWHM of
              the Gaussian is determined by the fwhm parameter, which was also
              used for the object finding.  Note that this assumes the spatial
              pixel (spaxel) sampling is the same in both x and y.

            - ``'whitelight'``: Use the whitelight image to determine a
              non-parametric spatial profile. The whitelight image is smoothed
              with a Gaussian kernel of width 0.5*sigma, where sigma is the
              standard deviation (fwhm/2.35) corresponding to the fwhm
              parameter. 

    spectrograph : str, pypeit.spectrographs.spectrograph.Spectrograph, optional
        The spectrograph used to take the data.
    show_qa : bool, optional
        If True, the function will display alignment QA  images in ginga.

    Returns
    -------
    sobjs : :class:`~pypeit.specobjs.SpecObjs`
        SpecObjs object containing the extracted spectrum
    spec2dobj : :class:`~pypeit.spec2dobj.Spec2DObj`
        Spec2DObj object containing a psuedo 2D spectrum for visualization purposes
    wl_img : :class:`numpy.ndarray`
        A whitelight image of the input cube (of type :class:`numpy.ndarray`)
        which is the average flux over the set of pixels in the wavelength range
        specified by wavemin and wavemax that are not masked by the badpixel
        mask cube or the sigma clipping mask.
    wl_ivar : :class:`numpy.ndarray`
        The inverse variance of the whitelight image.
    wl_gpm : :class:`numpy.ndarray`
        A good pixel mask for the whitelight image. A value of True indicates a
        good pixel.
    """
    if whitelight_range is None:
        whitelight_range = [np.min(wave), np.max(wave)]

    fwhm2sigma = 1.0 / (2 * np.sqrt(2 * np.log(2)))
    # Load the spectrograph
    _spectrograph = load_spectrograph(spectrograph) if isinstance(spectrograph, str) else spectrograph

    # Generate a spec1d object to hold the extracted spectrum
    log.info("Initialising a PypeIt SpecObj spec1d file")
    sobj = specobj.SpecObj(_spectrograph.pypeline, "DET01", SLITID=0)
    sobj.RA = wcscube.wcs.crval[0]
    sobj.DEC = wcscube.wcs.crval[1]
    sobj.SLITID = 0

    # Convert from counts/s/Ang/arcsec**2 to counts. The sensitivity function expects counts as input
    numwave, numyy, numxx = flxcube.shape
    dspat_x = np.abs(wcscube.wcs.cdelt[0] * wcscube.wcs.cunit[0].to(units.arcsec))
    dspat_y = np.abs(wcscube.wcs.cdelt[1] * wcscube.wcs.cunit[1].to(units.arcsec))
    arcsecSQ = dspat_x * dspat_y
    dspat = np.sqrt(arcsecSQ)
    if fluxed:
        # The datacube is flux calibrated, in units of 10^-17 erg/s/cm**2/Ang/arcsec**2
        # Scale the flux and ivar cubes to be in units of erg/s/cm**2/Ang
        unitscale = arcsecSQ
    else:
        # Scale the flux and ivar cubes to be in units of counts. pypeit_sensfunc expects counts as input
        deltawave = wcscube.wcs.cdelt[2]*wcscube.wcs.cunit[2].to(units.Angstrom)
        unitscale = exptime * deltawave * arcsecSQ


    # Apply the relevant scaling
    _flxcube = flxcube * unitscale
    _ivarcube = ivarcube / unitscale**2
    _gpmcube = np.logical_not(bpmcube)

    # Calculate the variance cube
    _varcube = utils.inverse(_ivarcube)

    # Generate a whitelight image, and fit a 2D Gaussian to estimate centroid and width
    log.info("Making white light image")
    wl_img, wl_ivar, wl_gpm = make_whitelight_fromcube(
        _flxcube, _ivarcube, _gpmcube, wave=wave, wavemin=whitelight_range[0],
        wavemax=whitelight_range[1]
    )
    popt, pcov, model, init_obj_position, flux_opt, sigma_opt = fitGaussian2D(
        wl_img, ivar=wl_ivar, gpm=wl_gpm, init_obj_position=manual_position,
        fwhm=fwhm/dspat, nsigma=snr_thresh, norm=False, pixelscale=dspat
    )
    _, xpos_gauss, ypos_gauss, sigma_x_gauss, sigma_y_gauss, theta_gauss, _ = popt
    gaussian_position = xpos_gauss, ypos_gauss
    
    # Object location for extraction 
    if manual_position is not None:
        xobj, yobj = manual_position
    else: 
        xobj, yobj = gaussian_position
    
    # Setup the coordinates of the mask
    y = np.linspace(0, numyy - 1, numyy * subpixel)
    x = np.linspace(0, numxx - 1, numxx * subpixel)
    yy, xx = np.meshgrid(y, x, indexing='ij')
    
    # Set the radius of the extraction boxcar for the sky determination
    if boxcar_radius is None:
        nsig = 4  # 4 sigma should be far enough... Note: percentage enclosed for 2D Gaussian = 1-np.exp(-0.5 * nsig**2)
        wid = nsig * max(popt[3], popt[4])
    else:
        # Set the user-defined radius
        wid = boxcar_radius / np.sqrt(arcsecSQ)
    # Set the radius of the extraction boxcar for the sky determination
    log.info(f"Using a boxcar radius of {wid*dspat:0.2f} arcsec")
    widsky = 2 * wid
    
    # Generate a mask
    log.info("Generating an object mask")
    newshape = (numyy * subpixel, numxx * subpixel)
    mask = np.zeros(newshape)
    ww = np.where((np.sqrt((xx - xobj) ** 2 + (yy - yobj) ** 2) < wid))
    mask[ww] = 1
    mask = utils.rebinND(mask, (numyy, numxx)).reshape(1, numyy, numxx)

    # Generate a sky mask
    log.info("Generating a sky mask")
    newshape = (numyy * subpixel, numxx * subpixel)
    smask = np.zeros(newshape)
    ww = np.where((np.sqrt((xx - xobj) ** 2 + (yy - yobj) ** 2) < widsky))
    smask[ww] = 1
    smask = utils.rebinND(smask, (numyy, numxx)).reshape(1, numyy, numxx)
    # Subtract off the object mask region, so that we just have an annulus around the object
    smask -= mask

    if skysub_resid:
        log.info("Subtracting the residual sky")
        # Subtract the residual sky from the datacube
        skymask = np.logical_not(bpmcube) * smask
        skycube = _flxcube * skymask
        skyspec = skycube.sum(axis=(1,2))
        nrmsky = skymask.sum(axis=(1,2))
        skyspec *= utils.inverse(nrmsky)
        _flxcube -= skyspec.reshape((numwave, 1, 1))
        # Now subtract the residual sky from the white light image
        sky_val = np.sum(wl_img[None, :, :] * smask) / np.sum(smask)
        wl_img -= sky_val
    else: 
        log.info("The residual sky will not be subtracted")
        skyspec = np.zeros(numwave)

    log.info("Extracting a boxcar spectrum of datacube")
    # Construct an image that contains the fraction of flux included in the
    # boxcar extraction at each wavelength interval
    norm_flux = wl_img[None,:,:] * mask
    norm_flux /= np.sum(norm_flux)
    # Extract boxcar
    cntmask = np.logical_not(bpmcube) * mask  # Good pixels within the masked region around the standard star
    flxscl = (norm_flux * cntmask).sum(axis=(1,2))  # This accounts for the flux that is missing due to masked pixels
    scimask = _flxcube * cntmask
    varmask = _varcube * cntmask**2
    nrmcnt = utils.inverse(flxscl)
    box_flux = scimask.sum(axis=(1,2)) * nrmcnt
    box_var = varmask.sum(axis=(1,2)) * nrmcnt**2
    box_gpm = flxscl > 1/3  # Good pixels are those where at least one-third of the standard star flux is measured

    # Store the BOXCAR extraction information
    sobj.BOX_R_PIX = wid  # Size of boxcar radius in pixels
    sobj.BOX_WAVE = wave.astype(float)
    if fluxed:
        sobj.BOX_FLAM = box_flux
        sobj.BOX_FLAM_SIG = np.sqrt(box_var)
        sobj.BOX_FLAM_IVAR = utils.inverse(box_var)
    else:
        sobj.BOX_COUNTS = box_flux
        sobj.BOX_COUNTS_SIG = np.sqrt(box_var)
        sobj.BOX_COUNTS_IVAR = utils.inverse(box_var)
        sobj.BOX_COUNTS_SKY = skyspec  # This is not the real sky, it is the residual sky. The datacube is presumed to be sky subtracted
    sobj.BOX_MASK = box_gpm
    sobj.S2N = np.median(box_flux * np.sqrt(utils.inverse(box_var)))

    # Now do the OPTIMAL extraction
    log.info("Extracting an optimal spectrum of datacube")
    # First, we need to rearrange the datacube and inverse variance cube into a 2D array.
    # The 3D -> 2D conversion is done so that there is a spectral and spatial dimension,
    # and the brightest white light pixel is transformed to be at the centre column of the 2D
    # array. Then, the second brightest white light pixel is transformed to be next to the centre
    # column of the 2D array, and so on. This is done so that the optimal extraction algorithm
    # can be applied.

    # Setup the coordinates
    x = np.linspace(0, wl_img.shape[1] - 1, wl_img.shape[1])
    y = np.linspace(0, wl_img.shape[0] - 1, wl_img.shape[0])
    xx, yy = np.meshgrid(x, y, indexing='xy')

    if opt_prof_method == 'user_gauss':
        log.info("Optimal extraction with user_gauss method:")
        log.info("------------------------------------------------")
        log.info(f"User provided FWHM: {fwhm:.2f} arcsec")
        log.info("------------------------------------------------")
        # Generate a Gaussian kernel
        fwhm_pix = fwhm/dspat
        intflux = 1
        sigma_x, sigma_y = fwhm_pix*fwhm2sigma, fwhm_pix*fwhm2sigma
        theta, offset, = 0.0, 0.0
        optkern = gaussian2D(
            (xx, yy), intflux, xobj, yobj, sigma_x, sigma_y, theta, offset).reshape(wl_img.shape)
        # Normalise the kernel
        optkern /= np.sum(optkern)
    elif opt_prof_method == 'fit_gauss':
        # Compute the Gaussian model without an offset to use as the optimal extraction profile
        popt_no_offset = popt.copy()
        popt_no_offset[-1] = 0.0
        gauss_profile = np.clip(gaussian2D((xx, yy), *popt_no_offset).reshape(wl_img.shape), 0.0, None)
        optkern = gauss_profile/np.sum(gauss_profile)
        # Print out the properties of the Gaussian
        log.info("Optimal extraction with fit_gauss method:")
        log.info("--------------------------------")
        log.info(f"FWHM_x: {sigma_x_gauss*dspat/fwhm2sigma:.2f} arcsec")
        log.info(f"FWHM_y: {sigma_y_gauss*dspat/fwhm2sigma:.2f} arcsec")
        log.info("--------------------------------")        
    elif opt_prof_method == 'whitelight':
        log.info("Optimal extraction with whitelight method: using whitelight image as a non-parametric spatial profile")
        sigma = fwhm/dspat*fwhm2sigma
        smoothed_wl_img = ndimage.gaussian_filter(wl_img, sigma=0.5*sigma, mode='constant', cval=0.0)
        # Create an apodization window using the coordinates and the specified center
        radius = np.sqrt((xx - xobj)**2 + (yy - yobj)**2)
        apodization_window = np.exp(-radius**2 / (2 * (5*sigma)**2))
        # Apply the apodization window to the smoothed image
        apo_smooth_wl_img = smoothed_wl_img * apodization_window
        optkern = apo_smooth_wl_img/np.sum(apo_smooth_wl_img)


    optkern_masked = optkern * mask[0,:,:]
    # Normalise the white light image
    optkern_masked /= np.sum(optkern_masked)
    asrt = np.argsort(optkern_masked, axis=None)
    # Need to ensure that the number of pixels in the object profile is even
    if asrt.size % 2 != 0:
        # Remove the pixel with the lowest kernel weight.
        # It should be a zero value (given the mask), so it doesn't matter if we remove it
        asrt = asrt[1:]
    # Now sort the indices of the pixels in the object profile
    tmp = asrt.reshape((asrt.size//2, 2))
    objprof_idx = np.append(tmp[:,0], tmp[::-1,1])
    objprof = optkern_masked[np.unravel_index(objprof_idx, optkern.shape)]

    # Now slice the datacube and inverse variance cube into a 2D array
    spec, spat = np.meshgrid(np.arange(numwave), objprof_idx, indexing='ij')
    spatspl = np.apply_along_axis(np.unravel_index, 1, spat, optkern.shape)
    # Now slice the datacube and corresponding cubes/vectors into a series of 2D arrays
    numspat = objprof_idx.size
    flxslice = (spec, spatspl[:,0,:], spatspl[:,1,:])
    flxcube2d = _flxcube[flxslice]
    ivarcube2d = _ivarcube[flxslice]
    gpmcube2d = np.logical_not(bpmcube[flxslice])
    waveimg = wave.reshape((numwave,1)).repeat(numspat, axis=1)
    skyimg = np.zeros((numwave, numspat))  # Note, the residual sky has already been subtracted off _flxcube
    oprof = objprof.reshape((1,numspat)).repeat(numwave, axis=0)
    thismask = np.ones_like(flxcube2d, dtype=bool)

    # Now do the optimal extraction
    sobj.TRACE_SPAT = np.full(numwave, numspat/2.0)
    sobj.extract_optimal(
        flxcube2d, ivarcube2d, gpmcube2d, waveimg, skyimg, thismask, oprof,
        min_frac_use=min_frac_use
    )

    # TODO :: The optimal extraction may suffer from residual DAR correction issues. This is because the
    #      :: object profile assumes that the white light image represents the true spatial profile of the
    #      :: object. One possibility is to fit a (linear?) model to the ratio of box/optimal extraction
    #      :: and then apply this model to the optimal extraction. This is a bit of a fudge.
    # Note that extract.extract_optimal() stores the optimal extraction in the
    # sobj.OPT_COUNTS, sobj.OPT_COUNTS_SIG, and sobj.OPT_COUNTS_IVAR attributes.
    # We need to store the fluxed extraction into the FLAM attributes (a slight fudge).
    if fluxed:
        sobj.OPT_FLAM = sobj.OPT_COUNTS
        sobj.OPT_FLAM_SIG = sobj.OPT_COUNTS_SIG
        sobj.OPT_FLAM_IVAR = sobj.OPT_COUNTS_IVAR

    # Make a specobjs object
    sobjs = specobjs.SpecObjs()
    sobjs.add_sobj(sobj)

    slit_left = np.full((numwave,1), 0.0)
    slit_righ = np.full((numwave,1), float(numspat))

    # TODO, fix hardwired PYP_SPEC
    det_container = _spectrograph.get_detector_par(1)
    slits = slittrace.SlitTraceSet(
        slit_left, slit_righ, _spectrograph.pypeline, detname=det_container.name, nspat=numspat,
        PYP_SPEC=_spectrograph.name, specmin=np.zeros(1), specmax=np.full(1, float(numwave)),
        maskdef_id=None, maskdef_objpos=None, maskdef_offset=None, maskdef_slitcen=None,
        maskdef_designtab=None
    )
    tilts = (waveimg - waveimg.min())/(waveimg.max() - waveimg.min())


    # Set the bit for pixels which were masked by the extraction.
    # For extractmask, True = Good, False = Bad
    bitmask = ImageBitMaskArray(flxcube2d.shape)
    bitmask.turn_on('BPM', select=np.logical_not(gpmcube2d))

    # Make a pseudo spec2d object with these outputs.
    spec2d = spec2dobj.Spec2DObj(
        sciimg=flxcube2d, ivarraw=ivarcube2d, skymodel=skyimg, bkg_redux_skymodel=None,
        objmodel=skyimg, ivarmodel=ivarcube2d, scaleimg=np.array([1.0], dtype=float),
        bpmmask=bitmask, detector=det_container, slits=slits, wavesol=None, waveimg=waveimg,
        tilts=tilts, sci_spat_flexure=None, sci_spec_flexure=None, vel_corr=None, vel_type=None,
        maskdef_designtab=None
    )

    if show_qa: 
        #  Show object finding QA
        whitelight_objfind_qa(wl_img, wl_ivar, wl_gpm, model, gaussian_position, init_obj_position,
                          manual_position=manual_position)
        # Show the extraction QA
        extract_chname = 'opt_prof_method:' + opt_prof_method
        viewer, ch_model = display.show_image(optkern_masked, chname=extract_chname, wcs_match=True,
                                                cuts=(0.0, np.max(optkern_masked)))

    # Return the specobjs object and the spec2d object
    return sobjs, spec2d, wl_img, wl_ivar, wl_gpm


def whitelight_objfind_qa(wl_img, wl_ivar, wl_gpm, gaussian_model, gaussian_position, init_obj_position, 
                          manual_position=None, channel_prefix=''):
    """
    Generate ginga QA for the white light image point source object finding. 

    Image coordinates are quoted as (x, y), matching Ginga and DS9 readouts.
    In numpy terms, if the image has shape (ny, nx), a position (x, y) refers
    to image[y, x].
    
    Parameters
    ----------
    wl_img : `numpy.ndarray`_
        The white light image
    wl_ivar : `numpy.ndarray`_
        The inverse variance of the white light image
    wl_gpm : `numpy.ndarray`_
        The good pixel mask of the white light image
    gaussian_model : `numpy.ndarray`_
        The 2D Gaussian model of the object from datacube.
    gaussian_position : tuple
        The object position in image (x, y) coordinates determined from the
        Gaussian fit to the object.
    init_obj_position : tuple
        The initial object position in image (x, y) coordinates determined
        from DAOStarFinder.
    manual_position : tuple, optional
        The manual extraction object position in image (x, y) coordinates.
    channel_prefix : str, optional
        The prefix to use for the channel name in ginga. Default is ''.
    """

    ny, nx = wl_img.shape
    mean, med, sigma = sigma_clipped_stats(wl_img[wl_gpm], sigma_lower=5.0, sigma_upper=5.0)
    cut_min = mean - 1.0 * sigma
    cut_max = mean + 5.0 * sigma
    viewer, ch_wl = display.show_image(
        wl_img, chname=channel_prefix + 'Whitelight', wcs_match=True, cuts=(cut_min, cut_max))
    mean_snr, med_snr, sigma_snr = sigma_clipped_stats((wl_img*np.sqrt(wl_ivar))[wl_gpm], 
                                                       sigma_lower=5.0, sigma_upper=5.0)
    cut_min_snr = mean_snr - 1.0 * sigma_snr
    cut_max_snr = mean_snr + 5.0 * sigma_snr
    viewer, ch_snr = display.show_image(
        wl_img*np.sqrt(wl_ivar), chname=channel_prefix + 'Whitelight S/N', wcs_match=True,
        cuts=(cut_min_snr, cut_max_snr))
    viewer, ch_model = display.show_image(
        gaussian_model, chname=channel_prefix + 'Gaussian Model',
        wcs_match=True, cuts=(cut_min, cut_max))

    # TODO Add WCS
    ch_list = [ch_wl, ch_model, ch_snr]
    for ich, ch in enumerate(ch_list):
        display.show_points(viewer, ch, [gaussian_position[1]], [gaussian_position[0]],
                            color='red',
                            legend='Gaussian           ; x={:.2f}, y={:.2f}'.format(gaussian_position[0],
                                                                                     gaussian_position[1]),
                            legend_spec=0.05*ny, legend_spat=0.5*nx)
        display.show_points(viewer, ch, [init_obj_position[1]], [init_obj_position[0]],
                            color='green',
                            legend='DAOStarFinder ; x={:.2f}, y={:.2f}'.format(init_obj_position[0],
                                                                               init_obj_position[1]),
                            legend_spec=0.10*ny, legend_spat=0.5*nx)
        if manual_position is not None:
            display.show_points(viewer, ch, [manual_position[1]], [manual_position[0]],
                            color='orange',
                            legend='Manual              ; x={:.2f}, y={:.2f}'.format(manual_position[0],
                                                                                     manual_position[1]),
                            legend_spec=0.15*ny, legend_spat=0.5*nx)
    


def make_good_skymask(slitimg, tilts):
    """
    Mask the spectral edges of each slit (i.e. the pixels near the ends of the
    detector in the spectral direction). Some extreme values of the tilts are
    only sampled with a small fraction of the pixels of the slit width. This
    leads to a bad extrapolation/determination of the sky model.

    Args:
        slitimg (`numpy.ndarray`_):
            An image of the slit indicating which slit each pixel belongs to
        tilts (`numpy.ndarray`_):
            Spectral tilts.

    Returns:
        `numpy.ndarray`_: A mask of the good sky pixels (True = good)
    """
    log.info("Masking edge pixels where the sky model is poor")
    # Initialise the GPM
    gpm = np.zeros(slitimg.shape, dtype=bool)
    # Find unique slits
    unq = np.unique(slitimg[slitimg>0])
    for uu in range(unq.size):
        # Find the x,y pixels in this slit
        ww = np.where((slitimg == unq[uu]) & (tilts != 0.0))
        # Mask the bottom pixels first
        wb = np.where(ww[0] == 0)[0]
        wt = np.where(ww[0] == np.max(ww[0]))[0]
        # Calculate the maximum tilt from the bottom row, and the miminum tilt from the top row
        maxtlt = np.max(tilts[0,  ww[1][wb]])
        mintlt = np.min(tilts[-1, ww[1][wt]])
        # Mask all values below this maximum
        gpm[ww] = (tilts[ww] >= maxtlt) & (tilts[ww] <= mintlt)  # The signs are correct here.
    return gpm


def get_output_filename(output_dir, fil, par_outfile, combine, native=False, idx=1):
    """
    Get the output filename of a datacube, given the input

    Parameters
    ----------
    output_dir : str
        The output directory to save the datacube.
    fil : str
        The spec2d filename.
    par_outfile : str
        The user-specified output filename (see cubepar['output_filename'])
    combine : bool
        Should the input frames be combined into a single datacube?
    native : bool, optional
        If the cube is saved at the native sampling of the spectrograph,
        add 'native' to the filename to indicate this.
    idx : int, optional
        Index of filename to be saved. Required if combine=False.

    Returns
    -------
    str
        The output filename to use.
    """
    # Perform a check on par_outfile
    par_outfile = '' if par_outfile is None else par_outfile
    # Prepare the native text
    native_txt = '_native' if native else ''
    if combine:
        if par_outfile == '':
            par_outfile = 'datacube.fits'
        # Check if we needs to append an extension
        return par_outfile if '.fits' in par_outfile else f'{par_outfile}.fits'
    if par_outfile == '':
        return fil.replace('spec2d_', f'spec3d{native_txt}_')
    # Finally, if nothing else, use the output filename as a prefix, and a numerical suffix
    return os.path.join(output_dir, os.path.splitext(par_outfile)[0] + f'{native_txt}_{idx:03}.fits')


def get_output_whitelight_filename(output_dir, outfile):
    """
    Given the output filename of a datacube, create an appropriate whitelight
    fits file name

    Parameters
    ----------
    output_dir : str
        The output directory to save the datacube.
    outfile : str
        The output filename used for the datacube.

    Returns
    -------
    A string containing the output filename to use for the whitelight image.
    """
    return os.path.join(output_dir, os.path.splitext(outfile)[0] + "_whitelight.fits")


def get_whitelight_pixels(all_wave, all_slitid, min_wl, max_wl):
    """
    Determine which pixels are included within the specified wavelength range

    Args:
        all_wave (`numpy.ndarray`_, list):
            List of `numpy.ndarray`_ wavelength images. The length of the list is the number of spec2d frames.
            Each element of the list contains a wavelength image that provides the wavelength at each pixel on
            the detector, with shape is (nspec, nspat).
        all_slitid (`numpy.ndarray`_, list):
            List of `numpy.ndarray`_ slitid images. The length of the list is the number of spec2d frames.
            Each element of the list contains a slitid image that provides the slit number at each pixel on
            the detector, with shape (nspec, nspat).
        min_wl (float):
            Minimum wavelength to consider
        max_wl (float):
            Maximum wavelength to consider

    Returns:
        :obj:`tuple`: The first element of the tuple is a list of `numpy.ndarray`_ slitid images
        (or a single `numpy.ndarray`_ slitid image if only one spec2d frame is provided),
        shape is (nspec, nspat), where a zero value corresponds to an excluded pixel
        (either outside the desired wavelength range, a bad pixel, a pixel not on the slit).
        All other pixels have a value equal to the slit number. The second element of the tuple
        is the wavelength difference between the maximum and minimum wavelength in the desired
        wavelength range.
    """
    # Check if lists or numpy arrays are input
    list_inputs = [all_wave, all_slitid]
    if all([isinstance(l, list) for l in list_inputs]):
        numframes = len(all_wave)
        if not all([len(l) == numframes for l in list_inputs]):
            raise PypeItError("All input lists must have the same length")
        # Store in the following variables
        _all_wave, _all_slitid = all_wave, all_slitid
    elif all([not isinstance(l, list) for l in list_inputs]):
        _all_wave, _all_slitid = [all_wave], [all_slitid]
        numframes = 1
    else:
        raise PypeItError("The input lists must either all be lists (of the same length) or all be numpy arrays")
    if max_wl < min_wl:
        raise PypeItError("The maximum wavelength must be greater than the minimum wavelength")
    # Initialise the output
    out_slitid = [np.zeros(_all_slitid[0].shape, dtype=int) for _ in range(numframes)]
    # Loop over all frames and find the pixels that are within the wavelength range
    if min_wl < max_wl:
        # Loop over files and determine which pixels are within the wavelength range
        for ff in range(numframes):
            ww = np.where((_all_wave[ff] > min_wl) & (_all_wave[ff] < max_wl))
            out_slitid[ff][ww] = _all_slitid[ff][ww]
    else:
        log.warning("Datacubes do not completely overlap in wavelength.")
        out_slitid = _all_slitid
        min_wl, max_wl = None, None
        for ff in range(numframes):
            this_wave = _all_wave[ff][_all_slitid[ff] > 0]
            tmp_min = np.min(this_wave)
            tmp_max = np.max(this_wave)
            if min_wl is None or tmp_min < min_wl:
                min_wl = tmp_min
            if max_wl is None or tmp_max > max_wl:
                max_wl = tmp_max
    # Determine the wavelength range
    wavediff = max_wl - min_wl
    # Need to return a single slitid image if only one frame, otherwise return a list of slitid images.
    # Also return the wavelength difference
    return out_slitid[0] if numframes == 1 else out_slitid, wavediff


def get_whitelight_range(wavemin, wavemax, wl_range):
    """
    Get the wavelength range to use for the white light images

    Parameters
    ----------
    wavemin : float
        Automatically determined minimum wavelength to use for making the white
        light image.
    wavemax : float
        Automatically determined maximum wavelength to use for making the white
        light image.
    wl_range : list
        Two element list containing the user-specified values to manually
        override the automated values determined by PypeIt.

    Returns
    -------
    wlrng : list
        A two element list containing the minimum and maximum wavelength to use
        for the white light images
    """
    wlrng = [wavemin, wavemax]
    if wl_range[0] is not None:
        if wl_range[0] < wavemin:
            log.warning(
                f"The user-specified minimum wavelength ({wl_range[0]:.2f}) to use for the white "
                f"light\nimages is lower than the recommended value ({wavemin:.2f}),\n"
                "which ensures that all spaxels cover the same wavelength range."
            )
        wlrng[0] = wl_range[0]
    if wl_range[1] is not None:
        if wl_range[1] > wavemax:
            log.warning(
                f"The user-specified maximum wavelength ({wl_range[1]:.2f}) to use for the white "
                "light\nimages is greater than the recommended value ({wavemax:.2f}),\n"
                "which ensures that all spaxels cover the same wavelength range."
            )
        wlrng[1] = wl_range[1]
    log.info("The white light images will cover the wavelength range: {0:.2f}A - {1:.2f}A".format(wlrng[0], wlrng[1]))
    return wlrng

def make_whitelight(
    output_wcs, flxcube, ivarcube, gpmcube, wave, output_dir, outfile, whitelight_range=None,
    overwrite=False
):
    """
    Generate a white light image using an input cube and write to a file.

    Parameters
    ----------
    output_wcs : :class:`astropy.wcs.WCS`
        Output world coordinate system.
    flxcube : :class:`numpy.ndarray`
        3D datacube (the final element contains the wavelength dimension).
    ivarcube : :class:`numpy.ndarray`
        3D inverse variance cube (the final element contains the wavelength
        dimension).
    gpmcube : :class:`numpy.ndarray`, bool
        3D good pixel mask cube (the final element contains the wavelength
        dimension).  A value of True indicates a good pixel.
    wave : :class:`numpy.ndarray`
        A 1D array containing the wavelength at each spectral coordinate of the
        datacube. The shape of the wavelength array is (nwave,).
    output_dir : str
        The output directory to save the datacube. 
    outfile : str
        The output filename for the datacube.        
    whitelight_range : list, optional
        A two element list that specifies the minimum and maximum wavelengths
        (in Angstroms) to use when constructing the white light image (format
        is: [min_wave, max_wave]). If None, the cube will be collapsed over the
        full wavelength range. If a list is provided and either element of the
        list is None, then the minimum/maximum wavelength range of that element
        will be set by the minimum/maximum wavelength of all_wave.
    overwrite : bool, optional
        Overwrite an existing files
    """

    whitelight_wcs = output_wcs.celestial
    # Check if the user requested a white light image
    if whitelight_range is not None:
        # Grab the WCS of the white light image
        # Determine the wavelength range of the whitelight image
        _whitelight_range = (wave[0] if whitelight_range[0] is None else whitelight_range[0],
                             wave[-1] if whitelight_range[1] is None else whitelight_range[1])
    else:
        _whitelight_range = (wave[0], wave[-1])

    log.info(
        f"White light image covers the wavelength range {_whitelight_range[0]:.2f} A - "
        f"{_whitelight_range[1]:.2f} A"
    )
    # Get the output filename for the white light image
    out_whitelight = get_output_whitelight_filename(output_dir, outfile)
    whitelight, ivar_whitelight, gpm_whitelight = make_whitelight_fromcube(
        flxcube, ivarcube, gpmcube, wave=wave, wavemin=_whitelight_range[0],
        wavemax=_whitelight_range[1]
    )
    log.info(f"Saving white light image as: {out_whitelight}")
    primary_hdu = fits.PrimaryHDU(whitelight, header=whitelight_wcs.to_header())
    # TODO: Primary HDUs should *NOT* have an EXTNAME keyword.  You can set the
    # name of the primary extension, I think (`primary_hdu.name =
    # 'WHITELIGHT'`), but EXTNAME is for extensions.
#    primary_hdu.header['EXTNAME'] = 'WHITELIGHT'
    ivar_hdu = fits.ImageHDU(ivar_whitelight, name='IVAR')
    gpm_hdu = fits.ImageHDU(gpm_whitelight.astype(np.uint8), name='GPM')

    hdul = fits.HDUList([primary_hdu, ivar_hdu, gpm_hdu])
    hdul.writeto(out_whitelight, overwrite=overwrite)


def make_whitelight_fromcube(cube, ivarcube, gpmcube, sigclip=5.0,
                             wave=None, wavemin=None, wavemax=None):
    """
    Generate a white light image using an input cube.

    Parameters
    ----------
    cube : :class:`numpy.ndarray`
        3D datacube (the first element contains the wavelength dimension)
    ivarcube : :class:`numpy.ndarray`
        3D inverse variance cube (the first element contains the wavelength dimension).
    gpmcube : :class:`numpy.ndarray`, bool
        3D good-pixel mask cube (the first element contains the wavelength
        dimension).  A value of True indicates a good pixel.
    sigclip : float, optional
        Flag outliers using sigma-clipping based on this sigma value (both above
        and blow the median).  If None, do not perform sigma clipping.
    wave : :class:`numpy.ndarray`, optional
        1D wavelength array. Only required if wavemin or wavemax are not None.
    wavemin : float, optional
        Minimum wavelength (same units as wave) to be included in the whitelight
        image.  You must provide wave as well if you want to reduce the
        wavelength range.
    wavemax : float, optional
        Maximum wavelength (same units as wave) to be included in the whitelight
        image.  You must provide wave as well if you want to reduce the
        wavelength range.

    Returns
    -------
    whitelight : :class:`numpy.ndarray`
        A whitelight image of the input cube, which is the average flux over the
        set of pixels in the wavelength range specified by wavemin and wavemax
        that are not masked by the badpixel mask cube or the sigma clipping
        mask.
    ivar_whitelight : :class:`numpy.ndarray`
        The inverse variance of the whitelight image.
    gpm_whitelight : :class:`numpy.ndarray`
        A good pixel mask for the whitelight image. A value of True indicates a
        good pixel.
    """
    # Make a wavelength cut, if requested
    if wavemin is not None or wavemax is not None:
        # Make some checks on the input
        if wave is None:
            raise PypeItError("wave variable must be supplied to create white light image with wavelength cuts")
        else:
            if wave.size != cube.shape[0]:
                raise PypeItError("wave variable should have the same length as the third axis of cube.")
        # assign wavemin & wavemax if one is not provided
        if wavemin is None:
            wavemin = np.min(wave)
        if wavemax is None:
            wavemax = np.max(wave)
        ww = np.where((wave >= wavemin) & (wave <= wavemax))[0]
        wmin, wmax = ww[0], ww[-1]+1
        cutcube = cube[wmin:wmax, :, :]
        cutivar = ivarcube[wmin:wmax, :, :]
        # Cut the bad pixel mask and convert it to a good pixel mask
        cutgpmcube = gpmcube[wmin:wmax, :, :]
    else:
        cutcube = cube.copy()
        cutivar = ivarcube.copy()
        cutgpmcube = gpmcube.copy()

    # Apply find_min_max_out
    data = np.ma.MaskedArray(cutcube, mask=np.logical_not(cutgpmcube))
    if sigclip is not None:
        sc = SigmaClip(sigma=sigclip, maxiters=25, cenfunc='median', stdfunc=utils.nan_mad_std)
        data_clipped, lower, upper = sc(data, axis=0, masked=True, return_bounds=True)
    gpm_sigclip = np.logical_not(data_clipped.mask)

    # Compute the average flux over the set of pixels that are not masked by gpm_sigclip
    npix_whitelight = np.sum(gpm_sigclip, axis=0)
    inv_npix_whitelight = utils.inverse(npix_whitelight)
    whitelight_sum = np.sum((cutcube*gpm_sigclip), axis=0)
    gpm_whitelight = npix_whitelight > 0
    whitelight = whitelight_sum*gpm_whitelight * inv_npix_whitelight

    # Compute the formal corresponding variance over the set of pixels that are not masked by
    # gpm_sigclip
    cut_var = utils.inverse(cutivar)
    var_sum_whitelight = np.sum((cut_var*gpm_sigclip), axis=0)
    var_whitelight = var_sum_whitelight * inv_npix_whitelight**2
    ivar_whitelight = utils.inverse(var_whitelight)*gpm_whitelight

    return whitelight, ivar_whitelight, gpm_whitelight


def load_imageWCS(filename, ext=0):
    """
    Load an image and return the image and the associated WCS.

    Args:
        filename (str):
            A fits filename of an image to be used when generating white light
            images. Note, the fits file must have a valid 3D WCS.
        ext (bool, optional):
            The extension that contains the image and WCS

    Returns:
        :obj:`tuple`: An `numpy.ndarray`_ with the 2D image data and a
        `astropy.wcs.WCS`_ with the image WCS.
    """
    imghdu = fits.open(filename)
    image = imghdu[ext].data
    imgwcs = wcs.WCS(imghdu[ext].header)
    # Return required info
    return image, imgwcs


def align_user_offsets(ifu_ra, ifu_dec, ra_offset, dec_offset):
    """
    Align the RA and DEC of all input frames, and then
    manually shift the cubes based on user-provided offsets.
    The offsets should be specified in arcseconds, and the
    ra_offset should include the cos(dec) factor.

    Args:
        ifu_ra (`numpy.ndarray`_):
            A list of RA values of the IFU (one value per frame)
        ifu_dec (`numpy.ndarray`_):
            A list of Dec values of the IFU (one value per frame)
        ra_offset (`numpy.ndarray`_):
            A list of RA offsets to be applied to the input pixel values (one value per frame).
            Note, the ra_offset MUST contain the cos(dec) factor. This is the number of degrees
            on the sky that represents the telescope offset.
        dec_offset (`numpy.ndarray`_):
            A list of Dec offsets to be applied to the input pixel values (one value per frame).
            This is the number of degrees on the sky that represents the telescope offset.

    Returns:
        A tuple containing a new set of RA and Dec offsets for each frame.
        Both arrays are of type `numpy.ndarray`_, and are in units of degrees.
    """
    # First, translate all coordinates to the coordinates of the first frame
    # Note: You do not need cos(dec) here, this just overrides the IFU coordinate centre of each frame
    #       The cos(dec) factor should be input by the user, and should be included in the self.opts['ra_offset']
    ref_shift_ra = ifu_ra[0] - ifu_ra
    ref_shift_dec = ifu_dec[0] - ifu_dec
    numfiles = len(ra_offset)
    out_ra_offsets = [0.0 for _ in range(numfiles)]
    out_dec_offsets = [0.0 for _ in range(numfiles)]
    for ff in range(numfiles):
        # Apply the shift
        out_ra_offsets[ff] = ref_shift_ra[ff] - ra_offset[ff]
        out_dec_offsets[ff] = ref_shift_dec[ff] - dec_offset[ff]
        log.info(
            f"Spatial shift of cube #{ff + 1}:\nRA, DEC (arcsec) = {ra_offset[ff]*3600.0:+0.3f} "
            f"E, {dec_offset[ff]*3600.0:+0.3f} N"
        )
    return out_ra_offsets, out_dec_offsets


def set_voxel_sampling(spatscale, specscale, dspat=None, dwv=None):
    """
    This function checks if the spatial and spectral scales of all frames are consistent.
    If the user has not specified either the spatial or spectral scales, they will be set here.

    Parameters
    ----------
    spatscale : `numpy.ndarray`_
        2D array, shape is (N, 2), listing the native spatial scales of N spec2d frames.
        spatscale[:,0] refers to the spatial pixel scale of each frame
        spatscale[:,1] refers to the slicer scale of each frame
        Each element of the array must be in degrees
    specscale : `numpy.ndarray`_
        1D array listing the native spectral scales of multiple frames. The length of this array should be equal
        to the number of frames you are using. Each element of the array must be in Angstrom
    dspat: :obj:`float`, optional
        Spatial scale to use as the voxel spatial sampling. If None, a new value will be derived based on the inputs
    dwv: :obj:`float`, optional
        Spectral scale to use as the voxel spectral sampling. If None, a new value will be derived based on the inputs

    Returns
    -------
    _dspat : :obj:`float`
        Spatial sampling
    _dwv : :obj:`float`
        Wavelength sampling
    """
    # Make sure all frames have consistent pixel scales
    ratio = (spatscale[:, 0] - spatscale[0, 0]) / spatscale[0, 0]
    if np.any(np.abs(ratio) > 1E-4):
        log.warning("The pixel scales of all input frames are not the same!")
        spatstr = ", ".join(["{0:.6f}".format(ss) for ss in spatscale[:,0]*3600.0])
        log.info("Pixel scales of all input frames:\n" + spatstr + "arcseconds")
    # Make sure all frames have consistent slicer scales
    ratio = (spatscale[:, 1] - spatscale[0, 1]) / spatscale[0, 1]
    if np.any(np.abs(ratio) > 1E-4):
        log.warning("The slicer scales of all input frames are not the same!")
        spatstr = ", ".join(["{0:.6f}".format(ss) for ss in spatscale[:,1]*3600.0])
        log.info("Slicer scales of all input frames:\n" + spatstr + "arcseconds")
    # Make sure all frames have consistent wavelength sampling
    ratio = (specscale - specscale[0]) / specscale[0]
    if np.any(np.abs(ratio) > 1E-2):
        log.warning("The wavelength samplings of the input frames are not the same!")
        specstr = ", ".join(["{0:.6f}".format(ss) for ss in specscale])
        log.info("Wavelength samplings of all input frames:\n" + specstr + "Angstrom")

    # If the user has not specified the spatial scale, then set it appropriately now to the largest spatial scale
    _dspat = np.max(spatscale) if dspat is None else dspat
    log.info("Adopting a square pixel spatial scale of {0:f} arcsec".format(3600.0 * _dspat))
    # If the user has not specified the spectral sampling, then set it now to the largest value
    _dwv = np.max(specscale) if dwv is None else dwv
    log.info("Adopting a wavelength sampling of {0:f} Angstrom".format(_dwv))
    return _dspat, _dwv


def check_inputs(list_inputs):
    """
    This function checks the inputs to several of the cube building routines, and makes sure they are all consistent.
    Often, this is to make check if all inputs are lists of the same length, or if all inputs are 2D `numpy.ndarray`_.
    The goal of the routine is to return a consistent set of lists of the input.

    Parameters
    ----------
    list_inputs : :obj:`list`
        A list of inputs to check.

    Returns
    -------
    list_inputs : :obj:`list`
        A list of inputs that have been checked for consistency.
    """
    if all([isinstance(l, list) for l in list_inputs]):
        # Several frames are being combined. Check the lists have the same length
        numframes = len(list_inputs[0])
        if not all([len(l) == numframes for l in list_inputs]):
            raise PypeItError("All input lists must have the same length")
        # The inputs are good, return as is
        return tuple(list_inputs)
    elif all([not isinstance(l, list) for l in list_inputs]):
        # Just a single frame - store as single element lists
        ret_list = ()
        for l in list_inputs:
            ret_list += ([l],)
        return ret_list
    else:
        raise PypeItError("The input arguments should all be of type 'list', or all not be of type 'list':")


def wcs_bounds(raImg, decImg, waveImg, slitid_img_gpm, ra_offsets=None, dec_offsets=None,
               ra_min=None, ra_max=None, dec_min=None, dec_max=None, wave_min=None, wave_max=None):
    """
    Calculate the bounds of the WCS and the expected edges of the voxels, based
    on user-specified parameters or the extremities of the data. 

    Parameters
    ----------
    raImg : (`numpy.ndarray`_, list):
        A list of 2D array containing the RA of each pixel, with shape (nspec, nspat)
    decImg : (`numpy.ndarray`_, list):
        A list of 2D array containing the Dec of each pixel, with shape (nspec, nspat)
    waveImg (`numpy.ndarray`_, list):
        A list of 2D array containing the wavelength of each pixel, with shape (nspec, nspat)
    slitid_img_gpm : (`numpy.ndarray`_, list):
        A list of 2D array containing the spat ID of each pixel, with shape (nspec, nspat).
        A value of 0 indicates that the pixel is not on a slit. All other values indicate the
        slit spatial ID.
    ra_offsets : list, optional
        A list of the RA offsets for each frame
    dec_offsets : list, optional
        A list of the Dec offsets for each frame
    ra_min : :obj:`float`, optional
        Minimum RA of the WCS
    ra_max : :obj:`float`, optional
        Maximum RA of the WCS
    dec_min : :obj:`float`, optional
        Minimum Dec of the WCS
    dec_max : :obj:`float`, optional
        Maximum Dec of the WCS
    wave_min : :obj:`float`, optional
        Minimum wavelength of the WCS
    wave_max : :obj:`float`, optional
        Maximum wavelength of the WCS

    Returns
    -------
    _ra_min : :obj:`float`
        Minimum RA of the WCS
    _ra_max : :obj:`float`
        Maximum RA of the WCS
    _dec_min : :obj:`float`
        Minimum Dec of the WCS
    _dec_max : :obj:`float`
        Maximum Dec of the WCS
    _wave_min : :obj:`float`
        Minimum wavelength of the WCS
    _wave_max : :obj:`float`
        Maximum wavelength of the WCS
    """
    # Check if the ra_offsets and dec_offsets are specified
    if ra_offsets is None or dec_offsets is None:
        if isinstance(raImg, list):
            ra_offsets = [0.0]*len(raImg)
            dec_offsets = [0.0]*len(raImg)
        else:
            ra_offsets = 0.0
            dec_offsets = 0.0
    # Check the inputs
    _raImg, _decImg, _waveImg, _slitid_img_gpm, _ra_offsets, _dec_offsets = \
        check_inputs([raImg, decImg, waveImg, slitid_img_gpm, ra_offsets, dec_offsets])
    numframes = len(_raImg)

    # Loop over the frames and get the bounds - start by setting the default values
    _ra_min, _ra_max = ra_min, ra_max
    _dec_min, _dec_max = dec_min, dec_max
    _wave_min, _wave_max = wave_min, wave_max
    for fr in range(numframes):
        # Only do calculations if the min/max inputs are not specified
        # Get the RA, Dec, and wavelength of the pixels on the slit
        if ra_min is None or ra_max is None:
            this_ra = _raImg[fr][_slitid_img_gpm[fr] > 0]
            tmp_min = np.min(this_ra) - _ra_offsets[fr]
            tmp_max = np.max(this_ra) - _ra_offsets[fr]
            if fr == 0 or tmp_min < _ra_min:
                _ra_min = tmp_min
            if fr == 0 or tmp_max > _ra_max:
                _ra_max = tmp_max
        if dec_min is None or dec_max is None:
            this_dec = _decImg[fr][_slitid_img_gpm[fr] > 0]
            tmp_min = np.min(this_dec) - _dec_offsets[fr]
            tmp_max = np.max(this_dec) - _dec_offsets[fr]
            if fr == 0 or tmp_min < _dec_min:
                _dec_min = tmp_min
            if fr == 0 or tmp_max > _dec_max:
                _dec_max = tmp_max
        if wave_min is None or wave_max is None:
            this_wave = _waveImg[fr][_slitid_img_gpm[fr] > 0]
            tmp_min = np.min(this_wave)
            tmp_max = np.max(this_wave)
            if fr == 0 or tmp_min < _wave_min:
                _wave_min = tmp_min
            if fr == 0 or tmp_max > _wave_max:
                _wave_max = tmp_max
    # Return the bounds
    return _ra_min, _ra_max, _dec_min, _dec_max, _wave_min, _wave_max


def create_wcs(raImg, decImg, waveImg, slitid_img_gpm, dspat, dwave,
               ra_offsets=None, dec_offsets=None,
               ra_min=None, ra_max=None, dec_min=None, dec_max=None, wave_min=None, wave_max=None,
               reference=None, collapse=False, equinox=2000.0, specname="PYP_SPEC"):
    """
    Create a WCS and the expected edges of the voxels, based on user-specified
    parameters or the extremities of the data.

    Parameters
    ----------
    raImg : (`numpy.ndarray`_, list):
        A list of 2D array containing the RA of each pixel, with shape (nspec, nspat)
    decImg : (`numpy.ndarray`_, list):
        A list of 2D array containing the Dec of each pixel, with shape (nspec, nspat)
    waveImg (`numpy.ndarray`_, list):
        A list of 2D array containing the wavelength of each pixel, with shape (nspec, nspat)
    slitid_img_gpm : (`numpy.ndarray`_, list):
        A list of 2D array containing the spat ID of each pixel, with shape (nspec, nspat).
        A value of 0 indicates that the pixel is not on a slit. All other values indicate the
        slit spatial ID.
    dspat : float
        Spatial size of each square voxel (in arcsec).
    dwave : float
        Linear wavelength step of each voxel (in Angstroms)
    ra_offsets : list, optional
        List of RA offsets for each frame (degrees)
    dec_offsets : list, optional
        List of Dec offsets for each frame (degrees)
    ra_min : float, optional
        Minimum RA of the WCS (degrees)
    ra_max : float, optional
        Maximum RA of the WCS (degrees)
    dec_min : float, optional
        Minimum Dec of the WCS (degrees)
    dec_max : float, optional
        Maximum Dec of the WCS (degrees)
    wave_min : float, optional
        Minimum wavelength of the WCS (degrees)
    wave_max : float, optional
        Maximum wavelength of the WCS (degrees)
    reference : str, optional
        Filename of a fits file that contains a WCS in the Primary HDU.
    collapse : bool, optional
        If True, the spectral dimension will be collapsed to a single channel
        (primarily for white light images)
    equinox : float, optional
        Equinox of the WCS
    specname : str, optional
        Name of the spectrograph

    Returns
    -------
    cubewcs : `astropy.wcs.WCS`_
        astropy WCS to be used for the combined cube
    voxedges : tuple
        A three element tuple (z, y, x) containing the bin edges in the x, y (RA, Dec spatial) and
        z (wavelength) dimensions
    reference_image : `numpy.ndarray`_
        The reference image to be used for the cross-correlation. Can be None.
    """
    # Setup the cube ranges
    _ra_min, _ra_max, _dec_min, _dec_max, _wave_min, _wave_max = \
        wcs_bounds(raImg, decImg, waveImg, slitid_img_gpm,
                   ra_offsets=ra_offsets, dec_offsets=dec_offsets,
                   ra_min=ra_min, ra_max=ra_max, dec_min=dec_min, dec_max=dec_max, wave_min=wave_min, wave_max=wave_max)

    # Grab cos(dec) for convenience. Use the average of the min and max dec
    cosdec = np.cos(0.5*(_dec_min+_dec_max) * np.pi / 180.0)

    # Number of voxels in each dimension
    numra = int((_ra_max - _ra_min) * cosdec / dspat)
    numdec = int((_dec_max - _dec_min) / dspat)
    numwav = int(np.round((_wave_max - _wave_min) / dwave))

    # If a white light WCS is being generated, make sure there's only 1 wavelength bin
    if collapse:
        dwave = _wave_max - _wave_min
        numwav = 1

    # Generate a master WCS to register all frames
    coord_min = [_ra_min, _dec_min, _wave_min]
    coord_dlt = [-dspat, dspat, dwave]

    # If a reference image is being used and a white light image is requested (collapse=True) update the celestial parts
    reference_image = None
    if reference is not None:
        # Load the requested reference image
        reference_image, imgwcs = load_imageWCS(reference)
        # Update the celestial WCS
        coord_min[:2] = imgwcs.wcs.crval
        coord_dlt[:2] = imgwcs.wcs.cdelt
        numdec, numra = reference_image.shape

    cubewcs = generate_WCS(coord_min, coord_dlt, numra, equinox=equinox, name=specname)
    log.info(
        f'\n{"-"*40}'
        "\nParameters of the WCS:"
        f"\nRA   min = {coord_min[0]}"
        f"\nDEC  min = {coord_min[1]}"
        f"\nWAVE min, max = {_wave_min}, {_wave_max}"
        f"\nSpaxel size = {3600.0 * dspat} arcsec"
        f"\nWavelength step = {dwave} A"
        f'\n{"-"*40}'
    )

    # Generate the output binning
    xbins = np.arange(1 + numra) - 0.5
    ybins = np.arange(1 + numdec) - 0.5
    spec_bins = np.arange(1 + numwav) - 0.5
    voxedges = (spec_bins, ybins, xbins)
    return cubewcs, voxedges, reference_image


def generate_WCS(crval, cdelt, numra, equinox=2000.0, name="PYP_SPEC"):
    """
    Generate a WCS that will cover all input spec2D files

    Args:
        crval (list):
            3 element list containing the [RA, DEC, WAVELENGTH] of
            the reference pixel
        cdelt (list):
            3 element list containing the delta values of the [RA,
            DEC, WAVELENGTH]
        numra (int):
            Number of RA values in the WCS. This is used to ensure
            that the convention of the WCS is so that North is up
            and East is to the left.
        equinox (float, optional):
            Equinox of the WCS

    Returns:
        `astropy.wcs.WCS`_ : astropy WCS to be used for the combined cube
    """
    # Create a new WCS object.
    log.info("Generating WCS")
    w = wcs.WCS(naxis=3)
    w.wcs.equinox = equinox
    w.wcs.name = name
    w.wcs.radesys = 'FK5'
    # Insert the coordinate frame
    w.wcs.cname = ['RA', 'DEC', 'Wavelength']
    w.wcs.cunit = [units.degree, units.degree, units.Angstrom]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN", "WAVE"]
    w.wcs.crval = crval  # RA, DEC, and wavelength zeropoints
    w.wcs.crpix = [numra, 0, 0]  # RA, DEC, and wavelength reference pixels
    #w.wcs.cd = np.array([[cdval[0], 0.0, 0.0], [0.0, cdval[1], 0.0], [0.0, 0.0, cdval[2]]])
    w.wcs.cdelt = cdelt
    w.wcs.lonpole = 180.0  # Native longitude of the Celestial pole
    w.wcs.latpole = 0.0  # Native latitude of the Celestial pole
    return w


def compute_weights_frompix(raImg, decImg, waveImg, sciImg, ivarImg, slitidImg, dspat, dwv, mnmx_wv, wghtsImg,
                            all_wcs, all_tilts, all_slits, all_align, all_dar, ra_offsets, dec_offsets,
                            ra_min=None, ra_max=None, dec_min=None, dec_max=None, wave_min=None, wave_max=None,
                            sn_smooth_npix=None, weight_method='auto',
                            reference_image=None, whitelight_range=None,
                            correct_dar=True, specname="PYPSPEC", init_obj_position=None, show_qa=False):
    r"""
    Calculate wavelength dependent optimal weights. The weighting is currently
    based on a relative :math:`(S/N)^2` at each wavelength. Note, this function
    first prepares a whitelight image, and then calls compute_weights() to
    determine the appropriate weights of each pixel.

    Parameters
    ----------

    raImg : `numpy.ndarray`_, list
        A list of 2D array containing the RA of each pixel, with shape (nspec, nspat)
    decImg : `numpy.ndarray`_, list
        A list of 2D array containing the Dec of each pixel, with shape (nspec, nspat)
    waveImg : `numpy.ndarray`_, list
        A list of 2D array containing the wavelength of each pixel, with shape (nspec, nspat)
    sciImg : `numpy.ndarray`_, list
        A list of 2D array containing the science image of each pixel, with shape (nspec, nspat)
    ivarImg : `numpy.ndarray`_, list
        A list of 2D array containing the inverse variance image of each pixel, with shape (nspec, nspat)
    slitidImg : `numpy.ndarray`_, list
        A list of 2D array containing the slit ID of each pixel, with shape (nspec, nspat)
    dspat : float
        The size of each spaxel on the sky (in degrees)
    dwv : float
        The size of each wavelength pixel (in Angstroms)
    mnmx_wv : `numpy.ndarray`_
        The minimum and maximum wavelengths of every slit and frame. The shape is (Nframes, Nslits, 2),
        The minimum and maximum wavelengths are stored in the [:,:,0] and [:,:,1] indices, respectively.
    wghtsImg : `numpy.ndarray`_, list
        A list of 2D array containing the weights of each pixel, with shape (nspec, nspat)
    all_wcs : `astropy.wcs.WCS`_, list
        A list of WCS objects, one for each frame.
    all_tilts : `numpy.ndarray`_, list
        2D wavelength tilts frame, or a list of tilt frames
    all_slits : :class:`~pypeit.slittrace.SlitTraceSet`, list
        Information stored about the slits, or a list of SlitTraceSet objects
    all_align : :class:`~pypeit.alignframe.AlignmentSplines`, list
        A Class containing the transformation between detector pixel
        coordinates and WCS pixel coordinates, or a list of Alignment
        Splines.
    all_dar : :class:`~pypeit.coadd3d.DARcorrection`, list
        A Class containing the DAR correction information, or a list of DARcorrection
        classes. If a list, it must be the same length as astrom_trans.
    ra_offsets : float, list
        RA offsets for each frame in units of degrees
    dec_offsets : float, list
        Dec offsets for each frame in units of degrees
    ra_min : float, optional
        Minimum RA of the WCS (degrees)
    ra_max : float, optional
        Maximum RA of the WCS (degrees)
    dec_min : float, optional
        Minimum Dec of the WCS (degrees)
    dec_max : float, optional
        Maximum Dec of the WCS (degrees)
    wave_min : float, optional
        Minimum wavelength of the WCS (degrees)
    wave_max : float, optional
        Maximum wavelength of the WCS (degrees)
    sn_smooth_npix : float, optional
        Number of pixels used for determining smoothly varying S/N ratio
        weights.  This is currently not required, since a relative weighting
        scheme with a polynomial fit is used to calculate the S/N weights.
    weight_method : `str`, optional
        Weight method to be used in :func:`~pypeit.coadd.sn_weights`.
        Options are ``'auto'``, ``'constant'``, ``'uniform'``, ``'wave_dependent'``, ``'relative'``, or
        ``'ivar'``. The default is ``'auto'``.  Behavior is as follows:

            - ``'auto'``: Use constant weights if rms_sn < 3.0, otherwise
                use wavelength dependent.

            - ``'constant'``: Constant weights based on rms_sn**2

            - ``'uniform'``: Uniform weighting.

            - ``'wave_dependent'``: Wavelength dependent weights will be
                used irrespective of the rms_sn ratio. This option will not
                work well at low S/N ratio although it is useful for objects
                where only a small fraction of the spectral coverage has high
                S/N ratio (like high-z quasars).

            - ``'relative'``: Calculate weights by fitting to the ratio of
                spectra? Note, relative weighting will only work well when
                there is at least one spectrum with a reasonable S/N, and a
                continuum.  RJC note - This argument may only be better when
                the object being used has a strong continuum + emission lines.
                The reference spectrum is assigned a value of 1 for all
                wavelengths, and the weights of all other spectra will be
                determined relative to the reference spectrum. This is
                particularly useful if you are dealing with highly variable
                spectra (e.g. emission lines) and require a precision better
                than ~1 per cent.

            - ``'ivar'``: Use inverse variance weighting. This is not well
                tested and should probably be deprecated.
    reference_image : `numpy.ndarray`_
        Reference image to use for the determination of the highest S/N spaxel in the image.
    correct_dar : bool, optional
        Correct for the differential atmospheric refraction.  Default is False.
    specname : str
        Name of the spectrograph
    init_obj_position : tuple, optional
        The initial guess for the object position in image (x, y) coordinates,
        matching Ginga or DS9 readouts. In numpy terms, if the image has shape
        (ny, nx), a position (x, y) refers to image[y, x]. If set, this value
        will be input into `fitGaussian2D` as the initial guess for the object
        position. The 2D Gaussian fit will then be performed with the position
        constrained to be within plus or minus fwhm/3 in x and y. If not set,
        the position will be determined by running DAOStarFinder on the image.
        Default is None.
    show_qa : bool, optional
        If True, show QA plots in ginga. 

    Returns
    -------
    weights : `numpy.ndarray`_ 
        a 1D array the same size as all_sci, containing relative wavelength
        dependent weights of each input pixel.
    """
    # Find the wavelength range where all frames overlap
    min_wl, max_wl = get_whitelight_range(np.max(mnmx_wv[:, :, 0]),  # The max blue wavelength
                                          np.min(mnmx_wv[:, :, 1]),  # The min red wavelength
                                          whitelight_range)  # The user-specified values (if any)
    # Get the good white light pixels
    slitid_img_gpm, wavediff = get_whitelight_pixels(waveImg, slitidImg, min_wl, max_wl)

    # Generate the WCS
    image_wcs, voxedge, reference_image = \
        create_wcs(raImg, decImg, waveImg, slitid_img_gpm, dspat, wavediff,
                   ra_offsets=ra_offsets, dec_offsets=dec_offsets,
                   ra_min=ra_min, ra_max=ra_max, dec_min=dec_min, dec_max=dec_max, wave_min=wave_min, wave_max=wave_max,
                   reference=reference_image, collapse=True, equinox=2000.0, specname=specname)

    # Generate the white light image
    # NOTE: hard-coding subpixel=1 in both directions for speed, and combining into a single image
    wl_full, wl_sig, wl_bpm = generate_image_subpixel(image_wcs, voxedge, sciImg, ivarImg, waveImg, slitid_img_gpm, wghtsImg,
                                      all_wcs, all_tilts, all_slits, all_align, all_dar, ra_offsets, dec_offsets,
                                      spec_subpixel=1, spat_subpixel=1, slice_subpixel=1, combine=True,
                                      correct_dar=correct_dar)

    # Compute the weights
    return compute_weights(raImg, decImg, waveImg, sciImg, ivarImg, slitidImg,
                           all_wcs, all_tilts, all_slits, all_align, all_dar, ra_offsets, dec_offsets,
                           wl_full, wl_sig, wl_bpm, dspat, dwv,
                           ra_min=ra_min, ra_max=ra_max, dec_min=dec_min, dec_max=dec_max, wave_min=wave_min,
                           sn_smooth_npix=sn_smooth_npix, weight_method=weight_method, correct_dar=correct_dar, 
                           init_obj_position=init_obj_position, show_qa=show_qa)

# TODO Refactor this, it should not be done this way, instead we should be computing the weights from the final aligned
# cubes, after sigma clipping is performed. See my notes in coadd3d.run()
def compute_weights(raImg, decImg, waveImg, sciImg, ivarImg, slitidImg,
                    all_wcs, all_tilts, all_slits, all_align, all_dar, ra_offsets, dec_offsets,
                    whitelight_img, whitelight_sigma, whitelight_bpm, dspat, dwv,
                    ra_min=None, ra_max=None, dec_min=None, dec_max=None, wave_min=None, wave_max=None,
                    sn_smooth_npix=None, weight_method='auto', correct_dar=True, fwhm=1.5, init_obj_position=None, show_qa=False):
    r"""
    Calculate wavelength dependent optimal weights. The weighting is currently
    based on a relative :math:`(S/N)^2` at each wavelength

    Parameters
    ----------

    raImg : `numpy.ndarray`_, list
        A list of 2D array containing the RA of each pixel, with shape (nspec, nspat)
    decImg : `numpy.ndarray`_, list
        A list of 2D array containing the Dec of each pixel, with shape (nspec, nspat)
    waveImg : `numpy.ndarray`_, list
        A list of 2D array containing the wavelength of each pixel, with shape (nspec, nspat)
    sciImg : `numpy.ndarray`_, list
        A list of 2D array containing the science image of each pixel, with shape (nspec, nspat)
    ivarImg : `numpy.ndarray`_, list
        A list of 2D array containing the inverse variance image of each pixel, with shape (nspec, nspat)
    slitidImg : `numpy.ndarray`_, list
        A list of 2D array containing the slit ID of each pixel, with shape (nspec, nspat)
    all_wcs : `astropy.wcs.WCS`_, list
        A list of WCS objects, one for each frame.
    all_tilts : `numpy.ndarray`_, list
        2D wavelength tilts frame, or a list of tilt frames
    all_slits : :class:`~pypeit.slittrace.SlitTraceSet`, list
        Information stored about the slits, or a list of SlitTraceSet objects
    all_align : :class:`~pypeit.alignframe.AlignmentSplines`, list
        A Class containing the transformation between detector pixel
        coordinates and WCS pixel coordinates, or a list of Alignment
        Splines.
    all_dar : :class:`~pypeit.coadd3d.DARcorrection`, list
        A Class containing the DAR correction information, or a list of DARcorrection
        classes. If a list, it must be the same length as astrom_trans.
    ra_offsets : float, list
        RA offsets for each frame in units of degrees
    dec_offsets : float, list
        Dec offsets for each frame in units of degrees
    whitelight_img : `numpy.ndarray`_
        A 2D array containing a white light image, that was created with the
        input ``all`` arrays.
    whitelight_sigma : `numpy.ndarray`_
        A 2D array containing the standard deviation of the white light image.
        Only used for the QA plot at present. 
    whitelight_bpm : `numpy.ndarray`_, bool
        A 2D array containing a bad pixel mask for the white light image.        
    dspat : float
        The size of each spaxel on the sky (in degrees)
    dwv : float
        The size of each wavelength pixel (in Angstroms)
    sn_smooth_npix : float, optional
        Number of pixels used for determining smoothly varying S/N ratio
        weights.  This is currently not required, since a relative weighting
        scheme with a polynomial fit is used to calculate the S/N weights.
    correct_dar : bool, optional
        Apply the DAR correction to the input data.  The default is True.
    weight_method : `str`, optional
        Weight method to be used in :func:`~pypeit.coadd.sn_weights`.
        Options are ``'auto'``, ``'constant'``, ``'uniform'``, ``'wave_dependent'``, ``'relative'``, or
        ``'ivar'``. The default is ``'auto'``.  Behavior is as follows:

            - ``'auto'``: Use constant weights if rms_sn < 3.0, otherwise
                use wavelength dependent.

            - ``'constant'``: Constant weights based on rms_sn**2

            - ``'uniform'``: Uniform weighting.

            - ``'wave_dependent'``: Wavelength dependent weights will be
                used irrespective of the rms_sn ratio. This option will not
                work well at low S/N ratio although it is useful for objects
                where only a small fraction of the spectral coverage has high
                S/N ratio (like high-z quasars).

            - ``'relative'``: Calculate weights by fitting to the ratio of
                spectra? Note, relative weighting will only work well when
                there is at least one spectrum with a reasonable S/N, and a
                continuum.  RJC note - This argument may only be better when
                the object being used has a strong continuum + emission lines.
                The reference spectrum is assigned a value of 1 for all
                wavelengths, and the weights of all other spectra will be
                determined relative to the reference spectrum. This is
                particularly useful if you are dealing with highly variable
                spectra (e.g. emission lines) and require a precision better
                than ~1 per cent.

            - ``'ivar'``: Use inverse variance weighting. This is not well
                tested and should probably be deprecated.

        fwhm : float, optional
            FWHM of the PSF in arcseconds. Use to determine the degree of smoothing of the whitelight image, the
            kernel size for the initial object finding, and the bounds of the parameters for the 2D Gaussian fit. 
            Default is 1.5 arcseconds.
        init_obj_position : tuple, optional
            The initial guess for the object position in image (x, y)
            coordinates, matching Ginga or DS9 readouts. In numpy terms, if the
            image has shape (ny, nx), a position (x, y) refers to image[y, x].
            If set, this value will be input into `fitGaussian2D` as the
            initial guess for the object position. The 2D Gaussian fit will
            then be performed with the position constrainted to be within plus
            or minus fwhm/3 in x and y. If not set, the position will be
            determined by running DAOStarFinder on the image. Default is None.
        show_qa : bool, optional
            If True, show the object detection QA plot in ginga. Default is False. 

    Returns
    -------
    all_wghts: list
        Either a 2D `numpy.ndarray`_ or a list of 2D `numpy.ndarray`_ arrays
        containing the optimal weights of each pixel for all frames, with shape
        (nspec, nspat).
    """
    log.info("Calculating the optimal weights of each pixel")
    # Check the inputs for combinations of lists or not, and then determine the number of frames
    _raImg, _decImg, _waveImg, _sciImg, _ivarImg, _slitidImg, \
        _all_wcs, _all_tilts, _all_slits, _all_align, _all_dar, _ra_offsets, _dec_offsets = \
            check_inputs([raImg, decImg, waveImg, sciImg, ivarImg, slitidImg,
                          all_wcs, all_tilts, all_slits, all_align, all_dar, ra_offsets, dec_offsets])
    numframes = len(_sciImg)

    # If there's only one frame, use uniform weighting
    if numframes == 1:
        log.warning("Only one frame provided.  Using uniform weighting.")
        return np.ones_like(sciImg)

    # Check the WCS bounds
    _ra_min, _ra_max, _dec_min, _dec_max, _wave_min, _wave_max = \
        wcs_bounds(_raImg, _decImg, _waveImg, _slitidImg, ra_offsets=_ra_offsets, dec_offsets=_dec_offsets,
                   ra_min=ra_min, ra_max=ra_max, dec_min=dec_min, dec_max=dec_max, wave_min=wave_min, wave_max=wave_max)

    # Find the location of the object with the highest S/N in the combined white light image
    _dspat = dspat*units.deg.to(units.arcsec)
    whitelight_ivar = utils.inverse(np.square(whitelight_sigma))
    popt, pcov, model, init_obj_position, flux_opt, sigma_opt = fitGaussian2D(
        whitelight_img, ivar=whitelight_ivar, gpm=np.logical_not(whitelight_bpm),
        fwhm=fwhm/_dspat, init_obj_position=init_obj_position, norm=False,
        pixelscale=_dspat)
    gaussian_position = popt[1], popt[2]
    if show_qa: 
        whitelight_objfind_qa(whitelight_img, utils.inverse(np.square(whitelight_sigma)), 
                                np.logical_not(whitelight_bpm), model, gaussian_position, 
                                init_obj_position, channel_prefix = f'Weights_')
    
    log.info(
        f"Highest S/N object located at spaxel (x, y) = {gaussian_position[0]:.2f}, "
        f"{gaussian_position[1]:.2f}"
    )

    # Make the bin edges to be at +/- 1 pixels around the maximum (i.e. summing 9 pixels total)
    numwav = int((_wave_max - _wave_min) / dwv)
    xbins = np.array([gaussian_position[0]-1, gaussian_position[0]+2]) - 0.5
    ybins = np.array([gaussian_position[1]-1, gaussian_position[1]+2]) - 0.5
    spec_bins = np.arange(1 + numwav) - 0.5
    bins = (spec_bins, ybins, xbins)

    # Grab cos(dec) for convenience. Use the average of the min and max dec.
    cosdec = np.cos(0.5 * (_dec_min + _dec_max) * np.pi / 180.0)
    # Number of spaxels in the RA direction
    numra = int((_ra_max - _ra_min) * cosdec / dspat)

    # Generate a 2D WCS to register all frames
    coord_min = [_ra_min, _dec_min, _wave_min]
    coord_dlt = [-dspat, dspat, dwv]
    whitelightWCS = generate_WCS(coord_min, coord_dlt, numra)
    wcs_scale = (1.0 * whitelightWCS.spectral.wcs.cunit[0]).to(units.Angstrom).value  # Ensures the WCS is in Angstroms

    # Extract the spectrum of the highest S/N object
    flux_stack = np.zeros((numwav, numframes))
    ivar_stack = np.zeros((numwav, numframes))
    for ff in range(numframes):
        log.info(f"Extracting spectrum of highest S/N detection from frame {ff + 1}/{numframes}")
        flxcube, sigcube, bpmcube, normcube, wave = \
            generate_cube_subpixel(whitelightWCS, bins, _sciImg[ff], _ivarImg[ff], _waveImg[ff],
                                   _slitidImg[ff], np.ones(_sciImg[ff].shape), _all_wcs[ff],
                                   _all_tilts[ff], _all_slits[ff], _all_align[ff], _all_dar[ff],
                                   _ra_offsets[ff], _dec_offsets[ff],
                                   spec_subpixel=1, spat_subpixel=1, slice_subpixel=1,
                                   skip_subpix_weights=True, correct_dar=correct_dar)
        # Store the flux and ivar spectra of the highest S/N object.
        # TODO :: This is the flux per spectral pixel, and not per detector pixel.  Is this correct?
        flux_stack[:, ff] = flxcube[:, 0, 0]
        ivar_stack[:, ff] = utils.inverse(sigcube[:, 0, 0])**2

    # Mask out any pixels that are zero in the flux or ivar stack
    mask_stack = (flux_stack != 0.0) & (ivar_stack != 0.0)
    # Obtain a wavelength of each pixel
    wcs_res = whitelightWCS.wcs_pix2world(np.vstack((np.zeros(numwav), np.zeros(numwav), np.arange(numwav))).T, 0)
    wcs_scale = (1.0 * whitelightWCS.wcs.cunit[2]).to_value(units.Angstrom)  # Ensures the WCS is in Angstroms
    wave_spec = wcs_scale * wcs_res[:, 2]
    # Compute the smoothing scale to use
    if sn_smooth_npix is None:
        sn_smooth_npix = int(np.round(0.1 * wave_spec.size))
    rms_sn, weights = coadd.sn_weights(utils.array_to_explist(flux_stack), utils.array_to_explist(ivar_stack), utils.array_to_explist(mask_stack),
                                       sn_smooth_npix=sn_smooth_npix, weight_method=weight_method, verbose=True)

    # Because we pass back a weights array, we need to interpolate to assign each detector pixel a weight
    all_wghts = numframes*[np.ones(_sciImg[0].shape)]
    for ff in range(numframes):
        ww = (slitidImg[ff] > 0)
        all_wghts[ff][ww] = interp1d(wave_spec, weights[ff], kind='cubic',
                                 bounds_error=False, fill_value="extrapolate")(waveImg[ff][ww])
    log.info("Optimal weighting complete")
    return all_wghts


def generate_image_subpixel(image_wcs, bins, sciImg, ivarImg, waveImg, slitid_img_gpm, wghtImg,
                            all_wcs, tilts, slits, astrom_trans, all_dar, ra_offset, dec_offset,
                            spec_subpixel=5, spat_subpixel=5, slice_subpixel=5, combine=False, correct_dar=True):
    """
    Generate a white light image from the input pixels

    Parameters
    ----------
    image_wcs : `astropy.wcs.WCS`_
        World coordinate system to use for the white light images.
    bins : tuple
        A 3-tuple (x,y,z) containing the histogram bin edges in x,y spatial
        and z wavelength coordinates
    sciImg : `numpy.ndarray`_, list
        A list of 2D science images, or a single 2D image containing the
        science data.
    ivarImg : `numpy.ndarray`_, list
        A list of 2D inverse variance images, or a single 2D image
        containing the inverse variance data.
    waveImg : `numpy.ndarray`_, list
        A list of 2D wavelength images, or a single 2D image containing the
        wavelength data.
    slitid_img_gpm : `numpy.ndarray`_, list
        A list of 2D slit ID images, or a single 2D image containing the
        slit ID data.
    wghtImg : `numpy.ndarray`_, list
        A list of 2D weight images, or a single 2D image containing the
        weight data.
    all_wcs : `astropy.wcs.WCS`_, list
        A list of WCS objects, or a single WCS object containing the WCS
        information of each image.
    tilts : `numpy.ndarray`_, list
        2D wavelength tilts frame, or a list of tilt frames (see all_idx)
    slits : :class:`~pypeit.slittrace.SlitTraceSet`, list
        Information stored about the slits, or a list of SlitTraceSet (see
        all_idx)
    astrom_trans : :class:`~pypeit.alignframe.AlignmentSplines`, list
        A Class containing the transformation between detector pixel
        coordinates and WCS pixel coordinates, or a list of Alignment
        Splines (see all_idx)
    all_dar : :class:`~pypeit.coadd3d.DARcorrection`, list
        A Class containing the DAR correction information, or a list of DARcorrection
        classes. If a list, it must be the same length as astrom_trans.
    ra_offset : :obj:`float`, list
        The RA offset to apply to each image, or a list of RA offsets.
    dec_offset : :obj:`float`, list
        The DEC offset to apply to each image, or a list of DEC offsets.
    spec_subpixel : :obj:`int`, optional
        What is the subpixellation factor in the spectral direction. Higher
        values give more reliable results, but note that the time required
        goes as (``spec_subpixel * spat_subpixel * slice_subpixel``). The
        default value is 5, which divides each detector pixel into 5 subpixels
        in the spectral direction.
    spat_subpixel : :obj:`int`, optional
        What is the subpixellation factor in the spatial direction. Higher
        values give more reliable results, but note that the time required
        goes as (``spec_subpixel * spat_subpixel * slice_subpixel``). The
        default value is 5, which divides each detector pixel into 5 subpixels
        in the spatial direction.
    slice_subpixel : :obj:`int`, optional
        What is the subpixellation factor in the slice direction. Higher
        values give more reliable results, but note that the time required
        goes as (``spec_subpixel * spat_subpixel * slice_subpixel``). The
        default value is 5, which divides each IFU slice into 5 subpixels
        in the slice direction.
    combine : :obj:`bool`, optional
        If True, all of the input frames will be combined into a single
        output. Otherwise, individual images will be generated.
    correct_dar : :obj:`bool`, optional
        If True, the DAR correction will be applied to the input images
        before generating the white light images. If False, the DAR
        correction will not be applied.

    Returns
    -------
    wl_imgs : `numpy.ndarray`_
        The white light images for all frames. If combine=True, this will be a single 2D image. 
        Otherwise, it will be a 3D array with dimensions (numra, numdec, numframes).
    sig_imgs : `numpy.ndarray`_
        The standard deviation images for all frames. If combine=True, this will be a single 2D image.
        Otherwise, it will be a 3D array with dimensions (numra, numdec, numframes).
    bpm_imgs : `numpy.ndarray`_, type=bool
        The bad pixel mask images for all frames. If combine=True, this will be a single 2D image.
        Otherwise, it will be a 3D array with dimensions (numra, numdec, numframes).
    """
    # Perform some checks on the input -- note, more complete checks are performed in subpixellate()
    _sciImg, _ivarImg, _waveImg, _slitid_img_gpm, _wghtImg, _all_wcs, _tilts, _slits, _astrom_trans, _all_dar, _ra_offset, _dec_offset = \
        check_inputs([sciImg, ivarImg, waveImg, slitid_img_gpm, wghtImg, all_wcs, tilts, slits, astrom_trans, all_dar, ra_offset, dec_offset])

    # Generate the white light images
    if combine:
        # Subpixellate
        img, sigimg, bpmimg, _ = subpixellate(image_wcs, bins, _sciImg, _ivarImg, _waveImg, _slitid_img_gpm, _wghtImg,
                                 _all_wcs, _tilts, _slits, _astrom_trans, _all_dar, _ra_offset, _dec_offset,
                                 spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel, slice_subpixel=slice_subpixel,
                                 skip_subpix_weights=True, correct_dar=correct_dar)
        return img[0, :, :], sigimg[0, :, :], bpmimg[0, :, :]
    else:
        # Prepare the array of white light images to be stored
        numframes = len(_sciImg)
        numra = bins[2].size - 1
        numdec = bins[1].size - 1
        all_wl_imgs = np.zeros((numdec, numra, numframes))
        all_sig_imgs = np.zeros((numdec, numra, numframes))
        all_bpm_imgs = np.zeros((numdec, numra, numframes), dtype=bool)
        # Loop through all frames and generate white light images
        for fr in range(numframes):
            log.info(f"Creating image {fr + 1}/{numframes}")
            # Subpixellate
            img, sigimg, bpmimg, _ = subpixellate(image_wcs, bins, _sciImg[fr], _ivarImg[fr], _waveImg[fr], _slitid_img_gpm[fr], _wghtImg[fr],
                                     _all_wcs[fr], _tilts[fr], _slits[fr], _astrom_trans[fr], _all_dar[fr], _ra_offset[fr], _dec_offset[fr],
                                     spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel, slice_subpixel=slice_subpixel,
                                     skip_subpix_weights=True, correct_dar=correct_dar)
            all_wl_imgs[:, :, fr] = img[0, :, :]
            all_sig_imgs[:, :, fr] = sigimg[0, :, :]
            all_bpm_imgs[:, :, fr] = bpmimg[0, :, :]
            
        # Return the constructed white light images
        return all_wl_imgs, all_sig_imgs, all_bpm_imgs


def generate_cube_subpixel(
    output_wcs, bins, sciImg, ivarImg, waveImg, slitid_img_gpm, wghtImg, all_wcs, tilts, slits,
    astrom_trans, all_dar, ra_offset, dec_offset, spec_subpixel=5, spat_subpixel=5,
    slice_subpixel=5, skip_subpix_weights=False, correct_dar=True
):
    """
    Save a datacube using the subpixel algorithm. Refer to the
    :func:`~pypeit.core.datcube.subpixellate` docstring for further details
    about this algorithm.
 
    Parameters
    ----------
    output_wcs : :class:`astropy.wcs.WCS`
        Output world coordinate system.
    bins : tuple
        A 3-tuple (z,y,x) containing the histogram bin edges in x,y spatial and
        z wavelength coordinates
    sciImg : :class:`numpy.ndarray`, list
        A list of 2D array containing the counts of each pixel. If a list, the
        shape of each numpy array is (nspec, nspat).
    ivarImg : :class:`numpy.ndarray`, list
        A list of 2D array containing the inverse variance of each pixel. If a
        list, the shape of each numpy array is (nspec, nspat).
    waveImg : :class:`numpy.ndarray`, list
        A list of 2D array containing the wavelength of each pixel. If a list,
        the shape of each numpy array is (nspec, nspat).
    slitid_img_gpm : :class:`numpy.ndarray`, list
        A list of 2D array containing the slitmask of each pixel. If a list, the
        shape of each numpy array is (nspec, nspat).  A zero value indicates
        that a pixel is either not on a slit or it is a bad pixel.  All other
        values are the slit spatial ID number.
    wghtImg : :class:`numpy.ndarray`, list
        A list of 2D array containing the weights of each pixel to be used in
        the combination. If a list, the shape of each numpy array is (nspec,
        nspat).
    all_wcs : :class:`astropy.wcs.WCS`, list
        A list of :class:`astropy.wcs.WCS` objects, one for each spec2d file.
    tilts : list
        A list of :class:`numpy.ndarray` objects, one for each spec2d file,
        containing the tilts of each pixel. The shape of each numpy array is
        (nspec, nspat).
    slits : :class:`~pypeit.slittrace.SlitTraceSet`, list
        A list of :class:`~pypeit.slittrace.SlitTraceSet` objects, one for each
        spec2d file, containing the properties of the slit for each spec2d file
    astrom_trans : :class:`~pypeit.alignframe.AlignmentSplines`, list
        A Class containing the transformation between detector pixel coordinates
        and WCS pixel coordinates, or a list of Alignment Splines (see all_idx)
    all_dar : :class:`~pypeit.coadd3d.DARcorrection`, list
        A class containing the DAR correction information, or a list of
        :class:`~pypeit.coadd3d.DARcorrection` classes. If a list, it must be
        the same length as ``astrom_trans``.
    ra_offset : float, list
        A float or list of floats containing the RA offset of each spec2d file.
    dec_offset : float, list
        A float or list of floats containing the DEC offset of each spec2d file.
    spec_subpixel : int, optional
        What is the subpixellation factor in the spectral direction? Higher
        values give more reliable results, but note that the time required goes
        as (``spec_subpixel * spat_subpixel``). The default value is 5, which
        divides each detector pixel into 5 subpixels in the spectral direction.
    spat_subpixel : int, optional
        What is the subpixellation factor in the spatial direction? Higher
        values give more reliable results, but note that the time required goes
        as (``spec_subpixel * spat_subpixel``). The default value is 5, which
        divides each detector pixel into 5 subpixels in the spatial direction.
    slice_subpixel : int, optional
        What is the subpixellation factor in the slice direction? Higher values
        give more reliable results, but note that the time required goes as
        (``slice_subpixel``). The default value is 5, which divides each IFU
        slice into 5 subslices in the slice direction.
    skip_subpix_weights : bool, optional
        If True, the computationally expensive step to calculate the
        subpixellation weights will be skipped. If set the True, note that the
        variance cubes returned will not be accurate. However, if you are not
        interested in the variance cubes, this can save a lot of time, and this
        is an example where you might consider setting this variable to True.
        The flux datacube is unaffected by this variable.  The default is False.
    correct_dar : bool, optional
        If True, the DAR correction will be applied to the datacube. If the DAR
        correction is not available, the datacube will not be corrected.

    Returns
    -------
    flxcube : :class:`numpy.ndarray`
        The datacube generated from the subpixellated inputs. The shape of the
        datacube is (nwave, nspat1, nspat2).
    sigcube : :class:`numpy.ndarray`
        The error cube (standard deviation).  Shape matches
        ``flxcube``.
    bpmcube : :class:`numpy.ndarray`
        The bad-pixel mask cube. Shape matches ``flxcube``.
    normcube : :class:`numpy.ndarray`
        A cube indicating the occupation number of a given pixel; i.e., the the
        number of input pixels contributing to each output pixel.  This is
        roughly the number of exposures; however, this will depend on the
        subsampling, relative offsets, and masking applied to each individual
        input cube.
    wave : :class:`numpy.ndarray`
        1D array containing the wavelength at each spectral coordinate of the
        datacube. The shape of the wavelength array is (nwave,).
    """
 
    # Subpixellate
    flxcube, sigcube, bpmcube, normcube = subpixellate(
        output_wcs, bins, sciImg, ivarImg, waveImg, slitid_img_gpm, wghtImg, all_wcs, tilts,
        slits, astrom_trans, all_dar, ra_offset, dec_offset, spec_subpixel=spec_subpixel,
        spat_subpixel=spat_subpixel, slice_subpixel=slice_subpixel,
        skip_subpix_weights=skip_subpix_weights, correct_dar=correct_dar
    )

    # Get wavelength of each pixel
    nspec = flxcube.shape[0]
    wcs_scale = (1.0*output_wcs.spectral.wcs.cunit[0]).to(units.Angstrom).value  # Ensures the WCS is in Angstroms
    wave = wcs_scale * output_wcs.spectral.wcs_pix2world(np.arange(nspec), 0)[0]

    return flxcube, sigcube, bpmcube, normcube, wave

# DEVELOPER NOTES: RJC is working towards making subpixellate a faster routine, and sometimes uses this decorator
# find out the bottlenecks and how to speed things up. Please leave this decorator in for the time-being, uncommented.
# @line_profiler.profile
def subpixellate(
    output_wcs, bins, sciImg, ivarImg, waveImg, slitid_img_gpm, wghtImg, all_wcs, tilts, slits,
    astrom_trans, all_dar, ra_offset, dec_offset, spec_subpixel=5, spat_subpixel=5,
    slice_subpixel=5, skip_subpix_weights=False, correct_dar=True, verbose=False
):
    r"""
    Subpixellate the input data into a datacube. This algorithm splits each
    detector pixel into multiple subpixels and each IFU slice into multiple
    subslices.  Then, the algorithm assigns each subdivided detector pixel to a
    voxel. For example, if ``spec_subpixel = spat_subpixel = slice_subpixel =
    5``, then each detector pixel is divided into :math:`5^3=125` subpixels.
    Alternatively, when spec_subpixel = spat_subpixel = slice_subpixel = 1, this
    corresponds to the nearest grid point (NGP) algorithm.

    .. important::

        If ``spec_subpixel > 1`` or ``spat_subpixel > 1`` or ``slice_subpixel >
        1``, the errors will be correlated, and the covariance is not being
        tracked, so the errors will not be (quite) right. There is a tradeoff
        one has to make between sampling and better looking cubes, versus no
        sampling and better behaved errors.

    Parameters
    ----------
    output_wcs : :class:`astropy.wcs.WCS`
        Output world coordinate system.
    bins : tuple
        A 3-tuple (z,y,x) containing the histogram bin edges in x,y spatial and
        z wavelength coordinates
    sciImg : :class:`numpy.ndarray`, list
        A list of 2D array containing the counts of each pixel. The shape of
        each 2D array is (nspec, nspat).
    ivarImg : :class:`numpy.ndarray`, list
        A list of 2D array containing the inverse variance of each pixel. The
        shape of each 2D array is (nspec, nspat).
    waveImg : :class:`numpy.ndarray`, list
        A list of 2D array containing the wavelength of each pixel. The shape of
        each 2D array is (nspec, nspat).
    slitid_img_gpm : :class:`numpy.ndarray`, list
        A list of 2D array containing the slitmask of each pixel. The shape of
        each 2D array is (nspec, nspat).  A zero value indicates that a pixel is
        either not on a slit or it is a bad pixel.  All other values are the
        slit spatial ID number.
    wghtImg : :class:`numpy.ndarray`, list
        A list of 2D array containing the weights of each pixel to be used in
        the combination. The shape of each 2D array is (nspec, nspat).
    all_wcs : :class:`astropy.wcs.WCS`, list
        A list of :class:`astropy.wcs.WCS` objects, one for each spec2d file.
    tilts : list
        A list of :class:`numpy.ndarray` objects, one for each spec2d file,
        containing the tilts of each pixel. The shape of each 2D array is
        (nspec, nspat).
    slits : :class:`~pypeit.slittrace.SlitTraceSet`, list
        A list of :class:`~pypeit.slittrace.SlitTraceSet` objects, one for each
        spec2d file, containing the properties of the slit for each spec2d file.
    astrom_trans : :class:`~pypeit.alignframe.AlignmentSplines`, list
        A Class containing the transformation between detector pixel coordinates
        and WCS pixel coordinates, or a list of Alignment Splines (see all_idx).
    all_dar : :class:`~pypeit.coadd3d.DARcorrection`, list
        A Class containing the DAR correction information, or a list of
        :class:`~pypeit.coadd3d.DARcorrection` classes. If a list, it must be
        the same length as astrom_trans.
    ra_offset : float, list
        A float or list of floats containing the RA offset of each spec2d file
        relative to the first spec2d file
    dec_offset : float, list
        A float or list of floats containing the DEC offset of each spec2d file
        relative to the first spec2d file
    spec_subpixel : int, optional
        What is the subpixellation factor in the spectral direction? Higher
        values give more reliable results, but note that the time required goes
        as (``spec_subpixel * spat_subpixel``). The default value is 5, which
        divides each detector pixel into 5 subpixels in the spectral direction.
    spat_subpixel : int, optional
        What is the subpixellation factor in the spatial direction? Higher
        values give more reliable results, but note that the time required goes
        as (``spec_subpixel * spat_subpixel``). The default value is 5, which
        divides each detector pixel into 5 subpixels in the spatial direction.
    slice_subpixel : int, optional
        What is the subpixellation factor in the slice direction? Higher values
        give more reliable results, but note that the time required goes as
        (``slice_subpixel``). The default value is 5, which divides each IFU
        slice into 5 subslices in the slice direction.
    skip_subpix_weights : bool, optional
        If True, the computationally expensive step to calculate the
        subpixellation weights will be skipped. If set the True, note that the
        variance cubes returned will not be accurate. However, if you are not
        interested in the variance cubes, this can save a lot of time, and this
        is an example where you might consider setting this variable to True.
        The flux datacube is unaffected by this variable.  The default is False.
    correct_dar : bool, optional
        If True, the DAR correction will be applied to the datacube. The default
        is True.
    verbose : bool, optional
        If True, the code will print out more information. The default is False.

    Returns
    -------
    flxcube : :class:`numpy.ndarray`
        The datacube generated from the subpixellated inputs. The shape of the
        datacube is (nwave, nspat1, nspat2).
    sigcube : :class:`numpy.ndarray`
        The error cube (standard deviation).  Shape matches
        ``flxcube``.
    bpmcube : :class:`numpy.ndarray`
        The bad-pixel mask cube. Shape matches ``flxcube``.
    normcube : :class:`numpy.ndarray`
        A cube indicating the occupation number of a given pixel; i.e., the the
        number of input pixels contributing to each output pixel.  This is
        roughly the number of exposures; however, this will depend on the
        subsampling, relative offsets, and masking applied to each individual
        input cube.
    """
    # Check the inputs for combinations of lists or not
    _sciImg, _ivarImg, _waveImg, _gpmImg, _wghtImg, _all_wcs, _tilts, _slits, _astrom_trans, _all_dar, _ra_offset, _dec_offset = \
        check_inputs([sciImg, ivarImg, waveImg, slitid_img_gpm, wghtImg, all_wcs, tilts, slits, astrom_trans, all_dar, ra_offset, dec_offset])
    numframes = len(_sciImg)

    # Prepare the output arrays
    outshape = (bins[0].size-1, bins[1].size-1, bins[2].size-1)
    binrng = np.array([[bins[0][0], bins[0][-1]], [bins[1][0], bins[1][-1]], [bins[2][0], bins[2][-1]]])
    voxscale = outshape / (binrng[:, 1] - binrng[:, 0])  # shape (3,)
    voxoffset = binrng[:, 0] * voxscale  # shape (3,)
    flxcube, varcube, normcube = np.zeros(outshape), np.zeros(outshape), np.zeros(outshape)
    # Divide each pixel into subpixels
    spec_offs = np.arange(0.5/spec_subpixel, 1, 1/spec_subpixel) - 0.5  # -0.5 is to offset from the centre of each pixel.
    spat_offs = np.arange(0.5/spat_subpixel, 1, 1/spat_subpixel) - 0.5  # -0.5 is to offset from the centre of each pixel.
    slice_offs = np.arange(0.5/slice_subpixel, 1, 1/slice_subpixel) - 0.5  # -0.5 is to offset from the centre of each slice.
    spat_x, spec_y = np.meshgrid(spat_offs, spec_offs)
    num_subpixels = spec_subpixel * spat_subpixel  # Number of subpixels (spat & spec) per detector pixel
    num_all_subpixels = num_subpixels * slice_subpixel  # Number of subpixels, including slice subpixels
    # Loop through all exposures
    for fr in range(numframes):
        onslit_gpm = _gpmImg[fr]
        this_onslit_gpm = onslit_gpm > 0
        this_specpos, this_spatpos = np.where(this_onslit_gpm)
        this_spatid = onslit_gpm[this_onslit_gpm]

        # Extract tilts and slits for convenience
        this_tilts = _tilts[fr]
        this_slits = _slits[fr]
        this_wcs = _all_wcs[fr]
        this_astrom_trans = _astrom_trans[fr]
        this_wght_subpix = _wghtImg[fr][this_onslit_gpm]
        this_sci = _sciImg[fr][this_onslit_gpm]
        this_var = utils.inverse(_ivarImg[fr][this_onslit_gpm])
        this_wav = _waveImg[fr][this_onslit_gpm]
        slshape = (this_slits.nspec, this_slits.nspat, slice_subpixel)
        raimg_slc = np.zeros(slshape, dtype=float)
        decimg_slc = np.zeros(slshape, dtype=float)
        for ss in range(slice_subpixel):
            # Generate an RA/Dec image for this subslice
            raimg_slc[:, :, ss], decimg_slc[:, :, ss], _ = this_slits.get_radec_image(
                this_wcs, this_astrom_trans, this_tilts, slice_offset=slice_offs[ss]
            )
        # Loop through all slits
        for sl, spatid in enumerate(this_slits.spat_id):
            if verbose:
                if numframes == 1:
                    log.info(f"Resampling slit {sl + 1}/{this_slits.nslits}")
                else:
                    log.info(f"Resampling slit {sl + 1}/{this_slits.nslits} of frame {fr + 1}/{numframes}")
            # Find the pixels on this slit
            this_sl = np.where(this_spatid == spatid)
            wpix = (this_specpos[this_sl], this_spatpos[this_sl])
            # Create an array to index each subpixel
            numpix = wpix[0].size
            if numpix == 0:
                # Slit is masked or has no good pixels - skip
                continue
            # Generate a spline between spectral pixel position and wavelength
            yspl = this_tilts[wpix] * (this_slits.nspec - 1)
            tiltpos = np.add.outer(yspl, spec_y).flatten()
            wspl = this_wav[this_sl]
            asrt = np.argsort(yspl, kind='stable')
            wave_spl = interp1d(yspl[asrt], wspl[asrt], kind='linear', bounds_error=False, fill_value='extrapolate')
            # Calculate the wavelength at each subpixel
            this_wave_subpix = wave_spl(tiltpos)
            # Calculate the DAR correction at each sub pixel
            ra_corr, dec_corr = 0.0, 0.0
            if correct_dar:
                # NOTE :: This routine needs the wavelengths to be expressed in Angstroms
                ra_corr, dec_corr = _all_dar[fr].correction(this_wave_subpix)
            # Calculate spatial and spectral positions of the subpixels
            spat_xx = np.add.outer(wpix[1], spat_x.flatten()).flatten()
            spec_yy = np.add.outer(wpix[0], spec_y.flatten()).flatten()
            # Transform this to spatial location
            spatpos_subpix = _astrom_trans[fr].transform(sl, spat_xx, spec_yy)
            spatpos = _astrom_trans[fr].transform(sl, wpix[1], wpix[0])
            ssrt = np.argsort(spatpos, kind='stable')
            # Initialize the voxel coordinates for each spec2D pixel
            vox_coord = np.full((numpix, num_all_subpixels, 3), -1, dtype=float)
            # Loop over the subslices
            for ss in range(slice_subpixel):
                if verbose and slice_subpixel > 1: 
                    # Only print this if there are multiple subslices
                    log.info(f"Resampling subslice {ss+1}/{slice_subpixel}")
                # Select the RA/Dec image for this subslice
                this_ra = raimg_slc[:,:,ss][this_onslit_gpm]
                this_dec = decimg_slc[:,:,ss][this_onslit_gpm]
                # Interpolate the RA/Dec over the subpixel spatial positions
                tmp_ra = this_ra[this_sl]
                tmp_dec = this_dec[this_sl]
                # Evaluate the RA/Dec at the subpixel spatial positions
                this_ra_int = utils.linear_interpolate_extrapolate(spatpos_subpix, spatpos[ssrt], tmp_ra[ssrt])
                this_dec_int = utils.linear_interpolate_extrapolate(spatpos_subpix, spatpos[ssrt], tmp_dec[ssrt])
                # Now apply the DAR correction and any user-supplied offsets
                this_ra_int += ra_corr - _ra_offset[fr]
                this_dec_int += dec_corr - _dec_offset[fr]
                # Convert world coordinates to voxel coordinates, then histogram
                sslo = ss * num_subpixels
                sshi = (ss + 1) * num_subpixels
                vox_coord[:,sslo:sshi,:] = output_wcs.wcs_world2pix(np.vstack((this_ra_int, this_dec_int, this_wave_subpix * 1.0E-10)).T, 0).reshape(numpix, num_subpixels, 3)[:,:,::-1]
            # Convert the voxel coordinates to a bin index
            if num_all_subpixels == 1 or skip_subpix_weights:
                subpix_wght = 1.0
            else:
                if verbose: 
                    log.info("Preparing subpixel weights")
                vox_index = np.floor(vox_coord * voxscale - voxoffset).astype(int)
                # Convert to a unique index
                vox_index = np.dot(vox_index, np.array([1, outshape[0], outshape[0]*outshape[1]]))
                # Calculate the number of repeated indices for each subpixel - this is the subpixel weights
                subpix_wght = utils.occurrences_sorted(vox_index)
            # Reshape the voxel coordinates
            vox_coord = vox_coord.reshape(numpix * num_all_subpixels, 3)
            # Use the "fast histogram" algorithm, that assumes regular bin spacing
            flxcube += histogramdd(vox_coord, bins=outshape, range=binrng, weights=np.repeat(this_sci[this_sl] * this_wght_subpix[this_sl], num_all_subpixels) * subpix_wght)
            varcube += histogramdd(vox_coord, bins=outshape, range=binrng, weights=np.repeat(this_var[this_sl] * this_wght_subpix[this_sl]**2, num_all_subpixels) * subpix_wght**2)  # NOTE :: This was changed from subpix_wght**3 to subpix_wght**2 by RJC on 2024-12-18
            normcube += histogramdd(vox_coord, bins=outshape, range=binrng, weights=np.repeat(this_wght_subpix[this_sl], num_all_subpixels) * subpix_wght)

    # Normalise the datacube and variance cube
    nc_inverse = utils.inverse(normcube)
    flxcube *= nc_inverse
    varcube *= nc_inverse**2
    bpmcube = normcube == 0

    # Return the datacube, variance cube and bad pixel cube
    return flxcube, np.sqrt(varcube), bpmcube, normcube



def make_whitelight_fromcube_old(cube, bpmcube, wave=None, wavemin=None, wavemax=None):
    """
    Generate a white light image using an input cube.

    Args:
        cube (`numpy.ndarray`_):
            3D datacube (the final element contains the wavelength dimension)
        bpmcube (`numpy.ndarray`_, bool):
            3D bad pixel mask cube (the final element contains the wavelength dimension).
            A value of True indicates a bad pixel.
        wave (`numpy.ndarray`_, optional):
            1D wavelength array. Only required if wavemin or wavemax are not
            None.
        wavemin (float, optional):
            Minimum wavelength (same units as wave) to be included in the
            whitelight image.  You must provide wave as well if you want to
            reduce the wavelength range.
        wavemax (float, optional):
            Maximum wavelength (same units as wave) to be included in the
            whitelight image.  You must provide wave as well if you want to
            reduce the wavelength range.

    Returns:
        A whitelight image of the input cube (of type `numpy.ndarray`_).
    """
    # Make a wavelength cut, if requested
    if wavemin is not None or wavemax is not None:
        # Make some checks on the input
        if wave is None:
            raise PypeItError(
                "wave variable must be supplied to create white light image with wavelength cuts"
            )
        else:
            if wave.size != cube.shape[2]:
                raise PypeItError(
                    "wave variable should have the same length as the third axis of cube."
                )
        # assign wavemin & wavemax if one is not provided
        if wavemin is None:
            wavemin = np.min(wave)
        if wavemax is None:
            wavemax = np.max(wave)
        ww = np.where((wave >= wavemin) & (wave <= wavemax))[0]
        wmin, wmax = ww[0], ww[-1]+1
        cutcube = cube[:, :, wmin:wmax]
        # Cut the bad pixel mask and convert it to a good pixel mask
        cutgpmcube = np.logical_not(bpmcube[:, :, wmin:wmax])
    else:
        cutcube = cube.copy()
        cutgpmcube = np.logical_not(bpmcube)
    # Now sum along the wavelength axis
    nrmval = np.sum(cutgpmcube, axis=2)
    nrmval[nrmval == 0] = 1.0
    return np.sum(cutcube*cutgpmcube, axis=2) / nrmval
