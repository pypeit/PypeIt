"""
Build a datacube from Binospec IFU spec1d files.

Unlike the general ``pypeit_coadd_datacube`` script designed for slicer-based
IFUs, this script handles the fiber-fed Binospec IFU by:

1. Reading extracted 1D fiber spectra from spec1d files (already
   sky-subtracted by the pipeline)
2. Mapping 320 science fibers per side to sky positions
3. Combining both detectors (640 science fibers total)
4. Interpolating scattered fiber positions onto a regular spatial grid

Each input file produces a separate output datacube.

.. include:: ../include/links.rst
"""

from __future__ import annotations

import argparse
from typing import TYPE_CHECKING

import numpy as np

from pypeit.scripts import scriptbase

if TYPE_CHECKING:
    from pypeit.spectrographs.spectrograph import Spectrograph


class BinospecIFUCube(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width: int | None = None) -> argparse.ArgumentParser:
        parser = super().get_parser(
            description='Build a datacube from Binospec IFU spec1d files.',
            width=width,
            default_log_file=True)
        parser.add_argument('files', type=str, nargs='+',
                            help='One or more PypeIt spec1d files, or a text '
                                 'file listing them (one per line)')
        parser.add_argument('-o', '--output', type=str, default=None,
                            help='Output FITS filename (only valid for a single '
                                 'input file; default: auto-generated)')
        parser.add_argument('--spatial_scale', type=float, default=0.27,
                            help='Output spatial pixel scale in arcsec (default: 0.27)')
        parser.add_argument('--boxcar', default=False, action='store_true',
                            help='Use boxcar (BOX) extraction columns instead '
                                 'of the default optimal (OPT) columns')
        parser.add_argument('--method', type=str, default='linear',
                            choices=['nearest', 'linear', 'cubic'],
                            help='Spatial interpolation method (default: linear)')
        return parser

    @classmethod
    def main(cls, args: argparse.Namespace) -> None:
        import os

        from astropy.io import fits

        from pypeit import log, PypeItError
        from pypeit.spectrographs.util import load_spectrograph

        cls.init_log(args)

        # ----------------------------------------------------------------
        # Parse input files
        # ----------------------------------------------------------------
        input_files: list[str] = []
        for f in args.files:
            if f.endswith('.txt'):
                with open(f, 'r') as fh:
                    for line in fh:
                        line = line.strip()
                        if line and not line.startswith('#'):
                            input_files.append(line)
            else:
                input_files.append(f)

        if len(input_files) == 0:
            raise PypeItError("No input files provided.")

        if args.output is not None and len(input_files) > 1:
            raise PypeItError("--output can only be used with a single input file.")

        # Only spec1d files are supported
        for f in input_files:
            if not os.path.basename(f).startswith('spec1d'):
                raise PypeItError(
                    f"Only spec1d files are supported; got {os.path.basename(f)}")

        log.info(f"Processing {len(input_files)} spec1d file(s)")

        # Load spectrograph and fiber layout (shared across all files)
        with fits.open(input_files[0]) as hdu:
            spectrograph = load_spectrograph(hdu[0].header['PYP_SPEC'])

        targetx, targety = spectrograph.load_sky_layout()

        # ----------------------------------------------------------------
        # Process each file into a separate datacube
        # ----------------------------------------------------------------
        for input_file in input_files:
            log.info(f"Building datacube for {os.path.basename(input_file)}")
            _build_cube_from_spec1d(input_file, args, spectrograph,
                                    targetx, targety)


def _build_cube_from_spec1d(spec1d_file: str, args: argparse.Namespace,
                            spectrograph: Spectrograph,
                            targetx: np.ndarray, targety: np.ndarray) -> None:
    """Build a single datacube from one spec1d file.

    Reads the already-extracted 1D fiber spectra from a PypeIt spec1d file
    and builds a datacube.  Uses OPT (optimal) extraction by default, or
    BOX (boxcar) extraction if ``--boxcar`` is specified.  Sky is already
    subtracted by the pipeline, so no additional sky subtraction is
    performed here.

    Parameters
    ----------
    spec1d_file : str
        Path to the spec1d FITS file.
    args : `argparse.Namespace`_
        Parsed command-line arguments.
    spectrograph : :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
        Spectrograph instance.
    targetx : `numpy.ndarray`_
        Fiber x positions on sky (arcsec).
    targety : `numpy.ndarray`_
        Fiber y positions on sky (arcsec).
    """
    import os

    import numpy as np
    from astropy.io import fits

    from pypeit import log, PypeItError
    from pypeit.specobjs import SpecObjs

    sobjs = SpecObjs.from_fitsfile(spec1d_file)
    if sobjs.nobj == 0:
        log.warning(f"No objects in {os.path.basename(spec1d_file)}, skipping")
        return

    # Choose extraction type
    prefix = 'BOX' if args.boxcar else 'OPT'
    log.info(f"  Using {prefix} extraction from spec1d")

    # ------------------------------------------------------------------
    # Step 1: Organize fiber spectra by detector
    # ------------------------------------------------------------------
    det_fiber_data = {}

    for det_name in ['DET01', 'DET02']:
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
                    f"{os.path.basename(spec1d_file)}, skipping")
        return

    # ------------------------------------------------------------------
    # Build the datacube
    # ------------------------------------------------------------------
    with fits.open(spec1d_file) as hdu:
        raw_hdr = hdu[0].header

    _build_cube_common(det_fiber_data, args, spectrograph,
                       targetx, targety, raw_hdr, spec1d_file)


def _build_cube_common(det_fiber_data: dict, args: argparse.Namespace,
                       spectrograph: Spectrograph,
                       targetx: np.ndarray, targety: np.ndarray,
                       raw_hdr, input_file: str) -> None:
    """Shared steps for building a datacube from fiber spectra.

    The input fiber spectra are assumed to be already sky-subtracted
    (as produced by the PypeIt pipeline in the spec1d files).

    Parameters
    ----------
    det_fiber_data : dict
        Per-detector fiber data. Keys are detector names (e.g. ``'DET01'``),
        values are dicts with keys ``'flux'``, ``'ivar'``, ``'wave'``,
        ``'fiber_meta'``.
    args : `argparse.Namespace`_
        Parsed command-line arguments.
    spectrograph : :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
        Spectrograph instance.
    targetx : `numpy.ndarray`_
        Fiber x positions on sky (arcsec).
    targety : `numpy.ndarray`_
        Fiber y positions on sky (arcsec).
    raw_hdr : `astropy.io.fits.Header`_
        Primary header from the input file.
    input_file : str
        Path to the input file (for generating the output filename).
    """
    import os

    import numpy as np
    from astropy import units
    from astropy.io import fits
    from astropy import wcs
    from scipy.interpolate import griddata
    from scipy.spatial import QhullError

    from pypeit import log

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

    # Resample each fiber onto the common wavelength grid.  Shared with the
    # 1D fiber extractor: interpolates variance (not inverse variance) so the
    # propagated uncertainty is statistically correct.
    from pypeit.core import datacube
    for data in det_fiber_data.values():
        nfibers = data['flux'].shape[0]
        flux_resamp = np.zeros((nfibers, n_wave))
        ivar_resamp = np.zeros((nfibers, n_wave))

        for i in range(nfibers):
            flux_resamp[i], ivar_resamp[i], _ = datacube.resample_spec_to_grid(
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
    # IDL scaling: positions / 3.0, then offset and scale
    scl = args.spatial_scale
    x_scaled = fiber_x / 3.0
    y_scaled = fiber_y / 3.0

    # Compute grid dimensions from scaled positions
    x_min, x_max = np.min(x_scaled), np.max(x_scaled)
    y_min, y_max = np.min(y_scaled), np.max(y_scaled)

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
    points = np.column_stack([x_scaled, y_scaled])
    cube = np.zeros((nx, ny, n_wave), dtype=np.float32)
    var_cube = np.zeros((nx, ny, n_wave), dtype=np.float32)

    log.info(f"Interpolating {n_wave} wavelength slices using "
              f"method='{args.method}'...")
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
                method=args.method, fill_value=0.0)
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
                    method=args.method, fill_value=0.0)
            except QhullError:
                var_cube[:, :, k] = griddata(
                    points[good], var_slice[good], (grid_x, grid_y),
                    method='nearest', fill_value=0.0)

    # ------------------------------------------------------------------
    # Step 7: Build WCS and write output
    # ------------------------------------------------------------------
    # Pointing and celestial CD matrix (shared convention with the 1D
    # fiber extractor; see MMTBINOSPECSpectrograph.ifu_sky_wcs).
    coord, cd = spectrograph.ifu_sky_wcs(raw_hdr, scl)
    (cd11, cd12), (cd21, cd22) = cd

    w = wcs.WCS(naxis=3)
    w.wcs.equinox = raw_hdr.get('EQUINOX', 2000.0)
    w.wcs.name = 'Binospec IFU'
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

    # Build output FITS
    hdr = w.to_header()
    hdr['INSTRUME'] = 'BINOSPEC'
    hdr['TELESCOP'] = 'MMT'
    hdr['IFUMODE'] = 'FIBER'
    hdr['NFIBERS'] = (n_sci_fibers, 'Number of science fibers')
    hdr['SPATSCL'] = (scl, 'Spatial pixel scale [arcsec]')
    hdr['WAVEMIN'] = (wave_grid[0], 'Minimum wavelength [Angstrom]')
    hdr['WAVEMAX'] = (wave_grid[-1], 'Maximum wavelength [Angstrom]')
    hdr['WAVESTP'] = (dwv, 'Wavelength step [Angstrom]')
    hdr['INTERP'] = (args.method, 'Spatial interpolation method')

    # Copy useful keywords from raw header
    for key in ['OBJECT', 'EXPTIME', 'DATE-OBS', 'DISPERSE', 'FILTER']:
        if key in raw_hdr:
            hdr[key] = raw_hdr[key]

    # Output filename
    if args.output is not None:
        outfile = args.output
    else:
        base = os.path.splitext(os.path.basename(input_file))[0]
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
