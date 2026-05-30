"""Interactive 1D fiber spectrum extractor for the MMT Binospec IFU.

Loads a single PypeIt ``spec1d`` file (covering both detector sides),
displays the IFU as a hexagonal flux map, and writes the user's combined
fiber selection as a :class:`~pypeit.onespec.OneSpec` FITS file.
"""
from __future__ import annotations

import argparse
from typing import TYPE_CHECKING

import numpy as np

from pypeit import PypeItError
from pypeit.scripts import scriptbase

if TYPE_CHECKING:
    from astropy.io.fits import Header

# Bright optical sky lines (vacuum, Angstrom).  Used only for the hex-display
# integrated flux so that residual sky-subtraction errors do not bias the
# per-fiber flux estimate.  Extraction itself uses the full wavelength range.
_SKY_LINE_WAVELENGTHS: tuple[float, ...] = (
    5577.34,                       # [OI]
    5890.00, 5896.00,              # Na I D
    6300.30, 6363.78,              # [OI]
    7340.0, 7821.0, 8344.0, 8430.0, 8761.0, 8867.0,  # OH band centers (approximate)
)
_SKY_LINE_MASK_HALFWIDTH: float = 10.0  # Angstrom


def _sky_line_mask(wave: np.ndarray) -> np.ndarray:
    """Return a boolean mask that is ``True`` inside any sky-line region.

    Parameters
    ----------
    wave : `numpy.ndarray`_
        1D array of wavelengths in Angstrom.

    Returns
    -------
    `numpy.ndarray`_
        1D boolean array, same shape as ``wave``, ``True`` where the
        wavelength falls within ``_SKY_LINE_MASK_HALFWIDTH`` of any line
        in ``_SKY_LINE_WAVELENGTHS``.
    """
    mask = np.zeros_like(wave, dtype=bool)
    for lam in _SKY_LINE_WAVELENGTHS:
        mask |= np.abs(wave - lam) < _SKY_LINE_MASK_HALFWIDTH
    return mask


def _project_to_sky(x_arcsec: np.ndarray, y_arcsec: np.ndarray,
                    raw_hdr: 'Header') -> tuple[np.ndarray, np.ndarray]:
    """Project instrument-frame fiber offsets to sky coordinates.

    Uses the same TAN WCS convention as
    :class:`~pypeit.scripts.binospec_ifu_cube.BinospecIFUCube` so that the
    extracted hex view aligns with the matching datacube.

    Parameters
    ----------
    x_arcsec, y_arcsec : `numpy.ndarray`_
        1D arrays of fiber offsets in instrument-frame arcseconds.
    raw_hdr : `astropy.io.fits.Header`_
        Primary header from the spec1d file, providing ``RA``, ``DEC``,
        and (optionally) ``POSANG``.

    Returns
    -------
    ra, dec : `numpy.ndarray`_
        Fiber RA/Dec in degrees, same shape as the inputs.
    """
    from astropy import units as u
    from astropy import wcs as astropy_wcs
    from astropy.coordinates import SkyCoord

    raval = raw_hdr.get('RA', 0.0)
    decval = raw_hdr.get('DEC', 0.0)
    try:
        coord = SkyCoord(raval, decval, unit=(u.hourangle, u.deg))
    except Exception:
        coord = SkyCoord(raval, decval, unit=(u.deg, u.deg))

    posang = raw_hdr.get('POSANG', 0.0)
    crota = np.radians(-posang)

    # 1 deg per pixel: x_arcsec/3600 becomes a "pixel" offset in the WCS.
    cdelt1 = -1.0 / 3600.0  # RA decreases with +x (east is +RA, +x is east)
    cdelt2 = 1.0 / 3600.0   # Dec increases with +y

    cd11 = cdelt1 * np.cos(crota)
    cd12 = abs(cdelt2) * np.sign(cdelt1) * np.sin(crota)
    cd21 = -abs(cdelt1) * np.sign(cdelt2) * np.sin(crota)
    cd22 = cdelt2 * np.cos(crota)

    w = astropy_wcs.WCS(naxis=2)
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.crval = [coord.ra.degree, coord.dec.degree]
    # crpix is FITS 1-indexed; pixel_to_world is 0-indexed.  Setting crpix=1
    # makes pixel (0, 0) == reference pixel == crval.
    w.wcs.crpix = [1.0, 1.0]
    w.wcs.cd = np.array([[cd11, cd12], [cd21, cd22]])
    w.wcs.lonpole = 180.0
    w.wcs.latpole = 0.0

    sky = w.pixel_to_world(np.asarray(x_arcsec, dtype=float),
                            np.asarray(y_arcsec, dtype=float))
    return np.atleast_1d(sky.ra.degree), np.atleast_1d(sky.dec.degree)


def _compute_fiber_fluxes(waves: list[np.ndarray],
                          fluxes: list[np.ndarray],
                          wave_min: float, wave_max: float) -> np.ndarray:
    """Per-fiber integrated flux with sky-line masking.

    Each fiber is integrated over its **own** native wavelength grid, with
    pixels in any sky-line region (see ``_SKY_LINE_WAVELENGTHS``) excluded.
    Used to colour the hex display; the actual extracted spectrum
    (:func:`_resample_and_combine`) is unaffected.

    Parameters
    ----------
    waves, fluxes : list of `numpy.ndarray`_
        Per-fiber native wavelengths and fluxes (parallel lists).
    wave_min, wave_max : float
        Slider-defined integration range, Angstrom.

    Returns
    -------
    `numpy.ndarray`_
        1D array of integrated fluxes, one per fiber.
    """
    out = np.zeros(len(waves), dtype=float)
    for i, (w, f) in enumerate(zip(waves, fluxes)):
        in_range = (w >= wave_min) & (w <= wave_max)
        sky = _sky_line_mask(w)
        mask = in_range & ~sky
        out[i] = np.nansum(f[mask])
    return out


def _resample_and_combine(waves: list[np.ndarray],
                          fluxes: list[np.ndarray],
                          ivars: list[np.ndarray]
                          ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Resample selected fibers onto a common grid and combine.

    Parameters
    ----------
    waves, fluxes, ivars : list of `numpy.ndarray`_
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
    dwv = float(np.median(diffs))
    n_wave = int(round((wave_max - wave_min) / dwv)) + 1
    wave_out = np.linspace(wave_min, wave_max, n_wave)

    n_fib = len(valid)
    flux_resamp = np.zeros((n_fib, n_wave))
    var_resamp = np.zeros((n_fib, n_wave))
    have_flux = np.zeros((n_fib, n_wave), dtype=bool)
    have_var = np.zeros((n_fib, n_wave), dtype=bool)

    for i, (w, f, iv) in enumerate(valid):
        good = (w > 0) & np.isfinite(f)
        if good.sum() < 2:
            continue
        srt = np.argsort(w[good])
        w_s = w[good][srt]
        f_s = f[good][srt]
        iv_s = iv[good][srt] if iv is not None else np.zeros_like(w_s)
        in_range = (wave_out >= w_s[0]) & (wave_out <= w_s[-1])
        flux_resamp[i, in_range] = np.interp(wave_out[in_range], w_s, f_s)
        have_flux[i, in_range] = True
        # Variance combine: build var on native grid only where ivar > 0,
        # then interp onto common grid.
        var_native = np.where(iv_s > 0, 1.0 / np.maximum(iv_s, 1e-300), 0.0)
        var_resamp[i, in_range] = np.interp(wave_out[in_range], w_s,
                                            var_native)
        have_var[i, in_range] = (var_resamp[i, in_range] > 0)

    n_good = have_flux.sum(axis=0)
    n_total = n_fib
    rescale = n_total / np.maximum(n_good, 1)
    flux_out = np.nansum(np.where(have_flux, flux_resamp, 0.0),
                         axis=0) * rescale

    var_out = np.nansum(np.where(have_var, var_resamp, 0.0), axis=0)
    ivar_out = np.where(var_out > 0, 1.0 / var_out, 0.0)

    return wave_out, flux_out, ivar_out


def _load_fibers(sobjs, spectrograph, targetx: np.ndarray,
                 targety: np.ndarray, prefix: str = 'OPT') -> list[dict]:
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
    prefix : {'OPT', 'BOX'}
        Extraction column prefix (optimal or boxcar).

    Returns
    -------
    list of dict
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

    fibers: list[dict] = []
    for det_name in ['DET01', 'DET02']:
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


def _write_onespec(wave: np.ndarray, flux: np.ndarray, ivar: np.ndarray,
                   raw_hdr: 'Header', pyp_spec: str, outfile: str) -> None:
    """Write a combined fiber spectrum as a OneSpec FITS file.

    Parameters
    ----------
    wave, flux, ivar : `numpy.ndarray`_
        1D arrays of equal length (output wavelength grid, summed flux,
        inverse variance).
    raw_hdr : `astropy.io.fits.Header`_
        Primary header of the source spec1d, copied through to the output.
    pyp_spec : str
        Spectrograph short name, written as ``PYP_SPEC`` so downstream
        PypeIt scripts can re-load the file.
    outfile : str
        Output FITS path.
    """
    from pypeit.onespec import OneSpec

    # `wave_grid_mid` is the (optional) uniformly-spaced grid used by 1D
    # coaddition; this writer is not a coadd, so leave it as None.
    one = OneSpec(wave=wave, wave_grid_mid=None, flux=flux, ivar=ivar,
                  PYP_SPEC=pyp_spec)
    # Do NOT set `one.head0` before `to_file` -- that triggers
    # `spectrograph.subheader_for_spec` which expects extra metadata keys
    # (e.g. target tables) not present in arbitrary input headers.  Pass
    # the raw header via `primary_hdr=` to copy keys through.
    one.to_file(outfile, primary_hdr=raw_hdr, overwrite=True)


class BinospecIFUExtract(scriptbase.ScriptBase):
    """Interactive fiber-spectrum extractor for the MMT Binospec IFU."""

    @classmethod
    def get_parser(cls, width: int | None = None) -> argparse.ArgumentParser:
        parser = super().get_parser(
            description='Interactive 1D fiber spectrum extractor for the '
                        'MMT Binospec IFU.',
            width=width,
            default_log_file=True)
        parser.add_argument('spec1d_file', type=str,
                            help='PypeIt spec1d FITS file (covering both '
                                 'detector sides)')
        parser.add_argument('-o', '--output', type=str, default=None,
                            help='Output OneSpec FITS filename '
                                 '(default: spec1d_<base>.fits → '
                                 'extract1d_<base>.fits)')
        parser.add_argument('--boxcar', default=False, action='store_true',
                            help='Use boxcar (BOX) extraction columns '
                                 'instead of optimal (OPT)')
        return parser

    @classmethod
    def main(cls, args: argparse.Namespace) -> None:
        import os

        from astropy.io import fits

        from pypeit import log
        from pypeit.specobjs import SpecObjs
        from pypeit.spectrographs.util import load_spectrograph

        cls.init_log(args)

        if not os.path.basename(args.spec1d_file).startswith('spec1d_'):
            raise PypeItError(
                f"Only spec1d files are supported; got "
                f"{os.path.basename(args.spec1d_file)}")

        log.info(f"Loading {args.spec1d_file}")
        sobjs = SpecObjs.from_fitsfile(args.spec1d_file)
        if sobjs.nobj == 0:
            raise PypeItError(f"No objects in {args.spec1d_file}")

        with fits.open(args.spec1d_file) as hdu:
            raw_hdr = hdu[0].header.copy()

        spectrograph = load_spectrograph(raw_hdr['PYP_SPEC'])
        if spectrograph.name != 'mmt_binospec_ifu':
            raise PypeItError(
                f"This script requires an mmt_binospec_ifu spec1d; "
                f"got PYP_SPEC={spectrograph.name!r}")
        targetx, targety = spectrograph.load_sky_layout()

        prefix = 'BOX' if args.boxcar else 'OPT'
        log.info(f"  Using {prefix} extraction")

        fibers = _load_fibers(sobjs, spectrograph, targetx, targety,
                              prefix=prefix)
        log.info(f"  Loaded {len(fibers)} science fibers")

        # Project each fiber's instrument-frame x,y to sky.
        x = np.array([f['x'] for f in fibers])
        y = np.array([f['y'] for f in fibers])
        ra, dec = _project_to_sky(x, y, raw_hdr)
        for f, r, d in zip(fibers, ra, dec):
            f['ra'] = float(r)
            f['dec'] = float(d)

        # Output filename
        if args.output is not None:
            outfile = args.output
        else:
            base = os.path.splitext(os.path.basename(args.spec1d_file))[0]
            if 'spec1d_' in base:
                base = base.replace('spec1d_', 'extract1d_')
            outfile = base + '.fits'

        gui = _ExtractGUI(fibers, raw_hdr, spectrograph.name, outfile)
        gui.show()


class _ExtractGUI:
    """Matplotlib widget GUI for selecting fibers and extracting spectra."""

    def __init__(self, fibers: list[dict], raw_hdr: 'Header', pyp_spec: str,
                 outfile: str):
        # Local imports; matplotlib is heavy and we only need it when the GUI
        # is actually launched (not on --help or module import).
        import matplotlib.pyplot as plt
        from matplotlib.collections import PatchCollection
        from matplotlib.patches import Circle, RegularPolygon, Rectangle
        from matplotlib.widgets import Button, RadioButtons, RangeSlider

        self._plt = plt
        self._PatchCollection = PatchCollection
        self._RegularPolygon = RegularPolygon
        self._Rectangle = Rectangle
        self._Circle = Circle
        self._Button = Button
        self._RadioButtons = RadioButtons
        self._RangeSlider = RangeSlider

        self.fibers = fibers
        self.raw_hdr = raw_hdr
        self.pyp_spec = pyp_spec
        self.outfile = outfile

        self.waves = [f['wave'] for f in fibers]
        self.fluxes = [f['flux'] for f in fibers]
        self.ivars = [f['ivar'] for f in fibers]
        self.ra = np.array([f['ra'] for f in fibers])
        self.dec = np.array([f['dec'] for f in fibers])
        self.dec_center = float(np.median(self.dec))
        self.cos_dec = np.cos(np.radians(self.dec_center))
        self.ra_scaled = self.ra * self.cos_dec
        self.nspec = len(fibers)

        # Slider extents.  Filter to fibers with at least one positive
        # wavelength so an all-zero/masked fiber does not crash np.min.
        valid_waves = [w[w > 0] for w in self.waves if np.any(w > 0)]
        if not valid_waves:
            raise PypeItError("No valid wavelength data in any fiber")
        self.wave_min = float(np.min([w.min() for w in valid_waves]))
        self.wave_max = float(np.max([w.max() for w in valid_waves]))
        self.current_wave_min = self.wave_min
        self.current_wave_max = self.wave_max

        self.fiber_fluxes = _compute_fiber_fluxes(
            self.waves, self.fluxes, self.wave_min, self.wave_max)

        # Selection state
        self.y_scale = 'linear'
        self.selection_shape = 'rectangle'
        self.selection_start = None
        self.selection_patch = None
        self.selected_region = None
        self.clicked_fibers: set[int] = set()
        self.fiber_mask: np.ndarray | None = None
        self.extracted_wave: np.ndarray | None = None
        self.extracted_flux: np.ndarray | None = None
        self.extracted_ivar: np.ndarray | None = None

        self._create_gui()

    def show(self) -> None:
        """Display the interactive GUI."""
        self._plt.show()

    def _create_gui(self) -> None:
        """Build all GUI axes, widgets, and event connections."""
        self.fig = self._plt.figure(figsize=(16, 8))
        self.ax_image = self._plt.subplot2grid((4, 3), (0, 0), rowspan=4,
                                               colspan=2)
        self.ax_spectrum = self._plt.subplot2grid((4, 3), (0, 2), rowspan=2,
                                                  colspan=1)
        self.ax_controls = self._plt.subplot2grid((4, 3), (2, 2), colspan=1)
        self.ax_coords = self._plt.subplot2grid((4, 3), (3, 2), colspan=1)
        self.ax_controls.axis('off')
        self.ax_coords.axis('off')

        self._create_hexagon_plot()

        self.ax_spectrum.set_xlabel('Wavelength (Å)')
        self.ax_spectrum.set_ylabel('Flux')
        self.ax_spectrum.set_title('Extracted Spectrum')
        self.ax_spectrum.grid(True, alpha=0.3)

        # Buttons
        self.btn_extract = self._Button(
            self._plt.axes([0.77, 0.34, 0.15, 0.05]), 'Extract Spectrum')
        self.btn_extract.on_clicked(self._on_extract)
        self.btn_reset = self._Button(
            self._plt.axes([0.77, 0.28, 0.15, 0.05]), 'Reset Display')
        self.btn_reset.on_clicked(self._on_reset)
        self.btn_save = self._Button(
            self._plt.axes([0.77, 0.22, 0.15, 0.05]), 'Save Spectrum')
        self.btn_save.on_clicked(self._on_save)

        # Mode radio
        self.radio = self._RadioButtons(
            self._plt.axes([0.68, 0.22, 0.08, 0.10]),
            ('Rectangle', 'Circle', 'Click Fibers'))
        self.radio.on_clicked(self._on_shape_change)

        # Y-axis scale radio for the extracted-spectrum panel
        self.yscale_radio = self._RadioButtons(
            self._plt.axes([0.77, 0.13, 0.08, 0.08]),
            ('Linear', 'Log'))
        self.yscale_radio.on_clicked(self._on_yscale_change)

        # Wavelength range slider
        self.wave_slider = self._RangeSlider(
            self._plt.axes([0.1, 0.05, 0.45, 0.03]), 'Wavelength (Å)',
            self.wave_min, self.wave_max,
            valinit=(self.wave_min, self.wave_max),
            valstep=((self.wave_max - self.wave_min) / 100) or None)
        self.wave_slider.on_changed(self._on_wave_change)

        # Mouse events
        self.fig.canvas.mpl_connect('button_press_event', self._on_press)
        self.fig.canvas.mpl_connect('motion_notify_event', self._on_motion)
        self.fig.canvas.mpl_connect('button_release_event', self._on_release)

        self._plt.tight_layout()
        self._plt.subplots_adjust(bottom=0.10)

    def _create_hexagon_plot(self) -> None:
        """Populate ``ax_image`` with a hexagonal fiber array flux map."""
        from matplotlib.ticker import FuncFormatter

        # Hexagon size from median nearest-neighbour separation
        distances = []
        for i in range(min(50, self.nspec)):
            dra = self.ra_scaled - self.ra_scaled[i]
            ddec = self.dec - self.dec[i]
            d = np.sqrt(dra**2 + ddec**2)
            d = d[d > 0]
            if d.size:
                distances.append(d.min())
        sep = float(np.median(distances)) if distances else 1.0 / 3600.0
        hex_radius = sep * 0.58

        valid = self.fiber_fluxes[np.isfinite(self.fiber_fluxes)
                                  & (self.fiber_fluxes != 0)]
        vmin = float(np.percentile(valid, 5)) if valid.size else 0.0
        vmax = float(np.percentile(valid, 95)) if valid.size else 1.0

        patches = [self._RegularPolygon((self.ra_scaled[i], self.dec[i]),
                                        numVertices=6, radius=hex_radius,
                                        orientation=np.pi / 2)
                   for i in range(self.nspec)]
        self.hex_collection = self._PatchCollection(patches, cmap='viridis',
                                                    edgecolors='gray',
                                                    linewidths=0.5)
        self.hex_collection.set_array(self.fiber_fluxes)
        self.hex_collection.set_clim(vmin, vmax)
        self.ax_image.add_collection(self.hex_collection)

        # nspec >= 1 is guaranteed by _load_fibers raising PypeItError on
        # empty input; otherwise ra_scaled.min() would crash here.
        pad = 5 * hex_radius
        self.ax_image.set_xlim(self.ra_scaled.min() - pad,
                               self.ra_scaled.max() + pad)
        self.ax_image.set_ylim(self.dec.min() - pad, self.dec.max() + pad)
        self.ax_image.set_aspect('equal')
        self.ax_image.xaxis.set_major_formatter(FuncFormatter(
            lambda x, pos: f'{x / self.cos_dec:.5f}'))
        self.ax_image.set_xlabel('RA (deg)')
        self.ax_image.set_ylabel('Dec (deg)')
        self.ax_image.set_title('Fiber Array (drag to select region)')
        self._plt.colorbar(self.hex_collection, ax=self.ax_image,
                           label='Integrated Flux')

    # ------------------------------------------------------------------
    # Callbacks
    # ------------------------------------------------------------------

    def _on_wave_change(self, val) -> None:
        wave_min, wave_max = val
        self.current_wave_min = float(wave_min)
        self.current_wave_max = float(wave_max)
        self.fiber_fluxes = _compute_fiber_fluxes(
            self.waves, self.fluxes,
            self.current_wave_min, self.current_wave_max)
        valid = self.fiber_fluxes[np.isfinite(self.fiber_fluxes)
                                   & (self.fiber_fluxes != 0)]
        if valid.size:
            self.hex_collection.set_clim(np.percentile(valid, 5),
                                          np.percentile(valid, 95))
        self.hex_collection.set_array(self.fiber_fluxes)
        if self.extracted_wave is not None:
            self._apply_spectrum_xlim()
        self.fig.canvas.draw()

    def _apply_spectrum_xlim(self) -> None:
        """Sync the spectrum axis to the current slider window."""
        if self.extracted_wave is None:
            return
        lo, hi = self.current_wave_min, self.current_wave_max
        self.ax_spectrum.set_xlim(lo, hi)
        in_window = ((self.extracted_wave >= lo)
                     & (self.extracted_wave <= hi))
        flux = self.extracted_flux[in_window]
        flux = flux[np.isfinite(flux)]
        if self.y_scale == 'log':
            # Log axis cannot represent non-positive flux; ignore it when
            # setting the limits (it is also masked from the plotted data).
            flux = flux[flux > 0]
            if flux.size:
                fmin, fmax = float(flux.min()), float(flux.max())
                self.ax_spectrum.set_ylim(fmin / 1.5, fmax * 1.5)
        elif flux.size:
            fmin, fmax = float(flux.min()), float(flux.max())
            pad = 0.05 * (fmax - fmin) if fmax > fmin else max(abs(fmax), 1.0)
            self.ax_spectrum.set_ylim(fmin - pad, fmax + pad)

    def _selected_indices(self) -> np.ndarray:
        mask = self._get_fiber_mask()
        for i in self.clicked_fibers:
            mask[i] = True
        return np.where(mask)[0]

    def _on_extract(self, event) -> None:
        from pypeit import log

        idxs = self._selected_indices()
        if idxs.size == 0:
            log.info("No fibers selected")
            return
        log.info(f"Extracting from {idxs.size} fibers")
        try:
            wave, flux, ivar = _resample_and_combine(
                [self.waves[i] for i in idxs],
                [self.fluxes[i] for i in idxs],
                [self.ivars[i] for i in idxs])
        except PypeItError as e:
            log.warning(str(e))
            return
        self.extracted_wave = wave
        self.extracted_flux = flux
        self.extracted_ivar = ivar
        self.fiber_mask = np.zeros(self.nspec, dtype=bool)
        self.fiber_mask[idxs] = True

        self._draw_spectrum()
        self.fig.canvas.draw()

    def _draw_spectrum(self) -> None:
        """Redraw the extracted spectrum honoring the current y-axis scale.

        On a log scale, non-positive flux values (e.g. sky-subtraction
        residuals) are masked so they neither appear in the trace nor
        distort the axis limits.
        """
        if self.extracted_wave is None:
            return
        wave = self.extracted_wave
        flux = self.extracted_flux
        ivar = self.extracted_ivar

        # Clamp tiny-ivar pixels so the envelope doesn't blow up the y-axis.
        sigma_raw = np.where(ivar > 0,
                             1.0 / np.sqrt(np.maximum(ivar, 1e-300)),
                             np.nan)
        med_sig = float(np.nanmedian(sigma_raw)) if np.any(ivar > 0) else 0.0
        sigma = np.where(np.isfinite(sigma_raw)
                          & (sigma_raw < 10.0 * med_sig),
                          sigma_raw, 0.0)

        plot_flux = flux
        lower = flux - sigma
        upper = flux + sigma
        if self.y_scale == 'log':
            # Mask non-positive flux: a log axis cannot represent it.  The
            # lower error envelope is clamped to a small positive value so
            # matplotlib does not drop the band where flux > 0 but
            # flux - sigma <= 0.
            bad = ~(flux > 0)
            plot_flux = np.where(bad, np.nan, flux)
            lower = np.where(bad, np.nan, np.maximum(lower, 1e-300))
            upper = np.where(bad, np.nan, upper)

        self.ax_spectrum.clear()
        self.ax_spectrum.set_yscale(self.y_scale)
        self.ax_spectrum.plot(wave, plot_flux, lw=0.8)
        self.ax_spectrum.fill_between(wave, lower, upper, alpha=0.3)
        self.ax_spectrum.set_xlabel('Wavelength (Å)')
        self.ax_spectrum.set_ylabel('Flux')
        n_fib = int(self.fiber_mask.sum()) if self.fiber_mask is not None else 0
        self.ax_spectrum.set_title(f'Extracted Spectrum ({n_fib} fibers)')
        self.ax_spectrum.grid(True, alpha=0.3)
        self._apply_spectrum_xlim()

    def _on_yscale_change(self, label: str) -> None:
        self.y_scale = label.lower()
        if self.extracted_wave is None:
            self.ax_spectrum.set_yscale(self.y_scale)
        else:
            self._draw_spectrum()
        self.fig.canvas.draw()

    def _on_save(self, event) -> None:
        from pypeit import log

        if self.extracted_wave is None:
            log.info("Nothing to save -- run Extract Spectrum first")
            return
        _write_onespec(self.extracted_wave, self.extracted_flux,
                       self.extracted_ivar, self.raw_hdr, self.pyp_spec,
                       self.outfile)
        log.info(f"Wrote {self.outfile}")

    def _on_reset(self, event) -> None:
        self.selected_region = None
        self.clicked_fibers.clear()
        self.selection_start = None
        if self.selection_patch:
            self.selection_patch.remove()
            self.selection_patch = None

        # Restore wavelength range to full extent.
        self.current_wave_min = self.wave_min
        self.current_wave_max = self.wave_max
        self.wave_slider.set_val((self.wave_min, self.wave_max))
        self.fiber_fluxes = _compute_fiber_fluxes(
            self.waves, self.fluxes, self.wave_min, self.wave_max)
        valid = self.fiber_fluxes[np.isfinite(self.fiber_fluxes)
                                   & (self.fiber_fluxes != 0)]
        vmin = float(np.percentile(valid, 5)) if valid.size else 0.0
        vmax = float(np.percentile(valid, 95)) if valid.size else 1.0
        self.hex_collection.set_array(self.fiber_fluxes)
        self.hex_collection.set_clim(vmin, vmax)

        # Reset highlighting.
        self.hex_collection.set_edgecolors('gray')
        self.hex_collection.set_linewidths(0.5)

        # Clear spectrum panel and extracted data.
        self.ax_spectrum.clear()
        self.ax_spectrum.set_yscale(self.y_scale)
        self.ax_spectrum.set_xlabel('Wavelength (Å)')
        self.ax_spectrum.set_ylabel('Flux')
        self.ax_spectrum.set_title('Extracted Spectrum')
        self.ax_spectrum.grid(True, alpha=0.3)
        self.extracted_wave = None
        self.extracted_flux = None
        self.extracted_ivar = None
        self.fiber_mask = None

        self.fig.canvas.draw()

    def _on_shape_change(self, label: str) -> None:
        self.selection_shape = label.lower()
        if self.selection_patch:
            self.selection_patch.remove()
            self.selection_patch = None
        if self.selection_shape == 'click fibers':
            self.ax_image.set_title('Fiber Array (click hexagons to add/remove)')
        else:
            self.ax_image.set_title('Fiber Array (drag to select region)')
        self.fig.canvas.draw()

    def _on_press(self, event) -> None:
        if event.inaxes != self.ax_image:
            return
        if self.selection_shape == 'click fibers':
            dra = self.ra_scaled - event.xdata
            ddec = self.dec - event.ydata
            d = np.sqrt(dra**2 + ddec**2)
            i = int(np.argmin(d))
            if d[i] < 3.0 / 3600.0:
                if i in self.clicked_fibers:
                    self.clicked_fibers.remove(i)
                else:
                    self.clicked_fibers.add(i)
                self._update_selection_display()
            return
        self.selection_start = (event.xdata, event.ydata)
        if self.selection_patch:
            self.selection_patch.remove()
            self.selection_patch = None

    def _on_motion(self, event) -> None:
        if (self.selection_shape == 'click fibers'
                or self.selection_start is None
                or event.inaxes != self.ax_image):
            return
        if self.selection_patch:
            self.selection_patch.remove()
            self.selection_patch = None
        x0, y0 = self.selection_start
        x1, y1 = event.xdata, event.ydata
        if self.selection_shape == 'rectangle':
            self.selection_patch = self._Rectangle(
                (x0, y0), x1 - x0, y1 - y0,
                fill=False, edgecolor='cyan', linewidth=2, linestyle='--')
        else:
            r = np.sqrt((x1 - x0)**2 + (y1 - y0)**2)
            self.selection_patch = self._Circle(
                (x0, y0), r, fill=False, edgecolor='cyan',
                linewidth=2, linestyle='--')
        self.ax_image.add_patch(self.selection_patch)
        # Live preview: temporarily compute the would-be `selected_region`
        # so `_update_selection_display` lights up the fibers under the
        # current rubber-band.  Restore the prior region so release is the
        # only place that commits.
        old_region = self.selected_region
        if self.selection_shape == 'rectangle':
            x0_true, x1_true = x0 / self.cos_dec, x1 / self.cos_dec
            self.selected_region = {
                'type': 'rectangle',
                'ra_min': min(x0_true, x1_true),
                'ra_max': max(x0_true, x1_true),
                'dec_min': min(y0, y1),
                'dec_max': max(y0, y1),
            }
        else:
            r_scaled = np.sqrt((x1 - x0)**2 + (y1 - y0)**2)
            self.selected_region = {
                'type': 'circle',
                'ra_center': x0 / self.cos_dec,
                'dec_center': y0,
                'radius_scaled': r_scaled,
            }
        self._update_selection_display()
        self.selected_region = old_region

    def _on_release(self, event) -> None:
        if (self.selection_shape == 'click fibers'
                or self.selection_start is None
                or event.inaxes != self.ax_image):
            self.selection_start = None
            if self.selection_patch:
                self.selection_patch.remove()
                self.selection_patch = None
                self.fig.canvas.draw()
            return
        x0_s, y0 = self.selection_start
        x1_s, y1 = event.xdata, event.ydata
        x0 = x0_s / self.cos_dec
        x1 = x1_s / self.cos_dec
        if self.selection_shape == 'rectangle':
            self.selected_region = {
                'type': 'rectangle',
                'ra_min': min(x0, x1), 'ra_max': max(x0, x1),
                'dec_min': min(y0, y1), 'dec_max': max(y0, y1),
            }
        else:
            r_scaled = np.sqrt((x1_s - x0_s)**2 + (y1 - y0)**2)
            self.selected_region = {
                'type': 'circle',
                'ra_center': x0, 'dec_center': y0,
                'radius_scaled': r_scaled,
            }
        self.selection_start = None
        self._update_selection_display()

    def _get_fiber_mask(self) -> np.ndarray:
        """Return a boolean mask of fibers inside ``self.selected_region``.

        Returns
        -------
        `numpy.ndarray`_
            1D boolean array of length ``self.nspec``, ``True`` for fibers
            whose sky position falls within the currently drawn selection
            region.  Returns all-``False`` if no region is defined.
        """
        mask = np.zeros(self.nspec, dtype=bool)
        r = self.selected_region
        if r is None:
            return mask
        if r['type'] == 'rectangle':
            mask = ((self.ra >= r['ra_min']) & (self.ra <= r['ra_max'])
                    & (self.dec >= r['dec_min']) & (self.dec <= r['dec_max']))
        else:  # circle (radius is in scaled coordinates)
            dra = self.ra_scaled - (r['ra_center'] * self.cos_dec)
            ddec = self.dec - r['dec_center']
            mask = (dra**2 + ddec**2) <= r['radius_scaled']**2
        return mask

    def _update_selection_display(self) -> None:
        """Redraw hexagon edges to highlight selected fibers.

        Fibers selected by the drawn region (``_get_fiber_mask``) **or** by
        individual clicks (``self.clicked_fibers``) are drawn with a thick
        cyan edge; all others revert to a thin gray edge.
        """
        region = self._get_fiber_mask()
        edges = []
        widths = []
        for i in range(self.nspec):
            if region[i] or i in self.clicked_fibers:
                edges.append('cyan')
                widths.append(3.0)
            else:
                edges.append('gray')
                widths.append(0.5)
        self.hex_collection.set_edgecolors(edges)
        self.hex_collection.set_linewidths(widths)
        self.fig.canvas.draw()
