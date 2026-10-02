"""
Matplotlib GUI for interactively selecting fibers from a fiber-fed IFU
spec1d file and extracting a combined 1D spectrum.

.. include:: ../include/links.rst
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle, RegularPolygon, Rectangle
from matplotlib.ticker import FuncFormatter
from matplotlib.widgets import Button, CheckButtons, RadioButtons, RangeSlider

from pypeit import log, PypeItError
from pypeit.core import datacube


# Bright optical sky lines (vacuum, Angstrom).  Used only for the integrated
# per-fiber flux shown in the fiber-IFU extractor display, so that residual
# sky-subtraction errors do not bias the per-fiber flux estimate.  Extraction
# itself uses the full wavelength range.
#
# TODO: Replace this hand-rolled list with PypeIt's built-in sky-line
#       resources (e.g. ``data/sky_spec/sky_single_mg.dat`` or the
#       ``OH_*_lines.dat`` line lists; cf. the matching TODO in pypeit/qa.py),
#       and expose the mask half-width as a GUI control rather than a module
#       constant.
_SKY_LINE_WAVELENGTHS = (
    5577.34,                       # [OI]
    5890.00, 5896.00,              # Na I D
    6300.30, 6363.78,              # [OI]
    7340.0, 7821.0, 8344.0, 8430.0, 8761.0, 8867.0,  # OH band centers (approx)
)
_SKY_LINE_MASK_HALFWIDTH = 10.0    # Angstrom


def sky_line_mask(wave):
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


def compute_fiber_fluxes(waves, fluxes, wave_min, wave_max):
    """Per-fiber integrated flux with sky-line masking.

    Each fiber is integrated over its **own** native wavelength grid, with
    pixels in any sky-line region (see ``_SKY_LINE_WAVELENGTHS``) excluded.
    Used to colour the fiber-IFU extractor display; the actual extracted
    spectrum (:func:`pypeit.core.datacube.resample_and_combine`) is
    unaffected.

    Parameters
    ----------
    waves, fluxes : :obj:`list` of `numpy.ndarray`_
        Per-fiber native wavelengths and fluxes (parallel lists).
    wave_min, wave_max : :obj:`float`
        Integration range, Angstrom.

    Returns
    -------
    `numpy.ndarray`_
        1D array of integrated fluxes, one per fiber.
    """
    out = np.zeros(len(waves), dtype=float)
    for i, (w, f) in enumerate(zip(waves, fluxes)):
        in_range = (w >= wave_min) & (w <= wave_max)
        sky = sky_line_mask(w)
        mask = in_range & ~sky
        out[i] = np.nansum(f[mask])
    return out


class BinospecIFUExtractGUI:
    """Matplotlib widget GUI for selecting fibers and extracting spectra.

    Displays the IFU as a hexagonal flux map, lets the user select fibers by
    drawing a region or clicking, and writes the combined selection as a
    :class:`~pypeit.onespec.OneSpec` FITS file.

    Parameters
    ----------
    fibers : :obj:`list` of :obj:`dict`
        Per-fiber records (see :func:`pypeit.core.datacube.load_fibers`), each
        carrying ``wave``, ``flux``, ``ivar``, ``ra``, and ``dec``.
    raw_hdr : `astropy.io.fits.Header`_
        Primary header of the source spec1d, copied to the output.
    pyp_spec : :obj:`str`
        Spectrograph short name, written as ``PYP_SPEC`` on the output.
    outfile : :obj:`str`
        Default output OneSpec FITS path.
    """

    def __init__(self, fibers, raw_hdr, pyp_spec, outfile):
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

        self.fiber_fluxes = compute_fiber_fluxes(
            self.waves, self.fluxes, self.wave_min, self.wave_max)

        # Selection state
        self.y_scale = 'linear'
        self.mask_sky = False
        self.selection_shape = 'rectangle'
        self.selection_start = None
        self.selection_patch = None
        self.selected_region = None
        self.clicked_fibers = set()
        self.fiber_mask = None
        self.extracted_wave = None
        self.extracted_flux = None
        self.extracted_ivar = None

        self._create_gui()

    def show(self):
        """Display the interactive GUI."""
        plt.show()

    def _create_gui(self):
        """Build all GUI axes, widgets, and event connections."""
        self.fig = plt.figure(figsize=(16, 8))
        self.ax_image = plt.subplot2grid((4, 3), (0, 0), rowspan=4,
                                         colspan=2)
        self.ax_spectrum = plt.subplot2grid((4, 3), (0, 2), rowspan=2,
                                            colspan=1)
        self.ax_controls = plt.subplot2grid((4, 3), (2, 2), colspan=1)
        self.ax_coords = plt.subplot2grid((4, 3), (3, 2), colspan=1)
        self.ax_controls.axis('off')
        self.ax_coords.axis('off')

        self._create_hexagon_plot()

        self.ax_spectrum.set_xlabel('Wavelength (Å)')
        self.ax_spectrum.set_ylabel('Flux')
        self.ax_spectrum.set_title('Extracted Spectrum')
        self.ax_spectrum.grid(True, alpha=0.3)

        # Buttons
        self.btn_extract = Button(
            plt.axes([0.77, 0.34, 0.15, 0.05]), 'Extract Spectrum')
        self.btn_extract.on_clicked(self._on_extract)
        self.btn_reset = Button(
            plt.axes([0.77, 0.28, 0.15, 0.05]), 'Reset Display')
        self.btn_reset.on_clicked(self._on_reset)
        self.btn_save = Button(
            plt.axes([0.77, 0.22, 0.15, 0.05]), 'Save Spectrum')
        self.btn_save.on_clicked(self._on_save)

        # Mode radio
        self.radio = RadioButtons(
            plt.axes([0.68, 0.22, 0.08, 0.10]),
            ('Rectangle', 'Circle', 'Click Fibers'))
        self.radio.on_clicked(self._on_shape_change)

        # Y-axis scale radio for the extracted-spectrum panel
        self.yscale_radio = RadioButtons(
            plt.axes([0.77, 0.13, 0.08, 0.08]),
            ('Linear', 'Log'))
        self.yscale_radio.on_clicked(self._on_yscale_change)

        # Sky-line masking toggle for the extracted-spectrum panel
        self.skymask_check = CheckButtons(
            plt.axes([0.86, 0.13, 0.12, 0.08]),
            ['Mask Sky Lines'], [self.mask_sky])
        self.skymask_check.on_clicked(self._on_skymask_change)

        # Wavelength range slider
        self.wave_slider = RangeSlider(
            plt.axes([0.1, 0.05, 0.45, 0.03]), 'Wavelength (Å)',
            self.wave_min, self.wave_max,
            valinit=(self.wave_min, self.wave_max),
            valstep=((self.wave_max - self.wave_min) / 100) or None)
        self.wave_slider.on_changed(self._on_wave_change)

        # Mouse events
        self.fig.canvas.mpl_connect('button_press_event', self._on_press)
        self.fig.canvas.mpl_connect('motion_notify_event', self._on_motion)
        self.fig.canvas.mpl_connect('button_release_event', self._on_release)

        plt.tight_layout()
        plt.subplots_adjust(bottom=0.10)

    def _create_hexagon_plot(self):
        """Populate ``ax_image`` with a hexagonal fiber array flux map."""
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

        patches = [RegularPolygon((self.ra_scaled[i], self.dec[i]),
                                  numVertices=6, radius=hex_radius,
                                  orientation=np.pi / 2)
                   for i in range(self.nspec)]
        self.hex_collection = PatchCollection(patches, cmap='viridis',
                                              edgecolors='gray',
                                              linewidths=0.5)
        self.hex_collection.set_array(self.fiber_fluxes)
        self.hex_collection.set_clim(vmin, vmax)
        self.ax_image.add_collection(self.hex_collection)

        # nspec >= 1 is guaranteed by load_fibers raising PypeItError on
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
        plt.colorbar(self.hex_collection, ax=self.ax_image,
                     label='Integrated Flux')

    # ------------------------------------------------------------------
    # Callbacks
    # ------------------------------------------------------------------

    def _on_wave_change(self, val):
        wave_min, wave_max = val
        self.current_wave_min = float(wave_min)
        self.current_wave_max = float(wave_max)
        self.fiber_fluxes = compute_fiber_fluxes(
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

    def _apply_spectrum_xlim(self):
        """Sync the spectrum axis to the current slider window."""
        if self.extracted_wave is None:
            return
        lo, hi = self.current_wave_min, self.current_wave_max
        self.ax_spectrum.set_xlim(lo, hi)
        in_window = ((self.extracted_wave >= lo)
                     & (self.extracted_wave <= hi))
        if self.mask_sky:
            in_window &= ~sky_line_mask(self.extracted_wave)
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

    def _selected_indices(self):
        mask = self._get_fiber_mask()
        for i in self.clicked_fibers:
            mask[i] = True
        return np.where(mask)[0]

    def _on_extract(self, event):
        idxs = self._selected_indices()
        if idxs.size == 0:
            log.info("No fibers selected")
            return
        log.info(f"Extracting from {idxs.size} fibers")
        try:
            wave, flux, ivar = datacube.resample_and_combine(
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

    def _draw_spectrum(self):
        """Redraw the extracted spectrum honoring the current display toggles.

        On a log scale, non-positive flux values (e.g. sky-subtraction
        residuals) are masked so they neither appear in the trace nor
        distort the axis limits.  When sky-line masking is enabled, pixels
        within ``_SKY_LINE_MASK_HALFWIDTH`` of a line in
        ``_SKY_LINE_WAVELENGTHS`` are likewise masked.
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

        if self.mask_sky:
            # Blank out sky-line regions so residuals there don't clutter
            # the displayed trace.  The saved spectrum is unaffected.
            sky = sky_line_mask(wave)
            plot_flux = np.where(sky, np.nan, plot_flux)
            lower = np.where(sky, np.nan, lower)
            upper = np.where(sky, np.nan, upper)

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

    def _on_yscale_change(self, label):
        self.y_scale = label.lower()
        if self.extracted_wave is None:
            self.ax_spectrum.set_yscale(self.y_scale)
        else:
            self._draw_spectrum()
        self.fig.canvas.draw()

    def _on_skymask_change(self, label):
        self.mask_sky = bool(self.skymask_check.get_status()[0])
        if self.extracted_wave is not None:
            self._draw_spectrum()
            self.fig.canvas.draw()

    def _on_save(self, event):
        if self.extracted_wave is None:
            log.info("Nothing to save -- run Extract Spectrum first")
            return
        outfile = self._default_save_path()
        datacube.write_onespec(self.extracted_wave, self.extracted_flux,
                               self.extracted_ivar, self.raw_hdr,
                               self.pyp_spec, outfile)
        log.info(f"Wrote {outfile}")

    def _default_save_path(self):
        """Return a non-clobbering default output path.

        Starts from ``self.outfile`` (``..._<base>.fits``) and inserts a
        two-digit counter before the extension (``..._<base>_00.fits``),
        incrementing the counter until a path that does not already exist on
        disk is found.

        Returns
        -------
        :obj:`str`
            The first ``..._<NN>.fits`` path (NN starting at ``00``) that does
            not yet exist.
        """
        path = Path(self.outfile)
        n = 0
        while True:
            candidate = path.with_name(f"{path.stem}_{n:02d}{path.suffix}")
            if not candidate.exists():
                return str(candidate)
            n += 1

    def _on_reset(self, event):
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
        self.fiber_fluxes = compute_fiber_fluxes(
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

    def _on_shape_change(self, label):
        self.selection_shape = label.lower()
        if self.selection_patch:
            self.selection_patch.remove()
            self.selection_patch = None
        if self.selection_shape == 'click fibers':
            self.ax_image.set_title('Fiber Array (click hexagons to add/remove)')
        else:
            self.ax_image.set_title('Fiber Array (drag to select region)')
        self.fig.canvas.draw()

    def _on_press(self, event):
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

    def _on_motion(self, event):
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
            self.selection_patch = Rectangle(
                (x0, y0), x1 - x0, y1 - y0,
                fill=False, edgecolor='cyan', linewidth=2, linestyle='--')
        else:
            r = np.sqrt((x1 - x0)**2 + (y1 - y0)**2)
            self.selection_patch = Circle(
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

    def _on_release(self, event):
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

    def _get_fiber_mask(self):
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

    def _update_selection_display(self):
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
