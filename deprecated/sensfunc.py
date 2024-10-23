def compute_blaze(self, wave, trace_spec, trace_spat, flatfile, box_radius=10.0,
                  min_blaze_value=1e-3, debug=False):
    """
    Compute the blaze function from a flat field image.

    Args:
        wave (`numpy.ndarray`_):
            Wavelength array. Shape = (nspec, norddet)
        trace_spec (`numpy.ndarray`_):
            Spectral pixels for the trace of the spectrum. Shape = (nspec, norddet)
        trace_spat (`numpy.ndarray`_):
            Spatial pixels for the trace of the spectrum. Shape = (nspec, norddet)
        flatfile (:obj:`str`):
            Filename for the flat field calibration image
        box_radius (:obj:`float`, optional):
            Radius of the boxcar extraction region used to extract the blaze function in pixels
        min_blaze_value (:obj:`float`, optional):
            Minimum value of the blaze function. Values below this are clipped and set to this value. Default=1e-3
        debug (:obj:`bool`, optional):
            Show plots useful for debugging. Default=False

    Returns:
        `numpy.ndarray`_: The log10 blaze function. Shape = (nspec, norddet)
        if norddet > 1, else shape = (nspec,)
    """
    flatImages = flatfield.FlatImages.from_file(flatfile, chk_version=self.chk_version)

    pixelflat_raw = flatImages.pixelflat_raw
    pixelflat_norm = flatImages.pixelflat_norm
    pixelflat_proc, flat_bpm = flat.flatfield(pixelflat_raw, pixelflat_norm)

    flux_box = moment1d(pixelflat_proc * np.logical_not(flat_bpm), trace_spat, 2 * box_radius, row=trace_spec)[0]

    pixtot = moment1d(pixelflat_proc * 0 + 1.0, trace_spat, 2 * box_radius, row=trace_spec)[0]
    pixmsk = moment1d(flat_bpm, trace_spat, 2 * box_radius, row=trace_spec)[0]

    mask_box = (pixmsk != pixtot) & np.isfinite(wave) & (wave > 0.0)

    # TODO This is ugly and redundant with spec_atleast_2d, but the order of operations compels me to do it this way
    blaze_function = (np.clip(flux_box * mask_box, 1e-3, 1e9)).reshape(-1, 1) \
        if flux_box.ndim == 1 else flux_box * mask_box
    wave_debug = wave.reshape(-1, 1) if wave.ndim == 1 else wave
    log10_blaze_function = np.zeros_like(blaze_function)
    norddet = log10_blaze_function.shape[1]
    for iorddet in range(norddet):
        blaze_function_smooth = utils.fast_running_median(blaze_function[:, iorddet], 5)
        blaze_function_norm = blaze_function_smooth / blaze_function_smooth.max()
        log10_blaze_function[:, iorddet] = np.log10(np.clip(blaze_function_norm, min_blaze_value, None))
        if debug:
            plt.plot(wave_debug[:, iorddet], log10_blaze_function[:, iorddet])
    if debug:
        plt.show()

    # TODO It would probably better to just return an array of shape (nspec, norddet) even if norddet = 1, i.e.
    # to get rid of this .squeeze()
    return log10_blaze_function.squeeze()
