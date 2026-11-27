VISTA/4MOST LRS-B red
====================

The experimental ``vista_4most_lrsb_red`` adapter uses the ``Fiber``
pipeline for the unbinned red detector of the second low-resolution
spectrograph (LRS-B). Setup identifies science, bias, daytime flat, and
simultaneous ThAr files from their ESO headers. The adapter maps the slit
positions to the science and sky fibers from the exposure's ``FIBINFO``
extension.

The commissioning example uses an external trace NPZ, selected through
``PYPEIT_4MOST_TRACE``, and imports a ThAr/Fabry–Pérot wavelength
solution and relative twilight throughput into the reduction's
calibration directory. PypeIt processes the raw science detector image,
models the sky from the assigned sky fibers, and writes boxcar and
optimal spectra. The optimal extraction uses spatial profiles from a
processed daytime lamp flat.

The current comparison covers ESO OB 5148227. Its median PypeIt to joint
continuum ratio is 0.950 in the 844–886 nm band. The offset persists
when both extractors use the same detector pixels, sky model, and
throughput factors. The commissioning run uses relative electron
counts and leaves the spectra in the observed wavelength frame.
