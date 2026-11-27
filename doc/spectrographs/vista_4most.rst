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

Commissioning aperture calibration
----------------------------------

An optional grey normalization anchors the optimally extracted counts to
the boxcar continuum and preserves the optimal spectral shape. Set
``PYPEIT_4MOST_APERTURE_SCALE`` to the positive wide-aperture-to-boxcar
ratio measured from an independent training exposure. The commissioning
example used ``1.013349531169947``, the median ratio for five stars in
ESO OB 5148227 with a 4.75-pixel aperture.

The adapter measures the boxcar-to-optimal ratio from common valid rows
in 8440--8860 Angstrom in vacuum. It applies the normalization to
science fibers with at least 900 median electrons in a two-pixel-radius
detector aperture and no neighboring trace within eight pixels whose
core exceeds 10 percent of the target's core. The boxcar extraction and
spectral masks retain their original values. The primary FITS header
records the training scale as ``APERSCL`` and the reduction log lists
the effective factor for each corrected fiber.

The normalization has been tested in 8440--8860 Angstrom on bright,
isolated commissioning targets. On 24 stars in the held-out OB 5148408,
the median continuum residual against the wide aperture went from
``-4.35%`` to ``-0.04%``. The propagated count uncertainties are conditional
on the normalization. They do not include a shared normalization
uncertainty across wavelengths.
