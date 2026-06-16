.. include:: ../include/links.rst

.. _focas:

**************
Subaru FOCAS
**************

Overview
========

This file summarizes several instrument specific
settings that are related to the Subaru/FOCAS spectrograph.

.. warning::

    PypeIt currently *cannot* reduce images produced by FOCAS
    in imaging mode or those taken with certain obsolete readout
    configurations. All image-handling assumes FOCAS images have been
    taken in spectroscopic mode with standard CCD readout parameters.
    PypeIt handles files that do not meet these criteria in two ways:

        - When running :ref:`pypeit_setup`, any frames not in
          spectroscopic mode or with non-standard readout will be ignored
          and should not appear in your :ref:`pypeit_file`.

        - If you add frames to the :ref:`pypeit_file` that are not in
          spectroscopic mode or have incompatible readout configurations,
          the method used to read the FOCAS files will fault.

At present,the code has only been testd on DET-ID=2 data, 
i.e. the second of the two detectors, which is where the 
primary target lands if one uses the longslit.

Deviations
==========

The default changes to the PypeIt parameters specific to FOCAS data are listed
here: :ref:`instr_par-subaru_focas`.  *You do not have to add these changes to
your PypeIt reduction file!*  This is just a listing of how the parameters used
for Subaru/FOCAS differ from the defaults listed in the preceding tables on
that page.

These are tuned to the standard calibration
set taken with FOCAS.

.. _focas_detector:

DETECTOR
========

FOCAS uses two 2048×4096 MIT/LL CCD detectors. Unlike other multi-detector instruments,
FOCAS does not use mosaic construction. The detector is oriented such that the
dispersion direction is along the 4096-pixel axis (y-direction) and the spatial
direction is along the 2048-pixel axis (x-direction).

The detector configuration is automatically handled by PypeIt with:

.. code-block:: ini

    [rdx]
        spectrograph = subaru_focas

This is the default configuration and typically does not need to be modified
by the user. 

Calibrations
============

Edge Tracing
------------

FOCAS long-slit and multi-object spectroscopy requires careful slit edge tracing.
The default `edge_thresh` parameter works well for most FOCAS configurations,
but may need adjustment for:

- Very short exposures
- High-resolution modes (VPH gratings)

If slit edges are not being detected properly, try:

.. code-block:: ini

    [calibrations]
        [[slitedges]]
            edge_thresh = 20

For multi-object spectroscopy (MOS), you may need to adjust additional parameters:

.. code-block:: ini

    [calibrations]
        [[slitedges]]
            edge_thresh = 20
            trace_thresh = 10
            minimum_slit_length = 6

Multi-Object Spectroscopy (MOS)
-------------------------------

FOCAS supports multi-object spectroscopy through slit masks. Unlike DEIMOS,
FOCAS does not currently support automatic matching to slit-mask design files.
Users must rely on the automatic slit detection and may need to manually
verify slit assignments.

For MOS observations, consider these parameters:

.. code-block:: ini

    [calibrations]
        [[slitedges]]
            sync_predict = nearest
            bound_detector = True
            minimum_slit_length = 8

Wavelength Calibration
----------------------

FOCAS wavelength calibration depends on the grating and comparison lamp used.
Common lamp combinations include:

- **HgCd lamps**: Standard for most gratings
- **Ne+Ar lamps**: Used for some high-resolution configurations
- **ThAr lamps**: Available for high-precision radial velocity work

The default lamp selection works for most cases, but can be overridden:

.. code-block:: ini

    [calibrations]
        [[wavelengths]]
            lamps = HgI, CdI


**Common Grating Configurations:**

- **150B, 300B, 600B**: Blue gratings, use HgCd lamps
- **150R, 300R, 600R**: Red gratings, use HgCd or Ne+Ar lamps  
- **VPH gratings**: High-resolution, may require ThAr lamps

Flat Fielding
-------------

FOCAS flat fielding generally works well with default parameters.
However, for certain configurations you may encounter:

- **Fringing in red wavelengths** (>7000Å): Common with red gratings
- **Illumination gradients**: Particularly in wide-slit modes
- **Scattered light**: Can affect faint object spectroscopy

To address fringing in red spectra:

.. code-block:: ini

    [calibrations]
        [[flatfield]]
            slit_illum_finecorr = False
            tweak_slits_thresh = 0.9

For severe illumination gradients:

.. code-block:: ini

    [calibrations]
        [[flatfield]]
            slit_illum_finecorr = True
            slit_illum_ref_idx = 0

Flexure
-------

FOCAS generally has good mechanical stability, but flexure corrections
may be needed for:

- **Large telescope movements** between calibrations and science
- **High-precision radial velocity work**

Enable flexure correction with:

.. code-block:: ini

    [flexure]
        spec_method = boxcar
        spec_maxshift = 20

For high-precision work, consider using sky lines for flexure correction:

.. code-block:: ini

    [flexure]
        spec_method = slitcen
        spec_maxshift = 10



Flux Calibration
================

FOCAS flux calibration requires standard star observations.
Recommended standards include:

- **Spectrophotometric standards**: BD+28d4211, HD 19445, etc.
- **White dwarf standards**: G191-B2B, GD 153, etc.

Standard star configuration:

.. code-block:: ini

    [fluxcalib]
        extinct_correct = True
        extinct_file = atm_extinct_subaru_maunakea.dat


Troubleshooting
===============

**Common Issues:**

1. **Wavelength solution fails**
   
   - Check lamp selection matches your grating
   - Verify arc lamp exposure times are adequate
   - Consider manual lamp line identification

2. **Poor flat fielding**
   
   - Check for saturated flat field exposures
   - Verify flat lamp matches science configuration

3. **Slit edges not detected**
   
   - Lower edge_thresh parameter
   - Verify adequate flat field signal


Additional Reading
==================

Here are additional docs related to Subaru/FOCAS:

.. toctree::
   :maxdepth: 1

   ../dev/add_missing_obj
   ../dev/fluxing