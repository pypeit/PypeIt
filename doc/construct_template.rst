======================================
Constructing a New Wavelength Template
======================================

.. index:: construct_template

Overview
========

This doc describes how to make the wavelength calibration used in
`Automated Algorithms`_

See :doc:`wave_calib` for a discussion of various algorithms and
see :doc:`master_wvcalib` for a discussion of the
main outputs and good/bad examples.

Finding Data
============

To build a new template, the requirements are wavelength calibrations for the
arc spectra over the full wavelength range needed. This may take multiple
exposures of arc lamps to cover the full range.

Using the Holy Grail Algorithm
==============================

Once you have a set of different arc lamp exposures that cover the desired range,
run through the data using the option of the Holy Grail to generate a
master wavelength calibration and select the slits where the technique is
successful. See :doc:`master_wvcalib` for an example QA plot that illustrates
a good wavelength solution.

With a complete set of correctly calibrated arc line data, identify wavelength
regions, with no gaps, that cover the full range in wavelength desired.

In situations where the Holy Grail fails, one may resort to

    pypeit_identify

Creating the Template
=====================

Once all of the slits and master wavelength calibration files are
selected, use the recipe in 'templ_keck_deimos.py' to turn the list of slits,
wavelength ranges, and master wavecalib files into a final template.

The input looks like the following:

    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
     ifiles=ifiles, det_cut=det_cut, chk=True,
     normalize=False, lowredux=False,
     subtract_conti=True, overwrite=overwrite,
     shift_wave=True)

This produces a file called 'outroot' that contains the template. It also
produces a plot of the final product.
