
.. include:: ../include/links.rst

.. _construct_template:

======================================
Constructing a New Wavelength Template
======================================


Overview
========

This doc describes how to create a wavelength template used to perform
wavelength calibration with the :ref:`wvcalib-algorithms`.

See :doc:`wave_calib` for a discussion of various algorithms and
see :doc:`wvcalib` for a discussion of the
main outputs and good/bad examples.

Finding Data
============

To build a new template, the requirements are wavelength calibrations for the
arc spectra over the full wavelength range needed. This may take multiple
exposures of arc lamps to cover the full range.

Using the Holy Grail Algorithm
==============================

Once you have a set of different arc lamp exposures that cover the desired range,
run through the data using the option of the :ref:`wvcalib-holygrail` algorithm to generate a
wavelength calibration and select the slits where the technique is
successful. See :doc:`wvcalib` for an example QA plot that illustrates
a good wavelength solution.

With a complete set of correctly calibrated arc line data, identify wavelength
regions, with no gaps, that cover the full range in wavelength desired.

In situations where the Holy Grail fails, one may resort to :ref:`pypeit_identify`.

If a calibrated arc line spectrum that covers the desired wavelength range is
already available, an ascii file with 2 columns, ``wavelength`` and ``flux``,
can be used instead of running the reduction with the :ref:`wvcalib-holygrail`
algorithm.

Creating the Template
=====================

Once all of the slits and wavelength calibration files are selected, use
:func:`~pypeit.core.wavecal.templates.build_template` to turn the
list of slits, wavelength ranges, and :doc:`wvcalib` files into a final
template.

For longevity, it is a good idea to install the `PypeIt Development Suite`_ (see
:ref:`here<dev-suite>`) and combine your solution with all PypeIt wavelength
solutions using the same instrument and grating. If you used the
:ref:`pypeit_identify` utility, it will output a file called ``wvcalib.fits``.
You should put this file in the appropriate instrument folder in the PypeIt
Development Suite templates directory:

.. code-block:: console

    $PYPEIT_DEV/dev_algorithms/wavelengths/template_files/

you should also rename this file to match the current formatting:

.. code-block:: console

    mv wvcalib.fits ${TELESCOPE}_${INSTRUMENT}_${GRATING}_${WAVECEN}.fits

Now, open the
``pypeit/core/wavecal/spectrographs/templ_TELESCOPE_INSTRUMENT.py`` file, and
either (1) add a new setup if yours doesn't already exist or (2) edit a setup if
you want to add your solution to a pre-existing setup.  The latter is only
usually warranted if you are extending the wavelength coverage of a previous
setup. The key information you will need is the spectral binning and the slit
``spatid``.  Both of these numbers are printed on the command line when you
complete the fitting with the :ref:`pypeit_identify` utility.  Once this is
setup, simply execute ``python templates.py`` from the directory that contains
``templ_TELESCOPE_INSTRUMENT.py``.

This will automatically put your solution into the ``reid_arxiv`` directory.
Add all of these files to git, and submit a PR for both PypeIt and the PypeIt
Development Suite.

Alternatively, if you want to construct a template file yourself (with the
possible disadvantage that your solution cannot be stitched together with
solutions from the same spectrograph+grating in the future), you can write a
script where the input looks like the following:

.. code-block:: python

    from pypeit.core.wavecal import templates

    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot, ifiles=ifiles,
                             det_cut=det_cut, chk=True, normalize=False, lowredux=False,
                             subtract_conti=True, overwrite=overwrite, shift_wave=True)

See :func:`~pypeit.core.wavecal.templates.build_template` for a description of
all the parameters, and ``pypeit/core/wavecal/spectrographs`` for examples of
the use of this function.

This produces a file called ``outroot`` that contains the template. The templates are saved in
``pypeit/data/arc_lines/reid_arxiv``. It also produces a plot of the final product.

