.. include:: ../include/links.rst

.. _apf_howto:

=================
APF-Levy HOWTO
=================

Overview
========

This document provides a guide to using the APF-Levy spectrograph with PypeIt.
It covers the specific configurations, calibration procedures, and data
reduction steps necessary.

The final extractions cover 55 orders in total, with the first order, physical
order 122, located at 3800 Angstroms. The reddest order, 68, is located at 6850
Angstroms.

Setup
=====

Organize the data
-----------------

PypeIt reduces all 3" and 8" long slits separately. The reductions share,
however, the WideFlat frames which are acquired using the 8" long slit with a 2"
width. These files are used to estimate the pixel flat images. Therefore, there
are additional, manual steps that the user must take. 

It is best to separate science and ThAr frames taken with the 3" deckers from
data taken with the 8" deckers, this is not required, but it does make the
reduction process easier. NarrowFlat images, flat field images taken with the 3"
long slit, are used only for reducing the 3" long slit science data.

Run pypeit_setup
----------------

After running :ref:`pypeit_setup` with a specific configuration - the ``-c A``
flag, there will be a pypeit file created in the appropriate directory for that
configuration. This file will have the name of the configuration, e.g., 
``apf_levy_A.pypeit``. If you mix your 3" and 8" data, you will also have
``apf_levy_B.pypeit`` if you type ``-c B``.  The number of configurations can be
determined by reading the file ``setup_files/apf_levy.sorted``. 

.. note::

    If you mix the data, you will get an error when running pypeit_setup.  

For 8" data, the ``WideFlat`` images should be used for ``trace`` frames as well
as ``pixflat`` frames, no editing is required.

For the 3" data, you need to edit the pypeit file. By default the WideFlat
images will be listed as ``pixflat,trace`` and the labeled must be updated to
just ``pixflat``. As there are 100 WideFlat images taken per night, editing by
hand can be time consuming. Use a search and replace function or a script to
automate the process. Below are two examples of how to do this:

.. code-block:: bash

    sed -i 's/pixflat,trace/pixflat/' apf_levy_A.pypeit
    perl -pi -e 's/pixflat,trace/pixflat/' apf_levy_A.pypeit

For the 3" data, use the ``NarrowFlat`` images for ``trace`` frames. These will
be identified automatically.

Non-standard binning
--------------------

The pipeline has not been tested for non-standard (non 1x1) 
binning. Please contact the developers if there are issues.


Main Run
========

Once the :doc:`../pypeit_file` is ready, the main call is
simply:

.. code-block:: bash

    run_pypeit apf_levy_A.pypeit -o

The ``-o`` indicates that any existing output files should be overwritten.  As
there are none, it is superfluous but we recommend (almost) always using it.

The :doc:`../running` doc describes the process in some
more detail.

This will take a while, as there are 55 orders per exposure which will be
traced, wavelength calibrated, flat fielded, and extracted. 

By default the 3" data has no local sky subtraction, but this can be enabled.
The 8" data has local sky subtraction enabled by default. To change this for
the 3" data add to the PypeIt file:

.. code-block:: ini

    [reduction]
        [[skysub]]
            no_local_sky = False

The other setting that is unusual for the 3" slit is that the optimal
extraction always uses a Gaussian profile instead of modeling the profile
with a B-spline. To restore normal behavior, add to the PypeIt file:

.. code-block:: ini

    [reduction]
        [[extraction]]
            sn_gauss = 4

The narrow length slit means that the model covers the full lenth of the slit,
so local sky subtraction is a challenge and the B-spline model often includes
the sky background in the fit.

Finally, the default box radius is 4 pixels, or 1.728 arc-seconds. This can be
reset using the boxcar_radius parameter in the extraction section of the PypeIt file:

.. code-block:: ini

    [reduction]
        [[extraction]]
            boxcar_radius = 1.296 # 3 pixels

Inspecting Files
================

Calibrations
------------

The first set are :doc:`../calibrations/calibrations`.
The APF Levy spectrograph has no moving parts that will
will change where the spectra land on the detector. 
This fixed format means that the trace and wavelength
calibrations are stable. 

Slit Edges
++++++++++

PypeIt will map the slit edges using the trace frames.
The orders locations of the orders to be extracted are
pre-defined in the spectrograph file inside PypeIt. This
means that, even if the order is not detected in the trace
frames, the order will still be extracted. 

Wavelengths
+++++++++++

One should inspect the :doc:`../qa` for the wavelength
calibration.  These are PNGs in the QA/PNG/ folder.

Note:  there are multiple files generated for every slit.
When the reduction is complete, you may prefer to scan
through them by opening the HTML file under ``QA/``.

The final wavelength solution is a two dimensional fit 
with pixel along the order one axis and the order itself
being the second. 

Remember, the default calibration is in vacuum wavelengths.

Spectra
-------

The code will generate 2D and 1D spectra outputs.  One per 
science frame, located in the ``Science/`` folder.

One can inspect the one dimensional spectra 
with :ref:`pypeit_show_1dspec`.

Using that tool, you can examine the individual orders. 
