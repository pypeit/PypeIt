
.. include:: ../include/links.rst

.. _build_archived_sensfuncs:

How to build archived sensitivity functions
===========================================

1. Make sure you have up to date versions of the `PypeIt
   <https://github.com/pypeit/PypeIt>`_ and `PypeIt Development Suite`_ 
   repositories from GitHub, and that PypeIt :ref:`installing` is complete.

2. Download the reduced slitless source images used to create Keck DEIMOS sensitivity functions from
   the `PypeIt dev-suite Google Drive`_. The files are located under ``DEIMOS_Dev/Throughput/throughput_gdw/data_products/``.

3. Run the ``create_deimos_sensfuncs.py`` script. The below examples assumes the files
   were downloaded to a ``data_products`` directory and will create the archival 
   sensitivity function files under ``sensfunc_files``.

   .. code-block:: bash

        mkdir sensfunc_files
        PypeIt-development-suite/sensfunc_archive/create_deimos_sensfuncs.py all data_products sensfunc_files

4. Copy the results to the PypeIt repository.

   .. code-block:: bash

        cd sensfunc_files
        cp keck_deimos_*.fits PypeIt/pypeit/data/sensfuncs/

How to build the scripts that generate the archived sensitivity functions
=========================================================================

The archived sensitivity functions for DEIMOS were created from raw data
gathered from scripted standard star observations. These are slitless
observations that PypeIt cannot reduce on its own, so they are reduced using
IDL scripts from Greg Wirth. 

Examples of this data can be found in the `PypeIt dev-suite Google Drive`_ Google Drive.

+----------------------------------------------------------------+-------------------------------------------------------------------------------+
| ``DEIMOS_Dev/Throughput/YYYY``                                 | Contains raw data files from that year.                                       |
+----------------------------------------------------------------+-------------------------------------------------------------------------------+
| ``DEIMOS_Dev/Throughput/throughput_gdw/data_products/extract`` | Contains the reduced data from IDL used to generate the sensitivity functions.|
+----------------------------------------------------------------+-------------------------------------------------------------------------------+

Creating the IDL reduced output
-------------------------------
TODO: Does brad want to document this?


Converting the IDL reduced files to spec1d
------------------------------------------
In the PypeIt-development-suite repo,  the ``batch_convert_throughput_to_spec1d.py``
script can create PypeIt spec1d files from the IDL output.

For example, if the files from Google Drive was downloaded to a
``data_products`` directory. The below would create a ``spec1d_files`` directory with
spec1d files created from the IDL output:

   .. code-block:: bash

      python $PYPEIT_DEV/dev_algorithms/fluxing/batch_convert_throughput_to_spec1d.py data_products/extract spec1d_files pad

Creating the individual sensitivity functions
---------------------------------------------
The sensitivity functions can be crated with a pypeit_sensfunc call. For DEIMOS we
used the following input file:

   .. code-block:: ini

      [sensfunc]
         multi_spec_det = 3,7
         algorithm = "IR"
         extrap_blu = 0.
         extrap_red = 0.


The sensitivity functions can be created in bulk using the ``run_sensfunc_on_all_spec1d.py`` script.
Assuming the sensfunc file above is in a file named ``deimos_arxiv_sensfunc``, the below
command will create individual sensfunc files for each matching spec1d file and place
them in the sens_files directory:

   .. code-block:: bash

      python $PYPEIT_DEV/dev_algorithms/fluxing/run_sensfunc_on_all_spec1d.py spec1d_files/600ZD/ *2023jan17*.fits sens_files/600ZD deimos_arxiv_sensfunc

Combining multiple sensivity functions
======================================

The archived DEIMOS sensitivity functions are created by combining 
multiple sensivity functions from different detectors at different 
grating tilts to get a wider wavelength coverage. The steps in this
process are typically:

   1 Choose sensitivity functions for grating. The functions should have good wavelength coverage,
     and low airmass.

   2 Start with the bluest sensitivity function and translate the next function up/down to roughly match.
     One or both function may be truncated to fit.

   3 If the gap between sensitivity function isn't smooth (for example on a detector boundary) 
     a polynomial fit can be used to join them.

Helper functions in TODO/stitcutils.py. See jupyter notebook todo for a detailed example of how the sesnitivity
functions for the 600ZD DEIMOS grating were stitched together.

