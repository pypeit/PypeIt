
.. include:: ../include/links.rst

.. _build_archived_sensfuncs:


Creating archival sensitivity functions
=======================================

Archival sensitivity functions for instruments can be created by combining 
multiple sensivity functions from different exposures at different 
grating tilts to get a wider wavelength coverage. The steps in this
process are typically:

#. Choose sensitivity functions for grating. The functions should have good wavelength coverage,
   and low airmass.

#. Start with the bluest sensitivity function and translate the next function up/down to roughly match.
   One or both function may be truncated to fit.

#. If the gap between sensitivity function isn't smooth (for example on a detector boundary) 
   a polynomial fit can be used to join them.

There are helper functions in $PYPEIT_DEV/sensfunc_archive/stitcutils.py intended to make this easier. See
``$PYPEIT_DEV/sensfunc_archive/create_deimos_sensfuncs.py`` and ``$PYPEIT_DEV/sensfunc_archive/stitchdeimos.py`` for examples of how this was done
for DEIMOS.

The DEIMOS archived sensitivity functions
-----------------------------------------

*Disclaimer*:
The DEIMOS archival sensitivity functions do not provide an absolute flux 
calibration.  Instead, they are only intended to remove the instrumental 
response, providing a relative flux calibration up to some unknown 
normalization.


Creating DEIMOS stitched sensitivity functions for every grating
****************************************************************

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
        $PYPEIT_DEV/sensfunc_archive/create_deimos_sensfuncs.py all data_products sensfunc_files

   The above creates spec1d files and sensfunc files from the reduced slitless source images. It then
   stitches those sensfunc files together to form sensfuncs that cover a wider wavelength range.  This results
   in one stitched archival sensfunc for each grating (600ZD, 830G, 900ZD, 1200B, and 1200G).
   
   To stitch together already created sensfunc files use the ``--reuse`` option. The pre-stitch sensfunc files used to 
   create the archival files are in the Google Drive under ``DEIMOS_Dev/sensfunc_stitch/pre_stitch_sensfuncs/``

4. Copy the results to the PypeIt repository.

   .. code-block:: bash

        cd sensfunc_files
        cp keck_deimos_*.fits PypeIt/pypeit/data/sensfuncs/

How each DEIMOS sensitivity functions was built
***********************************************

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


Converting the IDL reduced files to spec1d
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the PypeIt-development-suite repo,  the ``batch_convert_throughput_to_spec1d.py``
script can create PypeIt spec1d files from the IDL output.

By default, this script splits the IDL output into two detectors (DET03 and DET07). Because the IDL output isn't always
even in length, it is neccessary to pass "pad" to the script. Alternatively the script can generate spec1d files by treating the IDL output as a 
mosaic ('MSC03') by passing "nosplit":

For example, if the files from Google Drive was downloaded to a
``data_products`` directory. The below would create a ``spec1d_files`` directory with
spec1d files created from the IDL output:

   .. code-block:: bash

      # Create two detector spec1d files
      python $PYPEIT_DEV/dev_algorithms/fluxing/batch_convert_throughput_to_spec1d.py data_products/extract spec1d_files pad

      # Create a mosaic spec1d file
      python $PYPEIT_DEV/dev_algorithms/fluxing/batch_convert_throughput_to_spec1d.py data_products/extract spec1d_files nosplit

The current DEIMOS archival sensfuncs were created by splitting the output into two detectors using the "pad" parameter.

Creating the individual sensitivity functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The individual DEIMOS sensitivity functions were created with the following ``pypeit_sensfunc`` input file (See :ref:`sensitivity_file`).

   .. code-block:: ini

      [sensfunc]
         multi_spec_det = 3,7
         algorithm = "IR"
         extrap_blu = 0.
         extrap_red = 0.


The sensitivity functions can be created in bulk using the ``run_sensfunc_on_all_spec1d.py`` script.
Assuming the sensfunc input file above is named ``deimos_arxiv_sensfunc.sens``, the below
example will create individual sensfunc files for each matching 600ZD spec1d file and place
them in the ``sens_files`` directory:

   .. code-block:: bash

      python $PYPEIT_DEV/dev_algorithms/fluxing/run_sensfunc_on_all_spec1d.py spec1d_files/600ZD/ *2023jan17*.fits sens_files/600ZD deimos_arxiv_sensfunc.sens

