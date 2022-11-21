
.. include:: ../include/links.rst

.. _build_archived_sensfuncs:

How to build archived sensitivity functions
===========================================

1. Make sure you have up to date versions of the `PypeIt
   <https://github.com/pypeit/PypeIt>`_ and `PypeIt Development Suite`_ 
   repositories from GitHub, and that PypeIt :ref:`installing` is complete.

2. Download the reduced slitless source images used to create Keck DEIMOS sensitivity functions from
   the `PypeIt dev-suite Google Drive`_. The files are located under ``DEIMOS_Dev/Throughput/throughput_gdw/data_products/extract``.

   .. code-block:: bash

        rclone -P copy remote:/PypeIt-development-suite/DEIMOS_Dev/Throughput/throughput_gdw/data_products/extract/ source_files/

3. Run the ``create_deimos_sensfuncs.py`` script.

   .. code-block:: bash

        mkdir sensfunc_files
        PypeIt-development-suite/sensfunc_archive/create_deimos_sensfuncs.py all source_files sensfunc_files

4. Copy the results to the PypeIt repository.

   .. code-block:: bash

        cd sensfunc_files
        cp keck_deimos_*.fits PypeIt/pypeit/data/sensfuncs/


