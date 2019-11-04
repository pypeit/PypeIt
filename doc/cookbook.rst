===============
PypeIt Cookbook
===============

This document gives an overview on
how to run PypeIt, i.e. minimal detail is provided.
Notes on :doc:`installing` are found elsewhere.

We now also provide a set of Slides that provide a more
visual step-by-step.  Find them here at
the `PYPEIT HOWTO <https://tinyurl.com/pypeit-howto>`_
These should be considered to contain
the most up-to-date information.

The following outlines the standard steps for running
PypeIt on a batch of data.  There are alternate ways to
run these steps, but non-experts should adhere to the
following approach.

Outline
+++++++

Here is the basic outline of the work flow.  The
following is for one instrument in one working directory.

1. Organize/Prepare your data

  - Identify folder(s) with raw images
  - The raw images can be gzip compressed although the Python FITS reader works much more slowly on gzipped files
  - We will refer to that folder as RAWDIR

2. Run the :ref:`pypeit_setup` *without* the --custom option to handle instrument :doc:`setup`.

   Inputs are the path to the raw data with the data prefix (e.g. lrisb) and then
   one of the PypeIt-approved :doc:`instruments` (e.g. keck_lris_blue, shane_kast_red).
   Here is an example::

    pypeit_setup -r /full_path/RAWDIR/lrisb  -s keck_lris_blue

   This does the following:

 - Generates a setup_files/ folder that holds two files
 - Generates a dummy PypeIt reduction file within the folder [ignore it]
 - Generates a .sorted file which lists files sorted by setup

 You should scan the output WARNING messages for insufficient calibration files (e.g. missing arc frames)

3. Inspect the :ref:`sorted-file` to confirm the expected instrument configuration(s)

  - If needed, add more files to your RAWDIR
  - If you do, repeat Step 2 above

4. Run :ref:`pypeit_setup` *with* the --custom option

  This produces one folder per setup and a custom :doc:`pypeit_file`.
  Here is an example of the call::

    pypeit_setup -r /full_path/RAWDIR/lrisb  -s keck_lris_blue -c=all

  This generates one folder per setup and a unique :doc:`pypeit_file` file in each folder.


5. Prepare the custom :doc:`pypeit_file` for reducing a given setup

  - Enter one of the setup folders (e.g. kast_lris_blue_A)
  - Modify the custom :doc:`pypeit_file` as needed

    - trim/add calibration files
    - edit frametypes
    - Modify user-defined execution parameters

6. Run the reduction (described in :doc:`running`)

  - :ref:`run-pypeit` PypeIt_file
  - Hope for the best...  :)

7. Examine QA (:doc:`qa`)

  - When an exposure is fully reduced, a QA file (PDF) is generated in the QA folder
  - Examine the output as described in the :doc:`qa` documentation

8. Examine spectra
  - Examine the extracted 1D spectra with :ref:`pypeit-1dspec`
  - Examine the extracted 2D spectra with :ref:`pypeit-2dspec`

9.  Flux

10. Coadd (see :doc:`coadding`)

11. Repeat steps 5-10 for additional setups, as desired




