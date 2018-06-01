.. highlight:: rest

==============
PYPIT Cookbook
==============

This document gives an overview on
how to run PYPIT, i.e. minimal detail is provided.
Notes on :doc:`installing` are found elsewhere.

The following outlines the standard steps for running
PYPIT on a batch of data.  There are alternate ways to
run these steps, but non-experts should adhere to the
following approach.

Outline
+++++++

Here is the basic outline of the work flow.  The
following is for one instrument in one working directory.

1. Prepare your data

  - Identify folder(s) with raw images
  - The raw images can be gzip compressed although the Python FITS reader works much more slowly on gzipped files
  - We will refer to that folder as RAWDIR

2. Run the :ref:`pypit_setup` *without* the --custom option to handle instrument :doc:`setup`.

   Inputs are the path to the raw data with the data prefix (e.g. lrisb) and then
   one of the PYPIT-approved :doc:`instruments` (e.g. keck_lris_blue, shane_kast_red).
   Here is an example::

    pypit_setup /full_path/RAWDIR/lrisb  keck_lris_blue

   This does the following:

 - Generates a setup_files/ folder that holds a series of files
 - Generates the instrument PYPIT reduction file [not used further]
 - Generates the instrument .setups file (see :doc:`setup`)
 - Generates a .sorted file which lists files sorted by setup

 You should scan the output WARNING messages for insufficient calibration files (e.g. missing arc frames)

3. Inspect the :ref:`setups-file` and :ref:`sorted-file` to confirm the instrument configurations

  - If needed, add more files to your RAWDIR
  - If you do, repeat Step 2 above

4. Run :ref:`pypit_setup` *with* the --custom option

  This produces one folder per setup and a custom :doc:`pypit_file`.
  Here is an example of the call::

    pypit_setup /full_path/RAWDIR/lrisb  keck_lris_blue --custom


5. Prepare the custom :doc:`pypit_file` for reducing a given setup

  - Enter one of the setup folders (e.g. kast_lris_blue_setup_A)
  - Modify the custom :doc:`pypit_file` as needed

    - trim/add calibration files
    - edit frametypes

6. Run :ref:`run-calcheck` on the custom PYPIT file(s) (described in :doc:`calcheck`)

  - Modify the spect block in the PYPIT file to specify calibrations
  - Inspect the .calibs file for your PYPIT file.
  - Confirm calibration, science and standard frames
  - Further customize your PYPIT file, as needed

7. Run the reduction (described in :doc:`running`)

  - :ref:`run-pypit` PYPIT_file
  - Hope for the best...  :)

8. Examine QA (:doc:`qa`)

  - When an exposure is fully reduced, a QA file (PDF) is generated in the QA folder
  - Examine the output as described in the :doc:`qa` documentation

9. Examine spectra
  - Examine the extracted 1D spectra with :ref:`pypit-1dspec`
  - Examine the extracted 2D spectra with :ref:`pypit-2dspec`

10. Coadd (see :doc:`coadding`)

11. Repeat steps 4-8 for additional setups, as desired




