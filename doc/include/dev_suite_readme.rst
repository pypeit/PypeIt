The Python Spectroscopic Data Reduction Pipeline (PypeIt) development suite
===========================================================================

This repository provides data and input files used for testing `PypeIt
<https://github.com/pypeit/PypeIt>`__.

Installation
------------

To install, run:

.. code-block:: console

    git clone https://github.com/pypeit/PypeIt-development-suite.git

This pulls the versioned files into a directory called
``PypeIt-development-suite``.  Both the testing script in this repo and
the tests that are included with the PypeIt code expect an environmental
variable that points to this directory.  For example,

 - If you use c-shell, add the following to your ``~/.cshrc``:

   .. code-block:: console
   
    setenv PYPEIT_DEV ${HOME}/PypeIt-development-suite

 - If you use bash, add the following to your ``~/.bashrc`` or
   ``~/.bash_profile``:

   .. code-block:: console
   
    export PYPEIT_DEV=${HOME}/PypeIt-development-suite

Data Access
-----------

Given its volume, this repo does not contain the raw data.  Instead the
raw test data are hosted by in a 
`Google Drive <https://drive.google.com/drive/folders/1oh19siB1-F0jjmY-F_jr73eA-TQYEiFW?usp=sharing>`__.

The easiest way to access the
development suite data is download and use `Google Drive
for desktop <https://support.google.com/googleone/answer/10838124?visit_id=637915333936129509-3533094830&rd=1>`__.  Google
Drive for desktop is a Dropbox-like application that syncs your Google Drive
with a local directory on your machine.  

Using Google Drive for desktop, the you will be able to
access the development suite data at the path:

``/Volumes/GoogleDrive/My\ Drive/PypeIt-development-suite/``

If ``PypeIt-development-suite`` does not appear under ``My Drive``, go to the Google Drive web interface, find it under ``Shared with me``,
right click it, and select ``Add shortcut to drive``.

If you cannot or would rather
not use Google Drive for desktop, you can simply download the appropriate files directly using the Google Drive web interface (although this can be a bit onerous and does not keep in sync with the remote Team Drive). Alternately you can use `rclone <https://rclone.org/>`__ to manually copy or sync with the Google Drive.


The Google Drive contains two directories that should be accessible for
testing PypeIt (see below): ``RAW_DATA``, and ``CALIBS``.

  - If you download the data directly from the Google Drive, place them in
    your ``PypeIt-development-suite`` directory.  **Make sure that you *do
    not* add these directories to the repo!**

  - If you're using Google File Stream, add symlinks to your
    ``PypeIt-development-suite`` directory as follows (be sure to include
    the \ in the My\ Drive otherwise the space in "My Drive" will
    cause problems):

    .. code-block:: console

        cd $PYPEIT_DEV
        ln -s /Volumes/GoogleDrive/My\ Drive/PypeIt-development-suite/RAW_DATA  RAW_DATA
        ln -s /Volumes/GoogleDrive/My\ Drive/PypeIt-development-suite/CALIBS  CALIBS

Testing PypeIt
--------------

This suite includes several examples of the PypeIt data reduction process 
for some of the supported instruments.

If you have `PypeIt <https://github.com/pypeit/PypeIt>`__ installed, you
should be able to run the tests in this directory using the
``pypeit_test`` script. It has one required argument, the type of tests to run.

.. code-block:: console

    $ $PYPEIT_DEV/pypeit_test all

The above runs all pypeit tests, including the unit tests in the PypeIt 
installation. To run a subset of tests you can specify one or more test type. For example:

.. code-block:: console

    # Run the non-pytest tests, similar to develop in prior versions
    $ $PYPEIT_DEV/pypeit_test reduce afterburn ql

    # Run reduce tests followed by the vet pytest tests
    $ $PYPEIT_DEV/pypeit_test reduce vet

All of the supported test test types are shown in the table below.

============= ==============================================================================================================
Test Type     Description
============= ==============================================================================================================
pypeit_tests  Runs the pytest tests installed with PypeIt in pypeit/tests. These tests are self contained and can run in CI.
unit          Runs the pytest tests installed in $PYPEIT_DEV/unit. These tests require the RAW_DATA directory.
reduce        Runs the reduction tests that call run_pypeit directly.
afterburn     Runs the PypeIt tests that directly call PypeIt post-reduction  tools (e.g. flux calibration, coadding, etc).
ql            Runs the Quick Look tests.
vet           Runs the pytest tests that verify the results from earlier PypeIt tests.
all           Runs all of the above, in the order listed above.
list          This does not run any tests, instead it lists all of the supported instruments and setups. (See below).
============= ==============================================================================================================

Running Unit and Vet tests separately
-------------------------------------

The unit and vet tests can also be run directly using pytest. For example:

.. code-block:: console

    $ cd $PYPEIT_DEV

    # Run all dev-suite unit tests
    $ pytest unit_tests  

    # Run all dev-suite vet tests
    $ pytest vet_tests   

    # Run the script tests in both unit_tests and vet_tests
    $ pytest unit_tests/test_scripts.py vet_tests/test_scripts.py

See the `pytest docs <https://docs.pytest.org/>`__ for more information on running pytest.

Selecting test setups and instruments to test
---------------------------------------------

The ``-i`` and ``-s`` options to ``pypeit_test`` can be used to select multiple
instruments and setups to test. Setups can be specified in conjunction with
an instrument or by using a ``/`` between instrument and setup name. For example:

.. code-block:: console

    # Run all of pytest tests and all of the shane_kast_blue and shane_kast_red tests
    cd $PYPEIT_DEV
    $ ./pypeit_test all -i shane_kast_blue shane_kast_red

    # Run one reduce test from shane_kast_blue and shane_kast_red respectively
    $ ./pypeit_test reduce -i shane_kast_blue shane_kast_red -s 600_4310_d55 600_7500_d57

    # Run the same tests as above, using the / syntax
    $ ./pypeit_test reduce -s shane_kast_blue/600_4310_d55 shane_kast_red/600_7500_d57

Run ``pypeit_test list`` to see a list of all supported instruments and setups.

Test Reports
------------

A test report can be generated using the ``-r`` option, for example:

.. code-block:: console

    $PYPEIT_DEV/pypeit_test all -r test_report.txt

The contents of the report contain complete pytest output and additional information about the test setups run. Below is an example of some of the output:

.. code-block:: console

    PypeIt Unit Tests Results:
    -------------------------
    ============================= test session starts ==============================
    platform linux -- Python 3.8.5, pytest-6.2.4, py-1.10.0, pluggy-0.13.1 -- /root/miniconda3/bin/python
    cachedir: .pytest_cache
    hypothesis profile 'default' -> database=DirectoryBasedExampleDatabase('/tmp/REDUX_OUT/.hypothesis/examples')
    rootdir: /PypeIt, configfile: setup.cfg
    plugins: hypothesis-6.14.2, arraydiff-0.3, astropy-header-0.1.2, cov-2.12.1, doctestplus-0.10.0, filter-subpackage-0.1.1, openfiles-0.5.0, remotedata-0.3.2
    collecting ... collected 193 items

    ../../PypeIt/pypeit/tests/test_alignments.py::test_alignments PASSED     [  0%]
    ../../PypeIt/pypeit/tests/test_arc.py::test_detect_lines PASSED          [  1%]
    ../../PypeIt/pypeit/tests/test_archive.py::test_archive_meta PASSED      [  1%]
    ../../PypeIt/pypeit/tests/test_archive.py::test_archive_dir PASSED       [  2%]

    ...

    Reduced data for the following setups:
        shane_kast_blue/452_3306_d57
        shane_kast_blue/600_4310_d55
        shane_kast_blue/830_3460_d46

    Ran tests in 4 parallel processes

    -------------------------
    Test Setup: shane_kast_blue/830_3460_d46

    -------------------------
    Directories:
            Raw data: /tmp/PypeIt-development-suite/RAW_DATA/shane_kast_blue/830_3460_d46
        PypeIt output: /tmp/REDUX_OUT/shane_kast_blue/830_3460_d46
    Files:
        .pypeit file: None
    Std .pypeit file: None
    Tests:
    ----
    shane_kast_blue/830_3460_d46 pypeit  Result: --- PASSED

    Logfile:    /tmp/REDUX_OUT/shane_kast_blue/830_3460_d46/shane_kast_blue_830_3460_d46.test.log
    Process Id: 1188
    Start time: Thu Jan  5 07:20:05 2023
    End time:   Thu Jan  5 07:21:18 2023
    Duration:   0:01:13.575440
    Mem Usage:  1100517376
    Command:    run_pypeit /tmp/REDUX_OUT/shane_kast_blue/830_3460_d46/shane_kast_blue_830_3460_d46.pypeit -o

    Error Messages:

    End of Log:
    [INFO]    :: run_pypeit.py 111 main() - Generating QA HTML
    Wrote: /tmp/REDUX_OUT/shane_kast_blue/830_3460_d46/QA/MF_A.html
    Wrote: /tmp/REDUX_OUT/shane_kast_blue/830_3460_d46/QA/MF_A.html

    ...

    Test Summary
    --------------------------------------------------------
    --- PYTEST PYPEIT UNIT TESTS PASSED  205 passed, 137 warnings in 212.65s (0:03:32) ---
    --- PYTEST UNIT TESTS PASSED  118 passed, 260 warnings in 956.59s (0:15:56) ---
    --- PYTEST VET TESTS PASSED  29 passed, 160 warnings in 813.29s (0:13:33) ---
    --- PYPEIT DEVELOPMENT SUITE PASSED 155/155 TESTS  ---
    Coverage results:
    TOTAL                                              41653  12057    71%

    Testing Started at 2023-01-04T21:35:25.696424
    Testing Completed at 2023-01-05T07:35:30.909205
    Total Time: 10:00:05.212781

Code coverage
-------------

The dev suite can also collect coverage data for ``PypeIt`` using
`Coverage <https://coverage.readthedocs.io/>`__ . To do this add
``--coverage <coverage report file>`` to the ``pypeit_test`` command:

.. code-block:: console

    $ cd $PYPEIT_DEV
    $ ./pypeit_test all --coverage coverage_report_file.txt

The coverage report contains a file by file list of the coverage information, including missed lines. It ends with a summary of the total code coverage.
The unit tests, and deprecated sections of the ``PypeIt`` code base are omitted.
For example:

.. code-block:: console

    Name                                                               Stmts   Miss  Cover   Missing
    ------------------------------------------------------------------------------------------------
    /home/dusty/work/PypeIt/pypeit/__init__.py                            20      5    75%   45-49
    ...
    /home/dusty/work/PypeIt/pypeit/wavemodel.py                          330    286    13%   50-79, 108-137, 149-152, 195-211, 227-233, 247-257, 315-418, 459-535, 568-602, 627-648, 686-715, 765-788, 841-863
    /home/dusty/work/PypeIt/pypeit/wavetilts.py                          240    104    57%   90, 204-205, 297, 311-319, 435-525, 562-581, 593, 598-600, 611-614, 627-630, 635, 661-681, 688, 714-717, 722-729
    ------------------------------------------------------------------------------------------------
    TOTAL                                                              41785  22139    47%

Performance Statistics
----------------------

The dev suite collects performance information about how fast each test
runs and how much memory it uses. This can be output to a CSV file using
the ``-csv`` option.

.. code-block:: console

    $ cd $PYPEIT_DEV
    $ ./pypeit_test all --csv performance.csv

The output contains the duration of each test and its memory usage. For example (spaces added for clarity):

.. code-block:: console

    Setup,                       Test Type,        Start Time,                 End Time,                   Duration(s), Memory Usage (bytes), Duration (D:H:M:S), Memory Usage (MiB)
    gemini_gmos/GS_HAM_R400_860, pypeit,           2023-01-06 13:18:44.050069, 2023-01-06 14:42:15.791558, 5011.741489, 5697548288,           1:23:31.741489,     5433.60546875
    gemini_gmos/GS_HAM_R400_860, pypeit standards, 2023-01-06 14:42:15.791808, 2023-01-06 14:45:43.122017, 207.330209,  2606161920,           0:03:27.330209,     2485.4296875
    gemini_gmos/GS_HAM_R400_860, pypeit_sensfunc,  2023-01-06 14:45:43.124685, 2023-01-06 14:47:07.308025, 84.18334,    272601088,            0:01:24.183340,     259.97265625
    gemini_gmos/GS_HAM_R400_860, pypeit_flux_setup,2023-01-06 14:47:07.308275, 2023-01-06 14:47:08.978942, 1.670667,    249856,               0:00:01.670667,     0.23828125
    gemini_gmos/GS_HAM_R400_860, pypeit_flux,      2023-01-06 14:47:08.979198, 2023-01-06 14:47:12.503334, 3.524136,    210182144,            0:00:03.524136,     200.4453125

Parallel Testing
----------------

The development suite currently takes over 12 hours to run. This can be
sped up by running parallel tests:

.. code-block:: console

    ./pypeit_test -t 2 all

The number of threads that can be run depends on the amount of memory
available. Based on testing, the memory requirements of the devsuite are:

+------------+--------------+
| Number of  | Memory Usage |
|  Threads   |              |
+============+==============+
|        1   |  20 GiB      |
+------------+--------------+
|        2   |  40 GiB      |
+------------+--------------+
|        4   |  58 GiB      |
+------------+--------------+
|        6   |  78 GiB      |
+------------+--------------+
|        8   |  98 GiB      |
+------------+--------------+
|       12   | 133 GiB      |
+------------+--------------+
|       16   | 148 GiB      |
+------------+--------------+

For systems with more virtual CPUs than physical CPU
cores (i.e. Hyperthreading) the number of threads should not exceed the
number of physical cores, or else there could be a performance hit as
threads compete for resources.  

To keep all cpus active as long as possible ``pypeit_test`` runs the
slowest tests first. To do this it needs the ``test_priority_list`` file
which contains a list of all the test setups ordered from slowest to
fastest. This file is re-written everytime a run of the full test suite
passes, and should be kept up to date by periodically pushing it to git.

The pytest portion of the dev-suite currently cannot be run in parallel.

Headless Testing
----------------

Some of the tests in the dev-suite will start GUI applications. To run
in a headless environment where this isn't possible, QT must still be
installed.  To do so, first install ``PypeIt`` with QT5 support, as
documented `here
<https://pypeit.readthedocs.io/en/latest/installing.html>`__. Next
install the correct QT package for the OS.

=================  ================
OS Distribution    Package name
=================  ================
Ubuntu 21.04       qt5-default
Ubuntu 22.04       qtbase5-dev
Centos 7           qt5-qtbase-devel
=================  ================

Finally, set the ``QT_QPA_PLATFORM`` environment variable to
``offscreen``. There is a convenience file named
``source_headless_test.sh`` that will do this. For example:

.. code-block:: console

    source $PYPEIT_DEV/source_headless_test.sh

Running in Nautilus
-------------------

The dev-suite can be run in the `Nautilus cluster <https://ucsd-prp.gitlab.io/>`__.
To generate the YAML for a dev-suite job, use ``gen_kube_devsuite``.  If needed,
a specific branch of both the PypeIt repository and the Pypeit-development-suite
repository can be chosen using ``-p`` and ``-d`` respectively.  These default to
``develop`` if not specified. The YAML file can then be ran using ``kubectl``.

For example to run using the pypeit_branch on the PypeIt repo and the
devsuite_branch on the PypeIt-development-suite repo:

.. code-block:: console

    $ cd $PYPEIT_DEV/nauitilus
    $ ./gen_kube_devsuite devsuite-job-name devsuite_job_file.yml -p pypeit_branch -d devsuite_branch
    $ kubectl create -f devsuite_job_file.yml

The results of the dev-suite are copied to Nautilus S3 under
``s3://pypeit/Reports/``. And can be retrieved using the AWS CLI as
follows:

.. code-block:: console

    aws --endpoint  https://s3-west.nrp-nautilus.io s3 cp s3://pypeit/Reports/devsuite-job-name.report .

``rclone`` can also be used access the Nautilus S3 storage. When
configuring use ``https://s3-west.nrp-nautilus.io`` as the endpoint.

``gen_kube_devsuite`` has additional code for generating coverage
information and the test priority list. If ``--coverage`` and
``--priority_list`` are used, these files are also copied to S3:

.. code-block:: console

    $ ./gen_kube_devsuite coverage-job-name coverage_job_file.yml --coverage 
    $ ./gen_kube_devsuite priority-job-name priority_job_file.yml --priority_list
    $ kubectl create -f coverage_job_file.yml
    $ kubectl create -f priority_job_file.yml
    $ # Wait several hours 
    $ export ENDPOINT=https://s3-west.nrp-nautilus.io 
    $ aws --endpoint $ENDPOINT s3 cp s3://pypeit/Reports coverage-job-name.report .
    $ aws --endpoint $ENDPOINT s3 cp s3://pypeit/Reports/priority-job-name.report .
    $ aws --endpoint $ENDPOINT s3 cp s3://pypeit/Reports/coverage-job-name.coverage.report .
    $ aws --endpoint $ENDPOINT s3 cp s3://pypeit/Reports/priority-job-name.test_priority_list .

Notice that ``--coverage`` can affect the performance of tests, so it's best
not to run it and ``--priority_list`` together.

To monitor a test in Nautilus as it is running, the logs can be tailed:

.. code-block:: console

    $ kubectl get pods

    NAME                         READY    STATUS    RESTARTS
    devsuite-job-name--1-fpjxw   1/1      RUNNING   0

    $ kubectl logs -f devsuite-job-name--1-fpjxw

Additionally they can be monitored with the `Nautilus Grafana page <https://grafana.nrp-nautilus.io/?orgId=1>`__.

By default ``gen_kube_devsuite`` creates a job using a default container with PypeIt 
pre-installed. It also supports running with different python versions by
selecting a different container. For example:

.. code-block:: console

    $ ./gen_kube_devsuite devsuite-python3.11-job devsuite-python3.11-job.yml --container python3.11

Any of the standard python images in docker hub can be used. To use a different container the full 
download path must be given. For example:

.. code-block:: console

    $ ./gen_kube_devsuite devsuite-ubuntu-job devsuite-ubuntu-job.yml --container docker.io/library/ubuntu:22.04


Additional Options
------------------

.. code-block:: console

    $ $PYPEIT_DEV/pypeit_test -h
    usage: pypeit_test [-h] [-o OUTPUTDIR] [-i INSTRUMENTS [INSTRUMENTS ...]]
                   [-s SETUPS [SETUPS ...]] [--debug] [-p] [-m] [-t THREADS]
                   [-q] [-v] [--coverage COVERAGE] [-r REPORT] [-c CSV] [-w]
                   tests [tests ...]

    Run pypeit tests on a set of instruments. Typical call for testing pypeit when
    developing new code is `./pypeit_test all`. Execution requires you to have a
    PYPEIT_DEV environmental variable, pointing to the top-level directory of the
    dev-suite repository (typically the location of this script). Raw data for
    testing is expected to be at ${PYPEIT_DEV}/RAW_DATA. To run all tests for the
    supported instruments, use 'all'. To only run the basic reductions, use
    'reduce'. To only run the tests that use the results of the reductions, use
    'afterburn''. Use 'list' to view all supported setups.

    positional arguments:
    tests                 Which test types to run. Options are: pypeit_tests,
                            unit, reduce, afterburn, ql, vet, or all. Use list to
                            show all supported instruments and setups.

    optional arguments:
    -h, --help            show this help message and exit
    -o OUTPUTDIR, --outputdir OUTPUTDIR
                            Output folder. (default: REDUX_OUT)
    -i INSTRUMENTS [INSTRUMENTS ...], --instruments INSTRUMENTS [INSTRUMENTS ...]
                            One or more instruments to run tests for. Use
                            "pypeit_test list" to see all supported instruments.
                            (default: None)
    -s SETUPS [SETUPS ...], --setups SETUPS [SETUPS ...]
                            One or more setups to run tests for. Use "pypeit_test
                            list" to see all supported setups. (default: None)
    --debug               Debug using only blue setups (default: False)
    -p, --prep_only       Only prepare to execute run_pypeit, but do not
                            actually run it. (default: False)
    -m, --do_not_reuse_calibs
                            run pypeit without using any existing calibrations
                            (default: False)
    -t THREADS, --threads THREADS
                            Run THREADS number of parallel tests. (default: 1)
    -q, --quiet           Supress all output to stdout. If -r is not a given, a
                            report file will be written to
                            <outputdir>/pypeit_test_results.txt (default: False)
    -v, --verbose         Output additional detailed information while running
                            the tests and output a detailed report at the end of
                            testing. This has no effect if -q is given (default:
                            False)
    --coverage COVERAGE   Collect code coverage information. and write it to the
                            given file. (default: None)
    -r REPORT, --report REPORT
                            Write a detailed test report to REPORT. (default:
                            None)
    -c CSV, --csv CSV     Write performance numbers to a CSV file. (default:
                            None)
    -w, --show_warnings   Show warnings when running unit tests and vet tests.
                            (default: False)

