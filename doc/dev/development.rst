.. include:: ../include/links.rst

.. _development:

PypeIt Development Procedures and Guidelines
============================================

We encourage anyone to help us develop the ``PypeIt`` code base to better
suit your needs and to improve its algorithms.  If you do so, please
follow this list of procedures and guidelines.  In particular, please
note our :ref:`codeconduct`.

Installation
------------

If you plan to develop/alter the ``PypeIt`` code directly, you should
install the software via git.  You can fork the GitHub `repo`_, develop
within your own fork, and then perform a pull request from that fork.
Or you can clone the `repo`_ and work with it directly:

.. code-block:: console

    git clone https://github.com/pypeit/PypeIt
    cd PypeIt
    pip install -e .

It's important that you use the ``-e`` option when you run ``pip`` so that changes
made to your local installation are immediately available when you re-execute the code.
This will install the core dependencies, but not optional ones like PyQT, PySide2, or
any of the dependencies required for testing or building docs. To install everything you
might possibly need, do:

.. code-block:: console

    pip install -e ".[test,docs,pyside2,pyqt5,shapely]"

There is also a short-cut for this:

.. code-block:: console

    pip install -e ".[dev]"

.. warning::
    The quotes in the above installation lines are shell dependent.  For
    example, you need the quotes with zshell, but the quotes cause an error in
    bash.

It is very highly recommended to do your development install into its own clean
environment. See :ref:`installing` for examples of how to do this with either ``conda``
or ``virtualenv``.

This isn't required, but to simplify some of the commands below, I
assume that you have an environmental variable that points to your
installation of the ``PypeIt`` repository.  E.g., add the following to your
``~/.zshrc`` or ``~/.bash_profile``:

.. code-block:: console

    export PYPEIT_DIR=$HOME/Work/packages/PypeIt

Branches
--------

``PypeIt`` maintains two persistent branches:

 * ``release``: This is the primary stable version of the code.  Nominally, this
   is the version of the code that is installed when using `pip`_ and its the
   version that is used to produce latest `documentation`_.  Pull requests to
   this branch are only done before tagging a new release of the code or to
   perform critical bug hotfixes.  The release schedule is decided by the core
   developers.

 * ``develop``: This is the main development version of the code.  It should be
   stable enough to use, but it may contain experimental, unsupported code that
   is work in progress.

When editing the code, please create a new branch stemming from the
``develop`` branch:

.. code-block:: bash

    cd $PYPEIT_DIR
    git checkout develop
    git pull
    git checkout -b my_new_feature

Development Principles and Communication
----------------------------------------

The main thing to keep in mind when developing for ``PypeIt`` is that its
primary use is as an end-to-end reduction pipeline.  This has a few
implications that one should keep in mind:

 * By default, the execution of ``run_pypeit`` should continue either
   until a critical error is raised or the reduction is complete.  **No
   direct interaction with the code should be required at any point**,
   unless explicitly requested by the user.  ``PypeIt`` does have some
   interactive components, but most of these are executed via separate
   scripts.

 * Any input needed from the user for your feature should be provided by
   a parameter (preferred) or as a command-line argument.

 * When developing and debugging, you may need to interact with the code
   using `pdb`_ or `IPython.embed`_; however, these instances should be
   removed before performing a pull request.

 * The success or failure of any given procedure must be assessed via
   automatically generated quality-assessment figures (preferred) or via
   scripts that interact with the primary output files.

 * If your development includes adding a **new spectrograph** to the list of
   spectrograph data that ``PypeIt`` can reduce, see advice at :ref:`new_spec`.

 * If your development includes adding a **new executable script**, see advice
   at :ref:`new_script`.

Feature development in ``PypeIt`` is unlikely to be fully independent of
other development activities.  Your feature will likely depend on or
influence the outcome of other modules during the data-reduction
process.  This leads to a few important guidelines:

 * Make sure that your branch is always up-to-date with the ``develop``
   branch.  E.g.:

   .. code-block:: bash

        cd $PYPEIT_DIR
        git checkout develop
        git pull
        git checkout my_new_feature
        git merge --no-ff develop

 * Consider the effects of simultaneous development efforts on your work
   and vice versa.  For example, if you're working on a specific module
   of the code that depends on the result/datamodel of the
   wavelength-calibration module, you should communicate this and find
   out if someone else is developing that module and how/if they're
   changing it.  Depending on the scale of those changes, development
   priorities may need to be worked out to minimize merge conflicts and
   the need to immediately rework/refactor new code.

 * When you're ready to, you can submit a PR at any time, but the core
   development team will need to discuss the merge order to ensure a smooth
   process.

Our primary means of **communication** for development is the `PypeIt
developers Slack <https://pypeit.slack.com>`_.  Contact `X Prochaska`_ for
access.

Testing the Code
----------------

``PypeIt`` has two main methods for testing and verifying the code base,
unit tests and a dedicated development suite.

.. _dev-suite:

Development Suite
~~~~~~~~~~~~~~~~~

We have compiled a large suite of data from all instruments supported by
``PypeIt``.  These data are used to test that ``PypeIt`` is successful for *all*
instruments *anytime* new features are developed.  For access to the shared
Google TeamDrive, please contact `X Prochaska`_.

To test PypeIt using the data from the Google Drive:

 * Clone the `PypeIt-development-suite`_ repository:

   .. code-block:: bash

        git clone https://github.com/pypeit/PypeIt-development-suite.git

 * Download/sync the Drive to the repository.  The Drive ``CALIBS`` and
   ``RAW_DATA`` directories should be accessible.  For syncing, consider using
   `rclone`_.

   .. warning::

        The ``RAW_DATA`` directory currently contains about 43 Gb of data, and
        running the develop test below produces about 61 Gb (outdated) of
        reduced data.

 * Run the test suite on the setups designated for development purposes:

   .. code-block:: bash

        cd PypeIt-development-suite
        ./pypeit_test develop

   .. warning::

        The current script executes 89 (outdated) tests.  These tests are mostly
        reductions of example data, but they also include fluxing, flexure, and
        coadding tests.  The execution time is system dependent, but you should
        expect it to take approx. 14 hours (outdated).

   The test suite can be run in parallel using the -t option.

   .. code-block:: bash

        ./pypeit_test -t 2 develop

   The above runs tests using two processes run from parallel threads. The number of
   threads that can be used reliably depends on the amount of memory available. The
   below table shows the worst case amount of memory used per number of threads.

   +--------------+---------------+
   | # of threads | Memory used GB|
   +==============+===============+
   | 1            | 16            |
   +--------------+---------------+
   | 2            | 32            |
   +--------------+---------------+
   | 3            | 46            |
   +--------------+---------------+
   | 4            | 60            |
   +--------------+---------------+
   | 5            | 72            |
   +--------------+---------------+
   | 6            | 80            |
   +--------------+---------------+
   | 7            | 88            |
   +--------------+---------------+
   | 8            | 94            |
   +--------------+---------------+

To add new tests to the development suite

    #. Add the new data to shared Google Drive under ``RAW_DATA``. The tests are
       organized into setup directories under a directory named for the
       instrument.

    #. Add new a pypeit to the `PypeIt-development-suite`_ repo under
       ``pypeit_files``. The file name be lower case and named after the
       instrument and setup, for example: ``keck_deimos_1200g_m_7750.pypeit``.

    #. If desired, add any files for ``pypeit_sensfunc``, ``pypeit_flux_calib``,
       ``pypeit_coadd_1dspec``, ``pypeit_coadd_2dspec`` to the
       `PypeIt-development-suite`_ repo under ``sensfunc_files``,
       ``fluxing_files``, ``coadd1d_files``, ``coadd2d_files``, respectively.

    #. Edit ``test_setups.py`` in the `PypeIt-development-suite`_ under
       ``test_scripts``. Follow the instructions at the top of that file.

    #. Run the full development test suite to completion. Once all tests pass,
       the ``test_priority_file`` will be updated with the new test. This file
       tells the test scripts what order to run the tests in for optimum CPU
       utilization.  Commit ``test_priority_list`` and any other files added to
       the repository and submit a pull request.

.. _unit-tests:

Unit Tests
~~~~~~~~~~

Unit tests are located in the ``$PYPEIT_DIR/pypeit/tests`` directory.  To run
them, make sure you have `pytest`_ installed and then:

.. code-block:: bash

    cd $PYPEIT_DIR
    pytest

If some tests fail, you can run an individual test, e.g. ``test_wvcalib.py``
with

.. code-block:: bash

    cd $PYPEIT_DIR
    pytest -s pypeit/tests/test_wvcalib.py

Note that the "-s" option allows you to insert interactive debugging commands
into the test, here ``test_wvcalib.py``, to help determine why the test is
failing.

.. warning::

    Running these tests generates some files that should be ignored.  **Please
    do not add these test files to the repository.**  We're in the process of
    including some automatic clean-up in the testing functions.

Note also that the use of `pytest`_ requires the test dependencies to be
installed, e.g. via ``pip install -e .[test]``. It is also possible, and often
preferable, to run tests within their own isolated environments using `tox
<https://tox.readthedocs.io/en/latest/>`_. This provides the capability to
easily run tests against different versions of the various dependencies,
including different versions python. The available ``tox`` environments are
defined in ``$PYPEIT_DIR/tox.ini`` and can be listed by running ``tox -a``. To
run tests against the default dependencies using the default python, do:

.. code-block:: bash

    cd $PYPEIT_DIR
    tox -e test

To specify a python version, do something like:

.. code-block:: bash

    cd $PYPEIT_DIR
    tox -e py38-test

To test against, for example, the ``main`` branch for ``astropy`` on GitHub, you
can do:

.. code-block:: bash

    cd $PYPEIT_DIR
    tox -e py38-test-astropydev

Similar ``dev`` dependencies are configured for ``numpy``, ``ginga``, and
``linetools``, as well.

Some of the unit tests require the `Development Suite`_ and/or a set of "cooked"
data products with the expected result produced by testing ``PypeIt`` on
specific data.  The unit tests will skip the appropriate testing function if the
`Development Suite`_ path is not defined or if there is no ``Cooked/`` directory
in the relevant path.

For the unit tests to take advantage of the development-suite data,
``PYPEIT_DEV`` must be an environmental variable that points to the root
directory with the development suite data.  For example, include the following
line in your ``~/.zshrc`` or ``~/.bash_profile`` file:

.. code-block:: bash

    export PYPEIT_DEV=$HOME/Work/packages/PypeIt-development-suite

For unit tests that use the "cooked" data, ``PypeIt`` must find a directory
called ``$PYPEIT_DEV/Cooked/``.

Workflow
--------

A typical ``PypeIt`` development workflow is as follows:

 * Create a new branch stemming from the ``develop`` branch:

   .. code-block:: bash

        cd $PYPEIT_DIR
        git checkout develop
        git pull
        git checkout -b my_new_feature

 * Develop and debug the feature

 * Run the unit tests, fix any failures, add new tests that test your new
   feature(s), and/or modify the tests to accommodate your new feature:

   .. code-block:: bash

        cd $PYPEIT_DIR
        pytest

   or preferably:

   .. code-block:: bash

        cd $PYPEIT_DIR
        tox -e test

   The tests should be run so that they have access to the `Development
   Suite`_ (so that it can, e.g., test loading data), but this first
   round of tests can/should be run without the "cooked" output (e.g.,
   either delete or move the ``$PYPEIT_DEV/Cooked/`` directory).

 * Run the `Development Suite`_ and fix any failures:

   .. code-block:: bash

        cd $PYPEIT_DEV
        ./pypeit_test develop

 * (If needed) Build the cooked tar file (e.g., replace ``x.xx.x`` by
   incrementing the current version):

   .. code-block:: bash

        cd $PYPEIT_DEV
        ./build_cooked x.xx.x

 * Rerun and debug the tests (being sure to edit
   ``$PYPEIT_DIR/pypeit/tests/test_cooked.py`` to accept your cooked
   version):

   .. code-block:: bash

        cd $PYPEIT_DIR
        tox -e test

 * Edit ``$PYPEIT_DIR/CHANGES.rst`` to reflect your key developments and
   update the API `documentation`_.

   .. code-block:: bash

        cd $PYPEIT_DIR
        ./update_docs

 * Make sure all your edits are committed and pushed to the remote
   repository:

   .. code-block:: bash

        cd $PYPEIT_DIR
        git add -u
        git commit -m 'final prep for PR'
        git push

 * `Submit a Pull Request (PR)
   <https://github.com/pypeit/PypeIt/compare>`_. Unless otherwise
   requested, all PRs should be submitted to the ``develop`` branch.

Pull Request Acceptance Requirements
------------------------------------

Once you've submitted a pull request, we'll review your PR and provide
comments on the code.  The minimum requirements for acceptance of a PR
are as follows:

 * If your PR introduces a new instrument (see :ref:`new_spec`) that ``PypeIt``
   is to support for the long term, this instrument *must* be added to the
   `Development Suite`_.  That means raw data should be added to the Google
   Drive and a relevant test should be added to the ``$PYPEIT_DEV/pypeit_test``
   script (via a PR to the `PypeIt-development-suite`_) such that the new
   instrument is included in list of instruments tested by executing
   ``./pypeit_test develop``).

 * The Code Checks run by GitHub (see the Checks tab of the PR) on the remote
   repository must pass.

 * You have to post a successful report resulting from your execution of
   both the `Unit Tests`_ and the `Development Suite`_.  You should also
   have uploaded the "cooked" data to the TeamDrive; e.g.:

   .. code-block:: bash

        rclone copyto Cooked_pypeit_dev_vx.xx.x.tar.gz gdv:Cooked_pypeit_dev_vx.xx.x.tar.gz

   .. figure:: ../figures/tests_success.png

        Example posting of successful tests.

   For hotfixes, these tests can be circumvented at the discretion of
   the core developers in the cases where the hotfix is obviously
   correct.

 * All new methods and classes must be at least minimally documented.
   "Minimally documented" means that each method has a docstring that
   gives at least (1) a one sentence description of the purpose of the
   method, (2) a complete list of the required and optional arguments
   and their meaning, (3) a description of the returned objects, if
   there are any.  Documentation is expected to adhere to `Sphinx`_
   syntax; i.e., the docstrings should be `reStructuredText`_.  We
   accept both `Google-format docstrings`_ and `Numpy-format
   docstrings`_, but the former is preferred.

 * The docstrings for any changes to existing methods that were altered
   must have been modified so that they are up-to-date and accurate.

 * Spurious commented code used for debugging or testing is fine, but
   please let us know if you want it to be kept by adding a relevant
   comment, something like ``# TODO: Keep this around for now``, at the
   beginning of the commented block.  Otherwise, we're likely to remove
   the commented code when we come across it.

 * "Unsupported code," that is code that is experimental and still work
   in progress, should be minimized as much as is reasonable.  The
   relevant code block should be clearly marked as experimental or WIP,
   and it should not be executed by the main ``PypeIt`` executable,
   ``run_pypeit``.

 * At least two reviewers must accept the code.

Tagging Protocol
----------------

The core development team will regularly tag "release" versions of the
repository.  Tagging a release version of the code is triggered anytime the
development branch of the code is merged into the ``release`` branch.  The
tagging process is as follows:

 * At regular ``PypeIt`` telecons or over the ``PypeIt`` developers Slack, the
   core development team will decide to merge the ``develop`` branch into
   ``release``.

 * A branch is created off of ``develop`` (typically called ``staged``)
   and then a `PR <https://github.com/pypeit/PypeIt/compare>`_ is issued
   to merge ``staged`` into ``release``.  This merge must meet the same
   `Pull Request Acceptance Requirements`_ when merging new branches
   into ``develop``.  However, typically the unit and development-suite
   tests do not need to be run because ``develop`` has already passed
   these tests and ``staged`` will not make any substantive changes to
   the code.

   .. code-block:: bash

        cd $PYPEIT_DIR
        git checkout develop
        git pull
        git checkout -b staged

 * **Before being merged into master**, the code is tagged as follows:

    .. code-block:: bash

        cd $PYPEIT_DIR

        # Edit CHANGES.rst to reflect the ending of a development phase
        vi CHANGES.rst
        git add -u
        git commit -m 'tag CHANGES'

        # Create a tag of the form X.Y.Z (using 1.4.2 here as an example).
        # The current autogenerated version is found in pypeit/version.py
        # and the tag used here should be everything in front of the 'dev'
        # string in the current version.
        git tag 1.4.2

        # Push the changes and the new tag
        git push
        git push --tags

 * The PR is accepted and ``staged`` is merged into ``release``.  Hotfix
   branches are deleted, but the ``develop`` branch should not be.

 * The tag is released for `pip`_ installation.

    .. code-block:: bash

        # Make sure you have the most recent version of twine installed
        pip install twine --upgrade
        # Construct the pip distribution
        python setup.py sdist bdist_wheel
        # Test the upload
        twine upload --repository pypeit-test dist/*
        # Upload, this time it's for keeps
        twine upload --repository pypeit dist/*

    For the uploading, you need a ``~/.pypirc`` file that looks like this:

    .. code-block:: ini

        [distutils]
        index-servers =
            pypeit
            pypeit-test

        [pypeit]
        repository: https://upload.pypi.org/legacy/
        username = pypeit
        password = [ask for this]

        [pypeit-test]
        repository: https://test.pypi.org/legacy/
        username = pypeit
        password = [ask for this]

 * After a new version is uploaded to `pip`_, a new PR is automatically
   generated for the `conda-forge/pypeit-feedstock
   <https://github.com/conda-forge/pypeit-feedstock>`_ repository.  Follow the
   commit checklist there before merging that PR, which will trigger the release
   of a new pypeit package on conda-forge. For more information on how to
   manually update the conda-forge pypeit package, see :ref:`conda_forge`.

 * Now we need to advance the version of the code to a new development
   version and merge ``develop`` with the new ``release``:

    .. code-block:: bash

        cd $PYPEIT_DIR

        # Branch from release
        git checkout release
        git pull
        git checkout -b devup

        # Edit CHANGES.rst to begin the new dev section
        vi CHANGES.rst
        git add -u
        git commit -m 'CHANGES'

   The addition of this new commit will cause ``setuptools_scm`` to
   automatically increment the version based on the last tag that was pushed.
   This will be of the form ``{next_version}.dev{distance}+{scm letter}{revision
   hash}``.  See the `setuptools_scm documentation
   <https://github.com/pypa/setuptools_scm>`_ for more details.

 * Finally, a quick PR is issued that pulls ``devup`` into ``develop``.
   All of the `Pull Request Acceptance Requirements`_ should already be
   satisfied, meaning that the PR should be quickly accepted and merged.

DOI
---

If we wish to generate a new DOI for the code, it is as simple
as

 * Generate a `new release on GitHub
   <https://help.github.com/en/github/administering-a-repository/about-releases>`_.

 * Update the DOI in the README.md


----

This document was developed and mutually agreed upon by: Kyle Westfall,
J. Xavier Prochaska, Joseph Hennawi.

*Last Modified: 13 Jun 2021*

----


Additional Developer Links
--------------------------

.. Should move these rst files into the dev directory!

Here are some developer-specific docs:

.. toctree::
   :maxdepth: 1

   ../flow
   ../internals
   ../metadata
   ../new_script
   reports
   conda_forge
   fluxing
