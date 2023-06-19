.. include:: ../include/links.rst

.. _development:

Development Procedures and Guidelines
=====================================

We encourage anyone to help us develop the ``PypeIt`` code base to better
suit your needs and to improve its algorithms.  If you do so, please
follow this list of procedures and guidelines.  In particular, please
note our :ref:`codeconduct`.

Installation
------------

If you plan to develop/alter the ``PypeIt`` code directly, you should install
the software via git, instead of pip or conda.  See :ref:`developer_install`.

This isn't required, but to simplify some of the commands below, I assume that
you have an environmental variable that points to your installation of the
``PypeIt`` repository.  E.g., add something like the following to your
``~/.zshrc`` or ``~/.bash_profile``:

.. code-block:: console

    export PYPEIT_DIR=$HOME/Work/packages/PypeIt

Branches
--------

``PypeIt`` maintains two persistent branches:

 * ``release``: This is the primary stable version of the code.  Modulo any
   recent hotfixes, this is the version of the code that is installed when
   using `pip`_ and its the version that is used to produce the latest
   `documentation`_.  Pull requests to this branch are only done before tagging
   a new release of the code or to perform critical bug hotfixes.  The release
   schedule is irregular and decided by the core developers.

 * ``develop``: This is the main development version of the code.  It should be
   stable enough to use, but it may contain experimental, unsupported code that
   is work in progress.

When editing the code, please create a new branch stemming from the ``develop``
branch.  You should also pull and merge in the most recent version of the
``release`` branch to make sure your new branch includes any very recent
hot-fixes.  On the command line, you can do this as follows:

.. code-block:: bash

    cd $PYPEIT_DIR
    git checkout release
    git pull
    git checkout develop
    git pull
    git checkout -b my_new_feature
    git merge --no-ff release

.. note::

    In terms of the merge with the release branch, beware that you may need to
    start a new section of the ``CHANGES.rst`` file to reflect a jump in the
    version number.  This should only be necessary if your branch is the first
    one after a new tag is released.

Development Principles and Communication
----------------------------------------

The main thing to keep in mind when developing for ``PypeIt`` is that its
primary use is as an end-to-end reduction pipeline.  This has a few
implications:

 * By default, the execution of ``run_pypeit`` should continue either
   until a critical error is raised or the reduction is complete.  **No direct
   interaction with the code should be required at any point**.  ``PypeIt`` does
   have some interactive components, but these are executed via separate
   scripts.

 * Any input needed from the user for your feature should be provided by
   :ref:`parameters` (preferred) or as a command-line argument.

 * When developing and debugging, you may need to interact with the code
   using `pdb`_ or `IPython.embed`_; however, these instances should be
   removed before performing a pull request.

 * The success or failure of any given procedure must be assessed via
   automatically generated quality-assessment figures (preferred) or via
   scripts that interact with the primary output files.

 * See :doc:`here <new_spectrograph>` for guidance on 
   adding a **new spectrograph** to the list of spectrograph data that
   ``PypeIt`` can reduce.

 * If your development includes adding a **new executable script**, see advice
   at :ref:`new_script`.

Feature development in ``PypeIt`` is unlikely to be fully independent of
other development activities.  Your feature will likely depend on or
influence the outcome of other modules during the data-reduction
process.  This leads to a few important guidelines:

 * Make sure that your branch is always up-to-date with the ``develop`` *and*
   ``release`` branches.  E.g.:

   .. code-block:: bash

        cd $PYPEIT_DIR
        git checkout release
        git pull
        git checkout develop
        git pull
        git checkout my_new_feature
        git merge --no-ff release
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

Our primary means of **communication** for development is the `PypeIt developers
Slack <https://pypeit.slack.com>`_ and a biweekly telecon.  Contact `X
Prochaska`_ for Slack access and/or the relevant Zoom link.

Testing the Code
----------------

``PypeIt`` performs extensive testing using the :ref:`dev-suite`; follow that
link for more details on executing the tests.  What follows describes how to add
new tests.

.. _dev-suite-tests:

Development Suite
~~~~~~~~~~~~~~~~~

To add new tests to the development suite

    #. Add the new data to shared Google Drive under ``RAW_DATA``. The tests are
       organized into setup directories under a directory named for the
       instrument.

    #. Add a new :ref:`pypeit_file` specific to this data to the `PypeIt
       Development Suite`_ repo under ``pypeit_files``. The file name must be
       lower case and named after the instrument and setup, for example:
       ``keck_deimos_1200g_m_7750.pypeit``.

    #. If desired, add any files for ``pypeit_sensfunc``, ``pypeit_flux_calib``,
       ``pypeit_coadd_1dspec``, ``pypeit_coadd_2dspec`` to the
       `PypeIt Development Suite`_ repo under ``sensfunc_files``,
       ``fluxing_files``, ``coadd1d_files``, ``coadd2d_files``, respectively.

    #. Edit ``test_setups.py`` in the `PypeIt Development Suite`_ under
       ``test_scripts``. Follow the instructions at the top of that file.

    #. Run the full development test suite to completion. Once all tests pass,
       the ``test_priority_file`` will be updated with the new test. This file
       tells the test scripts what order to run the tests in for optimum CPU
       utilization.  Commit ``test_priority_list`` and any other files added to
       the dev-suite repository and submit a pull request.

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
    do not add these test files to the repository.**  We try to include clean-up
    as part of the tests, but these are not always caught.

Note also that the use of `pytest`_ requires the test dependencies to be
installed, e.g. via ``pip install -e .[test]``. It is also possible, and often
preferable, to run tests within their own isolated environments using `tox
<https://tox.readthedocs.io/en/latest/>`_. This provides the capability to
easily run tests against different versions of the various dependencies,
including different python versions. The available ``tox`` environments are
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

Unit tests included in the main PypeIt repo should *not* require large data
files.  Some files are kept in the repo for this purpose (see the
``pypeit/tests/files`` directory), but they should be minimized to keep the size
of the package distribution manageable.  Unit tests that require input data
files should instead be added to the `PypeIt Development Suite`_.

Workflow
--------

A typical ``PypeIt`` development workflow is as follows:

 * Create a new branch stemming from the ``develop`` branch (hot-fixes should
   instead branch from ``release``):

   .. code-block:: bash

        cd $PYPEIT_DIR
        git checkout release
        git pull
        git checkout develop
        git pull
        git checkout -b my_new_feature
        git merge --no-ff release

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

 * Run the `Development Suite`_ and fix any failures:

   .. code-block:: bash

        cd $PYPEIT_DEV
        ./pypeit_test develop

   .. warning::

        The `Development Suite`_ is *extensive* and takes significant computing
        resources and time.  The PypeIt development team consistently executes
        these tests using cloud computing.  We recommend you ensure that your
        pypeit branch successfully runs on either a specific instrument of
        interest or ``shane_kast_blue`` first, and then someone on the PypeIt
        development team can execute the tests in the cloud.  From the top-level
        directory of the `Development Suite`_, you can run all tests for
        ``shane_kast_blue`` as follows:

        .. code-block:: bash

            ./pypeit_test all -i shane_kast_blue

 * Edit ``$PYPEIT_DIR/CHANGES.rst`` to reflect your key developments and
   update the `documentation`_.  You can compile the docs using the
   ``update_docs`` script (see below), which is just a simple convenience script
   for executing ``make clean ; make html`` in the ``doc`` directory.

   .. code-block:: bash

        cd $PYPEIT_DIR
        ./update_docs

   *Any* warnings in the sphinx build of the docs *must* be fixed.  If you're
   having difficulty getting the right sphinx/rst incantation, ping the
   documentation channel in the `PypeIt Developers Slack
   <https://pypeit.slack.com>`__.

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


.. note::

   The addition of new commits causes ``setuptools_scm`` to automatically
   increment the version based on the last tag that was pushed.  This will be of
   the form ``{next_version}.dev{distance}+{scm letter}{revision hash}``.  See
   the `setuptools_scm documentation <https://github.com/pypa/setuptools_scm>`_
   for more details.

Pull Request Acceptance Requirements
------------------------------------

Once you've submitted a pull request, we'll review your PR and provide
comments on the code.  The minimum requirements for acceptance of a PR
are as follows:

 * If your PR introduces a new instrument (see :ref:`new_spec`) that ``PypeIt``
   is to support for the long term, this instrument *must* be added to the
   `Development Suite`_.  That means raw data should be added to the Google
   Drive (see :ref:`here <dev-suite>`) and relevant tests should be added to
   the ``$PYPEIT_DEV/pypeit_test`` script (via a PR to the `PypeIt Development
   Suite`_) such that the new instrument is included in list of instruments
   tested by executing ``./pypeit_test develop``.

 * The CI tests run by GitHub (see the Checks tab of the PR) on the remote
   repository must pass.

 * You (or someone running the tests on your behalf) must post a successful
   report resulting from your execution of the `Development Suite`_, which
   should look like this:

   .. figure:: ../figures/tests_success.png

        Example posting of successful tests.

   For hotfixes, these tests can be circumvented at the discretion of
   the core developers in the cases where the hotfix is obviously
   correct.

 * All new methods and classes must be at least minimally documented.
   "Minimally documented" means that each method has a docstring that
   gives at least:
   
    #. a one sentence description of the purpose of the method,

    #. a complete list of the required and optional arguments and their meaning,
    
    #. a description of the returned objects, if there are any.
   
   Documentation is expected to adhere to `Sphinx`_ syntax; i.e., the docstrings
   should be `reStructuredText`_.  We accept both `Google-format docstrings`_
   and `Numpy-format docstrings`_, but the former is preferred.

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

 * At biweekly ``PypeIt`` telecons or over the ``PypeIt`` developers Slack, the
   core development team will decide to merge the ``develop`` branch into
   ``release``.

 * A branch is created off of ``develop`` (typically called ``staged``)
   and then a `PR <https://github.com/pypeit/PypeIt/compare>`_ is issued
   to merge ``staged`` into ``release``.  This merge must meet the same
   `Pull Request Acceptance Requirements`_ when merging new branches
   into ``develop``.  However, typically the tests do not need to be run because
   ``develop`` has already passed these tests and ``staged`` will not make any
   substantive changes to the code.

 * **Before being merged into release**, the code is tagged as follows:

    .. code-block:: bash

        cd $PYPEIT_DIR

        # Edit CHANGES.rst to reflect the ending of a development phase
        vi CHANGES.rst
        git add -u
        git commit -m 'tag CHANGES'

        # Create a tag of the form X.Y.Z (using 1.12.1 here as an example).
        # The current autogenerated version is found in pypeit/version.py.
        git tag 1.12.1

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

*Last Modified: 28 Feb 2023*

----


Additional Developer Links
--------------------------

Here are some developer-specific docs:

.. toctree::
   :maxdepth: 1

   metadata
   new_script
   new_spectrograph
   reports
   conda_forge
   build_archived_sensfuncs
   fluxing


