.. include:: ../include/links.rst

.. _development:

=====================================
Development Procedures and Guidelines
=====================================

.. contents:: Contents
    :depth: 1
    :local:

----

We welcome and encourage community development of PypeIt!  All code contributors
are expected to follow our :ref:`codeconduct` and these development guidelines.

.. _developer_install:

Installing and working with your fork
=====================================

If you plan to develop the PypeIt code, you should install the code from a fork
of the main respository, as described below.

You may also need/want to create a fork of the :ref:`dev-suite`!  This process
proceeds similarly to the main code repository, **except that you do not need to
pip install the dev-suite repository.**

Setup a clean python environment
--------------------------------

This is the same as described for a user installation; see :ref:`environment`.

Fork the code
-------------

Except for a few maintainers, code development should be done in forks of the
code base.  To fork the code, go to `GitHub
<https://github.com/pypeit/PypeIt>`__ and click the "Fork" button at the upper
right.  Do not change the name of the repository (``PypeIt``).

When you fork the code, you will be asked if you want to only copy the
``release`` branch (default).  If you uncheck this box, you will get *all* of
the current branches, which is almost certainly not something you want to do.
See :ref:`below <checkout_remote_branch>` for one way to add the ``develop``
branch to your repo.

Finally, you are **strongly** encouraged to add branch protection rules for both
the ``release`` and ``develop`` branches in your fork (these rules are *not*
inherited from the main repository).

Follow a similar process for creating a fork of the :ref:`dev-suite`.

Local Install
-------------

To install from source, first create a local clone of your fork (replace
``{github_handle}`` with the GitHub handle used to create your fork):

.. code-block:: console

    git clone https://github.com/{github_handle}/PypeIt.git

This will create a ``PypeIt`` directory at the location the command above was
executed.

To install the dev-suite, perform a similar operation (from your fork):

.. code-block:: console

    git clone https://github.com/{github_handle}/PypeIt-development-suite.git

Install the main PypeIt code, including its package dependencies, using ``pip``
(even if you're in a conda environment):

.. code-block:: console

    cd PypeIt
    pip install -e ".[dev]"

This (specifically the ``-e`` option) creates an "editable" installation, which
means that any changes you make in the repository directory tree will take
immediate effect the next time the code is imported. Including the ``[dev]`` set
of optional dependencies ensures that all of the tools you need to test and
build PypeIt are installed. (Note that you may or may not need the quotes above
depending on your shell, and that you should avoid cutting and pasting these
commands into a terminal window.)

You *do not* need to execute a similar installation for the dev-suite.

Connect the main repository to your fork
----------------------------------------

To ensure you are up-to-date with the main repository, you should add it as a
remote:

.. code-block:: console

    % git remote -v
    origin	https://github.com/{github_handle}/PypeIt.git (fetch)
    origin	https://github.com/{github_handle}/PypeIt.git (push)
    % git remote add upstream https://github.com/pypeit/PypeIt.git
    % git remote -v
    origin	https://github.com/{github_handle}/PypeIt.git (fetch)
    origin	https://github.com/{github_handle}/PypeIt.git (push)
    upstream	https://github.com/pypeit/PypeIt.git (fetch)
    upstream	https://github.com/pypeit/PypeIt.git (push)

For the dev-suite:

.. code-block:: console

    % git remote -v
    origin	https://github.com/{github_handle}/PypeIt-development-suite.git (fetch)
    origin	https://github.com/{github_handle}/PypeIt-development-suite.git (push)
    % git remote add upstream https://github.com/pypeit/PypeIt-development-suite.git
    % git remote -v
    origin	https://github.com/{github_handle}/PypeIt-development-suite.git (fetch)
    origin	https://github.com/{github_handle}/PypeIt-development-suite.git (push)
    upstream	https://github.com/pypeit/PypeIt-development-suite.git (fetch)
    upstream	https://github.com/pypeit/PypeIt-development-suite.git (push)

Updating your fork
------------------

To update your fork so that it will recognize changes in the main repository,
run:

.. code-block:: console

    % git fetch upstream

Do this often!  Note, however, that this **does not** update the content of your
fork to match the changes in the upstream repo.

.. _checkout_remote_branch:

Checkout a remote branch and add it to your fork
------------------------------------------------

The most obvious application of this is for your first checkout of the
``develop`` branch into your fork.  However, you can do this with any branch
(and with any remote).  The steps are (1) update your repo, (2) checkout the
upstream branch, (3) push it to your fork (and change its remote tracking):

.. code-block:: console

    % git fetch upstream
    % git checkout -b develop upstream/develop
    % git push -u origin develop

Merging changes made to an upstream branch
------------------------------------------

The ``fetch`` command does *not* merge changes made to the upstream repo.  You
have to do that explicitly using either merge:

.. code-block:: console

    % git fetch upstream
    % git checkout develop
    % git merge upstream/develop
    % git push

or rebase:

.. code-block:: console

    % git fetch upstream
    % git checkout feature_branch
    % git rebase upstream/develop
    % git push --force

In the rebase example note two things:

- This does not update your local ``develop`` branch because we were not on our
  local ``develop`` branch when we executed the rebase command.  This *only*
  updates the ``feature_branch`` with the changes made to its original base
  branch, ``develop`` in this case.

- Rebasing often requires forcing the push to your fork, using the ``--force``
  option, because it alters the commit history.  If you try pushing without the
  force argument, you will likely get a warning that your branch is out of date
  and you need to pull the origin branch.  **Do not do this.**  Pulling in your
  remote branch after rebasing will lead to a number of conflicts and, likely,
  lost work.

Collaborating on a branch
-------------------------

If you're working with another developer (who has a different fork) on a single
branch, you can add their fork as another remote and proceed similarly to how
you work with the main repo:

.. code-block:: console

    % git remote add collab https://github.com/{collaborator}/PypeIt.git
    % git remote -v
    origin	https://github.com/{github_handle}/PypeIt.git (fetch)
    origin	https://github.com/{github_handle}/PypeIt.git (push)
    upstream	https://github.com/pypeit/PypeIt.git (fetch)
    upstream	https://github.com/pypeit/PypeIt.git (push)
    collab	https://github.com/{collaborator}/PypeIt.git (fetch)
    collab	https://github.com/{collaborator}/PypeIt.git (push)

Then you can both work on a ``feature_branch`` while keeping your branches in sync:

.. code-block:: console

    % git fetch collab
    % git checkout feature_branch
    % git merge collab/feature_branch
    % git push

When you're done, one of you can issue the PR from one of your forks into the
main repo.  If you want to keep your remote list clean, you can remove the
remote connection after finishing the collaborative work:

.. code-block:: console

    % git remote remove collab

.. _persistent-branches:

Persistent and development branches
===================================

PypeIt maintains two persistent branches:

 * ``release``: This is the primary stable version of the code.  Modulo any
   very recent :ref:`hotfixes`, this is closest to the most recently tagged and
   released version.  Pull requests to this branch are only done before tagging
   a new release of the code or to perform critical bug hotfixes.  The release
   schedule is discussed during our bi-weekly development meetings.

 * ``develop``: This is the main development version of the code.  It should be
   stable enough to use, but it may contain experimental, unsupported code that
   is work in progress.

When editing the code, please create a new branch stemming from the ``develop``
branch.  You should also pull and merge in the most recent version of the
``release`` branch to make sure your new branch includes any very recent
hotfixes.  On the command line, you can do this as follows from your ``PypeIt``
directory:

.. code-block:: bash

    # Fetch upstream changes
    git fetch upstream
    # Switch to your local release branch
    git checkout release
    # Update it
    git merge upstream/release
    # Push the updated branch to your fork
    git push
    # Switch to your local develop branch
    git checkout develop
    # Update it
    git merge upstream/develop
    # Push the updated branch to your fork
    git push
    # Create a new branch from your local develop branch
    git checkout -b my_new_feature
    # Make sure develop is up-to-date with release
    git merge release
    # Push the new branch to your fork
    git push -u origin my_new_feature

.. note::

    In terms of the merge with the release branch, beware that you may need to
    start a new release version doc that reflects the jump in the version
    number.  This should only be necessary if your branch is the first one after
    a new tag is released.  See :ref:`changelog`.

Development Principles and Communication
========================================

The main thing to keep in mind when developing for PypeIt is that its
primary use is as an end-to-end reduction pipeline.  This has a few
implications:

 * By default, the execution of ``run_pypeit`` should continue either until a
   critical error is raised or the reduction is complete.  **No direct
   interaction with the code should be required at any point**.  PypeIt does
   have some interactive components, but these are executed only if specifically
   requested by command-line arguments or via separate scripts.

 * Any input needed from the user for your feature should be provided by
   :ref:`parameters` (preferred) or as a command-line argument.

 * When developing and debugging, you may need to interact with the code using
   `IPython.embed`_; however, these instances should be removed before
   performing a pull request.

 * The success or failure of any given procedure must be assessed via
   automatically generated quality-assessment figures (preferred) or via
   scripts that interact with the primary output files.

 * See :doc:`here <new_spectrograph>` for guidance on 
   adding a **new spectrograph** to the list of spectroscopic data that
   PypeIt can reduce.

 * If your development includes adding a **new executable script**, see advice
   at :ref:`new_script`.

Feature development in PypeIt is unlikely to be fully independent of
other development activities.  Your feature will likely depend on or
influence the outcome of other modules during the data-reduction
process.  This leads to a few important guidelines:

 * Make sure that your branch is always up-to-date with the ``develop`` *and*
   ``release`` branches.  E.g.:

   .. code-block:: bash

        git fetch upstream
        git checkout release
        git merge upstream/release
        git checkout develop
        git merge upstream/develop
        git checkout my_new_feature
        git merge release
        git merge develop

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


.. _changelog:

Logging changes
---------------

It is important to log changes made to the code in a way that other developers
and eventually users can interpret.  In the past we have done this using the
single ``CHANGES.rst`` file; however, we now have version specific change logs
in the ``doc/releases`` directory.  In terms of development guidelines:

- Changes made to the code should be logged in the relevant development log.
  For example, all changes made *after* version 1.14.0 will be logged in a
  ``doc/release/1.14.1dev.rst`` file.  If the relevant file doesn't exist when
  you submit your PR, create it.

- Changes are expected to fall under a small set of broad categories, like
  improvements to performance for specific instruments, minor bug fixes, or
  datamodel changes (see previous release docs for examples).  When including
  your change, add it below the relevant heading; if no relevant heading
  exists, add a new one.

- Hotfixes merged directly to the ``release`` branch should *also be added to
  the relevant development log*.  I.e., these changes are not part of the
  released tag, even if they are in the "release" branch.  Again, if the
  relevant file doesn't exist when you perform the hotfix, create it in a way
  that it will get merged with the identical doc in the ``develop`` branch.  See
  :ref:`hotfixes`.

- When tagging, the development log will be renamed to the new tag version, and
  a new file should be created for the next development phase.  See
  :ref:`tagging`.

.. _hotfixes:

Hotfixes
--------

There may be bugs in the ``release`` version of the code that are not caught by
the tests, but significantly impact some users.  Fixing these issues leads to a
patch release of the code, following this procedure:

 * Checkout the ``release`` version of the code.

 * Create a new branch from the ``release`` version (*not* the ``develop``
   version).

 * Implement the hotfix.

 * Create a release doc that lists the hotfix.  For example, if the current
   development doc is ``doc/release/1.14.1dev.rst``: *create* the new file
   ``doc/release/1.14.1.rst``, edit it to include a description of the hotfix,
   add it to the repo, and rename the existing development doc to indicate the
   version increment (``doc/release/1.14.2dev.rst``).  Also update the
   ``doc/whatsnew.rst`` file to include the new release files.

 * Issue a PR and follow the :ref:`tagging`.  Following the example above, the
   new tag would be ``1.14.1``.

Testing
=======

PypeIt performs extensive testing using the :ref:`dev-suite`; follow that link
for more details on executing the tests.  Below, we describe how to add new
tests.

.. _dev-suite-tests:

Development Suite
-----------------

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

    #. Edit the ``test_scripts/test_setups.py`` file in the `PypeIt Development
       Suite`_ to include the new setup among the tests to perform.  Follow the
       instructions at the top of that file.

    #. Run the full development test suite to completion. Once all tests pass,
       the ``test_priority_file`` will be updated with the new test. This file
       tells the test scripts what order to run the tests in for optimum CPU
       utilization.  Commit ``test_priority_list`` and any other files added to
       the dev-suite repository and submit a pull request.

The :ref:`dev-suite` also contains unit tests that require use of data in the
``RAW_DATA`` directory and "vet" tests that are set of unit tests that require
output files from PypeIt scripts.  The former typically test simple
functionality of the PypeIt code, whereas the latter (vet tests) check the
results of the PypeIt scripts against the expected performance/result.

.. _unit-tests:

Unit Tests (GitHub CI)
----------------------

Unit tests performed by GitHub continuous integration (CI) are located in the
``$PYPEIT_DIR/pypeit/tests`` directory.  To run them, make sure you have
`pytest`_ installed (this should be true if you followed the developer
installation procedure) and then run, from your ``PypeIt`` directory:

.. code-block:: bash

    pytest

If some tests fail, you can run an individual test, e.g. ``test_wvcalib.py``
with

.. code-block:: bash

    pytest -s pypeit/tests/test_wvcalib.py

Note that the "-s" option allows you to insert interactive debugging commands
into the test, here ``test_wvcalib.py``, to help determine why the test is
failing.

.. warning::

    Running these tests may generate files that should be ignored.  **Please do
    not add these test files to the repository.**  We try to include clean-up as
    part of the tests, but these are not always caught.

Note also that the use of `pytest`_ requires the test dependencies to be
installed. It is also possible, and often preferable, to run tests within their
own isolated environments using `tox <https://tox.readthedocs.io/en/latest/>`_.
This provides the capability to easily run tests against different versions of
the various dependencies, including different python versions. The available
``tox`` environments are defined in the ``tox.ini`` file and can be listed by
running ``tox -a``. To run tests against the default dependencies using the
default python, do:

.. code-block:: bash

    tox -e test

To specify a python version, do something like:

.. code-block:: bash

    tox -e py312-test

To test against, for example, the ``main`` branch for ``astropy`` on GitHub, you
can do:

.. code-block:: bash

    cd $PYPEIT_DIR
    tox -e py312-test-astropydev

Similar ``dev`` dependencies are configured for ``numpy`` and ``ginga``, as
well.

Unit tests included in the main PypeIt repo should *not* require large data
files.  Some files are kept in the repo for this purpose (see the
``pypeit/data/tests`` directory), but they should be minimized to keep the size
of the repository manageable (these test files are *not* included in the package
distribution).  In general, unit tests that require input data files should
instead be added to the :ref:`dev-suite`.

Workflow
========

A typical PypeIt development workflow is as follows:

 * Create a new branch stemming from the ``develop`` branch (:ref:`hotfixes`
   should instead branch from ``release``); see :ref:`persistent-branches`.

 * Develop and debug the feature

 * Run the unit tests, fix any failures, add new tests that test your new
   feature(s), and/or modify the tests to accommodate your new feature.

 * Run the `Development Suite`_ and fix any failures.

   .. warning::

        The :ref:`dev-suite` is *extensive* and takes significant computing
        resources and time.  The PypeIt development team consistently executes
        these tests using cloud computing.  We recommend you ensure that your
        PypeIt branch successfully runs on either a specific instrument of
        interest or ``shane_kast_blue`` first, and then someone on the PypeIt
        development team can execute the tests in the cloud.  From the top-level
        directory of the :ref:`dev-suite`, you can run all tests for
        ``shane_kast_blue`` as follows, from the ``PypeIt-development-suite``
        directory:

        .. code-block:: bash

            ./pypeit_test all -i shane_kast_blue

 * Edit the relevant development log (e.g., ``doc/release/2.1.0dev.rst``) to
   include your key developments (see :ref:`changelog`).
 
 * Update the `documentation`_.  There is a convenience script to rebuild the
   documentation: In the top-level directory of the repo, run

   .. code-block:: bash

        ./update_docs

   Some documentation files are build in an automated way to keep them current
   as the code changes.  To successfully run one of these scripts, 
   ``$PYPEIT_DEV`` must be an environmental variable that points to the
   :ref:`dev-suite` such that code can access the ``RAW_DATA`` directory.  If
   the environmental variable is not defined or the ``RAW_DATA`` directory does
   not exist, the documentation updating script will skip these automated
   updates and issue a warning (in yellow).
   
   *Any* warnings from sphinx (in red or pink) *must* be fixed (i.e.,
   "warnings" should be considered the same as "errors" in this context).  The
   ``update_docs`` script will provide a summary at the end of its execution,
   letting you know if the build was successful or not.  If you're having
   difficulty getting the right sphinx/rst incantation, ping the documentation
   channel in the `PypeIt Developers Slack <https://pypeit.slack.com>`__.  Also
   note that, even if no warnings are issued, it's useful to check that the
   documentation formats as you expect.  After building the docs, you can open
   the ``doc/_build/html/index.html`` file to view and navigate through the
   documentation in its entirety.

   When debugging the docs, you might also find it useful to skip the scripts
   that build some of automated files (those that do not require access to the
   :ref:`dev-suite`).  You can do this by running

   .. code-block:: bash

        cd doc ; make htmlonly ; cd ../

 * Make sure all your edits are committed and pushed to your fork:

   .. code-block:: bash

        git add -u
        git commit -m 'final prep for PR'
        git push

 * `Submit a Pull Request (PR)
   <https://github.com/pypeit/PypeIt/compare>`_ from your fork to the main
   repository. Unless otherwise requested, all PRs should be submitted to the
   ``develop`` branch.

.. note::

   The addition of new commits causes ``setuptools_scm`` to automatically
   increment the version based on the last tag that was pushed.  This will be of
   the form ``{next_version}.dev{distance}+{scm letter}{revision hash}``.  See
   the `setuptools_scm documentation <https://github.com/pypa/setuptools_scm>`_
   for more details.

Pull Request Acceptance Requirements
====================================

Once you've submitted a pull request, two developers will review your PR and
provide comments on the code.  The minimum requirements for acceptance of a PR
are as follows:

 * If your PR introduces a new instrument (see :ref:`new_spec`) that PypeIt
   is to support for the long term, this instrument *must* be added to the
   :ref:`dev-suite`.  That means raw data should be added to the Google Drive
   and relevant tests should be added to the
   ``PypeIt-development-suite/pypeit_test`` script (via a PR to the
   :ref:`dev-suite` repo) such that the new instrument is included in the list
   of instruments tested by the testing script (``pypeit_test``).

 * The CI tests run by GitHub (see the Checks tab of the PR) on the remote
   repository must pass.

 * You (or someone running the tests on your behalf) must post a successful
   report resulting from your execution of the :ref:`dev-suite`, which should
   look something like this:

   .. figure:: ../figures/tests_success.png

        Example posting of successful tests.

   For :ref:`hotfixes`, these tests can be circumvented at the discretion of the
   maintainers in the cases where the hotfix is obviously correct.

 * All new methods and classes must be at least minimally documented.
   "Minimally documented" means that each method has a docstring that
   gives at least:
   
    #. a one sentence description of the purpose of the method,

    #. a complete list of the required and optional arguments and their meaning,
    
    #. a description of the returned objects, if there are any.
   
   Documentation is expected to adhere to `Sphinx`_ syntax; i.e., the docstrings
   should be `reStructuredText`_.  We accept both `Google-format docstrings`_
   and `Numpy-format docstrings`_ (preferred).

 * The docstrings for any changes to existing methods that were altered
   must have been modified so that they are up-to-date and accurate.

 * The documentation must be successfully recompiled.

 * Spurious commented code used for debugging or testing is fine, but
   please let us know if you want it to be kept by adding a relevant
   comment, something like ``# TODO: Keep this around for now``, at the
   beginning of the commented block.  Otherwise, we're likely to remove
   the commented code when we come across it.

 * "Unsupported code," that is code that is experimental and still work
   in progress, should be minimized as much as is reasonable.  The
   relevant code block should be clearly marked as experimental or WIP,
   and it should not be executed by the main PypeIt executable,
   ``run_pypeit``.

 * At least two reviewers must accept the code.

.. _tagging:

Tagging Protocol
================

The core development team will regularly tag "release" versions of the
repository.  Tagging a release version of the code is triggered anytime the
development branch of the code or a hotfix is merged into the ``release``
branch.  **Only a maintainer can tag the code.*** The tagging process is as
follows:

 * At biweekly PypeIt telecons or over the PypeIt developers Slack, the
   core development team will decide to merge the ``develop`` branch into
   ``release``.

 * A branch is created off of ``develop`` (typically called ``staged``) and then
   a `PR <https://github.com/pypeit/PypeIt/compare>`_ is issued to merge
   ``staged`` into ``release``.  This ``release...staged`` PR must meet the same
   `Pull Request Acceptance Requirements`_ when merging new branches into
   ``develop``.  Code review is expected to be limited (because all code changes
   will have been reviewed before pulling into ``develop``), but the result of
   the dev-suite tests must be shown and approved.  The reason for creating the
   new branch, instead of a direct ``release...develop`` PR, is to allow for the
   following updates to ``staged`` before merging (``develop`` is a protected
   branch and cannot be directly edited):

        * Fix any test failures.  As necessary, an accompanying :ref:`dev-suite`
          PR may be issued that includes test fixes required code changes.  If
          no code changes are required, a :ref:`dev-suite` PR should be issued
          that merges its ``develop`` branch directly into its ``main`` branch.

        * Make any final updates to the development log, and rename the log to
          the new tagged version (e.g., move ``1.14.1dev.rst`` to either
          ``1.14.1.rst`` or ``1.15.0.rst``).  The ``doc/whatsnew.rst`` should
          also be updated to reflect the file name change.

        * Update the list of supported versions in the ``SECURITY.md`` file.

        * Update the documentation by executing ``./update_docs``, add any
          updated files, and correct all errors/warnings.

 * Once the ``release`` branch and the :ref:`dev-suite` ``main`` branch are
   updated, the dev-suite tests are re-run using these two branches.  These
   tests must pass before tagging.  Once they pass, the code is tagged as
   follows:
          
    .. code-block:: bash

        # Create a tag of the form X.Y.Z (using 1.14.0 here as an example).
        # The current autogenerated version is found in pypeit/version.py.
        git checkout release
        git pull
        git tag 1.14.0

        # Push the new tag
        git push --tags

   Similarly, a matching tag is executed for the dev-suite code (these tags only
   exist for versions 1.15 and later).

 * The tag of the ``pypeit`` code-base (not the dev-suite) is released for
   `pip`_ installation.

    .. code-block:: bash

        git checkout 1.14.0
        # Make sure you have the most recent version of twine installed
        pip install twine --upgrade
        pip install build --upgrade
        # Construct the pip distribution
        python -m build --sdist --wheel .
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

DOI
---

If we wish to generate a new DOI for the code, it is as simple as

 * Generate a `new release on GitHub
   <https://help.github.com/en/github/administering-a-repository/about-releases>`_.

 * Update the DOI in the ``README.rst``


----

This document was developed and mutually agreed upon by: Kyle Westfall,
J. Xavier Prochaska, Joseph Hennawi.

*Last Modified: 16 Apr 2025*

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


