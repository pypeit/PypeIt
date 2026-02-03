
.. include:: include/links.rst

.. _installing:

============
Installation
============

This page provides detailed instructions for installing PypeIt.  For additional
troubleshooting, please consult the PypeIt :ref:`community` and/or `submit an
issue <https://github.com/pypeit/PypeIt/issues>`__ on GitHub.

.. contents:: Table of Contents
    :depth: 1
    :local:

----

.. _user:

User Installation
=================

.. _environment:

Setup a clean python environment
--------------------------------

PypeIt is available from the `Python Package Index <https://pypi.org/project/pypeit/>`_
(PyPI) and is installed via ``pip``.  This process also installs and/or upgrades
PypeIt's :ref:`dependencies`, and for this reason you should always
first set up a clean python environment in which to install PypeIt.  This mitigates
any possible dependency conflicts with other packages you use.

You can set up a new python environment using either `conda`_:

.. code-block:: console

    conda create -n pypeit python=3.11
    conda activate pypeit

or `venv`_:

.. code-block:: console

    python -m venv pypeit
    source pypeit/bin/activate

See the `Managing Environments with Conda
<https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`__
and/or the `venv documentation <https://docs.python.org/3/library/venv.html>`__ for
more details.  The `conda`_ installation method described below creates an
environment for you.

.. _installing-pip:

Install via ``pip``
-------------------

To install the latest release of PypeIt and its required dependencies, within your
virtual environment execute:

.. code-block:: console

    pip install pypeit

.. _optional-dependencies:

Optional Dependencies
^^^^^^^^^^^^^^^^^^^^^

PypeIt has a few optional dependencies that improve and/or expand functionality.

    - If you are generating datacubes (and performing an astrometric
      correction), you will also need the `scikit-image`_ package. It can be
      installed by including it in the optional dependencies, *e.g.*:

      .. code-block:: console

        pip install "pypeit[scikit-image]"

    - To take advantage of an interface that allows you to ingest PypeIt outputs
      into its ``Spectrum1D`` and ``SpectrumList`` objects (see
      :ref:`spec-1d-output`), you can include `specutils`_ in the installation
      like so:

      .. code-block:: console

        pip install "pypeit[specutils]"

.. note::

    Whether or not it is correct syntax to use the quotes in the commands above
    depends on your shell.  The above commands are specific to ZShell, whereas
    you don't need the quotes in Bash.  But, in any case, you should avoid
    copying these commands from your browser since the unicode for quotation
    marks may not be correct, leading to errors when they are directly pasted
    into a terminal window.

Development or pre-release versions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The developers may occasionally ask users to jump on specific branches of the
code to help them more quickly have their issue addressed.  **This is done on a
shared-risk basis!**  I.e., development versions are intrinsically more unstable
and, by definition, transient.

If you find yourself working in branches frequently, you should instead perform
a :ref:`developer_install`.

To install development versions of PypeIt directly from `GitHub
<https://github.com/pypeit/PypeIt>`_, use ``pip``.  First, **as always**, make
sure you are working in a python environment (this also works within a conda
environment).  Then run:

.. code-block:: console

    pip install --upgrade "pypeit[dev]@git+https://github.com/pypeit/PypeIt.git"

This will install the default branch, ``release``.  To install, e.g.,
the ``develop`` branch, run:

.. code-block:: console

    pip install --upgrade "pypeit[dev]@git+https://github.com/pypeit/PypeIt.git@develop"

To install a different branch, use the branch name instead of ``develop`` in the
command above.  Commit hashes, tag names, or git refs can also be specified; see
the `VCS Support documentation
<https://pip.pypa.io/en/stable/cli/pip_install/>`_ for details and examples.

Note that these installations include all the development dependencies (as
indicated by the ``pypeit[dev]`` part of the command).

Install via ``conda``
---------------------

We provide a YAML file that can be used to setup a virtual environment using
`conda`_.  This file creates a conda environment called ``pypeit``, and then
installs ``pypeit``, all of the required dependencies, and the optional
dependency ``specutils`` via ``pip``.  To use this:

    #. Download `environment.yml
       <https://raw.githubusercontent.com/pypeit/PypeIt/release/environment.yml>`__.

    #. Create the conda environment and install ``pypeit`` into it:

        .. code-block:: console

            conda env create -f environment.yml

    #. Activate it:

        .. code-block:: console

            conda activate pypeit

    #. Verify that the new environment was installed correctly and contains ``pypeit``:

        .. code-block:: console

            conda list

        Most of the packages listed will show as coming from the ``pypi`` channel, the rest
        from ``conda-forge``.

This environment should now be ready to use, and it will contain the latest
official ``pypeit`` release.

.. _upgrade:

Upgrading to a new version
--------------------------

Since either installation method above ultimately uses ``pip`` to install
PypeIt, upgrading the package should simply be a matter of executing:

.. code-block:: console

    pip install pypeit --upgrade

If this causes problems (*e.g.*, a new PypeIt script is unavailable or you
encounter script errors), first try uninstalling (``pip uninstall pypeit``) and
then reinstalling.  Also note that not all PypeIt versions are
backwards-compatible; see :ref:`versioning`.

.. _m1_macs:

Installation on Apple Silicon-based Macs
----------------------------------------

Both the `pip`_ and `conda`_ installation methods should be successful for Macs
that use Apple Silicon processors.  The full Anaconda installers also now
include support for Apple Silicon.  Solutions/Recommendations/Feedback for these
installation options are welcome; please `Submit an issue`_.

.. _install_windows:

Installation on Windows
-----------------------

Generally speaking, we encounter most installation issues for Windows users.

An alternative for running under Windows is to install the `Windows Subsystem
for Linux (WSL) <https://learn.microsoft.com/en-us/windows/wsl/install>`_.  This
in effect allows you to run PypeIt on a Windows computer as if it was running
under Linux.

For a pure Windows installation, our recommended installation procedure is as
follows.  Please post questions on our Users Slack if you have difficulties!

#. Download `Python for Windows <https://www.python.org/downloads/windows/>`_.

#. Run the installer.

    * Make sure "Add python.exe to Path" or "Add Python to environment
      variables" is selected before installing.

    * If you have Admin privileges click "Disable path length limit" after the
      installation succeeds.

#. Create a virtual environment as in `Setup a clean python environment
   <environment_>`__ and install PypeIt as described above.

If running ``python`` on Windows brings up a window for the Microsoft Store you
may want to change the application alias.  This is under ``Settings -> Apps ->
App execution aliases`` on Windows 10 and ``Settings -> Apps -> Advanced app
settings -> App execution aliases`` on Windows 11. Disable the ``App Installer``
options for the ``python.exe`` and ``python3.exe`` executables.

----

.. _versioning:

Versioning
==========

Traditionally, we provided no guarantee that *any* PypeIt version was backwards
compatible.  However, as of version 2.0.0, PypeIt uses `Semantic Versioning
<https://packaging.python.org/en/latest/discussions/versioning/>`__.  This
approach uses three version categories --- *major*, *minor*, and *patch* ---
where releases that increment the *major* version number are *not* backwards
compatible.  We expect the most common reason for incrementing the *major*
version number will be because of a backwards-incompatible change to either the
input configuration files (like the ``*.pypeit`` file) or the data models of the
primary output products.

.. important::

    When possible, we **always** recommend you use the most recent version of
    PypeIt and reprocess data as necessary.  The code is always improving, not
    just in functionality but also in robustness of data reduction and
    processing.

Beyond this, we emphasize two important considerations regarding PypeIt versioning:

- **Backwards-incompatible changes to datamodels can break simple viewing
  scripts.**  For example, ``pypeit_show_2dspec`` may fault when trying to view
  ``spec2d*`` files produced using a version of PypeIt that is not backwards
  compatible with your current version.  You can always maintain multiple python
  environments with different PypeIt versions installed or reprocess data with
  your currently installed PypeIt version.

- **Cached files are version-specific.**  Every time you upgrade PypeIt, we
  recommend deleting your existing cache and starting fresh!  See
  :ref:`view-cache`.  The only caveat to this is if you are actively using
  multiple versions of PypeIt (in different environments), meaning you will
  still be using old versions of the cached files.  Otherwise, you will end up
  with multiple versions of the same file on disk.  **Importantly**, the code
  also considers local files you have installed (using, e.g.,
  ``pypeit_install_linelist``) to be version specific.  If you have installed
  such files, you will need to re-install them *after* upgrading.

  If you have locally installed files, your upgrade may look something like
  this:

  .. code-block:: console

    # Check the cache contents
    pypeit_clean_cache -l
    # Delete everything
    pypeit_clean_cache --clear
    # Upgrade pypeit
    pip install pypeit --upgrade
    # Reinstall your local line lists
    pypeit_install_linelist /path/to/my/linelists/*_lines.dat

  .. note::

    If you find particular data files useful for your reductions, please
    consider issuing a PR to include them in the PypeIt repository.  This helps
    the community, and it means you'll avoid these upgrading complications.

----

.. _data_installation:

Additional Data and the PypeIt Cache
====================================

To limit the disk-space required for installation, most of PypeIt's static data
files are either not kept in the GitHub repository or distributed via `pip`_.
PypeIt uses the generalized cache system `provided by Astropy
<https://docs.astropy.org/en/stable/utils/data.html>`__ to interface with the
remote data, which maintains copies of the data files in a user-writeable
location that is independent of the PypeIt installation.  For most users, this
will be ``~/.pypeit/cache``, but the exact location can be set directly using
Astropy's `configuration system
<https://docs.astropy.org/en/stable/config/index.html#astropy-config>`__.  By
default, PypeIt will download necessary files at runtime if they are not already
cached.  Regardless of their location, remote or local, PypeIt essentially
organizes all its reference data into subdirectories of the ``pypeit/data``
directory in your package installation.  The following table gives the reference
name, subdirectory, and remote host for data in this directory tree:

.. _data_dir:

.. include:: include/data_dir.rst

Although most cached files are hosted on GitHub, a few particularly large files
are hosted on Amazon S3 cloud storage.  Note that a host of ``...`` means that
the files should be distributed with your package installation for these
directories.

As stated above, PypeIt will download remote files and store them in your cache
as they're needed to reduce your data.  I.e., you should mostly be able to
ignore the fact that the relevant files are remote, as long as you're running
the reductions while connected to the internet.  However, if you're preparing to
run a set of reductions and you would prefer to pre-load data that you expect to
need, we provide a few specific scripts for interacting with the cache, as
described below.

.. _view-cache:

Viewing/Removing Files in the Cache
-----------------------------------

The ``pypeit_clean_cache`` script allows you to view and/or delete files in the
cache.  To list the cache contents, use the ``-l`` option:

.. code-block:: console

    % pypeit_clean_cache -l
           HOST               BRANCH               SUBDIR FILE
         github               1.15.1                tests gemini_gnirs_32_1_spat_fit.npz
         github               1.15.1            sensfuncs keck_deimos_600ZD_sensfunc.fits
       s3_cloud                  ...   telluric/atm_grids TellPCA_3000_26000_R10000.fits
         github               1.15.1                tests solution_arrays.npz

Note that the files hosted on GitHub will be specific to a branch or version of
PypeIt.  **Every time you upgrade PypeIt, we recommend deleting your existing
cache and starting fresh!**

**Local files** that have been installed into the cache (e.g., using
``pypeit_install_linelist``) will appear as being hosted on GitHub and be
specific to the version of the code used to install it.  When you install local
files, keep two things in mind:

#. The current cache system *does not* keep track of the original on-disk
   location of these files.  When you install these local files into the cache,
   the original file will remain (as long as you don't move/delete it yourself),
   and they will not be removed by ``pypeit_clean_cache``.

#. However, as far as the cache is concerned, these files are specific to a
   given PypeIt version.  This means **you'll need to re-install them** when you
   upgrade PypeIt; otherwise, PypeIt will not recognize their existence in the
   cache.  We discuss upgrading :ref:`above<upgrade>`.

Some example uses for removing files include:

 - To remove your entire cache: ``pypeit_clean_cache --remove_all``.

 - To remove cached files for a specific version: ``pypeit_clean_cache -p 1.15.0``

 - To remove a specific file: ``pypeit_clean_cache -p gemini_gnirs_32_1_spat_fit.npz``

Pre-loading Cache Data
----------------------

Because a fresh install of PypeIt does not contain all of the ancillary data that
might be required for data reduction, users planning to run the pipeline without an
internet connection will need to cache the necessary data files ahead of time.  To ease
this process, we provide the ``pypeit_cache_github_data`` script.  For example, to
download the needed files for the ``keck_deimos`` spectrograph, you would execute:

.. code-block:: console

    $ pypeit_cache_github_data keck_deimos

(Alternatively, you can get all of the cached files hosted on GitHub by
performing a developer installation, if you prefer).  Once cached, the data will
be accessed by PypeIt without requiring an internet connection.  By default,
this script also downloads any files it finds that are *not* specific to a given
spectrograph.  (Unlike previous versions) This script does *not* download any
files hosted in ``s3_cloud`` (see the table :ref:`above<data_dir>`); instead,
use the scripts below.

.. note::

    Beware of rate limits imposed by GitHub.  If you run into this, try setting
    up an access token and export it as the ``GITHUB_TOKEN`` environmental
    variable; see
    `here <https://docs.github.com/en/rest/using-the-rest-api/rate-limits-for-the-rest-api>`__.

.. _install_atmosphere:

Atmospheric Models
------------------

Calculation of the sensitivity functions and general fitting of telluric
absorption uses a PCA decomposition of a massive grid of model atmosphere
spectra across many different observatories, with reference files that are only
a few MB in size. Earlier PypeIt versions used pre-computed model grids which
were much larger (several GB) and observatory-specific; note that the new PCA
models are explicitly designed to be observatory-agnostic, and typically only
differ in their intrinsic spectral resolution (except for experimental
Keck/HIRES and Keck/NIRSPEC [high-resolution mode] model files) to improve
computational efficiency.

.. note::

    Instruments that anticipate needing a telluric correction have the filename
    of a PCA model with sufficient resolution already included in the ``telgridfile``
    `TelluricPar keyword
    <https://pypeit.readthedocs.io/en/latest/pypeit_par.html#telluricpar-keywords>`__.
    While the ``pypeit_install_telluric`` script can be used to download the files
    ahead of time, they are light enough that this is no longer necessary.

It is still possible to use the legacy atmospheric grid files by modifying the
``teltype`` `TelluricPar keyword
<https://pypeit.readthedocs.io/en/latest/pypeit_par.html#telluricpar-keywords>`__.
Only a handful of the original observatory-specific grids are still available for
download on the `PypeIt dev-suite Google Drive`_. Users may select any of the files
in the Google Drive for their telluric correction, download them separately, then
install them using the ``--local`` option to ``pypeit_install_telluric``. For access
to other legacy grids please leave a message in the telluric channel on the
`PypeIt Users Slack <https://pypeit-users.slack.com>`__.

User-provided atmospheric extinction files and wavelength-calibration line lists
--------------------------------------------------------------------------------

As needed to improve their data reduction, users can "install" their own
atmospheric extinction files and/or wavelength-calibration line lists.  PypeIt
manages these *local* files within its cache system.  To install such files, use
the ``pypeit_install_extinctfile`` or ``pypeit_install_linelist`` script,
respectively; see :ref:`install_scripts`, :ref:`extinct-file`, and
:ref:`user_linelists`.

If you find specific files are generally useful/important to your data
reduction, we encourage you to submit a GitHub pull-request so that these files
can be included in the PypeIt repository.

.. important::

    Because PypeIt uses the cache system to manage the local files, it will
    associate each file with the version of the code used to install it in the
    cache.  Every time you upgrade your PypeIt version, you should delete the
    local files from the cache (this will not remove the local file itself) and
    re-install them using the upgraded version of PypeIt.  See :ref:`view-cache`
    and :ref:`upgrade`.

Quick-Look Calibration Files
----------------------------

.. note::

    We continue to work on cleaner installation solutions for these data
    products.  In the meantime, note that you will likely need to re-run the
    data-specific installation scripts described below every time you upgrade
    your installation.

Some of the quick-look reductions provided by PypeIt require canned calibration
files to speed up the data-reduction process, as appropriate for a quick-look
result.  These files are hosted in the ``QL_CALIB`` directory in the `PypeIt
dev-suite Google Drive`_.

To install the quick-look calibration files:

    #. Right-click on the ``QL_CALIB`` folder in the `PypeIt dev-suite
       Google Drive`_ and select the "Download" option from the drop-down menu.
       This will download a zip file containing the full directory contents.
       Its current size (as of 22 July 2021) is about 35 MB.

    #. Run the ``pypeit_install_ql_calibs`` script.  E.g.:

        .. code-block:: console

            $ pypeit_install_ql_calibs --zip ~/Downloads/QL_CALIB-20210722T162355Z-001.zip --odir my_path

The ``pypeit_install_ql_calibs`` script will unzip the downloaded file in the
``my_path`` directory and create a symlink to the extracted directory in the
``pypeit/data/`` directory of your PypeIt installation.  The script can
automatically delete the zip file using the ``--rmzip`` option.  If you already
have the ``QL_CALIB`` directory, you can also use the script to simply create
the symlink using the ``--ql_path`` option.

.. warning::

    The installation script simply creates symlinks to the downloaded data.
    This means that if you move the original data, the symlinks will become
    broken **and you will need to rerun the installation script.**

.. _devsuite-raw-data:

Raw Data
--------

Example raw data for all supported spectrographs are used in extensive testing
of the code base during development; see :ref:`dev-suite`.  General users should
not need access to these data; however, they may be useful for learning how to
use PypeIt before running it on your own data from the same instrument.
These data are stored in the ``RAW_DATA`` directory in the `PypeIt dev-suite
Google Drive`_, divided into subdirectories for each instrument and instrument
setup.  See also the `PypeIt Development Suite`_ GitHub repository, which
includes a :doc:`pypeit_file` for each instrument and setup used during
development testing.

----

.. _notes:

Important Package Notes
=======================

.. _interactive:

Interactive Tools
-----------------

Interactive tools in PypeIt are generally built using the `QT
<https://www.qt.io/>`_ windowing toolkit. The ``qtpy`` package is used to
provide an abstract interface to the two most widely used QT bindings for
Python (see :ref:`dependencies`):

* `pyqt <https://riverbankcomputing.com/software/pyqt/intro>`_
* `PySide <https://wiki.qt.io/Qt_for_Python>`_

At least one of those bindings must be installed for the interactive GUIs to
work. By default PypeIt will install ``pyqt6``. Other backends can be used
by installing them manually via ``pip`` or ``conda`` and then setting the ``QT_API``
environment variable. See the `QtPy documentation <https://github.com/spyder-ide/qtpy>`_
for more details.

ginga Plugins
-------------

PypeIt requires the ``ginga`` viewer and uses at least one ``ginga`` plugin
to enable specific display functionality. No special considerations are needed
to have these plugins installed; however, you can check that they're enabled by
running the following script with the following result::

    $ pypeit_chk_plugins
    [INFO]    :: All required plugins found: SlitWavelength, Spec1dView

If the check is unsuccessful, you will see an error message listing
the missing plugins. If you have a problem, please `submit an issue
<https://github.com/pypeit/PypeIt/issues>`__.

----

.. _dependencies:

Package Dependencies
====================

All PypeIt dependencies are installed along with the installation of
PypeIt itself.  Beware this means that packages in your current environment
may be updated depending on the PypeIt version requirements, which is why we
recommend you :ref:`environment` for PypeIt!  The current version
requirements for both users and developers are:

.. include:: include/dependencies_table.rst

Dependency Caveats
------------------

Some users have run into the following complications when installing the
PypeIt dependencies.  If you run into any more, please `submit an issue
<https://github.com/pypeit/PypeIt/issues>`__.

.. TODO: IS THIS FIRST ITEM STILL TRUE?

- At the moment, an implicit dependency on QT bindings remains because of our
  dependence on ``linetools``.

----

.. _test_installation:

Test Your Installation
======================

Tagged versions of PypeIt are extensively tested before distribution.
However, it is worth testing that your installation has been successful, as
follows.

The most basic tests that PypeIt has been properly installed is to get the
help dialog for one of its main executables.  I.e., from a terminal widow, type:

.. code-block:: console

    run_pypeit -h

A second basic test is to try to import PypeIt from within a python session.
For example:

.. code-block:: console

    python
    >>> import pypeit

**To ensure that your installation of ``pyqt6`` works**,
you can try to use ``pypeit_show_1dspec`` on one of the test files distributed
with the package.  Below is a zshell command-line incantation (it's likely the
same in bash) that will locate a test spec1D file and attempt to use
:ref:`pypeit_show_1dspec` to show it:

.. code-block:: console

    python -c "from pypeit import dataPaths; print(dataPaths.tests.get_file_path('spec1d_r153-J0025-0312_KASTr_20150123T025323.850.fits'))" | xargs -I {} pypeit_show_1dspec {}

If ``pyqt6`` or another Qt backend is correctly installed, this should show a test
spectrum from the Shane/KAST spectrograph at Lick Observatory.

----

.. _install_troubleshoot:

Troubleshooting
===============

If you have trouble installing PypeIt, you're encouraged to join our `PypeIt
Users Slack <https://pypeit-users.slack.com>`__ using `this invitation
link <invite_>`_ and post your issue to the ``#installing`` channel.  Below is an
incomplete list of issues that users have reported in the past.  In addition to
posting to the Users Slack if your issue isn't among those listed below, *please
let us know if these suggestions do not work for you.*

----

**I am trying to install PypeIt for the first time and it fails!**:  The root
problem of this can be system dependent:

 - First, *always* make sure you install the code into a fresh environment.

 - If you're on Windows, make sure you follow the :ref:`install_windows`
   instructions.

 - Occasionally, the installation may fail because of incompatible dependencies.
   This may be because of recent releases of one of PypeIt's dependencies; i.e.,
   updates to packages since the most recent PypeIt release.  **Please let us
   know if this happens!**  Once notified, we will try to issue a new release
   asap that corrects the incompatibility.  In the short-term, we may ask you to
   install old versions of packages that we know work.

----

**I am trying to upgrade PypeIt and it fails!**:  First try uninstalling your
current PypeIt version:

.. code-block:: bash

    pip uninstall pypeit

Then reinstall it.  If that also fails, try creating a fresh environment and
reinstalling PypeIt in that new environment.

----

**The installation process succeeded, but the code is faulting!**:  This could
be for a few reasons:

 - Recall that PypeIt isn't necessarily backwards compatible.  If you've
   upgraded PypeIt and tried to use it with data that was reduced by a previous
   version, the fault may because of changes between versions.  You will either
   need to revert to your previous version or reprocess the data.

 - This may be because of dependency changes.  A tell-tale signature of this is
   if you get errors associate with missing or unknown keywords or arguments.
   This is may be because of recent releases of one of PypeIt's dependencies;
   i.e., updates to packages since the most recent PypeIt release.  **Please let
   us know if this happens!**  Once notified, we will try to issue a new release
   asap that corrects the incompatibility.  In the short-term, we may ask you to
   install old versions of packages that we know work.

----

**The installation process succeeded and the code completes without faulting,
but the output looks wrong!**:  This could happen for any number of reasons.
*We always welcome reports of failures!*  Either `submit an issue
<https://github.com/pypeit/PypeIt/issues>`__ or report it on the PypeIt Users
Slack.  However, here are a few things to note and/or try:

 - Make sure you have checked your calibrations; see :ref:`calibrations`.  The
   issue may be related to a parameter that you can change.

 - If you don't see any ``spec1d`` files in your ``Science`` folder, this is
   likely because the code didn't find any objects; see :ref:`object_finding`.

 - If you've recently upgraded the code, this may be related to changes in
   dependencies that the developers didn't catch.  PypeIt performs *a lot* of
   testing before issuing a new release, but does not have complete test
   coverage and performance validation.  This means silent failures are the most
   difficult to catch.

 - And, of course, the code will have bugs.  If you find one, the more
   information you provide the developers, the easier it will be for us to track
   down the issue.  Valuable information includes your OS, OS version, python
   version, PypeIt version, and the full Traceback provided with the error.  QA
   plots and ``ginga`` screen grabs that illustrate the issue are also very
   helpful!
