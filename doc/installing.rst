
.. include:: include/links.rst

.. _installing:

============
Installation
============

.. warning::

    **Apple Silicon users** who are having issues installing PypeIt may need
    to set up an environment configured for x86-64.  See :ref:`m1_macs` for detailed
    steps.  If that also fails, please `Submit an issue`_ and/or reach out to
    our user :ref:`community`.

Below, we provide detailed instructions for installing PypeIt.  For
troubleshooting, please consult the PypeIt :ref:`community` and/or `submit
an issue <https://github.com/pypeit/PypeIt/issues>`__ on GitHub.

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
PypeIt's :ref:`dependencies`, and for this reason, we highly (!!) recommend you
first set up a clean python environment in which to install PypeIt.  This mitigates
any possible dependency conflicts with other packages you use.

You can set up a new python environment using either `conda`_:

.. code-block:: console

    conda create -n pypeit python=3.11
    conda activate pypeit

or `virtualenv`_:

.. code-block:: console

    virtualenv pypeit
    source pypeit/bin/activate

See the `Managing Environments with Conda
<https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_
and/or `Virtualenv documentation <https://virtualenv.pypa.io/en/latest/>`_
for more details. See also `virtualenvwrapper
<https://virtualenvwrapper.readthedocs.io/en/latest/>`_ as an option for more
easily managing `virtualenv`_ environments. The `conda`_ installation method described below
creates an environment for you.

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

Install via ``conda``
---------------------

`conda`_ is a popular and widely-used package and environment manager.  We
provide a YAML file that can be used to set up the virtual environment for
you.  This file creates a conda environment called ``pypeit``, and then
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

        Most of the packages listed will show as coming from the ``pypi`` channel.

This environment should now be ready to use and contain the latest official
``pypeit`` release.

.. _upgrade:

Upgrading to a new version of PypeIt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since either installation method above ultimately uses ``pip`` to install
PypeIt, upgrading the package should simply be a matter of executing:

.. code-block:: console

    pip install pypeit --upgrade

If this causes problems (*e.g.*, a new PypeIt script is unavailable or you
encounter script errors), first try uninstalling (``pip uninstall pypeit``) and
then reinstalling.  There are two important things to keep in mind when
upgrading:

 - **PypeIt datamodels are not necessarily backwards-compatible.**  This means
   that, *e.g.*, ``pypeit_show_2dspec`` may fault when trying to view
   ``spec2d*`` files produced with your existing PypeIt version after upgrading
   to a new version.  **The best approach is to always re-reduce data you're
   still working with anytime you update PypeIt.**

 - As opposed to earlier versions of PypeIt, the cached files are now
   version-specific.  **Every time you upgrade pypeit, we recommend deleting
   your existing cache and starting fresh!**  The only caveat to this is if you
   are actively using multiple versions of PypeIt, meaning you will still be
   using old versions of the cached files.  Otherwise, you will end up with
   multiple versions of the same file on disk.  **Importantly**, the code also
   considers local files you have installed (using, e.g.,
   ``pypeit_install_linelist``) to be version specific.  If you have installed
   such files, you will need to re-install them *after* upgrading.

If you have locally installed files, your upgrade may look something like this:

.. code-block:: console

    pypeit_clean_cache --remove_all
    pip install pypeit --upgrade
    pypeit_install_linelist /path/to/my/linelists/*_lines.dat

.. note::

    If you find particular data files useful for your reductions, please
    consider issuing a PR to include them file in the PypeIt repository.  This
    helps the community, and it means you'll avoid these upgrading
    complications.
    
.. _m1_macs:

User Installation on Apple Silicon-based Macs
---------------------------------------------

Both the `pip`_ and `conda`_ installation methods should be successful for Macs
that uses Apple Silicon processors.  The full Anaconda installers also now
include support for Apple Silicon.

If the above does not work, you may need to set up a virtual environment configured
for x86-64:

    #. Install `miniconda <https://docs.conda.io/en/main/miniconda.html>`_.
    #. Use ``conda`` to create an environment configured for x86-64:

        .. code-block:: console

            conda create -n pypeit -y
            conda activate pypeit
            conda config --env --set subdir osx-64

    #. Install whichever version of Python you want (*e.g.*):

        .. code-block:: console

            conda install python=3.11

    #. Install PypeIt via ``pip`` as above.

Solutions/Recommendations/Feedback for these installation options are welcome;
please `Submit an issue`_.

User Installation on Windows
---------------------------------------------

#. Download `Python for Windows <https://www.python.org/downloads/windows/>`_.

#. Run the installer.
   
    * Make sure "Add python.exe to Path" or "Add Python to environment variables" is selected before installing.
    * If you have Admin privileges click "Disable path length limit" after the installation succeeds.

#. Downloand and run the `Visual Studio build tools <https://visualstudio.microsoft.com/visual-cpp-build-tools/>`_ installer.

    * Only "Desktop Development with C++" needs to be checked.
    * Click install

#. Create a virtual environment as in `Setup a clean python environment <environment_>`__ and install PypeIt as described above.

If running ``python`` on Windows brings up a window for the Microsoft Store you may want to change the application alias.
This is under ``Settings -> Apps -> App execution aliases`` on Windows 10 and ``Settings -> Apps -> Advanced app settings -> App execution aliases``
on Windows 11. Disable the ``App Installer`` options for the ``python.exe`` and ``python3.exe`` executables.
 
An alternative for running under Windows is to install the `Windows Subsystem for Linux (WSL) <https://learn.microsoft.com/en-us/windows/wsl/install>`_.
This in effect allows you to run PypeIt under Linux under Windows.

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
``astropy``'s `configuration system
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
all files should be distributed with your package installation for these
directories.

As stated above, PypeIt will download remote files and store them in your cache
as they're needed to reduce your data.  I.e., you should mostly be able to
ignore the fact that the relevant files are remote, as long as you're running
the reductions while connected to the internet.  However, if you're preparing to
run a set of reductions and you would prefer to pre-load data that you expect to
need, we provide a few specific scripts for interacting the cache, as described
below.

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
PypeIt.  **Every time you upgrade pypeit, we recommend deleting your existing
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

 - To remove cached files for a specific version: ``pypeit_clean_cache -v 1.15.0``

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

Atmospheric Model Grids
-----------------------

Calculation of the sensitivity functions and general fitting of telluric
absorption uses a grid of model atmosphere spectra.  Files with these
pre-computed model grids are large (few GB).  Recent updates allow for the use
of a PCA decomposition of these models to provided smaller (few MB) reference
files and faster performance; however, both methods are still functional.

For the larger, atmospheric-grid files, note that we provide model spectra for
atmospheric conditions specific to an observatory; however, a model grid is not
provided for all observatories with spectrographs supported by PypeIt.  If you
do not find models for your observatory, you can use the Maunakea model as an
approximation. It includes a large grid of different parameters and should be
good enough for most purposes.

.. note::

    Instruments that anticipate needing a telluric grid have its filename
    already included in the ``telgridfile`` `TelluricPar keyword
    <https://pypeit.readthedocs.io/en/latest/pypeit_par.html#telluricpar-keywords>`__.
    The needed model grid will download automatically when required by the code,
    but given the size of these files and your downlink speed, this may take
    some time.

To install an atmospheric grid or PCA file independent of a reduction, use the
``pypeit_install_telluric`` script, calling the filename of the grid required.
For example, if you need the file ``TelFit_MaunaKea_3100_26100_R20000.fits``,
you would execute:

.. code-block:: console

    $ pypeit_install_telluric TelFit_MaunaKea_3100_26100_R20000.fits

The downloaded file will exist in the PypeIt cache.  Unlike files hosted on
``github``, these files will persist through upgrades of your installation.  To
force the update of a telluric model grid file to the latest version, simply run
``pypeit_install_telluric`` with the ``--force_update`` option.

If you require a telluric grid that is not presently hosted in the cloud, the
code will instruct you to download the file separately from the `PypeIt
dev-suite Google Drive`_.  Users may select any of the files in the Google Drive
for their telluric correction, download them separately, then install them using
the ``--local`` option to ``pypeit_install_telluric``.

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
    cache.  Every time you upgrade your pypeit version, you should delete the
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

Interactive tools in PypeIt are built using the `QT
<https://www.qt.io/>`_ windowing toolkit. The ``qtpy`` package is used to
provide an abstract interface to the two most widely used QT bindings for
Python (see :ref:`dependencies`):

* `pyqt <https://riverbankcomputing.com/software/pyqt/intro>`_
* `PySide <https://wiki.qt.io/Qt_for_Python>`_

At least one of those bindings must be installed for the interactive GUIs to
work. By default ``pypeit`` will install ``pyqt6``. Other backends can be used
by installing them manually via ``pip`` or ``conda`` and then setting the ``QT_API``
environment variable. See the `QtPy documentation <https://github.com/spyder-ide/qtpy>`_
for more details.

C code
------

Significant speed gains in PypeIt can be enabled via compilation of the C
versions of the b-spline fitting code. Compilation of the C code should happen
automatically when you install PypeIt.  However, you can check that the C
code was compiled successfully by running the ``pypeit_c_enabled`` script. What
you should see is:

.. code-block:: console

    $ pypeit_c_enabled
    Successfully imported bspline C utilities.

If no message is printed, the C code could not be imported.

Some notes if you have problems installing the C code:

    - the code will still run successfully by falling back to slower,
      pure-python implementations

    - to successfully compile the C code, you may need to update ``gcc`` and/or
      ``Xcode`` for Mac users

    - for some Mac users, you may also need to update your OS if you're using a
      particularly old version (*e.g.*, 10.10 Yosemite)

Some of the C code uses `OpenMP <https://www.openmp.org/>`_ to parallelize loops
and take advantage of multiple cores/threads. This support is transparent and the code
will work single-threaded if OpenMP is not available. GCC supports OpenMP
out of the box, however the ``clang`` compiler that Apple's XCode provides does not. So
for optimal performance on Apple hardware, you will want to install GCC via ``homebrew``
or ``macports`` and specify its use when installing ``pypeit``. For example, if you installed
GCC via ``homebrew``, you would get ``pypeit`` to use it by doing, for example:

.. code-block:: console

    $ export CC=/opt/homebrew/bin/gcc
    $ pip install pypeit

Basically, ``pypeit`` checks the ``CC`` environment variable for what compiler to use so configure
that as needed to use your desired compiler. The ``pypeit_c_enabled`` script can be used to check if
your compiler has OpenMP support.

ginga Plugins
-------------

PypeIt requires the ``ginga`` viewer and uses at least one ``ginga`` plugin
to enable specific display functionality. No special considerations are needed
to have these plugins installed; however, you can check that they're enabled by
running the following script with the following result::

    $ pypeit_chk_plugins
    [INFO]    :: All required plugins found: SlitWavelength

If the check is unsuccessful, you will see an error message listing
the missing plugins. If you have a problem, please `submit an issue
<https://github.com/pypeit/PypeIt/issues>`__.

----

.. _dependencies:

Package Dependencies
====================

All PypeIt dependencies are installed along with the installation of
PypeIt itself.  Beware this means that packages in your current environment
may be updated depending on the PypeIt version requirements (which is why we
recommend you :ref:`environment` for PypeIt).  The current version
requirements for both users and developers are:

.. include:: include/dependencies_table.rst

Dependency Caveats
------------------

Some users have run into the following complications when installing the
PypeIt dependencies.  If you run into any more, please `submit an issue
<https://github.com/pypeit/PypeIt/issues>`__.

.. TODO: IS THIS FIRST ITEM STILL TRUE?

- At the moment, an implicit dependency on QT bindings remains because of our dependence on ``linetools``.

----

.. _developer_install:

Developer Installation
======================

We, of course, welcome and encourage community development of PypeIt.
Please see our :ref:`codeconduct` and the :ref:`development`.

Developer install via ``pip``
-----------------------------

Install pre-release or development versions of PypeIt directly from `GitHub
<https://github.com/pypeit/PypeIt>`_ using ``pip`` as follows. If you already
have a ``pypeit`` environment set up, run:

.. code-block:: console

    pip install --upgrade "git+https://github.com/pypeit/PypeIt#egg=pypeit"

If you're installing in a clean environment, be sure to include the optional
dependencies as well:

.. code-block:: console

    pip install --upgrade "git+https://github.com/pypeit/PypeIt#egg=pypeit"

These commands will install the default branch, ``release``. You can also
specify a different branch, such as the main ``develop`` branch:

.. code-block:: console

    pip install --upgrade "git+https://github.com/pypeit/PypeIt.git@develop#egg=pypeit"

Commit hashes, tag names, or git refs can also be specified; see the `VCS
Support documentation
<https://pip.pypa.io/en/stable/reference/pip_install/#vcs-support>`_ for details
and examples.

Developer install from source
-----------------------------

Developers doing code development will likely want to set up an "editable"
install that points to a locally checked out copy of the GitHub repository.  We
highly recommended using ``pip`` to install the repository and to
:ref:`environment` for code development.

To install from source (after setting up the python environment), first clone
(your fork of) the repository:

.. code-block:: console

    git clone https://github.com/pypeit/PypeIt.git

Then install the code, include the development dependencies:

.. code-block:: console

    cd PypeIt
    pip install -e ".[dev]"

An "editable" install means that any changes you make in the repository
directory tree will become immediately available the next time the code is
imported. Including the ``[dev]`` set of optional dependencies ensures that all
of the tools you need to test and build PypeIt are installed. (Again, note that
you may or may not need the quotes above depending on your shell, and that you
should avoid cutting and pasting these commands into a terminal window.)

Finally, you may want to add lines to your relevant shell configuration file
(*e.g.*, ``.zshrc`` or ``.bashrc``) that activate the relevant environment
whenever you start up a new shell.  For example:

.. code-block:: console

    conda activate pypeit

Otherwise you will need to always type this command at the terminal prompt to
activate the pypeit environment.

----

.. _test_installation:

Test Your Installation
======================

Tagged versions of PypeIt are extensively tested before distribution.
However, it is worth testing that your installation has been successful, as
follows.

User Tests
----------

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

    python -c "from importlib import resources; print(resources.files('pypeit') / 'tests/files/spec1d_r153-J0025-0312_KASTr_20150123T025323.850.fits')" | xargs -I {} pypeit_show_1dspec {}

If ``pyqt6`` or another Qt backend is correctly installed, this should show a test
spectrum from the Shane/KAST spectrograph.

Developer Tests
---------------

If you performed a developer installation by cloning the repository into a local
directory (*e.g.*, ``~/PypeIt``), you can run the standard unit tests within the
PypeIt environment by executing:

.. code-block:: console

    cd ~/PypeIt
    pytest

.. TODO: WHO AMONG THE CORE DEVELOPERS USE TOX?  WE USE IT FOR CI TESTS, BUT
.. SHOULD WE BE RECOMMENDING IT FOR USERS?

To test within isolated environments and against different versions of various
dependencies, we recommend using ``tox``:

.. code-block:: console

    cd PypeIt
    tox -e test

or, *e.g.*:

.. code-block:: console

    tox -e test-astropydev

Run ``tox -a`` to see a list of available test environments.

In either case, over 100 tests should pass, nearly 100 will be skipped and none
should fail. The skipped tests only run if the PypeIt development is installed
and configured; see :ref:`dev-suite`.


