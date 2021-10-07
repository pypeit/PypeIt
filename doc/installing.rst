
.. include:: include/links.rst

.. _installing:

============
Installation
============

.. DO WE HAVE A RELEVANT LINK FOR THE PYPEIT USERS SLACK?

Below, we provide detailed instructions for installing ``PypeIt``.  For
troubleshooting, please consult our ``PypeIt`` user community via our PypeIt
Users Slack and/or `submit an issue <https://github.com/pypeit/PypeIt/issues>`__
on GitHub.

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

Both methods discussed below for installing ``PypeIt`` (via `pip`_ or `conda`_)
also install or upgrade its :ref:`dependencies`.  For this reason, we highly
(!!) recommended you first set up a clean python environment in which to install
``PypeIt``.  This mitigates any possible dependency conflicts with other
packages you use.

You can setup a new python environment using `virtualenv`_:

.. code-block:: console

    virtualenv pypeit
    source pypeit/bin/activate

or `conda`_:

.. code-block:: console

    conda create -n pypeit python=3.8
    conda activate pypeit

See the `Virtualenv documentation <https://virtualenv.pypa.io/en/latest/>`_
and/or `Managing Environments with Conda
<https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_
for more details. See also `virtualenvwrapper
<https://virtualenvwrapper.readthedocs.io/en/latest/>`_ as an option for more
easily managing `virtualenv`_ environments.

Install via ``pip`` (recommended)
---------------------------------

To install the latest release of ``PypeIt`` and its required dependencies, execute
either

.. code-block:: console

    pip install "pypeit[pyqt5]"

to select the ``PyQT5`` QT bindings or

.. code-block:: console

    pip install "pypeit[pyside2]"

to select ``PySide2``; see :ref:`interactive`.

If you are generating datacubes (and performing an astrometric correction), you
will also need the `scikit-image`_ package. It can be installed by including it
in the optional dependencies, e.g.:

.. code-block:: console

    pip install "pypeit[pyside2,scikit-image]"

.. note::

    Whether or not it is correct syntax to use the quotes in the commands above
    depends on your shell.  The above commands are specific to ZShell, whereas
    you don't need the quotes in Bash.  But, in any case, you should avoid
    copying these commands from your browser since the unicode for quotation
    marks may not be correct, leading to errors when they are directly pasted
    into a terminal window.

Install via ``conda``
---------------------

`conda`_ is a popular and widely-used package and environment manager. We
provide a yaml file that can be used to setup a conda environment called
``pypeit`` that contains all of the required dependencies.  To use this:

    #. Download `environment.yml
       <https://github.com/pypeit/PypeIt/blob/release/environment.yml>`__.

    #. Create the conda environment:

        .. code-block:: console

            conda env create -f environment.yml

    #. Activate it:

        .. code-block:: console

            conda activate pypeit

    #. Verify that the new environment was installed correctly:

        .. code-block:: console

            conda env list

    #. Install latest ``pypeit`` via ``pip`` as above or perform a developer
       install as below or install the latest release from ``conda-forge``:

        .. code-block:: console

            conda install -c conda-forge pypeit

    #. Install the preferred QT binding either via ``pip`` as above or via
       conda:

        .. code-block:: console

            conda install -c conda-forge pyqt

       or

        .. code-block:: console

            conda install -c conda-forge pyside2

Upgrading to a new version
--------------------------

Upgrading ``PypeIt`` should simply be a matter of executing:

.. code-block:: console

    pip install pypeit --upgrade

or 

.. code-block:: console

    conda install pypeit --upgrade

depending on how you first installed ``PypeIt``.  If this causes problems (e.g.,
a new ``PypeIt`` script is unavailable or you encounter script errors), first
try uninstalling (e.g., ``pip uninstall pypeit``) and then reinstalling.

----

.. _data_installation:

Additional Data
===============

Some data used by ``PypeIt`` is not kept in the GitHub repository or distributed
via `pip`_ because of its large size.  These include:

 - Raw data included in our development suite used for extensive testing of the code base,
 - Atmospheric model grids used for telluric correction and flux calibration, and
 - Canned data-reduction products used by quick-look scripts.

These are all located in our open-access `PypeIt dev-suite Google Drive`_.

.. note::

    We continue to work on cleaner installation solutions for these data
    products, particularly for the telluric grids and the quick-look master
    files.  In the meantime, note that you will likely need to re-run the
    data-specific installation scripts described below every time you upgrade
    your installation (via `pip`_ or `conda`_).

Raw Data
--------

Example raw data for all supported spectrographs are used in extensive testing
of the code base during development; see :ref:`dev-suite`.  General users should
not need access to these data; however, they may be useful for learning how to
use ``PypeIt`` before running it on your own data from the same instrument.
These data are stored in the ``RAW_DATA`` directory in the `PypeIt dev-suite
Google Drive`_, divided into subdirectories for each instrument and instrument
setup.  See also the `PypeIt-development-suite`_ GitHub repository, which
includes a :doc:`pypeit_file` for each instrument and setup used during
development testing.

Atmospheric Model Grids
-----------------------

Calculation of the sensitivity functions for IR instruments and general fitting
of telluric absorption uses a grid of model atmosphere spectra.  These model
grids are provided in the ``Telluric`` directory in the `PypeIt dev-suite Google
Drive`_ and range in size from 3.5-7.7 GB.  Each file provides model spectra for
atmospheric conditions specific to an observatory; however, a model grid is not
provided for all observatories with spectrographs supported by ``PypeIt``.  If
you do not find models for your observatory, you can use the Maunakea model as
an approximation. It includes a large grid of different parameters and should be
good enough for most purposes.

To install the model grids:

    #. Download the relevant file(s) from the Telluric directory in the `PypeIt
       dev-suite Google Drive`_.

    #. Run the ``pypeit_install_telluric`` script.  For example, if you've
       downloaded ``TelFit_MaunaKea_3100_26100_R200000.fits`` to ``my_path``,
       then you would execute:

        .. code-block:: console

            pypeit_install_telluric --path my_path

       or

        .. code-block:: console

            cd my_path
            pypeit_install_telluric

The ``pypeit_install_telluric`` script simply searches for all the ``TelFit*``
files in the path it is provided and creates symlinks to those files in the
``pypeit/data/telluric/atm_grids`` directory in your installation of ``PypeIt``.

.. warning::

    The installation script simply creates symlinks to the downloaded data.
    This means that if you move the original data, the symlinks will become
    broken **and you will need to rerun the installation script.**


Quick-look Master Files
-----------------------

Some of the quick-look reductions provided by ``PypeIt`` require canned master
files to speed up the data-reduction process, as appropriate for a quick-look
result.  These files are hosted in the ``QL_MASTERS`` directory in the `PypeIt
dev-suite Google Drive`_.

To install the quick-look master files:

    #. Right-click on the ``QL_MASTERS`` folder in the `PypeIt dev-suite
       Google Drive`_ and select the "Download" option from the drop-down menu.
       This will download a zip file containing the full directory contents.
       Its current size (as of 22 July 2021) is about 35 MB.

    #. Run the ``pypeit_install_ql_masters`` script.  E.g.:

        .. code-block:: console

            pypeit_install_ql_masters --zip ~/Downloads/QL_MASTERS-20210722T162355Z-001.zip --odir my_path

The ``pypeit_install_ql_masters`` script will unzip the downloaded file in the
``my_path`` directory and create a symlink to the extracted directory in the
``pypeit/data/`` directory of your ``PypeIt`` installation.  The script can
automatically delete the zip file using the ``--rmzip`` option.  If you already
have the ``QL_MASTERS`` directory, you can also use the script to simply create
the symlink using the ``--ql_path`` option.

.. warning::

    The installation script simply creates symlinks to the downloaded data.
    This means that if you move the original data, the symlinks will become
    broken **and you will need to rerun the installation script.**

----

.. _notes:

Important Package Notes
=======================

.. _interactive:

Interactive Tools
-----------------

Interactive tools in ``PypeIt`` are built using the `QT
<https://www.qt.io/>`_ windowing toolkit. The ``qtpy`` package is used to
provide an abstract interface to the two most widely used QT bindings for
Python (see :ref:`dependencies`):

* `pyqt5 <https://riverbankcomputing.com/software/pyqt/intro>`_ 
* `PySide2 <https://wiki.qt.io/Qt_for_Python>`_ 

At least one of those bindings must be installed for the interative GUIs to
work. **Do not install both!**  These two packages do not play nicely together.
We strongly recommend that you use ``pyqt5``, unless you are attracted to the
more flexible licensing that ``PySide2`` provides.  ``PySide2`` can occasionally
cause GUIs to crash because of conflicts with other packages in your environment
that use ``pyqt5`` (all the more reason to isolate your ``PypeIt`` installation
in its own environment).

C code
------

Significant speed gains in ``PypeIt`` can be enabled via compilation of the C
versions of the b-spline fitting code. Compilation of the C code should happen
automatically when you install ``PypeIt``.  However, you can check that the C
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
      particularly old version (e.g., 10.10 Yosemite)

ginga Plugins
-------------

``PypeIt`` requires the ``ginga`` viewer and uses at least one ``ginga`` plugin
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

All ``PypeIt`` dependencies are installed along with the installation of
``PypeIt`` itself.  Beware this means that packages in your current environment
may be updated depending on the ``PypeIt`` version requirements (which is why we
recommend you :ref:`environment` for ``PypeIt``).  The current version
requirements for both users and developers are:

.. include:: include/dependencies_table.rst

Dependency Caveats
------------------

Some users have run into the following complications when installing the
``PypeIt`` dependencies.  If you run into any more, please `submit an issue
<https://github.com/pypeit/PypeIt/issues>`__.

.. IS THIS FIRST ITEM STILL TRUE?

- At the moment, an implicit dependency on QT bindings remains (either PyQT5 or
  PySide2) because of our dependence on ``linetools``.

- Note that ``shapely`` is provided as an optional dependency, but is only
  currently used by one method that calculates the spaxel area for KCWI output
  datacubes.

----

.. _developer_install:

Developer Installation
======================

We, of course, welcome and encourage community development of ``PypeIt``.
Please see our :ref:`codeconduct` and the :ref:`development`.

Install via ``pip``
-------------------

Install pre-release or development versions of ``PypeIt`` directly from `GitHub
<https://github.com/pypeit/PypeIt>`_ using ``pip`` as follows. If you already
have a ``pypeit`` environment setup, run:

.. code-block:: console

    pip install --upgrade "git+https://github.com/pypeit/PypeIt#egg=pypeit"

If you're installing in a clean environment, be sure to include the optional
dependencies as well:

.. code-block:: console

    pip install --upgrade "git+https://github.com/pypeit/PypeIt#egg=pypeit[pyqt5,shapely]"

These commands will install the default branch, ``release``. You can also
specify a different branch, such as the main ``develop`` branch:

.. code-block:: console

    pip install --upgrade "git+https://github.com/pypeit/PypeIt.git@develop#egg=pypeit[pyqt5,shapely]"

Commit hashes, tag names, or git refs can also be specified; see the `VCS
Support documentation
<https://pip.pypa.io/en/stable/reference/pip_install/#vcs-support>`_ for details
and examples.

Install from source
-------------------

Developers doing code development will likely want to set up an "editable"
install that points to a locally checked out copy of the GitHub repository.  We 
highly recommended using ``pip`` to install the repository and to
:ref:`environment` for code development.

To install from source (after setting up the python enviroment), first clone the
repository:

.. code-block:: console

    git clone https://github.com/pypeit/PypeIt.git

Then install the code, include the development dependencies:

.. code-block:: console

    cd PypeIt
    pip install -e ".[dev,pyqt5]"

An "editable" install means that any changes you make in the repository
directory tree will become immediately available the next time the code is
imported. Including the ``[dev]`` set of optional dependencies ensures that all
of the tools you need to test and build ``PypeIt`` are installed. The ``pyqt5``
installation option instructs the script to use the PyQt5 Qt backend. (Again,
note that you may or may not need the quotes above depending on your shell, and
that you should avoid cutting and pasting these commands into a terminal
window.)

Finally, you may want to add lines to your relevant shell configuration file
(e.g., ``.zshrc`` or ``.bashrc``) that activate the relevant environment
whenever you start up a new shell.  For example:

.. code-block:: console

    conda activate pypeit

Otherwise you will need to always type this command at the terminal prompt to
activate the pypeit environment.

----

.. _test_installation:

Test Your Installation
======================

Tagged versions of ``PypeIt`` are extensively tested before distribution.
However, it is worth testing that your installation has been successful, as
follows.

User Tests
----------

The most basic tests that ``PypeIt`` has been properly installed is to get the
help dialog for one of its main executables.  I.e., from a terminal widow, type:

.. code-block:: console

    run_pypeit -h

A second basic test is to try to import ``PypeIt`` from within a python session.
For example:

.. code-block:: console

    python
    >>> import pypeit

Developer Tests
---------------

If you performed a developer installation by cloning the repository into a local
directory (e.g., ``~/PypeIt``), you can run the standard unit tests within the
``PypeIt`` environment by executing:

.. code-block:: console

    cd ~/PypeIt
    pytest

.. WHO AMONG THE CORE DEVELOPERS USE TOX?  WE USE IT FOR CI TESTS, BUT SHOULD WE BE RECOMMENDING IT FOR USERS?

To test within isolated environments and against different versions of various
dependencies, we recommend using ``tox``:

.. code-block:: console

    cd PypeIt
    tox -e test

or, e.g.:

.. code-block:: console

    tox -e test-astropydev

Run ``tox -a`` to see a list of available test environemts.

In either case, over 100 tests should pass, nearly 100 will be skipped and none
should fail. The skipped tests only run if the PypeIt development is installed
and configured; see :ref:`dev-suite`.


