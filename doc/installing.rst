============
Installation
============

This document describes how to install ``PypeIt`` for both users and developers.

Supported Platforms
===================
PypeIt is only supported for CPython versions 3.7+. PypeIt may work on other versions of Python,
however no guarantees or support will be provided.

Installing with ``conda`` or ``pip``
====================================
The easiest way to install PypeIt is to use ``conda``. Conda is part of the
`Anaconda Python Distribution <https://www.anaconda.com/products/individual>`_,
though we recommend using the Miniconda installer for faster installation and
smaller downloads.

If you use Anaconda on your system, ``conda`` is the recommended way to install
PypeIt. Otherwise, use ``pip``.

Once you've completed the installation steps, you should
:ref:`test_installation`.

Installing with ``conda``
-------------------------
Make sure you have either `Anaconda <https://www.anaconda.com/products/individual>`_,
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_, or
`Miniforge <https://github.com/conda-forge/miniforge#download>`_ installed.
Use this command to install PypeIt::

    conda install -c conda-forge pypeit

----

Installing with ``pip``
-----------------------

  #. First, you need to ensure that PypeIt's dependencies are installed. Download `requirements.txt <https://github.com/pypeit/PypeIt/blob/master/pypeit/requirements.txt>`__.

  #. Install the dependencies with this command::

        pip install -r requirements.txt

  #. Then use this command to install PypeIt::

        pip install pypeit

This has been known to fail on some systems (and we're working to fix
the issue). If you have problems, instead try::

    pip install git+https://github.com/pypeit/PypeIt.git

If that also fails, please `submit an issue
<https://github.com/pypeit/PypeIt/issues>`__.

.. _dependencies:

Installing Dependencies
=======================

Installing ``PypeIt`` (except from ``conda``) will not automatically ensure that all the
dependencies (and their appropriate versions) are installed and
up-to-date. Below we provide two ways of ensuring that the relevant
dependencies are available.

The package and version requirements for ``PypeIt`` are:

* `python <http://www.python.org/>`_ version 3.7 or later
* `numpy <http://www.numpy.org/>`_ version 1.18.0 or later
* `scipy <http://www.scipy.org/>`_ version 1.4 or later
* `astropy <http://www.astropy.org/>`_ version 4.0 or later
* `matplotlib <http://matplotlib.org/>`_  version 3.1 or later
* `scikit-learn <https://scikit-learn.org/stable/>`_ -- version 0.20 or later
* `numba <https://numba.pydata.org/>`_ version 0.39.0 or later
* `QtPy <https://pypi.org/project/qtpy>`_ version 1.9 or later
  * one of PySide2 or PyQt5
* `pyyaml <https://pyyaml.org/>`_ -- version 5.1
* `configobj <https://pypi.org/project/configobj/>`_ -- version 5.0.6 or later
* `IPython <https://ipython.org>`_ -- version 7.2.0 or later
* `requests <https://requests.readthedocs.io/en/master/>`_ -- version 2.23 or later
* `packaging <https://pypi.org/project/packaging/>`_ -- version 19.0 or later
* `extension-helpers <https://pypi.org/project/extension-helpers/>`_ -- version 0.1 or later
* `ginga <https://pypi.org/project/ginga/>`_ -- version 3.0 or later
* `linetools <https://pypi.org/project/linetools/>`_ -- version 0.3.dev2231 or later
* `pytest <https://pypi.org/project/pytest/>`_ -- version 3.0.7 or later
* `shapely <https://pypi.org/project/Shapely/>`_ -- version 1.7 or later; optional, **required for KCWI only**

Developer-only items
--------------------

If you are developing, you may need the following packages for building the documentation:

* `sphinx <https://www.sphinx-doc.org/en/master/>`_ -- version 4.0 or later
* sphinx-automodapi
* sphinx_rtd_theme

How to install dependencies
---------------------------

#. Download `requirements.txt <https://github.com/pypeit/PypeIt/blob/master/pypeit/requirements.txt>`__.

#. Install the dependencies using ``conda``::

        conda install -y -c conda-forge --file requirements.txt

  or ``pip``::

        pip install -r requirements.txt

Note that this is a "system-wide" installation, and will
replace/upgrade any current versions of the packages you already have
installed.

Dependency Caveats
------------------

Some users have run into the following complications when installing the
``PypeIt`` dependencies:

 - Installation of ``numba`` has been known to fault because of an issue with
   the LLVM compiler. For one particular case, the solution was to revert to
   llvm version 9.0.1 using `Homebrew <https://brew.sh/>`_ and then add
   environmental variables to your shell rc that point to the reverted
   directory structure.

 - Note that ``shapely`` is listed as an optional dependency, but is only
   currently used by one method that calculates the spaxel area for KCWI
   output datacubes.

 - Note that you will likely have either PyQt5 or PySide2 installed and PypeIt requires
   one of these two packages installed for any of the GUIs to work, in particular
   ``pypeit_show_2dspec`` and ``pypeit_show_1dspec``. If you have neither installed, and
   both ``pypeit_show_2dspec`` and ``pypeit_show_1dspec`` crash, then you should install
   one of PySide2 or PyQt5 (``pip install pyside2`` or ``pip install pyqt5`` or
   ``conda install pyside2`` or ``conda install pyqt``).

 - For the developer-only (``Sphinx``) packages, download
   `requirements_doc.txt <https://github.com/pypeit/PypeIt/blob/master/requirements_doc.txt>`_
   and install with ``conda install -y -c conda-forge --file requirements_doc.txt`` or
   ``pip install -r requirements_doc.txt``.

----

Install from the git source (Advanced)
======================================

If ``conda`` and ``pip`` are unsuccessful, or if you are planning to use any of the
``PypeIt`` development branches, then you should install directly
from GitHub.

 #. Clone the repository and navigate to the PypeIt directory::

        git clone https://github.com/pypeit/PypeIt.git
        cd PypeIt

 #. Install dependencies

  - into your current environment with ONE of the following commands (you may need to also install
    PySide2/PyQt5 if they are not already installed as per note above)::

        pip install -r pypeit/requirements.txt

        conda install -y -c conda-forge --file requirements.txt

  - or into a new environment::

        conda env create -f environment.yml
        conda activate pypeit
        # install PySide2/PyQt5 from conda/pip using using ONE of the following lines
        conda install -c conda-forge pyside2
        conda install pyqt
        pip install pyside2
        pip install pyqt5

 #. Install PypeIt::

        pip install -e .

Installing the code this way ensures that virtually all changes to files in
the ``PypeIt`` directory take immediate effect the next time you
import the code.

----

Compiled Code and Plug-ins
==========================

C code
------

Significant speed gains in ``PypeIt`` can be enabled via compilation
of the C code version of the b-spline fitting code. Compilation of
the C code should happen automatically when you execute ``pip
install`` or ``pip install -e .``. When installing from ``conda``,
the C compilation has already been done for you. You can check that the C
code was compiled successfully by running the ``pypeit_c_enabled``
script. What you should see is::

    $ pypeit_c_enabled
    Successfully imported bspline C utilities.

If no message is printed, the C code could not be imported.

Some notes if you have problems installing the C code:

    - the code will still run successfully by falling back to slower,
      pure-python implementations
    - to successfully compile the C code, you may need to update
      `gcc` and/or `Xcode` for Mac users
    - for some Mac users, you may also need to update your OS if
      you're using a particularly old version (e.g., 10.10 Yosemite)

ginga Plugins
-------------

``PypeIt`` now (as of version 1.0.7dev) requires the ``ginga`` viewer
and uses at least one ``ginga`` plugin to enable specific display
functionality. No special considerations are needed to have these
plugins installed; however, you can check that they're enabled by
running the following script with the following result::

    $ pypeit_chk_plugins
    [INFO]    :: All required plugins found: SlitWavelength

If the check is unsuccessful, you will see an error message listing
the missing plugins. If you have a problem, please `submit an issue
<https://github.com/pypeit/PypeIt/issues>`__.

----

.. _test_installation:

Test Your Installation
======================

In order to assess whether ``PypeIt`` has been properly installed, we
suggest you run the following tests:

1. Ensure ``run_pypeit`` works
------------------------------

Go to a directory outside of the ``PypeIt`` directory (e.g. your home
directory) and run the main executable. E.g.,::

	cd
	run_pypeit -h

This should fail if any of the requirements are not satisfied; see
:ref:`dependencies`.


2. Run the ``PypeIt`` unit tests
--------------------------------

If you cloned the repo (i.e., you did *not* use `pip` or `conda`) then you can
run the standard tests by doing::

    cd PypeIt
    python setup.py test

Or, alternatively::

    cd PypeIt/pypeit/tests
    python -m pytest . -W ignore

Over 100 tests should pass, nearly 100 will be skipped (unless
you are a developer) and none should fail.

----

Developers
==========

For developers, see :ref:`development`.

Also, test scripts for development purposes are available at the
`PypeIt Development Suite <https://github.com/pypeit/PypeIt-development-suite>`_.

