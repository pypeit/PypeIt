============
Installation
============

This document describes how to install ``PypeIt`` for both users and developers.

----

.. _dependencies:

Installing Dependencies
=======================

Installing ``PypeIt`` will not automatically ensure that all the
dependencies (and their appropriate versions) are installed and
up-to-date. Below we provide two ways of ensuring that the relevant
dependencies are available.

The package and version requirements for ``PypeIt`` are:

* `python <http://www.python.org/>`_ version 3.7 or later
* `numpy <http://www.numpy.org/>`_ version 1.18.0 or later
* `astropy <http://www.astropy.org/>`_ version 4.0 or later
* `scipy <http://www.scipy.org/>`_ version 1.4 or later
* `matplotlib <http://matplotlib.org/>`_  version 3.1 or later
* `numba <https://numba.pydata.org/>`_ version 0.39.0 or later
* `PySide2 <https://wiki.qt.io/Qt_for_Python>`_ version 5
* `pyyaml <https://pyyaml.org/>`_ -- version 5.1
* `configobj <https://pypi.org/project/configobj/>`_ -- version 5.0.6 or later
* `scikit-learn <https://scikit-learn.org/stable/>`_ -- version 0.20 or later
* `IPython <https://ipython.org>`_ -- version 7.2.0 or later
* `ginga <https://pypi.org/project/ginga/>`_ -- version 3.0 or later
* `requests <https://requests.readthedocs.io/en/master/>`_ -- version 2.23 or later
* `packaging <https://pypi.org/project/packaging/>`_ -- version 19.0 or later
* `linetools <https://pypi.org/project/linetools/>`_ -- version 0.2 or later
* `extension_helpers <https://pypi.org/project/extension-helpers/>`_ -- version 0.1 or later
* `shapely <https://pypi.org/project/Shapely/>`_ -- version 1.7 or later; optional, **required for KCWI only**
* `pytest <https://pypi.org/project/pytest/>`_ -- version 3.0.7 or later; optional, developers only

Developer-only items
--------------------

If you are developing, you may need the following packages:

* `sphinx <https://www.sphinx-doc.org/en/master/>`_ -- version 4.0 or later
* sphinx_automodapi (pip install only)
* sphinx_rtd_theme (pip install only)

Create a conda environment (recommended)
----------------------------------------

We highly recommend using `Anaconda <https://www.anaconda.com/>`_ as
a package and environment manager. We provide a yaml file that can be
used to setup a conda environment called ``PypeIt``.  To use this:

 #. Download `environment.yml <https://github.com/pypeit/PypeIt/blob/master/environment.yml>`__.

 #. Create the conda environment::

        conda env create -f environment.yml

 #. Activate it::

        conda activate pypeit

 #. Verify that the new environment was installed correctly::

        conda env list

Install via ``pip`` 
-------------------

To install the dependencies using `pip <https://pypi.org/project/pip/>`_:

 #. Download `requirements.txt <https://github.com/pypeit/PypeIt/blob/master/pypeit/requirements.txt>`__.

 #. Install the dependencies::

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
 
 - At the moment, an implicit dependency on PyQt5 remains (in addition to
   PySide2) because of our dependence on ``linetools``.

 - Note that ``shapely`` is listed as an optional dependency, but is only
   currently used by one method that calculates the spaxel area for KCWI
   output datacubes.

 - For the developer-only (``Sphinx``) packages, download
   `requirements_doc.txt <https://github.com/pypeit/PypeIt/blob/master/requirements_doc.txt>`_
   and install with ``pip install -r requirements_doc.txt``.

----

Installing PypeIt
=================

Please read all of the text in this sub-section before choosing which
of the two methods described below for how you wish to install
``PypeIt``. Once you've completed the installation steps, you should
:ref:`test_installation`.

Install using pip
-----------------

If you are not using code on the edge of development, then
we recommend that you install ``PypeIt`` with ``pip``::

    pip install pypeit

This has been known to fail on some systems (and we're working to fix
the issue). If you have problems, instead try::

    pip install git+https://github.com/pypeit/PypeIt.git

If that also fails, please `submit an issue
<https://github.com/pypeit/PypeIt/issues>`__.

Install from the git source
---------------------------

If ``pip`` is unsuccessful or if you are planning to use any of the
``PypeIt`` development branches, then you should install directly
from GitHub.

 #. Clone the repository::

        git clone https://github.com/pypeit/PypeIt.git

 #. This will create a ``PypeIt`` directory in your current path. To install::

        cd PypeIt
        python setup.py develop

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
install`` or ``python setup.py develop``. You can check that the C
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

If you cloned the repo (i.e., you did *not* use `pip`) then you can
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

