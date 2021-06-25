
.. _installing:

============
Installation
============

This document describes how to install ``PypeIt`` for both users and developers.

----

.. _dependencies:

Package Dependencies
=======================

The package and version requirements for ``PypeIt`` currently are:

* `python <http://www.python.org/>`_ version 3.7 or later
* `numpy <http://www.numpy.org/>`_ version 1.18.0 or later
* `astropy <http://www.astropy.org/>`_ version 4.0 or later
* `scipy <http://www.scipy.org/>`_ version 1.4 or later
* `matplotlib <http://matplotlib.org/>`_  version 3.1 or later
* `numba <https://numba.pydata.org/>`_ version 0.39.0 or later
* `pyyaml <https://pyyaml.org/>`_ -- version 5.1
* `configobj <https://pypi.org/project/configobj/>`_ -- version 5.0.6 or later
* `scikit-learn <https://scikit-learn.org/stable/>`_ -- version 0.20 or later
* `IPython <https://ipython.org>`_ -- version 7.2.0 or later
* `ginga <https://pypi.org/project/ginga/>`_ -- version 3.0 or later
* `requests <https://requests.readthedocs.io/en/master/>`_ -- version 2.23 or later
* `packaging <https://pypi.org/project/packaging/>`_ -- version 19.0 or later
* `linetools <https://pypi.org/project/linetools/>`_ -- version 0.2 or later
* `extension_helpers <https://pypi.org/project/extension-helpers/>`_ -- version 0.1 or later
* `qtpy <https://github.com/spyder-ide/qtpy>`_ -- version 1.9 or later

The interactive tools in ``PypeIt`` are built using the `QT <https://www.qt.io/>`_ windowing toolkit. The ``qtpy`` package
is used to provide an abstract interface to the two most widely used QT bindings for Python:

* `pyqt5 <https://riverbankcomputing.com/software/pyqt/intro>`_ -- version 5
* `PySide2 <https://wiki.qt.io/Qt_for_Python>`_ -- version 5

At least one of those bindings must be installed for the interative GUIs to work. DO NOT INSTALL BOTH, as these
two packages do not play nicely together. We strongly recommend that you go with pyqt5, unless you are attracted
to the more flexible licensing that PySide2 provides.  PySide2 can occasionally cause GUIs to crash because
of conflicts with other packages in your environment that use pyqt5.

Developer-only items
--------------------

If you are developing, you will need the following packages:

* `pytest <https://pypi.org/project/pytest/>`_ -- version 3.0.7 or later
* `tox <https://tox.readthedocs.io/en/latest/>`_

Building the documentation requires the following extra dependencies:

* `sphinx <https://www.sphinx-doc.org/en/master/>`_ -- version 4.0 or later
* sphinx_automodapi (pip install only)
* sphinx_rtd_theme (pip install only)

Installing via ``pip`` (recommended)
------------------------------------

The recommended method for installing ``PypeIt`` and its dependencies, both required and optional,
is via `pip <https://pypi.org/project/pip/>`_. It is very highly recommended to first set up a clean environment
in which to install ``PypeIt`` so that possible dependency conflicts can be avoided. This can be done via ``virtualenv``:

        virtualenv pypeit
        source pypeit/bin/activate

or via ``conda``:

        conda create -n pypeit python=3.8
        conda activate pypeit

See the `Virtualenv documentation <https://virtualenv.pypa.io/en/latest/>`_ and/or `Managing Environments with Conda
<https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_ for more details. See also
`virtualenvwrapper <https://virtualenvwrapper.readthedocs.io/en/latest/>`_ as an option for more easily managing
``virtualenv`` virtual environments.

To install the latest release of ``PypeIt`` and its required dependencies, do either:

        pip install "pypeit[pyqt5]"

to select the ``PyQT5`` QT bindings or:

        pip install "pypeit[pyside2]"

to select ``PySide2``.

If you are using KWCI, you may also need the ``shapely`` package. It can be installed by
including it in the optional dependencies, e.g.:

        pip install "pypeit[pyside2,shapely]"

Installing via ``conda``
------------------------

`Conda <https://docs.conda.io/projects/conda/en/latest/index.html>`_ is
a popular and widely-used package and environment manager. We provide a yaml file that can be
used to setup a conda environment called ``pypeit`` that contains all of the required dependencies.  To use this:

 #. Download `environment.yml <https://github.com/pypeit/PypeIt/blob/master/environment.yml>`__.

 #. Create the conda environment::

        conda env create -f environment.yml

 #. Activate it::

        conda activate pypeit

 #. Verify that the new environment was installed correctly::

        conda env list

 #. Install latest ``pypeit`` via ``pip`` as above or perform a developer install as below or install the latest
 release from ``conda-forge``::

        conda install -c conda-forge pypeit

#. Install the preferred QT binding either via ``pip`` as above or via conda::

        conda install -c conda-forge pyqt

or

        conda install -c conda-forge pyside2

Developer Install via ``pip``
-----------------------------

It is also possible to install pre-release or development versions of ``PypeIt`` directly from `GitHub <https://github.com/pypeit/PypeIt>`_
using ``pip``. If you already have a ``pypeit`` environment setup, do:

        pip install --upgrade "git+https://github.com/pypeit/PypeIt#egg=pypeit"

If you're installing in a clean environment, be sure to include the optional dependencies as well:

        pip install --upgrade "git+https://github.com/pypeit/PypeIt#egg=pypeit[pyqt5,shapely]"

Those commands will install the default branch, ``master``. You can also specify the branch you wish to use:

        pip install --upgrade "git+https://github.com/pypeit/PypeIt.git@develop#egg=pypeit[pyqt5,shapely]"

Commit hashes, tag names, or git refs can also be specified. See the `VCS Support documentation
<https://pip.pypa.io/en/stable/reference/pip_install/#vcs-support>`_ for details and examples.


Developer Install from Source
-----------------------------

Developers doing code development will likely want to set up an "editable" install that points to a locally checked out
copy of the GitHub repository. It is highly recommened to use ``pip`` for this as well so that the dependencies can be managed.
It is also recommended to install all optional dependencies within the environment used for ``PypeIt`` development. First,
we highly recommend setting up a clean environment in which to install ``PypeIt`` to avoid conflicts, as described
above. Using ``conda`` you would execute the commands::

        conda create -n pypeit python=3.8
        conda activate pypeit

Second, clone the repository::

        git clone https://github.com/pypeit/PypeIt.git

Then perform the install, by entering the PypeIt directory and running the installation script via pip::

        cd PypeIt
        pip install -e ".[dev,pyqt5]"

An "editable" install means that any changes you make in that code tree will become immediately available the next
time the code is imported. Including the ``[dev]`` set of optional dependencies ensures that all of the tools you
need to test and build ``PypeIt`` are installed. The ``pyqt5`` intructs the script to use the PyQt5 Qt backend. (Note that
you may or may not need the quotes above depending on your shell, and that you should avoid cutting and pasting these
commands since the quotation marks may not paste correctly). Finally, you may want to add::

        conda activate pypeit

to your .bashrc or .tcshrc in order to activate thge pypeit python environment when you launch a shell via the terminal.
Otherwise you will need to always type this command at the terminal prompt to activate the pypeit environment.

If any of this fails, please `submit an issue
<https://github.com/pypeit/PypeIt/issues>`__.


Dependency Caveats
------------------

Some users have run into the following complications when installing the
``PypeIt`` dependencies:

 - Installation of ``numba`` has been known to fault because of an issue with
   the LLVM compiler. For one particular case, the solution was to revert to
   llvm version 9.0.1 using `Homebrew <https://brew.sh/>`_ and then add
   environmental variables to your shell rc that point to the reverted
   directory structure. ``numba`` also does not yet officially support
   Python 3.9.

 - At the moment, an implicit dependency on QT bindings remains (either PyQT5 or
   PySide2) because of our dependence on ``linetools``.

 - Note that ``shapely`` is provided as an optional dependency, but is only
   currently used by one method that calculates the spaxel area for KCWI
   output datacubes.

----

Compiled Code and Plug-ins
==========================

C code
------

Significant speed gains in ``PypeIt`` can be enabled via compilation
of the C code version of the b-spline fitting code. Compilation of
the C code should happen automatically when you execute ``pip
install`` or ``pip install -e .``. You can check that the C
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

If you cloned the repo (i.e., you did *not* use ``pip`` or ``conda``), then you can
run the standard tests within the environment you created for ``PypeIt`` development by doing::

    cd PypeIt
    pytest

To test your development work within isolated environments and against different versions of
various dependencies, it is recommended to use ``tox``::

    cd PypeIt
    tox -e test

or, e.g.::

    tox -e test-astropydev

Run ``tox -a`` to see a list of available test environemts.

In either case, over 100 tests should pass, nearly 100 will be skipped and none should fail. The skipped
tests only run if the PypeIt development is installed and configured.


.. THERE NEEDS TO BE A "DATA" SECTION HERE THAT INSTRUCTS USERS ON HOW TO
.. DOWNLOAD AND INSTALL, E.G., THE TELLURIC GRIDS, QSO MODELS, STANDARD STARS,
.. ETC.!!

----

Developers
==========

For developers, see :ref:`development`.

Also, test scripts for development purposes are available at the
`PypeIt Development Suite <https://github.com/pypeit/PypeIt-development-suite>`_.

