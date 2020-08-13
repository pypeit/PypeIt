============
Installation
============

This document describes how to install PypeIt for both users and developers.

Installing Dependencies
=======================

Installing pypeit will not automatically ensure that all the
dependencies (and their appropriate versions) are installed and
up-to-date.  This must be done manually.

The required packages for pypeit are listed in the `requirements.txt
file
<https://github.com/pypeit/PypeIt/blob/master/pypeit/requirements.txt>`__.
If you download the file, you can install/update all the dependencies
by executing::

    pip install -r requirements.txt

However, we highly recommend that you use `Anaconda
<https://www.continuum.io/downloads/>`_, as described below. The
package and version requirements for pypeit are:

* `python <http://www.python.org/>`_ version 3.7 or later
* `numpy <http://www.numpy.org/>`_ version 1.18.0 or later
* `astropy <http://www.astropy.org/>`_ version 4.0 or later
* `scipy <http://www.scipy.org/>`_ version 1.4 or later
* `matplotlib <http://matplotlib.org/>`_  version 3.1 or later
* `numba <https://numba.pydata.org/>`_ version 0.39.0 or later (optional - some speed ups, inc. wavecal)
* `PyQT5 <https://wiki.python.org/moin/PyQt/>`_ version 5 (needed for linetools)
* `pyyaml <https://pyyaml.org/>`_ -- version 5.1
* `configobj <https://pypi.org/project/configobj/>`_ -- version 5.0.6 or later
* `scikit-learn <https://scikit-learn.org/stable/>`_ -- version 0.20 or later
* `IPython <https://ipython.org>`_ -- version 7.2.0 or later
* `extension_helpers <https://pypi.org/project/extension-helpers/>`_ -- version 0.1 or later
* `ginga <https://pypi.org/project/ginga/>`_ -- version 3.0 or later
* `packaging <https://pypi.org/project/packaging/>`_ -- version 19.0 or later
* `linetools <https://pypi.org/project/linetools/>`_ -- version 0.2 or later (see also below)
* `pytest <https://pypi.org/project/pytest/>`_ -- version 3.0.7 or later

With Anaconda, you can check the presence of these packages with::

	conda list "^python$|numpy|astropy$|scipy$|matplotlib|numba|PyQT|ginga|pyyaml|h5py"

If the packages have been installed, this command should print out
all the packages and their version numbers.

If any of the packages are out-of-date, they can be updated with a
command like::

	conda update scipy

The only exception (true for ginga?) is that you *must* use ``pip
install`` to install ``linetools``::

    pip install linetools

Installing PypeIt
=================

Please read all of the text in this sub-section before choosing how you
wish to install ``PypeIt``.

There are two methods you can use to install ``PypeIt``:

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

A simple check of the installation can be done by trying to call the
main ``PypeIt`` executable to provide the usage help::

    run_pypeit -h

This should fail if any of the requirements are not satisfied. As
stated above, you can download the `requirements.txt file
<https://github.com/pypeit/PypeIt/blob/master/pypeit/requirements.txt>`__
to install/update all of the dependencies::

    pip install -r requirements.txt

Install from the git source
---------------------------

If ``pip`` is unsuccessful or if you are planning to use any of the
``PypeIt`` development branches, then you should install directly
from GitHub.

First, clone the repository::

    git clone https://github.com/pypeit/PypeIt.git

This will create a ``PypeIt/`` directory in your current path. Then,
install with::

    cd PypeIt
    python setup.py develop

(Installing this way ensures that virtually all changes to files in
the ``PypeIt/`` directory take immediate effect the next time you
import the code.)

Again, do a test run as above to check you have all the requirements.

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

Testing the Installation
========================

In order to assess whether ``PypeIt`` has been properly installed, we
suggest you run the following tests:

1. Ensure ``run_pypeit`` works
------------------------------

Go to a directory outside of the ``PypeIt`` directory (e.g. your home
directory) and run the main executable. E.g.,::

	cd
	run_pypeit -h


2. Run the ``PypeIt`` unit tests
--------------------------------

If you cloned the repo (i.e., you did *not* use `pip`) then you can
run the standard tests by doing::

    cd PypeIt
    python setup.py test

Or, alternatively::

    cd PypeIt/pypeit/tests
    python -m pytest . -W ignore


3. Try running ``PypeIt`` on our extensive development test suite
-----------------------------------------------------------------

**This should only be done for those wishing to develop the code.**
Regardless, ask for help if you really want to do this.

We have provided a suite of tests that you can download and run via
this repo: `TestSuite
<https://github.com/pypeit/PypeIt-development-suite>`_

Install it by cloning the GitHub repository (do **not** install this
in the ``PypeIt`` source directory tree)::

	git clone https://github.com/pypeit/PypeIt-development-suite.git

To run the test::

	cd PypeIt-development-suite
	./pypeit_test develop

.. note::

	``pypeit_test`` can also, e.g., take the argument ``kast``
	instead of ``develop`` to only test data from the Shane Kast
	spectrograph.

The test takes a while to run but should run without issue if all the
packages have been properly installed.

Developers
==========

For developers, see :doc:`development`.

