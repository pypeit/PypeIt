.. highlight:: rest

****************
Installing PypeIt
****************

This document will describe how to install PypeIt.

Installing Dependencies
=======================
Though we have tried to keep the number of dependencies low,
there are a few packages that need to be installed (various python packages,
GSL, and linetools).

In general, we recommend that you use Anaconda for the majority
of these installations.

Detailed installation instructions are presented below:

Python Dependencies
-------------------

PypeIt depends on the following list of Python packages. 

We recommend that you use `Anaconda <https://www.continuum.io/downloads/>`_ to install and/or update these packages.

* `python <http://www.python.org/>`_ versions 2.7, or 3.5 or later (2.7 will be deprecated)
* `numpy <http://www.numpy.org/>`_ version 1.13 or later
* `astropy <http://www.astropy.org/>`_ version 2.0.5 or later
* `scipy <http://www.scipy.org/>`_ version 0.18.1 or later
* `matplotlib <http://matplotlib.org/>`_  version 1.5.3 or later
* `numba <https://numba.pydata.org/>`_ version 0.39.0 or later (optional - some speed ups, inc. wavecal)
* `PyQT5 <https://wiki.python.org/moin/PyQt/>`_ version 5 (needed for linetools)
* `h5py <https://www.h5py.org/>`_ version 2.6 (for data I/O)
*  yaml -- On Python 3 (at least), you may need to install pyyaml
* `future <https://pypi.python.org/pypi/future/0.6.0>`_ version 0.16 or later
*  configobj -- version 5.0.6 or later

These packages need to be installed by cloning from GitHub:

* `ginga <https://github.com/profxj/ginga>`_ JXP's fork of Ginga

If you are using Anaconda, you can check the presence of these packages with::

	conda list "^python$|numpy|astropy$|scipy$|matplotlib|numba|PyQT|ginga|yaml|h5py"

If the packages have been installed, this command should print out all the packages and their version numbers.

If any of the packages are out of date, they can be updated with a command like::

	conda update scipy

For Ginga, it is currently necessary that you install the fork maintained by JXP::

    git clone https://github.com/profxj/ginga
    cd ginga
    python setup.py install

We hope to make a plug-in for PypeIt instead in the future.

Installing Linetools
--------------------
The latest version of `linetools <https://github.com/linetools/linetools/>`_ is
also required for PypeIt.
Linetools is a package designed for the analysis of 1-D spectra.
The installation steps for linetools are provided
`here <http://linetools.readthedocs.io/en/latest/install.html/>`_.
Do **not** use the pip install.

According to the linetools documentation page, "If you wish to have
full functionality of the GUIs and are using MacOSX, then you probably
need to change your backend from macosx to TkAgg in the matplotlibrc file."


GSL
---

GSL installation
++++++++++++++++

The package complies Cython code that links to gsl routines.
These must be installed on your system prior to PypeIt installation.
We recommend that if you need to install GSL that you use Anaconda,
e.g.::

    conda install -c https://conda.anaconda.org/asmeurer gsl

You are also required to point the ENVIRONMENTAL variable
GSL_PATH to the path above the lib/ and include/ directories
You may determine this path with::

    gsl-config --prefix

It is possible you will also need to set the
LD_LIBRARY_PATH environmental variable to the gsl lib directory,
e.g.::

    export LD_LIBRARY_PATH=/u/xavier/anaconda/lib

.. _GSLELCAPITAN:

GSL on Mac OSX El Capitan
+++++++++++++++++++++++++
.. warning::

	**The above method for installing GSL with Anaconda will not work
	if you are using Mac OSX El Capitan!**

The Mac OSX El Capitan operating system introduced
"Sytem Integrity Protection" (SIP), which restricts root access to as well
as the creation of symlinks in SIP-protected folders (ex: /usr, /bin etc).
The /Users folder, where Anaconda generally installs packages,
is also SIP-protected. This means that the relative paths produced by
some of our Cython code are interfered with by SIP and will cause PypeIt to crash.

Here are some hacks to make the anaconda installation work as
well as some alternate installation instructions:

**1) Replace relative paths in compiled Cython files with full path** 
::

	 #in this example, GSL is installed in '/Users/USERNAME/anaconda/lib/'
	 cd PypeIt/pypeit/
	 install_name_tool -change "@rpath/./libgsl.0.dylib" "/Users/USERNAME/anaconda/lib/libgsl.0.dylib" arcyextract.so
	 install_name_tool -change "@rpath/./libgslcblas.0.dylib" "/Users/USERNAME/anaconda/lib/libgslcblas.0.dylib" arcyextract.so
	 install_name_tool -change "@rpath/./libgsl.0.dylib" "/Users/USERNAME/anaconda/lib/libgsl.0.dylib" arcytrace.so
	 install_name_tool -change "@rpath/./libgslcblas.0.dylib" "/Users/USERNAME/anaconda/lib/libgslcblas.0.dylib" arcytrace.so
	 install_name_tool -change "@rpath/./libgsl.0.dylib" "/Users/USERNAME/anaconda/lib/libgsl.0.dylib" arcycomb.so
	 install_name_tool -change "@rpath/./libgslcblas.0.dylib" "/Users/USERNAME/anaconda/lib/libgslcblas.0.dylib" arcycomb.so
	 install_name_tool -change "@rpath/./libgsl.0.dylib" "/Users/USERNAME/anaconda/lib/libgsl.0.dylib" arcyproc.so
	 install_name_tool -change "@rpath/./libgslcblas.0.dylib" "/Users/USERNAME/anaconda/lib/libgslcblas.0.dylib" arcyproc.so
	 install_name_tool -change "@rpath/./libgsl.0.dylib" "/Users/USERNAME/anaconda/lib/libgsl.0.dylib" arcyutils.so
	 install_name_tool -change "@rpath/./libgslcblas.0.dylib" "/Users/USERNAME/anaconda/lib/libgslcblas.0.dylib" arcyutils.so
	 install_name_tool -change "@rpath/./libgsl.0.dylib" "/Users/USERNAME/anaconda/lib/libgsl.0.dylib" arcyarc.so
	 install_name_tool -change "@rpath/./libgslcblas.0.dylib" "/Users/USERNAME/anaconda/lib/libgslcblas.0.dylib" arcyarc.so
	 

**2) Disable System Integrity Protection**

This is a last resort solution and we do not
recommend it due to security concerns. Instructions for how
to do this can be
found `here <https://www.quora.com/How-do-I-turn-off-the-rootless-in-OS-X-El-Capitan-10-11/>`_.


**3) Install GSL with Homebrew instead of Anaconda**

Since `Homebrew <http://brew.sh/>`_ installs programs in /usr/local , which is not SIP protected, this should work without additional hacks.::

  brew install gsl

in which case the ``GSL_PATH`` variable should be set to ``/usr/local/Cellar/gsl/1.16/``, where ``1.16`` might have to
be replaced with whatever version number you have installed.

Since Homebrew installs programs in /usr/local , which is not
SIP protected, this should work without additional hacks.


Installing PypeIt
================

We recommend that you grab the code from github::

	#go to the directory where you would like to install PypeIt.
	git clone https://github.com/pypeit/pypeit.git

From there, you can build and install either with install or develop, e.g.::

	cd PypeIt
	python setup.py develop

or::

	cd PypeIt
	python setup.py install

This should compile all the necessary Cython files, etc.

If your python installation requires root access, you'll need to use sudo with the "-E" option to pass environment variables.

	sudo -E python setup.py develop


Tests
=====
In order to assess whether PypeIt has been properly installed,
we suggest you run the following tests:

1. Ensure run_pypeit works
-------------------------
Go to a directory outside of the PypeIt directory (e.g. your home directory),
then type run_pypeit.::

	cd
	run_pypeit


2. Run the PypeIt unit tests
---------------------------

Enter the PypeIt directory and do::

    python setup.py test


3. Try the test suite
---------------------
We have provided a suite of tests that you can download and run via this Repo:
`TestSuite <https://github.com/pypeit/PypeIt-development-suite>`_

It can be installed as follows::

	# we suggest installing this in the directory above PypeIt
	git clone https://github.com/pypeit/PypeIt-development-suite.git

To run the test::

	cd PypeIt-development-suite
	./pypeit_test all

.. note::

	pypeit_test can also take the argument kast instead of all. 


The test takes a while to run but should run without issue if all the packages have been properly installed. 


**If you installed GSL with anaconda, a common error from running ./pypeit_test all is:**

|[BUG]     :: There appears to be a bug on Line 7 of arproc.py with error:

| dlopen(/Users/USERNAME/software/PypeIt/pypeit/arcyextract.so, 2): Library not loaded: @rpath/./libgsl.0.dylib

| Referenced from: /Users/USERNAME/software/PypeIt/pypeit/arcyextract.so


**To fix this bug:**

a) Make sure GSL_PATH and LD_LIBRARY_PATH are defined in your .bashrc or .tcshrc file and that the appropriate rc file has been sourced

b) If that does not work, check out :ref:`GSLELCAPITAN`.
