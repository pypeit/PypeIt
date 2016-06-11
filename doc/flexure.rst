.. highlight:: rest

****************
Installing PYPIT
****************

This document will describe how flexure corrections
are performed in PYPIT.

Overview
========

The general approach is to compare the sky model
from the observation with an archived sky model.  By
default, the Paranal sky spectrum is used, as derived
from the SDSS codes.

We have also found that the
Mauna Kea sky is very dark in the blue and also include
specific sky models for LRIS.


-------------------

PYPIT depends on the following list of Python packages. 

We recommend that you use `Anaconda <https://www.continuum.io/downloads/>`_ to install and/or update these packages.

* `python <http://www.python.org/>`_ versions 2.7, or 3.3 or later
* `numpy <http://www.numpy.org/>`_ version 1.10 or later
* `astropy <http://www.astropy.org/>`_ version 1.0 or later
* `scipy <http://www.scipy.org/>`_ version 0.17 or later
* `matplotlib <http://matplotlib.org/>`_  version 1.4 or later
* `PyQT4 <https://wiki.python.org/moin/PyQt/>`_ version 4 (needed for linetools)
* `Ginga <https://ginga.readthedocs.io/en/latest/>`_ latest version (highly recommended)

If you are using Anaconda, you can check the presence of these packages with::

	conda list "^python$|numpy|astropy$|scipy$|matplotlib|PyQT|ginga"

If the packages have been installed, this command should print out all the packages and their version numbers.  

If any of these packages are missing you can install them with a command like::

	conda install PyQT

If any of the packages are out of date, they can be updated with a command like::

	conda update scipy

Installing Linetools
--------------------
The latest version of `Linetools <https://github.com/linetools/linetools/>`_ is also required for PYPIT. 
Linetools is a package designed for the analysis of 1-D spectra. The installation steps for linetools are provided below but much more detailed documentation is provided `here <http://linetools.readthedocs.io/en/latest/install.html/>`_. 

Linetools can be installed from github or using pip. 
To install a development version of linetools from the source::

	#go to the directory where you would like to install linetools
	git clone https://github.com/linetools/linetools.git
	cd linetools
	python setup.py develop

To install linetools using pip::

	pip install linetools

According to the linetools documentation page, "If you wish to have full functionality of the GUIs and are using MacOSX, then you probably need to change your backend from macosx to TkAgg in the matplotlibrc file."

If you'd like to make sure that linetools has been properly installed, you can run the following test::

	python -c 'import linetools; linetools.test()'

The pip installed version will result in quite a few errors. If you cloned the github directory, there will be 1 error from not having specutils installed. This test is old and specutils is no longer required for linetools. 

GSL
---

The package complies Cython code that links to gsl routines.
These must be installed on your system prior to PYPIT installation.
We recommend that if you need to install GSL that you use Anaconda,
e.g.::

    conda install -c https://conda.anaconda.org/asmeurer gsl

Beware:  multiple installations of libgsl is likely to cause a problem..

You are also required to point the ENVIRONMENTAL variable
GSL_PATH to the path above the lib/ and include/ directories and the LD_LIBRARY_PATH variable to the lib/ directory.

You should be able to determine that path with::

    gsl-config --prefix

In my case, anaconda installs packages in /Users/USERNAME/anaconda 

You then need to edit your .bashrc, .tcshrc, or .cshrc file with this path.
If you are using bash (and the path mentioned above as an example)::

	export GSL_PATH='/Users/USERNAME/anaconda'
	export LD_LIBRARY_PATH='/Users/USERNAME/anaconda/lib/'

If instead you are using tcsh or csh::

	setenv GSL_PATH "/Users/USERNAME/anaconda"
	setenv LD_LIBRARY_PATH "/Users/USERNAME/anaconda/lib/"

.. _GSLELCAPITAN:

GSL installation on Mac OSX El Capitan
++++++++++++++++++++++++++++++++++++++
.. warning::

	**The above method for installing GSL with Anaconda will not work if you are using Mac OSX El Capitan!**

The Mac OSX El Capitan operating system introduced "Sytem Integrity Protection" (SIP), which restricts root access to as well as the creation of symlinks in SIP-protected folders (ex: /usr, /bin etc). The /Users folder, where Anaconda generally installs packages, is also SIP-protected. This means that the relative paths produced by some of our Cython code are interfered with by SIP and will cause PYPIT to crash. 


Here are some hacks to make the anaconda installation work as well as some alternate installation instructions:

**1) Replace relative paths in compiled Cython files with full path** 
::

	 #in this example, GSL is installed in '/Users/USERNAME/anaconda/lib/'
	 cd PYPIT/pypit/
	 install_name_tool -change "@rpath/./libgsl.0.dylib" "/Users/cleibler/anaconda/lib/libgsl.0.dylib" arcyextract.so
	 install_name_tool -change "@rpath/./libgslcblas.0.dylib" "/Users/USERNAME/anaconda/lib/libgslcblas.0.dylib" arcyextract.so
	 

**2) Disable System Integrity Protection**

This is a last resort solution and we do not recommend it due to security concerns. Instructions for how to do this can be found `here <https://www.quora.com/How-do-I-turn-off-the-rootless-in-OS-X-El-Capitan-10-11/>`_. 


**3) Install GSL with Homebrew instead of Anaconda**

Since Homebrew installs programs in /usr/local , which is not SIP protected, this should work without additional hacks. 
ASHER CAN YOU WRITE UP SOME INSTRUCTIONS FOR THE HOMEBREW INSTALLATION OF GSL?


Installing PYPIT
================

We recommend that you grab the code from github::

	#go to the directory where you would like to install PYPIT.
	git clone https://github.com/PYPIT/PYPIT.git

From there, you can build and install either with install or develop, e.g.::

	cd PYPIT
	python setup.py develop

or::

	cd PYPIT
	python setup.py install

This should compile all the necessary Cython files, etc.

Tests
=====
In order to assess whether PYPIT has been properly installed, we suggest you run the following tests:

1. Ensure run_pypit works
----------------------------
Go to a directory outside of the PYPIT directory (e.g. your home directory), then type run_pypit. 
::

	cd
	run_pypit


2. Try the test suite
---------------------
We have provided a suite of tests that you can download and run via this Repo:
`TestSuite <https://github.com/PYPIT/PYPIT-development-suite>`_

It can be installed as follows::

	# we suggest installing this in the directory above PYPIT
	git clone https://github.com/PYPIT/PYPIT-development-suite.git

To run the test::

	cd PYPIT-development-suite
	./pypit_test all

.. note::

	pypit_test can also take the argument kast instead of all. 


The test takes a while to run but should run without issue if all the packages have been properly installed. 


**If you installed GSL with anaconda, a common error from running ./pypit_test all is:**

|[BUG]     :: There appears to be a bug on Line 7 of arproc.py with error:

| dlopen(/Users/USERNAME/software/PYPIT/pypit/arcyextract.so, 2): Library not loaded: @rpath/./libgsl.0.dylib

| Referenced from: /Users/USERNAME/software/PYPIT/pypit/arcyextract.so


**To fix this bug:**

a) Make sure GSL_PATH and LD_LIBRARY_PATH are defined in your .bashrc or .tcshrc file and that the appropriate rc file has been sourced

b) If that does not work, check out :ref:`GSLELCAPITAN`.