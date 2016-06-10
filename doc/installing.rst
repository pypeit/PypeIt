.. highlight:: rest

****************
Installing PYPIT
****************

This document will describe how to install PYPIT.

Installing Dependencies
=======================
Though we have tried to keep the number of dependencies low, there are a few packages that need to be installed (various python packages, GSL, and linetools).

In general, we recommend that you use Anaconda for the majority of these installations. 

Detailed installation instructions are presented below:

Python Dependencies
-------------------

PYPIT depends on the following list of Python packages. 

We recommend that you use `Anaconda <https://www.continuum.io/downloads/>`_ to install and/or update these packages.

* `python <http://www.python.org/>`_ versions 2.7, or 3.3 or later
* `numpy <http://www.numpy.org/>`_ version 1.10 or later
* `astropy <http://www.astropy.org/>`_ version 1.0 or later
* `scipy <http://www.scipy.org/>`_ version 0.17 or later
* `matplotlib <http://matplotlib.org/>`_  version 1.4 or later
* `PyQT4 <https://wiki.python.org/moin/PyQt/>`_ version 4 (needed for linetools)

If you are using Anaconda, you can check the presence of these packages with::

	conda list "^python$|numpy|astropy$|scipy$|matplotlib|PyQT"

If the packages have been installed, this command should print out all the packages and their version numbers.  

If any of these packages are missing you can install them with a command like::

	conda install PyQT

If any of the packages are out of date, they can be updated with a command like::

	conda update scipy

Installing Linetools
--------------------
The latest version of `Linetools <https://github.com/linetools/linetools/>`_ is also required for PYPIT. 
Linetools is a package designed for the analysis of 1-D spectra. The installation steps for linetools are provided below but much more detailed documentation is provided `here <http://linetools.readthedocs.io/en/latest/install.html/>`_. 

In addition to the python packages mentioned above, linetools requires the installation of `specutils <https://github.com/astropy/specutils/>`_. Specutils can either be cloned directly from github or installed via pip.

To install specutils with pip, use the command::

	pip install --no-deps specutils

Linetools can be installed from github or using pip. 
To install a development version of linetools from the source::

	mkdir linetools
	cd linetools
	git clone https://github.com/linetools/linetools.git
	python setup.py develop

To install linetools using pip::

	pip install linetools

According to the linetools documentation page, "If you wish to have full functionality of the GUIs and are using MacOSX, then you probably need to change your backend from macosx to TkAgg in the matplotlibrc file."

If you'd like to make sure that linetools has been properly installed, you can run the following test::

	python -c 'import linetools; linetools.test()'

If the test reports no errors, then linetools was installed without issue. 

GSL
---

The package complies Cython code that links to gsl routines.
These must be installed on your system prior to PYPIT installation.
We recommend that if you need to install GSL that you use Anaconda,
e.g.::

    conda install -c https://conda.anaconda.org/asmeurer gsl

Beware:  multiple installations of libgsl is likely to cause a problem..

You are also required to point the ENVIRONMENTAL varible
GSL_PATH to the path above the lib/ and include/ files.

You should be able to determine that path with::

    gsl-config --prefix

In my case, anaconda installs packages in /Users/USERNAME/anaconda 

You then need to edit your .bashrc, .tcshrc, or .cshrc file with this path.
If you are using bash (and the path mentioned above as an example)::

	export GSL_PATH='/Users/USERNAME/anaconda'

If instead you are using tcsh or csh::

	setenv GSL_PATH "/Users/USERNAME/anaconda"

Installing PYPIT
================

We recommend that you grab the code from github::

    mkdir PYPIT
	cd PYPIT
	git clone https://github.com/PYPIT/PYPIT.git

From there, you can build and install either with install or develop, e.g.::

    python setup.py develop

or::

	python setup.py install

This should compile all the necessary Cython files, etc.

**If you will be editing/contributing to the code:** 

Note that the default branch you will be in once you clone PYPIT is the master branch. You may want to switch to a different branch by using a command like::

	 git checkout OTHERBRANCHNAME

Check which branch you are in with the command::

	git branch

(an asterisk will appear next to the branch that is currently checked out)

Tests
=====
In order to assess whether PYPIT has been properly installed, we suggest you run the following tests:

1. Ensure run_pypit works
----------------------------
Go to a directory outside of the PYPIT directory (e.g. your home directory), then type run_pypit. 
::
	cd
	run_pypit

Note: This may not work if you are checked into a branch that is not the master branch (e.g. new_docs). 

If this step is failing, try::

	git checkout master

Then try again. 


2. Try the test suite
---------------------
We 
::
    python setup.py test

There is also a suite of tests that you can download and
run via this Repo:
`TestSuite <https://github.com/PYPIT/PYPIT-development-suite>`_
