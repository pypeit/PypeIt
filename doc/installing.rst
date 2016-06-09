.. highlight:: rest

****************
Installing PYPIT
****************

This document will describe how to install PYPIT.
Camille will be editing this.

GSL
===

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
    export GSL_PATH=/usr/local   # Or whatever

Python Dependencies
===================

Here are the list of Python packages required.

* `python <http://www.python.org/>`_ versions 2.7, or 3.3 or later
* `numpy <http://www.numpy.org/>`_ version 1.10 or later
* `astropy`_ version 1.0 or later
* `scipy <http://www.scipy.org/>`_ version 0.17 or later
* `matplotlib <http://matplotlib.org/>`_  version 1.4 or later
* `linetools <https://github.com/linetools/linetools>`_  latest version

Source Code
===========

You are recommended to grab the code from github::

    git clone https://github.com/PYPIT/PYPIT.git

From there, you can build and install either with install or develop, e.g.::

    python setup.py develop

This should generate the Cython files, etc.

Tests
=====

It is worth your while try::

    python setup.py test

There is also a suite of tests that you can download and
run via this Repo:
`TestSuite <https://github.com/PYPIT/PYPIT-development-suite>`_
