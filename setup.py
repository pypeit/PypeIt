#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# NOTE: The configuration for the package, including the name, version, and
# other information are set in the setup.cfg file.

import os
import sys

from setuptools import setup

from extension_helpers import get_extensions

# First provide helpful messages if contributors try and run legacy commands
# for tests or docs.

TEST_HELP = """
Note: running tests via 'python setup.py test' is now deprecated. The recommended method
is to run:

    tox -e test-alldeps

The Python version can also be specified, e.g.:

    tox -e py38-test-alldeps

You can list all available environments by doing:

    tox -a

If you don't already have tox installed, you can install it by doing:

    pip install tox

If you want to run all or part of the test suite within an existing environment,
you can use pytest directly:

    pip install -e .[dev]
    pytest

For more information, see:

  http://docs.astropy.org/en/latest/development/testguide.html#running-tests
"""

if 'test' in sys.argv:
    print(TEST_HELP)
    sys.exit(1)

VERSION_TEMPLATE = """
# Note that we need to fall back to the hard-coded version if either
# setuptools_scm can't be imported or setuptools_scm can't determine the
# version, so we catch the generic 'Exception'.
try:
    from setuptools_scm import get_version
    version = get_version(root='..', relative_to=__file__)
except Exception:
    version = '{version}'
""".lstrip()

setup(use_scm_version={'write_to': os.path.join('pypeit', 'version.py'),
                       'write_to_template': VERSION_TEMPLATE},
      ext_modules=get_extensions())
