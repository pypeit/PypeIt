"""
Version checking.
"""

from __future__ import absolute_import, division, print_function

import pkg_resources
import os

requirements_file = pkg_resources.resource_filename('pypit', 'requirements.txt')
install_requires = [line.strip().replace('==', '>=') for line in open(requirements_file)
                    if not line.strip().startswith('#') and line.strip() != '']
for requirement in install_requires:
    pkg, version = requirement.split('>=')
    try:
        pv = pkg_resources.get_distribution(pkg).version
    except pkg_resources.DistributionNotFound:
        raise ImportError("Package: {:s} not installed!".format(pkg))
    else:
        if pkg_resources.parse_version(pv) < pkg_resources.parse_version(version):
            print("Version of package {:s} = {:s}".format(pkg, pv))
            raise ImportError("You need version >= {:s}".format(version))

'''
from distutils.version import LooseVersion

import warnings

try:
    from xastropy.xutils import xdebug as debugger
except ImportError:
    import pdb as debugger

# check these
import numpy as np
import scipy
try:
    import ginga
except ImportError:
    warnings.warn("Ginga is not installed.  You may wish to do so.")
import h5py
import astropy
import matplotlib
import yaml
import linetools
import future


class VersionError(Exception):
    pass

minimum_versions = {'scipy': '0.17.0'}


def version_check():
    """
    Raises an error if there is a mismatched dependency.
    """
    # loop through dependencies and versions
    for dep, ver in minimum_versions.items():
        if LooseVersion(globals()[dep].__version__) < LooseVersion(ver):
            raise VersionError('Update ' + dep + ' to at least version ' + ver + '!')
'''
