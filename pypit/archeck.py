"""
Version checking.
"""

from __future__ import absolute_import, division, print_function

from distutils.version import LooseVersion

import warnings
from pypit import armsgs

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

class VersionError(Exception):
    pass

minimum_versions = {'scipy': '0.17.0'}

def version_check():
    '''
    Raises an error if there is a mismatched dependency.
    '''
    # loop through dependencies and versions
    for dep, ver in minimum_versions.items():
        if LooseVersion(globals()[dep].__version__) < LooseVersion(ver):
            raise VersionError('Update ' + dep + ' to at least version ' + ver + '!')
