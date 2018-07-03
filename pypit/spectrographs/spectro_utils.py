""" Utilities for spectograph codes
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)
import numpy as np

from astropy.io import fits

from pypit import msgs
from pypit.spectrographs import deimos
from pypit.spectrographs import lris
from pypit import arparse

from pypit import ardebug as debugger


def load_spec_class(spectrograph=None, data_file=None):
    if spectrograph is not None:
