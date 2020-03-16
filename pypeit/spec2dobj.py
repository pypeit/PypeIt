"""
Module for the Spec2DObj class

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import copy
import inspect
from IPython import embed

import numpy as np

from scipy import interpolate

from astropy import units
from astropy.table import Table

from linetools.spectra import xspectrum1d

from pypeit import msgs
from pypeit.core import parse
from pypeit.core import flux_calib
from pypeit import utils
from pypeit import datamodel
from pypeit.images import detector_container

naming_model = {}
for skey in ['SPAT', 'SLIT', 'DET', 'SCI','OBJ', 'ORDER']:
    naming_model[skey.lower()] = skey

class Spec2DObj(datamodel.DataContainer):
    """Class to handle 2D spectral image outputs of PypeIt

    One generates one of these Objects for each detector in the exposure.

    Args:

    Attributes:
        See datamodel

    """
    version = '1.0.0'

    datamodel = {
        'sciiimg': dict(otype=np.ndarray, atype=np.floating, desc='2D processed science image'),
        'ivarraw': dict(otype=np.ndarray, atype=np.floating, desc='2D processed inverse variance image'),
        'skymodel': dict(otype=np.ndarray, atype=np.floating, desc='2D sky model image'),
        'objmodel': dict(otype=np.ndarray, atype=np.floating, desc='2D object model image'),
        'ivarmodel': dict(otype=np.ndarray, atype=np.floating, desc='2D ivar model image'),
        'mask': dict(otype=np.ndarray, atype=np.int32, desc='2D mask image'),
        'detector': dict(otype=detector_container.DetectorContainer, desc='Detector DataContainer'),
        'det': dict(otype=int, desc='Detector index'),
    }

    def __init__(self, det, sciimg, ivarraw, skymodel, objmodel, ivarmodel, mask, detector):

        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=_d)

    def _init_internals(self):
        pass


class AllSpec2DObj(dict):
    """
    Simple object to hold Spec2DObj objects
    and I/O

    Restrict keys to be type int or 'meta'
    and items to be `Spec2DObj`_

    """

    def __init__(self):
        self['meta'] = {}

    def __setitem__(self, item, value):
        # Check item
        if not isinstance(item, int) and item != 'meta':
            raise KeyError('Key must be an integer, i.e. detector number or "meta"')
        # Check value
        if isinstance(item, int):
            assert isinstance(value, Spec2DObj), 'Item must be a Spec2DObj'
        self.__dict__[item] = value
