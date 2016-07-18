'''
Module for dealing with multi-object spectroscopy slitmasks
'''
from __future__ import (print_function, absolute_import,
                        division, unicode_literals)

import numpy as np
from astropy.io import fits as pyfits

from pypit import armsgs

try:
    from xastropy.xutils import xdebug as debugger
except:
    import pdb as debugger

# Logging
msgs = armsgs.get_logger()

class Slit(object):
    '''
    Slit with mask coordinates
    '''

    def __init__(self, left_edge, right_edge, name):
        '''
        Parameters
        ----------
        left_edge : float, 
        right_edge : 
        name : str, name of slit
        '''
        self.left_edge = left_edge
        self.right_edge = right_edge
        self.name = name

    def __repr__(self):
        return '<' + self.__name__ + ': ' + self.name + '>'

        
class Slitmask(object):
    '''
    Generic slitmask class, should be sub-classed for specific instruments.
    '''

    def __init__(self):
        self.mask_name = None
        self.mask_pa = None
        self.slits = []

    def __len__(self):
        return len(self.slits)
    
    def __repr__(self):
        s = ('<' + self.__name__ + ': ' + self.mask_name +  ' (' +
             str(len(self)) + ' slits)>')

class DEIMOS_slitmask(Slitmask):
    pass

class LRIS_slitmask(Slitmask):
    pass
