# encoding: utf-8
"""
Define the telescopes parameters used by Pypit.
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

from pypeit.par.pypeitpar import TelescopePar

#TODO: Remove 'Par' from class name?

class KeckTelescopePar(TelescopePar):
    def __init__(self):
        super(KeckTelescopePar, self).__init__(name='KECK',
                                               longitude=155.47833,
                                               latitude=19.82833,
                                               elevation=4160.0)

class ShaneTelescopePar(TelescopePar):
    def __init__(self):
        super(ShaneTelescopePar, self).__init__(name='SHANE',
                                                longitude=121.6428,
                                                latitude=37.3413889,
                                                elevation=1283.0)
                                
class WHTTelescopePar(TelescopePar):
    def __init__(self):
        super(WHTTelescopePar, self).__init__(name='WHT',
                                              longitude=17.8947,
                                              latitude=26.7636,
                                              elevation=2396.0)
                                
class APFTelescopePar(TelescopePar):
    def __init__(self):
        super(APFTelescopePar, self).__init__(name='APF',
                                              longitude=121.642778,
                                              latitude=37.34138889,
                                              elevation=1283.0)
                                
class TNGTelescopePar(TelescopePar):
    def __init__(self):
        super(TNGTelescopePar, self).__init__(name='TNG',
                                              longitude=17.88906,
                                              latitude=28.754,
                                              elevation=2387.2)

class VLTTelescopePar(TelescopePar):
    def __init__(self):
        super(VLTTelescopePar, self).__init__(name='VLT',
                                               longitude=70.404830556,
                                               latitude=-24.6271666666,
                                               elevation=2635.43)
