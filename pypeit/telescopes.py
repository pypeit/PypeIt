# encoding: utf-8
"""
Define the telescopes parameters used by Pypit.
"""
from pypeit.par.pypeitpar import TelescopePar

#TODO: Remove 'Par' from class name?

class GeminiNTelescopePar(TelescopePar):
    def __init__(self):
        super(GeminiNTelescopePar, self).__init__(name='GEMINI-N',
                                               longitude=155.47833,
                                               latitude=19.82833,
                                               elevation=4160.0)

class KeckTelescopePar(TelescopePar):
    def __init__(self):
        super(KeckTelescopePar, self).__init__(name='KECK',
                                               longitude=155.47833,
                                               latitude=19.82833,
                                               elevation=4160.0,
                                               fratio=15,
                                               diameter=10)

class MagellanTelescopePar(TelescopePar):
    def __init__(self):
        super(MagellanTelescopePar, self).__init__(name='MAGELLAN',
                                               longitude=70.6858,
                                               latitude=-29.0283,
                                               elevation=2516.0,
                                               diameter=6.5)

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
class GeminiSTelescopePar(TelescopePar):
    def __init__(self):
        super(GeminiSTelescopePar, self).__init__(name='GEMINI-S',
                                              longitude=70.8062,              # Longitude of the telescope (NOTE: West should correspond to positive longitudes)
                                              latitude=-30.1691,              # Latitude of the telescope
                                              elevation=2200.0)               # Elevation of the telescope (in m)

class LBTTelescopePar(TelescopePar):
    def __init__(self):
        super(LBTTelescopePar, self).__init__(name='LBT',
                                               longitude=109.889064,
                                               latitude=32.701308,
                                               elevation=3221.0)

class MMTTelescopePar(TelescopePar):
    def __init__(self):
        super(MMTTelescopePar, self).__init__(name='MMT',
                                               longitude=110.885,
                                               latitude=31.6883,
                                               elevation=2616.0)
