# encoding: utf-8
"""
Define the telescopes parameters used by Pypit.
"""
from pypeit.par.pypeitpar import TelescopePar

#TODO: Remove 'Par' from class name?

class GTCTelescopePar(TelescopePar):
    def __init__(self):
        super(GTCTelescopePar, self).__init__(name='GTC',
                                               longitude=17.877,
                                               latitude=28.762,
                                               elevation=2348.0,
                                               eff_aperture=73.0)


# eff_aperture of Keck, Shane from xidl
class KeckTelescopePar(TelescopePar):
    def __init__(self):
        super(KeckTelescopePar, self).__init__(name='KECK',
                                               longitude=155.47833,
                                               latitude=19.82833,
                                               elevation=4160.0,
                                               fratio=15,
                                               diameter=10,
                                               eff_aperture=72.3674)


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
                                                elevation=1283.0,
                                                diameter=3.05,
                                                eff_aperture=6.3617)

                                
class WHTTelescopePar(TelescopePar):
    def __init__(self):
        super(WHTTelescopePar, self).__init__(name='WHT',
                                              longitude=17.8947,
                                              latitude=26.7636,
                                              elevation=2396.0,
                                              diameter=4.2)
                                
class APFTelescopePar(TelescopePar):
    def __init__(self):
        super(APFTelescopePar, self).__init__(name='APF',
                                              longitude=121.642778,
                                              latitude=37.34138889,
                                              elevation=1283.0,
                                              diameter=2.4)
                                
class TNGTelescopePar(TelescopePar):
    def __init__(self):
        super(TNGTelescopePar, self).__init__(name='TNG',
                                              longitude=17.88906,
                                              latitude=28.754,
                                              elevation=2387.2,
                                              diameter=3.58)

class VLTTelescopePar(TelescopePar):
    def __init__(self):
        super(VLTTelescopePar, self).__init__(name='VLT',
                                              longitude=70.404830556,
                                              latitude=-24.6271666666,
                                              elevation=2635.43,
                                              diameter=8.2,
                                              eff_aperture=51.2)
# VLT aperture from https://www.eso.org/observing/etc/doc/formulabook/node15.html
# This seems unrealistic given that pi(8.2^2)/4 = 52.81

class NTTTelescopePar(TelescopePar):
    def __init__(self):
        super(NTTTelescopePar, self).__init__(name='NTT',
                                              longitude=289.2700,
                                              latitude=-29.2567,
                                              elevation=2375,
                                              diameter=3.58)

class GeminiNTelescopePar(TelescopePar):
    def __init__(self):
        super(GeminiNTelescopePar, self).__init__(name='GEMINI-N',
                                               longitude=155.47833,
                                               latitude=19.82833,
                                               elevation=4160.0,
                                               diameter=8.1)
class GeminiSTelescopePar(TelescopePar):
    def __init__(self):
        super(GeminiSTelescopePar, self).__init__(name='GEMINI-S',
                                              longitude=70.7367,              # Longitude of the telescope (NOTE: West should correspond to positive longitudes)
                                              latitude=-30.24075,              # Latitude of the telescope
                                              elevation=2750.0,               # Elevation of the telescope (in m)
                                              diameter=8.1)

class SOARTelescopePar(TelescopePar):
    def __init__(self):
        super(SOARTelescopePar, self).__init__(name='SOAR',
                                              longitude=70.7336,              # Longitude of the telescope (NOTE: West should correspond to positive longitudes)
                                              latitude=-30.2379,              # Latitude of the telescope
                                              elevation=2713.0,               # Elevation of the telescope (in m)
                                              diameter=4.1)                   # Ignores central obscuration

class LBTTelescopePar(TelescopePar):
    def __init__(self):
        super(LBTTelescopePar, self).__init__(name='LBT',
                                              longitude=109.889064,
                                              latitude=32.701308,
                                              elevation=3221.0,
                                              diameter=8.4)

class KPNOTelescopePar(TelescopePar):
    def __init__(self):
        super(KPNOTelescopePar, self).__init__(name='KPNO',
                                              longitude=111.616111,
                                              latitude=31.9516666,
                                              elevation=2098.0,
                                              diameter=4.0,
                                              eff_aperture=11.2)
# KPNO from https://en.wikipedia.org/wiki/Nicholas_U._Mayall_Telescope

class MMTTelescopePar(TelescopePar):
    def __init__(self):
        super(MMTTelescopePar, self).__init__(name='MMT',
                                              longitude=110.885,
                                              latitude=31.6883,
                                              elevation=2616.0,
                                              diameter=6.5)

class NOTTelescopePar(TelescopePar):
    def __init__(self):
        super(NOTTelescopePar, self).__init__(name='NOT',
                                              longitude=17.88432979,
                                              latitude=28.7543303,
                                              elevation=2465.5,
                                              diameter=2.56)

class P200TelescopePar(TelescopePar):
    def __init__(self):
        super(P200TelescopePar, self).__init__(name='P200',
                                               longitude=116.86489,
                                               latitude=33.35631,
                                               elevation=1713.,
                                               diameter=5.1)

class BokTelescopePar(TelescopePar):
    def __init__(self):
        super(BokTelescopePar, self).__init__(name='BOK',
                                              longitude=111.6004,
                                              latitude=31.9629,
                                              elevation=2071.1)

class LDTTelescopePar(TelescopePar):
    def __init__(self):
        super(LDTTelescopePar, self).__init__(name='LDT',
                                              longitude=111.4223,
                                              latitude=34.7443,
                                              elevation=2361.0,
                                              fratio=6.1,
                                              diameter=4.3,
                                              eff_aperture=49.5)

