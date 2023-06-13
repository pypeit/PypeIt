# encoding: utf-8
"""
Define the telescopes parameters used by Pypit.

NOTE: Longitudes are measured increasing to the east, so west longitudes are negative.
"""
from pypeit.par.pypeitpar import TelescopePar
from astropy.coordinates import EarthLocation
from astropy import units

#TODO: Remove 'Par' from class name?

class GTCTelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('Roque de los Muchachos')
        super(GTCTelescopePar, self).__init__(name='GTC',
                                               longitude=loc.lon.to(units.deg).value,
                                               latitude=loc.lat.to(units.deg).value,
                                               elevation=loc.height.to(units.m).value,
                                               eff_aperture=73.0)


# eff_aperture of Keck, Shane from xidl
class KeckTelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('W. M. Keck Observatory')
        super(KeckTelescopePar, self).__init__(name='KECK',
                                               longitude=loc.lon.to(units.deg).value,
                                               latitude=loc.lat.to(units.deg).value,
                                               elevation=loc.height.to(units.m).value,
                                               fratio=15,
                                               diameter=10,
                                               eff_aperture=72.3674)


class MagellanTelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('Las Campanas Observatory')
        super(MagellanTelescopePar, self).__init__(name='MAGELLAN',
                                               longitude=loc.lon.to(units.deg).value,
                                               latitude=loc.lat.to(units.deg).value,
                                               elevation=loc.height.to(units.m).value,
                                               diameter=6.5)

class ShaneTelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('Lick Observatory')
        super(ShaneTelescopePar, self).__init__(name='SHANE',
                                                longitude=loc.lon.to(units.deg).value,
                                                latitude=loc.lat.to(units.deg).value,
                                                elevation=loc.height.to(units.m).value,
                                                diameter=3.05,
                                                eff_aperture=6.3617)

                                
class WHTTelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('Roque de los Muchachos')
        super(WHTTelescopePar, self).__init__(name='WHT',
                                              longitude=loc.lon.to(units.deg).value,
                                              latitude=loc.lat.to(units.deg).value,
                                              elevation=loc.height.to(units.m).value,
                                              diameter=4.2)
                                
class APFTelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('Lick Observatory')
        super(APFTelescopePar, self).__init__(name='APF',
                                              longitude=loc.lon.to(units.deg).value,
                                              latitude=loc.lat.to(units.deg).value,
                                              elevation=loc.height.to(units.m).value,
                                              diameter=2.4)
                                
class TNGTelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('Roque de los Muchachos')
        super(TNGTelescopePar, self).__init__(name='TNG',
                                              longitude=loc.lon.to(units.deg).value,
                                              latitude=loc.lat.to(units.deg).value,
                                              elevation=loc.height.to(units.m).value,
                                              diameter=3.58)

class VLTTelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('Paranal Observatory')
        super(VLTTelescopePar, self).__init__(name='VLT',
                                              longitude=loc.lon.to(units.deg).value,
                                              latitude=loc.lat.to(units.deg).value,
                                              elevation=loc.height.to(units.m).value,
                                              diameter=8.2,
                                              eff_aperture=51.2)
# VLT aperture from https://www.eso.org/observing/etc/doc/formulabook/node15.html
# This seems unrealistic given that pi(8.2^2)/4 = 52.81

class NTTTelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('La Silla Observatory')
        super(NTTTelescopePar, self).__init__(name='NTT',
                                              longitude=loc.lon.to(units.deg).value,
                                              latitude=loc.lat.to(units.deg).value,
                                              elevation=loc.height.to(units.m).value,
                                              diameter=3.58)

class GeminiNTelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('gemini_north')
        super(GeminiNTelescopePar, self).__init__(name='GEMINI-N',
                                               longitude=loc.lon.to(units.deg).value,
                                               latitude=loc.lat.to(units.deg).value,
                                               elevation=loc.height.to(units.m).value,
                                               diameter=8.1)
class GeminiSTelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('gemini_south')
        super(GeminiSTelescopePar, self).__init__(name='GEMINI-S',
                                                  longitude=loc.lon.to(units.deg).value,
                                                  latitude=loc.lat.to(units.deg).value,
                                                  elevation=loc.height.to(units.m).value,
                                                  diameter=8.1)

class SOARTelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('Cerro Pachon')
        super(SOARTelescopePar, self).__init__(name='SOAR',
                                               longitude=loc.lon.to(units.deg).value,
                                               latitude=loc.lat.to(units.deg).value,
                                               elevation=loc.height.to(units.m).value,
                                               diameter=4.1)                   # Ignores central obscuration

class LBTTelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('Large Binocular Telescope')
        super(LBTTelescopePar, self).__init__(name='LBT',
                                              longitude=loc.lon.to(units.deg).value,
                                              latitude=loc.lat.to(units.deg).value,
                                              elevation=loc.height.to(units.m).value,
                                              diameter=8.4)

class KPNOTelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('Kitt Peak National Observatory')
        super(KPNOTelescopePar, self).__init__(name='KPNO',
                                               longitude=loc.lon.to(units.deg).value,
                                               latitude=loc.lat.to(units.deg).value,
                                               elevation=loc.height.to(units.m).value,
                                               diameter=4.0,
                                               eff_aperture=11.2)
# KPNO from https://en.wikipedia.org/wiki/Nicholas_U._Mayall_Telescope

class MMTTelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('Whipple Observatory')
        super(MMTTelescopePar, self).__init__(name='MMT',
                                              longitude=loc.lon.to(units.deg).value,
                                              latitude=loc.lat.to(units.deg).value,
                                              elevation=loc.height.to(units.m).value,
                                              diameter=6.5)

class NOTTelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('Roque de los Muchachos')
        super(NOTTelescopePar, self).__init__(name='NOT',
                                              longitude=loc.lon.to(units.deg).value,
                                              latitude=loc.lat.to(units.deg).value,
                                              elevation=loc.height.to(units.m).value,
                                              diameter=2.56)

class P200TelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('Palomar')
        super(P200TelescopePar, self).__init__(name='P200',
                                               longitude=loc.lon.to(units.deg).value,
                                               latitude=loc.lat.to(units.deg).value,
                                               elevation=loc.height.to(units.m).value,
                                               diameter=5.1)

class BokTelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('Kitt Peak')
        super(BokTelescopePar, self).__init__(name='BOK',
                                              longitude=loc.lon.to(units.deg).value,
                                              latitude=loc.lat.to(units.deg).value,
                                              elevation=loc.height.to(units.m).value)

class LDTTelescopePar(TelescopePar):
    def __init__(self):
        loc = EarthLocation.of_site('Discovery Channel Telescope')
        super(LDTTelescopePar, self).__init__(name='LDT',
                                              longitude=loc.lon.to(units.deg).value,
                                              latitude=loc.lat.to(units.deg).value,
                                              elevation=loc.height.to(units.m).value,
                                              fratio=6.18,
                                              diameter=4.25,
                                              eff_aperture=12.37)


# TODO provisional values
class JWSTTelescopePar(TelescopePar):
    def __init__(self):
        super(JWSTTelescopePar, self).__init__(name='JWST',
                                              longitude=0.0,
                                              latitude=0.0,
                                              elevation=0.0,
                                              diameter=6.5)
