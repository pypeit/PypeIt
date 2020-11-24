
from pypeit.spectrographs import spectrograph

from pypeit.spectrographs import gemini_flamingos
from pypeit.spectrographs import gemini_gmos
from pypeit.spectrographs import gemini_gnirs
from pypeit.spectrographs import keck_deimos
from pypeit.spectrographs import keck_hires
from pypeit.spectrographs import keck_kcwi
from pypeit.spectrographs import keck_lris
from pypeit.spectrographs import keck_mosfire
from pypeit.spectrographs import keck_nires
from pypeit.spectrographs import keck_nirspec
from pypeit.spectrographs import lbt_luci
from pypeit.spectrographs import lbt_mods
from pypeit.spectrographs import magellan_fire
from pypeit.spectrographs import magellan_mage
from pypeit.spectrographs import mdm_osmos
from pypeit.spectrographs import mmt_binospec
from pypeit.spectrographs import mmt_mmirs
from pypeit.spectrographs import not_alfosc
from pypeit.spectrographs import p200_dbsp
from pypeit.spectrographs import p200_tspec
from pypeit.spectrographs import shane_kast
from pypeit.spectrographs import tng_dolores
from pypeit.spectrographs import vlt_fors
from pypeit.spectrographs import vlt_xshooter
from pypeit.spectrographs import wht_isis

# Build the names of the supported spectrographs
import numpy as np

def spectrograph_subclasses():
    subc = np.array(spectrograph.Spectrograph.__subclasses__())
    nss = np.array([len(c.__subclasses__()) for c in subc])
    while np.any(nss > 0):
        add_subc = np.empty(0, dtype=object)
        keep = np.ones(subc.size, dtype=bool)
        for i in range(subc.size):
            if nss[i] == 0:
                continue
            keep[i] = False
            add_subc = np.append(add_subc, subc[i].__subclasses__())
        subc = np.append(subc[keep], add_subc)
        nss = np.array([len(c.__subclasses__()) for c in subc])
    return subc

supported_spectrographs = np.sort([c.name for c in spectrograph_subclasses()]).tolist()
