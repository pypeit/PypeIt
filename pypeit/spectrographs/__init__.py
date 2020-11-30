
from pypeit.spectrographs import spectrograph

# The import of all the spectrograph modules here is what enables the dynamic
# compiling of all the available spectrographs below
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
from pypeit.spectrographs import mmt_bluechannel
from pypeit.spectrographs import mmt_mmirs
from pypeit.spectrographs import not_alfosc
from pypeit.spectrographs import p200_dbsp
from pypeit.spectrographs import p200_tspec
from pypeit.spectrographs import shane_kast
from pypeit.spectrographs import tng_dolores
from pypeit.spectrographs import vlt_fors
from pypeit.spectrographs import vlt_xshooter
from pypeit.spectrographs import wht_isis

# Build the list of names for the available spectrographs
import numpy as np

def all_subclasses(cls):
    """
    Thanks to:
    https://stackoverflow.com/questions/3862310/how-to-find-all-the-subclasses-of-a-class-given-its-name
    """
    return set(cls.__subclasses__()).union(
        [s for c in cls.__subclasses__() for s in all_subclasses(c)])

def spectrograph_classes():
    # Recursively collect all subclasses
    spec_c = np.array(list(all_subclasses(spectrograph.Spectrograph)))
    # Select spectrograph classes with a defined name; spectrographs without a
    # name are either undefined or a base class.
    spec_c = spec_c[[c.name is not None for c in spec_c]]
    # Construct a dictionary with the spectrograph name and class
    srt = np.argsort(np.array([c.name for c in spec_c]))
    return dict([ (c.name,c) for c in spec_c[srt]])

available_spectrographs = list(spectrograph_classes().keys())

