
from pypeit.spectrographs.util import spectrograph_classes
import numpy

cls = spectrograph_classes()
keys = []
for k, v in cls.items():
    keys += v().configuration_keys()

print('Unique configuration keys: ', numpy.unique(keys))

