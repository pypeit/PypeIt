
from pypeit import spectrographs
import numpy

cls = spectrographs.spectrograph_classes()
keys = []
for k, v in cls.items():
    keys += v().configuration_keys()

print('Unique configuration keys: ', numpy.unique(keys))

