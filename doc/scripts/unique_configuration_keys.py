
from pypeit import spectrographs
import numpy as np

cls = spectrographs.spectrograph_classes()
keys = []
for k, v in cls.items():
    keys += v().configuration_keys()

print('Unique configuration keys: ', np.unique(keys))

