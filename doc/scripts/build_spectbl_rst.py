"""
Dynamically build table listing available spectrographs.
"""

import os
import time
from pkg_resources import resource_filename

from IPython import embed

import numpy

from pypeit.utils import to_string, string_table
from pypeit.spectrographs import spectrograph_classes

#-----------------------------------------------------------------------------

def write_spec_table(path):
    ofile = os.path.join(path, 'spectrographs_table.rst')

    spec = spectrograph_classes()
    nspec = len(spec.keys())

    data_table = numpy.empty((nspec+1, 7), dtype=object)
    data_table[0,:] = ['``PypeIt`` Name', '``PypeIt`` Class', 'Telescope', 'Camera',
                       'Pipeline Approach', 'Supported', 'Comments']
    for i,cls in enumerate(spec.values()):
        data_table[i+1,0] = cls.name
        data_table[i+1,1] = ':class:`~' + cls.__module__ + '.' + cls.__name__ + '`'
        data_table[i+1,2] = cls.telescope['name']
        data_table[i+1,3] = cls.camera
        data_table[i+1,4] = cls.pypeline
        data_table[i+1,5] = to_string(cls.supported)
        data_table[i+1,6] = '' if cls.comment is None else cls.comment

    lines = string_table(data_table, delimeter='rst')
    with open(ofile, 'w') as f:
        f.write(lines)
    print('Wrote: {}'.format(ofile))


if __name__ == '__main__':
    t = time.perf_counter()

    pypeit_root = os.path.dirname(resource_filename('pypeit', ''))
    path = os.path.join(pypeit_root, 'doc', 'include')
    if not os.path.isdir(path):
        os.makedirs(path)

    write_spec_table(path)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



