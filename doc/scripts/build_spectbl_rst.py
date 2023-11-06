"""
Dynamically build table listing available spectrographs.
"""

from importlib import resources
import time

import numpy

from pypeit.utils import to_string, string_table
from pypeit.spectrographs import spectrograph_classes

from IPython import embed

#-----------------------------------------------------------------------------

def write_spec_table(path):
    ofile = path / 'spectrographs_table.rst'

    spec = spectrograph_classes()
    nspec = len(spec.keys())

    data_table = numpy.empty((nspec+1, 9), dtype=object)
    data_table[0,:] = ['``PypeIt`` Name', '``PypeIt`` Class', 'Telescope', 'Camera', 'URL',
                       'Pipeline', 'Supported', 'QL Tested', 'Comments']
    for i,cls in enumerate(spec.values()):
        data_table[i+1,0] = cls.name
        data_table[i+1,1] = ':class:`~' + cls.__module__ + '.' + cls.__name__ + '`'
        data_table[i+1,2] = cls.telescope['name']
        data_table[i+1,3] = cls.camera
        data_table[i+1,4] = '' if cls.url is None else f'`Link <{cls.url}>`__'
        data_table[i+1,5] = cls.pypeline
        data_table[i+1,6] = to_string(cls.supported)
        data_table[i+1,7] = to_string(cls.ql_supported)
        data_table[i+1,8] = '' if cls.comment is None else cls.comment

    lines = string_table(data_table, delimeter='rst')
    with open(ofile, 'w') as f:
        f.write(lines)
    print('Wrote: {}'.format(ofile))


if __name__ == '__main__':
    t = time.perf_counter()

    output_root = resources.files('pypeit').parent / 'doc' / 'include'
    if not output_root.is_dir():
        output_root.mkdir(parents=True)

    write_spec_table(output_root)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



