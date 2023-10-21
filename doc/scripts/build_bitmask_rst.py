"""
Dynamically build tables for the bitmasks.
"""

from importlib import resources
import time

import numpy

from pypeit.utils import to_string, string_table

from IPython import embed

#-----------------------------------------------------------------------------

def write_bitmask_table(obj, path):
    bm = obj()
    ofile = path / f'{obj.__name__.lower()}_table.rst'

    data_table = numpy.empty((bm.nbits+1, 4), dtype=object)
    data_table[0,:] = ['Bit Name', 'Bit Number', 'Decimal Value', 'Description']
    for i,k in enumerate(bm.bits.keys()):
        data_table[i+1,0] = k
        data_table[i+1,1] = to_string(bm.bits[k])
        data_table[i+1,2] = to_string(int(2**bm.bits[k]))
        data_table[i+1,3] = to_string(bm.descr[bm.bits[k]])

    lines = string_table(data_table, delimeter='rst')
    with open(ofile, 'w') as f:
        f.write(lines)
    print('Wrote: {}'.format(ofile))


if __name__ == '__main__':
    t = time.perf_counter()

    path = resources.files('pypeit').parent / 'doc' / 'include'
    if not path.is_dir():
        path.mkdir(parents=True)

    from pypeit.images.imagebitmask import ImageBitMask
    write_bitmask_table(ImageBitMask, path)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



