"""
Dynamically build table listing available spectrographs.
"""

from importlib import resources
import time

import numpy as np
from astropy import table

from pypeit.utils import to_string, string_table
from pypeit import data

from IPython import embed


def write_table(path):
    ofile = path / 'calspec_table.rst'

    # get table data
    star_file = data.Paths.standards / 'calspec' / "calspec_info.txt"
    star_tbl = table.Table.read(star_file, comment='#', format='ascii')

    data_table = np.empty((len(star_tbl)+1, len(star_tbl.keys())), dtype=object)
    data_table[0,:] = star_tbl.keys()

    for i, row in enumerate(star_tbl):
        data_table[i+1,:] = list(row)

    lines = string_table(data_table, delimeter='rst')
    with open(ofile, 'w') as f:
        f.write(lines)
    print('Wrote: {}'.format(ofile))


if __name__ == '__main__':
    t = time.perf_counter()

    output_root = resources.files('pypeit').parent / 'doc' / 'include'
    if not output_root.is_dir():
        output_root.mkdir(parents=True)

    write_table(output_root)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



