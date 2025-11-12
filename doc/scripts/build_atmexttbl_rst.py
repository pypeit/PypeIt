"""
Dynamically build table listing available standard stars.
"""

from importlib import resources

import numpy as np
from astropy import table

from pypeit.utils import string_table
from pypeit import dataPaths

from IPython import embed


def write_tables(path):

    # Output file
    ofile = path / 'atmext_table.rst'

    # get table data
    file = dataPaths.extinction.get_file_path(f'extinction_curves.txt')
    tbl = table.Table.read(file, comment='#', format='ascii')

    # create the table for the rst file
    data_table = np.empty((len(tbl)+1, len(tbl.keys())), dtype=object)
    data_table[0,:] = tbl.keys()

    for i, row in enumerate(tbl):
        data_table[i+1,:] = [str(r) for r in row]

    lines = string_table(data_table, delimeter='rst')
    with open(ofile, 'w') as f:
        f.write(lines)
    print(f'Wrote: {ofile}')



if __name__ == '__main__':
    output_root = resources.files('pypeit').parent / 'doc' / 'include'
    if not output_root.is_dir():
        output_root.mkdir(parents=True)

    write_tables(output_root)
