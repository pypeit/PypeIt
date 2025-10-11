"""
Dynamically build table listing available standard stars.
"""

from importlib import resources

import numpy as np
from astropy import table

from pypeit.utils import string_table
from pypeit import dataPaths
from pypeit.core import standard

from IPython import embed


def write_tables(path):

    archives = list(standard.get_archive_sets()) + ['blackbody']
    files = [dataPaths.standards.get_file_path(f'{archive}/{archive}_info.txt') 
             for archive in archives]
    archives += ['kurucz93']
    files += [dataPaths.standards.get_file_path('kurucz93/schmidt-kaler_table.txt')]

    for archive, f in zip(archives,files):
        # Output file
        ofile = path / f'{archive}_table.rst'

        # get table data
        tbl = table.Table.read(f, comment='#', format='ascii')

        # create the table for the rst file
        data_table = np.empty((len(tbl)+1, len(tbl.keys())), dtype=object)
        data_table[0,:] = tbl.keys()

        for i, row in enumerate(tbl):
            data_table[i+1,:] = [str(v) for v in row]

        lines = string_table(data_table, delimeter='rst')
        with open(ofile, 'w') as f:
            f.write(lines)
        print(f'Wrote: {ofile}')


if __name__ == '__main__':
    output_root = resources.files('pypeit').parent / 'doc' / 'include'
    if not output_root.is_dir():
        output_root.mkdir(parents=True)

    write_tables(output_root)
