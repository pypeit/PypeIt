"""
Dynamically build table listing available standard stars.
"""

from importlib import resources

import astropy.table

from pypeit import dataPaths

from IPython import embed


def write_tables(path):

    # Output file
    ofile = path / 'atmext_table.rst'

    # get table data
    file = dataPaths.extinction.get_file_path(f'extinction_curves.txt')
    tbl = astropy.table.Table.read(file, comment='#', format='ascii')
    tbl.meta = {}

    tbl.write(ofile, format="ascii.rst", overwrite=True)
    print(f'Wrote: {ofile}')


if __name__ == '__main__':
    output_root = resources.files('pypeit').parent / 'doc' / 'include'
    if not output_root.is_dir():
        output_root.mkdir(parents=True)

    write_tables(output_root)
