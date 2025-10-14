"""
Construct an rst table listing the data directories and whether or not they're
expected to be hosted by the cache system.
"""
from importlib import resources

import astropy.table

from pypeit.pypeitdata import PypeItDataPaths

from IPython import embed

def write_table(ofile):
    """
    Write a reStructureText file with a table listing the directories in
    pypeit/data, their ``pypeit.dataPaths`` attribute, and where the data is
    hosted.

    See ../installing.rst for more description.
    """
    data_table = []
    # Use the class dictionary to fill the table
    for key, meta in PypeItDataPaths.defined_paths.items():

        data_table.append(
            {
                'Reference': key,
                'Subdirectory': meta['path'],
                'Host': '...' if meta['host'] is None else meta['host'],
            }
        )

    tbl = astropy.table.Table(data_table)
    tbl.sort('Reference')

    tbl.write(ofile, format="ascii.rst", overwrite=True)
    print(f'Wrote: {ofile}')


def main():
    output_root = resources.files('pypeit').parent / 'doc' / 'include'
    if not output_root.is_dir():
        raise NotADirectoryError(f'{output_root} does not exist!')

    ofile = output_root / 'data_dir.rst'
    write_table(ofile)


if __name__ == '__main__':
    main()

