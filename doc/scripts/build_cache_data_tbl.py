"""
Construct an rst table listing the data directories and whether or not they're
expected to be hosted by the cache system.
"""
from importlib import resources

import numpy as np

from pypeit.pypeitdata import PypeItDataPaths
from pypeit.utils import string_table

from IPython import embed

def write_table(ofile):
    """
    Write a reStructureText file with a table listing the directories in
    pypeit/data, their ``pypeit.dataPaths`` attribute, and where the data is
    hosted.

    See ../installing.rst for more description.
    """
    # Table columns
    data_table = [['Reference', 'Subdirectory', 'Host']]
    # Use the class dictionary to fill the table
    for key, meta in PypeItDataPaths.defined_paths.items():
        data_table += [[key, meta['path'], '...' if meta['host'] is None else meta['host']]]

    # Convert to a numpy array and sort by the attribute name
    data_table = np.atleast_1d(data_table)
    srt = np.argsort(data_table[:,0])
    data_table = data_table[srt]

    # Convert the 2D array of values into text lines of an RST table
    lines = string_table(data_table, delimeter='rst')
    # Print the table
    with open(ofile, 'w') as f:
        f.write(lines)
    # Report
    print('Wrote: {}'.format(ofile))


def main():
    output_root = resources.files('pypeit').parent / 'doc' / 'include'
    if not output_root.is_dir():
        raise NotADirectoryError(f'{output_root} does not exist!')

    ofile = output_root / 'data_dir.rst'
    write_table(ofile)


if __name__ == '__main__':
    main()

