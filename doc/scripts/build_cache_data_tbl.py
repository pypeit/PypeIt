"""
Construct an rst table with the detector properties
"""
from importlib import resources

import numpy as np

from pypeit.pypeitdata import PypeItDataPaths
from pypeit.utils import string_table

from IPython import embed

def write_table(ofile):

    data_table = [['Reference', 'Subdirectory', 'Host']]
    for key, meta in PypeItDataPaths.defined_paths.items():
        data_table += [[key, meta['path'], '...' if meta['host'] is None else meta['host']]]

    data_table = np.atleast_1d(data_table)
    srt = np.argsort(data_table[:,0])
    data_table = data_table[srt]

    lines = string_table(data_table, delimeter='rst')
    with open(ofile, 'w') as f:
        f.write(lines)
    print('Wrote: {}'.format(ofile))


def main():
    output_root = resources.files('pypeit').parent / 'doc' / 'include'
    if not output_root.is_dir():
        raise NotADirectoryError(f'{output_root} does not exist!')

    ofile = output_root / 'data_dir.rst'
    write_table(ofile)


if __name__ == '__main__':
    main()

