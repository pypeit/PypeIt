"""
Construct an rst table with the dependencies
"""

from IPython import embed
import os
import configparser
from pkg_resources import resource_filename

import numpy

from pypeit.utils import string_table


def write_dependency_table(setup_file, path):
    ofile = os.path.join(path, 'dependencies_table.rst')

    setup = configparser.ConfigParser()
    setup.read(setup_file)

    user_requires = numpy.sort(setup['options']['install_requires'].split('\n')[1:])
    dev_requires = numpy.sort(setup['options.extras_require']['dev'].split('\n')[1:])
    required_python = setup['options']['python_requires']

    data_table = numpy.empty((3, 2), dtype=object)
    data_table[0,:] = ['Python Version', f'``{required_python}``']
    data_table[1,:] = ['Required for users', ', '.join([f'``{u}``' for u in user_requires])]
    data_table[2,:] = ['Required for developers', ', '.join([f'``{d}``' for d in dev_requires])]

    lines = string_table(data_table, delimeter='rst', has_header=False)
    with open(ofile, 'w') as f:
        f.write(lines)
    print('Wrote: {}'.format(ofile))


def main():
    output_root = os.path.join(os.path.split(os.path.abspath(resource_filename('pypeit', '')))[0],
                               'doc', 'include')
    if not os.path.isdir(output_root):
        raise NotADirectoryError(f'{output_root} does not exist!')

    setup_file = os.path.join(os.path.split(os.path.abspath(resource_filename('pypeit', '')))[0],
                               'setup.cfg')
    if not os.path.isfile(setup_file):
        raise FileNotFoundError(f'{setup_file} does not exist!')
    write_dependency_table(setup_file, output_root)


if __name__ == '__main__':
    main()

