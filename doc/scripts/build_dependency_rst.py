"""
Construct an rst table with the dependencies
"""

import configparser
from importlib import resources

import numpy

from pypeit.utils import string_table

from IPython import embed


def write_dependency_table(setup_file, path):
    ofile = path / 'dependencies_table.rst'

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
    pypeit_root = resources.files('pypeit').parent
    output_root = pypeit_root / 'doc' / 'include'
    if not output_root.is_dir():
        raise NotADirectoryError(f'{output_root} does not exist!')

    setup_file = pypeit_root / 'setup.cfg'
    if not setup_file.is_file():
        raise FileNotFoundError(f'{setup_file} does not exist!')
    write_dependency_table(setup_file, output_root)


if __name__ == '__main__':
    main()

