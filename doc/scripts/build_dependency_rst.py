"""
Construct an rst table with the dependencies
"""

import tomllib
from importlib import resources

import numpy

from pypeit.utils import string_table

from IPython import embed


def write_dependency_table(project_file, path):
    ofile = path / 'dependencies_table.rst'

    with open(project_file, 'rb') as f:
        project = tomllib.load(f)

    user_requires = numpy.sort(project['project']['dependencies']) 
    dev_requires = numpy.sort(project['project']['optional-dependencies']['dev'])
    required_python = project['project']['requires-python']

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

    project_file = pypeit_root / 'pyproject.toml'
    if not project_file.is_file():
        raise FileNotFoundError(f'{project_file} does not exist!')
    write_dependency_table(project_file, output_root)


if __name__ == '__main__':
    main()

