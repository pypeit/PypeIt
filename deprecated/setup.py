# !usr/bin/env python
# -*- coding: utf-8 -*-
#
# Licensed under a 3-clause BSD license.

import sys
import os
import glob

from setuptools import setup, find_packages

from extension_helpers import get_extensions

def get_data_files():
    """ Build the list of data files to include.  """
    data_files = []

    # Walk through the data directory, adding all files
    data_generator = os.walk('pypeit/data')
    for path, directories, files in data_generator:
        for f in files:
            data_path = '/'.join(path.split('/')[1:])
            data_files.append(os.path.join(data_path, f))

    # Add pipeline and spectrograph settings
    settings = glob.glob('pypeit/settings/settings.*')
    settings = ['/'.join(path.split('/')[1:]) for path in settings]
    data_files.extend(settings)

    return data_files


def get_scripts():
    """ Grab all the scripts in the bin directory.  """
    scripts = []
    if os.path.isdir('bin'):
        scripts = [ fname for fname in glob.glob(os.path.join('bin', '*'))
                                if not os.path.basename(fname).endswith('.rst') ]
    return scripts


def get_requirements():
    """ Get the requirements from a system file.  """
    name = 'pypeit/requirements.txt'

    requirements_file = os.path.join(os.path.dirname(__file__), name)
    install_requires = [line.strip().replace('==', '>=') for line in open(requirements_file)
                        if not line.strip().startswith('#') and line.strip() != '']
    return install_requires


NAME = 'pypeit'
# do not use x.x.x-dev.  things complain.  instead use x.x.xdev
VERSION = '1.3.1dev'
RELEASE = 'dev' not in VERSION

# To enable pypeit ginga global plugin
entry_points = """
[ginga.rv.plugins]
SlitWavelength = pypeit.display:setup_SlitWavelength
"""

def run_setup(data_files, scripts, packages, install_requires):

    # TODO: Are any/all of the *'d keyword arguments needed? I.e., what
    # are the default values?

    setup(name=NAME,
          provides=NAME,                                                # *
          version=VERSION,
          license='BSD3',
          description='PypeIt Spectroscopic Reduction',
          long_description=open('README.md').read(),
          author='PypeIt Collaboration',
          author_email='pypeit@ucolick.org',
          keywords='pypeit PypeIt astronomy Keck UCO Lick data reduction',
          url='https://github.com/pypeit/PypeIt',
          packages=packages,
          package_data={'pypeit': data_files, '': ['*.rst', '*.txt']},
          python_requires='>=3.7',
          include_package_data=True,
          scripts=scripts,
          install_requires=install_requires,
          requires=[ 'Python (>3.7.0)' ],                               # *
          zip_safe=False,                                               # *
          use_2to3=False,                                               # *
          setup_requires=[ 'pytest-runner' ],
          tests_require=[ 'pytest' ],
          ext_modules=get_extensions(),
          entry_points=entry_points,
          classifiers=[
              'Development Status :: 4 - Beta',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: BSD License',
              'Natural Language :: English',
              'Operating System :: OS Independent',
              'Programming Language :: Python',
              'Programming Language :: Python :: 3.7',
              'Topic :: Documentation :: Sphinx',
              'Topic :: Scientific/Engineering :: Astronomy',
              'Topic :: Software Development :: Libraries :: Python Modules',
              'Topic :: Software Development :: User Interfaces'
          ])

#-----------------------------------------------------------------------
if __name__ == '__main__':

    # Compile the data files to include
    data_files = get_data_files()
    # Compile the scripts in the bin/ directory
    scripts = get_scripts()
    # Get the packages to include
    packages = find_packages()
    # Collate the dependencies based on the system text file
    install_requires = get_requirements()
    install_requires = []  # Remove this line to enforce actual installation
    # Run setup from setuptools
    run_setup(data_files, scripts, packages, install_requires)


