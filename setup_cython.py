#!/usr/bin/env python
#
# python setup_cython.py build_ext --inplace
# Licensed under a 3-clause BSD style license - see LICENSE.rst
# conda install -c https://conda.anaconda.org/asmeurer gsl
from __future__ import absolute_import, division, print_function
#
# Standard imports
#
import glob
import numpy, os
from Cython.Distutils import build_ext
from distutils.extension import Extension
#
# setuptools' sdist command ignores MANIFEST.in
#
#from distutils.command.sdist import sdist as DistutilsSdist
from setuptools import setup
#
# DESI support code.
#
#from desiutil.setup import DesiTest, DesiVersion, get_version
#
# Begin setup
#
setup_keywords = dict()
#
# THESE SETTINGS NEED TO BE CHANGED FOR EVERY PRODUCT.
#
setup_keywords['name'] = 'pypit'
setup_keywords['description'] = 'PYPIT Spectroscopic Reduction'
setup_keywords['author'] = 'PYPIT Collaboration'
setup_keywords['author_email'] = 'pypit@ucolick.org'
setup_keywords['license'] = 'BSD'
setup_keywords['url'] = 'https://github.com/pypit/pypit'
#
# END OF SETTINGS THAT NEED TO BE CHANGED.
#
setup_keywords['version'] = '0.6.dev0' #get_version(setup_keywords['name'])
#
# Use README.rst as long_description.
#
setup_keywords['long_description'] = ''
if os.path.exists('README.md'):
    with open('README.md') as readme:
        setup_keywords['long_description'] = readme.read()
#
# Set other keywords for the setup function.  These are automated, & should
# be left alone unless you are an expert.
#
# Treat everything in bin/ except *.rst as a script to be installed.
#
if os.path.isdir('bin'):
    setup_keywords['scripts'] = [fname for fname in glob.glob(os.path.join('bin', '*'))
        if not os.path.basename(fname).endswith('.rst')]
setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['requires'] = ['Python (>2.7.0)']
# setup_keywords['install_requires'] = ['Python (>2.7.0)']
setup_keywords['zip_safe'] = False
setup_keywords['use_2to3'] = True
setup_keywords['packages'] = ['pypit'] #find_packages('pypit')
#setup_keywords['package_dir'] = {'':''}
#setup_keywords['cmdclass'] = {'version': DesiVersion, 'test': DesiTest, 'sdist': DistutilsSdist}
#setup_keywords['test_suite']='{name}.test.{name}_test_suite.{name}_test_suite'.format(**setup_keywords)

# Cython

include_gsl_dir = os.getenv('GSL_PATH')+'/include/'
lib_gsl_dir = os.getenv('GSL_PATH')+'/lib/'
pyx_files = glob.glob('pypit/*.pyx')
setup_keywords['ext_modules']=[]
for pyx_file in pyx_files:
    pyx_split = pyx_file.split('.')
    # Generate Extension
    ext = Extension(pyx_split[0], [pyx_file],
        include_dirs=[numpy.get_include(),
                    include_gsl_dir],
        library_dirs=[lib_gsl_dir],
        libraries=["gsl","gslcblas"]
    )
    # Append
    setup_keywords['ext_modules'].append(ext)
setup_keywords['cmdclass']={'build_ext': build_ext}
import pdb
pdb.set_trace()
#
# Autogenerate command-line scripts.
#
# setup_keywords['entry_points'] = {'console_scripts':['desiInstall = desiutil.install.main:main']}

#
# Add internal data directories.
#
setup_keywords['package_data'] = {'pypit': ['data/extinction/*',
                                            'data/arc_lines/*',
                                            'data/standards/*',
                                            'data/sky_spec/*',]}
#
# Run setup command.
#
setup(**setup_keywords)

