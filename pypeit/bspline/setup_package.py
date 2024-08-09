import os
import sys

from setuptools import Extension
from extension_helpers import add_openmp_flags_if_available


C_BSPLINE_PKGDIR = os.path.relpath(os.path.dirname(__file__))

SRC_FILES = [os.path.join(C_BSPLINE_PKGDIR, filename)
              for filename in ['src/bspline.c']]

extra_compile_args=['-UNDEBUG', '-ffast-math']
if not sys.platform.startswith('win'):
    extra_compile_args.append('-fPIC')


def get_extensions():
    extension = Extension(name='pypeit.bspline._bspline', sources=SRC_FILES,
                          extra_compile_args=extra_compile_args, language='c',                      
                          export_symbols=['bspline_model', 'solution_arrays',
                                          'cholesky_band', 'cholesky_solve',
                                          'intrv'])
    
    # extension_helpers will check for opnmp support by trying to build
    # some test code with the appropriate flag. openmp provides a big performance boost, but some
    # systems, notably apple's version of clang that xcode provides, don't support it out of the box.

    add_openmp_flags_if_available(extension)
    return [extension]
