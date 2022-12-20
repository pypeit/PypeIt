import os
import sys
from distutils.extension import Extension

C_BSPLINE_PKGDIR = os.path.relpath(os.path.dirname(__file__))

SRC_FILES = [os.path.join(C_BSPLINE_PKGDIR, filename)
              for filename in ['src/bspline.c']]

extra_compile_args=['-UNDEBUG']
if not sys.platform.startswith('win'):
    extra_compile_args.append('-fPIC')

def get_extensions():
    return [Extension(name='pypeit.bspline._bspline', sources=SRC_FILES,
                      extra_compile_args=extra_compile_args, language='c',
                      export_symbols=['bspline_model', 'solution_arrays',
                                      'cholesky_band', 'cholesky_solve',
                                      'intrv'])]
