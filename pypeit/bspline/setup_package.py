import os
import sys
import tempfile, subprocess, shutil

from distutils.extension import Extension

# the most reliable way to check for openmp support in the C compiler is to try to build
# some test code with the -fopenmp flag. openmp provides a big performance boost, but some
# systems, notably apple's version of clang that xcode provides, don't support it out of the box.

# see http://openmp.org/wp/openmp-compilers/
omp_test = \
r"""
#include <omp.h>
#include <stdio.h>
int main() {
#pragma omp parallel
printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
}
"""

def check_for_openmp():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    if 'CC' in os.environ:
        c_compiler = os.environ['CC']
    else:
        c_compiler = 'gcc'

    filename = r'test.c'
    with open(filename, 'w') as file:
        file.write(omp_test)
    with open(os.devnull, 'w') as fnull:
        result = subprocess.call([c_compiler, '-fopenmp', filename],
                                 stdout=fnull, stderr=fnull)

    os.chdir(curdir)
    # clean up test code
    shutil.rmtree(tmpdir)

    # return code from compiler process is 0 if it completed successfully
    if result == 0:
        return True
    else:
        return False

C_BSPLINE_PKGDIR = os.path.relpath(os.path.dirname(__file__))

SRC_FILES = [os.path.join(C_BSPLINE_PKGDIR, filename)
              for filename in ['src/bspline.c']]

extra_compile_args=['-UNDEBUG', '-ffast-math']
if not sys.platform.startswith('win'):
    extra_compile_args.append('-fPIC')

if check_for_openmp():
    extra_compile_args.append('-fopenmp')
    extra_link_args = ['-fopenmp']
else:
    extra_link_args = []

def get_extensions():
    return [Extension(name='pypeit.bspline._bspline', sources=SRC_FILES,
                      extra_compile_args=extra_compile_args, language='c',
                      extra_link_args=extra_link_args,
                      export_symbols=['bspline_model', 'solution_arrays',
                                      'cholesky_band', 'cholesky_solve',
                                      'intrv'])]
