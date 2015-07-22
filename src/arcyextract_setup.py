from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy, os

include_gsl_dir = os.getenv('GSL_PATH')+'/include/'
lib_gsl_dir = os.getenv('GSL_PATH')+'/lib/'

ext = Extension("arcyextract", ["arcyextract.pyx"],
    include_dirs=[numpy.get_include(),
                include_gsl_dir],
    library_dirs=[lib_gsl_dir],
    libraries=["gsl","gslcblas"]
)

setup(ext_modules=[ext],
        cmdclass = {'build_ext': build_ext})