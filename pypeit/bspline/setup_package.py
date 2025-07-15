import os
import sys

from setuptools import Extension
from extension_helpers import add_openmp_flags_if_available


C_BSPLINE_PKGDIR = os.path.relpath(os.path.dirname(__file__))

SRC_FILES = [os.path.join(C_BSPLINE_PKGDIR, filename) for filename in ["src/bspline.c"]]

extra_compile_args = ["-UNDEBUG", "-O3", "-ffast-math"]
extra_link_args = []

if not sys.platform.startswith("win"):
    extra_compile_args.append("-fPIC")

if "CONDA_PREFIX" in os.environ:
    # If we are in a conda environment, we need to add the include and lib directories
    # so that the libraries from the llvm-openmp package can be found.
    extra_compile_args.append("-I" + os.environ["CONDA_PREFIX"] + "/include")
    extra_link_args.append("-L" + os.environ["CONDA_PREFIX"] + "/lib")

if "darwin" in sys.platform:
    # On macOS, we need manually add the -Xclang and -fopenmp flags to unlock OpenMP.
    extra_compile_args.extend(["-Xclang", "-fopenmp"])

def get_extensions():
    extension = Extension(
        name="pypeit.bspline._bspline",
        sources=SRC_FILES,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        language="c",
        export_symbols=[
            "bspline_model",
            "solution_arrays",
            "cholesky_band",
            "cholesky_solve",
            "intrv",
        ],
    )

    # extension_helpers will check for opnmp support by trying to build
    # some test code with the appropriate flag. openmp provides a big performance boost, but some
    # systems, notably apple's version of clang that xcode provides, don't support it out of the box.

    add_openmp_flags_if_available(extension)
    return [extension]
