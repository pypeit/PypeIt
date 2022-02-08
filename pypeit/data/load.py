# -*- coding: utf-8 -*-
"""
Data locations for built-in PypeIt data files

.. include:: ../include/links.rst
"""
import os
from pkg_resources import resource_filename

import numpy as np

from astropy.table import Table

from pypeit import io
from pypeit import msgs


# Hardwired Paths
# TODO: There is probably a better way to package these...
TELGRID_PATH = os.path.join(resource_filename('pypeit', 'data'), 'telluric', 'atm_grids')
ARCLINES_DIR = os.path.join(resource_filename('pypeit', 'data'), 'arc_lines')
REID_ARXIV_PATH = os.path.join(ARCLINES_DIR, 'reid_arxiv')
LINE_PATH = os.path.join(ARCLINES_DIR, 'lists')
NIST_PATH = os.path.join(ARCLINES_DIR, 'NIST')
ARC_PLOT_PATH = os.path.join(ARCLINES_DIR, 'plots')


def load_telluric_grid(filename):

    # Check for existance of file parameter
    if not filename:
        msgs.error("No file specified for telluric correction.  "
                   "See https://pypeit.readthedocs.io/en/latest/telluric.html")

    # Add standard data path to the filename, as contained in default pypeit pars
    file_with_path = os.path.join(TELGRID_PATH, filename)

    # Check for existance of file
    if not os.path.isfile(file_with_path):
        msgs.error(f"File {file_with_path} is not on your disk.  "
                   "You likely need to download the Telluric files.  "
                   "See https://pypeit.readthedocs.io/en/release/installing.html#atmospheric-model-grids")

    return io.fits_open(file_with_path)


def load_thar_spec():
    return io.fits_open(ARCLINES_DIR+'thar_spec_MM201006.fits')
