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


class Paths():
    """Paths List of hardwired paths within the pypeit.data module

    [extended_summary]
    """

    # Class Attributes -- Hardwired Paths
    data = resource_filename('pypeit', 'data')

    # Telluric Corrections
    telgrid = os.path.join(data, 'telluric', 'atm_grids')
    tel_model = os.path.join(data, 'telluric', 'models')

    # Wavelength Calibrations
    arclines = os.path.join(data, 'arc_lines')
    reid_arxiv = os.path.join(data, 'arc_lines', 'reid_arxiv')
    linelist = os.path.join(data, 'arc_lines', 'lists')
    nist = os.path.join(data, 'arc_lines', 'NIST')
    arc_plot = os.path.join(data, 'arc_lines', 'plots')

    # Flux Calibrations
    standards = os.path.join(data, 'standards')
    extinction = os.path.join(data, 'extinction')
    skisim = os.path.join(data, 'skisim')
    filters = os.path.join(data, 'filters')

    # Other
    sky_spec = os.path.join(data, 'sky_spec')
    static_calibs = os.path.join(data, 'static_calibs')

    def __init__(self):
        pass



# Loading Functions for Particular File Types ================================#
def load_telluric_grid(filename):

    # Check for existance of file parameter
    if not filename:
        msgs.error("No file specified for telluric correction.  "
                   "See https://pypeit.readthedocs.io/en/latest/telluric.html")

    # Add standard data path to the filename, as contained in default pypeit pars
    file_with_path = os.path.join(Paths.telgrid, filename)

    # Check for existance of file
    if not os.path.isfile(file_with_path):
        msgs.error(f"File {file_with_path} is not on your disk.  "
                   "You likely need to download the Telluric files.  "
                   "See https://pypeit.readthedocs.io/en/release/installing.html#atmospheric-model-grids")

    return io.fits_open(file_with_path)


def load_thar_spec():
    return io.fits_open(os.path.join(Paths.arclines, 'thar_spec_MM201006.fits'))
