#  Class for organinizing PYPIT setup
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np

#from importlib import reload

from astropy.io import fits
from astropy.table import hstack

from linetools import utils as ltu

from pypit import msgs
from pypit import ardebug as debugger
from pypit import arload
from pypit import arsort
from pypit import ginga

from scipy import ndimage

# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)


class SetupClass(object):
    """Class to handle setup

    Parameters
    ----------

    Attributes
    ----------
    """
    def __init__(self, settings, fitstbl=None):

        # Required parameters
        self.settings = settings

        # Other parameters
        self.fitstbl = fitstbl

        # Attributes

    def build_fitstbl(self, file_list):
        self.fitstbl, updates = arload.load_headers(file_list, self.settings.spect,
                                                    self.settings.argflag)
        self.fitstbl.sort('time')
        return self.fitstbl

    def match_to_science(self):
        self.fitstbl = arsort.new_match_to_science(self.fitstbl,
                                         self.settings.spect,
                                         self.settings.argflag)
        return self.fitstbl

    def sort_data(self):
        filetypes = arsort.new_sort_data(self.fitstbl,
                                        self.settings.spect,
                                        self.settings.argflag)

        # hstack me -- Might over-write self.fitstbl here
        self.fitstbl = hstack([self.fitstbl, filetypes])

        #
        return self.fitstbl

    def write(self, outfile, overwrite=True):
        """
        Write to FITS

        Parameters
        ----------
        outfile : str
        """
        tmp = self.fitstbl.copy()
        #tmp.remove_column('utc')
        tmp.write(outfile, overwrite=overwrite)

    def run(self, armlsd=True, ignore_orders=False, add_user_slits=None):
        """ Main driver for tracing slits.

          Code flow
           1.  Determine approximate slit edges (left, right)
             1b.    Trim down to one pixel per edge per row [seems wasteful, but ok]
           2.  Give edges ID numbers + stitch together partial edges (match_edges)
             2b.   first maxgap option -- NOT recommended
           3.  Assign slits (left, right) ::  Deep algorithm
           4.  For ARMLSD
              -- Trace crude the edges
              -- Do a multi-slit sync to pair up left/right edges
           5.  Remove short slits -- Not recommended for ARMLSD
           6.  Fit left/right slits
           7.  Synchronize
           8.  Extrapolate into blank regions (PCA)

        Parameters
        ----------
        armlsd : bool (optional)
          Running longslit or multi-slit?
        ignore_orders : bool (optional)
          Perform ignore_orders algorithm (recommended only for echelle data)
        add_user_slits : list of lists
          List of 2 element lists, each an [xleft, xright] pair specifying a slit edge
          These are specified at mstrace.shape[0]//2

        Returns
        -------
        lcen : ndarray
          Left edge traces
        rcen  : ndarray
          Right edge traces
        extrapord
        """
        # Specify a single slit?
        if len(self.settings['trace']['slits']['single']) > 0:  # Single slit
            self._edgearr_single_slit()
            self.user_set = True
        else:  # Generate the edgearr from the input trace image
            self._edgearr_from_binarr()
            self.user_set = False


    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: >'.format(self.__class__.__name__)
        return txt



