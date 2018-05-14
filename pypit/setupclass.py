#  Class for organinizing PYPIT setup
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np

#from importlib import reload

from astropy.io import fits

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

    def sort_data(self):
        filetypes = arsort.new_sort_data(self.fitstbl,
                                        self.settings.spect,
                                        self.settings.argflag)
        return filetypes

    def show(self, attr, display='ginga'):
        """
        Display an image or spectrum in TraceSlits

        Parameters
        ----------
        attr : str
          'edges' -- Show the mstrace image and the edges
          'edgearr' -- Show the edgearr image
          'siglev' -- Show the Sobolev image
        display : str (optional)
          'ginga' -- Display to an RC Ginga
        """
        if attr == 'edges':
            viewer, ch = ginga.show_image(self.mstrace)
            if self.lcen is not None:
                ginga.show_slits(viewer, ch, self.lcen, self.rcen, np.arange(self.lcen.shape[1]) + 1, pstep=50)
        elif attr == 'edgearr':
            # TODO -- Figure out how to set the cut levels
            debugger.show_image(self.edgearr)
        elif attr == 'siglev':
            # TODO -- Figure out how to set the cut levels
            debugger.show_image(self.siglev)

    def write(self, root):
        """
        Write the main pieces of TraceSlits to the hard drive
          FITS -- mstrace and other images
          JSON -- steps, settings, ts_dict

        Parameters
        ----------
        root : str
          Path+root name for the output files
        """

        # Images
        outfile = root+'.fits'
        hdu = fits.PrimaryHDU(self.mstrace)
        hdu.name = 'MSTRACE'
        hdulist = [hdu]
        if self.edgearr is not None:
            hdue = fits.ImageHDU(self.edgearr)
            hdue.name = 'EDGEARR'
            hdulist.append(hdue)
        if self.siglev is not None:
            hdus = fits.ImageHDU(self.siglev)
            hdus.name = 'SIGLEV'
            hdulist.append(hdus)
        hdup = fits.ImageHDU(self.pixlocn)
        hdup.name = 'PIXLOCN'
        hdulist.append(hdup)
        if self.input_binbpx:  # User inputted
            hdub = fits.ImageHDU(self.binbpx)
            hdub.name = 'BINBPX'
            hdulist.append(hdub)
        if self.lcen is not None:
            hdulf = fits.ImageHDU(self.lcen)
            hdulf.name = 'LCEN'
            hdulist.append(hdulf)
            hdurt = fits.ImageHDU(self.rcen)
            hdurt.name = 'RCEN'
            hdulist.append(hdurt)

        # Write
        hdul = fits.HDUList(hdulist)
        hdul.writeto(outfile, overwrite=True)
        msgs.info("Writing TraceSlit arrays to {:s}".format(outfile))

        # dict of steps, settings and more
        out_dict = {}
        out_dict['settings'] = self.settings
        if self.tc_dict is not None:
            out_dict['tc_dict'] = self.tc_dict
        out_dict['steps'] = self.steps
        # Clean+Write
        outfile = root+'.json'
        clean_dict = ltu.jsonify(out_dict)
        ltu.savejson(outfile, clean_dict, overwrite=True, easy_to_read=True)
        msgs.info("Writing TraceSlit dict to {:s}".format(outfile))

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



