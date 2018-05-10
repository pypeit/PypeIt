# Module for guiding Slit/Order tracing
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np

from importlib import reload

from astropy.io import fits

from linetools import utils as ltu

from pypit import msgs
from pypit import ardebug as debugger
from pypit import arload
from pypit import arparse
from pypit import arutils
from pypit import ginga

from scipy import ndimage

# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

# Place these here or elsewhere?
#  Wherever they be, they need to be defined, described, etc.
default_settings = dict(run={'spectrograph': 'UNKNOWN'},
                        detector={'numamplifiers': 1,
                                  'dispaxis': 0,  # Spectra aligned with columns
                                  'dataext': None})


class ProcessImages(object):
    """Base class to guide image loading+processing

    Parameters
    ----------
    file_list : list
      List of filenames
    det : int, optional
      Detector index, starts at 1
    settings : dict, optional
      Settings for trace slits
    user_settings : dict, optional
      Allow for user to over-ride individual internal/default settings
      without providing a full settings dict

    Attributes
    ----------
    images : list
    stack : ndarray
    """
    def __init__(self, file_list, det=1, settings=None, user_settings=None):

        # Parameters
        self.file_list = file_list
        self.det = det
        self.dnum = arparse.get_dnum(self.det)

        # Settings
        if settings is None:
            self.settings = default_settings.copy()
        else:
            # The copy allows up to update settings with user settings without changing the original
            self.settings = settings.copy()
        if user_settings is not None:
            # This only works to replace entire dicts
            #    Hopefully parsets will work more cleanly..
            self.settings.update(user_settings)

        # Main (possible) outputs
        self.stack = None
        self.steps = []

        # Key Internals
        self.images = None
        self.headers = None

    @property
    def nfiles(self):
        return len(self.file_list)

    @property
    def nloaded(self):
        if self.images is None:
            return 0
        else:
            return len(self.images)

    def load(self):
        self.images = []  # Zeros out any previous load
        self.headers = []  # Zeros out any previous load
        for ifile in self.file_list:
            img, head = arload.load_raw_frame(self.settings['run']['spectrograph'], ifile, self.det,
                                        dataext=self.settings['detector']['dataext'],
                                        disp_dir=self.settings['detector']['dispaxis'])
            # Save
            self.images.append(img)
            self.headers.append(head)
        # Step
        self.steps.append(inspect.stack()[0][3])

    def show(self, attr, display='ginga'):
        if attr == 'edges':
            viewer, ch = ginga.show_image(self.mstrace)
            ginga.show_slits(viewer, ch, self.lcen, self.rcen, np.arange(self.lcen.shape[1]) + 1, pstep=50)
        elif attr == 'edgearr':
            # TODO -- Figure out how to set the cut levels
            debugger.show_image(self.edgearr)
        elif attr == 'siglev':
            # TODO -- Figure out how to set the cut levels
            debugger.show_image(self.siglev)

    def write(self, root):

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
            hdus = fits.ImageHDU(self.edgearr)
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

        # dict
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

        # Assign a number to each edge 'grouping'
        self._match_edges()

        # Add in a single left/right edge?
        self._add_left_right()

        # If slits are set as "close" by the user, take the absolute value
        # of the detections and ignore the left/right edge detections
        #  Use of maxgap is NOT RECOMMENDED
        if self.settings['trace']['slits']['maxgap'] is not None:
            self._maxgap_prep()

        # Assign edges
        self._assign_edges()

        # Handle close edges (as desired by the user)
        #  JXP does not recommend using this method for multislit
        if self.settings['trace']['slits']['maxgap'] is not None:
            self._maxgap_close()

        # Final left/right edgearr fussing (as needed)
        if not self.user_set:
            self._final_left_right()

        # Trace crude me
        #   -- Mainly to deal with duplicates and improve the traces
        #   -- Developed for ARMLSD not ARMED
        if armlsd:
            self._mslit_tcrude()

        # Synchronize and add in edges
        if armlsd:
            self._mslit_sync()

        # Add user input slits
        if add_user_slits is not None:
            self.add_user_slits(add_user_slits)

        # Ignore orders/slits on the edge of the detector when they run off
        #    Recommended for Echelle only
        if ignore_orders:
            self._ignore_orders()

        # Fit edges
        self.set_lrminx()
        self._fit_edges('left')
        self._fit_edges('right')

        # Are we done, e.g. longslit?
        #   Check if no further work is needed (i.e. there only exists one order)
        if self.chk_for_longslit():
            return self.lcen, self.rcen, np.zeros(1, dtype=np.bool)

        # Synchronize
        #   For multi-silt, mslit_sync will have done most of the work already..
        self._synchronize()

        # PCA?
        #  Whether or not a PCA is performed, lcen and rcen are generated for the first time
        self._pca()

        # Remove any slits that are completely off the detector
        #   Also remove short slits here for multi-slit and long-slit (aligntment stars)
        self.trim_slits(usefracpix=armlsd)

        # Illustrate where the orders fall on the detector (physical units)
        if msgs._debug['trace']:
            self.show('edges')
            debugger.set_trace()

        # Finish
        return self.lcen, self.rcen, self.extrapord

    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: >'.format(self.__class__.__name__)
        return txt


