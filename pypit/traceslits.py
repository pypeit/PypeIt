# Module for guiding Slit/Order tracing
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np

from importlib import reload

from astropy.io import fits

from linetools import utils as ltu

from pypit import msgs
from pypit import ardebug as debugger
from pypit import artrace
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
default_settings = dict(trace={'slits': {'single': [],
                               'function': 'legendre',
                               'polyorder': 3,
                               'diffpolyorder': 2,
                               'fracignore': 0.01,
                               'medrep': 0,
                               'number': -1,
                               'maxgap': None,
                               'sigdetect': 20.,
                               'pca': {'params': [3,2,1,0,0,0], 'type': 'pixel',
                                       'extrapolate': {'pos': 0, 'neg':0}},
                               'sobel': {'mode': 'nearest'}}})

class TraceSlits(object):
    """Class to guide slit/order tracing

    Parameters:
    ----------

    Attributes:
    ----------
    """
    @classmethod
    def from_files(cls, root):
        # FITS
        fits_file = root+'.fits'
        hdul = fits.open(fits_file)
        names = [ihdul.name for ihdul in hdul]

        mstrace = hdul[names.index('MSTRACE')].data
        pixlocn = hdul[names.index('PIXLOCN')].data
        if 'BINBPX' in names:
            binbpx = hdul[names.index('BINBPX')].data
            msgs.info("Loading BPM from {:s}".format(fits_file))
        else:
            binbpx = None

        # JSON
        json_file = root+'.json'
        ts_dict = ltu.loadjson(json_file)
        slf = cls(mstrace, pixlocn, binbpx=binbpx, settings=ts_dict['settings'])

        # Fill in a bit more
        slf.steps = ts_dict['steps']
        if 'LCEN' in names:
            slf.lcen = hdul[names.index('LCEN')].data
            slf.rcen = hdul[names.index('RCEN')].data
            msgs.info("Loading LCEN, RCEN from {:s}".format(fits_file))
        if 'EDGEARR' in names:
            slf.edgearr = hdul[names.index('EDGEARR')].data
            msgs.info("Loading EDGEARR from {:s}".format(fits_file))
        if 'SIGLEV' in names:
            slf.siglev = hdul[names.index('SIGLEV')].data
            msgs.info("Loading SIGLEV from {:s}".format(fits_file))
        slf.tc_dict = ts_dict['tc_dict']

        # Return
        return slf

    def __init__(self, mstrace, pixlocn, binbpx=None, settings=None, det=None, ednum=100000):
        # TODO -- Remove pixlocn as a required item

        # Required attributes
        self.mstrace = mstrace
        if settings is None:
            self.settings = default_settings.copy()
        else:
            self.settings = settings
        if binbpx is None: # Bad pixel array
            self.binbpx = np.zeros_like(mstrace)
            self.input_binbpx = False # For writing
        else:
            self.binbpx = binbpx
            self.input_binbpx = True
        self.pixlocn = pixlocn

        # Optional attributes
        self.det = det

        # Generate binarr
        self.make_binarr()

        # Main outputs
        self.tc_dict = None
        self.edgearr = None
        self.steps = []
        self.lcen = None
        self.rcen = None
        self.extrapord = None

        # Key Internals
        self.ednum = ednum
        self.siglev = None
        self.user_set = None
        self.lcnt = None
        self.rcnt = None
        self.lmin = None
        self.lmax = None
        self.rmin = None
        self.rmax = None

    def make_binarr(self):
        # Generate first edgearr from mstrace or user-supplied
        #  Only filter in the spectral dimension, not spatial!
        self.binarr = ndimage.uniform_filter(self.mstrace, size=(3, 1), mode='mirror')

    def _edgearr_from_binarr(self):
        self.siglev, self.edgearr = artrace.edgearr_from_binarr(self.binarr, self.binbpx,
                                                                medrep=self.settings['trace']['slits']['medrep'],
                                                                sobel_mode=self.settings['trace']['slits']['sobel']['mode'],
                                                                sigdetect=self.settings['trace']['slits']['sigdetect'],
                                                                number_slits=self.settings['trace']['slits']['number'])
        self.user_set = False
        # Step
        self.steps.append(inspect.stack()[0][3])

    def _edgearr_single_slit(self):
        #  Note this is different from add_user_slits (which is handled below)
        #  This trace slits single option is likely to be deprecated
        iledge, iredge = (self.det-1)*2, (self.det-1)*2+1
        ledge = self.settings['trace']['slits']['single'][iledge]
        redge = self.settings['trace']['slits']['single'][iredge]
        self.edgearr = artrace.edgearr_from_user(self.mstrace.shape, ledge, redge, self.det)
        self.user_set = True
        self.siglev = None
        # Step
        self.steps.append(inspect.stack()[0][3])

    def _add_left_right(self):
        self.edgearr, self.lcnt, self.rcnt = artrace.edgearr_add_left_right(
            self.edgearr, self.binarr, self.binbpx, self.lcnt, self.rcnt, self.ednum)
        # Step
        self.steps.append(inspect.stack()[0][3])

    def _add_user_slits(self, add_user_slits, run_to_finish=False):
        # Reset (if needed) -- For running after PYPIT took a first pass
        self.reset_edgearr_ednum()
        # Add user input slits
        self.edgearr = artrace.add_user_edges(self.edgearr, self.siglev, self.tc_dict, add_user_slits)
        # Finish
        if run_to_finish:
            self.set_lrminx()
            self._fit_edges('left')
            self._fit_edges('right')
            self._synchronize()
            self._pca()
            self.trim_slits()
        # Step
        self.steps.append(inspect.stack()[0][3])

    def _assign_edges(self):

        # Assign left edges
        msgs.info("Assigning left slit edges")
        if self.lcnt == 1:
            self.edgearr[np.where(self.edgearr <= -2*self.ednum)] = -self.ednum
        else:
            artrace.assign_slits(self.binarr, self.edgearr, lor=-1, settings=self.settings)
        # Assign right edges
        msgs.info("Assigning right slit edges")
        if self.rcnt == 1:
            self.edgearr[np.where(self.edgearr >= 2*self.ednum)] = self.ednum
        else:
            artrace.assign_slits(self.binarr, self.edgearr, lor=+1, settings=self.settings)
        # Steps
        self.steps.append(inspect.stack()[0][3])

    def chk_for_longslit(self):
        # Are we done, e.g. longslit?
        #   Check if no further work is needed (i.e. there only exists one order)
        if (self.lmax+1-self.lmin == 1) and (self.rmax+1-self.rmin == 1):
            plxbin = self.pixlocn[:, :, 0].copy()
            minvf, maxvf = plxbin[0, 0], plxbin[-1, 0]
            # Just a single order has been identified (i.e. probably longslit)
            msgs.info("Only one slit was identified. Should be a longslit.")
            xint = self.pixlocn[:, 0, 0]
            # Finish
            self.lcen = np.zeros((self.mstrace.shape[0], 1))
            self.rcen = np.zeros((self.mstrace.shape[0], 1))
            self.lcen[:, 0] = arutils.func_val(self.lcoeff[:, 0], xint,
                                                  self.settings['trace']['slits']['function'],
                                             minv=minvf, maxv=maxvf)
            self.rcen[:, 0] = arutils.func_val(self.rcoeff[:, 0], xint,
                                                  self.settings['trace']['slits']['function'],
                                             minv=minvf, maxv=maxvf)
            return True
        else:
            return False

    def _final_left_right(self):
        # Final left/right edgearr fussing (as needed)
        self.edgearr, self.lcnt, self.rcnt = artrace.edgearr_final_left_right(
            self.edgearr, self.ednum, self.siglev)
        # Steps
        self.steps.append(inspect.stack()[0][3])

    def _fit_edges(self, side):
        # Setup for fitting
        plxbin = self.pixlocn[:, :, 0].copy()
        plybin = self.pixlocn[:, :, 1].copy()

        # Fit
        if side == 'left':
            self.lcoeff, self.lnmbrarr, self.ldiffarr, self.lwghtarr = artrace.fit_edges(
                self.edgearr, self.lmin, self.lmax, plxbin, plybin,
                left=True, polyorder=self.settings['trace']['slits']['polyorder'],
                function=self.settings['trace']['slits']['function'])
        else:
            self.rcoeff, self.rnmbrarr, self.rdiffarr, self.rwghtarr = artrace.fit_edges(
                self.edgearr, self.rmin, self.rmax, plxbin, plybin,
                left=False, polyorder=self.settings['trace']['slits']['polyorder'],
                function=self.settings['trace']['slits']['function'])

        # Steps
        self.steps.append(inspect.stack()[0][3]+'_{:s}'.format(side))

    def _ignore_orders(self):
        # Ignore orders/slits on the edge of the detector when they run off
        #    Recommended for Echelle only
        self.edgearr, self.lmin, self.lmax, self.rmin, self.rmax = artrace.edgearr_ignore_orders(
            self.edgearr, self.settings['trace']['slits']['fracignore'])
        # Steps
        self.steps.append(inspect.stack()[0][3])

    def _match_edges(self):
        # Assign a number to each edge 'grouping'
        __edgearr = self.edgearr.copy()
        self.lcnt, self.rcnt = artrace.new_match_edges(__edgearr, self.ednum)
        self.edgearr = __edgearr
        # Sanity check (unlikely we will ever hit this)
        if self.lcnt >= self.ednum or self.rcnt >= self.ednum:
            msgs.error("Found more edges than allowed by ednum. Set ednum to a larger number.")
        # Step
        self.steps.append(inspect.stack()[0][3])

    def _maxgap_prep(self):
        self.edgearrcp = self.edgearr.copy()
        self.edgearr[np.where(self.edgearr < 0)] += 1 + np.max(self.edgearr) - np.min(self.edgearr)
        # Step
        self.steps.append(inspect.stack()[0][3])

    def _maxgap_close(self):
        # Handle close edges (as desired by the user)
        #  JXP does not recommend using this method for multislit
        self.edgearr = artrace.edgearr_close_slits(self.binarr, self.edgearr,
                                              self.edgearrcp, self.ednum, self.settings)
        # Step
        self.steps.append(inspect.stack()[0][3])

    def _mslit_sync(self, debug=False):
        if debug:
            reload(artrace)
        #
        self.edgearr = artrace.mslit_sync(self.edgearr, self.tc_dict, self.ednum)
        # Step
        self.steps.append(inspect.stack()[0][3])

    def _mslit_tcrude(self):
        # Trace crude me
        self.edgearr, self.tc_dict = artrace.tcrude_edgearr(self.edgearr, self.siglev, self.ednum)
        # Step
        self.steps.append(inspect.stack()[0][3])

    def _pca(self):
        if self.settings['trace']['slits']['pca']['type'] == 'order':
            self._pca_order_slit_edges()
        elif self.settings['trace']['slits']['pca']['type'] == 'pixel':
            self._pca_pixel_slit_edges()
        else: # No PCA
            self.set_lrcen()

    def _pca_order_slit_edges(self):
        plxbin = self.pixlocn[:, :, 0].copy()
        self.lcen, self.rcen, self.extrapord = artrace.pca_order_slit_edges(self.binarr, self.edgearr,
                                                                    self.lcent, self.rcent, self.gord,
                                                                    self.lcoeff, self.rcoeff, plxbin,
                                                                    self.slitcen, self.pixlocn, self.settings)
        # Step
        self.steps.append(inspect.stack()[0][3])

    def _pca_pixel_slit_edges(self):
        plxbin = self.pixlocn[:, :, 0].copy()
        self.lcen, self.rcen, self.extrapord = artrace.pca_pixel_slit_edges(self.binarr,
                                                                            self.edgearr, self.lcoeff, self.rcoeff,
                                                                            self.ldiffarr, self.rdiffarr, self.lnmbrarr,
                                                                            self.rnmbrarr, self.lwghtarr, self.rwghtarr, self.lcent,
                                                                            self.rcent, plxbin, self.settings)
        # Step
        self.steps.append(inspect.stack()[0][3])

    def reset_edgearr_ednum(self):
        if np.max(self.edgearr) < self.ednum:
            neg = np.where(self.edgearr < 0)
            self.edgearr[neg] -= (self.ednum - 1)
            pos = np.where(self.edgearr > 0)
            self.edgearr[pos] += (self.ednum - 1)

    def set_lrminx(self):
        ww = np.where(self.edgearr < 0)
        self.lmin, self.lmax = -np.max(self.edgearr[ww]), -np.min(self.edgearr[ww])  # min/max are switched because of the negative signs
        ww = np.where(self.edgearr > 0)
        self.rmin, self.rmax = np.min(self.edgearr[ww]), np.max(self.edgearr[ww])  # min/max are switched because of the negative signs

    def set_lrcen(self):
        allord = np.arange(self.lcent.shape[0])
        maskord = np.where((np.all(self.lcent, axis=1) == False) | (np.all(self.rcent, axis=1) == False))[0]
        ww = np.where(np.in1d(allord, maskord) == False)[0]
        self.lcen = self.lcent[ww, :].T.copy()
        self.rcen = self.rcent[ww, :].T.copy()
        self.extrapord = np.zeros(self.lcen.shape[1], dtype=np.bool)

    def _synchronize(self):
        plxbin = self.pixlocn[:, :, 0].copy()
        msgs.info("Synchronizing left and right slit traces")
        self.lcent, self.rcent, self.gord, self.lcoeff, self.ldiffarr, self.lnmbrarr, self.lwghtarr, self.rcoeff, self.rdiffarr, self.rnmbrarr, self.rwghtarr = artrace.synchronize_edges(
            self.binarr, self.edgearr, plxbin, self.lmin, self.lmax, self.lcoeff, self.rmin, self.rcoeff,
            self.lnmbrarr, self.ldiffarr, self.lwghtarr, self.rnmbrarr, self.rdiffarr, self.rwghtarr, self.settings)
        self.slitcen = 0.5*(self.lcent+self.rcent).T
        # Step
        self.steps.append(inspect.stack()[0][3])

    def trim_slits(self, usefracpix=True):
        nslit = self.lcen.shape[1]
        mask = np.zeros(nslit)
        fracpix = int(self.settings['trace']['slits']['fracignore']*self.mstrace.shape[1])
        for o in range(nslit):
            if np.min(self.lcen[:, o]) > self.mstrace.shape[1]:
                mask[o] = 1
                msgs.info("Slit {0:d} is off the detector - ignoring this slit".format(o+1))
            elif np.max(self.rcen[:, o]) < 0:
                mask[o] = 1
                msgs.info("Slit {0:d} is off the detector - ignoring this slit".format(o + 1))
            if usefracpix:
                if np.median(self.rcen[:,o]-self.lcen[:,o]) < fracpix:
                    mask[o] = 1
                    msgs.info("Slit {0:d} is less than fracignore - ignoring this slit".format(o + 1))
        # Trim
        wok = np.where(mask == 0)[0]
        self.lcen = self.lcen[:, wok]
        self.rcen = self.rcen[:, wok]
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


    def __repr__(self):
        # Generate sets string
        return 'blah'
        #return ('<SpecObjExp: {:s} == Setup {:s} Object at {:g} in Slit at {:g} with det={:s}, scidx={:d} and objtype={:s}>'.format(
        #        self.idx, self.config, self.xobj, self.slitcen, sdet, self.scidx, self.objtype))


def run(mstrace, pixlocn, det=None, settings=None,
                       binbpx=None, ignore_orders=False,
                       add_user_slits=None, armlsd=True):
    """ Main driver for tracing slits.
    May be renamed to trace_slits

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
    mstrace : ndarray
    pixlocn : ndarray
    det : int (optional)
    settings : dict (optional)
    ednum : int (optional)
       A large dummy number used for slit edge assignment.
       ednum should be larger than the number of edges detected
    binbpx : ndarray
      Bad pixel mask
    det : int (optional)
      Required for single slilt specification
    ignore_orders : bool (optional)
      Perform ignore_orders algorithm (recommended only for echelle data)
    add_user_slits : list of lists
      List of 2 element lists, each an [xleft, xright] pair specifying a slit edge
      These are specified at mstrace.shape[0]//2
    armlsd : bool (optional)
      Running with ARMLSD ?

    Returns
    -------
    lcen : ndarray
      Left edge traces
    rcen  : ndarray
      Right edge traces
    extrapord
    """
    # Init TraceSlits
    tslits = TraceSlits(mstrace, pixlocn, binbpx=binbpx, settings=settings, det=det)

    # Specify a single slit?
    if len(settings['trace']['slits']['single']) > 0:  # Single slit
        tslits._edgearr_single_slit()
    else:  # Generate the edgearr from the input trace image
        tslits._edgearr_from_binarr()

    # Assign a number to each edge 'grouping'
    tslits._match_edges()

    # Add in a single left/right edge?
    tslits._add_left_right()

    # If slits are set as "close" by the user, take the absolute value
    # of the detections and ignore the left/right edge detections
    #  Use of maxgap is NOT RECOMMENDED
    if settings['trace']['slits']['maxgap'] is not None:
        tslits._maxgap_prep()

    # Assign edges
    tslits._assign_edges()

    # Handle close edges (as desired by the user)
    #  JXP does not recommend using this method for multislit
    if settings['trace']['slits']['maxgap'] is not None:
        tslits._maxgap_close()

    # Final left/right edgearr fussing (as needed)
    if not tslits.user_set:
        tslits._final_left_right()

    # Trace crude me
    #   -- Mainly to deal with duplicates and improve the traces
    #   -- Developed for ARMLSD not ARMED
    if armlsd:
        tslits._mslit_tcrude()

    # Synchronize and add in edges
    if armlsd:
        tslits._mslit_sync()

    # Add user input slits
    if add_user_slits is not None:
        tslits._add_user_slits(add_user_slits)

    # Ignore orders/slits on the edge of the detector when they run off
    #    Recommended for Echelle only
    if ignore_orders:
        tslits._ignore_orders()

    # Fit edges
    tslits.set_lrminx()
    tslits._fit_edges('left')
    tslits._fit_edges('right')

    # Are we done, e.g. longslit?
    #   Check if no further work is needed (i.e. there only exists one order)
    if tslits.chk_for_longslit():
        return tslits.lcen, tslits.rcen, np.zeros(1, dtype=np.bool), tslits

    # Synchronize
    #   For multi-silt, mslit_sync will have done most of the work already..
    debugger.set_trace()
    tslits._synchronize()

    # PCA?
    tslits._pca()

    # Remove any slits that are completely off the detector
    #   Also remove short slits here for multi-slit and long-slit (aligntment stars)
    tslits.trim_slits(usefracpix=armlsd)

    # Illustrate where the orders fall on the detector (physical units)
    '''
    if msgs._debug['trace']:
        viewer, ch = ginga.show_image(mstrace)
        ginga.show_slits(viewer, ch, lcen, rcen, np.arange(nslit) + 1)
        debugger.set_trace()
    '''

    # Finish
    return tslits.lcen, tslits.rcen, tslits.extrapord, tslits
