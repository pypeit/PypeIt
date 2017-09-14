'''
Implements DEIMOS-specific functions, including reading in slitmask design files.
'''
from __future__ import absolute_import, division, print_function

import glob
import numpy as np
# from astropy.io import fits
import astropy.io.fits as pyfits

from pypit import armsgs
from pypit.arparse import load_sections
from pypit import ardebug as debugger
#from IPython import embed

try:
    basestring
except NameError:  # For Python 3
    basestring = str


# Logging
msgs = armsgs.get_logger()


def read_deimos(raw_file):
    """
    Read a raw DEIMOS data frame (one or more detectors)
    Packed in a multi-extension HDU
    Based on pypit.arlris.read_lris...
       Based on readmhdufits.pro

    Parameters
    ----------
    raw_file : str
      Filename
    det : int
      detector index starting at 1
    trim : bool, optional
      Trim the image?

    Returns
    -------
    array : ndarray
      Combined image
    header : FITS header
    sections : list
      List of datasec, oscansec, ampsec sections
    """

    # Check for file; allow for extra .gz, etc. suffix
    fil = glob.glob(raw_file + '*')
    if len(fil) != 1:
        msgs.error("Found {:d} files matching {:s}".format(len(fil)))
    # Read
    try:
        msgs.info("Reading DEIMOS file: {:s}".format(fil[0]))
    except AttributeError:
        print("Reading DEIMOS file: {:s}".format(fil[0]))

    hdu = pyfits.open(fil[0])
    head0 = hdu[0].header

    # Get post, pre-pix values
    precol = head0['PRECOL']
    postpix = head0['POSTPIX']
    preline = head0['PRELINE']
    postline = head0['POSTLINE']

    # Setup for datasec, oscansec
    dsec = []
    osec = []

    # get the x and y binning factors...
    binning = head0['BINNING']
    xbin, ybin = [int(ibin) for ibin in binning.split(',')]

    datsec = hdu[7].header['DATASEC']
    detsec = hdu[7].header['DETSEC']
    x1_dat, x2_dat, y1_dat, y2_dat = np.array(load_sections(datsec)).flatten()
    x1_det, x2_det, y1_det, y2_det = np.array(load_sections(detsec)).flatten()

    # This rotates the image to be increasing wavelength to the top
    data = np.rot90((hdu[7].data).T,k=2)
    nx=data.shape[0]
    ny=data.shape[1]

    dsec = '[{:d}:{:d},{:d}:{:d}]'.format(postpix+1, nx-precol, y1_dat, y2_dat)  # Eliminate lines
    osec = '[{:d}:{:d},{:d}:{:d}]'.format(1, postpix, y1_dat, y2_dat)  # Eliminate lines
    return data, head0, (dsec, osec)


