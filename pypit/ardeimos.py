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
    detlsize = head0['DETLSIZE']
    x0, x_npix, y0, y_npix = np.array(load_sections(detlsize)).flatten()

    # Create final image
    image = np.zeros((x_npix,y_npix+4*postpix))

    # Setup for datasec, oscansec
    dsec = []
    osec = []


    # get the x and y binning factors...
    binning = head0['BINNING']
    if binning != '1,1':
        msgs.error("This binning for DEIMOS might not work.  But it might..")

    xbin, ybin = [int(ibin) for ibin in binning.split(',')]

    nchip = 8
    ii = 2048
    jj = 4096

    for tt in range(nchip):
        data, oscan = deimos_read_1chip(hdu, tt+1)


        #if n_elements(nobias) eq 0 then nobias = 0
        # Fill the image up

        # y indices
        if tt < 4:
            y1, y2 = 0, jj
        else:
            y1, y2 = jj, 2*jj
        o_y1, o_y2 = y1, y2

        # x
        x1, x2 = (tt%4)*ii, (tt%4 + 1)*ii
        o_x1 = 4*ii + (tt%4)*postpix
        o_x2 = o_x1 + postpix

        # Fill
        image[y1:y2, x1:x2] = data
        image[o_y1:o_y2, o_x1:o_x2] = oscan

        # Sections
        idsec = '[{:d}:{:d},{:d}:{:d}]'.format(y1, y2, x1, x2)
        iosec = '[{:d}:{:d},{:d}:{:d}]'.format(o_y1, o_y2, o_x1, o_x2)
        dsec.append(idsec)
        osec.append(iosec)
    # Return
    return image, head0, (dsec,osec)

def deimos_read_1chip(hdu,chipno):

    # Extract datasec from header
    datsec = hdu[chipno].header['DATASEC']
    detsec = hdu[chipno].header['DETSEC']
    postpix = hdu[0].header['POSTPIX']
    precol = hdu[0].header['PRECOL']

    x1_dat, x2_dat, y1_dat, y2_dat = np.array(load_sections(datsec)).flatten()
    x1_det, x2_det, y1_det, y2_det = np.array(load_sections(detsec)).flatten()

    # This rotates the image to be increasing wavelength to the top
    #data = np.rot90((hdu[chipno].data).T, k=2)
    #nx=data.shape[0]
    #ny=data.shape[1]


    # Science data
    fullimage = hdu[chipno].data
    data = fullimage[x1_dat:x2_dat,y1_dat:y2_dat]

    # Overscan
    oscan = fullimage[:,y2_dat:]

    # Flip as needed
    if x1_det > x2_det:
        data = np.flipud(data)
        oscan = np.flipud(oscan)
    if y1_det > y2_det:
        data = np.fliplr(data)
        oscan = np.fliplr(oscan)


    return data, oscan
