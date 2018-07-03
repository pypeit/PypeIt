'''
Implements DEIMOS-specific functions, including reading in slitmask design files.
'''
from __future__ import absolute_import, division, print_function

import glob
import numpy as np
from astropy.io import fits

from pypit import msgs
from pypit import arparse
from pypit.spectrographs import spectroclass

from pypit import ardebug as debugger


# Logging
#msgs = armsgs.get_logger()

class DEIMOSSpectrograph(spectroclass.Spectrograph):

    def __init__(self):

        spectroclass.Spectrograph.__init__(self)
        self.spectrograph = 'keck_lris'  # Note this is a base of keck_lris_red and keck_lris_blue; we might have to sub-class each

    def load_raw_img_head(self, raw_file, det=1, **null_kwargs):
        raw_img, head0, _ = read_deimos(raw_file, det=det)

        return raw_img, head0

    def get_datasec(self, filename, det, settings_det):

        datasec, oscansec, naxis0, naxis1 = [], [], 0, 0
        temp, head0, secs = read_deimos(filename, det)
        for kk in range(settings_det['numamplifiers']):
            datasec.append(arparse.load_sections(secs[0][kk], fmt_iraf=False))
            oscansec.append(arparse.load_sections(secs[1][kk], fmt_iraf=False))

        # Need naxis0, naxis1 too
        naxis0 = temp.shape[0]
        naxis1 = temp.shape[1]

        # Return
        return datasec, oscansec, naxis0, naxis1


def read_deimos(raw_file, det=None):
    """
    Read a raw DEIMOS data frame (one or more detectors)
    Packed in a multi-extension HDU
    Based on pypit.arlris.read_lris...
       Based on readmhdufits.pro

    Parameters
    ----------
    raw_file : str
      Filename

    Returns
    -------
    array : ndarray
      Combined image
    header : FITS header
    sections : tuple
      List of datasec, oscansec sections
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

    hdu = fits.open(fil[0])
    head0 = hdu[0].header

    # Get post, pre-pix values
    precol = head0['PRECOL']
    postpix = head0['POSTPIX']
    preline = head0['PRELINE']
    postline = head0['POSTLINE']
    detlsize = head0['DETLSIZE']
    x0, x_npix, y0, y_npix = np.array(arparse.load_sections(detlsize)).flatten()

    # Create final image
    if det is None:
        image = np.zeros((x_npix,y_npix+4*postpix))

    # Setup for datasec, oscansec
    dsec = []
    osec = []

    # get the x and y binning factors...
    binning = head0['BINNING']
    if binning != '1,1':
        msgs.error("This binning for DEIMOS might not work.  But it might..")

    xbin, ybin = [int(ibin) for ibin in binning.split(',')]

    # DEIMOS detectors
    nchip = 8


    if det is None:
        chips = range(nchip)
    else:
        chips = [det-1] # Indexing starts at 0 here
    # Loop
    for tt in chips:
        data, oscan = deimos_read_1chip(hdu, tt+1)


        #if n_elements(nobias) eq 0 then nobias = 0


        # One detector??
        if det is not None:
            image = np.zeros((data.shape[0],data.shape[1]+oscan.shape[1]))

        # Indexing
        x1, x2, y1, y2, o_x1, o_x2, o_y1, o_y2 = indexing(tt, postpix, det=det)

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


def indexing(itt, postpix, det=None):
    """
    Some annoying book-keeping for instrument placement.

    Parameters
    ----------
    itt : int
    postpix : int
    det : int, optional

    Returns
    -------

    """
    # Deal with single chip
    if det is not None:
        tt = 0
    else:
        tt = itt
    ii = 2048
    jj = 4096
    # y indices
    if tt < 4:
        y1, y2 = 0, jj
    else:
        y1, y2 = jj, 2*jj
    o_y1, o_y2 = y1, y2

    # x
    x1, x2 = (tt%4)*ii, (tt%4 + 1)*ii
    if det is None:
        o_x1 = 4*ii + (tt%4)*postpix
    else:
        o_x1 = ii + (tt%4)*postpix
    o_x2 = o_x1 + postpix

    # Return
    return x1, x2, y1, y2, o_x1, o_x2, o_y1, o_y2

def deimos_read_1chip(hdu,chipno):
    """ Read one of the DEIMOS detectors

    Parameters
    ----------
    hdu : HDUList
    chipno : int

    Returns
    -------
    data : ndarray
    oscan : ndarray
    """

    # Extract datasec from header
    datsec = hdu[chipno].header['DATASEC']
    detsec = hdu[chipno].header['DETSEC']
    postpix = hdu[0].header['POSTPIX']
    precol = hdu[0].header['PRECOL']

    x1_dat, x2_dat, y1_dat, y2_dat = np.array(arparse.load_sections(datsec)).flatten()
    x1_det, x2_det, y1_det, y2_det = np.array(arparse.load_sections(detsec)).flatten()

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

    # Return
    return data, oscan

def bpm(det):
    """ Generate a BPM for DEIMOS
    Currently assumes 1x1 binning

    Parameters
    ----------
    det : int

    Returns
    -------
    bpix : ndarray
      0 = ok; 1. = Mask

    """
    bpix = np.zeros((4096, 2048))
    if det == 1:
        bpix[:,1052:1054] = 1.
    elif det == 2:
        bpix[:,0:4] = 1.
        bpix[:,376:380] = 1.
        bpix[:,2047] = 1.
    elif det == 3:
        bpix[:,851] = 1.
    elif det == 4:
        bpix[:,0:4] = 1.
        bpix[:,997:998] = 1.
    elif det == 5:
        bpix[:,129] = 1.
    elif det == 7:
        bpix[:,426:428] = 1.
    elif det == 8:
        bpix[:,931] = 1.
        bpix[:,933] = 1.

    return bpix


