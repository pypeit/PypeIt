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

class KeckDEIMOSSpectrograph(spectroclass.Spectrograph):
    """
    Child to handle Keck/DEIMOS specific code
    """

    def __init__(self):

        # Get it started
        spectroclass.Spectrograph.__init__(self)
        self.spectrograph = 'keck_deimos'

    def load_raw_img_head(self, raw_file, det=None, **null_kwargs):
        """
        Wrapper to the raw image reader for DEIMOS

        Args:
            raw_file:  str, filename
            det: int, REQUIRED
              Desired detector
            **null_kwargs:
              Captured and never used

        Returns:
            raw_img: ndarray
              Raw image;  likely unsigned int
            head0: Header

        """
        raw_img, head0, _ = read_deimos(raw_file, det=det)

        return raw_img, head0

    def get_datasec(self, filename, det, settings_det):
        """
        Load up the datasec and oscansec and also naxis0 and naxis1

        Args:
            filename: str
              data filename
            det: int
              Detector specification
            settings_det: ParSet
              numamplifiers

        Returns:
            datasec: list
            oscansec: list
            naxis0: int
            naxis1: int
        """

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

    def bpm(self, det=None, **null_kwargs):
        """ Generate a BPM for DEIMOS
        Currently assumes 1x1 binning

        Parameters
        ----------
        det : int, REQUIRED
        **null_kwargs:
            Captured and never used

        Returns
        -------
        bpix : ndarray
          0 = ok; 1 = Mask

        """
        bpix = np.zeros((4096, 2048), dtype=int)
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

    def setup_arcparam(self, arcparam, disperser=None, fitstbl=None, arc_idx=None,
                       msarc_shape=None, **null_kwargs):
        """

        Args:
            arcparam:
            disperser:
            fitstbl:
            arc_idx:
            msarc_shape:
            binspectral:
            **null_kwargs:

        Returns:

        """
        arcparam['wv_cen'] = fitstbl['dispangle'][arc_idx]
        # TODO -- Should set according to the lamps that were on
        arcparam['lamps'] = ['ArI','NeI','KrI','XeI']
        if disperser == '830G': # Blaze 8640
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=0.47 # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0]
            arcparam['wvmnx'][0] = 550.
            arcparam['wvmnx'][1] = 11000.
            arcparam['min_ampl'] = 3000.  # Lines tend to be very strong
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))

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


