""" Module for LRIS specific codes
"""
from __future__ import absolute_import, division, print_function

try:
    basestring
except NameError:  # For Python 3
    basestring = str

import glob

import numpy as np
from astropy.io import fits

from pypit import msgs
from pypit import arparse
from pypit import ardebug as debugger
from pypit.spectrographs import spectrograph
from ..par.pypitpar import DetectorPar
from .. import telescopes

class GeminiGMOSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Gemini/GMOS specific code
    """

    def __init__(self):

        # Get it started
        super(GeminiGMOSSpectrograph, self).__init__()
        self.spectrograph = 'gemini_gmos'
        self.telescope = telescopes.GeminiTelescopePar()

    def load_raw_img_head(self, raw_file, det=None, **null_kwargs):
        """
        Wrapper to the raw image reader for LRIS

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
        raw_img, head0, _ = read_gmos(raw_file, det=det)

        return raw_img, head0

    '''
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
        temp, head0, secs = read_gmos(filename, det)
        for kk in range(settings_det['numamplifiers']):
            datasec.append(arparse.load_sections(secs[0][kk], fmt_iraf=False))
            oscansec.append(arparse.load_sections(secs[1][kk], fmt_iraf=False))

        # Need naxis0, naxis1 too
        naxis0 = temp.shape[0]
        naxis1 = temp.shape[1]

        # Return
        return datasec, oscansec, naxis0, naxis1
    '''

    def get_image_shape(self, filename=None, det=None, **null_kwargs):
        """
        Overrides :class:`Spectrograph.get_image_shape` for LRIS images.

        Must always provide a file.
        """
        # Cannot be determined without file
        if filename is None:
            raise ValueError('Must provide a file to determine the shape of an LRIS image.')

        # Use a file
        self._check_detector()
        self.naxis = (self.load_raw_frame(filename, det=det)[0]).shape
        return self.naxis

class GeminiGMOSSSpectrograph(GeminiGMOSSpectrograph):
    """
    Child to handle Keck/LRISb specific code
    """
    def __init__(self):

        # Get it started
        super(GeminiGMOSSSpectrograph, self).__init__()
        self.spectrograph = 'gemini_gmos_south'
        self.camera = 'GMOS-S'
        self.detector = [
            # Detector 1
            DetectorPar(dataext         = 1,  # Not sure this is used
                        dispaxis        = 1,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.080,
                        darkcurr        = 0.0,
                        saturation      = 65535.,
                        nonlinear       = 0.86,
                        numamplifiers   = 4,
                        gain            = [1.83]*4,
                        ronoise         = [3.98]*4,
                        suffix          = '_01'
                        ),
            # Detector 2
            DetectorPar(dataext         = 2,  # Not sure this is used
                        dispaxis        = 1,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.080,
                        darkcurr        = 0.0,
                        saturation      = 65535.,
                        nonlinear       = 0.86,
                        numamplifiers   = 4,
                        gain            = [1.83]*4,
                        ronoise         = [3.98]*4,
                        suffix          = '_02'
                        ),
            # Detector 3
            DetectorPar(dataext         = 3,  # Not sure this is used
                        dispaxis        = 1,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.080,
                        darkcurr        = 0.0,
                        saturation      = 65535.,
                        nonlinear       = 0.86,
                        numamplifiers   = 4,
                        gain            = [1.83]*4,
                        ronoise         = [3.98]*4,
                        suffix          = '_03'
                        ),
        ]

    def setup_arcparam(self, arcparam, disperser=None, **null_kwargs):
        """
        Setup the arc parameters

        Args:
            arcparam: dict
            disperser: str, REQUIRED
            **null_kwargs:
              Captured and never used

        Returns:
            arcparam is modified in place

        """
        arcparam['lamps'] = ['CuI', 'ArI']
        if 'R150' in disperser:
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=0.63 # Ang per pixel (unbinned)
            arcparam['b1']= 4.54698031e-04
            arcparam['b2']= -6.86414978e-09
            arcparam['wvmnx'][1] = 6000.
            arcparam['wv_cen'] = 4000.
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))


def read_gmos(raw_file, det=1):
    """

    Parameters
    ----------
    raw_file : str
      Filename
    det : int, optional
      Detector number; Default = 1

    Returns
    -------
    array : ndarray
      Combined image 
    header : FITS header
    sections : list
      List of datasec, oscansec, ampsec sections
    """

    # Check for file; allow for extra .gz, etc. suffix
    fil = glob.glob(raw_file+'*') 
    if len(fil) != 1:
        msgs.error("Found {:d} files matching {:s}".format(len(fil)))

    # Read
    msgs.info("Reading LRIS file: {:s}".format(fil[0]))
    hdu = fits.open(fil[0])
    head0 = hdu[0].header
    head1 = hdu[1].header

    # Setup for datasec, oscansec
    dsec = []
    osec = []

    # get the x and y binning factors...
    binning = head1['CCDSUM']
    xbin, ybin = [int(ibin) for ibin in binning.split(' ')]

    # First read over the header info to determine the size of the output array...
    n_ext = len(hdu)-1  # Number of extensions (usually 12)

    datasec = head1['DATASEC']
    x1, x2, y1, y2 = np.array(arparse.load_sections(datasec, fmt_iraf=False)).flatten()
    biassec = head1['BIASSEC']
    b1, b2, b3, b4 = np.array(arparse.load_sections(biassec, fmt_iraf=False)).flatten()
    nxb = b2-b1 + 1

    # determine the output array size...
    nx = x2*4 + nxb*4
    ny = y2-y1+1

    # Deal with detectors
    if det in [1,2,3,4]:
        n_ext = n_ext // 3
        det_idx = np.arange(n_ext, dtype=np.int) + (det-1)*n_ext
        ndet = 1
        order = range((det-1)*4+1,(det-1)*4+5)
    elif det is 'all':
        debugger.set_trace() # NEED TO SET THIS UP
        ndet = 2
        det_idx = np.arange(n_ext).astype(int)
    else:
        raise ValueError('Bad value for det')

    # allocate output array...
    array = np.zeros( (nx, ny) )

    # insert extensions into master image...
    for kk, i in enumerate(order):

        # grab complete extension...
        data, overscan, datasec, biassec, x1, x2 = gemini_read_amp(hdu, i+1)
                            #, linebias=linebias, nobias=nobias, $
                            #x1=x1, x2=x2, y1=y1, y2=y2, gaindata=gaindata)
        # insert components into output array...
        inx = data.shape[0]
        xs = inx*kk
        xe = xs + inx

        # insert data...
        # Data section
        section = '[:,{:d}:{:d}]'.format(xs, xe)  # Eliminate lines
        dsec.append(section)
        array[xs:xe, :] = data

        #; insert postdata...
        xs = nx - n_ext*nxb + kk*nxb
        xe = xs + nxb
        section = '[:,{:d}:{:d}]'.format(xs, xe)
        osec.append(section)
        array[xs:xe, :] = overscan

    # make sure BZERO is a valid integer for IRAF
    obzero = head1['BZERO']
    #head0['O_BZERO'] = obzero
    head0['BZERO'] = 32768-obzero

    # Return, transposing array back to goofy Python indexing
    return array, head0, (dsec, osec)


def gemini_read_amp(inp, ext):
    """
    Read one amplifier of an LRIS multi-extension FITS image

    Parameters
    ----------
    inp: tuple 
      (str,int) filename, extension
      (hdu,int) FITS hdu, extension

    Returns
    -------
    data
    predata
    postdata
    x1
    y1

    ;------------------------------------------------------------------------
    function lris_read_amp, filename, ext, $
      linebias=linebias, nobias=nobias, $
      predata=predata, postdata=postdata, header=header, $
      x1=x1, x2=x2, y1=y1, y2=y2, GAINDATA=gaindata
    ;------------------------------------------------------------------------
    ; Read one amp from LRIS mHDU image
    ;------------------------------------------------------------------------
    """
    # Parse input
    if isinstance(inp, basestring):
        hdu = fits.open(inp)
    else:
        hdu = inp

    head1 = hdu[1].header

    # Deal with binning
    binning = head1['CCDSUM']
    xbin, ybin = [int(ibin) for ibin in binning.split(' ')]

    # get entire extension...
    temp = hdu[ext].data.transpose()
    tsize = temp.shape
    nxt = tsize[0]

    # parse the DETSEC keyword to determine the size of the array.
    header = hdu[ext].header
    detsec = header['DETSEC']
    x1, x2, y1, y2 = np.array(arparse.load_sections(detsec, fmt_iraf=False)).flatten()

    # parse the DATASEC keyword to determine the size of the science region (unbinned)
    datasec = header['DATASEC']
    xdata1, xdata2, ydata1, ydata2 = np.array(arparse.load_sections(datasec, fmt_iraf=False)).flatten()

    # grab the components...
    data = temp[xdata1-1:xdata2,:]

    # Overscan
    biassec = header['BIASSEC']
    xdata1, xdata2, ydata1, ydata2 = np.array(arparse.load_sections(biassec, fmt_iraf=False)).flatten()
    overscan = temp[xdata1-1:xdata2,:]

    # Return
    return data, overscan, datasec, biassec, x1, x2



