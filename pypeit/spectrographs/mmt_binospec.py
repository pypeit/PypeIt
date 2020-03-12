""" Module for MMT/BINOSPEC specific codes
"""
import glob
import numpy as np
from astropy.io import fits

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.core import parse

## FW ToDo: I subtract the overscan when reading data rather than pasrse it to procimg.py since procimg does not
##    support subtracting overscan in both x and y-axis

class MMTBINOSPECSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle MMT/BINOSPEC specific code
    """
    def __init__(self):
        # Get it started
        super(MMTBINOSPECSpectrograph, self).__init__()
        self.spectrograph = 'mmt_binospec'
        self.telescope = telescopes.MMTTelescopePar()
        self.camera = 'BINOSPEC'
        self.numhead = 11
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 1,
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.24,
                            darkcurr        = 3.0, ##ToDO: To Be update
                            saturation      = 65535.,
                            nonlinear       = 0.95,  #ToDO: To Be update
                            numamplifiers   = 4,
                            gain            = [1.085,1.046,1.042,0.975],
                            ronoise         = [3.2]*4,
                            suffix          = '_01'
                            ),
                # Detector 2
                pypeitpar.DetectorPar(
                            dataext         = 2,
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.24,
                            darkcurr        = 3.0, ##ToDO: To Be update
                            saturation      = 65535.,
                            nonlinear       = 0.95, #ToDO: To Be update
                            numamplifiers   = 4,
                            gain            = [1.028,1.163,1.047,1.045],
                            ronoise         = [3.2]*4,
                            suffix          = '_02'
                )]

    def init_meta(self):
        """
        Generate the meta data dict
        Note that the children can add to this

        Returns:
            self.meta: dict (generated in place)

        """
        meta = {}
        # Required (core)
        meta['ra'] = dict(ext=1, card='RA')
        meta['dec'] = dict(ext=1, card='DEC')
        meta['target'] = dict(ext=1, card='OBJECT')
        meta['decker'] = dict(ext=1, card=None, default='default')
        meta['binning'] = dict(card=None, compound=True)

        meta['mjd'] = dict(ext=1, card='MJD')
        meta['exptime'] = dict(ext=1, card='EXPTIME')
        meta['airmass'] = dict(ext=1, card='AIRMASS')
        # Extras for config and frametyping
        meta['dispname'] = dict(ext=1, card='DISPERS1')
        meta['idname'] = dict(ext=1, card='IMAGETYP')

        # used for arclamp
        meta['lampstat01'] = dict(ext=1, card='HENEAR')
        # used for flatlamp, SCRN is actually telescope status
        meta['lampstat02'] = dict(ext=1, card='SCRN')

        # Ingest
        self.meta = meta


    def compound_meta(self, headarr, meta_key):
        """

        Args:
            headarr: list
            meta_key: str

        Returns:
            value

        """
        if meta_key == 'binning':
            binspatial, binspec = parse.parse_binning(headarr[1]['CCDSUM'])
            binning = parse.binning2string(binspec, binspatial)
            return binning

    def default_pypeit_par(self):
        """
        Set default parameters for MMT/BINOSPEC reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'mmt_binospec'
        # Frame numbers
        par['calibrations']['standardframe']['number'] = 0
        par['calibrations']['biasframe']['number'] = 0
        par['calibrations']['pixelflatframe']['number'] = 5
        par['calibrations']['traceframe']['number'] = 5
        par['calibrations']['arcframe']['number'] = 5
        par['calibrations']['arcframe']['process']['overscan'] ='median'
        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.5
        par['calibrations']['wavelengths']['sigdetect'] = 20.
        par['calibrations']['wavelengths']['fwhm']= 5.0
        par['calibrations']['wavelengths']['lamps'] = ['ArI', 'ArII']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'

        # Tilt and slit parameters
        par['calibrations']['tilts']['tracethresh'] =  10.0
        par['calibrations']['tilts']['spat_order'] = 6
        par['calibrations']['tilts']['spec_order'] = 6
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Flats
        par['calibrations']['flatfield']['illumflatten'] = True

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0
        ## Do not perform global sky subtraction for standard stars
        ## FW: The slit is too wide and the target usually close to the center where is the boundary
        ## two detectors. global_sky will cause some problem sometime.
        par['reduce']['skysub']['global_sky_std']  = False

        # Flexure
        par['flexure']['method'] = 'skip'

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] ='nothing'
        #par['scienceframe']['process']['overscan'] ='median'

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 100]
        par['calibrations']['arcframe']['exprng'] = [20, None]
        par['calibrations']['darkframe']['exprng'] = [20, None]
        par['scienceframe']['exprng'] = [20, None]
        return par

    def configuration_keys(self):
        return ['dispname']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype == 'science':
            return good_exp & (fitstbl['lampstat01'] == 'off') & (fitstbl['lampstat02'] == 'stowed') & (fitstbl['exptime'] > 100.0)
        if ftype == 'standard':
            return good_exp & (fitstbl['lampstat01'] == 'off') & (fitstbl['lampstat02'] == 'stowed') & (fitstbl['exptime'] <= 100.0)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['lampstat01'] == 'on')
        if ftype in ['pixelflat', 'trace']:
            return good_exp & (fitstbl['lampstat01'] == 'off') & (fitstbl['lampstat02'] == 'deployed')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


    def get_rawimage(self, raw_file, det):
        """
        Load up the raw image and generate a few other bits and pieces
        that are key for image processing

        Args:
            raw_file (str):
            det (int):

        Returns:
            tuple:
                raw_img (np.ndarray) -- Raw image for this detector
                hdu (astropy.io.fits.HDUList)
                exptime (float)
                rawdatasec_img (np.ndarray)
                oscansec_img (np.ndarray)

        """
        # Check for file; allow for extra .gz, etc. suffix
        fil = glob.glob(raw_file + '*')
        if len(fil) != 1:
            msgs.error("Found {:d} files matching {:s}".format(len(fil)))

        # Read
        msgs.info("Reading BINOSPEC file: {:s}".format(fil[0]))
        hdu = fits.open(fil[0])
        head0 = hdu[0].header
        head1 = hdu[1].header

        # TOdO Store these parameters in the DetectorPar.
        # Number of amplifiers (could pull from DetectorPar but this avoids needing the spectrograph, e.g. view_fits)
        numamp = 4

        # get the x and y binning factors...
        binning = head1['CCDSUM']
        xbin, ybin = [int(ibin) for ibin in binning.split(' ')]

        # First read over the header info to determine the size of the output array...
        datasec = head1['DATASEC']
        x1, x2, y1, y2 = np.array(parse.load_sections(datasec, fmt_iraf=False)).flatten()
        nxb = x1 - 1

        # determine the output array size...
        nx = (x2 - x1 + 1) * int(numamp/2) + nxb * int(numamp/2)
        ny = (y2 - y1 + 1) * int(numamp/2)

        #datasize = head1['DETSIZE']
        #_, nx, _, ny = np.array(parse.load_sections(datasize, fmt_iraf=False)).flatten()

        # allocate output array...
        array = np.zeros((nx, ny))
        rawdatasec_img = np.zeros_like(array, dtype=int)
        oscansec_img = np.zeros_like(array, dtype=int)

        if det == 1:  # A DETECTOR
            order = range(1, 5, 1)
        elif det == 2:  # B DETECTOR
            order = range(5, 9, 1)

        # insert extensions into master image...
        for kk, jj in enumerate(order):
            # grab complete extension...
            data, overscan, datasec, biassec = binospec_read_amp(hdu, jj)

            # insert components into output array...
            inx = data.shape[0]
            xs = inx * kk
            xe = xs + inx

            iny = data.shape[1]
            ys = iny * kk
            yn = ys + iny

            b1, b2, b3, b4 = np.array(parse.load_sections(biassec, fmt_iraf=False)).flatten()

            if kk == 0:
                array[b2:inx+b2,:iny] = data #* self.detector[det-1]['gain'][kk]
                rawdatasec_img[b2:inx+b2,:iny] = kk + 1
                array[:b2,:iny] = overscan
                oscansec_img[2:b2,:iny] = kk + 1
            elif kk == 1:
                array[b2+inx:2*inx+b2,:iny] = np.flipud(data) #* self.detector[det-1]['gain'][kk]
                rawdatasec_img[b2+inx:2*inx+b2:,:iny] = kk + 1
                array[2*inx+b2:,:iny] = overscan
                oscansec_img[2*inx+b2:,:iny] = kk + 1
            elif kk == 2:
                array[b2+inx:2*inx+b2,iny:] = np.fliplr(np.flipud(data)) #* self.detector[det-1]['gain'][kk]
                rawdatasec_img[b2+inx:2*inx+b2,iny:] = kk + 1
                array[2*inx+b2:, iny:] = overscan
                oscansec_img[2*inx+b2:, iny:] = kk + 1
            elif kk == 3:
                array[b2:inx+b2,iny:] = np.fliplr(data) #* self.detector[det-1]['gain'][kk]
                rawdatasec_img[b2:inx+b2,iny:] = kk + 1
                array[:b2,iny:] = overscan
                oscansec_img[2:b2,iny:] = kk + 1

        # Need the exposure time
        exptime = hdu[self.meta['exptime']['ext']].header[self.meta['exptime']['card']]
        # Return, transposing array back to orient the overscan properly
        return np.fliplr(np.flipud(array)), hdu, exptime, np.fliplr(np.flipud(rawdatasec_img)), np.fliplr(np.flipud(oscansec_img))


def binospec_read_amp(inp, ext):
    """
    Read one amplifier of an MMT BINOSPEC multi-extension FITS image

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
    if isinstance(inp, str):
        hdu = fits.open(inp)
    else:
        hdu = inp
    # get entire extension...
    temp = hdu[ext].data.transpose()
    tsize = temp.shape
    nxt = tsize[0]
    nyt = tsize[1]

    # parse the DETSEC keyword to determine the size of the array.
    header = hdu[ext].header
    #detsec = header['DETSEC']
    #x1, x2, y1, y2 = np.array(parse.load_sections(detsec, fmt_iraf=False)).flatten()

    # parse the DATASEC keyword to determine the size of the science region (unbinned)
    datasec = header['DATASEC']
    xdata1, xdata2, ydata1, ydata2 = np.array(parse.load_sections(datasec, fmt_iraf=False)).flatten()
    datasec = '[{:}:{:},{:}:{:}]'.format(xdata1 - 1, xdata2, ydata1-1, ydata2)

    # Overscan X-axis
    # Since pypeit can only subtract overscan along one axis, I'm subtract the overscan here using median method.
    if xdata1 > 1:
        overscanx = temp[0:xdata1-1, :]
        overscanx_vec = np.median(overscanx, axis=0)
        temp = temp - overscanx_vec[None,:]
    data = temp[xdata1 - 1:xdata2, ydata1 -1 : ydata2]

    # Overscan Y-axis
    if ydata2<nyt:
        os1, os2 = ydata2, nyt-1
        overscany = temp[xdata1 - 1:xdata2, ydata2:os2]
        overscany_vec = np.median(overscany, axis=1)
        data = data -  overscany_vec[:,None]

    biassec = '[0:{:},{:}:{:}]'.format(xdata1-1, ydata1-1, ydata2)
    xos1, xos2, yos1, yos2 = np.array(parse.load_sections(biassec, fmt_iraf=False)).flatten()
    overscan = np.zeros_like(temp[xos1:xos2, yos1:yos2]) # Give a zero fake overscan at the edge of each amplifiers

    return data, overscan, datasec, biassec


'''
    # Subtract overscan along y-axis
    if ydata2<nyt-2:
        oscany = np.median(temp[:, ydata2:nyt],axis=1)
        #from scipy import signal
        #ossub = signal.savgol_filter(oscany, 65, 5)
        ossub = oscany.copy()
        temp = temp - ossub[:, None]
    # grab the components...
    data = temp[xdata1 - 1:xdata2, ydata1 -1 : ydata2]

    # Overscan
    biassec = '[0:{:},{:}:{:}]'.format(xdata1-1, ydata1-1, ydata2)
    xdata1, xdata2, ydata1, ydata2 = np.array(parse.load_sections(biassec, fmt_iraf=False)).flatten()
    from IPython import embed
    embed()
    ## ToDO: Figure out the real overscan along x-axis
    overscan = temp[xdata1:xdata2, ydata1:ydata2]
    #overscan =np.zeros_like(temp[xdata1:xdata2, ydata1:ydata2]) + np.median(temp[xdata1:xdata2,:np.max([nyt-ydata2,50])])

    # Return
    return data, overscan, datasec, biassec


'''