""" Module for LRIS specific codes
"""
import glob
import os
import numpy as np
from astropy.io import fits

from pypeit import msgs
from pypeit.spectrographs import spectrograph
from ..par.pypeitpar import DetectorPar
from .. import telescopes
from pypeit.core import framematch
from pypeit.core import parse
from pypeit.par import pypeitpar

from IPython import embed

class GeminiGMOSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Gemini/GMOS specific code
    """

    def __init__(self):

        # Get it started
        super(GeminiGMOSSpectrograph, self).__init__()
        self.timeunit = 'isot'  # Synthesizes date+time

    def init_meta(self):
        """
        Generate the meta data dictionary.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(ext=0, card='MASKNAME')
        self.meta['binning'] = dict(card=None, compound=True)
        # TODO: Can we define the card here that compound meta uses to
        # set the binning?  Would be better to have all header cards
        # collected in this function...
#        self.meta['binning'] = dict(ext=1, card='CCDSUM')

        self.meta['mjd'] = dict(ext=0, card='OBSEPOCH')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='GRATING')
        self.meta['dispangle'] = dict(ext=0, card='CENTWAVE', rtol=1e-5)
        self.meta['dichroic'] = dict(ext=0, card='FILTER1')

        self.meta['datasec'] = dict(ext=1, card='DATASEC')

    def configuration_keys(self):
        """
        Extra keys for defining the configuration

        Returns:

        """
        return ['datasec']


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

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype == 'science':
            return good_exp & (fitstbl['target'] != 'CuAr') & (fitstbl['target'] != 'GCALflat') & (fitstbl['target'] != 'Bias')
            #& (fitstbl['idname'] == 'OBJECT')
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['target'] == 'CuAr')#& (fitstbl['idname'] == 'ARC')
        if ftype in ['pixelflat', 'trace']:
            return good_exp & (fitstbl['target'] == 'GCALflat')#& (fitstbl['idname'] == 'FLAT')
        if ftype == 'bias':
            return good_exp & (fitstbl['target'] == 'Bias')#& (fitstbl['idname'] == 'BIAS')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Keck LRISb reductions.
        """
        par = pypeitpar.PypeItPar()
        par['calibrations']['slitedges']['edge_thresh'] = 20.
        par['calibrations']['slitedges']['fit_order'] = 3

        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.40  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect'] = 5.  # Doesn't work for reddest chip
        par['calibrations']['wavelengths']['lamps'] = ['CuI', 'ArI', 'ArII']
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['nsnippet'] = 1  # 3 detectors splitting is already a lot

        par['calibrations']['tilts']['tracethresh'] = 10.  # Deals with faint CuAr lines

        #   IF YOU CHANGE THIS, YOU WILL NEED TO DEAL WITH THE OVERSCAN GOING ALONG ROWS
        #for key in par['calibrations'].keys():
        #    if 'frame' in key:
        #        par['calibrations'][key]['process']['overscan'] = 'median'

        # Overscan subtract the images
        #par['calibrations']['biasframe']['useframe'] = 'overscan'

        # Alter the method used to combine pixel flats
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['pixelflatframe']['process']['sig_lohi'] = [10.,10.]

        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        # Always correct for flexure, starting with default parameters
        par['flexure']['method'] = 'boxcar'

        # Set the default exposure time ranges for the frame typing
        #par['scienceframe']['exprng'] = [30, None]

        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        .. todo::
            Document the changes made!

        Args:
            scifile (str):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`pypeit.par.parset.ParSet`: The PypeIt paramter set
            adjusted for configuration specific parameter values.
        """
        par = self.__class__.default_pypeit_par() if inp_par is None else inp_par

        headarr = self.get_headarr(scifile)

        # Turn PCA off for long slits
        if 'arcsec' in self.get_meta_value(headarr, 'decker'):
            par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        return par

    def configuration_keys(self):
        """
        Return the metadata keys that defines a unique instrument
        configuration.

        This list is used by :class:`pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:

            list: List of keywords of data pulled from meta
        """
        cfg_keys = super(GeminiGMOSSpectrograph, self).configuration_keys()
        # Add grating tilt
        return cfg_keys+['dispangle']

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
        msgs.info("Reading GMOS file: {:s}".format(fil[0]))
        hdu = fits.open(fil[0])
        head0 = hdu[0].header
        head1 = hdu[1].header

        # Number of amplifiers (could pull from DetectorPar but this avoids needing the spectrograph, e.g. view_fits)
        numamp = (len(hdu) - 1) // 3

        # get the x and y binning factors...
        binning = head1['CCDSUM']
        xbin, ybin = [int(ibin) for ibin in binning.split(' ')]

        # First read over the header info to determine the size of the output array...
        datasec = head1['DATASEC']
        x1, x2, y1, y2 = np.array(parse.load_sections(datasec, fmt_iraf=False)).flatten()
        biassec = head1['BIASSEC']
        b1, b2, b3, b4 = np.array(parse.load_sections(biassec, fmt_iraf=False)).flatten()
        nxb = b2 - b1 + 1

        # determine the output array size...
        nx = (x2 - x1 + 1) * numamp + nxb * numamp
        ny = y2 - y1 + 1

        # allocate output array...
        array = np.zeros((nx, ny))
        rawdatasec_img = np.zeros_like(array, dtype=int)
        oscansec_img = np.zeros_like(array, dtype=int)

        if numamp == 2:  # E2V
            if det == 1:  # BLUEST DETECTOR
                order = range(6, 4, -1)
            elif det == 2:  # NEXT
                order = range(3, 5)
            elif det == 3:  # REDDEST DETECTOR
                order = range(1, 3)
        elif numamp == 4:  # Hamamatsu
            if det == 1:  # BLUEST DETECTOR
                order = range(12, 8, -1)
            elif det == 2:  # BLUEST DETECTOR
                order = range(8, 4, -1)
            elif det == 3:  # BLUEST DETECTOR
                order = range(4, 0, -1)
        else:
            embed()

        # insert extensions into master image...
        for kk, jj in enumerate(order):
            # grab complete extension...
            data, overscan, datasec, biassec, x1, x2 = gemini_read_amp(hdu, jj)
            # insert components into output array...
            inx = data.shape[0]
            xs = inx * kk
            xe = xs + inx

            # insert data...
            # Data section
            #section = '[{:d}:{:d},:]'.format(xs * xbin, xe * xbin)  # Eliminate lines
            #dsec.append(section)
            array[xs:xe, :] = np.flipud(data)
            rawdatasec_img[xs:xe, :] = kk+1

            # ; insert postdata...
            xs = nx - numamp * nxb + kk * nxb
            xe = xs + nxb

            #osection = '[{:d}:{:d},:]'.format(xs * xbin, xe * xbin)  # TRANSPOSED FOR WHAT COMES
            #osec.append(osection)
            array[xs:xe, :] = overscan
            oscansec_img[xs:xe, :] = kk+1

        # Need the exposure time
        exptime = hdu[self.meta['exptime']['ext']].header[self.meta['exptime']['card']]
        # Return, transposing array back to orient the overscan properly
        return array.T, hdu, exptime, rawdatasec_img.T, oscansec_img.T


class GeminiGMOSSHamSpectrograph(GeminiGMOSSpectrograph):
    """
    Child to handle Gemini/GMOS-S instrument with Hamamatsu detector
    """
    def __init__(self):

        # Get it started
        super(GeminiGMOSSHamSpectrograph, self).__init__()
        self.spectrograph = 'gemini_gmos_south_ham'
        self.camera = 'GMOS-S'
        self.telescope = telescopes.GeminiSTelescopePar()

        self.detector = [  #  Hamamatsu (since 2014)
            # Detector 1
            DetectorPar(dataext         = 1,  # Not sure this is used
                        specaxis        = 1,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.080,
                        darkcurr        = 0.0,
                        saturation      = 129000.,
                        nonlinear       = 0.95,
                        numamplifiers   = 4,
                        gain            = [1.83]*4,
                        ronoise         = [3.98]*4,
                        suffix          = '_01'
                        ),
            # Detector 2
            DetectorPar(dataext         = 2,  # Not sure this is used
                        specaxis        = 1,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.080,
                        darkcurr        = 0.0,
                        saturation      = 123000.,
                        nonlinear       = 0.95,
                        numamplifiers   = 4,
                        gain            = [1.83]*4,
                        ronoise         = [3.98]*4,
                        suffix          = '_02'
                        ),
            # Detector 3
            DetectorPar(dataext         = 3,  # Not sure this is used
                        specaxis        = 1,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.080,
                        darkcurr        = 0.0,
                        saturation      = 125000.,
                        nonlinear       = 0.95,
                        numamplifiers   = 4,
                        gain            = [1.83]*4,
                        ronoise         = [3.98]*4,
                        suffix          = '_03'
                        ),
        ]
        self.numhead = 13

    def bpm(self, filename, det, shape=None):
        """ Generate a BPM

        Parameters
        ----------
        det : int, REQUIRED
        **null_kwargs:
           Captured and never used

        Returns
        -------
        badpix : ndarray

        """
        # Get the empty bpm: force is always True
        bpm_img = self.empty_bpm(filename, det, shape=shape)

        if det == 1:
            msgs.info("Using hard-coded BPM for det=1 on GMOSs")

            # TODO: Fix this
            # Get the binning
            hdu = fits.open(filename)
            binning = hdu[1].header['CCDSUM']
            hdu.close()

            # Apply the mask
            xbin = int(binning.split(' ')[0])
            badc = 616//xbin
            bpm_img[badc,:] = 1
        elif det == 2:
            msgs.info("Using hard-coded BPM for det=2 on GMOSs")

            # Get the binning
            hdu = fits.open(filename)
            binning = hdu[1].header['CCDSUM']
            hdu.close()

            # Apply the mask
            xbin = int(binning.split(' ')[0])
            if xbin != 2:
                embed()
            # Up high
            badr = (898*2)//xbin # Transposed
            bpm_img[badr:badr+(8*2)//xbin,:] = 1
            # Down low
            badr = (161*2)//xbin # Transposed
            bpm_img[badr,:] = 1
        elif det == 3:
            msgs.info("Using hard-coded BPM for det=3 on GMOSs")

            # Get the binning
            hdu = fits.open(filename)
            binning = hdu[1].header['CCDSUM']
            hdu.close()

            # Apply the mask
            xbin = int(binning.split(' ')[0])
            if xbin != 2:
                embed()
            badr = (281*2)//xbin # Transposed
            bpm_img[badr:badr+(2*2)//xbin,:] = 1

        return bpm_img

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        .. todo::
            Document the changes made!

        Args:
            scifile (str):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`pypeit.par.parset.ParSet`: The PypeIt paramter set
            adjusted for configuration specific parameter values.
        """
        # Start with instrument wide
        par = super(GeminiGMOSSHamSpectrograph, self).config_specific_par(scifile, inp_par=inp_par)

        if self.get_meta_value(scifile, 'dispname')[0:4] == 'R400':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gmos_r400_ham.fits'
        elif self.get_meta_value(scifile, 'dispname')[0:4] == 'B600':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gmos_b600_ham.fits'
        #
        return par


class GeminiGMOSNSpectrograph(GeminiGMOSSpectrograph):
    """
    Child to handle Gemini/GMOS-N instrument
    """
    def __init__(self):

        # Get it started
        super(GeminiGMOSNSpectrograph, self).__init__()
        self.telescope = telescopes.GeminiNTelescopePar()
        self.camera = 'GMOS-N'


class GeminiGMOSNHamSpectrograph(GeminiGMOSNSpectrograph):
    """
    Child to handle Gemini/GMOS-N instrument with Hamamatsu detector
    Used since February 2017
    """
    def __init__(self):

        # Get it started
        super(GeminiGMOSNHamSpectrograph, self).__init__()
        self.spectrograph = 'gemini_gmos_north_ham'

        self.detector = [  #  Hamamatsu
            # Detector 1
            DetectorPar(dataext         = 1,  # Not sure this is used
                        specaxis        = 1,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.0807,
                        darkcurr        = 0.0,
                        saturation      = 129000.,
                        nonlinear       = 0.95,
                        numamplifiers   = 4,
                        gain            = [1.63]*4,
                        ronoise         = [4.14]*4,
                        suffix          = '_01'
                        ),
            # Detector 2
            DetectorPar(dataext         = 2,  # Not sure this is used
                        specaxis        = 1,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.0807,
                        darkcurr        = 0.0,
                        saturation      = 123000.,
                        nonlinear       = 0.95,
                        numamplifiers   = 4,
                        gain            = [1.63]*4,
                        ronoise         = [4.14]*4,
                        suffix          = '_02'
                        ),
            # Detector 3
            DetectorPar(dataext         = 3,  # Not sure this is used
                        specaxis        = 1,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.0807,
                        darkcurr        = 0.0,
                        saturation      = 125000.,
                        nonlinear       = 0.95,
                        numamplifiers   = 4,
                        gain            = [1.63]*4,
                        ronoise         = [4.14]*4,
                        suffix          = '_03'
                        ),
        ]
        self.numhead = 13

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        .. todo::
            Document the changes made!

        Args:
            scifile (str):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`pypeit.par.parset.ParSet`: The PypeIt paramter set
            adjusted for configuration specific parameter values.
        """
        # Start with instrument wide
        par = super(GeminiGMOSNHamSpectrograph, self).config_specific_par(scifile, inp_par=inp_par)

        if self.get_meta_value(scifile, 'dispname')[0:4] == 'R400':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gmos_r400_ham.fits'
        elif self.get_meta_value(scifile, 'dispname')[0:4] == 'B600':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gmos_b600_ham.fits'
        #
        return par

class GeminiGMOSNE2VSpectrograph(GeminiGMOSNSpectrograph):
    """
    Child to handle Gemini/GMOS-N instrument with E2V detector
    Used until February 2017
    """
    def __init__(self):

        # Get it started
        super(GeminiGMOSNE2VSpectrograph, self).__init__()

        self.spectrograph = 'gemini_gmos_north_e2v'

        self.detector = [  #  E2V
            # Detector 1
            DetectorPar(dataext         = 1,  # Not sure this is used
                        specaxis        = 1,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.0728,  # arcsec per pixel
                        darkcurr        = 0.0,
                        saturation      = 110900.,
                        nonlinear       = 0.95,
                        numamplifiers   = 2,
                        gain            = [2.27]*2,
                        ronoise         = [3.32]*2,
                        suffix          = '_01'
                        ),
            # Detector 2
            DetectorPar(dataext         = 2,  # Not sure this is used
                        specaxis        = 1,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.0728,
                        darkcurr        = 0.0,
                        saturation      = 115500.,
                        nonlinear       = 0.95,
                        numamplifiers   = 2,
                        gain            = [2.27]*2,
                        ronoise         = [3.32]*2,
                        suffix          = '_02'
                        ),
            # Detector 3
            DetectorPar(dataext         = 3,  # Not sure this is used
                        specaxis        = 1,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.0728,
                        darkcurr        = 0.0,
                        saturation      = 116700.,
                        nonlinear       = 0.95,
                        numamplifiers   = 2,
                        gain            = [2.27]*2,
                        ronoise         = [3.32]*2,
                        suffix          = '_03'
                        ),
        ]
        self.numhead = 7

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        .. todo::
            Document the changes made!

        Args:
            scifile (str):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`pypeit.par.parset.ParSet`: The PypeIt paramter set
            adjusted for configuration specific parameter values.
        """
        # Start with instrument wide
        par = super(GeminiGMOSNE2VSpectrograph, self).config_specific_par(scifile, inp_par=inp_par)

        if self.get_meta_value(scifile, 'dispname')[0:4] == 'R400':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gmos_r400_e2v.fits'
        #
        return par

def gemini_read_amp(inp, ext):
    """
    Read one amplifier of an Gemini GMOS multi-extension FITS image

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

    # parse the DETSEC keyword to determine the size of the array.
    header = hdu[ext].header
    detsec = header['DETSEC']
    x1, x2, y1, y2 = np.array(parse.load_sections(detsec, fmt_iraf=False)).flatten()

    # parse the DATASEC keyword to determine the size of the science region (unbinned)
    datasec = header['DATASEC']
    xdata1, xdata2, ydata1, ydata2 = np.array(parse.load_sections(datasec, fmt_iraf=False)).flatten()

    # grab the components...
    data = temp[xdata1-1:xdata2,:]

    # Overscan
    biassec = header['BIASSEC']
    xdata1, xdata2, ydata1, ydata2 = np.array(parse.load_sections(biassec, fmt_iraf=False)).flatten()
    overscan = temp[xdata1-1:xdata2,:]

    # Return
    return data, overscan, datasec, biassec, x1, x2



