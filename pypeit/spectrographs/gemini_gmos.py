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
        # Set wave tilts order
        par['calibrations']['slits']['sigdetect'] = 20.
        par['calibrations']['slits']['trace_npoly'] = 3
        # TODO: No longer a parameter
#        par['calibrations']['slits']['fracignore'] = 0.02
#        par['calibrations']['slits']['pcapar'] = [3,2,1,0]

        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.40  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect'] = 5.  # Doesn't work for reddest chip
        par['calibrations']['wavelengths']['lamps'] = ['CuI', 'ArI', 'ArII']
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['nsnippet'] = 1  # 3 detectors splitting is already a lot

        # Overscan subtract the images
        #par['calibrations']['biasframe']['useframe'] = 'overscan'

        # Alter the method used to combine pixel flats
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['pixelflatframe']['process']['sig_lohi'] = [10.,10.]

        # Scienceimage default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()
        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()

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
        par = self.default_pypeit_par() if inp_par is None else inp_par
        # TODO: Should we allow the user to override these?
        if self.get_meta_value(scifile, 'dispname')[0:4] == 'R400':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gmos_r400.fits'
        #
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

    def load_raw_frame(self, raw_file, det=None, **null_kwargs):
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
            hdu: HDUList

        """
        raw_img, hdu, _ = read_gmos(raw_file, det=det)

        return raw_img, hdu

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

    def get_image_section(self, inp=None, det=1, section='datasec'):
        """
        Return a string representation of a slice defining a section of
        the detector image.

        Overwrites base class function to use :func:`read_gmos` to get
        the image sections.

        .. todo::
            - It feels really ineffiecient to just get the image section
              using the full :func:`read_gmos`.  Can we parse that
              function into something that can give you the image
              section directly?

        This is done separately for the data section and the overscan
        section in case one is defined as a header keyword and the other
        is defined directly.

        Args:
            inp (:obj:`str`, `astropy.io.fits.Header`_, optional):
                String providing the file name to read, or the relevant
                header object.  Default is None, meaning that the
                detector attribute must provide the image section
                itself, not the header keyword.
            det (:obj:`int`, optional):
                1-indexed detector number.
            section (:obj:`str`, optional):
                The section to return.  Should be either 'datasec' or
                'oscansec', according to the
                :class:`pypeitpar.DetectorPar` keywords.

        Returns:
            tuple: Returns three objects: (1) A list of string
            representations for the image sections, one string per
            amplifier.  The sections are *always* returned in PypeIt
            order: spectral then spatial.  (2) Boolean indicating if the
            slices are one indexed.  (3) Boolean indicating if the
            slices should include the last pixel.  The latter two are
            always returned as True following the FITS convention.
        """
        # Read the file
        if inp is None:
            msgs.error('Must provide Gemini GMOS file to get image section.')
        elif not os.path.isfile(inp):
            msgs.error('File {0} does not exist!'.format(inp))
        temp, head0, secs = read_gmos(inp, det=det)
        if section == 'datasec':
            return secs[0], False, False
        elif section == 'oscansec':
            # Need to flip these
            # TODO: What does the above mean?
            return secs[1], False, False
        else:
            raise ValueError('Unrecognized keyword: {0}'.format(section))



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
                        specaxis        = 0,  # FW: this should be 0 for gmos_ham
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
                        specaxis        = 0,
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
                        specaxis        = 0,
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

    def bpm(self, shape=None, filename=None, det=None, **null_kwargs):
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
        self.empty_bpm(shape=shape, filename=filename, det=det)

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
            self.bpm_img[badc,:] = 1
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
            self.bpm_img[badr:badr+(8*2)//xbin,:] = 1
            # Down low
            badr = (161*2)//xbin # Transposed
            self.bpm_img[badr,:] = 1
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
            self.bpm_img[badr:badr+(2*2)//xbin,:] = 1

        return self.bpm_img



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
    """
    def __init__(self):

        # Get it started
        super(GeminiGMOSNHamSpectrograph, self).__init__()
        self.spectrograph = 'gemini_gmos_north_ham'

        self.detector = [  #  Hamamatsu (since 2011)
            # Detector 1
            DetectorPar(dataext         = 1,  # Not sure this is used
                        specaxis        = 0,  # FW: this should be 0 for gmos_ham
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
                        specaxis        = 0,
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
                        specaxis        = 0,
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


class GeminiGMOSNE2VSpectrograph(GeminiGMOSNSpectrograph):
    """
    Child to handle Gemini/GMOS-N instrument with E2V detector
    """
    def __init__(self):

        # Get it started
        super(GeminiGMOSNE2VSpectrograph, self).__init__()

        self.spectrograph = 'gemini_gmos_north_e2v'

        self.detector = [  #  E2V
            # Detector 1
            DetectorPar(dataext         = 1,  # Not sure this is used
                        specaxis        = 0,  # I think this is ignored
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
                        specaxis        = 0,
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
                        specaxis        = 0,
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

    def init_meta(self):
        """
        Generate the meta data dictionary.
        """
        super(GeminiGMOSNE2VSpectrograph, self).init_meta()
        self.meta['exptime'] = dict(ext=0, card='EXPOSURE')


# TODO: Put this inside the class or abstract it.
def read_gmos(raw_file, det=1):
    """
    Read the GMOS data file

    Parameters
    ----------
    raw_file : str
      Filename
    detector_par : ParSet
      Needed for numamplifiers if not other things
    det : int, optional
      Detector number; Default = 1

    Returns
    -------
    array : ndarray
      Combined image 
    hdu : FITS HDUList
    sections : list
      List of datasec, oscansec, ampsec sections
    """

    # Check for file; allow for extra .gz, etc. suffix
    fil = glob.glob(raw_file+'*') 
    if len(fil) != 1:
        msgs.error("Found {:d} files matching {:s}".format(len(fil)))

    # Read
    msgs.info("Reading GMOS file: {:s}".format(fil[0]))
    hdu = fits.open(fil[0])
    head0 = hdu[0].header
    head1 = hdu[1].header

    # Number of amplifiers (could pull from DetectorPar but this avoids needing the spectrograph, e.g. view_fits)
    numamp = (len(hdu)-1)//3

    # Setup for datasec, oscansec
    dsec = []
    osec = []

    # get the x and y binning factors...
    binning = head1['CCDSUM']
    xbin, ybin = [int(ibin) for ibin in binning.split(' ')]

    # First read over the header info to determine the size of the output array...
    datasec = head1['DATASEC']
    x1, x2, y1, y2 = np.array(parse.load_sections(datasec, fmt_iraf=False)).flatten()
    biassec = head1['BIASSEC']
    b1, b2, b3, b4 = np.array(parse.load_sections(biassec, fmt_iraf=False)).flatten()
    nxb = b2-b1 + 1

    # determine the output array size...
    nx = (x2-x1+1)*numamp + nxb*numamp
    ny = y2-y1+1

    # allocate output array...
    array = np.zeros( (nx, ny) )

    if numamp == 2:
        if det == 1: # BLUEST DETECTOR
            order = range(6,4,-1)
        elif det == 2: # BLUEST DETECTOR
            order = range(3,5)
        elif det == 3: # BLUEST DETECTOR
            order = range(1,3)
    elif numamp == 4:
        if det == 1: # BLUEST DETECTOR
            order = range(12,8,-1)
        elif det == 2: # BLUEST DETECTOR
            order = range(8,4,-1)
        elif det == 3: # BLUEST DETECTOR
            order = range(4,0,-1)
    else:
        embed()

    # insert extensions into master image...
    for kk, jj in enumerate(order):

        # grab complete extension...
        data, overscan, datasec, biassec, x1, x2 = gemini_read_amp(hdu, jj)
                            #, linebias=linebias, nobias=nobias, $
                            #x1=x1, x2=x2, y1=y1, y2=y2, gaindata=gaindata)
        # insert components into output array...
        inx = data.shape[0]
        xs = inx*kk
        xe = xs + inx

        # insert data...
        # Data section
        #section = '[:,{:d}:{:d}]'.format(xs, xe)  # Eliminate lines
        section = '[{:d}:{:d},:]'.format(xs*xbin, xe*xbin)  # Eliminate lines
        dsec.append(section)
        array[xs:xe, :] = np.flipud(data)

        #; insert postdata...
        xs = nx - numamp*nxb + kk*nxb
        xe = xs + nxb
        #debugger.set_trace()
        #section = '[:,{:d}:{:d}]'.format(xs, xe)
        osection = '[{:d}:{:d},:]'.format(xs*xbin, xe*xbin)  # TRANSPOSED FOR WHAT COMES
        osec.append(osection)
        array[xs:xe, :] = overscan


    # make sure BZERO is a valid integer for IRAF
    obzero = head1['BZERO']
    #head0['O_BZERO'] = obzero
    head0['BZERO'] = 32768-obzero

    # Return, transposing array back to goofy Python indexing
    return array, hdu, (dsec, osec)


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



