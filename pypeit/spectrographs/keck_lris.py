""" Module for LRIS specific codes
"""
from __future__ import absolute_import, division, print_function

import glob

import numpy as np
from astropy.io import fits

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph

from pypeit import debugger

class KeckLRISSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/LRIS specific code
    """
    def __init__(self):
        # Get it started
        super(KeckLRISSpectrograph, self).__init__()
        self.spectrograph = 'keck_lris_base'
        self.telescope = telescopes.KeckTelescopePar()

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Keck LRISr reductions.
        """
        par = pypeitpar.PypeItPar()
        # Set wave tilts order
        par['calibrations']['slits']['sigdetect'] = 30.
        par['calibrations']['slits']['pcapar'] = [3,2,1,0]
        # 1D wavelengths
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20  # Might be grism dependent
        # Always sky subtract, starting with default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()
        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['pixelflatframe']['exprng'] = [None, 30]
        par['calibrations']['traceframe']['exprng'] = [None, 30]
        par['scienceframe']['exprng'] = [29, None]
        return par

    def header_keys(self):
        """
        Return a dictionary with the header keywords to read from the
        fits file.

        Returns:
            dict: A nested dictionary with the header keywords to read.
            The first level gives the extension to read and the second
            level gives the common name for header values that is passed
            on to the PypeItMetaData object.
        """
        hdr_keys = {}
        hdr_keys[0] = {}
        hdr_keys[1] = {}
        hdr_keys[2] = {}
        hdr_keys[3] = {}
        hdr_keys[4] = {}

        # Copied over defaults
        hdr_keys[0]['idname'] = 'OBSTYPE'
        hdr_keys[0]['time'] = 'MJD-OBS'
        #hdr_keys[0]['date'] = 'DATE'
        hdr_keys[0]['utc'] = 'UTC'
        hdr_keys[0]['ut'] = 'UT'
        hdr_keys[0]['ra'] = 'RA'
        hdr_keys[0]['dec'] = 'DEC'
        hdr_keys[0]['airmass'] = 'AIRMASS'
        hdr_keys[0]['binning'] = 'BINNING'
        hdr_keys[0]['decker'] = 'SLITNAME'
        hdr_keys[0]['dichroic'] = 'DICHNAME'

        hdr_keys[0]['target'] = 'TARGNAME'
        hdr_keys[0]['exptime'] = 'ELAPTIME'
        hdr_keys[0]['hatch'] = 'TRAPDOOR'
        hdr_keys[0]['dispname'] = 'GRANAME'
        hdr_keys[0]['dispangle'] = 'GRANGLE'
        hdr_keys[0]['wavecen'] = 'WAVELEN'
        hdr_keys[0]['spectrograph'] = 'INSTRUME'
        hdr_keys[1]['NAXIS01'] = 'NAXIS'
        hdr_keys[2]['NAXIS02'] = 'NAXIS'
        hdr_keys[3]['NAXIS03'] = 'NAXIS'
        hdr_keys[4]['NAXIS04'] = 'NAXIS'
        hdr_keys[1]['CCDGEOM'] = 'CCDGEOM'
        hdr_keys[1]['CCDNAME01'] = 'CCDNAME'
        hdr_keys[3]['CCDNAME02'] = 'CCDNAME'

        lamp_names = ['MERCURY', 'NEON', 'ARGON', 'CADMIUM', 'ZINC', 'KRYPTON', 'XENON',
                      'FEARGON', 'DEUTERI', 'FLAMP1', 'FLAMP2', 'HALOGEN']
        for kk,lamp_name in enumerate(lamp_names):
            hdr_keys[0]['lampstat{:02d}'.format(kk+1)] = lamp_name

        return hdr_keys

    def metadata_keys(self):
        return super(KeckLRISSpectrograph, self).metadata_keys() \
                    + ['binning', 'dichroic', 'dispangle']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype == 'science':
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == 'open')
        if ftype == 'bias':
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == 'closed')
        if ftype == 'pixelflat' or ftype == 'trace':
            # Flats and trace frames are typed together
            return good_exp & self.lamps(fitstbl, 'dome') & (fitstbl['hatch'] == 'open')
        if ftype == 'pinhole' or ftype == 'dark':
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype == 'arc':
            return good_exp & self.lamps(fitstbl, 'arcs') & (fitstbl['hatch'] == 'closed')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)
  
    def lamps(self, fitstbl, status):
        """
        Check the lamp status.

        Args:
            fitstbl (:obj:`astropy.table.Table`):
                The table with the fits header meta data.
            status (:obj:`str`):
                The status to check.  Can be `off`, `arcs`, or `dome`.
        
        Returns:
            numpy.ndarray: A boolean array selecting fits files that
            meet the selected lamp status.

        Raises:
            ValueError:
                Raised if the status is not one of the valid options.
        """
        if status == 'off':
            # Check if all are off
            return np.all(np.array([ (fitstbl[k] == 'off') | (fitstbl[k] == 'None')
                                        for k in fitstbl.keys() if 'lampstat' in k]), axis=0)
        if status == 'arcs':
            # Check if any arc lamps are on
            arc_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(1,9) ]
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys()
                                            if k in arc_lamp_stat]), axis=0)
        if status == 'dome':
            # Check if any dome lamps are on
            dome_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(9,13) ]
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys()
                                            if k in dome_lamp_stat]), axis=0)
        raise ValueError('No implementation for status = {0}'.format(status))
        
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
        raw_img, head0, _ = read_lris(raw_file, det=det)

        return raw_img, head0

    def get_image_section(self, filename, det, section='datasec'):
        """
        Return a string representation of a slice defining a section of
        the detector image.

        Overwrites base class function to use :func:`read_lris` to get
        the image sections.

        .. todo::
            - It feels really ineffiecient to just get the image section
              using the full :func:`read_lris`.  Can we parse that
              function into something that can give you the image
              section directly?

        This is done separately for the data section and the overscan
        section in case one is defined as a header keyword and the other
        is defined directly.
        
        Args:
            filename (str):
                data filename
            det (int):
                Detector number
            section (:obj:`str`, optional):
                The section to return.  Should be either datasec or
                oscansec, according to the :class:`DetectorPar`
                keywords.

        Returns:
            list, bool: A list of string representations for the image
            sections, one string per amplifier, followed by three
            booleans: if the slices are one indexed, if the slices
            should include the last pixel, and if the slice should have
            their order transposed.
        """
        # Read the file
        temp, head0, secs = read_lris(filename, det)
        if section == 'datasec':
            return secs[0], False, False, False
        elif section == 'oscansec':
            return secs[1], False, False, False
        else:
            raise ValueError('Unrecognized keyword: {0}'.format(section))

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

    def get_match_criteria(self):
        match_criteria = {}
        for key in framematch.FrameTypeBitMask().keys():
            match_criteria[key] = {}
        #
        match_criteria['standard']['match'] = {}
        match_criteria['standard']['match']['dispname'] = ''
        match_criteria['standard']['match']['dichroic'] = ''
        match_criteria['standard']['match']['binning'] = ''
        match_criteria['standard']['match']['decker'] = ''
        # Bias
        match_criteria['bias']['match'] = {}
        match_criteria['bias']['match']['binning'] = ''
        # Pixelflat
        match_criteria['pixelflat']['match'] = match_criteria['standard']['match'].copy()
        # Traceflat
        match_criteria['trace']['match'] = match_criteria['standard']['match'].copy()
        # Arc
        match_criteria['arc']['match'] = match_criteria['standard']['match'].copy()

        # Return
        return match_criteria


class KeckLRISBSpectrograph(KeckLRISSpectrograph):
    """
    Child to handle Keck/LRISb specific code
    """
    def __init__(self):
        # Get it started
        super(KeckLRISBSpectrograph, self).__init__()
        self.spectrograph = 'keck_lris_blue'
        self.camera = 'LRISb'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 1,
                            dispaxis        = 0,
                            dispflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.135,
                            darkcurr        = 0.0,
                            saturation      = 65535.,
                            nonlinear       = 0.86,
                            numamplifiers   = 2,
                            gain            = [1.55, 1.56],
                            ronoise         = [3.9, 4.2],
                            datasec         = ['',''],      # These are provided by read_lris
                            oscansec        = ['',''],
                            suffix          = '_01blue'
                            ),
                #Detector 2
                pypeitpar.DetectorPar(
                            dataext         = 2,
                            dispaxis        = 0,
                            dispflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.135,
                            darkcurr        = 0.,
                            saturation      = 65535.,
                            nonlinear       = 0.86,
                            numamplifiers   = 2,
                            gain            = [1.63, 1.70],
                            ronoise         = [3.6, 3.6],
                            datasec         = ['',''],      # These are provided by read_lris
                            oscansec        = ['',''],
                            suffix          = '_02blue'
                            )]
        self.numhead = 5
        # Uses default timeunit
        # Uses default primary_hdrext
        self.sky_file = 'sky_LRISb_600.fits'

    def default_pypeit_par(self):
        """
        Set default parameters for Keck LRISr reductions.
        """
        par = KeckLRISSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'keck_lris_blue'

        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect'] = 10.0
        #par['calibrations']['wavelengths']['lowest_nsig'] = 10.0
        par['calibrations']['wavelengths']['lamps'] = ['NeI', 'ArI', 'CdI', 'KrI', 'XeI', 'ZnI', 'HgI']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 1


        return par

    def check_headers(self, headers):
        """
        Check headers match expectations for an LRISb exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = { '0.INSTRUME': 'LRISBLUE',
                               '1.NAXIS': 2,
                               '2.NAXIS': 2,
                               '3.NAXIS': 2,
                               '4.NAXIS': 2,
                             '1.CCDGEOM': 'e2v (Marconi) CCD44-82',
                             '1.CCDNAME': '00151-14-1' }
        super(KeckLRISBSpectrograph, self).check_headers(headers, expected_values=expected_values)

    def header_keys(self):
        hdr_keys = super(KeckLRISBSpectrograph, self).header_keys()
        hdr_keys[0]['filter1'] = 'BLUFILT'
        return hdr_keys


class KeckLRISRSpectrograph(KeckLRISSpectrograph):
    """
    Child to handle Keck/LRISr specific code
    """
    def __init__(self):
        # Get it started
        super(KeckLRISRSpectrograph, self).__init__()
        self.spectrograph = 'keck_lris_red'
        self.camera = 'LRISr'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         =1,
                            dispaxis        =0,
                            dispflip        = False,
                            xgap            =0.,
                            ygap            =0.,
                            ysize           =1.,
                            platescale      =0.135,
                            darkcurr        =0.0,
                            saturation      =65535.,
                            nonlinear       =0.76,
                            numamplifiers   =2,
                            gain            =[1.255, 1.18],
                            ronoise         =[4.64, 4.76],
                            datasec         = ['',''],      # These are provided by read_lris
                            oscansec        = ['',''],
                            suffix          ='_01red'
                            ),
                #Detector 2
                pypeitpar.DetectorPar(
                            dataext         =2,
                            dispaxis        =0,
                            dispflip        = False,
                            xgap            =0.,
                            ygap            =0.,
                            ysize           =1.,
                            platescale      =0.135,
                            darkcurr        =0.,
                            saturation      =65535., 
                            nonlinear       =0.76,
                            numamplifiers   =2,
                            gain            =[1.191, 1.162],
                            ronoise         =[4.54, 4.62],
                            datasec         = ['',''],      # These are provided by read_lris
                            oscansec        = ['',''],
                            suffix          ='_02red'
                            )]
        self.numhead = 5
        # Uses default timeunit
        # Uses default primary_hdrext
        # TODO why isn't there a sky file set here?
        # self.sky_file ?

    def default_pypeit_par(self):
        """
        Set default parameters for Keck LRISr reductions.
        """
        par = KeckLRISSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'keck_lris_red'
        #
        par['calibrations']['slits']['sigdetect'] = 50.

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['NeI', 'ArI', 'CdI', 'KrI', 'XeI', 'ZnI', 'HgI']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['sigdetect'] = 10.0
        # reidentification stuff
        #par['calibrations']['wavelengths']['method'] = 'reidentify'
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_red_400_8500_d560.json'


        return par

    def get_lacosmics_par(self,proc_par,binning='1x1'):
        par = self.default_pypeit_par()
        default_sigclip = par['scienceframe']['process']['sigclip']
        default_objlim = par['scienceframe']['process']['objlim']
        # Check whether the user has changed the parameters.
        if (proc_par['sigclip'] == default_sigclip) and (proc_par['objlim'] == default_objlim):
            # Unbinned LRISr needs very aggressive LACosmics parameters.
            if binning is '1x1':
                sigclip = 3.0
                objlim = 0.5
            else:
                sigclip = 5.0
                objlim = 5.0
        else:
            sigclip = proc_par['sigclip']
            objlim = proc_par['objlim']
        return sigclip, objlim

    def check_headers(self, headers):
        """
        Check headers match expectations for an LRISr exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = { '0.INSTRUME': 'LRIS',
                               '1.NAXIS': 2,
                               '2.NAXIS': 2,
                               '3.NAXIS': 2,
                               '4.NAXIS': 2,
                             '1.CCDGEOM': 'LBNL Thick High-Resistivity',
                             '1.CCDNAME': '19-3',
                             '3.CCDNAME': '19-2' }
        super(KeckLRISRSpectrograph, self).check_headers(headers, expected_values=expected_values)

    def header_keys(self):
        hdr_keys = super(KeckLRISRSpectrograph, self).header_keys()
        hdr_keys[0]['filter1'] = 'REDFILT'
        return hdr_keys
            
    def bpm(self, filename=None, det=None, **null_kwargs):
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
        self.empty_bpm(filename=filename, det=det)
        
        # Only defined for det=2
        if det == 2:
            msgs.info("Using hard-coded BPM for det=2 on LRISr")

            # Get the binning
            hdu = fits.open(filename)
            binning = hdu[0].header['BINNING']
            hdu.close()

            # Apply the mask
            xbin = int(binning.split(',')[0])
            badc = 16//xbin
            self.bpm_img[:, 0:badc] = 1

        return self.bpm_img


def read_lris(raw_file, det=None, TRIM=False):
    """
    Read a raw LRIS data frame (one or more detectors)
    Packed in a multi-extension HDU
    Based on readmhdufits.pro

    Parameters
    ----------
    raw_file : str
      Filename
    det : int, optional
      Detector number; Default = both
    TRIM : bool, optional
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
    fil = glob.glob(raw_file+'*') 
    if len(fil) != 1:
        msgs.error("Found {:d} files matching {:s}".format(len(fil)))

    # Read
    msgs.info("Reading LRIS file: {:s}".format(fil[0]))
    hdu = fits.open(fil[0])
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

    # First read over the header info to determine the size of the output array...
    n_ext = len(hdu)-1  # Number of extensions (usually 4)
    xcol = []
    xmax = 0
    ymax = 0
    xmin = 10000
    ymin = 10000
    for i in np.arange(1, n_ext+1):
        theader = hdu[i].header
        detsec = theader['DETSEC']
        if detsec != '0':
            # parse the DETSEC keyword to determine the size of the array.
            x1, x2, y1, y2 = np.array(parse.load_sections(detsec, fmt_iraf=False)).flatten()

            # find the range of detector space occupied by the data
            # [xmin:xmax,ymin:ymax]
            xt = max(x2, x1)
            xmax = max(xt, xmax)
            yt =  max(y2, y1)
            ymax = max(yt, ymax)

            # find the min size of the array
            xt = min(x1, x2)
            xmin = min(xmin, xt)
            yt = min(y1, y2)
            ymin = min(ymin, yt)
            # Save
            xcol.append(xt)

    # determine the output array size...
    nx = xmax - xmin + 1
    ny = ymax - ymin + 1

    # change size for binning...
    nx = nx // xbin
    ny = ny // ybin

    # Update PRECOL and POSTPIX
    precol = precol // xbin
    postpix = postpix // xbin

    # Deal with detectors
    if det in [1,2]:
        nx = nx // 2
        n_ext = n_ext // 2
        det_idx = np.arange(n_ext, dtype=np.int) + (det-1)*n_ext
        ndet = 1
    elif det is None:
        ndet = 2
        det_idx = np.arange(n_ext).astype(int)
    else:
        raise ValueError('Bad value for det')

    # change size for pre/postscan...
    if not TRIM:
        nx += n_ext*(precol+postpix)
        ny += preline + postline

    # allocate output array...
    array = np.zeros( (nx, ny) )
    order = np.argsort(np.array(xcol))

    # insert extensions into master image...
    for kk, i in enumerate(order[det_idx]):

        # grab complete extension...
        data, predata, postdata, x1, y1 = lris_read_amp(hdu, i+1)
                            #, linebias=linebias, nobias=nobias, $
                            #x1=x1, x2=x2, y1=y1, y2=y2, gaindata=gaindata)
        # insert components into output array...
        if not TRIM:
            # insert predata...
            buf = predata.shape
            nxpre = buf[0]
            xs = kk*precol
            xe = xs + nxpre
            '''
            if keyword_set(VERBOSITY) then begin
                section = '['+stringify(xs)+':'+stringify(xe)+',*]'
                message, 'inserting extension '+stringify(i)+ $
                         ' predata  in '+section, /info
            endif 
            '''
            array[xs:xe, :] = predata

            # insert data...
            buf = data.shape
            nxdata = buf[0]
            nydata = buf[1]
            xs = n_ext*precol + kk*nxdata #(x1-xmin)/xbin
            xe = xs + nxdata
            # Data section
            section = '[{:d}:{:d},{:d}:{:d}]'.format(preline,nydata-postline, xs, xe)  # Eliminate lines
            dsec.append(section)
            #print('data',xs,xe)
            array[xs:xe, :] = data   # Include postlines

            #; insert postdata...
            buf = postdata.shape
            nxpost = buf[0]
            xs = nx - n_ext*postpix + kk*postpix
            xe = xs + nxpost 
            section = '[:,{:d}:{:d}]'.format(xs, xe)
            osec.append(section)
            '''
            if keyword_set(VERBOSITY) then begin
                section = '['+stringify(xs)+':'+stringify(xe)+',*]'
                message, 'inserting extension '+stringify(i)+ $
                         ' postdata in '+section, /info
            endif 
            '''
            array[xs:xe, :] = postdata
        else:
            buf = data.shape
            nxdata = buf[0]
            nydata = buf[1]

            xs = (x1-xmin)//xbin
            xe = xs + nxdata 
            ys = (y1-ymin)//ybin
            ye = ys + nydata - postline

            yin1 = preline
            yin2 = nydata - postline 

            '''
            if keyword_set(VERBOSITY) then begin
                section = '['+stringify(xs)+':'+stringify(xe)+ $
                          ','+stringify(ys)+':'+stringify(ye)+']'
                message, 'inserting extension '+stringify(i)+ $
                         ' data     in '+section, /info
            endif 
            '''
            array[xs:xe, ys:ye] = data[:, yin1:yin2]

    # make sure BZERO is a valid integer for IRAF
    obzero = head0['BZERO']
    head0['O_BZERO'] = obzero
    head0['BZERO'] = 32768-obzero

    # Return, transposing array back to goofy Python indexing
    return array.T, head0, (dsec, osec)


def lris_read_amp(inp, ext):
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
    if isinstance(inp, str):
        hdu = fits.open(inp)
    else:
        hdu = inp

    # Get the pre and post pix values
    # for LRIS red POSTLINE = 20, POSTPIX = 80, PRELINE = 0, PRECOL = 12
    head0 = hdu[0].header
    precol = head0['precol']
    postpix = head0['postpix']

    # Deal with binning
    binning = head0['BINNING']
    xbin, ybin = [int(ibin) for ibin in binning.split(',')]
    precol = precol//xbin
    postpix = postpix//xbin

    # get entire extension...
    temp = hdu[ext].data.transpose() # Silly Python nrow,ncol formatting
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
    predata = temp[0:precol, :]
    # datasec appears to have the x value for the keywords that are zero
    # based. This is only true in the image header extensions
    # not true in the main header.  They also appear inconsistent between
    # LRISr and LRISb!
    #data     = temp[xdata1-1:xdata2-1,*]
    #data = temp[xdata1:xdata2+1, :]
    if (xdata1-1) != precol:
        msgs.error("Something wrong in LRIS datasec or precol")
    xshape = 1024 // xbin
    if (xshape+precol+postpix) != temp.shape[0]:
        msgs.error("Wrong size for in LRIS detector somewhere.  Funny binning?")
    data = temp[precol:precol+xshape,:]
    postdata = temp[nxt-postpix:nxt, :]

    # flip in X as needed...
    if x1 > x2:
        xt = x2
        x2 = x1
        x1 = xt
        data = np.flipud(data) #reverse(temporary(data),1)

    # flip in Y as needed...
    if y1 > y2:
        yt = y2
        y2 = y1
        y1 = yt
        data = np.fliplr(data)
        predata = np.fliplr(predata)
        postdata = np.fliplr(postdata)

    '''
    #; correct gain if requested...
    if keyword_set(GAINDATA) then begin
        gain = gainvalue( gaindata, header)
        data = FLOAT(temporary(data)) * gain
        predata = FLOAT(temporary(predata)) * gain
        postdata = FLOAT(temporary(postdata)) * gain
    endif
    '''

    '''
    ;; optional bias subtraction...
    if ~ keyword_set(NOBIAS) then begin
        if keyword_set( LINEBIAS) then begin
            ;; compute a bias for each line...
            bias = median( postdata, dim=1)

            ;; subtract for data...
            buf = size(data)
            nx = buf[1]
            ny = buf[2]
            data2 = fltarr(nx,ny)
            for i=0,nx-1 do begin
                data2[i,*] = float(data[i,*]) - bias
            endfor 
            data = data2
        endif else begin
            ;; compute a scalar bias....
            bias = median( postdata)
            data -= bias
        endelse
    endif
    '''

    return data, predata, postdata, x1, y1


'''
def bpm(slf, camera, fitsdict, det):
    """  Wrapper for core_bpm
    Will likely be deprecated

    Parameters
    ----------
    slf
    camera
    fitsdict
    det

    Returns
    -------
    badpix : ndarray

    """
    sidx = slf._idx_sci[0]
    # Binning
    xbin, ybin = [int(ii) for ii in fitsdict['binning'][sidx].split(',')]
    return core_bpm(xbin, ybin, camera, det)
'''



def convert_lowredux_pixelflat(infil, outfil):
    """ Convert LowRedux pixelflat to PYPIT format
    Returns
    -------

    """
    # Read
    hdu = fits.open(infil)
    data = hdu[0].data

    #
    prihdu = fits.PrimaryHDU()
    hdus = [prihdu]
    prihdu.header['FRAMETYP'] = 'pixelflat'

    # Detector 1
    img1 = data[:,:data.shape[1]//2]
    hdu = fits.ImageHDU(img1)
    hdu.name = 'DET1'
    prihdu.header['EXT0001'] = 'DET1-pixelflat'
    hdus.append(hdu)

    # Detector 2
    img2 = data[:,data.shape[1]//2:]
    hdu = fits.ImageHDU(img2)
    hdu.name = 'DET2'
    prihdu.header['EXT0002'] = 'DET2-pixelflat'
    hdus.append(hdu)

    # Finish
    hdulist = fits.HDUList(hdus)
    hdulist.writeto(outfil, clobber=True)
    print('Wrote {:s}'.format(outfil))

