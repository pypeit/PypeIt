""" Module for LRIS specific codes
"""
from __future__ import absolute_import, division, print_function

import glob

import numpy as np
from astropy.io import fits

from pypeit import msgs
from pypeit.core import parse
from pypeit.par.pypeitpar import DetectorPar
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit import telescopes
from pypeit.core import fsort

from pypeit import debugger


class KeckLRISSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/LRIS specific code
    """
    def __init__(self):
        # Get it started
        super(KeckLRISSpectrograph, self).__init__()
        self.spectrograph = 'keck_lris'
        self.telescope = telescopes.KeckTelescopePar()

    def lris_header_keys(self):
#        def_keys = self.default_header_keys()

        hdr_keys = {}
        hdr_keys[0] = {}
        hdr_keys[1] = {}
        hdr_keys[2] = {}
        hdr_keys[3] = {}
        hdr_keys[4] = {}

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

        hdr_keys[0]['lampstat01'] = 'MERCURY'
        hdr_keys[0]['lampstat02'] = 'NEON'
        hdr_keys[0]['lampstat03'] = 'ARGON'
        hdr_keys[0]['lampstat04'] = 'CADMIUM'
        hdr_keys[0]['lampstat05'] = 'ZINC'
        hdr_keys[0]['lampstat06'] = 'KRYPTON'
        hdr_keys[0]['lampstat07'] = 'XENON'
        hdr_keys[0]['lampstat08'] = 'FEARGON'
        hdr_keys[0]['lampstat09'] = 'DEUTERIUM'
        hdr_keys[0]['lampstat10'] = 'FLAMP1'
        hdr_keys[0]['lampstat11'] = 'FLAMP2'

        return hdr_keys

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
        for key in fsort.ftype_list:
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
                DetectorPar(dataext         = 1,
                            dispaxis        = 0,
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
                DetectorPar(dataext         = 2,
                            dispaxis        = 0,
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
                            )
            ]
        self.numhead = 5
        # Uses default timeunit
        # Uses default primary_hdrext
        self.sky_file = 'sky_LRISb_600.fits'

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Keck LRISb reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'keck_lris_blue'
        # Use the ARMS pipeline
        par['rdx']['pipeline'] = 'ARMS'
        # Set wave tilts order
        par['calibrations']['slits']['sigdetect'] = 30.
        par['calibrations']['slits']['pcapar'] = [3,2,1,0]
        # Always sky subtract, starting with default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()
        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        return par

    def check_header(self, headers):
        """Validate elements of the header."""
        chk_dict = {}
        # chk_dict is 1-indexed!
        chk_dict[2] = {}
        # THIS CHECK IS A MUST! It performs a standard check to make sure the data are 2D.
        chk_dict[2]['NAXIS'] = 2
        # Check the CCD name
        chk_dict[2]['CCDGEOM'] = 'e2v (Marconi) CCD44-82'
        chk_dict[2]['CCDNAME'] = '00151-14-1'
        return chk_dict

    def header_keys(self):
        head_keys = self.lris_header_keys()
        head_keys[0]['filter1'] = 'BLUFILT'
        return head_keys

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
        arcparam['lamps'] = ['NeI', 'ArI', 'CdI', 'KrI', 'XeI', 'ZnI', 'HgI']
        if disperser == '600/4000':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=0.63 # Ang per pixel (unbinned)
            arcparam['b1']= 4.54698031e-04
            arcparam['b2']= -6.86414978e-09
            arcparam['wvmnx'][1] = 6000.
            arcparam['wv_cen'] = 4000.
            arcparam['min_ampl'] = 1000.0
        elif disperser == '400/3400':
            pass
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=1.02
            arcparam['b1']= 2.72694493e-04
            arcparam['b2']= -5.30717321e-09
            arcparam['wvmnx'][1] = 6000.
            arcparam['min_ampl'] = 1000.0
        elif disperser == '300/5000':
            arcparam['n_first'] = 2
            arcparam['wv_cen'] = 4500.
            arcparam['disp'] = 1.43
            arcparam['min_ampl'] = 1000.0
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))

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
                DetectorPar(dataext         =1,
                            dispaxis        =0,
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
                DetectorPar(dataext         =2,
                            dispaxis        =0,
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
                            )
            ]
        self.numhead = 5
        # Uses default timeunit
        # Uses default primary_hdrext
        # TODO why isn't there a sky file set here?
        # self.sky_file ?

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Keck LRISr reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'keck_lris_red'
        # Use the ARMS pipeline
        par['rdx']['pipeline'] = 'ARMS'
        # Set wave tilts order
        par['calibrations']['slits']['sigdetect'] = 30.
        par['calibrations']['slits']['pcapar'] = [3,2,1,0]
        # Always sky subtract, starting with default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()
        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        return par


    def header_keys(self):
        head_keys = self.lris_header_keys()
#        head_keys[0]['filter1'] = 'BLUFILT'
        return head_keys

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

    def setup_arcparam(self, arcparam, disperser=None, fitstbl=None, arc_idx=None,
                       msarc_shape=None, binspectral=None, **null_kwargs):
        """
        Setup the arc parameters

        Args:
            arcparam: dict
            disperser: str, REQUIRED

        Returns:
            arcparam is modified in place

        """
        #arcparam['wv_cen'] = fitstbl['wavecen'][arc_idx]
        # Should set according to the lamps that were on
        arcparam['lamps'] = ['ArI','NeI','HgI','KrI','XeI']
        '''
        if disperser == '600/7500':
            arcparam['n_first']=3 # Too much curvature for 1st order
            arcparam['disp']=0.80 # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][1] = 11000.
        elif disperser == '600/10000':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=0.80 # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][1] = 12000.
        elif disperser == '400/8500':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=1.19 # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][1] = 11000.
            arcparam['min_ampl'] = 3000.  # Lines tend to be very strong
            arcparam['nsig_rej_final'] = 5.
        elif disperser == '900/5500':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=0.53 # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][1] = 7000.
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))
        '''

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
