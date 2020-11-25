"""
Implements HIRES-specific functions, including reading in slitmask design
files.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import glob

from IPython import embed

import numpy as np

from scipy import interpolate

from pypeit import msgs
from pypeit import telescopes
from pypeit import io
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph


class KECKHIRESSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle KECK/HIRES specific code.

    This spectrograph is not yet supported.
    """
    ndet = 1
    telescope = telescopes.KeckTelescopePar()
    pypeline = 'Echelle'

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()
        # Correct for flexure using the default approach
        par['flexure'] = pypeitpar.FlexurePar()
        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(ext=0, card='DECKNAME')
        self.meta['binning'] = dict(ext=0, card='BINNING')

        self.meta['mjd'] = dict(ext=0, card='MJD')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        self.meta['dispname'] = dict(ext=0, card='ECHNAME')
        # Extras for config and frametyping
#        self.meta['echangl'] = dict(ext=0, card='ECHANGL')
#        self.meta['xdangl'] = dict(ext=0, card='XDANGL')

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.

        Args:
            ftype (:obj:`str`):
                Type of frame to check. Must be a valid frame type; see
                frame-type :ref:`frame_type_defs`.
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for one or more frames to check.
            exprng (:obj:`list`, optional):
                Range in the allowed exposure time for a frame of type
                ``ftype``. See
                :func:`pypeit.core.framematch.check_frame_exptime`.

        Returns:
            `numpy.ndarray`_: Boolean array with the flags selecting the
            exposures in ``fitstbl`` that are ``ftype`` type frames.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        # TODO: Allow for 'sky' frame type, for now include sky in
        # 'science' category
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'Object')
        if ftype == 'standard':
            return good_exp & ((fitstbl['idname'] == 'Std') | (fitstbl['idname'] == 'Object'))
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'Bias')
        if ftype == 'dark':
            return good_exp & (fitstbl['idname'] == 'Dark')
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & ((fitstbl['idname'] == 'Flat') | (fitstbl['idname'] == 'IntFlat'))
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'Line')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

# TODO: This spectrograph is way out of date
#    def load_raw_img_head(self, raw_file, det=None, **null_kwargs):
#        """
#        Wrapper to the raw image reader for HIRES
#
#        Args:
#            raw_file (:obj:`str`):
#                filename
#            det (:obj:`int`, optional):
#              Desired detector.  Despite default value, cannot be
#              ``None`` (todo: set a sensible default).
#            **null_kwargs:
#              Captured and never used
#
#        Returns:
#            tuple: Raw image and header
#
#        """
#        raw_img, head0, _ = read_hires(raw_file, det=det)
#
#        return raw_img, head0
#
#    def get_image_section(self, inp=None, det=1, section='datasec'):
#        """
#        Return a string representation of a slice defining a section of
#        the detector image.
#
#        Overwrites base class function to use :func:`read_hires` to get
#        the image sections.
#
#        .. todo ::
#            - It is really ineffiecient.  Can we parse
#              :func:`read_hires` into something that can give you the
#              image section directly?
#
#        This is done separately for the data section and the overscan
#        section in case one is defined as a header keyword and the other
#        is defined directly.
#
#        Args:
#            inp (:obj:`str`, `astropy.io.fits.Header`_, optional):
#                String providing the file name to read, or the relevant
#                header object.  Default is None, meaning that the
#                detector attribute must provide the image section
#                itself, not the header keyword.
#            det (:obj:`int`, optional):
#                1-indexed detector number.
#            section (:obj:`str`, optional):
#                The section to return.  Should be either 'datasec' or
#                'oscansec', according to the
#                :class:`pypeitpar.DetectorPar` keywords.
#
#        Returns:
#            tuple: Returns three objects: (1) A list of string
#            representations for the image sections, one string per
#            amplifier.  The sections are *always* returned in PypeIt
#            order: spectral then spatial.  (2) Boolean indicating if the
#            slices are one indexed.  (3) Boolean indicating if the
#            slices should include the last pixel.  The latter two are
#            always returned as True following the FITS convention.
#        """
#        # Read the file
#        if inp is None:
#            msgs.error('Must provide Keck HIRES file to get image section.')
#        elif not os.path.isfile(inp):
#            msgs.error('File {0} does not exist!'.format(inp))
#        temp, head0, secs = read_hires(inp, det)
#        if section == 'datasec':
#            return secs[0], False, False
#        elif section == 'oscansec':
#            return secs[1], False, False
#        else:
#            raise ValueError('Unrecognized keyword: {0}'.format(section))
#
#     def get_datasec_img(self, filename, det=1, force=True):
#         """
#         Create an image identifying the amplifier used to read each pixel.
#
#         Args:
#             filename (str):
#                 Name of the file from which to read the image size.
#             det (:obj:`int`, optional):
#                 Detector number (1-indexed)
#             force (:obj:`bool`, optional):
#                 Force the image to be remade
#
#         Returns:
#             `numpy.ndarray`: Integer array identifying the amplifier
#             used to read each pixel.
#         """
#         if self.datasec_img is None or force:
#             # Check the detector is defined
#             self._check_detector()
#             # Get the image shape
#             raw_naxis = self.get_raw_image_shape(filename, det=det)
#
#             # Binning is not required because read_hires accounts for it
# #            binning = self.get_meta_value(filename, 'binning')
#
#             data_sections, one_indexed, include_end, transpose \
#                     = self.get_image_section(filename, det, section='datasec')
#
#             # Initialize the image (0 means no amplifier)
#             self.datasec_img = np.zeros(raw_naxis, dtype=int)
#             for i in range(self.detector[det-1]['numamplifiers']):
#                 # Convert the data section from a string to a slice
#                 datasec = parse.sec2slice(data_sections[i], one_indexed=one_indexed,
#                                           include_end=include_end, require_dim=2,
#                                           transpose=transpose) #, binning=binning)
#                 # Assign the amplifier
#                 self.datasec_img[datasec] = i+1
#         return self.datasec_img


class KECKHIRESRSpectrograph(KECKHIRESSpectrograph):
    """
    Child to handle KECK/HIRES-R specific code

    .. warning::
        Spectrograph not yet supported
    """
    name = 'keck_hires_red'
    camera = 'HIRES_R'

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        # Adjustments to slit and tilts for NIR
        par['calibrations']['slitedges']['edge_thresh'] = 600.
        par['calibrations']['slitedges']['fit_order'] = 5
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['left_right_pca'] = True

        par['calibrations']['tilts']['tracethresh'] = 20
        # Bias
        par['calibrations']['biasframe']['useframe'] = 'bias'

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['ThAr']
        #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['rms_threshold'] = 0.25
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        # Reidentification parameters
        #par['calibrations']['wavelengths']['method'] = 'reidentify'
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_xshooter_nir.json'
        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] ='nothing'
        par['calibrations']['standardframe']['exprng'] = [None, 600]
        par['scienceframe']['exprng'] = [600, None]

        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        super().init_meta()
        self.meta['decker'] = dict(ext=0, card='DECKNAME')


def indexing(itt, postpix, det=None,xbin=None,ybin=None):
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
    ii = int(np.round(2048/xbin))
    jj = int(np.round(4096/ybin))
    # y indices
    y1, y2 = 0, jj
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

def hires_read_1chip(hdu,chipno):
    """ Read one of the HIRES detectors

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

    x1_dat, x2_dat, y1_dat, y2_dat = np.array(parse.load_sections(datsec)).flatten()
    x1_det, x2_det, y1_det, y2_det = np.array(parse.load_sections(detsec)).flatten()

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

def read_hires(raw_file, det=None):
    """
    Read a raw HIRES data frame (one or more detectors).

    Data are unpacked from the multi-extension HDU.  Function is
    based :func:`pypeit.spectrographs.keck_lris.read_lris`, which
    was based on the IDL procedure ``readmhdufits.pro``.
    
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
        msgs.error('Found {0} files matching {1}'.format(len(fil), raw_file + '*'))
    # Read
    try:
        msgs.info("Reading HIRES file: {:s}".format(fil[0]))
    except AttributeError:
        print("Reading HIRES file: {:s}".format(fil[0]))

    hdu = io.fits_open(fil[0])
    head0 = hdu[0].header

    # Get post, pre-pix values
    precol = head0['PRECOL']
    postpix = head0['POSTPIX']
    preline = head0['PRELINE']
    postline = head0['POSTLINE']
    detlsize = head0['DETLSIZE']
    x0, x_npix, y0, y_npix = np.array(parse.load_sections(detlsize)).flatten()

    # Create final image
    if det is None:
        image = np.zeros((x_npix,y_npix+4*postpix))

    # Setup for datasec, oscansec
    dsec = []
    osec = []

    # get the x and y binning factors...
    binning = head0['BINNING']
    if binning != '3,1':
        msgs.warn("This binning for HIRES might not work.  But it might..")

    xbin, ybin = [int(ibin) for ibin in binning.split(',')]

    # HIRES detectors
    nchip = 3


    if det is None:
        chips = range(nchip)
    else:
        chips = [det-1] # Indexing starts at 0 here
    # Loop
    for tt in chips:
        data, oscan = hires_read_1chip(hdu, tt+1)

        # One detector??
        if det is not None:
            image = np.zeros((data.shape[0],data.shape[1]+oscan.shape[1]))

        # Indexing
        x1, x2, y1, y2, o_x1, o_x2, o_y1, o_y2 = indexing(tt, postpix, det=det,xbin=xbin,ybin=ybin)

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


