""" Module to define the Spectrograph class
"""
from __future__ import absolute_import, division, print_function

import os

from abc import ABCMeta
from pkg_resources import resource_filename

import numpy as np
from astropy.io import fits

from linetools.spectra import xspectrum1d

from .. import arparse
from .. import arpixels
from ..par.pypitpar import DetectorPar, TelescopePar
from pypit import msgs

try:
    basestring
except NameError:
    basestring = str

# TODO: Consider changing the name of this to Instrument
class Spectrograph(object):
    """
    Generic class for spectrograph-specific codes

    Attributes:
        spectrograph (str):
            The name of the spectrograph.  See
            :func:`pypit.spectrographs.util.valid_spectrographs` for the
            currently supported spectrographs.
        telescope (:class:`TelescopePar`):
            Parameters of the telescope that feeds this spectrograph.
        camera (str):
            Name of the spectrograph camera.
        detector (list):
            A list of instances of
            :class:`pypit.par.pypitpar.DetectorPar` with the parameters
            for each detector in the spectrograph
        naxis (tuple):
            A tuple with the lengths of the two axes for current
            detector image.
        datasec_img (:obj:`numpy.ndarray`):
            An image identifying the amplifier that reads each detector
            pixel.
        bpm_img (:obj:`numpy.ndarray`):
            The bad-pixel mask for the currently read detector.
        sky_file (str):
            A file with an archived sky spectrum, used for flexure
            corrections.
    """
    __metaclass__ = ABCMeta

    def __init__(self):
        self.spectrograph = 'generic'
        self.telescope = None
        self.camera = None
        self.detector = None
        self.naxis = None
        self.datasec_img = None
        self.bpm_img = None

        # Default time unit
        self.timeunit = 'mjd'

        # Default extension with the primary header data
        #   used by arsave.save_2d_images
        self.primary_hdrext = 0
        self.numhead = 0

        self.sky_file = None

    @staticmethod
    def default_sky_spectrum():
        """
        Return the path to the default sky spectrum: currently
        'pypit/data/sky_spec/paranal_sky.fits' in the pypit source
        distribution.
        """
        return os.path.join(resource_filename('pypit', 'data/sky_spec/'), 'paranal_sky.fits')

    def add_to_fitstbl(self, fitstbl):
        pass

    def _check_telescope(self):
        # Check the detector
        if self.telescope is None:
            raise ValueError('Must define the telescope used to take the observations.')
        if not isinstance(self.telescope, TelescopePar):
                raise TypeError('Telescope parameters must be one of those specified in'
                                'pypit.telescopes.')

    def _check_detector(self):
        # Check the detector
        if self.detector is None:
            raise ValueError('Must first define spectrograph detector parameters!')
        for d in self.detector:
            if not isinstance(d, DetectorPar):
                raise TypeError('Detector parameters must be specified using DetectorPar.')

    def load_raw_frame(self, raw_file, det=None):
        """
        Load the image (converted to np.float) and primary header of the input file

        The image is transposed, as needed, so that the spectral dimension
        runs along the columns

        Args:
            raw_file:  str, filename
            dataext: int, optional
              Extension in the FITS list for the data
            det: int, optional
              Desired detector

        Returns:
            img: ndarray
              Converted to np.float and transposed if necessary
            head0: Header

        """
        self._check_detector()
        _det = 1 if det is None else det

        # Load the raw image
        raw_img, head0 = self.load_raw_img_head(raw_file, dataext=self.detector[_det-1]['dataext'],
                                                det=_det)

        # Turn to float
        img = raw_img.astype(np.float)
        # Transpose?
        if self.detector[_det-1]['dispaxis'] == 1:
            img = img.T
        # Return
        return img, head0

    def load_raw_img_head(self, raw_file, dataext, **null_kwargs):
        """
        Generic raw image reader

        Args:
            raw_file: str
            dataext: int
            **null_kwargs:
              Captured and never used

        Returns:
            raw_img: ndarray
              Raw image;  likely unsigned int
            head0: Header

        """
        # Open and go
        hdulist = fits.open(raw_file)
        raw_img = hdulist[dataext].data
        head0 = hdulist[0].header
        # Return
        return raw_img, head0

    def get_image_section(self, filename, det, section='datasec'):
        """
        Return a string representation of a slice defining a section of
        the detector image.

        This default function of the base class tries to get the image
        section in two ways, first by assuming the image section
        defined by the detector is a header keyword, and then by just
        assuming the detector provides the image section directly.

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
        # Check the section is one of the detector keywords
        # TODO: Is this too broad?  I didn't want to specify 'datasec'
        # and 'oscansec' in case those change.
        if section not in self.detector[det-1].keys():
            raise ValueError('Unrecognized keyword: {0}'.format(section))

        # Check the detector is defined
        self._check_detector()

        # Get the data section
        try:
            # Try using the image sections as header keywords
            hdu = fits.open(filename)
            image_sections = [ hdu[self.detector[det-1]['dataext']].header[key] \
                                    for key in self.detector[det-1][section] ]
            # If this is successful, assume normal FITS header formatting
            one_indexed = True
            include_last = True
            transpose = True
        except:
            # Otherwise use the detector definition directly
            image_sections = self.detector[det-1][section]
            # and assume the string is python formatted.
            one_indexed = False
            include_last = False
            transpose = False

        return image_sections, one_indexed, include_last, transpose

    def get_datasec_img(self, filename, det=1, force=True):
        """
        Create an image identifying the amplifier used to read each
        pixel.

        .. todo::
            - I find 1-indexing to be highly annoying...
            - Check for overlapping amplifiers?
            - Consider renaming this datasec_ampid or something like
              that.  I.e., the image's main purpose is to tell you where
              the amplifiers are for the data section
          
        Args:
            filename (str):
                Name of the file from which to read the image size.
                TODO: Can this be obtained from the fitstbl?
            det (int):
                Detector number (1-indexed)
            force (:obj:`bool`, optional):
                Force the image to be remade

        Returns:
            `numpy.ndarray`: Integer array identifying the amplifier
            used to read each pixel.
        """
        if self.datasec_img is None or force:
            # Check the detector is defined
            self._check_detector()
            # Get the image shape
            self.get_image_shape(filename=filename, det=det)

            data_sections, one_indexed, include_end, transpose \
                    = self.get_image_section(filename, det, section='datasec')

            # Initialize the image (0 means no amplifier)
            self.datasec_img = np.zeros(self.naxis, dtype=int)
            for i in range(self.detector[det-1]['numamplifiers']):
                # Convert the data section from a string to a slice
                # TODO: Should consider doing this in DetectorPar
                datasec = arparse.sec2slice(data_sections[i], one_indexed=one_indexed,
                                            include_end=include_end, require_dim=2,
                                            transpose=transpose)
                # Assign the amplifier
                self.datasec_img[datasec] = i+1
#            self.datasec_img = arpixels.pix_to_amp(self.naxis[0], self.naxis[1], 
#                                                   self.detector[det-1]['datasec'],
#                                                   self.detector[det-1]['numamplifiers'])
        return self.datasec_img

    def get_image_shape(self, filename=None, det=None, force=True):
        """
        Get the shape of the image data for a given detector using a
        file.  :attr:`detector` must be defined.

        Fails if filename is None and the instance does not have a
        predefined :attr:`naxis`.  If the filename is None, always
        returns the predefined :attr:`naxis`.  If the filename is
        provided, the header of the associated detector is used.  If the
        the detector is set to None, the primary header of the file is
        used.
        
        Args:
            filename (:obj:`str`, optional):
                Name of the fits file with the header to use.
            det (:obj:`int`, optional):
                1-indexed number of the detector.  Default is None.  If
                None, the primary extension is used.  Otherwise the
                internal detector parameters are used to determine the
                extension to read.
            force (:obj:`bool`, optional):
                Force the image shape to be redetermined.
            null_kwargs (dict):
                Used to catch any extraneous keyword arguments.
        
        Returns:
            tuple: Tuple of two integers with the length of each image
            axes.

        Raises:
            ValueError:
                Raised if the image shape cannot be determined from the
                input and available attributes.
        """
        # Cannot be determined
        if (self.naxis is None or force) and filename is None:
            raise ValueError('Cannot determine image shape!  Must have NAXIS predefined or '
                             'provide a file to read.')

        # Return the predefined value
        if filename is None:
            return self.naxis

        # Use a file
        self._check_detector()
        self.naxis = (self.load_raw_frame(filename, det=det)[0]).shape
        return self.naxis

    def empty_bpm(self, shape=None, filename=None, det=1, force=True):
        """
        Generate a generic (empty) BPM.
        
        This requires a successful call to :func:`get_image_shape`.

        .. todo::
            Any reason this isn't returned as a boolean array?

        Args:
            shape: tuple, REQUIRED
            **null_kwargs:

        Returns:
            bpm: ndarray, int
              0=not masked; 1=masked

        """
        if self.bpm_img is None or force:
            if shape is None:
                # Check the detector is defined
                self._check_detector()
                # Get the image shape
                _shape = self.get_image_shape(filename=filename, det=det)
            else:
                _shape = shape
            self.bpm_img = np.zeros(_shape, dtype=np.int8)
        return self.bpm_img

    def bpm(self, shape=None, filename=None, det=1, force=True):
        """
        Generate a default bad-pixel mask.

        Currently identical to calling :func:`empty_bpm`.

        Args:
            shape: tuple, REQUIRED
            **null_kwargs:

        Returns:
            bpm: ndarray, int
              0=not masked; 1=masked

        """
        return self.empty_bpm(shape=shape, filename=filename, det=det, force=force)

    def default_header_keys(self):
        def_head_keys = {}
        def_head_keys[0] = {}
        def_head_keys[0]['target'] = 'OBJECT'     # Header keyword for the name given by the observer to a given frame
        def_head_keys[0]['idname'] = 'OBSTYPE'    # The keyword that identifies the frame type (i.e. bias, flat, etc.)
        def_head_keys[0]['time'] = 'MJD-OBS'      # The time stamp of the observation (i.e. decimal MJD)
        def_head_keys[0]['date'] = 'DATE'         # The UT date of the observation which is used for heliocentric (in the format YYYY-MM-DD  or  YYYY-MM-DDTHH:MM:SS.SS)
        def_head_keys[0]['ra'] = 'RA'             # Right Ascension of the target
        def_head_keys[0]['dec'] = 'DEC'           # Declination of the target
        def_head_keys[0]['airmass'] = 'AIRMASS'   # Airmass at start of observation
        def_head_keys[0]['binning'] = 'BINNING'   # Binning
        def_head_keys[0]['exptime'] = 'EXPTIME'   # Exposure time keyword
        def_head_keys[0]['decker'] = 'SLITNAME'
        def_head_keys[0]['dichroic'] = 'DICHNAME' # Dichroic name
        def_head_keys[0]['dispname'] = 'GRISNAME' # Grism name
        # Return
        return def_head_keys

    def header_keys(self):
        return self.default_header_keys()

    def get_headarr(self, filename, reduce_par):
        headarr = ['None' for k in range(self.numhead)]
        # Try to load em up
        try:
            for k in range(self.numhead):
                headarr[k] = fits.getheader(filename, ext=k)
        except:
            if reduce_par['setup']:
                msgs.warn("Bad header in extension {0:d} of file:".format(filename))
                msgs.warn("Proceeding on the hopes this was a calibration file, otherwise consider removing.")
            else:
                msgs.error("Error reading header from extension {0:d} of file:".format(filename))
        return headarr

    def check_headers(self, headers):
        pass

    def setup_arcparam(self, **null_kwargs):
        return None

    @property
    def ndet(self):
        """Return the number of detectors."""
        return 0 if self.detector is None else len(self.detector)

    def archive_sky_spectrum(self):
        """
        Load an archived sky spectrum based on :attr:`sky_file`.
        
        Returns:
            str, :class:`linetools.xspectrum1d.XSpectrum1D`: The name of
            the file and the instance of :class:`XSpectrum1D` with the
            spectrum data.

        Raises:
            FileNotFoundError:
                Raised if the file does not exist as written or in the
                pypit/data/sky_spec/ directory in the source
                distribution.
        """
        # No file was defined
        if self.sky_file is None:
            self.sky_file = Spectrograph.default_sky_spectrum()
            msgs.warn('Using default sky spectrum: {0}'.format(self.sky_file))

        if os.path.isfile(self.sky_file):
            # Found directly
            return self.sky_file, xspectrum1d.XSpectrum1D.from_file(self.sky_file)

        root = resource_filename('pypit', 'data/sky_spec/')
        _sky_file = os.path.join(root, self.sky_file)
        if os.path.isfile(_sky_file):
            # Found within the source distribution
            return self.sky_file, xspectrum1d.XSpectrum1D.from_file(_sky_file)

        # File could not be read
        raise FileNotFoundError('Could not find archive sky spectrum: {0} or {1}'.format(
                                    self.sky_file, _sky_file))


    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: '.format(self.__class__.__name__)
        txt += ' spectrograph={:s},'.format(self.spectrograph)
        txt += ' camera={:s}'.format(self.camera)
        txt += '>'
        return txt
