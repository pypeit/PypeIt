""" Module to define the Spectrograph class
"""
from __future__ import absolute_import, division, print_function

from abc import ABCMeta
from pkg_resources import resource_filename

import numpy as np
from astropy.io import fits

from linetools.spectra import xspectrum1d

from .. import arpixels
from ..par.pypitpar import DetectorPar, TelescopePar

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

        # Default sky spectrum
        self.sky_file = os.path.join(resource_filename('pypit', 'data/sky_spec/'),
                                     'paranal_sky.fits')

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

    def get_datasec(self, filename, det):
        """
        Return a string representation of the image slice defining the
        data section of the detector.
        
        This default function of the base class assumes that the
        detector attribute provides the header keyword with the data
        section in the extension appropriate to each detector.

        Args:
            filename (str):
                data filename
            det (int):
                Detector number

        Returns:
            list, bool: A list of the strings defining the data section
            for each amplifier and a flag that the data sections should
            be transposed to appropriately match python numpy arrays.

        """
        self._check_detector()
        hdu = fits.open(filename)
        return [ hdu[self.detector[det-1]['dataext']].header[key] \
                                for key in self.detector[det-1]['datasec'] ], True, True, True

    @staticmethod
    def parse_sec2slice(subarray, one_indexed=False, include_end=False, require_dim=None,
                        transpose=False):
        """
        Convert a string representation of an array subsection (slice)
        into a list of slice objects.

        .. todo::
            This is really a general utility that should go somewhere
            else.

        Args:
            subarray (str):
                The string to convert.  Should have the form of normal
                slice operation, e.g.::
                    
                    print(x[:10])
                    print(x[10:])
                    print(x[::-1])
                    print(x[:2,:,:])
                
                The string ignores whether or not the string has the
                brackets '[]', but the string must contain the
                appropriate ':' and ',' characters.
            one_indexed (:obj:`bool`, optional):
                The string should be interpreted as 1-indexed.  Default
                is to assume python indexing.
            include_end (:obj:`bool`, optional):
                **If** the end is defined, adjust the slice such that
                the last element is included.  Default is to exclude the
                last element as with normal python slicing.
            require_dim (:obj:`int`, optional):
                Test if the string indicates the slice along the proper
                number of dimensions.
            transpose (:obj:`bool`, optional):
                Transpose the order of the returned slices.  The
                following are equivalent::
                    
                    tslices = parse_sec2slice('[:10,10:]')[::-1]
                    tslices = parse_sec2slice('[:10,10:]', transpose=True)

        Returns:
            list: A list of slice objects, one per dimension of the
            prospective array.

        Raises:
            TypeError:
                Raised if the input `subarray` is not a string.
            ValueError:
                Raised if the string does not match the required
                dimensionality or if the string does not look like a
                slice.
        """
        # Check it's a string
        if not isinstance(subarray, basestring):
            raise TypeError('Can only parse string-based subarray sections.')
        # Remove brackets if they're included
        sections = subarray.strip('[]').split(',')
        # Check the dimensionality
        ndim = len(sections)
        if require_dim is not None and ndim != require_dim:
            raise ValueError('Number of slices ({0}) does not match '.format(ndim) + 
                             'required dimensions ({0}).'.format(require_dim))
        # Convert the slice of each dimension from a string to a slice
        # object
        slices = []
        for s in sections:
            # Must be able to find the colon
            if ':' not in s:
                raise ValueError('Unrecognized slice string: {0}'.format(s))
            # Initial conversion
            s = [ None if x == '' else int(x) for x in s.split(':') ]
            if len(s) < 3:
                # Include step
                s += [ None ]
            if one_indexed:
                # Decrement to convert from 1- to 0-indexing
                s = [ None if x is None else x-1 for x in s ]
            if include_end and s[1] is not None:
                # Increment to include last 
                s[1] += 1
            # Append the new slice
            slices += [slice(*s)]

        return slices[::-1] if transpose else slices

    def get_datasec_img(self, filename, det=1, force=True):
        """
        Create an image identifying the amplifier used to read each
        pixel.

        .. todo::
            - I find 1-indexing to be highly annoying...
            - Check for overlapping amplifiers?

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

            # Get the data sections
            one_indexed = False
            include_end = False
            transpose = False
            try:
                # First try the default approach
                data_sections, one_indexed, include_end, transpose \
                            = self.get_datasec(filename, det)
            except:
                # If that doesn't work, assume the detector definitions
                # are correct and appropriately defined
                # (Danger, Will Robinson! Danger!)
                data_sections = self.detector[det-1]['datasec']

            # Initialize the image (0 means no amplifier)
            self.datasec_img = np.zeros(self.naxis, dtype=int)
            for i in range(self.detector[det-1]['numamplifiers']):
                # Convert the data section from a string to a slice
                # TODO: Should consider doing this in DetectorPar
                datasec = Spectrograph.parse_sec2slice(data_sections[i], one_indexed=one_indexed,
                                                       include_end=include_end, require_dim=2)
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
            raise ValueError('No archived sky spectrum defined for spectrograph {0}'.format(
                                        self.spectrograph))
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

