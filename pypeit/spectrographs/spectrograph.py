""" Module to define the Spectrograph class
"""
from __future__ import absolute_import, division, print_function

import os
import warnings

from abc import ABCMeta
from pkg_resources import resource_filename

import numpy as np
from astropy.io import fits

from linetools.spectra import xspectrum1d

from pypeit import msgs
from pypeit.core import parse
from pypeit.par import pypeitpar
from pypeit.core import pixels
from pypeit.metadata import PypeItMetaData

from pypeit import debugger

class Spectrograph(object):
    """
    Generic class for spectrograph-specific codes

    Attributes:
        spectrograph (str):
            The name of the spectrograph.  See
            :func:`pypeit.spectrographs.util.valid_spectrographs` for the
            currently supported spectrographs.
        telescope (:class:`TelescopePar`):
            Parameters of the telescope that feeds this spectrograph.
        detector (list):
            A list of instances of
            :class:`pypeit.par.pypeitpar.DetectorPar` with the parameters
            for each detector in the spectrograph
        naxis (tuple):
            A tuple with the lengths of the two axes for current
            detector image; often trimmmed.
        raw_naxis (tuple):
            A tuple with the lengths of the two axes for untrimmed detector image.
        datasec_img (:obj:`numpy.ndarray`):
            An image identifying the amplifier that reads each detector
            pixel.
        bpm_img (:obj:`numpy.ndarray`):
            The bad-pixel mask for the currently read detector.
    """
    __metaclass__ = ABCMeta

    def __init__(self):
        self.spectrograph = 'base'
        self.telescope = None
        self.detector = None
        self.naxis = None
#        self.raw_naxis = None
        self.datasec_img = None
        self.bpm_img = None

        # Default time unit
        self.timeunit = 'mjd'

        # Default extension with the primary header data
        #   used by arsave.save_2d_images
        self.primary_hdrext = 0
        self.numhead = 0

        self.minexp = 0  # NEED TO TIE TO INSTRUMENT PAR INSTEAD

        # Init Calibrations Par
#        self._set_calib_par()

        # Init meta
        self.meta_data_model = PypeItMetaData.get_meta_data_model()
        self.init_meta()
        self.validate_metadata()

    @staticmethod
    def default_pypeit_par():
        return pypeitpar.PypeItPar()

    def config_specific_par(self, par, scifile):
        """
        Used to modify the ParSet from metadata
        drawn from the input file

        Args:
            par: ParSet
            scifile: str

        Returns:
            par

        """
        return par

    def _check_telescope(self):
        # Check the detector
        if self.telescope is None:
            raise ValueError('Must define the telescope used to take the observations.')
        if not isinstance(self.telescope, pypeitpar.TelescopePar):
                raise TypeError('Telescope parameters must be one of those specified in'
                                'pypeit.telescopes.')

    def _check_detector(self):
        # Check the detector
        if self.detector is None:
            raise ValueError('Must first define spectrograph detector parameters!')
        for d in self.detector:
            if not isinstance(d, pypeitpar.DetectorPar):
                raise TypeError('Detector parameters must be specified using DetectorPar.')

#    def _set_calib_par(self, user_supplied=None):
#        pass

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
        if self.detector[_det-1]['specaxis'] == 1:
            img = img.T
        if self.detector[_det-1]['specflip'] is True:
            img = np.flip(img, axis=0)
        if self.detector[_det-1]['spatflip'] is True:
            img = np.flip(img, axis=1)

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

    def get_image_section(self, inp=None, det=1, section='datasec'):
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
            inp (:obj:`str`, `astropy.io.fits.Header`, optional):
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
            A list of string representations for the image sections, one
            string per amplifier, followed by three booleans: if the
            slices are one indexed, if the slices should include the
            last pixel, and if the slice should have their order
            transposed.
        """
        # Check the section is one of the detector keywords
        if section not in self.detector[det-1].keys():
            raise ValueError('Unrecognized keyword: {0}'.format(section))

        # Check the detector is defined
        self._check_detector()

        # Get the data section
        try:
            # Parse inp
            if inp is None:
                # Force the call to the except block
                raise KeyError
            elif isinstance(inp, str):
                hdu = fits.open(inp)
                hdr = hdu[self.detector[det-1]['dataext']].header
            elif isinstance(inp, fits.Header):
                hdr = inp
            else:
                msgs.error('Input must be a filename or a fits.Header object.')

            # Try using the image sections as header keywords
            image_sections = [ hdr[key] for key in self.detector[det-1][section] ]
        except KeyError:
            # Expect a KeyError to be thrown if the section defined by
            # the detector attribute is actually the image section
            # string itself
            image_sections = self.detector[det-1][section]
            if not isinstance(image_sections, list):
                image_sections = [image_sections]

        # Always assume normal FITS header formatting
        one_indexed = True
        include_last = True
        transpose = self.detector[det-1]['specaxis'] == 0

        return image_sections, one_indexed, include_last, transpose

    def get_datasec_img(self, filename, det=1, force=True):
        """
        Create an image identifying the amplifier used to read each pixel.

        .. todo::
            - I find 1-indexing to be highly annoying...
            - Check for overlapping amplifiers?
            - Consider renaming this datasec_ampid or something like
              that.  I.e., the image's main purpose is to tell you where
              the amplifiers are for the data section
          
        Args:
            filename (str):
                Name of the file from which to read the image size.
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
            raw_naxis = self.get_raw_image_shape(filename, det=det)

            binning = self.get_meta_value(filename, 'binning')
#            binning = self.parse_binning(filename)

            data_sections, one_indexed, include_end, transpose \
                    = self.get_image_section(filename, det, section='datasec')

            # Initialize the image (0 means no amplifier)
            self.datasec_img = np.zeros(raw_naxis, dtype=int)
            for i in range(self.detector[det-1]['numamplifiers']):
                # Convert the data section from a string to a slice
                datasec = parse.sec2slice(data_sections[i], one_indexed=one_indexed,
                                          include_end=include_end, require_dim=2,
                                          transpose=transpose, binning=binning)
                # Assign the amplifier
                self.datasec_img[datasec] = i+1
        return self.datasec_img

    def get_raw_image_shape(self, filename, det=None, force=True):
        """
        Get the *untrimmed* shape of the image data for a given detector using a
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

        # Use a file
        self._check_detector()
        raw_naxis = (self.load_raw_frame(filename, det=det)[0]).shape
        return raw_naxis

    def empty_bpm(self, shape=None, filename=None, det=1):
        """
        Generate a generic (empty) BPM.
        If shape is None, this requires a successful call to :func:`get_image_shape`.

        .. todo::
            Any reason this isn't returned as a boolean array?

        Args:
            shape: tuple, REQUIRED

        Returns:
            bpm: ndarray, int
              0=not masked; 1=masked

        """

        # TODO The logic heere is botched. shape has to be set for the code not to crash so why is filename even
        # an argument???
        if shape is None:
            msgs.error("THIS IS NOT GOING TO WORK")
            # Check the detector is defined
            self._check_detector()
            # Get the image shape
            _shape = self.get_raw_image_shape(filename, det=det)
        else:
            _shape = shape
        # JFH I think all masks should be boolean aside from the bitmask.
        self.bpm_img = np.zeros(_shape, dtype=np.int8)
        # Return
        return self.bpm_img

    def bpm(self, shape=None, filename=None, det=1):
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
        ## JFH TODO Filename should not be an option since it can never be used in empty_bpm since the code crashes
        # before that.
        return self.empty_bpm(shape=shape, filename=filename, det=det)

    '''
    # TODO: (KBW) I've removed all the defaults.  Should maybe revisit
    # this
    def default_header_keys(self):
        def_head_keys = {}
        def_head_keys[0] = {}
        return def_head_keys

    def header_keys(self):
        return self.default_header_keys()
    '''

    def configuration_keys(self):
        """
        Return the metadata keys that defines a unique instrument
        configuration.

        This list is used by :class:`pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:

            list: List of keywords of data pulled from file headers and
            used to constuct the :class:`pypeit.metadata.PypeItMetaData`
            object.
        """
        return ['dispname', 'dichroic', 'decker']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file

        Returns:
            pypeit_keys: list

        """
        pypeit_keys = ['filename', 'frametype']
        # Core
        core_meta = PypeItMetaData.define_core_meta()
        pypeit_keys += list(core_meta.keys())  # Might wish to order these
        # Add in config_keys (if new)
        for key in self.configuration_keys():
            if key not in pypeit_keys:
                pypeit_keys.append(key)
        # Finish
        return pypeit_keys

    def compound_meta(self, headarr, meta_key):
        """
        Methods to generate meta in a more complex manner than simply
        reading from the header

        These are defined per spectrograph, as needed

        Args:
            headarr: list
              List of headers
            meta_key: str

        Returns:
            value:

        """
        return None

    def init_meta(self):
        """
        Define how meta values are dervied from the spectrograph files

        Returns:
            self.meta defined

        """
        self.meta = {}

    # TODO: Change this so that it uses one argument for ifile or
    # headarray? ala, get_image_section, etc?
    def get_meta_value(self, ifile, meta_key, headarr=None, required=False, ignore_bad_header=False,
                       usr_row=None):
        """
        Return meta data from a given file (or its array of headers)

        Args:
            ifile: str or None
              Input filename
            meta_key: str or list of str
            headarr: list, optional
              List of headers
            required: bool, optional
              Require the meta key to be returnable
            ignore_bad_header: bool, optional
              Over-ride required;  not recommended
            usr_row: Row
              Provides user supplied frametype (and other things not used)

        Returns:
            value: value or list of values

        """
        if headarr is None:
            headarr = self.get_headarr(ifile)

        # Loop?
        if isinstance(meta_key, list):
            values = []
            for mdict in meta_key:
                values.append(self.get_meta_value(mdict, headarr=headarr, required=required))
            #
            return values

        # Are we prepared to provide this meta data?
        if meta_key not in self.meta.keys():
            if required:
                msgs.error("Need to allow for meta_key={} in your meta data".format(meta_key))
            else:
                msgs.warn("Requested meta data does not exist...")
                return None
        # Is this not derivable?  If so, use the default
        #   or search for it as a compound method
        value = None
        if self.meta[meta_key]['card'] is None:
            if 'default' in self.meta[meta_key].keys():
                value = self.meta[meta_key]['default']
            elif 'compound' in self.meta[meta_key].keys():
                value = self.compound_meta(headarr, meta_key)
            else:
                msgs.error("Failed to load spectrograph value for meta: {}".format(meta_key))
        else:
            # Grab from the header, if we can
            try:
                value = headarr[self.meta[meta_key]['ext']][self.meta[meta_key]['card']]
            except (KeyError, TypeError):
                value = None

        if value is None:
            # Was this required?
            if required:
                kerror = True
                if not ignore_bad_header:
                    # Is this meta required for this frame type (Spectrograph specific)
                    if ('required_ftypes' in self.meta[meta_key]) and (usr_row is not None):
                        kerror = False
                        # Is it required?
                        for ftype in usr_row['frametype'].split(','):
                            if ftype in self.meta[meta_key]['required_ftypes']:
                                kerror = True
                    # Bomb out?
                    if kerror:
                        msgs.error('Required meta "{:s}" did not load!  You may have a corrupt header'.format(meta_key))
                else:
                    msgs.warn("Required card {:s} missing from your header of {:s}.  Proceeding with risk..".format(
                        self.meta[meta_key]['card'], ifile))
            return None

        # Deal with dtype (DO THIS HERE OR IN METADATA?  I'M TORN)
        if self.meta_data_model[meta_key]['dtype'] == str:
            value = str(value).strip()
        elif self.meta_data_model[meta_key]['dtype'] == int:
            value = int(value)
        elif self.meta_data_model[meta_key]['dtype'] == float:
            value = float(value)
        elif self.meta_data_model[meta_key]['dtype'] == tuple:
            assert isinstance(value, tuple)
        else:
            debugger.set_trace()
        # Return
        return value

    def validate_metadata(self):
        """
        Validates the meta definitions of the Spectrograph
        by making a series of comparisons to the meta data model
        definied in metadata.py

        Returns:

        """
        # Load up
        core_meta = PypeItMetaData.define_core_meta()
        meta_data_model = PypeItMetaData.get_meta_data_model()
        # Check core
        for key in core_meta:
            assert key in self.meta.keys(), 'key {:s} not defined in spectrograph meta!'.format(key)
        # Check for rtol for config keys that are type float
        for key in self.configuration_keys():
            if meta_data_model[key]['dtype'] in [float]:
                assert 'rtol' in self.meta[key].keys(), 'rtol not set for key {:s} not defined in spectrograph meta!'.format(key)
        # Now confirm all meta are in the data model
        for key in self.meta.keys():
            if key not in self.meta_data_model.keys():
                msgs.error("Meta data {:s} not in meta_data_model".format(key))

    def get_headarr(self, filename, strict=True):
        """
        Read the header data from all the extensions in the file.

        Args:
            filename (:obj:`str`):
                Name of the file to read.
            strict (:obj:`bool`, optional):
                Function will fault if :func:`fits.getheader` fails to
                read any of the headers.  Set to False to report a
                warning and continue.

        Returns:
            list: Returns a list of :attr:`numhead` :obj:`fits.Header`
            objects with the extension headers.
        """
        headarr = ['None']*self.numhead
        for k in range(self.numhead):
            try:
                headarr[k] = fits.getheader(filename, ext=k)
            except:
                if strict:
                    msgs.error("Header error in extension {0} in {1}.".format(k, filename))
                else:
                    msgs.warn('Bad header in extension {0} in {1}'.format(k, filename) 
                              + msgs.newline() + 'Proceeding on the hopes this was a '
                              + 'calibration file, otherwise consider removing.')
        return headarr

#    def get_match_criteria(self):
#        msgs.error("You need match criteria for your spectrograph.")

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        raise NotImplementedError('Frame typing not defined for {0}.'.format(self.spectrograph))

    def idname(self, ftype):
        """
        Return the `idname` for the selected frame type for this instrument.

        Args:
            ftype (str):
                File type, which should be one of the keys in
                :class:`pypeit.core.framematch.FrameTypeBitMask`.

        Returns:
            str: The value of `idname` that should be available in the
            `PypeItMetaData` instance that identifies frames of this
            type.
        """
        raise NotImplementedError('Header keyword with frame type not defined for {0}.'.format(
                                  self.spectrograph))

#    def parse_binning(self, inp, det=1, key='BINNING'):
#        """
#        Get the pixel binning for an image.
#
#        Args:
#            inp (:obj:`str`, `astropy.io.fits.Header`):
#                String providing the file name to read, or the relevant
#                header object.
#            det (:obj:`int`, optional):
#                1-indexed detector number.
#            key (:obj:`str`, optional):
#                Header key with the binning.  This is included as an
#                argument in the base-class implementation so that it can
#                be called by the derived classes for most cases, when
#                the binning is in a single keyword.
#
#        Returns:
#            str: String representation of the binning.  The ordering is
#            as provided in the header, regardless of which axis is
#            designated as the dispersion axis.  It is expected that this
#            be used with :func:`pypeit.core.parse.sec2slice` to setup
#            the data and overscane sections of the image data.
#
#        Raises:
#            PypeItError:
#                Raised if `inp` is not one of the accepted types.
#        """
#        # Get the header
#        # TODO: Read primary header by default instead?
#        if isinstance(inp, str):
#            hdu = fits.open(inp)
#            hdr = hdu[self.detector[det-1]['dataext']].header
#        elif isinstance(inp, fits.Header):
#            hdr = inp
#        else:
#            msgs.error('Input must be a filename or fits.Header object')
#
#        # Parse the keyword
#        binning = parse.parse_binning(hdr[key])
#
#        # Return comma-separated string
#        return ','.join([ str(b) for b in binning])

    @property
    def ndet(self):
        """Return the number of detectors."""
        return 0 if self.detector is None else len(self.detector)

    '''
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
                pypeit/data/sky_spec/ directory in the source
                distribution.
        """
        # No file was defined
        if self.sky_file is None:
            self.sky_file = Spectrograph.default_sky_spectrum()
            warnings.warn('Using default sky spectrum: {0}'.format(self.sky_file))

        if os.path.isfile(self.sky_file):
            # Found directly
            return self.sky_file, xspectrum1d.XSpectrum1D.from_file(self.sky_file)

        root = resource_filename('pypeit', 'data/sky_spec/')
        _sky_file = os.path.join(root, self.sky_file)
        if os.path.isfile(_sky_file):
            # Found within the source distribution
            return self.sky_file, xspectrum1d.XSpectrum1D.from_file(_sky_file)

        # File could not be read
        raise FileNotFoundError('Could not find archive sky spectrum: {0} or {1}'.format(
                                    self.sky_file, _sky_file))
    '''

    @property
    def pypeline(self):
        return 'MultiSlit'

    def mm_per_pix(self, det=1):
        """
        Return the spatial scale at the telescope focal plane in mm per
        pixel at the detector.

        The fratio and diameter of the telescope must be defined.

        Args:
            det (:obj:`int`, optional):
                Detector to use for the spectrograph platescale.

        Returns:
            float: The spatial scale at the telescope focal plane in mm
            per detector pixel scale.
        
        Raises:
            ValueError: 
                Raised if the telescope is undefined, any of the numbers
                needed for the calculation are not available, or the
                selected detector is out of range.
        """
        if det > self.ndet:
            raise ValueError('Selected detector out of range; det={0}..{1}.'.format(1,self.ndet))
        tel_platescale = None if self.telescope is None else self.telescope.platescale()
        if self.telescope is None or tel_platescale is None or \
                self.detector[det-1]['platescale'] is None:
            raise ValueError('Incomplete information to calculate mm per pixel.')

        return self.detector[det-1]['platescale']/tel_platescale

    def slit_minmax(self, nslits, binspectral=1):
        """
        Generic routine to determine the minimum and maximum spectral pixel for slitmasks. This functionality
        is most useful for echelle observations where the orders cutoff. This generic routine will be used for most
        slit spectrographs and simply sets the minimum and maximum spectral pixel for the slit to be -+ infinity.

        Args:
            binspectral:

        Returns:

        """
        spec_min = np.full(nslits, -np.inf)
        spec_max = np.full(nslits, np.inf)
        return spec_min, spec_max


    def slitmask(self, tslits_dict, pad = None):
        """
         Generic routine to construct a slitmask image from a tslits_dict. Children of this class can
         overload this function to implement instrument specific slitmask behavior, for example setting
         where the orders on an echelle spectrograph end

         Parameters
         -----------
         tslits_dict: dict
            Trace slits dictionary with slit boundary information

         Optional Parameters
         pad: int or float
            Padding of the slit boundaries

         Returns
         -------
         slitmask: ndarray int
            Image with -1 where there are no slits/orders, and an integer where there are slits/order with the integer
            indicating the slit number going from 0 to nslit-1 from left to right.

         """

        # These lines are always the same
        pad = tslits_dict['pad'] if pad is None else pad
        slitmask = pixels.slit_pixels(tslits_dict['lcen'],tslits_dict['rcen'],tslits_dict['nspat'], pad=pad)
        return slitmask

    # This routine is only for echelle spectrographs. It returns the plate scale order by order
    def order_platescale(self, binning=None):
        pass

    # This routine is only for echelle spectrographs. It returns the plate scale order by order
    def slit2order(self, slit):
        pass


    def __repr__(self):
        # Generate string
        txt = '<{:s}: '.format(self.__class__.__name__)
        txt += ' spectrograph={:s},'.format(self.spectrograph)
        txt += ' telescope={:s},'.format(self.telescope['name'])
        txt += '>'
        return txt

