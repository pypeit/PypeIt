"""
Defines the abstract `Spectrograph` class, which is the parent class for
all instruments served by PypeIt.

The key functionality of this base class and its derived classes are to
provide instrument-specific:
    - file I/O routines
    - detector properties (see
      :class:`pypeit.par.pypeitpar.DetectorPar`)
    - telescope properties (see
      :class:`pypeit.par.pypeitpar.TelescopePar`)
    - fits header keywords that are collated and injested into PypeIt's
      metadata table that it uses throughout the reduction
    - header keyword values to check to confirm a fits file has been
      taken with the selected instrument
    - default methods for automatically determining the type of each
      exposure that PypeIt was asked to reduce
    - header keywords to use when matching calibration frames to science
      frames
    - methods used to generate and/or read bad-pixel masks for an
      exposure
    - default parameters for PypeIt's algorithms
    - method to access an archival sky spectrum

.. _astropy.io.fits: http://docs.astropy.org/en/stable/io/fits/
.. _astropy.io.fits.Header: http://docs.astropy.org/en/stable/io/fits/api/headers.html
.. _numpy.ndarray: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html

"""
import os
import warnings

from abc import ABCMeta
from pkg_resources import resource_filename

import numpy as np
from astropy.io import fits

from linetools.spectra import xspectrum1d

from pypeit import msgs
from pypeit.core.wavecal import wvutils
from pypeit.core import parse
from pypeit.core import procimg
from pypeit.par import pypeitpar
from pypeit.core import pixels
from pypeit.metadata import PypeItMetaData

from IPython import embed

class Spectrograph(object):
    """
    Abstract class whose derived classes dictate instrument-specific
    behavior in PypeIt.

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
        rawdatasec_img (:obj:`numpy.ndarray`):
            An image identifying the amplifier that reads each detector pixel.
        oscansec_img (:obj:`numpy.ndarray`):
            An image identifying the amplifier that reads each detector pixel
    """
    __metaclass__ = ABCMeta

    def __init__(self):
        self.spectrograph = 'base'
        self.telescope = None
        self.detector = None
        self.naxis = None
#        self.raw_naxis = None
        self.rawdatasec_img = None
        self.oscansec_img = None

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

    def nonlinear_counts(self, det, datasec_img=None):
        """
        Return the counts at which the detector response becomes
        non-linear.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
        
        Returns:
            float or np.ndarray:
                Counts at which detector response becomes nonlinear.
                If datasec_img is provided, an image with the same shape
                is returned with the gain accounted for
        """
        nonlinear_counts = self.detector[det-1]['saturation']*self.detector[det-1]['nonlinear']
        # Generate an image, applying the gain?
        if datasec_img is not None:
            gain = np.atleast_1d(self.detector[det-1]['gain']).tolist()
            nonlinear_counts = nonlinear_counts * procimg.gain_frame(datasec_img, gain)
        # Return
        return nonlinear_counts

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.
        
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
        return self.default_pypeit_par() if inp_par is None else inp_par

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

    def load_raw_extras(self, hdu, det):
        """
        Method to load up other bits and pieces related to the raw image
        Bundling them avoids multiple calls on the image

        Args:
            raw_file (str):
            det (int):

        Returns:
            tuple:
                exptime (float)
                rawdatasec_img  (np.ndarray)
                oscansec_img (np.ndarray)
                binning_raw (str)
        """
        # Binning
        headarr = self.get_headarr(hdu)
        binning = self.get_meta_value(headarr, 'binning')
        if self.detector[det-1]['specaxis'] == 1:
            binning_raw = (',').join(binning.split(',')[::-1])
        else:
            binning_raw = binning

        # Exposure time
        exptime = self.get_meta_value(headarr, 'exptime')

        # Rawdatasec, oscansec images
        rawdatasec_img = self.get_pixel_img(hdu, det, 'datasec', binning_raw)
        oscansec_img = self.get_pixel_img(hdu, det, 'oscansec', binning_raw)

        # Return it all
        return exptime, rawdatasec_img, oscansec_img, binning_raw

    def load_raw_frame(self, raw_file, det=1):
        r"""
        Load the image and header for an exposure taken with this
        spectrograph.

        Args:
            raw_file (:obj:`str`):
                File with the image data.
            det (:obj:`int`, optional):
                1-indexed detector number.

        Returns:
            Returns an `numpy.ndarray`_ with the image data and
            `astropy.io.fits.HDUList`_ object with the image and header
            data, respectively.  The image data is always returned with
            floating-point type.
        """
        # Check the detector is defined
        self._check_detector()

        hdu = fits.open(raw_file)
        raw_img = hdu[self.detector[det-1]['dataext']].data
        # Turn to float
        img = raw_img.astype(float)

        # Return
        return img, hdu

    def raw_is_transposed(self, det=1):
        """
        Indicates that raw files read by `astropy.io.fits`_ yields an
        image with the spatial dimension along rows, meaning that the
        image must be transposed to match the uniform PypeIt format of
        the spectral dimension along rows.

        Args:
            det (:obj:`int`, optional):
                1-indexed detector number.

        Returns:
            :obj:`bool`: Flat that transpose is required.
        """
        return self.detector[det-1]['specaxis'] == 1

    def get_image_section(self, hdulist, det, section):
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
            hdulist (`astropy.io.fits.Header`_):
                Header list of the image
            det (:obj:`int`):
                1-indexed detector number.
            section (:obj:`str`):
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
        # Check the section is one of the detector keywords
        if section not in self.detector[det-1].keys():
            raise ValueError('Unrecognized keyword: {0}'.format(section))

        # Check the detector is defined
        self._check_detector()

        # Get the data section
        hdr = hdulist[self.detector[det-1]['dataext']].header
        # Try using the image sections as header keywords
        try:
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
        #transpose = self.detector[det-1]['specaxis'] == 0

        return image_sections, one_indexed, include_last

    '''
    # THIS WILL PROBABLY NEED TO COME BACK
    def get_datasec_img(self, filename, det):
        """
        Generate and return the datasec image in the PypeIt reference
        frame, e.g. trimmed + oriented

        Returns:
            np.ndarray

        """
        rdimg = self.get_rawdatasec_img(filename=filename, det=det)
        # Fuss
        rdimg = procimg.trim_frame(rdimg, rdimg < 1)
        dimg = self.orient_image(rdimg, det)
        # Return
        return dimg
    '''

    def get_pixel_img(self, hdulist, det, section, binning_raw):
        """
        Create an image identifying the amplifier used to read each pixel.
        This is in the *raw* data format

        .. todo::
            - I find 1-indexing to be highly annoying...
            - Check for overlapping amplifiers?
            - Consider renaming this datasec_ampid or something like
              that.  I.e., the image's main purpose is to tell you where
              the amplifiers are for the data section

        Note on data format
        --------------------
         binning_pypeit = the binning  in the PypeIt convention of (spec, spat)
         binning_raw = the binning in the format of the raw data.
         In other words: PypeIt requires spec to be the first dimension of the image as read into python. If the
         files are stored the other way with spat as the first dimension (as read into python), then the transpose
         flag manages this, which is basically the value of the self.detector[det-1]['specaxis'] above.
         (Note also that BTW the python convention of storing images is transposed relative to the fits convention
         and the datasec typically written to headers. However this flip is dealt with explicitly in the
         parse.spec2slice code and is NOT the transpose we are describing and flipping here).
         TODO Add a blurb on the PypeIt data model.

        Args:
            filename (str):
                Name of the file from which to read the image size.
            section (str):  'datasec' or 'oscansec'
            det (int):
                Detector number (1-indexed)

        Returns:
            `numpy.ndarray`: Integer array identifying the amplifier
            used to read each pixel.
        """
        # Check the detector is defined
        self._check_detector()
        # Get the image shape
        #raw_naxis = self.get_raw_image_shape(filename, det=det)
        #binning_raw = self.get_meta_value(filename, 'binning')

        raw_naxis = self.get_raw_image_shape(hdulist, det=det)

        data_sections, one_indexed, include_end \
                    = self.get_image_section(hdulist, det, section=section)

        # Initialize the image (0 means no amplifier)
        pix_img = np.zeros(raw_naxis, dtype=int)
        for i in range(self.detector[det-1]['numamplifiers']):
            # Convert the data section from a string to a slice
            datasec = parse.sec2slice(data_sections[i], one_indexed=one_indexed,
                                      include_end=include_end, require_dim=2,
                                      binning=binning_raw)
            # Assign the amplifier
            pix_img[datasec] = i+1

        return pix_img

    def get_raw_image_shape(self, hdulist, det=None):
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
            hdulist (:obj:`str`, optional):
                Name of the fits file with the header to use.
            det (:obj:`int`, optional):
                1-indexed number of the detector.  Default is None.  If
                None, the primary extension is used.  Otherwise the
                internal detector parameters are used to determine the
                extension to read.
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
        #hdu = fits.open(filename)
        header = hdulist[self.detector[det-1]['dataext']].header
        shape = (header['NAXIS2'], header['NAXIS1'])  # Usual Python vs. the world
        return shape

    def header_cards_for_spec(self):
        """
        Define the header cards to be written to spec1d (and maybe spec2d) files.
        These refer to the keys in the fitstbl

        Returns:
            list: Keys for header cards of spec1d

        """
        core_meta = PypeItMetaData.define_core_meta()
        header_cards = list(core_meta.keys())
        # Add a few more
        header_cards += ['filename']  # For fluxing
        return header_cards

    def orient_image(self, rawimage, det):
        """
        Orient the image into the PypeIt frame

        Args:
            rawimage (np.ndarray):
                Image in the raw frame
            det (int):
                Detector index

        Returns:
            np.ndarray:  Oriented image

        """
        image = rawimage.copy()
        # Transpose?
        if self.raw_is_transposed(det):
            image = image.T
        # Flip spectral axis?
        if self.detector[det-1]['specflip'] is True:
            image = np.flip(image, axis=0)
        # Flip spatial axis?
        if self.detector[det-1]['spatflip'] is True:
            image = np.flip(image, axis=1)
        return image

    def empty_bpm(self, filename, det):
        """
        Generate a generic (empty) bad-pixel mask.

        Even though they are both optional, either the precise shape for
        the image (`shape`) or an example file that can be read to get
        the shape (`filename` using :func:`get_image_shape`) *must* be
        provided.

        Args:
            filename (:obj:`str`):
                An example file to use to get the image shape.
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to 0.
        """
        # Load the raw frame
        hdu = fits.open(filename)
        exptime, rawdatasec_img, oscansec_img, _ = self.load_raw_extras(hdu, det)
        # Trim + reorient
        trim = procimg.trim_frame(rawdatasec_img, rawdatasec_img < 1)
        orient = self.orient_image(trim, det)

        # Generate
        bpm_img = np.zeros(orient.shape, dtype=np.int8)

        # Return
        return bpm_img

    def bpm(self, filename, det):
        """
        Generate a default bad-pixel mask.

        Currently identical to calling :func:`empty_bpm`.

        Even though they are both optional, either the precise shape for
        the image (`shape`) or an example file that can be read to get
        the shape (`filename` using :func:`get_image_shape`) *must* be
        provided.

        Args:
            filename (:obj:`str`):
                An example file to use to get the image shape.
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        return self.empty_bpm(filename, det)

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

    def get_meta_value(self, inp, meta_key, required=False, ignore_bad_header=False, usr_row=None):
        """
        Return meta data from a given file (or its array of headers)

        Args:
            inp (str or list):
              Input filename or headarr list
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
        if isinstance(inp, str):
            headarr = self.get_headarr(inp)
        else:
            headarr = inp

        # Loop?
        if isinstance(meta_key, list):
            values = []
            for mdict in meta_key:
                values.append(self.get_meta_value(headarr, mdict, required=required))
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
                    msgs.warn("Required card {:s} missing from your header.  Proceeding with risk..".format(
                        self.meta[meta_key]['card']))
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
            assert key in self.meta.keys(), \
                'key {:s} not defined in spectrograph meta!'.format(key)
        # Check for rtol for config keys that are type float
        for key in self.configuration_keys():
            if meta_data_model[key]['dtype'] in [float]:
                assert 'rtol' in self.meta[key].keys(), \
                    'rtol not set for key {:s} not defined in spectrograph meta!'.format(key)
        # Now confirm all meta are in the data model
        for key in self.meta.keys():
            if key not in self.meta_data_model.keys():
                msgs.error("Meta data {:s} not in meta_data_model".format(key))

    def get_headarr(self, inp, strict=True):
        """
        Read the header data from all the extensions in the file.

        Args:
            inp (:obj:`str` or hdulist):
                Name of the file to read or the hdulist
            strict (:obj:`bool`, optional):
                Function will fault if :func:`fits.getheader` fails to
                read any of the headers.  Set to False to report a
                warning and continue.

        Returns:
            list: Returns a list of :attr:`numhead` :obj:`fits.Header`
            objects with the extension headers.
        """
        # Faster to open the whole file and then assign the headers,
        # particularly for gzipped files (e.g., DEIMOS)
        if isinstance(inp, str):
            try:
                hdu = fits.open(inp)
            except:
                if strict:
                    msgs.error('Problem opening {0}.'.format(inp))
                else:
                    msgs.warn('Problem opening {0}.'.format(inp) + msgs.newline()
                              + 'Proceeding, but should consider removing this file!')
                    return ['None']*self.numhead
        else:
            hdu = inp
        return [hdu[k].header for k in range(self.numhead)]

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

    @property
    def ndet(self):
        """Return the number of detectors."""
        return 0 if self.detector is None else len(self.detector)

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


    def order_platescale(self, order_vec, binning=None):
        """
        This routine is only for echelle spectrographs. It returns the plate scale order by order

        Args:
            order_vec (np.ndarray):
            binning:

        Returns:
            np.ndarray

        """
        pass

    def order_vec(self, slit_spat_pos):
        """
        Convert an array of slit_spat_pos values to order numbers

        Args:
            slit_spat_pos (np.ndarray): Slit positions

        Returns:
            np.ndarray: Order numbers

        """
        order_vec = np.zeros(slit_spat_pos.size, dtype=int)
        for kk, ipos in enumerate(slit_spat_pos):
            order_vec[kk], indx= self.slit2order(ipos)
        # Return
        return order_vec

    @property
    def norders(self):
        pass

    @property
    def order_spat_pos(self):
        pass

    @property
    def orders(self):
        pass

    @property
    def spec_min_max(self):
        return None

    @property
    def dloglam(self):
        pass

    @property
    def loglam_minmax(self):
        pass

    def slit2order(self, slit_spat_pos):
        """
        This routine is only for fixed-format echelle spectrographs.
        It returns the order of the input slit based on its slit_pos

        Args:
            slit_spat_pos (float):  Slit position (spatial at 1/2 the way up)

        Returns:
            order, indx

            order (int): order number
            indx  (int): order index

        """
        indx = np.arange(self.norders)
        # Find closest
        try:
            iorder = [np.argmin(np.abs(slit-self.order_spat_pos)) for slit in slit_spat_pos]
        except TypeError:
            iorder = np.argmin(np.abs(slit_spat_pos-self.order_spat_pos))

        # Check
        if np.any(np.abs(self.order_spat_pos[iorder] - slit_spat_pos) > 0.05):
            msgs.warn("Bad echelle format for VLT-XSHOOTER or you are performing a 2-d coadd with different order locations."
                      "Returning order vector with the same number of orders you requested")
            iorder = np.arange(slit_spat_pos.size)
            return self.orders[iorder], indx[iorder]
        else:
            return self.orders[iorder], indx[iorder]


    def slit_minmax(self, slit_spat_pos, binspectral=1):
        """

        Args:
            slit_spat_pos (float or ndarray):
                normalized slit_spatial position as computed by trace_slits.slit_spat_pos
            binspectral (int): default=1
               spectral binning

        Returns:

        """

        if self.spec_min_max is None:
            try:
                nslit = len(slit_spat_pos)
            except TypeError:
                nslit = 1
            return np.vstack((np.asarray([-np.inf]*nslit), np.asarray([np.inf]*nslit)))

        else:
            try:
                iorder = [np.argmin(np.abs(slit-self.order_spat_pos)) for slit in slit_spat_pos]
            except TypeError:
                iorder = np.argmin(np.abs(slit_spat_pos-self.order_spat_pos))
            return self.spec_min_max[:, iorder]/binspectral


    def wavegrid(self, binning=None, midpoint=False,samp_fact=1.0):
        """
        Routine to generate a fixed wavelength grid in log_10 lambda. Mostly used by echelle spectrographs

        Args:
            binning:
            midpoint:
            samp_fact:

        Returns:

        """

        binspectral, binspatial = parse.parse_binning(binning)
        logmin, logmax = self.loglam_minmax
        loglam_grid = wvutils.wavegrid(logmin, logmax, self.dloglam*binspectral, samp_fact=samp_fact)
        if midpoint:
            loglam_grid = loglam_grid + self.dloglam*binspectral/samp_fact/2.0

        return np.power(10.0,loglam_grid)


    def __repr__(self):
        # Generate string
        txt = '<{:s}: '.format(self.__class__.__name__)
        txt += ' spectrograph={:s},'.format(self.spectrograph)
        txt += ' telescope={:s},'.format(self.telescope['name'])
        txt += '>'
        return txt

