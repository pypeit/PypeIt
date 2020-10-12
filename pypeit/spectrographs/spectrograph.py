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

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os
from copy import deepcopy
import warnings

from abc import ABCMeta
from pkg_resources import resource_filename

import numpy as np
from astropy.io import fits

from pypeit import msgs
from pypeit import utils
from pypeit.core.wavecal import wvutils
from pypeit.core import parse
from pypeit.core import procimg
from pypeit.core import meta
from pypeit.par import pypeitpar

from IPython import embed

# TODO: Create an EchelleSpectrograph derived class that holds all of
# the echelle specific methods.

class Spectrograph:
    # TODO: This docstring needs to be updated.
    """
    Abstract base class whose derived classes dictate
    instrument-specific behavior in PypeIt.

    Attributes:
        spectrograph (:obj:`str`):
            The name of the spectrograph. See
            :func:`pypeit.spectrographs.util.valid_spectrographs` for
            the currently supported spectrographs.
        telescope (:class:`TelescopePar`):
            Parameters of the telescope that feeds this spectrograph.
        detector (:obj:`list`):
            A list of instances of
            :class:`pypeit.par.pypeitpar.DetectorPar` with the
            parameters for each detector in the spectrograph
        naxis (:obj:`tuple`):
            A tuple with the lengths of the two axes for current
            detector image; often trimmmed.
        raw_naxis (tuple):
            A tuple with the lengths of the two axes for untrimmed
            detector image.
        rawdatasec_img (:obj:`numpy.ndarray`):
            An image identifying the amplifier that reads each detector
            pixel.
        oscansec_img (:obj:`numpy.ndarray`):
            An image identifying the amplifier that reads each detector
            pixel
        slitmask (:class:`pypeit.spectrographs.slitmask.SlitMask`):
            Provides slit and object coordinate data for an
            observation. Not necessarily populated for all
            spectrograph instantiations.
    """
    __metaclass__ = ABCMeta

    ndet = None

    def __init__(self):
        self.spectrograph = 'base'
        self.camera = 'base'
        self.telescope = None
        self.camera = None
        self.dispname = None
        self.detector = None
        self.naxis = None
        self.rawdatasec_img = None
        self.oscansec_img = None
        self.slitmask = None

        # Default time unit
        self.timeunit = 'mjd'

        # Default extension with the primary header data
        #   used by arsave.save_2d_images
        self.primary_hdrext = 0

        # Init meta
        self.meta_data_model = meta.get_meta_data_model()
        self.init_meta()
        self.validate_metadata()

        # Validate detector
        assert self.ndet > 0

    @staticmethod
    def default_pypeit_par():
        return pypeitpar.PypeItPar()

    def nonlinear_counts(self, detector_par, datasec_img=None, apply_gain=True):
        """
        Return the counts at which the detector response becomes
        non-linear.

        Default is to apply the gain, i.e. return this is counts not ADU

        Args:
            detector_par (:class:`pypeit.par.pypeitpar.DetectorPar`):
            datasec_img (np.ndarray, optional):
                If provided, nonlinear_counts is returned as an image.
                DO NOT USE THIS OPTION; IT IS NOT YET IMPLEMENTED
                DOWNSTREAM.
            apply_gain (bool, optional):
                Apply gain in the calculation, i.e. convert to counts
                If only a float is returned, (i.e. no datasec_img is provided)
                then the mean of the gains for all amplifiers is adopted

        Returns:
            float, np.ndarray: Counts at which detector response becomes
            nonlinear.  If datasec_img is provided, an image with the
            same shape is returned
        """
        # Deal with gain
        gain = np.atleast_1d(detector_par['gain']).tolist()
        if not apply_gain:  # Set to 1 if gain is not to be applied
            gain = [1. for item in gain]
        # Calculation without gain
        nonlinear_counts = detector_par['saturation']*detector_par['nonlinear']
        # Finish
        if datasec_img is not None:  # 2D image
            nonlinear_counts = nonlinear_counts * procimg.gain_frame(datasec_img, gain)
        else:  # float
            nonlinear_counts = nonlinear_counts * np.mean(gain)
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
            :class:`pypeit.par.parset.ParSet`: The PypeIt parameter set
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

    def raw_is_transposed(self, detector_par):
        """
        Indicates that raw files read by `astropy.io.fits`_ yields an
        image with the spatial dimension along rows, meaning that the
        image must be transposed to match the uniform PypeIt format of
        the spectral dimension along rows.

        Args:
            detector_par (:class:`pypeit.par.pypeitpar.DetectorPar`):

        Returns:
            :obj:`bool`: Flag that transpose is required.
        """
        return detector_par['specaxis'] == 1

    def parse_spec_header(self, header):
        """
        Parses an input header for key spec items

        Args:
            header (:class:`astropy.io.fits.Header`):

        Returns:
            dict:

        """
        spec_dict = {}
        #
        core_meta_keys = list(meta.define_core_meta().keys())
        core_meta_keys += ['filename']
        for key in core_meta_keys:
            if key.upper() in header.keys():
                spec_dict[key.upper()] = header[key.upper()]
        # Return
        return spec_dict

    def subheader_for_spec(self, row_fitstbl, raw_header, extra_header_cards=[],
                           allow_missing=False):
        """
        Generate a dict that will be added to the Header of spectra
        files generated by PypeIt,
        e.g.  :class:`pypeit.specobjs.SpecObjs`

        Args:
            row_fitstbl (:class:`astropy.table.Row` or :class:`astropy.io.fits.Header`):
            raw_header (:class:`astropy.io.fits.Header`):
            extra_header_cards (list, optional):
                Additional header cards to add, if present

        Returns:
            :obj:`dict`: -- Used to generate a Header or table downstream
        """
        subheader = {}

        core_meta = meta.define_core_meta()
        # Core
        for key in core_meta.keys():
            try:
                subheader[key] = (row_fitstbl[key], core_meta[key]['comment'])
            except KeyError:
                if not allow_missing:
                    msgs.error("Key: {} not present in your fitstbl/Header".format(key))
        # Add a few more
        for key in ['filename']:  # For fluxing
            subheader[key] = row_fitstbl[key]

        # The following are pulled from the original header, if available
        header_cards = ['INSTRUME', 'DETECTOR'] + extra_header_cards  # For specDB and more
        for card in header_cards:
             if card in raw_header.keys():
                 subheader[card] = raw_header[card]  # Self-assigned instrument name

        # Specify which pipeline created this file
        subheader['PYPELINE'] = self.pypeline
        subheader['PYP_SPEC'] = (self.spectrograph, 'PypeIt: Spectrograph name')

        # Observatory and Header supplied Instrument
        telescope = self.telescope
        subheader['TELESCOP'] = (telescope['name'], 'Telescope')
        subheader['LON-OBS'] = (telescope['longitude'], 'Telescope longitude')
        subheader['LAT-OBS'] = (telescope['latitude'], 'Telescope latitute')
        subheader['ALT-OBS'] = (telescope['elevation'], 'Telescope elevation')

        # Return
        return subheader

    def orient_image(self, detector_par, rawimage):
        """
        Orient the image into the PypeIt spec,spat configuration

        Args:
            detector_par (:class:`pypeit.par.pypeitpar.DetectorPar`):
            rawimage (np.ndarray):
                Image in the raw frame
            det (int):
                Detector index

        Returns:
            np.ndarray:  Oriented image

        """
        image = rawimage.copy()
        # Transpose?
        if self.raw_is_transposed(detector_par):
            image = image.T
        # Flip spectral axis?
        if detector_par['specflip'] is True:
            image = np.flip(image, axis=0)
        # Flip spatial axis?
        if detector_par['spatflip'] is True:
            image = np.flip(image, axis=1)
        return image

    ## TODO: JFH Are these bad pixel masks in the raw frame, or the flipped/transposed pypeit frame??
    def empty_bpm(self, filename, det, shape=None):
        """
        Generate a generic (empty) bad-pixel mask.

        Even though they are both optional, either the precise shape for
        the image (`shape`) or an example file that can be read to get
        the shape (`filename` using :func:`get_image_shape`) *must* be
        provided.

        Args:
            filename (:obj:`str` or None):
                An example file to use to get the image shape.
                If None, shape must be provided
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.
            shape (tuple, optional):
                Processed image shape
                Required if filename is None
                Ignored if filename is not None

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to 0.
        """
        # Load the raw frame
        if filename is not None:
            detector_par, _,  _, _, rawdatasec_img, _ = self.get_rawimage(filename, det)
            # Trim + reorient
            trim = procimg.trim_frame(rawdatasec_img, rawdatasec_img < 1)
            orient = self.orient_image(detector_par, trim)#, det)
            #
            shape = orient.shape
        else: # This is risky if you don't really know what you are doing!
            if shape is None:
                msgs.error("Must specify shape if filename is None")

        # Generate
        bpm_img = np.zeros(shape, dtype=np.int8)

        # Return
        return bpm_img

    def bpm_frombias(self, msbias, det, bpm_img):
        """
        Generate a bad-pixel mask from a master bias frame.

        Args:
            msbias (`numpy.ndarray`):
                Master bias frame used to identify bad pixels
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.
            bpm_img (`numpy.ndarray`):
                bad pixel mask

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to 0.
        """
        msgs.info("Generating a BPM for det={0:d} on {1:s}".format(det, self.camera))
        medval = np.median(msbias.image)
        madval = 1.4826 * np.median(np.abs(medval - msbias.image))
        ww = np.where(np.abs(msbias.image - medval) > 10.0 * madval)
        bpm_img[ww] = 1

        # Return
        return bpm_img

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask.

        Currently identical to calling :func:`empty_bpm`.

        Even though they are both optional, either the precise shape for
        the image (`shape`) or an example file that can be read to get
        the shape (`filename` using :func:`get_image_shape`) *must* be
        provided.

        Args:
            filename (:obj:`str` or None):
                An example file to use to get the image shape.
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.
            shape (tuple, optional):
                Processed image shape
                Required if filename is None
                Ignored if filename is not None
            msbias (`numpy.ndarray`, optional):
                Master bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Generate an empty BPM first
        bpm_img = self.empty_bpm(filename, det, shape=shape)

        # Fill in bad pixels if a master bias frame is provided
        if msbias is not None:
            bpm_img = self.bpm_frombias(msbias, det, bpm_img)

        return bpm_img

    def get_slitmask(self, filename):
        """
        Empty for base class.  See derived classes.
        """
        return None

    def mask_to_pixel_coordinates(self, x=None, y=None, wave=None, order=1, filename=None,
                                  corners=False):
        """
        Returns an error message if `mask_to_pixel_coordinates` crashed because `use_maskdesign`
        is set to True for a spectrograph that does not support it.
        """
        msgs.error('This spectrograph does not support the use of mask design. Set `use_maskdesign=False`')

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
        core_meta = meta.define_core_meta()
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

    def get_detector_par(self, hdu, det):
        pass

    def get_rawimage(self, raw_file, det):
        """
        Load up the raw image and generate a few other bits and pieces
        that are key for image processing

        Parameters
        ----------
        raw_file : :obj:`str`
            File to read
        det : :obj:`int`
            Detector to read

        Returns
        -------
        detector_par : :class:`pypeit.par.pypeitpar.DetectorPar`
        raw_img : `numpy.ndarray`_
            Raw image for this detector
        hdu : `astropy.io.fits.HDUList`_
            Opened fits file
        exptime : :obj:`float`
        rawdatasec_img : `numpy.ndarray`_
        oscansec_img : `numpy.ndarray`_

        """
        # Open
        hdu = fits.open(raw_file)

        # Grab the DetectorPar
        detector = self.get_detector_par(hdu, det)

        # Raw image
        raw_img = hdu[detector['dataext']].data.astype(float)
        # TODO -- Move to FLAMINGOS2 spectrograph
        # Raw data from some spectrograph (i.e. FLAMINGOS2) have an
        # addition extention, so I add the following two lines. It's
        # easier to change here than writing another get_rawimage
        # function in the spectrograph file.
        if raw_img.ndim == 3:
            raw_img = raw_img[0]

        # Extras
        headarr = self.get_headarr(hdu)

        # Exposure time (used by ProcessRawImage)
        exptime = self.get_meta_value(headarr, 'exptime')

        # Rawdatasec, oscansec images
        binning = self.get_meta_value(headarr, 'binning')
        if detector['specaxis'] == 1:
            binning_raw = (',').join(binning.split(',')[::-1])
        else:
            binning_raw = binning

        for section in ['datasec', 'oscansec']:

            # Get the data section
            # Try using the image sections as header keywords
            # TODO -- Deal with user windowing of the CCD (e.g. Kast red)
            #  Code like the following maybe useful
            #hdr = hdu[detector[det - 1]['dataext']].header
            #image_sections = [hdr[key] for key in detector[det - 1][section]]
            # Grab from Detector
            image_sections = detector[section]
            #if not isinstance(image_sections, list):
            #    image_sections = [image_sections]
            # Always assume normal FITS header formatting
            one_indexed = True
            include_last = True

            # Initialize the image (0 means no amplifier)
            pix_img = np.zeros(raw_img.shape, dtype=int)
            for i in range(detector['numamplifiers']):

                if image_sections is not None:  # and image_sections[i] is not None:
                    # Convert the data section from a string to a slice
                    datasec = parse.sec2slice(image_sections[i], one_indexed=one_indexed,
                                              include_end=include_last, require_dim=2,
                                              binning=binning_raw)
                    # Assign the amplifier
                    pix_img[datasec] = i+1
            # Finish
            if section == 'datasec':
                rawdatasec_img = pix_img.copy()
            else:
                oscansec_img = pix_img.copy()

        # Return
        return detector, raw_img, hdu, exptime, rawdatasec_img, oscansec_img

    def get_lamps_status(self, headarr):
        """
        Return a string containing the information on the lamp status

        Args:
            headarr (list of fits headers):
              list of headers

        Returns:
            str: A string that uniquely represents the lamp status
        """
        # Loop through all lamps and collect their status
        kk = 1
        lampstat = []
        while True:
            lampkey = 'lampstat{:02d}'.format(kk)
            if lampkey not in self.meta.keys():
                break
            # Pull value from header
            lampstat += self.get_meta_value(headarr, lampkey)
            kk += 1
        return "_".join(lampstat)

    def get_meta_value(self, inp, meta_key, required=False, ignore_bad_header=False,
                       usr_row=None, no_fussing=False):
        """
        Return meta data from a given file (or its array of headers)

        Args:
            inp (str or list):
              Input filename or headarr list
            meta_key (str or list of str):
            headarr (list, optional)
              List of headers
            required (bool, optional):
              Require the meta key to be returnable
            ignore_bad_header: bool, optional
              Over-ride required;  not recommended
            usr_row: Row
              Provides user supplied frametype (and other things not used)
            no_fussing (bool, optional):
                No type checking or anything.  Just pass back the first value retrieved
                Mainly for bound pairs of meta, e.g. ra/dec

        Returns:
            value: value or list of values

        """
        headarr = self.get_headarr(inp) if isinstance(inp, str) else inp

        # Loop?
        if isinstance(meta_key, list):
            return [self.get_meta_value(headarr, mdict, required=required) for mdict in meta_key]

        # Are we prepared to provide this meta data?
        if meta_key not in self.meta.keys():
            if required:
                msgs.error("Need to allow for meta_key={} in your meta data".format(meta_key))
            else:
                msgs.warn("Requested meta data for meta_key={} does not exist...".format(meta_key))
                return None

        # Check if this meta key is required
        if 'required' in self.meta[meta_key].keys():
            required = self.meta[meta_key]['required']

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

        # Return now?
        if no_fussing:
            return value

        # Deal with 'special' cases
        if meta_key in ['ra', 'dec'] and value is not None:
            # TODO: Can we get rid of the try/except here and instead get to the heart of the issue?
            try:
                ra, dec = meta.convert_radec(self.get_meta_value(headarr, 'ra', no_fussing=True),
                                    self.get_meta_value(headarr, 'dec', no_fussing=True))
            except:
                msgs.warn('Encounter invalid value of your coordinates. Give zeros for both RA and DEC')
                ra, dec = 0.0, 0.0
            value = ra if meta_key == 'ra' else dec

        # JFH Added this bit of code to deal with situations where the
        # header card is there but the wrong type, e.g. MJD-OBS =
        # 'null'
        try:
            if self.meta_data_model[meta_key]['dtype'] == str:
                retvalue = str(value).strip()
            elif self.meta_data_model[meta_key]['dtype'] == int:
                retvalue = int(value)
            elif self.meta_data_model[meta_key]['dtype'] == float:
                retvalue = float(value)
            elif self.meta_data_model[meta_key]['dtype'] == tuple:
                if not isinstance(value, tuple):
                    msgs.error('dtype for {0} is tuple, but value '.format(meta_key)
                               + 'provided is {0}.  Casting is not possible.'.format(type(value)))
                retvalue = value
            castable = True
        except:
            retvalue = None
            castable = False

        # JFH Added the typing to prevent a crash below when the header
        # value exists, but is the wrong type. This causes a crash
        # below when the value is cast.
        if value is None or not castable:
            # Was this required?
            if required:
                kerror = True
                if not ignore_bad_header:
                    # Is this meta required for this frame type (Spectrograph specific)
                    if ('required_ftypes' in self.meta[meta_key]) and (usr_row is not None):
                        kerror = False
                        # Is it required?
                        # TODO: Use numpy.isin ?
                        for ftype in usr_row['frametype'].split(','):
                            if ftype in self.meta[meta_key]['required_ftypes']:
                                kerror = True
                    # Bomb out?
                    if kerror:
                        embed(header=utils.embed_header())
                        msgs.error('Required meta "{0}" did not load!'.format(meta_key)
                                   + 'You may have a corrupt header.')
                else:
                    msgs.warn('Required card {0} missing '.format(self.meta[meta_key]['card'])
                              + 'from your header.  Proceeding with risk...')
            return None

        # Return
        return retvalue

    def get_wcs(self, hdr, slits, platescale, wave0, dwv):
        """Get the WCS for a frame

        Parameters
        ----------
        hdr : fits header
            The header of the raw frame. The information in this
            header will be extracted and returned as a WCS.
        slits : :class:`pypeit.slittrace.SlitTraceSet`
            Master slit edges
        platescale : float
            platescale of an unbinned pixel in arcsec/pixel (e.g. detector.platescale)
        wave0 : float
            wavelength zeropoint
        dwv : float
            delta wavelength per spectral pixel

        Returns
        -------
        astropy.wcs : An astropy WCS object.
        """
        msgs.warn("No WCS setup for spectrograph: {0:s}".format(self.spectrograph))
        return None

    def get_datacube_bins(self, slitlength, minmax, num_wave):
        """Calculate the bin edges to be used when making a datacube

        Parameters
        ----------
        slitlength : int
            Length of the slit in pixels
        minmax : `numpy.ndarray`_
            An array of size (nslits, 2), listing the minimum and maximum pixel
            locations on each slit relative to the reference location (usually
            the centre of the slit). This array is returned by the function
            `slittrace.SlitTraceSet.get_radec_image`_
        num_wave : int
            Number of wavelength steps = int(round((wavemax-wavemin)/delta_wave))

        Returns
        -------
        tuple : Three 1D numpy.ndarray providing the bins to use when constructing a histogram
                of the spec2d files. The elements are (x, y, lambda).
        """
        msgs.warn("No datacube setup for spectrograph: {0:s}".format(self.spectrograph))
        return None

    def validate_metadata(self):
        """
        Validates the meta definitions of the Spectrograph
        by making a series of comparisons to the meta data model
        definied in metadata.py
        """
        # Load up
        # TODO: Can we indicate if the metadata element is core instead
        # of having to call both of these?
        core_meta = meta.define_core_meta()
        # KBW: These should have already been defined to self
        #meta_data_model = meta.get_meta_data_model()

        # Check core
        core_keys = np.array(list(core_meta.keys()))
        indx = np.invert(np.isin(core_keys, list(self.meta.keys())))
        if np.any(indx):
            msgs.error('Required keys {0} not defined by spectrograph!'.format(core_keys[indx]))

        # Check for rtol for config keys that are type float
        config_keys = np.array(self.configuration_keys())
        indx = ['rtol' not in self.meta[key].keys() if self.meta_data_model[key]['dtype'] == float 
                    else False for key in config_keys]
        if np.any(indx):
            msgs.error('rtol not set for {0} keys in spectrograph meta!'.format(config_keys[indx]))

        # Now confirm all meta are in the data model
        meta_keys = np.array(list(self.meta.keys()))
        indx = np.invert(np.isin(meta_keys, list(self.meta_data_model.keys())))
        if np.any(indx):
            msgs.error('Meta data keys {0} not in metadata model'.format(meta_keys[indx]))

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
                    return ['None']*999 # self.numhead
        else:
            hdu = inp
        return [hdu[k].header for k in range(len(hdu))]

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
    def pypeline(self):
        return 'MultiSlit'

#    JXP says -- LEAVE THIS HERE FOR NOW. WE MAY NEED IT
#    def mm_per_pix(self, det=1):
#        """
#        Return the spatial scale at the telescope focal plane in mm per
#        pixel at the detector.
#
#        The fratio and diameter of the telescope must be defined.
#
#        Args:
#            det (:obj:`int`, optional):
#                Detector to use for the spectrograph platescale.
#
#        Returns:
#            float: The spatial scale at the telescope focal plane in mm
#            per detector pixel scale.
#
#        Raises:
#            ValueError:
#                Raised if the telescope is undefined, any of the numbers
#                needed for the calculation are not available, or the
#                selected detector is out of range.
#        """
#        if det > self.ndet:
#            raise ValueError('Selected detector out of range; det={0}..{1}.'.format(1,self.ndet))
#        tel_platescale = None if self.telescope is None else self.telescope.platescale()
#        if self.telescope is None or tel_platescale is None or \
#                self.detector[det-1]['platescale'] is None:
#            raise ValueError('Incomplete information to calculate mm per pixel.')
#
#        return self.detector[det-1]['platescale']/tel_platescale

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

    @property
    def norders(self):
        return None

    def check_disperser(self):
        """
        Ensure that the disperser is defined.
        """
        if self.dispname is None:
            msgs.error('Disperser used for observations is required.  Reinit with an example '
                       'science frame.')

    @property
    def order_spat_pos(self):
        return None

    @property
    def orders(self):
        return None

    @property
    def spec_min_max(self):
        return None

    @property
    def dloglam(self):
        return None

    @property
    def loglam_minmax(self):
        return None

    # TODO : This code needs serious work.  e.g. eliminate the try/except
    def slit_minmax(self, slit_spat_pos, binspectral=1):
        """

        Args:
            slit_spat_pos (float or ndarray):
                normalized slit_spatial position as computed by edgetrace.slit_spat_pos
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

