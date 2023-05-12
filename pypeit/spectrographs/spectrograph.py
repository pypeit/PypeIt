"""
Defines the abstract :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
class, which is the parent class for all instruments served by PypeIt.

The key functionality of this base class and its derived classes are to
provide instrument-specific:

    - file I/O routines
    - detector properties (see
      :class:`~pypeit.images.detector_container.DetectorContainer`)
    - telescope properties (see :class:`~pypeit.par.pypeitpar.TelescopePar`)
    - fits header keywords that are collated and injested into a metadata
      table that it uses throughout the reduction (see
      :class:`~pypeit.metadata.PypeItMetaData`)
    - header keyword values to check to confirm a fits file has been taken
      with the selected instrument
    - default methods for automatically determining the type of each exposure
      that PypeIt was asked to reduce
    - header keywords to use when matching calibration frames to science
      frames
    - methods used to generate and/or read bad-pixel masks for an exposure
    - default parameters for the reduction algorithms
    - methods to access an archival sky spectrum

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from abc import ABCMeta
from pathlib import Path

from IPython import embed

import numpy as np
from astropy.io import fits

from pypeit import msgs
from pypeit import io
from pypeit.core import parse
from pypeit.core import procimg
from pypeit.core import meta
from pypeit.par import pypeitpar
from pypeit.images.detector_container import DetectorContainer
from pypeit.images.mosaic import Mosaic


# TODO: Create an EchelleSpectrograph derived class that holds all of
# the echelle specific methods.

class Spectrograph:
    """
    Abstract base class for all instrument-specific behavior in PypeIt.

    Attributes:
        dispname (:obj:`str`):
            Name of the dispersing element.
        rawdatasec_img (`numpy.ndarray`_):
            An image identifying the amplifier that reads each detector
            pixel.
        oscansec_img (`numpy.ndarray`_):
            An image identifying the amplifier that reads each detector
            pixel
        slitmask (:class:`~pypeit.spectrographs.slitmask.SlitMask`):
            Provides slit and object coordinate data for an
            observation. Not necessarily populated for all
            spectrograph instantiations.
        primary_hdrext (:obj:`int`):
            0-indexed number of the extension in the raw frames with the
            primary header data.
        meta (:obj:`dict`):
            Instrument-specific metadata model, linking header information to
            metadata elements required by PypeIt.
    """
    __metaclass__ = ABCMeta

    ndet = None
    """
    Number of detectors for this instrument.
    """

    name = None
    """
    The name of the spectrograph. See :ref:`instruments` for the currently
    supported spectrographs.
    """

    telescope = None
    """
    Instance of :class:`~pypeit.par.pypeitpar.TelescopePar` providing
    telescope-specific metadata.
    """

    camera = None
    """
    Name of the spectrograph camera or arm.
    This is used by specdb, so use that naming convention
    """

    url = None
    """
    Reference url
    """

    header_name = None
    """
    Name of the spectrograph camera or arm from the Header.
    Usually the INSTRUME card.
    """

    pypeline = 'MultiSlit'
    """
    String used to select the general pipeline approach for this
    spectrograph.
    """

    ech_fixed_format = None
    """
    If an echelle spectrograph, this will be set to a boolean indicating whether it is a fixed format or tiltable 
    echelle. 
    """

    supported = False
    """
    Flag that PypeIt code base has been sufficiently tested with data
    from this spectrograph that it is officially supported by the development
    team.
    """

    ql_supported = False
    """
    Flag that PypeIt code base has been sufficiently tested with data
    from this spectrograph in quicklook mode
    that it is officially supported by the development team.
    """


    comment = None
    """
    A brief comment or description regarding PypeIt usage with this
    spectrograph.
    """

    meta_data_model = meta.get_meta_data_model()
    """
    Metadata model that is generic to all spectrographs.
    """

    allowed_extensions = None
    """
    Defines the allowed extensions for the input fits files.
    """

    def __init__(self):
        self.dispname = None
        self.rawdatasec_img = None
        self.oscansec_img = None
        self.slitmask = None

        # Extension with the primary header data
        self.primary_hdrext = 0

        # Generate and check the instrument-specific metadata definition
        self.init_meta()
        self.validate_metadata()
        if self.pypeline == 'Echelle' and self.ech_fixed_format is None:
            msgs.error('ech_fixed_format must be set for echelle spectrographs')

        # TODO: Is there a better way to do this?
        # Validate the instance by checking that the class has defined the
        # number of detectors
        assert self.ndet > 0

        # TODO: Add a call to _check_telescope here?

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = cls.name
        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            scifile (:obj:`str`):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        return self.__class__.default_pypeit_par() if inp_par is None else inp_par

    def update_edgetracepar(self, par):
        """
        This method is used in :func:`pypeit.edgetrace.EdgeTraceSet.maskdesign_matching`
        to update EdgeTraceSet parameters when the slitmask design matching is not feasible
        because too few slits are present in the detector.

        **This method is not defined for all spectrographs.**

        Args:
            par (:class:`pypeit.par.pypeitpar.EdgeTracePar`):
                The parameters used to guide slit tracing.

        Returns:
            :class:`pypeit.par.pypeitpar.EdgeTracePar`
            The modified parameters used to guide slit tracing.
        """

        return par

    # TODO: Include a script in doc/ that will automatically document the
    # parameters set here.
    @staticmethod
    def ql_par():
        """
        Return the list of parameters specific to quick-look reduction of
        science frames for this spectrograph.

        The following parameters are set by default for any spectrograph.
        Derived classes may override this function.

        Returns:
            :obj:`dict`: Dictionary with a form that matches
            :class:`~pypeit.par.pypeitpar.PypeItPar` with the parameters to be
            set.
        """
        return dict(
            # Calibration + misc
            rdx = dict(
                ignore_bad_headers = True
            ),
            baseprocess = dict(
                use_biasimage = False
            ),
            calibrations = dict(
                # This allows for science frames only and "cooked" calibrations
                raise_chk_error = False,
                flatfield = dict(
                    saturated_slits = 'mask'
                )
            ),
            scienceframe = dict(
                process = dict(
                    mask_cr = False
                )
            ),
            # Reduction parameters
            reduce = dict(
                extraction = dict(
                    skip_optimal = True
                ),
                findobj = dict(
                    skip_second_find = True
                )
            )
        )

    def _check_extensions(self, filename):
        """
        Check if this filename has an allowed extension

        Args:
            filename (:obj:`str`, `Path`_):
                Input raw fits filename
        """
        if self.allowed_extensions is not None:
            _filename = Path(filename).resolve()
            if _filename.suffix not in self.allowed_extensions:
                msgs.error(f'The input file ({_filename.name}) does not have a recognized '
                           f'extension ({_filename.suffix}).  The allowed extensions for '
                           f'{self.name} include {",".join(self.allowed_extensions)}.')

    def _check_telescope(self):
        """Check the derived class has properly defined the telescope."""
        if self.telescope is None:
            raise ValueError('Must define the telescope used to take the observations.')
        if not isinstance(self.telescope, pypeitpar.TelescopePar):
                raise TypeError('Telescope parameters must be one of those specified in'
                                'pypeit.telescopes.')

    def raw_is_transposed(self, detector_par):
        """
        Check if raw image files are transposed with respect to the
        PypeIt convention.

        Indicates that reading raw files with `astropy.io.fits`_ yields an
        image with the spatial dimension along the first axis of the 2D
        array. This means that the image must be transposed to match the
        PypeIt convention of having the spectral dimension along the
        first axis.

        Args:
            detector_par (:class:`~pypeit.images.detector_container.DetectorContainer`):
                Detector-specific metadata.

        Returns:
            :obj:`bool`: Flag that transpose is required.
        """
        return detector_par['specaxis'] == 1

    # TODO: This circumvents all the infrastructure we have for pulling
    # metadata from headers. Why aren't we using self.meta and
    # self.get_meta_value? See pypeit.metadata.PypeItMetaData._build()
    def parse_spec_header(self, header):
        """
        Parses an input header for key spectrograph items.

        Args:
            header (`astropy.io.fits.Header`_):
                Fits header read from a file.

        Returns:
            :obj:`dict`: Dictionary with the metadata read from ``header``.
        """
        spec_dict = {}
        # The keys in spec_dict should be the CORE metadata,
        #   spectrograph CONFIGURATION KEYS, and the FILENAME
        core_meta_keys = list(meta.define_core_meta().keys())
        core_meta_keys += self.configuration_keys()
        core_meta_keys += ['filename']
        for key in core_meta_keys:
            if key.upper() in header.keys():
                spec_dict[key.upper()] = header[key.upper()]
        # Return
        return spec_dict

    def subheader_for_spec(self, row_fitstbl, raw_header, extra_header_cards=None,
                           allow_missing=False):
        """
        Generate a dict that will be added to the Header of spectra files
        generated by PypeIt (e.g. :class:`~pypeit.specobjs.SpecObjs`).

        Args:
            row_fitstbl (dict-like):
                Typically an `astropy.table.Row`_ or
                `astropy.io.fits.Header`_ with keys defined by
                :func:`~pypeit.core.meta.define_core_meta`.
            raw_header (`astropy.io.fits.Header`_):
                Header that defines the instrument and detector, meaning that
                the header must contain the ``INSTRUME`` and ``DETECTOR``
                header cards. If provided, this must also contain the header
                cards provided by ``extra_header_cards``.
            extra_header_cards (:obj:`list`, optional):
                Additional header cards from ``raw_header`` to include in the
                output dictionary. Can be an empty list or None.
            allow_missing (:obj:`bool`, optional):
                Ignore any keywords returned by
                :func:`~pypeit.core.meta.define_core_meta` are not present in
                ``row_fitstbl``. Otherwise, raise ``PypeItError``.

        Returns:
            :obj:`dict`: Dictionary with data to include an output fits
            header file or table downstream.
        """
        subheader = {}

        core_meta = meta.define_core_meta()
        # Core Metadata Keys -- These must be present
        for key in core_meta.keys():
            try:
                subheader[key] = (row_fitstbl[key], core_meta[key]['comment'])
            except KeyError:
                if not allow_missing:
                    msgs.error(f"Core Meta Key: {key} not present in your fitstbl/Header")
        # Configuration Keys -- In addition to Core Meta,
        #   other Config-Specific values; optional
        for key in self.configuration_keys():
            if key not in subheader:
                try:
                    subheader[key] = row_fitstbl[key]
                except KeyError:
                    # If configuration_key is not in row_fitstbl, warn but move on
                    msgs.warn(f"Configuration Key: {key} not present in your fitstbl/Header")
        # Add a few more
        for key in ['filename']:  # For fluxing
            subheader[key] = row_fitstbl[key]

        # The following are pulled from the original header, if available
        header_cards = ['INSTRUME', 'DETECTOR', 'DATE-OBS'] + self.raw_header_cards()
        if extra_header_cards is not None:
            header_cards += extra_header_cards  # For specDB and more
        for card in header_cards:
             if card in raw_header.keys():
                 subheader[card] = raw_header[card]  # Self-assigned instrument name

        # Specify which pipeline created this file
        subheader['PYPELINE'] = self.pypeline
        subheader['PYP_SPEC'] = (self.name, 'PypeIt: Spectrograph name')

        # Observatory and Header supplied Instrument
        subheader['TELESCOP'] = (self.telescope['name'], 'Telescope')
        subheader['LON-OBS'] = (self.telescope['longitude'], 'Telescope longitude')
        subheader['LAT-OBS'] = (self.telescope['latitude'], 'Telescope latitute')
        subheader['ALT-OBS'] = (self.telescope['elevation'], 'Telescope elevation')

        # Return
        return subheader

    def orient_image(self, detector_par, rawimage):
        """
        Orient the image into the PypeIt configuration: (spectral,
        spatial).

        Args:
            detector_par (:class:`pypeit.images.detector_container.DetectorContainer`):
                Detector metadata.
            rawimage (`numpy.ndarray`_):
                Image from the raw frame

        Returns:
            `numpy.ndarray`_: Re-oriented image.
        """
        image = rawimage.copy()
        # Transpose?
        if self.raw_is_transposed(detector_par):
            image = image.T
        # Flip spectral axis?
        if detector_par['specflip']:
            image = np.flip(image, axis=0)
        # Flip spatial axis?
        if detector_par['spatflip']:
            image = np.flip(image, axis=1)
        return image

    # TODO: JFH Are these bad pixel masks in the raw frame, or the
    # flipped/transposed pypeit frame?? KBW: Does the new description of
    # "shape" answer this? (JXP please check I edited this correctly).
    def empty_bpm(self, filename, det, shape=None):
        """
        Generate a generic (empty) bad-pixel mask.

        Even though they are both optional, either the precise shape for the
        image (``shape``) or an example file that can be read to get the
        shape (``filename``) *must* be provided. In the latter, the file is
        read, trimmed, and re-oriented to get the output shape. If both
        ``shape`` and ``filename`` are provided, ``shape`` is ignored.

        This is the generic function provided in the base class meaning that all
        pixels are returned as being valid/unmasked.

        Args:
            filename (:obj:`str`):
                An example file to use to get the image shape. Can be None,
                but ``shape`` must be provided, if so. Note the overhead of
                this function is large if you ``filename``. You're better off
                providing ``shape``, if you know it.
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.
            shape (:obj:`tuple`, optional):
                Processed image shape. I.e., if the image for this instrument
                is re-oriented, trimmed, etc, this shape must be that of the
                re-oriented (trimmed, etc) image. This is required if
                ``filename`` is None, but ignored otherwise.

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set to 1 and
            an unmasked value set to 0. The shape of the returned image should
            be that of a trimmed and oriented PypeIt processed image. This
            function specifically is the generic method for the base class,
            meaning that all pixels are returned as unmasked (0s).
        """
        # TODO: I think shape should take precedence over filename, not the
        # other way around.  I.e., if we know the shape we want going in, why
        # are we still defaulting to reading, trimming, and re-orienting an
        # input image just to figure out that shape?

        # Load the raw frame
        if filename is None:
            _shape = shape
        else:
            detector_par, _,  _, _, rawdatasec_img, _ = self.get_rawimage(filename, det)
            # Trim + reorient
            trim = procimg.trim_frame(rawdatasec_img, rawdatasec_img < 1)
            orient = self.orient_image(detector_par, trim)#, det)
            _shape = orient.shape

        # Shape must be defined at this point.
        if _shape is None:
            msgs.error('Must specify shape if filename is None.')

        # Generate
        # TODO: Why isn't this a boolean array?
        return np.zeros(_shape, dtype=np.int8)

    def bpm_frombias(self, msbias, bpm_img, thresh=10.):
        """
        Generate a bad-pixel mask from a processed bias frame.

        Args:
            msbias (:class:`~pypeit.images.pypeitimage.PypeItImage`):
                Processed bias frame used to identify bad pixels.
            bpm_img (`numpy.ndarray`_):
                Zeroth-order bad pixel mask; i.e., generated using
                :func:`~pypeit.spectrographs.spectrograph.Spectrograph.empty_bpm`.
                **Must** be the same shape as ``msbias``.
            thresh (:obj:`float`, optional):
                The sigma threshold used to identify bad pixels.

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set to 1
            and an unmasked value set to 0. The shape of the returned image
            is the same as the provided ``msbias`` and ``bpm_img`` images.
        """
        # Check that the bias has the correct shape
        if msbias.image.shape != bpm_img.shape:
            msgs.error(f'Shape mismatch between processed bias {msbias.image.shape} and expected '
                       f'BPM {bpm_img.shape}.')
        # Setup
        nimg = 1 if bpm_img.ndim == 2 else bpm_img.shape[0]

        # NOTE: expand_dims does *not* copy the array, so we need to do so
        # explicitly here ...
        _bpm_img = np.expand_dims(bpm_img.copy(), 0) if nimg == 1 else bpm_img.copy()
        # ... but not here.  I.e., we don't want to change the input bpm_img,
        # and we only ever access the values of msbias.
        _bias = np.expand_dims(msbias.image, 0) if nimg == 1 else msbias.image

        # Treat each bias image separately
        for i in range(nimg):
            medval = np.median(_bias[i])
            madval = 1.4826 * np.median(np.absolute(_bias[i] - medval))
            _bpm_img[i,np.abs(_bias[i] - medval) > thresh * madval] = 1
        # Done
        return _bpm_img[0] if nimg == 1 else _bpm_img

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask.

        Even though they are both optional, either the precise shape for
        the image (``shape``) or an example file that can be read to get
        the shape (``filename`` using :func:`get_image_shape`) *must* be
        provided.

        Args:
            filename (:obj:`str`):
                An example file to use to get the image shape.  Can be None.
            det (:obj:`int`, :obj:`tuple`):
                1-indexed detector(s) to read.  An image mosaic is selected
                using a :obj:`tuple` with the detectors in the mosaic, which
                must be one of the allowed mosaics returned by
                :func:`allowed_mosaics`.
            shape (:obj:`tuple`, optional):
                Processed image shape.  If ``filename`` is None, this *must* be
                provided; otherwise, this is ignored.
            msbias (:class:`~pypeit.images.pypeitimage.PypeItImage`, optional):
                Processed bias frame.  If provided, it is used by
                :func:`~pypeit.spectrographs.spectrograph.Spectrograph.bpm_frombias`
                to identify bad pixels.

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set to 1 and
            an unmasked value set to 0.
        """
        # Validate the entered (list of) detector(s)
        nimg, _det = self.validate_det(det)
        # If using a mosaic, the shape of *all* processed images must be
        # identical.  Therefore, we only need to generate the empty BPM for one
        # of the detectors, like so:
        bpm_img = self.empty_bpm(filename, _det[0], shape=shape)
        # Repeat it if necessary
        if nimg > 1:
            bpm_img = np.tile(bpm_img, (nimg,1,1))

        if msbias is None:
            # Not using bias to identify bad pixels, so we're done
            return bpm_img

        msgs.info(f'Generating a BPM using bias for det={_det} for {self.name}')
        return self.bpm_frombias(msbias, bpm_img)

    def list_detectors(self, mosaic=False):
        """
        List the *names* of the detectors in this spectrograph.

        This is primarily used :func:`~pypeit.slittrace.average_maskdef_offset`
        to measure the mean offset between the measured and expected slit
        locations.  **This method is not defined for all spectrographs.**

        Detectors separated along the dispersion direction should be ordered
        along the first axis of the returned array.  For example, Keck/DEIMOS
        returns:

        .. code-block:: python

            dets = np.array([['DET01', 'DET02', 'DET03', 'DET04'],
                             ['DET05', 'DET06', 'DET07', 'DET08']])

        such that all the bluest detectors are in ``dets[0]``, and the slits
        found in detectors 1 and 5 are just from the blue and red counterparts
        of the same slit.

        Args:
            mosaic (:obj:`bool`, optional):
                Is this a mosaic reduction?
                It is used to determine how to list the detector, i.e., 'DET' or 'MSC'.

        Returns:
            `numpy.ndarray`_: The list of detectors in a `numpy.ndarray`_.  If
            the array is 2D, there are detectors separated along the dispersion
            axis.
        """
        return None

    def get_lamps(self, fitstbl):
        """
        Extract the list of arc lamps used from header.

        This method is not defined for all spectrographs. This base-class
        method raises an exception.
        """
        msgs.error('This spectrograph does not support the use of lamps list from header. '
                   'provide a list of lamps using the parameter `lamps` in WavelengthSolutionPar')

    def get_slitmask(self, filename):
        """
        Empty for base class.  See derived classes.
        """
        return None

    def mask_to_pixel_coordinates(self, x=None, y=None, wave=None, order=1, filename=None,
                                  corners=False):
        """
        Predict detector pixel coordinates for a given set of slit-mask
        coordinates.

        This method is not defined for all spectrographs. This base-class
        method raises an exception. This may be because ``use_maskdesign``
        has been set to True for a spectrograph that does not support it.
        """
        msgs.error('This spectrograph does not support the use of mask design. '
                   'Set `use_maskdesign=False`')

    def get_maskdef_slitedges(self, ccdnum=None, filename=None, debug=None, 
                              binning:str=None, trc_path:str=None):
        """
        Provides the slit edges positions predicted by the slitmask design.

        This method is not defined for all spectrographs. This base-class
        method raises an exception. This may be because ``use_maskdesign``
        has been set to True for a spectrograph that does not support it.

        Args:
            trc_path(str, optional): Path to the first trace file used to generate the trace flat
            binning(str, optional): spec,spat binning of the flat field image
            filename (:obj:`str` or :obj:`list`, optional): Name of the file holding the mask design info
                or the maskfile and wcs_file in that order
            debug (:obj:`bool`, optional): Debug
            ccdnum (:obj:`int`, optional): detector number
        """
        msgs.error('This spectrograph does not support the use of mask design. '
                   'Set `use_maskdesign=False`')

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        This list is used by :class:`~pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            :obj:`list`: List of keywords of data pulled from file headers
            and used to constuct the :class:`~pypeit.metadata.PypeItMetaData`
            object.
        """
        return ['dispname', 'dichroic', 'decker']

    def same_configuration(self, configs, check_keys=True):
        """
        Check if a set of instrument setup/configurations are all the same,
        within the tolerance set for this spectrograph for each configuration
        key.

        Args:
            configs (:obj:`dict`, array-like):
                A list or dictionary of configuration parameters.  Each list or
                dictionary element must be a dictionary with the values for the
                configuration-defining metadata parameters for this
                spectrograph.
            check_keys (:obj:`bool`, optional):
                Check that the keys in each configuration match the expected
                keys defined by :func:`configuration_keys`.  If False, the keys
                are set by the first configuration in ``configs``, *not* those
                returned by :func:`configuration_keys`.
        
        Returns:
            :obj:`bool`: Flag that all configurations in the list are the same.
        """
        cfg_id = list(configs.keys() if isinstance(configs, dict) else range(len(configs)))

        if check_keys:
            # The list of metadata keys used to define a unique configuration
            cfg_meta = self.configuration_keys()
            # Check that the relevant keys are in the first configuration
            for key in cfg_meta:
                if key not in configs[cfg_id[0]].keys():
                    msgs.error(f'Configuration {cfg_id[0]} missing required key, {key}.  Cannot '
                               'determine if configurations are the same!')
                if key not in self.meta.keys():
                    msgs.error(f'CODING ERROR: {key} is a configuration key but not defined in '
                               f'the metadata dictionary for {self.__class__.__name__}!')
        else:
            cfg_meta = configs[cfg_id[0]].keys()

        # Match against all of the other configurations
        for _cfg_id in cfg_id[1:]:
            matched = []
            for key in cfg_meta:
                if key not in configs[_cfg_id].keys():
                    msgs.error(f'Configuration {_cfg_id} missing required key, {key}.  Cannot '
                               'determine if configurations are the same!')
                # TODO: Instead check if 'rtol' exists and is not None?
                if isinstance(configs[cfg_id[0]][key], (float, np.floating)) \
                        and isinstance(configs[_cfg_id][key], (float, np.floating)):
                    # NOTE: No float-valued metadata can be 0!
                    matched += [np.abs(configs[cfg_id[0]][key]-configs[_cfg_id][key])
                                    / configs[cfg_id[0]][key] < self.meta[key]['rtol']]
                else:
                    matched += [np.all(configs[cfg_id[0]][key] == configs[_cfg_id][key])]
            if not np.all(matched):
                # We found a difference so return
                return False
        # Went through all configurations and didn't find any differences
        return True

    def raw_header_cards(self):
        """
        Return additional raw header cards to be propagated in
        downstream output files for configuration identification.

        The list of raw data FITS keywords should be those used to populate
        the :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.configuration_keys`
        or are used in :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.config_specific_par`
        for a particular spectrograph, if different from the name of the
        PypeIt metadata keyword.

        This list is used by :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.subheader_for_spec`
        to include additional FITS keywords in downstream output files.

        Returns:
            :obj:`list`: List of keywords from the raw data files that should
            be propagated in output files.
        """
        return []

    def modify_config(self, row, cfg):
        """
        Modify the configuration dictionary for a given frame. This method is
        used in :func:`~pypeit.metadata.PypeItMetaData.set_configurations` to
        modify in place the configuration requirement to assign a specific frame
        to the current setup.

        **This method is not defined for all spectrographs.**

        Args:
            row (`astropy.table.Row`_):
                The table row with the metadata for one frame.
            cfg (:obj:`dict`):
                Dictionary with metadata associated to a specific configuration.

        Returns:
            :obj:`dict`: modified dictionary with metadata associated to a
            specific configuration.
        """
        return cfg

    def valid_configuration_values(self):
        """
        Return a fixed set of valid values for any/all of the configuration
        keys.

        Method is undefined for the base class.

        Returns:
            :obj:`dict`: A dictionary with any/all of the configuration keys
            and their associated discrete set of valid values. If there are
            no restrictions on configuration values, None is returned.
        """
        pass

    # TODO: This feels like something that should be in the PypeItMetaData
    # class, not the spectrograph class.
    def vet_instrument(self, meta_tbl):
        """
        Confirm the metadata gathered for a set of measurements are all unique
        and from this spectrograph, according to the expected instrument name in
        the headers of its raw data files.

        This function *only* issues warnings; no exceptions are raised.

        Args:
            meta_tbl (`astropy.table.Table`_):
                Table with the meta data; see
                :class:`~pypeit.metadata.PypeItMetaData`.
        """
        if 'instrument' in meta_tbl.keys():
            if self.header_name is None:
                msgs.error('CODING ERROR: header_name is not defined for '
                           f'{self.__class__.__name__}!')
            # Check that there is only one instrument
            #  This could fail if one mixes is much older calibs
            indx = meta_tbl['instrument'].data != None
            instr_names = np.unique(meta_tbl['instrument'].data[indx])
            if len(instr_names) != 1:
                msgs.warn(f'More than one instrument in your dataset! {instr_names} \n'
                          'Proceed with great caution...')
            # Check the name
            if instr_names[0] != self.header_name:
                msgs.warn('The instrument name in the headers of the raw files does not match the '
                          f'expected one! Found {instr_names[0]}, expected {self.header_name}.  '
                          'You may have chosen the wrong PypeIt spectrograph name!')


    def config_independent_frames(self):
        """
        Define frame types that are independent of the fully defined
        instrument configuration.

        By default, bias and dark frames are considered independent of a
        configuration; however, at the moment, these frames can only be
        associated with a *single* configuration. That is, you cannot take
        afternoon biases, change the instrument configuration during the
        night, and then use the same biases for both configurations. See
        :func:`~pypeit.metadata.PypeItMetaData.set_configurations`.

        This method returns a dictionary where the keys of the dictionary are
        the list of configuration-independent frame types. The value of each
        dictionary element can be set to one or more metadata keys that can
        be used to assign each frame type to a given configuration group. See
        :func:`~pypeit.metadata.PypeItMetaData.set_configurations` and how it
        interprets the dictionary values, which can be None.

        Returns:
            :obj:`dict`: Dictionary where the keys are the frame types that
            are configuration-independent and the values are the metadata
            keywords that can be used to assign the frames to a configuration
            group.
        """
        return {'bias': None, 'dark': None}

    def get_comb_group(self, fitstbl):
        """
        Automatically assign combination groups and background images by parsing
        known dither patterns.

        This method is used in
        :func:`~pypeit.metadata.PypeItMetaData.set_combination_groups`, and
        directly modifies the ``comb_id`` and ``bkg_id`` columns in the provided
        table.

        **This method is not defined for all spectrographs.**  This base-class
        implementation simply returns the input table.

        Args:
            fitstbl(`astropy.table.Table`_):
                The table with the metadata for all the frames.

        Returns:
            `astropy.table.Table`_: modified fitstbl.
        """
        return fitstbl

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :func:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
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
        Methods to generate metadata requiring interpretation of the header
        data, instead of simply reading the value of a header card.

        Method is undefined in this base class.

        Args:
            headarr (:obj:`list`):
                List of `astropy.io.fits.Header`_ objects.
            meta_key (:obj:`str`):
                Metadata keyword to construct.

        Returns:
            object: Metadata value read from the header(s).
        """
        pass

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}

    def meta_key_map(self):
        """
        Print the mapping of the pypeit-specific metadata keywords to the
        header cards used for this spectrograph.

        .. note::
            Metadata keys with header cards that are None have no simple
            mapping between keyword and header card; their values are set by
            some combination of header keywords as defined by
            :func:`compound_meta` specific to each spectrograph.
        """
        meta_keys = list(self.meta.keys())
        meta_cards = [str(self.meta[key]['card']) for key in meta_keys]
        nk = max(12, max([len(key) for key in meta_keys]))
        nc = max(11, max([len(card) for card in meta_cards]))
        print('')
        print('{0}   {1}'.format('Metadata Key'.center(nk), 'Header Card'.center(nc)))
        print('-'*nk + '   ' + '-'*nc)
        for key, card in zip(meta_keys, meta_cards):
            print('{0}   {1}'.format(key.rjust(nk), card.rjust(nc)))
        print('')

    def get_detector_par(self, det, hdu=None):
        """
        Read/Set the detector metadata.

        This method is needed by some instruments that require the detector
        metadata to be interpreted from the output files. This method is
        undefined in the base class.
        """
        pass

    def get_mosaic_par(self):
        """
        Return the hard-coded parameters needed to construct detector mosaics.
        """
        pass

    @property
    def allowed_mosaics(self):
        """
        The list of allowed detector mosaics.

        For instruments with no allowed detector mosaics, this *must* be
        returned as an empty list.
        """
        return []

    def get_det_name(self, det):
        """
        Return a name for the detector or mosaic.

        This is a simple wrapper for
        :func:`pypeit.images.detector_container.DetectorContainer.get_name` or
        :func:`pypeit.images.mosaic.Mosaic.get_name`, depending on the type of
        ``det``.

        Args:
            det (:obj:`int`, :obj:`tuple`):
                The 1-indexed detector number(s).  If a tuple, it must include
                detectors designated as a viable mosaic for
                :attr:`spectrograph`; see
                :func:`~pypeit.spectrographs.spectrograph.Spectrograph.allowed_mosaics`.

        Returns:
            :obj:`str`: Name for the detector or mosaic.
        """
        if isinstance(det, tuple):
            # The "detector" is a mosaic.
            if det not in self.allowed_mosaics:
                msgs.error(f'{det} is not an allowed mosaic for {self.name}.')
            return Mosaic.get_name(self.allowed_mosaics.index(det)+1)

        # Single detector
        if det <= 0 or det > self.ndet:
            msgs.error(f'{det} is not a valid detector for {self.name}.')
        return DetectorContainer.get_name(det)

    def get_det_id(self, det):
        """
        Return an identifying index for the detector or mosaic.

        Args:
            det (:obj:`int`, :obj:`tuple`):
                The 1-indexed detector number(s).  If a tuple, it must include
                detectors designated as a viable mosaic for
                :attr:`spectrograph`; see
                :func:`~pypeit.spectrographs.spectrograph.Spectrograph.allowed_mosaics`.

        Returns:
            :obj:`int`: Unique index for the detector or mosaic.
        """
        if isinstance(det, tuple):
            # The "detector" is a mosaic.
            if det not in self.allowed_mosaics:
                msgs.error(f'{det} is not an allowed mosaic for {self.name}.')
            return self.allowed_mosaics.index(det)+1

        # Single detector
        if det <= 0 or det > self.ndet:
            msgs.error(f'{det} is not a valid detector for {self.name}.')
        return det

    def select_detectors(self, subset=None):
        """
        Vet and return a set of valid detectors or detector mosaics for this
        spectrograph.

        By default, the method returns a list selecting all the detectors
        individually.  A subset can be selected using one of the following
        methods:

            - If ``subset`` is a string, it is assumed that you're selecting a
              set of detector and spatial pixel coordinate combinations needed
              to reduce a single slit; see ``slitspatnum`` in the
              :ref:`pypeitpar`.  The string is parsed into the list of relevant
              detectors.

            - If ``subset`` is a list, integer, or tuple, it is parsed into a
              set of single detector or detector mosaics to reduce.

        Args:
            subset (:obj:`int`, :obj:`tuple`, :obj:`list`, :obj:`str`, optional):
                An object used select a subset of detectors to reduce.  See
                description above.  Note detectors are 1-indexed.

        Returns:
            :obj:`list`: List of unique detectors or detector mosaics to be
            reduced.

        Raises:
            PypeItError: Raised if any of the detectors or detector mosaics
            specified by ``subset`` are invalid.
        """
        if subset is None:
            return np.arange(1, self.ndet+1).tolist()

        if isinstance(subset, str):
            _subset = parse.parse_slitspatnum(subset)[0].tolist()
            # Convert detector to int/tuple
            new_dets = []
            for item in _subset:
                if 'DET' in item:
                    idx = np.where(self.list_detectors() == item)[0][0]
                    new_dets.append(idx+1)
                elif 'MSC' in item:
                    idx = np.where(self.list_detectors(mosaic=True) == item)[0][0]
                    new_dets.append(self.allowed_mosaics[idx])
            _subset = new_dets
        elif isinstance(subset, (int, tuple)):
            _subset = [subset]
        else:
            _subset = subset

        allowed = np.arange(1, self.ndet+1).tolist() + self.allowed_mosaics
        if any([s not in allowed for s in _subset]):
            msgs.error('Selected detectors or detector mosaics contain invalid values.')

        # Require the list contains unique items
        # DP: I had to modify this, because list(set(_subset)) was changing the order of the detectors
        unique_subset = []
        for item in _subset:
            if item not in unique_subset:
                unique_subset.append(item)

        return unique_subset

    @property
    def default_mosaic(self):
        """
        Return the default detector mosaic.

        For instruments with no allowed detector mosaics, this *must* be
        returned as None.
        """
        return None

    def validate_det(self, det):
        """
        Validate the detector(s) provided.

        Args:
            det (:obj:`int`, :obj:`tuple`):
                Selection of 1-indexed detector(s).  An image mosaic is selected
                using a :obj:`tuple` with the detectors in the mosaic, which
                must be one of the allowed mosaics returned by
                :func:`allowed_mosaics`.

        Returns:
            :obj:`tuple`: Returns the number of detectors and a tuple with the
            validated detectors.  If ``det`` is provided as a single integer,
            the latter is just a one-element tuple.  Otherwise, it is identical
            to the input ``det`` tuple, assuming the validation was successful.
        """
        if isinstance(det, tuple):
            if det not in self.allowed_mosaics:
                msgs.error(f'Selected detectors {det} are not an allowed mosaic for {self.name}.')
            return len(det), det
        if not isinstance(det, (int, np.integer)):
            msgs.error(f'Provided det must have type tuple or integer, not {type(det)}.')
        return 1, (det,)

    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces that are key
        for image processing.

        .. warning::

            - When reading multiple detectors for a mosaic, this function
              expects all detector arrays to have exactly the same shape.

        Parameters
        ----------
        raw_file : :obj:`str`, `Path`_
            File to read
        det : :obj:`int`, :obj:`tuple`
            1-indexed detector(s) to read.  An image mosaic is selected using a
            :obj:`tuple` with the detectors in the mosaic, which must be one of
            the allowed mosaics returned by :func:`allowed_mosaics`.

        Returns
        -------
        detector_par : :class:`~pypeit.images.detector_container.DetectorContainer`, :class:`~pypeit.images.mosaic.Mosaic`
            Detector metadata parameters for one or more detectors.
        raw_img : `numpy.ndarray`_
            Raw image for this detector.  Shape is 2D if a single detector image
            is read and 3D if multiple detectors are read.  E.g., the image from
            the first detector in the tuple is accessed using ``raw_img[0]``.
        hdu : `astropy.io.fits.HDUList`_
            Opened fits file
        exptime : :obj:`float`
            Exposure time *in seconds*.
        rawdatasec_img : `numpy.ndarray`_
            Data (Science) section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.  Shape
            is identical to ``raw_img``.
        oscansec_img : `numpy.ndarray`_
            Overscan section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.  Shape
            is identical to ``raw_img``.
        """
        # Check extension and then open
        self._check_extensions(raw_file)
        hdu = io.fits_open(raw_file)

        # Validate the entered (list of) detector(s)
        nimg, _det = self.validate_det(det)

        # Grab the detector or mosaic parameters
        mosaic = None if nimg == 1 else self.get_mosaic_par(det, hdu=hdu)
        detectors = [self.get_detector_par(det, hdu=hdu)] if nimg == 1 else mosaic.detectors

        # Grab metadata from the header
        # NOTE: These metadata must be *identical* for all images when reading a
        # mosaic
        headarr = self.get_headarr(hdu)

        # Exposure time (used by RawImage)
        # NOTE: This *must* be (converted to) seconds.
        exptime = self.get_meta_value(headarr, 'exptime')

        # Rawdatasec, oscansec images
        binning = self.get_meta_value(headarr, 'binning')
        # NOTE: This means that `specaxis` must be the same for all detectors in
        # a mosaic
        if detectors[0]['specaxis'] == 1:
            binning_raw = (',').join(binning.split(',')[::-1])
        else:
            binning_raw = binning

        raw_img = [None]*nimg
        rawdatasec_img = [None]*nimg
        oscansec_img = [None]*nimg
        for i in range(nimg):

            # Raw image
            raw_img[i] = hdu[detectors[i]['dataext']].data.astype(float)
            # Raw data from some spectrograph (i.e. FLAMINGOS2) have an addition
            # extention, so I add the following two lines. It's easier to change
            # here than writing another get_rawimage function in the
            # spectrograph file.
            # TODO: This feels dangerous, but not sure what to do about it...
            if raw_img[i].ndim != 2:
                raw_img[i] = np.squeeze(raw_img[i])
            if raw_img[i].ndim != 2:
                msgs.error(f"Raw images must be 2D; check extension {detectors[i]['dataext']} "
                           f"of {raw_file}.")

            for section in ['datasec', 'oscansec']:

                # Get the data section
                # Try using the image sections as header keywords
                # TODO -- Deal with user windowing of the CCD (e.g. Kast red)
                #  Code like the following maybe useful
                #hdr = hdu[detector[det - 1]['dataext']].header
                #image_sections = [hdr[key] for key in detector[det - 1][section]]
                # Grab from Detector
                image_sections = detectors[i][section]
                #if not isinstance(image_sections, list):
                #    image_sections = [image_sections]
                # Always assume normal FITS header formatting
                one_indexed = True
                include_last = True

                # Initialize the image (0 means no amplifier)
                pix_img = np.zeros(raw_img[i].shape, dtype=int)
                for j in range(detectors[i]['numamplifiers']):

                    if image_sections is not None:  # and image_sections[i] is not None:
                        # Convert the data section from a string to a slice
                        datasec = parse.sec2slice(image_sections[j], one_indexed=one_indexed,
                                                  include_end=include_last, require_dim=2,
                                                  binning=binning_raw)
                        # Assign the amplifier
                        pix_img[datasec] = j+1

                # Finish
                if section == 'datasec':
                    rawdatasec_img[i] = pix_img.copy()
                else:
                    oscansec_img[i] = pix_img.copy()

        if nimg == 1:
            # Return single image
            return detectors[0], raw_img[0], hdu, exptime, rawdatasec_img[0], oscansec_img[0]

        if any([img.shape != raw_img[0].shape for img in raw_img[1:]]):
            msgs.error('All raw images in a mosaic must have the same shape.')
        # Return all images for mosaic
        return mosaic, np.array(raw_img), hdu, exptime, np.array(rawdatasec_img), \
                np.array(oscansec_img)

    def get_lamps_status(self, headarr):
        """
        Return a string containing the information on the lamp status.

        Args:
            headarr (:obj:`list`):
                A list of 1 or more `astropy.io.fits.Header`_ objects.

        Returns:
            :obj:`str`: A string that uniquely represents the lamp status.
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

    def get_meta_value(self, inp, meta_key, required=False,
                       ignore_bad_header=False,
                       usr_row=None, no_fussing=False):
        """
        Return meta data from a given file (or its array of headers).

        Args:
            inp (:obj:`str`, `Path`_, `astropy.io.fits.Header`_, :obj:`list`):
                Input filename, an `astropy.io.fits.Header`_ object, or a list
                of `astropy.io.fits.Header`_ objects.  If None, function simply
                returns None without issuing any warnings/errors, unless
                ``required`` is True.
            meta_key (:obj:`str`, :obj:`list`):
                A (list of) string(s) with the keywords to read from the file
                header(s).
            required (:obj:`bool`, optional):
                The metadata are required and must be available. If it is not,
                the method will raise an exception.  Metadata requirements can
                be globally defined and/or frame-type specific for each
                spectrograph, which will override any value given here.  See the
                ``required`` and ``required_ftype`` keyword in
                ``self.meta[meta_key]``.
            ignore_bad_header (:obj:`bool`, optional):
                PypeIt expects certain metadata values to have specific
                datatypes. If the keyword finds the appropriate data but it
                cannot be cast to the correct datatype, this parameter
                determines whether or not the method raises an exception. If
                True, the incorrect type is ignored. It is recommended that this
                be False unless you know for sure that PypeIt can proceed
                appropriately.  This flag takes precedence over ``required``;
                i.e., if ``ignore_bad_header`` is True, ``required`` is ignored.
            usr_row (`astropy.table.Table`_, optional):
                A single row table with the user-supplied frametype. This is
                used to determine if the metadata value is required for each
                frametype. Must contain a columns called `frametype`;
                everything else is ignored.
            no_fussing (:obj:`bool`, optional):
                No type checking or anything. Just pass back the first value
                retrieved. This is mainly for bound pairs of meta, e.g.
                ra/dec.

        Returns:
            Value recovered for (each) keyword.  Can be None.
        """
        if isinstance(inp, (str, Path, fits.HDUList)):
            headarr = self.get_headarr(inp)
        elif inp is None or isinstance(inp, list):
            headarr = inp
        elif isinstance(inp, fits.Header):
            headarr = [inp]
        else:
            msgs.error(f'Unrecognized type for input: {type(inp)}')
        
        if headarr is None:
            if required:
                msgs.error(f'Unable to access required metadata value for {meta_key}.  Input is '
                           f'either a bad file or an invalid argument to get_meta_value: {inp}.')
            return None

        # Loop?
        if isinstance(meta_key, list):
            return [self.get_meta_value(headarr, key, required=required) for key in meta_key]

        # Are we prepared to provide this meta data?
        if meta_key not in self.meta.keys():
            if required:
                msgs.error("Need to allow for meta_key={} in your meta data".format(meta_key))
            else:
                msgs.warn("Requested meta data for meta_key={} does not exist...".format(meta_key))
                return None

        # Is this meta required for this frame type (Spectrograph specific)
        if 'required_ftypes' in self.meta[meta_key] and usr_row is not None:
            required = False
            for ftype in self.meta[meta_key]['required_ftypes']:
                if ftype in usr_row['frametype']:
                    required = True

        # Check if this meta key is required
        if 'required' in self.meta[meta_key].keys():
            required = self.meta[meta_key]['required']

        # Is this not derivable?  If so, use the default
        #   or search for it as a compound method
        value = None
        try:
            # TODO: Change so that 'card' isn't a required keyword?  ala:
            # if 'card' not in self.meta[meta_key].keys() or self.meta[meta_key]['card'] is None:
            if self.meta[meta_key]['card'] is None:
                if 'default' in self.meta[meta_key].keys():
                    value = self.meta[meta_key]['default']
                elif 'compound' in self.meta[meta_key].keys():
                    value = self.compound_meta(headarr, meta_key)
                else:
                    msgs.error(f"Failed to load spectrograph value for meta: {meta_key}")
            else:
                # Grab from the header, if we can
                value = headarr[self.meta[meta_key]['ext']][self.meta[meta_key]['card']]
        except (KeyError, TypeError) as e:
            if ignore_bad_header or not required:
                msgs.warn(f"Bad Header key ({meta_key}), but we'll try to continue on..")
            else:
                raise e

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
                        msgs.error('Required meta "{0}" did not load!'.format(meta_key)
                                   + 'You may have a corrupt header.')
                else:
                    msgs.warn('Required card {0} missing '.format(self.meta[meta_key]['card'])
                              + 'from your header.  Proceeding with risk...')
            return None

        # Return
        return retvalue

    def get_wcs(self, hdr, slits, platescale, wave0, dwv, spatial_scale=None):
        """
        Construct/Read a World-Coordinate System for a frame.

        This is undefined in the base class.

        Args:
            hdr (`astropy.io.fits.Header`_):
                The header of the raw frame. The information in this
                header will be extracted and returned as a WCS.
            slits (:class:`~pypeit.slittrace.SlitTraceSet`):
                Slit traces.
            platescale (:obj:`float`):
                The platescale of an unbinned pixel in arcsec/pixel (e.g.
                detector.platescale).
            wave0 (:obj:`float`):
                The wavelength zeropoint.
            dwv (:obj:`float`):
                Change in wavelength per spectral pixel.
            spatial_scale (:obj:`float`, None, optional):
                The spatial scale (units=arcsec/pixel) of the WCS to be used.
                This variable is fixed, and is independent of the binning.
                If spatial_scale is set, it will be used for the spatial size
                of the WCS and the platescale will be ignored. If None, then
                the platescale will be used.

        Returns:
            `astropy.wcs.wcs.WCS`_: The world-coordinate system.
        """
        msgs.warn("No WCS setup for spectrograph: {0:s}".format(self.name))
        return None

    def get_datacube_bins(self, slitlength, minmax, num_wave):
        r"""
        Calculate the bin edges to be used when making a datacube.

        Args:
            slitlength (:obj:`int`):
                Length of the slit in pixels
            minmax (`numpy.ndarray`_):
                An array with the minimum and maximum pixel locations on each
                slit relative to the reference location (usually the centre
                of the slit). Shape must be :math:`(N_{\rm slits},2)`, and is
                typically the array returned by
                :func:`~pypeit.slittrace.SlitTraceSet.get_radec_image`.
            num_wave (:obj:`int`):
                Number of wavelength steps.  Given by::
                    int(round((wavemax-wavemin)/delta_wave))

        Returns:
            :obj:`tuple`: Three 1D `numpy.ndarray`_ providing the bins to use
            when constructing a histogram of the spec2d files. The elements
            are :math:`(x,y,\lambda)`.
        """
        msgs.warn("No datacube setup for spectrograph: {0:s}".format(self.name))
        return None

    def fit_2d_det_response(self, det_resp, gpmask):
        r"""
        Perform a 2D model fit to the instrument-specific detector response.

        Args:
            det_resp(`numpy.ndarray`_):
                An image of the detector response.
            gpmask (`numpy.ndarray`_):
                Good pixel mask (True=good), the same shape as ff_struct.

        Returns:
            `numpy.ndarray`_: A model fit to the detector response.
        """
        msgs.warn("2D detector response is not implemented for spectrograph: {0:s}".format(self.name))
        return np.ones_like(det_resp)

    def validate_metadata(self):
        """
        Validates the definitions of the Spectrograph metadata by making a
        series of comparisons to the metadata model defined by
        :func:`pypeit.core.meta.define_core_meta` and :attr:`meta`.
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
            inp (:obj:`str`, `Path`_, `astropy.io.fits.HDUList`_):
                Name of the file to read or the previously opened HDU list.  If
                None, the function will simply return None.
            strict (:obj:`bool`, optional):
                Function will fault if :func:`fits.getheader` fails to read
                any of the headers. Set to False to report a warning and
                continue.

        Returns:
            :obj:`list`: A list of `astropy.io.fits.Header`_ objects with the
            extension headers.  If ``strict`` is False and ``inp`` is a file
            name to be opened, the function will return None if
            :func:`~pypeit.io.fits_open` faults for any reason.
        """
        if inp is None:
            return None

        # Faster to open the whole file and then assign the headers,
        # particularly for gzipped files (e.g., DEIMOS)
        if isinstance(inp, (str, Path)):
            self._check_extensions(inp)
            try:
                hdul = io.fits_open(inp)
            except:
                if strict:
                    msgs.error(f'Cannot open {inp}.')
                else:
                    msgs.warn(f'Cannot open {inp}.  Proceeding, but consider removing this file!')
                    return None
        elif isinstance(inp, (list, fits.HDUList)):
            # TODO: If a list, check that the list elements are HDUs?
            hdul = inp
        else:
            msgs.error(f'Input to get_headarr has incorrect type: {type(inp)}.')
        return [hdu.header for hdu in hdul]

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

        Raises:
            NotImplementedError:
                Raised by the base class to denote that any derived class has
                not been properly defined.
        """
        raise NotImplementedError('Frame typing not defined for {0}.'.format(self.name))

    def idname(self, ftype):
        """
        Return the ``idname`` for the selected frame type for this
        instrument.

        Args:
            ftype (:obj:`str`):
                Frame type, which should be one of the keys in
                :class:`~pypeit.core.framematch.FrameTypeBitMask`.

        Returns:
            :obj:`str`: The value of ``idname`` that should be available in
            the :class:`~pypeit.metadata.PypeItMetaData` instance that
            identifies frames of this type.

        Raises:
            NotImplementedError:
                Raised by the base class to denote that any derived class has
                not been properly defined.
        """
        raise NotImplementedError('Header keyword with frame type not defined for {0}.'.format(
                                  self.name))

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

    def get_echelle_angle_files(self):
        """ Pass back the files required
        to run the echelle method of wavecalib

        Returns:
            list: List of files
        """
        msgs.error(f'Echelle angle files not ready for {self.name}')

    def order_platescale(self, order_vec, binning=None):
        """
        Return the platescale for each echelle order.

        This routine is only defined for echelle spectrographs, and it is
        undefined in the base class.

        Args:
            order_vec (`numpy.ndarray`_):
                The vector providing the order numbers.
            binning (:obj:`str`, optional):
                The string defining the spectral and spatial binning.

        Returns:
            `numpy.ndarray`_: An array with the platescale for each order
            provided by ``order``.
        """
        pass

    @property
    def norders(self):
        """
        Number of orders for this spectograph. Should only defined for
        echelle spectrographs, and it is undefined for the base class.
        """
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
        """ Return the expected spatial position of each echelle order.

        This is for fixed-format echelle spectrographs (e.g. X-Shooter)
        And is measured 1/2 way up the chip (spectral)
        and in normalized units (0-1)

        Returns:
            `numpy.ndarray`_: An array with values. 
            The length of the provided array much 
            match self.norders
        """
        return None

    @property
    def orders(self):
        """
        Return the order number for each echelle order.

        Returns:
            `numpy.ndarray`_: An array with values. 
            Order number.  Must have lenght of self.norders
        """
        return None

    @property
    def spec_min_max(self):
        """
        Return the minimum and maximum spectral pixel expected for the
        spectral range of each order.
        """
        return None

    @property
    def dloglam(self):
        """
        Return the logarithmic step in wavelength for output spectra.
        """
        return None

    @property
    def loglam_minmax(self):
        """
        Return the base-10 logarithm of the first and last wavelength for
        ouput spectra.
        """
        return None

    # TODO : This code needs serious work.  e.g. eliminate the try/except
    def slit_minmax(self, slit_spat_pos, binspectral=1):
        """
        Adjust the minimum and maximum spectral pixel expected for the
        spectral range of each echelle order by accounting for the spectral
        binning.

        Args:
            slit_spat_pos (:obj:`float`, `numpy.ndarray`_):
                Spatial position of each slit/order normalized by the full
                spatial extent of the detector.
            binspectral (:obj:`int`, optional):
                Number of pixels binned in the spectral direction.

        Returns:
            `numpy.ndarray`_: The minimum and maximum (binned) pixel that
            includes the valid spectra of each slit/order.
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
    
    def spec1d_match_spectra(self, sobjs):
        """
        Match up slits in a :class:`~pypeit.specobjs.SpecObjs` object.

        This typically done across multiple detectors; see
        :func:`pypeit.specrographs.keck_deimos.spec1d_match_spectra`.

        Args:
            sobjs (:class:`~pypeit.specobjs.SpecObjs`):
                Spec1D objects

        Returns:
            :obj:`tuple`: Arrays that provide the indices of slits matched
            across multiple detectors.
        """
        msgs.error(f'Method to match slits across detectors not defined for {self.name}')


    def tweak_standard(self, wave_in, counts_in, counts_ivar_in, gpm_in, meta_table):
        """

        This routine is for performing instrument/disperser specific tweaks to standard stars so that sensitivity
        function fits will be well behaved. For example, masking second order light. For instruments that don't
        require such tweaks it will just return the inputs, but for isntruments that do this function is overloaded
        with a method that performs the tweaks.

        Parameters
        ----------
        wave_in: `numpy.ndarray`_
            Input standard star wavelengths (:obj:`float`, ``shape = (nspec,)``)
        counts_in: `numpy.ndarray`_
            Input standard star counts (:obj:`float`, ``shape = (nspec,)``)
        counts_ivar_in: `numpy.ndarray`_
            Input inverse variance of standard star counts (:obj:`float`, ``shape = (nspec,)``)
        gpm_in: `numpy.ndarray`_
            Input good pixel mask for standard (:obj:`bool`, ``shape = (nspec,)``)
        meta_table: :obj:`dict`
            Table containing meta data that is slupred from the :class:`~pypeit.specobjs.SpecObjs`
            object.  See :meth:`~pypeit.specobjs.SpecObjs.unpack_object` for the
            contents of this table.

        Returns
        -------
        wave_out: `numpy.ndarray`_
            Output standard star wavelengths (:obj:`float`, ``shape = (nspec,)``)
        counts_out: `numpy.ndarray`_
            Output standard star counts (:obj:`float`, ``shape = (nspec,)``)
        counts_ivar_out: `numpy.ndarray`_
            Output inverse variance of standard star counts (:obj:`float`, ``shape = (nspec,)``)
        gpm_out: `numpy.ndarray`_
            Output good pixel mask for standard (:obj:`bool`, ``shape = (nspec,)``)
        """
        return wave_in, counts_in, counts_ivar_in, gpm_in

    def calc_pattern_freq(self, frame, rawdatasec_img, oscansec_img, hdu):
        """
        Calculate the pattern frequency using the overscan region that covers
        the overscan and data sections. Using a larger range allows the
        frequency to be pinned down with high accuracy.

        Parameters
        ----------
        frame : `numpy.ndarray`_
            Raw data frame to be used to estimate the pattern frequency.
        rawdatasec_img : `numpy.ndarray`_
            Array the same shape as ``frame``, used as a mask to identify the
            data pixels (0 is no data, non-zero values indicate the amplifier
            number).
        oscansec_img : `numpy.ndarray`_
            Array the same shape as ``frame``, used as a mask to identify the
            overscan pixels (0 is no data, non-zero values indicate the
            amplifier number).
        hdu : `astropy.io.fits.HDUList`_
            Opened fits file.

        Returns
        -------
        patt_freqs : :obj:`list`
            List of pattern frequencies.
        """
        msgs.info("Pattern noise removal is not implemented for spectrograph {0:s}".format(self.name))
        return []


    def __repr__(self):
        """Return a string representation of the instance."""
        txt = '<{:s}: '.format(self.__class__.__name__)
        txt += ' spectrograph={:s},'.format(self.name)
        txt += ' telescope={:s},'.format(self.telescope['name'])
        txt += ' pypeline={:s},'.format(self.pypeline)
        txt += '>'
        return txt


