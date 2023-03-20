"""
Implements the calibration frame base class.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
from pathlib import Path

from IPython import embed

import numpy as np

from astropy.io import fits

from pypeit import msgs
from pypeit import datamodel
from pypeit import io

class CalibFrame(datamodel.DataContainer):
    """
    An abstract class for calibration frames.  The primary purpose of the class
    is to set the naming scheme for all processed calibration files.
    """

    calib_type = None
    """
    The type of the calibration frame, primarily used to set the name of the
    output file.
    """

    calib_file_format = 'fits'
    """
    The extension and file format of the output file.  Should be ``'fits'`` or
    ``'fits.gz'`` (for gzipped output).
    """

    # TODO: Add an astropy.Table with the metadata?    
    datamodel = {'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name')}
    """
    Default datamodel for any :class:`CalibFrame`.  Derived classes should
    instantiate their datamodels by first inheriting from the base class.  E.g.:

    .. code-block:: python

        class ArcFrame(CalibFrame):
            datamodel = {**CalibFrame.datamodel, ...}

    """

    internals = ['calib_id', 'calib_key', 'calib_dir']
    """
    Base class internals.  The :attr:`internals` of any derived class should
    also include these.  E.g.:

    .. code-block:: python

        class ArcFrame(CalibFrame):
            internals = CalibFrame.internals + ['arc_specific_internal']

    """

    def _validate(self):
        """
        Validation method that is executed every time a :class:`CalibFrame` is
        instantiated.

        Ensures:

            - :attr:`calib_type` and :attr:`datamodel` are defined, and

            - any members of :attr:`datamodel` of the base class are also
              members of the derived class.

        """
        if self.calib_type is None:
            msgs.error(f'CODING ERROR: Must define calib_type for {self.__class__.__name__}.')
        if self.datamodel is None:
            msgs.error(f'CODING ERROR: datamodel cannot be None for {self.__class__.__name__}.')
        for key in CalibFrame.datamodel.keys():
            if key not in self.keys():
                msgs.error(f'CODING ERROR: datamodel for {self.__class__.__name__} must inherit '
                           'all datamodel components from CalibFrame.datamodel.')

    def set_paths(self, odir, setup, calib_id, detname):
        """
        Set the internals necessary to construct the IO path for the calibration
        file.

        Nothing is returned; this function is used to set :attr:`calib_dir` and
        :attr:`calib_key`.

        Args:
            odir (:obj:`str`, `Path`_):
                Output directory for the processed calibration frames
            setup (:obj:`str`):
                The string identifier for the instrument setup/configuration;
                see :func:`~pypeit.metadata.PypeItMetaData.unique_configurations`.
            calib_id (:obj:`list`, :obj:`str`, :obj:`int`):
                Identifiers for one or more calibration groups for this
                calibration frame.  The individual string and the string
                elements of the list must be single numbers.  The only exception
                is the string 'all'.
            detname (:obj:`str`):
                The identifier used for the detector or detector mosaic for the
                relevant instrument; see
                :func:`~pypeit.spectrograph.spectrograph.Spectrograph.get_det_name`.
        """
        self.calib_dir = Path(odir).resolve()
        # TODO: Keep this, or throw an error if the directory doesn't exist instead.
        if not self.calib_dir.exists():
            self.calib_dir.mkdir(parents=True)
        # TODO: Would like to transition to using Path objects instead of strings!
        self.calib_dir = str(self.calib_dir)
        self.calib_id = CalibFrame.ingest_calib_id(calib_id)
        self.calib_key = self.construct_calib_key(setup, self.calib_id, detname)

    # NOTE: Only need to overload to_file because the only thing special about
    # CalibFrame is that the paths are pre-defined.
    def to_file(self, **kwargs):
        """
        Overrides the base-class function, forcing the naming convention.

        Args:
            **kwargs (optional):
                Passed directly to
                :func:`~pypeit.datamodel.DataContainer.to_file`.  One of the
                keyword arguments *cannot* be ``'filename'``.  The filenames for
                :class:`CalibFrame` objects are predefined.  If ``'hdr'`` is a
                provided keyword, the CalibFrame-specific keywords are added
                using :func:`build_header` and passed to on to 
                :func:`~pypeit.datamodel.DataContainer.to_file`.
        """
        if 'hdr' in kwargs:
            hdr = self.build_header(hdr=kwargs['hdr'])
            kwargs.pop('hdr')
        else:
            hdr = self.build_header()
        # TODO: This used to set overwrite=True, but overwrite should be part of
        # kwargs.  The default of DataContainer is overwrite=False, so I need to
        # make sure that the behavior is maintained.
        super().to_file(self.get_path(), hdr=hdr, **kwargs)

    @classmethod
    def from_hdu(cls, hdu, chk_version=True, **kwargs):
        """
        Instantiate the object from an HDU extension.

        Args:
            hdu (`astropy.io.fits.HDUList`_, `astropy.io.fits.ImageHDU`_, `astropy.io.fits.BinTableHDU`_):
                The HDU(s) with the data to use for instantiation.
            chk_version (:obj:`bool`, optional):
                If True, raise an error if the datamodel version or
                type check failed. If False, throw a warning only.
            **kwargs:
                Passed directly to :func:`~pypeit.datamodel.DataContainer._parse`.
        """
        d, dm_version_passed, dm_type_passed, parsed_hdus = cls._parse(hdu, **kwargs)
        # Check version and type?
        if not dm_type_passed:
            msgs.error(f'The HDU(s) cannot be parsed by a {cls.__name__} object!')
        if not dm_version_passed:
            _f = msgs.error if chk_version else msgs.warn
            _f(f'Current version of {cls.__name__} object in code ({cls.version}) '
               'does not match version used to write your HDU(s)!')

        # Instantiate
        self = cls.from_dict(d=d)

        # Calibration frame attributes
        # NOTE: If multiple HDUs are parsed, this assumes that the information
        # necessary to set all the calib internals is always in *every* header.
        # BEWARE!
        hdr_to_parse = hdu[parsed_hdus[0]].header
        # TODO: Consider making the next two lines a function so that they don't
        # need to be repeated in PypeItCalibrationImage.
        self.calib_key, self.calib_dir = CalibFrame.parse_key_dir(hdr_to_parse)
        self.calib_id = self.ingest_calib_id(hdr_to_parse['CALIBID'].split(','))
        # Check that the parsed HDUs are for the correct calibration frame type.
        # TODO: This is superfluous because the _parse statement above will
        # already check that the datamodel type passed.  I.e., we allow
        # circumventing the datamodel version, but we never allow circumventing
        # the datamodel *type*.
#        if 'CALIBTYP' in hdr_to_parse:
#            if hdr_to_parse['CALIBTYP'] != cls.calib_type:
#                msgs.error(f'Primary header of {_ifile} does not have the correct calibration '
#                            f'frame type!  Found {hdr_to_parse["CALIBTYP"]}, but expected '
#                            f'{cls.calib_type}.')

        return self

    @staticmethod
    def parse_key_dir(inp, from_filename=False):
        """
        Grab the identifying key and directory by parsing the input.

        Args:
            inp (:obj:`str`, `astropy.io.fits.Header`_):
                Either a filename or a Header of a FITS file
            from_filename (:obj:`bool`, optional):
                If True, ``inp`` must be a string providing the calibration file
                name, which must follow the expected naming convention.  If
                False, ``inp`` must be an `astropy.io.fits.Header`_ or a file
                from which a header can be read.

        Returns:
            :obj:`tuple`:  Two strings with the identifying key and directory of
            the processed calibration frame.
        """
        if from_filename:
            path = Path(inp).resolve()
            return '_'.join(path.name.split('.')[0].split('_')[1:]), str(path.parent)

        if isinstance(inp, str):
            with io.fits_open(inp) as hdu:
                ext = None
                for h in hdu:
                    if 'CALIBKEY' in h.header and 'CALIBDIR' in h.header:
                        ext = h.name
                        break
                if ext is None:
                    msgs.error(f'None of the headers in {inp} have CALIBKEY or CALIBDIR keywords!')
                return hdu[ext].header['CALIBKEY'], hdu[ext].header['CALIBDIR']

        if isinstance(inp, fits.Header):
            if 'CALIBKEY' not in inp or 'CALIBDIR' not in inp:
                msgs.error('Header does not include CALIBKEY or CALIBDIR.')
            return inp['CALIBKEY'], inp['CALIBDIR']

        msgs.error('Input object must have type str or astropy.io.fits.Header, not {type(inp)}.')

    @staticmethod
    def ingest_calib_id(calib_id):
        """
        Ingest the calibration group IDs, converting the input into a list of strings.

        Args:
            calib_id (:obj:`list`, :obj:`str`, :obj:`int`):
                Identifiers for one or more calibration groups for this
                calibration frame.  The individual string and the string
                elements of the list must be single numbers.  The only exception
                is the string 'all'.

        Returns:
            :obj:`list`: List of string representations of single calibration
            groups.
        """
        _calib_id = calib_id if isinstance(calib_id, list) else [calib_id]
        _calib_id = np.unique([str(c) for c in _calib_id])
        if 'all' in _calib_id and len(_calib_id) != 1:
            msgs.warn(f'Calibration groups set to {_calib_id}, resetting to simply "all".')
            _calib_id = np.array(['all'])
        for c in _calib_id:
            if c == 'all':
                continue
            try:
                _c = int(c)
            except ValueError:
                # TODO: Not sure this is strictly necessary
                msgs.error(f'Invalid calibration group {c}; must be convertible to an integer.')
        return _calib_id.tolist()

    @staticmethod
    def construct_calib_key(setup, calib_id, detname):
        """
        Construct the identifier used for a given set of calibrations.

        The identifier is the combination of the configuration, the calibration
        group, and the detector.  The configuration ID is the same as included
        in the configuration column (A, B, C, etc), the calibration group is the
        same as the calibration bit number, and the detector identifier.

        Args:
            setup (:obj:`str`):
                The string identifier for the instrument setup/configuration;
                see :func:`~pypeit.metadata.PypeItMetaData.unique_configurations`.
            calib_id (:obj:`list`, :obj:`str`, :obj:`int`):
                Identifiers for one or more calibration groups for this
                calibration frame.  The individual string and the string
                elements of the list must be single numbers.  The only exception
                is the string 'all'.
            detname (:obj:`str`):
                The identifier used for the detector or detector mosaic for the
                relevant instrument; see
                :func:`~pypeit.spectrograph.spectrograph.Spectrograph.get_det_name`.

        Returns:
            :obj:`str`: Calibration identifier.
        """
        _calib_id = CalibFrame.ingest_calib_id(calib_id)
        return f'{setup}_{"-".join(_calib_id)}_{detname}'

    @staticmethod
    def parse_calib_key(calib_key):
        """
        Given the calibration key identifier, parse its different components.

        To see how the key is constructed, see :func:`construct_calib_key`.

        Args:
            calib_key (:obj:`str`):
                The calibration key identifier to parse.

        Returns:
            :obj:`tuple`: The three components of the calibration key.
        """
        setup, calib_id, detname = calib_key.split('_')
        return setup, ','.join(calib_id.split('-')), detname

    @classmethod
    def construct_file_name(cls, calib_key, calib_dir=None):
        """
        Generate a calibration frame filename.

        Args:
            calib_key (:obj:`str`):
                String identifier of the calibration group.  See
                :func:`construct_calib_key`.
            calib_dir (:obj:`str`, `Path`_, optional):
                If provided, return the full path to the file given this
                directory.

        Returns:
            :obj:`str`: File path or file name
        """
        if None in [cls.calib_type, cls.calib_file_format]:
            msgs.error(f'CODING ERROR: {cls.__name__} does not have all '
                       'the attributes needed to construct its filename.')
        if calib_key is None:
            msgs.error('CODING ERROR: calib_key cannot be None when constructing the '
                       f'{cls.__name__} file name.')
        filename = f'{cls.calib_type}_{calib_key}.{cls.calib_file_format}'
        if calib_dir is None:
            return filename
        return str(Path(calib_dir).resolve() / filename)

    def get_path(self):
        """
        Return the path to the output file.

        This is a simple wrapper for the :func:`construct_file_name` classmethod
        that uses the existing values of :attr:`calib_key` and :attr`calib_dir`.

        Returns:
            :obj:`str`: File path or file name.  This is always the full path if
            :attr:`calib_dir` is defined.
        """
        return self.__class__.construct_file_name(self.calib_key, calib_dir=self.calib_dir)

    def build_header(self, hdr=None):
        """
        Initialize the calibration frame header.

        This builds a generic header that is written to all PypeIt calibration
        frames.

        Args:
            hdr (`astropy.io.fits.Header`, optional):
                Header object to update.  The object is modified in-place and
                also returned. If None, an empty header is instantiated, edited,
                and returned.

        Returns:
            `astropy.io.fits.Header`_: The initialized (or edited) fits header.
        """
        # Standard init
        _hdr = io.initialize_header(hdr)

        # Save the calibration frame type and key and version, in case the file
        # name is changed.
        _hdr['CALIBTYP'] = (self.calib_type, 'PypeIt: Calibration frame type')
        _hdr['CALIBVER'] = (self.version, 'PypeIt: Calibration datamodel version')
        if self.calib_dir is not None:
            _hdr['CALIBDIR'] = (self.calib_dir, 'PypeIt: Calibration file directory')
        if self.calib_key is not None:
            _hdr['CALIBKEY'] = (self.calib_key, 'PypeIt: Calibration key')
        if self.calib_id is not None:
            _hdr['CALIBID'] = (','.join(self.calib_id), 'PypeIt: Calibration key')

        # Return
        return _hdr


