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
from pypeit.pypmsgs import PypeItError
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

    # TODO: Add an astropy.Table into the base-class data model that includes
    # the subset of `fitstbl` with the metadata for the raw calibration frames?
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

        Nothing is returned; this function is used to set :attr:`calib_dir`,
        :attr:`calib_id`, and :attr:`calib_key`.

        Args:
            odir (:obj:`str`, `Path`_):
                Output directory for the processed calibration frames
            setup (:obj:`str`):
                The string identifier for the instrument setup/configuration;
                see :func:`~pypeit.metadata.PypeItMetaData.unique_configurations`.
            calib_id (:obj:`str`, :obj:`list`, :obj:`int`):
                Identifiers for one or more calibration groups for this
                calibration frame.  Strings (either as individually entered or
                as elements of a provided list) can be single or comma-separated
                integers.  Otherwise, all strings must be convertible to
                integers; the only exception is the string 'all'.
            detname (:obj:`str`):
                The identifier used for the detector or detector mosaic for the
                relevant instrument; see
                :func:`~pypeit.spectrograph.spectrograph.Spectrograph.get_det_name`.
        """
        self.calib_dir = Path(odir).resolve()
        # TODO: Keep this, or throw an error if the directory doesn't exist instead?
        if not self.calib_dir.exists():
            self.calib_dir.mkdir(parents=True)
        # TODO: Use Path object instead of string here?
        self.calib_dir = str(self.calib_dir)
        self.calib_id = CalibFrame.ingest_calib_id(calib_id)
        self.calib_key = self.construct_calib_key(setup, self.calib_id, detname)

    def copy_calib_internals(self, other):
        """
        Copy the internals from another :class:`CalibFrame` to this one.

        Args:
            other (:class:`CalibFrame`):
                Object to copy from.
        """
        for attr in CalibFrame.internals:
            setattr(self, attr, getattr(other, attr))

    # NOTE: Only need to overload to_file because the only thing special about
    # CalibFrame is that the paths are pre-defined.
    def to_file(self, file_path=None, overwrite=True, **kwargs):
        """
        Overrides the base-class function, forcing the naming convention.

        Args:
            file_path (:obj:`str`, `Path`_, optional):
                Full path for the file to be written.  This should be used
                *very* rarely.  The whole point of the :class:`CalibFrame` is to
                follow a deterministic I/O naming structure, and use of this
                option circumvents that.  You should instead typically use
                :func:`set_paths` so that the file path is defined
                automatically.
            overwrite (:obj:`bool`, optional):
                Flag to overwrite any existing files.  This overrides the
                default of the base class, meaning that anytime a calibration
                frame is written to disk it will overwrite any existing files by
                default!
            **kwargs (optional):
                Passed directly to
                :func:`~pypeit.datamodel.DataContainer.to_file`.
        """
        _file_path = self.get_path() if file_path is None else Path(file_path).resolve()
        super().to_file(_file_path, overwrite=overwrite, **kwargs)

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
        self.calib_keys_from_header(hdu[parsed_hdus[0]].header)
        return self

    def calib_keys_from_header(self, hdr):
        """
        (Attempt to) Fill the calibration keys based on the provided hdr.

        If successful, this sets the values for the calibration
        :attr:`internals`.

        Args:
            hdr (`astropy.io.fits.Header`_):
                Header to parse
        """
        try:
            self.calib_key, self.calib_dir = CalibFrame.parse_key_dir(hdr)
        except PypeItError as e:
            msgs.warn(f'{e}')
        if 'CALIBID' in hdr:
            self.calib_id = self.ingest_calib_id(hdr['CALIBID'])
        else:
            msgs.warn('Header does not have CALIBID card; cannot parse calibration IDs.')
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
                    msgs.error(f'None of the headers in {inp} have both CALIBKEY and CALIBDIR '
                               'keywords!')
                return hdu[ext].header['CALIBKEY'], hdu[ext].header['CALIBDIR']

        if isinstance(inp, fits.Header):
            if 'CALIBKEY' not in inp or 'CALIBDIR' not in inp:
                msgs.error('Header does not include CALIBKEY and/or CALIBDIR.')
            return inp['CALIBKEY'], inp['CALIBDIR']

        msgs.error(f'Input object must have type str or astropy.io.fits.Header, not {type(inp)}.')

    @staticmethod
    def ingest_calib_id(calib_id):
        """
        Ingest the calibration group IDs, converting the input into a list of strings.

        Args:
            calib_id (:obj:`str`, :obj:`list`, :obj:`int`):
                Identifiers for one or more calibration groups for this
                calibration frame.  Strings (either as individually entered or
                as elements of a provided list) can be single or comma-separated
                integers.  Otherwise, all strings must be convertible to
                integers; the only exception is the string 'all'.

        Returns:
            :obj:`list`: List of string representations of single calibration
            group integer identifiers.

        Examples:

            >>> CalibFrame.ingest_calib_id('all')
            ['all']
            >>> CalibFrame.ingest_calib_id(['all', 1])
            [WARNING] :: Calibration groups set to ['1' 'all'], resetting to simply "all".
            ['all']
            >>> CalibFrame.ingest_calib_id('1,2')
            ['1', '2']
            >>> CalibFrame.ingest_calib_id(['1,2', '5,8', '3'])
            ['1', '2', '3', '5', '8']
            >>> CalibFrame.ingest_calib_id([2, 1, 2])
            ['1', '2']

        """
        if isinstance(calib_id, str):
            _calib_id = calib_id.split(',')
        elif isinstance(calib_id, list):
            _calib_id = calib_id
        else:
            _calib_id = [calib_id]
        _calib_id = np.unique(np.concatenate([str(c).split(',') for c in _calib_id]))
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
        group(s), and the detector.  The configuration ID is the same as
        included in the configuration column (A, B, C, etc), the calibration
        group is a dash-separated list of the calibration group identifiers or
        "all", and the detector/mosaic identifier (e.g., DET01, MSC02) is set by
        the detector number or mosaic tuple (see
        :func:`~pypeit.spectrographs.spectrograph.Spectrograph.get_det_name`).

        Args:
            setup (:obj:`str`):
                The string identifier for the instrument setup/configuration;
                see :func:`~pypeit.metadata.PypeItMetaData.unique_configurations`.
            calib_id (:obj:`str`, :obj:`list`, :obj:`int`):
                Identifiers for one or more calibration groups for this
                calibration frame.  See :func:`ingest_calib_id`.
            detname (:obj:`str`):
                The identifier used for the detector or detector mosaic for the
                relevant instrument; see
                :func:`~pypeit.spectrograph.spectrograph.Spectrograph.get_det_name`.

        Returns:
            :obj:`str`: Calibration identifier.
        """
        return f'{setup}_{"-".join(CalibFrame.ingest_calib_id(calib_id))}_{detname}'

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
            :obj:`str`, `Path`_: File path if ``calib_dir`` is provided,
            otherwise the file name
        """
        if None in [cls.calib_type, cls.calib_file_format]:
            msgs.error(f'CODING ERROR: {cls.__name__} does not have all '
                       'the attributes needed to construct its filename.')
        if calib_key is None:
            msgs.error('CODING ERROR: calib_key cannot be None when constructing the '
                       f'{cls.__name__} file name.')
        filename = f'{cls.calib_type}_{calib_key}.{cls.calib_file_format}'
        return filename if calib_dir is None else Path(calib_dir).resolve() / filename

    def get_path(self):
        """
        Return the path to the output file.

        This is a simple wrapper for the :func:`construct_file_name` classmethod
        that uses the existing values of :attr:`calib_key` and :attr`calib_dir`.

        Returns:
            :obj:`str`, `Path`_: File path or file name.  This is always the
            full path if :attr:`calib_dir` is defined.
        """
        return self.__class__.construct_file_name(self.calib_key, calib_dir=self.calib_dir)

    def _base_header(self, hdr=None):
        """
        Override the base class method to add useful/identifying internals to
        the header.

        Args:
            hdr (`astropy.io.fits.Header`, optional):
                Header object to update.  The provided object is *not* edited,
                only copied.

        Returns:
            `astropy.io.fits.Header`_: The new/edited fits header.
        """
        _hdr = super()._base_header(hdr=hdr)
        _hdr['CALIBTYP'] = (self.calib_type, 'PypeIt: Calibration frame type')
        if self.calib_dir is not None:
            _hdr['CALIBDIR'] = (self.calib_dir, 'PypeIt: Calibration file directory')
        if self.calib_key is not None:
            _hdr['CALIBKEY'] = (self.calib_key, 'PypeIt: Calibration key')
        if self.calib_id is not None:
            _hdr['CALIBID'] = (','.join(self.calib_id), 'PypeIt: Calibration groups')
        return _hdr

    @classmethod
    def glob(cls, calib_dir, setup, calib_id, detname=None):
        """
        Search for calibration files.

        Args:
            calib_dir (:obj:`str`, `Path`_):
                Directory to search
            setup (:obj:`str`):
                The setup/configuration identifier (e.g., A, B, C, ...) of the calibrations
            calib_id (:obj:`str`, :obj:`int`):
                The *single* calibration group of the calibrations
            detname (:obj:`str`, optional):
                The identifier of the detector/mosaic of the calibrations.  If
                None, any relevant calibrations are returned.

        Returns:
            :obj:`list`: List of paths to applicable calibration files.  If no
            relevant files are found or if ``calib_dir`` is not an existing
            directory, None is returned.
        """
        # Check the path exists
        _calib_dir = Path(calib_dir).resolve()
        if not _calib_dir.exists():
            return None

        # Construct the search string
        search = f'{cls.calib_type}*{setup}*'
        if detname is not None:
            search += f'{detname}*'
        search += f'{cls.calib_file_format}'

        # Find all the relevant calibrations in the directory
        files = np.array(sorted(_calib_dir.glob(search)))
        if files.size == 0:
            return None

        # For the remaining files, find the ones that have applicable
        # calibration groups
        keep = np.ones(files.size, dtype=bool)
        for i, f in enumerate(files):
            _calib_id = cls.parse_calib_key(cls.parse_key_dir(str(f), from_filename=True)[0])[1]
            if _calib_id == 'all' or str(calib_id) in cls.ingest_calib_id(_calib_id):
                continue
            keep[i] = False

        # Return the applicable calibrations
        return files[keep].tolist() if any(keep) else None

