"""
Provides a class that handles the fits metadata required by PypeIt.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os
from pathlib import Path
import io
import string
from copy import deepcopy
import datetime

from IPython import embed

import numpy as np

from astropy import table, time

from pypeit import msgs
from pypeit import inputfiles
from pypeit.core import framematch
from pypeit.core import flux_calib
from pypeit.core import parse
from pypeit.core import meta
from pypeit.io import dict_to_lines
from pypeit.par import PypeItPar
from pypeit.bitmask import BitMask


# TODO: Turn this into a DataContainer
# Initially tried to subclass this from astropy.table.Table, but that
# proved too difficult.
class PypeItMetaData:
    """
    Provides a table and interface to the relevant fits file metadata
    used during the reduction.

    The content of the fits table is dictated by the header keywords
    specified for the provided spectrograph. It is expected that this
    table can be used to set the frame type of each file.

    The metadata is validated using checks specified by the provided
    spectrograph class.

    For the data table, one should typically provide either the file
    list from which to grab the data from the fits headers or the
    data directly. If neither are provided the table is instantiated
    without any data.

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph used to collect the data save to each file.
            The class is used to provide the header keyword data to
            include in the table and specify any validation checks.
        par (:obj:`pypeit.par.pypeitpar.PypeItPar`):
            PypeIt parameters used to set the code behavior.
        files (:obj:`str`, :obj:`list`, optional):
            The list of files to include in the table.
        data (table-like, optional):
            The data to include in the table.  The type can be anything
            allowed by the instantiation of
            :class:`astropy.table.Table`.
        usrdata (:obj:`astropy.table.Table`, optional):
            A user provided set of data used to supplement or overwrite
            metadata read from the file headers.  The table must have a
            `filename` column that is used to match to the metadata
            table generated within PypeIt.  **Note**: This is ignored if
            `data` is also provided.  This functionality is only used
            when building the metadata from the fits files.
        strict (:obj:`bool`, optional):
            Function will fault if there is a problem with the reading
            the header for any of the provided files; see
            :func:`~pypeit.spectrographs.spectrograph.get_headarr`.  Set
            to False to instead report a warning and continue.

    Attributes:
        spectrograph
            (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph used to collect the data save to each file.
            The class is used to provide the header keyword data to
            include in the table and specify any validation checks.
        par (:class:`pypeit.par.pypeitpar.PypeItPar`):
            PypeIt parameters used to set the code behavior.  If not
            provided, the default parameters specific to the provided
            spectrograph are used.
        configs (:obj:`dict`):
            A dictionary of the unique configurations identified.
        type_bitmask (:class:`pypeit.core.framematch.FrameTypeBitMask`):
            The bitmask used to set the frame type of each fits file.
        calib_bitmask (:class:`BitMask`):
            The bitmask used to keep track of the calibration group bits.
        table (:class:`astropy.table.Table`):
            The table with the relevant metadata for each fits file to
            use in the data reduction.
    """
    def __init__(self, spectrograph, par, files=None, data=None, usrdata=None, 
                 strict=True):

        if data is None and files is None:
            # Warn that table will be empty
            msgs.warn('Both data and files are None in the instantiation of PypeItMetaData.'
                      '  The table will be empty!')

        # Initialize internals
        self.spectrograph = spectrograph
        self.par = par
        if not isinstance(self.par, PypeItPar):
            raise TypeError('Input parameter set must be of type PypeItPar.')
        self.type_bitmask = framematch.FrameTypeBitMask()

        # Build table
        self.table = table.Table(data if files is None 
                                 else self._build(files, strict=strict, 
                                                  usrdata=usrdata))

        # Merge with user data, if present
        if usrdata is not None:
            self.merge(usrdata)

        # Impose types on specific columns
        self._impose_types(['comb_id', 'bkg_id', 'manual'], [int, int, str])

        # Initialize internal attributes
        self.configs = None
        self.calib_bitmask = None

        # Initialize columns that the user might add
        self.set_user_added_columns()
        # Validate instrument name
        self.spectrograph.vet_instrument(self.table)

    def _impose_types(self, columns, types):
        """
        Impose a set of types on certain columns.

        .. note::
            :attr:`table` is edited in place.

        Args:
            columns (:obj:`list`):
                List of column names
            types (:obj:`list`):
                List of types
        """
        for c,t in zip(columns, types):
            if c in self.keys():
                self.table[c] = self.table[c].astype(t)

    def _build(self, files, strict=True, usrdata=None):
        """
        Generate the fitstbl that will be at the heart of PypeItMetaData.

        Args:
            files (:obj:`str`, :obj:`list`):
                One or more files to use to build the table.
            strict (:obj:`bool`, optional):
                Function will fault if :func:`fits.getheader` fails to
                read any of the headers.  Set to False to report a
                warning and continue.
            usrdata (astropy.table.Table, optional):
                Parsed for frametype for a few instruments (e.g. VLT)
                where meta data may not be required

        Returns:
            dict: Dictionary with the data to assign to :attr:`table`.

        """
        # Allow for single files
        _files = files if hasattr(files, '__len__') else [files]

        # Build lists to fill
        data = {k:[] for k in self.spectrograph.meta.keys()}
        data['directory'] = ['None']*len(_files)
        data['filename'] = ['None']*len(_files)

        # Build the table
        for idx, ifile in enumerate(_files):
            _ifile = Path(ifile).resolve()
            # User data (for frame type)
            if usrdata is None:
                usr_row = None
            else:
                # TODO: This check should be done elsewhere
                # Check
                if _ifile.name != usrdata['filename'][idx]:
                    msgs.error('File name list does not match user-provided metadata table.  See '
                               'usrdata argument of instantiation of PypeItMetaData.')
                usr_row = usrdata[idx]

            # Add the directory and file name to the table
            data['directory'][idx] = str(_ifile.parent)
            data['filename'][idx] = _ifile.name
            if not data['directory'][idx]:
                data['directory'][idx] = '.'

            # Read the fits headers.  NOTE: If the file cannot be opened,
            # headarr will be None, and the subsequent loop over the meta keys
            # will fill the data dictionary with None values.
            msgs.info(f'Adding metadata for {data["filename"][idx]}')
            headarr = self.spectrograph.get_headarr(_ifile, strict=strict)

            # Grab Meta
            for meta_key in self.spectrograph.meta.keys():
                value = self.spectrograph.get_meta_value(headarr, meta_key, 
                                                         required=strict,
                                                         usr_row=usr_row, 
                        ignore_bad_header = (
                            self.par['rdx']['ignore_bad_headers'] or strict))
                if isinstance(value, str) and '#' in value:
                    value = value.replace('#', '')
                    msgs.warn('Removing troublesome # character from {0}.  Returning {1}.'.format(
                              meta_key, value))
                data[meta_key].append(value)

        # JFH Changed the below to not crash if some files have None in
        # their MJD. This is the desired behavior since if there are
        # empty or corrupt files we still want this to run.

        # Validate, print out a warning if there is problem
        try:
            time.Time(data['mjd'], format='mjd')
        except ValueError:
            mjd = np.asarray(data['mjd'])
            filenames = np.asarray(data['filename'])
            bad_files = filenames[mjd == None]
            # Print status message
            msg = f'Time invalid for {len(bad_files)} files.\nContinuing, but the following ' \
                  'frames either could not be opened, are empty, or have corrupt headers:\n'
            for file in bad_files:
                msg += f'    {file}\n'
            msgs.warn(msg)

        # Return
        return data

    # TODO:  In this implementation, slicing the PypeItMetaData object
    # will return an astropy.table.Table, not a PypeItMetaData object.
    def __getitem__(self, item):
        return self.table.__getitem__(item)

    def __setitem__(self, item, value):
        return self.table.__setitem__(item, value)

    def __len__(self):
        return self.table.__len__()

    def __repr__(self):
        return self.table._base_repr_(html=False,
                            descr_vals=['PypeItMetaData:\n',
                                        '              spectrograph={0}\n'.format(
                                                                    self.spectrograph.name),
                                        '              length={0}\n'.format(len(self))])

    def _repr_html_(self):
        return self.table._base_repr_(html=True, max_width=-1,
                            descr_vals=['PypeItMetaData: spectrograph={0}, length={1}\n'.format(
                                                    self.spectrograph.name, len(self))])

    @staticmethod
    def default_keys():
        return [ 'directory', 'filename', 'instrume' ]

    def keys(self):
        return self.table.keys()

    def sort(self, col):
        return self.table.sort(col)

    def merge(self, usrdata, match_type=True):
        """
        Use the provided table to supplement or overwrite the metadata.

        If the internal table already contains the column in `usrdata`,
        the function will try to match the data type of the `usrdata`
        column to the existing data type.  If it can't it will just add
        the column anyway, with the type in `usrdata`.  You can avoid
        this step by setting `match_type=False`.

        Args:
            usrdata (:obj:`astropy.table.Table`):
                A user provided set of data used to supplement or
                overwrite metadata read from the file headers.  The
                table must have a `filename` column that is used to
                match to the metadata table generated within PypeIt.
            match_type (:obj:`bool`, optional):
                Attempt to match the data type in `usrdata` to the type
                in the internal table.  See above.

        Raises:
            TypeError: 
                Raised if `usrdata` is not an `astropy.io.table.Table`
            KeyError:
                Raised if `filename` is not a key in the provided table.
        """
        meta_data_model = meta.get_meta_data_model()
        # Check the input
        if not isinstance(usrdata, table.Table):
            raise TypeError('Must provide an astropy.io.table.Table instance.')
        if 'filename' not in usrdata.keys():
            raise KeyError('The user-provided table must have \'filename\' column!')

        # Make sure the data are correctly ordered
        srt = [np.where(f == self.table['filename'])[0][0] for f in usrdata['filename']]

        # Convert types if possible
        existing_keys = list(set(self.table.keys()) & set(usrdata.keys()))
        radec_done = False
        if len(existing_keys) > 0 and match_type:
            for key in existing_keys:
                if len(self.table[key].shape) > 1:  # NOT ALLOWED!!
                    # TODO: This should be converted to an assert statement...
                    raise ValueError('CODING ERROR: Found high-dimensional column.')
                    #embed(header='372 of metadata')
                elif key in meta_data_model.keys(): # Is this meta data??
                    dtype = meta_data_model[key]['dtype']
                else:
                    dtype = self.table[key].dtype
                # Deal with None's properly
                nones = usrdata[key] == 'None'
                usrdata[key][nones] = None
                # Rest
                # Allow for str RA, DEC (backwards compatability)
                if key in ['ra', 'dec'] and not radec_done:
                    ras, decs = meta.convert_radec(usrdata['ra'][~nones].data,
                                                usrdata['dec'][~nones].data)
                    usrdata['ra'][~nones] = ras.astype(dtype)
                    usrdata['dec'][~nones] = decs.astype(dtype)
                    radec_done = True
                else:
                    usrdata[key][~nones] = usrdata[key][~nones].astype(dtype)

        # Include the user data in the table
        for key in usrdata.keys():
            self.table[key] = usrdata[key][srt]

    def remove_rows(self, rows, regroup=False):
        """
        Remove the provided rows from the data table.

        This edits the object directly, nothing is returned.

        Args:
            rows (:obj:`int`, array-like):
                One or more rows that should be *removed* from the datatable.
                This is passed directly to `astropy.table.Table.remove_rows`_;
                see astropy documentation to confirm allowed types.
            regroup (:obj:`bool`, optional):
                If True, reset the setup/configuration, calibration, and
                combination groups.
        """
        self.table.remove_rows(np.atleast_1d(rows))
        if regroup:
            for col in ['setup', 'calib', 'calibbit', 'comb_id', 'bkg_id']:
                if col in self.keys():
                    del self.table[col]
            self.set_configurations()
            self.set_calibration_groups()
            self.set_combination_groups()

    def finalize_usr_build(self, frametype, setup):
        """
        Finalize the build of the table based on user-provided data,
        typically pulled from the PypeIt file.

        This function:
            - sets the frame types based on the provided object
            - sets all the configurations to the provided `setup`
            - assigns all frames to a single calibration group, if the
              'calib' column does not exist
            - if the 'comb_id' column does not exist, this sets the
              combination groups to be either undefined or to be unique
              for each science or standard frame, see
              :func:`set_combination_groups`.

        .. note::
            This should only be run if all files are from a single
            instrument configuration.  :attr:`table` is modified
            in-place.

        See also: :func:`pypeit.pypeitsetup.PypeItSetup.run`.

        .. todo::
            - Why isn't frametype just in the user-provided data?  It
              may be (see get_frame_types) and I'm just not using it...

        Args:
            frametype (:obj:`dict`):
                A dictionary with the types designated by the user.  The
                file name and type are expected to be the key and value
                of the dictionary, respectively.  The number of keys
                therefore *must* match the number of files in
                :attr:`table`.  For frames that have multiple types, the
                types should be provided as a string with
                comma-separated types.
            setup (:obj:`str`):
                If the 'setup' columns does not exist, fill the
                configuration setup columns with this single identifier.
        """
        self.get_frame_types(user=frametype)
        # TODO: Add in a call to clean_configurations? I didn't add it
        # here, because this method is only called for a preconstructed
        # pypeit file, which should nominally follow an execution of
        # pypeit_setup. If the user edits back in a frame that has an
        # invalid key, at least for now the DEIMOS image reader will
        # fault.
        self.set_configurations(fill=setup)
        self.set_calibration_groups(default=True)
        self.set_combination_groups()

    def get_configuration(self, indx, cfg_keys=None, modified=False):
        """
        Return the configuration dictionary for a given frame.

        Args:
            indx (:obj:`int`):
                The index of the table row to use to construct the
                configuration.
            cfg_keys (:obj:`list`, optional):
                The list of metadata keys to use to construct the
                configuration.  If None, the `configuration_keys` of
                :attr:`spectrograph` is used.
            modified (:obj:`bool`, optional):
                Return the configuration as modified by the
                spectrograph-specific
                :func:`~pypeit.spectrographs.spectrograph.Spectrograph.modify_config`.

        Returns:
            dict: A dictionary with the metadata values from the
            selected row.
        """
        _cfg_keys = self.spectrograph.configuration_keys() if cfg_keys is None else cfg_keys
        cfg = {k:self.table[k][indx] for k in _cfg_keys}
        return self.spectrograph.modify_config(self.table[indx], cfg) if modified else cfg

    def construct_obstime(self, row):
        """
        Construct the MJD of when the frame was observed.

        .. todo::
            - Consolidate with :func:`convert_time` ?

        Args:
            row (:obj:`int`):
                The 0-indexed row of the frame.
        
        Returns:
            astropy.time.Time: The MJD of the observation.
        """
        return time.Time(self['mjd'][row], format='mjd')

    def construct_basename(self, row, obstime=None):
        """
        Construct the root name primarily for PypeIt file output.

        Args:
            row (:obj:`int`):
                The 0-indexed row of the frame.
            obstime (:class:`astropy.time.Time`, optional):
                The MJD of the observation.  If None, constructed using
                :func:`construct_obstime`.
        
        Returns:
            str: The root name for file output.
        """
        _obstime = self.construct_obstime(row) if obstime is None else obstime
        tiso = time.Time(_obstime, format='isot')
        dtime = datetime.datetime.strptime(tiso.value, '%Y-%m-%dT%H:%M:%S.%f')
        return '{0}-{1}_{2}_{3}{4}'.format(self['filename'][row].split('.fits')[0],
                                           self['target'][row].replace(" ", ""),
                                           self.spectrograph.camera,
                                           datetime.datetime.strftime(dtime, '%Y%m%dT'),
                                           tiso.value.split("T")[1].replace(':',''))

#    def get_setup(self, row, det=None, config_only=False):
#        """
#        Construct the setup dictionary.
#
#        .. todo::
#            - This is for backwards compatibility, but we should
#              consider reformatting it.  And it may be something to put
#              in the relevant spectrograph class.
#
#        Args:
#            row (:obj:`int`):
#                The 0-indexed row used to construct the setup.
#            det (:obj:`int`, optional):
#                The 1-indexed detector to include.  If None, all
#                detectors are included.
#            config_only (:obj:`bool`, optional):
#                Just return the dictionary with the configuration, don't
#                include the top-level designation of the configuration
#                itself.
#
#        Returns:
#            dict: The pypeit setup dictionary with the default format.
#
#        Raises:
#            PypeItError:
#                Raised if the 'setup' isn't been defined.
#        """
#        if 'setup' not in self.keys():
#            msgs.error('Cannot provide instrument setup without \'setup\' column; '
#                       'run set_configurations.')
#        dispname = 'none' if 'dispname' not in self.keys() else self['dispname'][row]
#        dispangle = 'none' if 'dispangle' not in self.keys() else self['dispangle'][row]
#        dichroic = 'none' if 'dichroic' not in self.keys() else self['dichroic'][row]
#        decker = 'none' if 'decker' not in self.keys() else self['decker'][row]
#        slitwid = 'none' if 'slitwid' not in self.keys() else self['slitwid'][row]
#        slitlen = 'none' if 'slitlen' not in self.keys() else self['slitlen'][row]
#        binning = '1,1' if 'binning' not in self.keys() else self['binning'][row]
#
#        skey = 'Setup {}'.format(self['setup'][row])
#        # Key names *must* match configuration_keys() for spectrographs
#        setup = {skey:
#                    {'--':
#                        {'disperser': {'dispname': dispname, 'dispangle':dispangle},
#                         'dichroic': dichroic,
#                         'slit': {'decker': decker, 'slitwid':slitwid, 'slitlen':slitlen},
#                         'binning': binning,  # PypeIt orientation binning of a science image
#                         }
#                     }
#                 }
#        #_det = np.arange(self.spectrograph.ndet)+1 if det is None else [det]
#        #for d in _det:
#        #    setup[skey][str(d).zfill(2)] \
#        #            = {'binning': binning, 'det': d,
#        #               'namp': self.spectrograph.detector[d-1]['numamplifiers']}
#        return setup[skey] if config_only else setup

    def get_configuration_names(self, ignore=None, return_index=False, configs=None):
        """
        Get the list of the unique configuration names.
        
        This provides just the list of setup identifiers ('A', 'B',
        etc.) and the row index where it first occurs.  This is
        different from :func:`unique_configurations` because the latter
        determines and provides the configurations themselves.

        This is mostly a convenience function for the writing routines.

        Args:
            ignore (:obj:`list`, optional):
                Ignore configurations in the provided list.
            return_index (:obj:`bool`, optional):
                Return row indices with the first occurence of these
                configurations.
            configs (:obj:`str`, :obj:`list`, optional):
                One or more strings used to select the configurations
                to include in the returned objects. If ``'all'``,
                pass back all configurations. Otherwise, only return
                the configurations matched to this provided string or
                list of strings (e.g., ['A','C']).

        Returns:
            :obj:`tuple`, `numpy.ndarray`_: The list of unique setup names.  A
            tuple is returned with a second `numpy.ndarray`_ object providing
            the indices of the first occurrence of these setups, if requested
            (using ``return_index``).

        Raises:
            PypeItError:
                Raised if the 'setup' isn't been defined.
        """
        if 'setup' not in self.keys():
            msgs.error('Cannot get setup names; run set_configurations.')

        # Unique configurations
        # NOTE: This annoyingly returns Column types, not np.arrays!  So need to
        # specify that the data from the column be used.
        setups, indx = np.unique(self['setup'].data, return_index=True)

        if ignore is not None:
            # Remove the selected configurations to ignore
            rm = np.logical_not(np.isin(setups, ignore))
            setups = setups[rm]
            indx = indx[rm]

        # Restrict
        _configs = None if configs is None else np.atleast_1d(configs)
        # TODO: Why do we need to specify 'all' here? Can't `configs is
        # None` mean that you want all the configurations? Or can we
        # make the default 'all'?
        if configs is not None and 'all' not in _configs:
            use = np.isin(setups, _configs)
            setups = setups[use]
            indx = indx[use]

        return (setups, indx) if return_index else setups

    def _get_cfgs(self, copy=False, rm_none=False):
        """
        Convenience method to return :attr:`configs` with possible
        alterations.

        This method *should not* be called by any method outside of
        this class; use :func:`unique_configurations` instead.

        Args:
            copy (:obj:`bool`, optional):
                Return a deep copy of :attr:`configs` instead of the
                object itself.
            rm_none (:obj:`bool`, optional):
                Remove any configurations set to 'None'. If copy is
                True, this is done *after* :attr:`configs` is copied
                to a new dictionary.

        Returns:
            :obj:`dict`: A nested dictionary, one dictionary per
            configuration with the associated metadata for each.
        """
        _cfg = deepcopy(self.configs) if copy else self.configs
        if rm_none and 'None' in _cfg.keys():
            del _cfg['None']
        return _cfg

    @staticmethod
    def maximum_number_of_configurations():
        nascii = len(string.ascii_uppercase)
        return nascii + nascii**2

    @staticmethod
    def configuration_generator(start=0):
        nascii = len(string.ascii_uppercase)
        return (string.ascii_uppercase[i] if i < nascii \
                    else string.ascii_uppercase[int((i-nascii)//nascii)] \
                            + string.ascii_uppercase[int((i-nascii) % nascii)]
                    for i in range(start, nascii + nascii**2))

    @property
    def n_configs(self):
        if self.configs is None:
            msgs.error('Configurations not defined by PypeItMetaData object.  Execute '
                       'unique_configurations first.')
        return len(list(self.configs.keys()))

    def unique_configurations(self, force=False, copy=False, rm_none=False):
        """
        Return the unique instrument configurations.

        If run before the ``'setup'`` column is initialized, this function
        determines the unique instrument configurations by finding
        unique combinations of the items in the metadata table listed by
        the spectrograph ``configuration_keys`` method.

        If run after the ``'setup'`` column has been set, this simply
        constructs the configuration dictionary using the unique
        configurations in that column.

        This is used to set the internal :attr:`configs`. If this
        attribute is not None, this function simply returns
        :attr:`config` (cf. ``force``).

        .. warning::

            Any frame types returned by the
            :func:`~pypeit.spectrographs.spectrograph.Spectrograph.config_independent_frames`
            method for :attr:`spectrograph` will be ignored in the
            construction of the unique configurations. If
            :func:`~pypeit.spectrographs.spectrograph.Spectrograph.config_independent_frames`
            does not return None and the frame types have not yet
            been defined (see :func:`get_frame_types`), this method
            will fault!

        Args:
            force (:obj:`bool`, optional):
                Force the configurations to be redetermined.  Otherwise
                the configurations are only determined if
                :attr:`configs` has not yet been defined.
            copy (:obj:`bool`, optional):
                Return a deep copy of :attr:`configs` instead of the
                object itself.
            rm_none (:obj:`bool`, optional):
                Remove any configurations set to 'None'. If copy is
                True, this is done *after* :attr:`configs` is copied
                to a new dictionary.

        Returns:
            :obj:`dict`: A nested dictionary, one dictionary per
            configuration with the associated metadata for each.

        Raises:
            PypeItError:
                Raised if there are list of frame types to ignore but
                the frame types have not been defined yet.
        """
        if self.configs is not None and not force:
            return self._get_cfgs(copy=copy, rm_none=rm_none)

        if 'setup' in self.keys():
            msgs.info('Setup column already set.  Finding unique configurations.')
            uniq, indx = np.unique(self['setup'], return_index=True)
            ignore = uniq == 'None'
            if np.sum(ignore) > 0:
                msgs.warn(f'Ignoring {np.sum(ignore)} frames with configuration set to None.')
            self.configs = {}
            for i in range(len(uniq)):
                if ignore[i]:
                    continue
                self.configs[uniq[i]] = self.get_configuration(indx[i])
            msgs.info('Found {0} unique configurations.'.format(len(self.configs)))
            return self._get_cfgs(copy=copy, rm_none=rm_none)

        msgs.info('Using metadata to determine unique configurations.')

        # sort self.table
        mjd = self.table['mjd'].copy()
        # Deal with possibly None mjds if there were corrupt header cards
        mjd[mjd == None] = -99999.0
        isort = np.argsort(mjd)
        self.table = self.table[isort]

        # If the frame types have been set, ignore anything listed in
        # the ignore_frames
        indx = np.arange(len(self))
        ignore_frames = self.spectrograph.config_independent_frames()
        if ignore_frames is not None:
            if 'frametype' not in self.keys():
                msgs.error('To ignore frames, types must have been defined; run get_frame_types.')
            ignore_frames = list(ignore_frames.keys())
            msgs.info('Unique configurations ignore frames with type: {0}'.format(ignore_frames))
            use = np.ones(len(self), dtype=bool)
            for ftype in ignore_frames:
                use &= np.logical_not(self.find_frames(ftype))
            indx = indx[use]
        if len(indx) == 0:
            msgs.error('No frames to use to define configurations!')

        # Instantiate the configuration generator
        cfg_gen = PypeItMetaData.configuration_generator()

        # TODO: Placeholder: Allow an empty set of configuration keys
        # meaning that the instrument setup has only one configuration.
        if len(self.spectrograph.configuration_keys()) == 0:
            self.configs = {}
            self.configs[next(cfg_gen)] = {}
            msgs.info('All files assumed to be from a single configuration.')
            return self._get_cfgs(copy=copy, rm_none=rm_none)

        # Use the first file to set the first unique configuration
        self.configs = {}
        self.configs[next(cfg_gen)] = self.get_configuration(indx[0])

        # Check if any of the other files show a different configuration.
        for i in indx[1:]:
            cfg = self.get_configuration(i)
            matched = False
            for c in self.configs.values():
                if self.spectrograph.same_configuration([c,cfg]):
                    matched = True
                    break
            if matched:
                # Matches an existing configuration, so move on to the next
                # frame
                continue

            try:
                # Get the next setup identifier
                setup = next(cfg_gen)
            except StopIteration:
                msgs.error('Cannot assign more configurations!  Either something went wrong'
                           'or you are trying to reduce data from more than '
                           f'{PypeItMetaData.maximum_number_of_configurations()} setups!')
            # Add the configuration
            self.configs[setup] = cfg

        msgs.info(f'Found {len(self.configs)} unique configuration(s).')
        return self._get_cfgs(copy=copy, rm_none=rm_none)

    def set_configurations(self, configs=None, force=False, fill=None):
        """
        Assign each frame to a configuration (setup) and include it
        in the metadata table.

        The internal table is edited *in place*. If the 'setup'
        column already exists, the configurations are **not** reset
        unless you call the function with ``force=True``.

        Args:
            configs (:obj:`dict`, optional):
                A nested dictionary, one dictionary per configuration
                with the associated values of the metadata associated
                with each configuration.  The metadata keywords in the
                dictionary should be the same as in the table, and the
                keywords used to set the configuration should be the
                same as returned by the spectrograph
                `configuration_keys` method.  The latter is not checked.
                If None, this is set by :func:`unique_configurations`. 
            force (:obj:`bool`, optional):
                Force the configurations to be reset.
            fill (:obj:`str`, optional):
                If the 'setup' column does not exist, fill the
                configuration setup columns with this single identifier.
                Ignores other inputs.

        Raises:
            PypeItError:
                Raised if none of the keywords in the provided
                configuration match with the metadata keywords. Also
                raised when some frames cannot be assigned to a
                configuration, the spectrograph defined frames that
                have been ignored in the determination of the unique
                configurations, but the frame types have not been set
                yet.
        """
        # Configurations have already been set
        if 'setup' in self.keys() and not force:
            return

        if 'setup' not in self.keys() and fill is not None:
            self['setup'] = fill
            return

        _configs = self.unique_configurations() if configs is None else configs
        for k, cfg in _configs.items():
            if len(set(cfg.keys()) - set(self.keys())) > 0:
                msgs.error('Configuration {0} defined using unavailable keywords!'.format(k))

        # define the column 'setup' in self.table
        nrows = len(self)
        col = table.Column(data=['None'] * nrows, name='setup', dtype='U25')
        self.table.add_column(col)
        is_science = self.find_frames('science')    # Science frames can only have one configuration
        for i in range(nrows):
            for d, cfg in _configs.items():
                # modify the configuration items only for specific frames. This is instrument dependent.
                mod_cfg = self.spectrograph.modify_config(self.table[i], cfg)
                this_cfg = self.get_configuration(i, modified=True)
                if self.spectrograph.same_configuration([this_cfg, mod_cfg], check_keys=False):
                    if d in self.table['setup'][i]:
                        continue
                    elif self.table['setup'][i] == 'None':
                        self.table['setup'][i] = d
                    elif not is_science[i]:
                        self.table['setup'][i] += f',{d}'

        # Check if any of the configurations are not set
        not_setup = self.table['setup'] == 'None'
        if not np.any(not_setup):
            # All are set, so we're done
            return

        # Some frame types may have been ignored
        ignore_frames = self.spectrograph.config_independent_frames()
        if ignore_frames is None:
            # Nope, we're still done
            return

        # At this point, we need the frame type to continue
        if 'frametype' not in self.keys():
            msgs.error('To account for ignored frames, types must have been defined; run '
                       'get_frame_types.')

        # For each configuration, determine if any of the frames with
        # the ignored frame types should be assigned to it:
        for cfg_key in _configs.keys():
            in_cfg = np.array([cfg_key in _setup for _setup in self.table['setup']])
            for ftype, metakey in ignore_frames.items():

                # TODO: For now, use this assert to check that the
                # metakey is either not set or a string
                assert metakey is None or isinstance(metakey, str), \
                    'CODING ERROR: metadata keywords set by config_indpendent_frames are not ' \
                    'correctly defined for {0}; values must be None or a string.'.format(
                        self.spectrograph.__class__.__name__)

                # Get the list of frames of this type without a
                # configuration
                indx = not_setup & self.find_frames(ftype)
                if not np.any(indx):
                    continue
                if metakey is None:
                    # No matching meta data defined, so just set all
                    # the frames to this (first) configuration
                    self.table['setup'][indx] = cfg_key
                    continue

                # Find the unique values of meta for this configuration
                uniq_meta = np.unique(self.table[metakey][in_cfg].data)
                # Warn the user that the matching meta values are not
                # unique for this configuration.
                if uniq_meta.size != 1:
                    msgs.warn('When setting the instrument configuration for {0} '.format(ftype)
                              + 'frames, configuration {0} does not have unique '.format(cfg_key)
                              + '{0} values.' .format(meta))
                # Find the frames of this type that match any of the
                # meta data values
                indx &= np.isin(self.table[metakey], uniq_meta)
                # assign
                new_cfg_key = np.full(len(self.table['setup'][indx]), 'None', dtype=object)
                for c in range(len(self.table['setup'][indx])):
                    if cfg_key in self.table['setup'][indx][c]:
                        new_cfg_key[c] = self.table['setup'][indx][c]
                    if self.table['setup'][indx][c] == 'None':
                        new_cfg_key[c] = cfg_key
                    else:
                        new_cfg_key[c] = self.table['setup'][indx][c] + ',{}'.format(cfg_key)
                self.table['setup'][indx] = new_cfg_key

    def clean_configurations(self):
        """
        Ensure that configuration-defining keywords all have values
        that will yield good PypeIt reductions. Any frames that do
        not are removed from :attr:`table`, meaning this method may
        modify that attribute directly.

        The valid values for configuration keys is set by
        :func:`~pypeit.spectrographs.spectrograph.Spectrograph.valid_configuration_values`.
        """
        cfg_limits = self.spectrograph.valid_configuration_values()
        if cfg_limits is None:
            # No values specified, so we're done
            return

        good = np.ones(len(self), dtype=bool)
        for key in cfg_limits.keys():
            # NOTE: For now, check that the configuration values were
            # correctly assigned in the spectrograph class definition.
            # This should probably go somewhere else or just removed.
            assert isinstance(cfg_limits[key], list), \
                'CODING ERROR: valid_configuration_values is not correctly defined ' \
                'for {0}; values must be a list.'.format(self.spectrograph.__class__.__name__)

            # Check that the metadata are valid for this column.
            indx = np.isin(self[key], cfg_limits[key])
            if not np.all(indx):
                msgs.warn('Found frames with invalid {0}.'.format(key))
            good &= indx

        if np.all(good):
            # All values good, so we're done
            return

        # Alert the user that some of the frames are going to be
        # removed
        msg = 'The following frames have configurations that cannot be reduced by PypeIt' \
              ' and will be removed from the metadata table (pypeit file):\n'
        indx = np.where(np.logical_not(good))[0]
        for i in indx:
            msg += '    {0}\n'.format(self['filename'][i])
        msgs.warn(msg)
        # And remove 'em
        self.table = self.table[good]

    def find_configuration(self, setup, index=False):
        """
        Find all frames associated with the provided setup/configuration.

        Args:
            setup (:obj:`str`):
                The setup/configuration to search on.
            index (:obj:`bool`, optional):
                Return an array of 0-indexed indices instead of a
                boolean array.

        Returns:
            numpy.ndarray: A boolean array, or an integer array if
            ``index=True``, with the table rows associated with the requested
            setup/configuration.
        """
        if 'setup' not in self.keys():
            msgs.error('Configurations not set; first execute self.unique_configurations.')

        # NOTE: frames can be associated with multiple setups (namely biases),
        # meaning that we have to split the string by any separating commas.
        # This has the added benefit that ``'A' in 'AA'`` is True (a bug in our
        # setup naming system in the albeit crazy case when someone is trying to
        # parse more than 26 setups), but ``'A' in 'AA'.split(',')`` is False.
        in_cfg = np.array([setup in _set.split(',') for _set in self.table['setup']])
        return np.where(in_cfg)[0] if index else in_cfg

    # TODO: Make some tests for this!
    def _set_calib_group_bits(self):
        """
        Set the calibration group bit based on the string values of the
        'calib' column.
        """
        # NOTE: This is a hack to ensure the type of the *elements* of the calib
        # column are all strings, but that the type of the column remains as
        # "object".  I'm calling this a hack because doing this is easier than
        # trying to track down everywhere calib is changed to values that may or
        # may not be integers instead of strings.
        self['calib'] = np.array([str(c) for c in self['calib']], dtype=object)
        # Collect and expand any lists
        group_names = np.unique(np.concatenate(
                        [s.split(',') for s in self['calib'] if s not in ['all', 'None']]))
        # Expand any ranges
        keep_group = np.ones(group_names.size, dtype=bool)
        added_groups = []
        for i, name in enumerate(group_names):
            if ':' in name:
                # Parse the range
                keep_group[i] = False
                added_groups += [str(n) for n in parse.str2list(name)]
        # Combine and find the unique *integer* identifiers
        group_names = np.unique(np.asarray(added_groups + 
                                            (group_names[keep_group]).tolist()).astype(int))

        # Define the bitmask and initialize the bits
        self.calib_bitmask = BitMask(group_names)
        self['calibbit'] = 0

        # Set the bits based on the selected groups
        for i in range(len(self)):
            if self['calib'][i] == 'None':
                continue
            if self['calib'][i] == 'all':
                self['calibbit'][i] = self.calib_bitmask.turn_on(self['calibbit'][i], group_names)
                continue
            self['calibbit'][i] = self.calib_bitmask.turn_on(self['calibbit'][i], 
                                                             parse.str2list(self['calib'][i]))

    def _check_calib_groups(self):
        """
        Check that the calibration groups are valid.

        This currently only checks that the science frames are
        associated with one calibration group.

        TODO: Is this appropriate for NIR data?

        """
        is_science = self.find_frames('science')
        for i in range(len(self)):
            if not is_science[i]:
                continue
            if len(self.calib_bitmask.flagged_bits(self['calibbit'][i])) > 1:
                msgs.error('Science frames can only be assigned to a single calibration group.')

    @property
    def n_calib_groups(self):
        """Return the number of calibration groups."""
        return None if self.calib_bitmask is None else self.calib_bitmask.nbits
                
    @property
    def calib_groups(self):
        """Return the calibration group identifiers."""
        return None if self.calib_bitmask is None else self.calib_bitmask.keys()
                
    def set_calibration_groups(self, global_frames=None, default=False, force=False):
        """
        Group calibration frames into sets.
        
        Requires the 'setup' column to have been defined.  For now this
        is a simple grouping of frames with the same configuration.

        .. todo::
            - Maintain a detailed description of the logic.

        The 'calib' column has a string type to make sure that it
        matches with what can be read from the pypeit file.  The
        'calibbit' column is actually what is used to determine the
        calibration group of each frame; see :attr:`calib_bitmask`.

        Args:
            global_frames (:obj:`list`, optional):
                A list of strings with the frame types to use in all
                calibration groups (e.g., ['bias', 'dark']).
            default (:obj:`bool`, optional):
                If the 'calib' column is not present, set a single
                calibration group *for all rows*.
            force (:obj:`bool`, optional):
                Force the calibration groups to be reconstructed if
                the 'calib' column already exists.

        Raises:
            PypeItError:
                Raised if 'setup' column is not defined, or if
                `global_frames` is provided but the frame types have not
                been defined yet.
        """
        # Set the default if requested and 'calib' doesn't exist yet
        if 'calib' not in self.keys() and default:
            self['calib'] = '0'
            # Make sure the calibbit column does not exist
            if 'calibbit' in self.keys():
                del self['calibbit']

        # Groups have already been set
        if 'calib' in self.keys() and 'calibbit' in self.keys() and not force:
            return

        # Groups have been set but the bits have not (likely because the
        # data was read from a pypeit file)
        if 'calib' in self.keys() and 'calibbit' not in self.keys() and not force:
            self._set_calib_group_bits()
            self._check_calib_groups()
            return

        # TODO: The rest of this just nominally sets the calibration
        # group based on the configuration.  This will change!

        # The configuration must be present to determine the calibration
        # group
        if 'setup' not in self.keys():
            msgs.error('Must have defined \'setup\' column first; try running set_configurations.')
        configs = np.unique(np.concatenate([_setup.split(',') for _setup in self['setup'].data])).tolist()
        if 'None' in configs:
            configs.remove('None')      # Ignore frames with undefined configurations
        n_cfg = len(configs)

        # TODO: Science frames can only have one calibration group

        # Assign everything from the same configuration to the same
        # calibration group; this needs to have dtype=object, otherwise
        # any changes to the strings will be truncated at 4 characters.
        self.table['calib'] = np.full(len(self), 'None', dtype=object)
        for i in range(n_cfg):
            in_cfg = np.array([configs[i] in _set for _set in self.table['setup']]) & (self['framebit'] > 0)
            icalibs = np.full(len(self['calib'][in_cfg]), 'None', dtype=object)
            for c in range(len(self['calib'][in_cfg])):
                if self['calib'][in_cfg][c] == 'None':
                    icalibs[c] = str(i)
                else:
                    icalibs[c] = self['calib'][in_cfg][c] + ',{}'.format(i)
            self['calib'][in_cfg] = icalibs

        # Allow some frame types to be used in all calibration groups
        # (like biases and darks)
        if global_frames is not None:
            if 'frametype' not in self.keys():
                msgs.error('To set global frames, types must have been defined; '
                           'run get_frame_types.')

            calibs = '0' if n_cfg == 1 else ','.join(np.arange(n_cfg).astype(str))
            for ftype in global_frames:
                indx = np.where(self.find_frames(ftype))[0]
                for i in indx:
                    self['calib'][i] = calibs

        # Set the bits based on the string representation of the groups
        self._set_calib_group_bits()
        # Check that the groups are valid
        self._check_calib_groups()

    def find_frames(self, ftype, calib_ID=None, index=False):
        """
        Find the rows with the associated frame type.

        If the index is provided, the frames must also be matched to the
        relevant science frame.

        Args:
            ftype (str):
                The frame type identifier.  See the keys for
                :class:`pypeit.core.framematch.FrameTypeBitMask`.  If
                set to the string 'None', this returns all frames
                without a known type.
            calib_ID (:obj:`int`, optional):
                Index of the calibration group that it must match.  If None,
                any row of the specified frame type is included.
            index (:obj:`bool`, optional):
                Return an array of 0-indexed indices instead of a
                boolean array.

        Returns:
            numpy.ndarray: A boolean array, or an integer array if
            index=True, with the rows that contain the frames of the
            requested type.  

        Raises:
            PypeItError:
                Raised if the `framebit` column is not set in the table.
        """
        if 'framebit' not in self.keys():
            msgs.error('Frame types are not set.  First run get_frame_types.')
        if ftype == 'None':
            return self['framebit'] == 0
        # Select frames
        indx = self.type_bitmask.flagged(self['framebit'], ftype)

        if calib_ID is not None:
            # Select frames in the same calibration group
            indx &= self.find_calib_group(calib_ID)

        # Return
        return np.where(indx)[0] if index else indx

    def find_frame_files(self, ftype, calib_ID=None):
        """
        Return the list of files with a given frame type.

        The frames must also match the science frame index, if it is
        provided.

        Args:
            ftype (str):
                The frame type identifier.  See the keys for
                :class:`pypeit.core.framematch.FrameTypeBitMask`.
            calib_ID (:obj:`int`, optional):
                Index of the calibration group that it must match.  If None,
                any row of the specified frame type is included.

        Returns:
            list: List of file paths that match the frame type and
            science frame ID, if the latter is provided.
        """
        return self.frame_paths(self.find_frames(ftype, calib_ID=calib_ID))

    def frame_paths(self, indx):
        """
        Return the full paths to one or more frames.

        Args:
            indx (:obj:`int`, array-like):
                One or more 0-indexed rows in the table with the frames
                to return.  Can be an array of indices or a boolean
                array of the correct length.
        Returns:
            list: List of the full paths of one or more frames.
        """
        if isinstance(indx, (int,np.integer)):
            return os.path.join(self['directory'][indx], self['filename'][indx])
        return [os.path.join(d,f) for d,f in zip(self['directory'][indx], self['filename'][indx])]

    def set_frame_types(self, type_bits, merge=True):
        """
        Set and return a Table with the frame types and bits.
        
        Args:
            type_bits (numpy.ndarray):
                Integer bitmask with the frame types.  The length must
                match the existing number of table rows.

            merge (:obj:`bool`, optional):
                Merge the types and bits into the existing table.  This
                will *overwrite* any existing columns.
        
        Returns:
            `astropy.table.Table`: Table with two columns, the frame
            type name and bits.
        """
        # Making Columns to pad string array
        ftype_colmA = table.Column(self.type_bitmask.type_names(type_bits), name='frametype')

        # KLUDGE ME
        #
        # TODO: It would be good to get around this.  Is it related to
        # this change?
        # http://docs.astropy.org/en/stable/table/access_table.html#bytestring-columns-in-python-3
        #
        # See also:
        #
        # http://docs.astropy.org/en/stable/api/astropy.table.Table.html#astropy.table.Table.convert_bytestring_to_unicode
        #
        # Or we can force type_names() in bitmask to always return the
        # correct type...
        if int(str(ftype_colmA.dtype)[2:]) < 9:
            ftype_colm = table.Column(self.type_bitmask.type_names(type_bits), dtype='U9',
                                      name='frametype')
        else:
            ftype_colm = ftype_colmA

        fbits_colm = table.Column(type_bits, name='framebit')
        t = table.Table([ftype_colm, fbits_colm])

        if merge:
            self['frametype'] = t['frametype']
            self['framebit'] = t['framebit']
        return t

    def edit_frame_type(self, indx, frame_type, append=False):
        """
        Edit the frame type by hand.

        Args:
            indx (:obj:`int`):
                The 0-indexed row in the table to edit
            frame_type (:obj:`str`, :obj:`list`):
                One or more frame types to append/overwrite.
            append (:obj:`bool`, optional):
                Append the frame type.  If False, all existing frame
                types are overwitten by the provided type.
        """
        if not append:
            self['framebit'][indx] = 0
        self['framebit'][indx] = self.type_bitmask.turn_on(self['framebit'][indx], flag=frame_type)
        self['frametype'][indx] = self.type_bitmask.type_names(self['framebit'][indx])

    def get_frame_types(self, flag_unknown=False, user=None, merge=True):
        """
        Generate a table of frame types from the input metadata object.

        .. todo::
            - Here's where we could add a SPIT option.
    
        Args:
            flag_unknown (:obj:`bool`, optional):
                Instead of crashing out if there are unidentified files,
                leave without a type and continue.
            user (:obj:`dict`, optional):
                A dictionary with the types designated by the user.  The
                file name and type are expected to be the key and value
                of the dictionary, respectively.  The number of keys
                therefore *must* match the number of files in
                :attr:`table`.  For frames that have multiple types, the
                types should be provided as a string with
                comma-separated types.
            merge (:obj:`bool`, optional):
                Merge the frame typing into the exiting table.

        Returns:
            :obj:`astropy.table.Table`: A Table with two columns, the
            type names and the type bits.  See
            :class:`pypeit.core.framematch.FrameTypeBitMask` for the
            allowed frame types.
        """
        # Checks
        if 'frametype' in self.keys() or 'framebit' in self.keys():
            msgs.warn('Removing existing frametype and framebit columns.')
        if 'frametype' in self.keys():
            del self.table['frametype']
        if 'framebit' in self.keys():
            del self.table['framebit']

        # Start
        msgs.info("Typing files")
        type_bits = np.zeros(len(self), dtype=self.type_bitmask.minimum_dtype())
    
        # Use the user-defined frame types from the input dictionary
        if user is not None:
            if len(user.keys()) != len(self):
                if len(np.unique(self['filename'].data)) != len(self):
                    msgs.error('Your pypeit file has duplicate filenames which is not allowed.')
                else:
                    msgs.error('The user-provided dictionary does not match table length.')
            msgs.info('Using user-provided frame types.')
            for ifile,ftypes in user.items():
                indx = self['filename'] == ifile
                try:
                    type_bits[indx] = self.type_bitmask.turn_on(type_bits[indx], flag=ftypes.split(','))
                except ValueError as err:
                    msgs.error(f'Improper frame type supplied!{msgs.newline()}'
                               f'{err}{msgs.newline()}'
                               'Check your PypeIt Reduction File')
            return self.set_frame_types(type_bits, merge=merge)
    
        # Loop over the frame types
        for i, ftype in enumerate(self.type_bitmask.keys()):
    
            # Include a combination of instrument-specific checks using
            # combinations of the full set of metadata
            exprng = self.par['scienceframe']['exprng'] if ftype == 'science' \
                else self.par['calibrations']['{0}frame'.format(ftype)]['exprng']
            # TODO: Use & or | ?  Using idname above gets overwritten by
            # this if the frames to meet the other checks in this call.
#            indx &= self.spectrograph.check_frame_type(ftype, self.table, exprng=exprng)
            indx = self.spectrograph.check_frame_type(ftype, self.table, exprng=exprng)
            # Turn on the relevant bits
            type_bits[indx] = self.type_bitmask.turn_on(type_bits[indx], flag=ftype)
    
        # Find the nearest standard star to each science frame
        # TODO: Should this be 'standard' or 'science' or both?
        if 'ra' not in self.keys() or 'dec' not in self.keys():
            msgs.warn('Cannot associate standard with science frames without sky coordinates.')
        else:
            # TODO: Do we want to do this here?
            indx = self.type_bitmask.flagged(type_bits, flag='standard')
            for b, f, ra, dec in zip(type_bits[indx], self['filename'][indx], self['ra'][indx],
                                     self['dec'][indx]):
                if ra == 'None' or dec == 'None':
                    msgs.warn('RA and DEC must not be None for file:' + msgs.newline() + f)
                    msgs.warn('The above file could be a twilight flat frame that was'
                              + msgs.newline() + 'missed by the automatic identification.')
                    b = self.type_bitmask.turn_off(b, flag='standard')
                    continue

                # If an object exists within 20 arcmins of a listed standard,
                # then it is probably a standard star
                foundstd = flux_calib.find_standard_file(ra, dec, check=True)
                b = self.type_bitmask.turn_off(b, flag='science' if foundstd else 'standard')
    
        # Find the files without any types
        indx = np.logical_not(self.type_bitmask.flagged(type_bits))
        if np.any(indx):
            msgs.info("Couldn't identify the following files:")
            for f in self['filename'][indx]:
                msgs.info(f)
            if not flag_unknown:
                msgs.error("Check these files before continuing")
    
        # Finish up (note that this is called above if user is not None!)
        msgs.info("Typing completed!")
        return self.set_frame_types(type_bits, merge=merge)

    def set_pypeit_cols(self, write_bkg_pairs=False, write_manual=False):
        """
        Generate the list of columns to be included in the fitstbl
        (nearly the complete list).

        Args:
            write_bkg_pairs (:obj:`bool`, optional):
                Add additional ``PypeIt`` columns for calib, comb_id
                and bkg_id
            write_manual (:obj:`bool`, optional):
                Add additional ``PypeIt`` columns for manual extraction


        Returns:
            `numpy.ndarray`_: Array of columns to be used in the fits
            table>
        """
        # Columns for output
        columns = self.spectrograph.pypeit_file_keys() + ['calib']

        extras = []

        # comb, bkg columns
        if write_bkg_pairs:
            extras += ['comb_id', 'bkg_id']
        # manual
        if write_manual:
            extras += ['manual']
        for key in extras:
            if key not in columns:
                columns += [key]

        # Take only those present
        output_cols = np.array(columns)
        return output_cols[np.isin(output_cols, self.keys())].tolist()

    def set_combination_groups(self, assign_objects=True):
        """
        Set combination groups.

        .. note::
            :attr:`table` is edited in place.

        This function can be used to initialize the combination group
        and background group columns, and/or to initialize the combination
        groups to the set of objects (science or standard frames) to a
        unique integer.

        If the 'comb_id' or 'bkg_id' columns do not exist, they're set
        to -1.  

        Args:
            assign_objects (:obj:`bool`, optional):
                If all of 'comb_id' values are less than 0 (meaning they're
                unassigned), the combination groups are set to be unique for
                each standard and science frame.  For some instruments (e.g.,
                Keck/NIRES), this will also parse known dither patterns and use
                them to set default difference-imaging groups.

        """
        if 'comb_id' not in self.keys():
            self['comb_id'] = -1
        if 'bkg_id' not in self.keys():
            self['bkg_id'] = -1

        # NOTE: Importantly, this if statement means that, if the user has
        # defined any non-negative combination IDs in their pypeit file, none of
        # this automated assignment logic is executed.
        if assign_objects and np.all(self['comb_id'] < 0):
            # find_frames will throw an exception if framebit is not
            # set...
            sci_std_idx = np.where(np.any([self.find_frames('science'),
                                           self.find_frames('standard')], axis=0))[0]
            self['comb_id'][sci_std_idx] = np.arange(len(sci_std_idx), dtype=int) + 1
            # update comb_id and bkg_id with dither pattern specific to the used instrument
            # if get_comb_group() is not defined in the relevant spectrograph self.table is unchanged
            self.table = self.spectrograph.get_comb_group(self.table)

            if 'calib' in self.keys():
                # Re-set the calibbit in case calib was changed by get_comb_group().
                # If calib was not changed, calibbit values will be unchanged
                self._set_calib_group_bits()
                # Check that the groups are valid
                self._check_calib_groups()

    def set_user_added_columns(self):
        """
        Set columns that the user *might* add

        .. note::
            :attr:`table` is edited in place.

        This function can be used to initialize columns
        that the user might add

        """
        if 'manual' not in self.keys():
            self['manual'] = ''

    def write_sorted(self, ofile, overwrite=True, ignore=None, 
                     write_bkg_pairs=False, write_manual=False):
        """
        Write the sorted file.

        The sorted file lists all the unique instrument configurations
        (setups) and the frames associated with each configuration.  The
        output data table is identical to the pypeit file output.

        .. todo::
            - This is for backwards compatibility, but we should
              consider reformatting/removing it.

        Args:
            ofile (:obj:`str`, `Path`_):
                Name for the output sorted file.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file with the same name.
            ignore (:obj:`list`, optional):
                Ignore configurations in the provided list.
            write_bkg_pairs (:obj:`bool`, optional):
                Add additional ``PypeIt`` columns for calib, comb_id
                and bkg_id
            write_manual (:obj:`bool`, optional):
                Add additional ``PypeIt`` columns for manual extraction

        Raises:
            PypeItError:
                Raised if the 'setup' isn't been defined.
        """
        if 'setup' not in self.keys():
            msgs.error('Cannot write sorted instrument configuration table without \'setup\' '
                       'column; run set_configurations.')

        _ofile = Path(ofile).resolve()    
        if _ofile.exists() and not overwrite:
            msgs.error(f'{_ofile} already exists.  Use ovewrite=True to overwrite.')

        # Grab output columns
        output_cols = self.set_pypeit_cols(write_bkg_pairs=write_bkg_pairs,
                                           write_manual=write_manual)

        cfgs = self.unique_configurations(copy=ignore is not None)
        if ignore is not None:
            for key in cfgs.keys():
                if key in ignore:
                    del cfgs[key]

        # Construct file
        ff = open(_ofile, 'w')
        for setup in cfgs.keys():
            # Get the subtable of frames taken in this configuration
            indx = np.array([setup in _set for _set in self['setup']])
            if not np.any(indx):
                continue
            subtbl = self.table[output_cols][indx]
            if 'calib' in output_cols:
                # calib can be a str with a list of values because in some cases (e.g. MOSFIRE) the same
                # calibration files are used for different setups. Here we update calib to have only the
                # value relevant for this setup.
                # find the calib value in this setup that is not a list (which is probably a science/standard)
                no_list = np.array([',' not in str(cc) for cc in subtbl['calib']])
                if np.any(no_list):
                    # assign the calib value in this setup that is not a list to frames that have calib as a list
                    subtbl['calib'][np.logical_not(no_list)] = subtbl['calib'][no_list][0]

            # Write the file
            ff.write('##########################################################\n')
            ff.write('Setup {:s}\n'.format(setup))
            ff.write('\n'.join(dict_to_lines(cfgs[setup], level=1)) + '\n')
            ff.write('#---------------------------------------------------------\n')
            mjd = subtbl['mjd'].copy()
            # Deal with possibly None mjds if there were corrupt header cards
            mjd[mjd == None] = -99999.0
            isort = np.argsort(mjd)
            subtbl = subtbl[isort]
            # This needs to match the format for writing file blocks in pypeit.inputfiles.InputFile
            subtbl.write(ff, format='ascii.fixed_width', bookend=False)
        ff.write('##end\n')
        ff.close()

    def write_pypeit(self, output_path=None, cfg_lines=None,
                     write_bkg_pairs=False, write_manual=False,
                     configs=None, config_subdir=True,
                     version_override=None, date_override=None):
        """
        Write a pypeit file in data-table format.

        The pypeit file is the main configuration file for PypeIt,
        configuring the control-flow and algorithmic parameters and
        listing the data files to read.  This function writes the
        columns selected by the
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.pypeit_file_keys`,
        which can be specific to each instrument.

        Args:
            output_path (:obj:`str`, optional):
                Root path for the output pypeit files. If None, set
                to current directory. If the output directory does
                not exist, it is created.
            cfg_lines (:obj:`list`, optional):
                The list of configuration lines to include in the file.
                If None are provided, the vanilla configuration is
                included.
            write_bkg_pairs (:obj:`bool`, optional):
                When constructing the
                :class:`pypeit.metadata.PypeItMetaData` object, include
                two columns called `comb_id` and `bkg_id` that identify
                object and background frame pairs.  
            write_manual (:obj:`bool`, optional):
                Add additional ``PypeIt`` columns for manual extraction
            configs (:obj:`str`, :obj:`list`, optional):
                One or more strings used to select the configurations
                to include in the returned objects. If ``'all'``,
                pass back all configurations. Otherwise, only return
                the configurations matched to this provided string or
                list of strings (e.g., ['A','C']). See
                :attr:`configs`.
            config_subdir (:obj:`bool`, optional):
                Flag to place the pypeit file in a subdirectory named for each
                configuration.  If True, the pypeit file is written to
                ``{spec}_{config}/{spec}_{config}.pypeit`` (e.g.,
                ``shane_kast_blue_A/shane_kast_blue_A.pypeit``).  If False, the
                pypeit file is placed directly in the ``output_path``.
            version_override (:obj:`str`, optional):
                Override the current version and use this one instead.  **For
                documentation purposes only!**
            date_override (:obj:`str`, optional):
                Override the current date and use this one instead.  **For
                documentation purposes only!**

        Raises:
            PypeItError:
                Raised if the 'setup' isn't defined and split is True.

        Returns:
            :obj:`list`: List of ``PypeIt`` files generated.
        """
        # Set output path
        if output_path is None:
            output_path = os.getcwd()

        # Find unique configurations, always ignoring any 'None'
        # configurations...
        cfg = self.unique_configurations(copy=True, rm_none=True)

        # Get the setups to write
        if configs is None or configs == 'all' or configs == ['all']:
            cfg_keys = list(cfg.keys())
        else:
            _configs = configs if isinstance(configs, list) else [configs]
            cfg_keys = [key for key in cfg.keys() if key in _configs]

        if len(cfg_keys) == 0:
            msgs.error('No setups to write!')

        # Grab output columns
        output_cols = self.set_pypeit_cols(write_bkg_pairs=write_bkg_pairs,
                                           write_manual=write_manual)

        # Write the pypeit files
        ofiles = [None]*len(cfg_keys)
        for j,setup in enumerate(cfg_keys):
            # Create the output directory
            root = '{0}_{1}'.format(self.spectrograph.name, setup)
            if config_subdir:
                odir = os.path.join(output_path, root)
                if not os.path.isdir(odir):
                    os.makedirs(odir)
            else:
                odir = output_path
            # Create the output file name
            ofiles[j] = os.path.join(odir, '{0}.pypeit'.format(root))

            # Setup dict
            setup_dict = {}
            setup_dict[f'Setup {setup}'] = {}
            for key in cfg[setup]:
                setup_dict[f'Setup {setup}'][key] = cfg[setup][key]
            
            # Get the paths
            in_cfg = np.array([setup in _set for _set in self.table['setup']])
            if not np.any(in_cfg):
                continue
            paths = np.unique(self['directory'][in_cfg]).tolist()

            # Get the data lines
            subtbl = self.table[output_cols][in_cfg]
            if 'calib' in output_cols:
                # calib can be a str with a list of values because in some cases (e.g. MOSFIRE) the same
                # calibration files are used for different setups. Here we update calib to have only the
                # value relevant for this setup.
                # find the calib value in this setup that is not a list (which is probably a science/standard)
                no_list = np.array([',' not in str(cc) for cc in subtbl['calib']])
                if np.any(no_list):
                    # assign the calib value in this setup that is not a list to frames that have calib as a list
                    subtbl['calib'][np.logical_not(no_list)] = subtbl['calib'][no_list][0]
            subtbl.sort(['frametype','filename'])
            #with io.StringIO() as ff:
            #    subtbl.write(ff, format='ascii.fixed_width')
            #    data_lines = ff.getvalue().split('\n')[:-1]

            # Config lines
            if cfg_lines is None:
                cfg_lines = ['[rdx]']
                cfg_lines += ['    spectrograph = {0}'.format(self.spectrograph.name)]

            # Instantiate a PypeItFile
            pypeItFile = inputfiles.PypeItFile(cfg_lines, paths, subtbl, setup_dict)
            # Write
            pypeItFile.write(ofiles[j], version_override=version_override,
                             date_override=date_override) 

        # Return
        return ofiles

    def write(self, output=None, rows=None, columns=None, sort_col=None, overwrite=False,
              header=None):
        """
        Write the metadata either to a file or to the screen.

        The method allows you to set the columns to print and which column to
        use for sorting.

        Args:
            output (:obj:`str`, optional):
                Output signature or file name. If None, the table contents
                are printed to the screen. If ``'table'``, the table that
                would have been printed/written to disk is returned.
                Otherwise, the string is interpreted as the name of an ascii
                file to which to write the table contents.
            rows (`numpy.ndarray`_, optional):
                A boolean vector selecting the rows of the table to write. If
                None, all rows are written. Shape must match the number of
                the rows in the table.
            columns (:obj:`str`, :obj:`list`, optional):
                A list of columns to include in the output file. Can be
                provided as a list directly or as a comma-separated string.
                If None or ``'all'``, all columns in are written; if
                ``'pypeit'``, the columns are the same as those included in
                the pypeit file. Each selected column must be a valid pypeit
                metadata keyword, specific to :attr:`spectrograph`.
                Additional valid keywords, depending on the processing level
                of the metadata table, are directory, filename, frametype,
                framebit, setup, calib, and calibbit.
            sort_col (:obj:`str`, optional):
                Name of the column to use for sorting the output. If
                None, the table is printed in its current state.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file; otherwise raise an
                exception.
            header (:obj:`str`, :obj:`list`, optional):
                One or more strings to write to the top of the file, on
                string per file line; ``# `` is added to the beginning of
                each string. Ignored if ``output`` does not specify an output
                file.

        Returns:
            `astropy.table.Table`: The table object that would have been
            written/printed if ``output == 'table'``. Otherwise, the method
            always returns None.

        Raises:
            ValueError:
                Raised if the columns to include are not valid, or if the
                column to use for sorting is not valid.
            FileExistsError:
                Raised if overwrite is False and the file exists.
        """
        # Check the file can be written (this is here because the spectrograph
        # needs to be defined first)
        ofile = None if output in [None, 'table'] else output
        if ofile is not None and os.path.isfile(ofile) and not overwrite:
            raise FileExistsError(f'{ofile} already exists; set flag to overwrite.')

        # Check the rows input
        if rows is not None and len(rows) != len(self.table):
            raise ValueError('Boolean vector selecting output rows has incorrect length.')

        # Get the columns to return
        if columns in [None, 'all']:
            tbl_cols = list(self.keys())
        elif columns == 'pypeit':
            tbl_cols = self.set_pypeit_cols(write_bkg_pairs=True)
        else:
            all_cols = list(self.keys())
            tbl_cols = columns if isinstance(columns, list) else columns.split(',')
            badcol = [col not in all_cols for col in tbl_cols]
            if np.any(badcol):
                raise ValueError('The following columns are not valid: {0}'.format(
                                 ', '.join(tbl_cols[badcol])))

        # Make sure the basic parameters are the first few columns; do them in
        # reverse order so I can always insert at the beginning of the list
        for col in ['framebit', 'frametype', 'filename', 'directory']:
            if col not in tbl_cols:
                continue
            indx = np.where([t == col for t in tbl_cols])[0][0]
            if indx != 0:
                tbl_cols.insert(0, tbl_cols.pop(indx))

        # Make sure the dithers and combination and background IDs are the last
        # few columns
        ncol = len(tbl_cols)
        for col in ['dithpat', 'dithpos', 'dithoff', 'calib', 'comb_id', 'bkg_id']:
            if col not in tbl_cols:
                continue
            indx = np.where([t == col for t in tbl_cols])[0][0]
            if indx != ncol-1:
                tbl_cols.insert(ncol-1, tbl_cols.pop(indx))

        # Copy the internal table so that it is unaltered
        output_tbl = self.table.copy()

        # Select the output rows if a vector was provided
        if rows is not None:
            output_tbl = output_tbl[rows]

        # Select and sort the data by a given column
        if sort_col is not None:
            if sort_col not in self.keys():
                raise ValueError(f'Cannot sort by {sort_col}.  Not a valid column.')
            # Ignore any NoneTypes
            indx = output_tbl[sort_col] != None
            is_None = np.logical_not(indx)
            srt = np.append(np.where(is_None)[0],
                            np.where(indx)[0][np.argsort(output_tbl[sort_col][indx].data)])
            output_tbl = output_tbl[tbl_cols][srt]
        else:
            output_tbl = output_tbl[tbl_cols]

        if output == 'table':
            # Instead of writing, just return the modified table
            return output_tbl

        # Always write the table in ascii format
        with io.StringIO() as ff:
            output_tbl.write(ff, format='ascii.fixed_width')
            data_lines = ff.getvalue().split('\n')[:-1]

        if ofile is None:
            # Output file not defined so just print it
            print('\n'.join(data_lines))
            return None

        # Write the output to an ascii file
        with open(ofile, 'w') as f:
            if header is not None:
                _header = header if isinstance(header, list) else [header]
                for h in _header:
                    f.write(f'# {h}\n')
            f.write('\n')
            f.write('\n'.join(data_lines))
            f.write('\n')

        # Just to be explicit that the method returns None when writing to a
        # file...
        return None

    def find_calib_group(self, grp):
        """
        Find all the frames associated with the provided calibration group.
        
        Args:
            grp (:obj:`int`):
                The calibration group integer.

        Returns:
            numpy.ndarray: Boolean array selecting those frames in the
            table included in the selected calibration group.

        Raises:
            PypeItError:
                Raised if the 'calibbit' column is not defined.
        """
        if 'calibbit' not in self.keys():
            msgs.error('Calibration groups are not set.  First run set_calibration_groups.')
        return self.calib_bitmask.flagged(self['calibbit'].data, grp)

    def find_frame_calib_groups(self, row):
        """
        Find the calibration groups associated with a specific frame.
        """
        return self.calib_bitmask.flagged_bits(self['calibbit'][row])


