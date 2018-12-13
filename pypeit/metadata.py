"""
Provides a class that handles the fits metadata required by PypeIt.
"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from __future__ import unicode_literals

import os
import io
import string

import numpy as np
import yaml

import datetime
from astropy import table, coordinates, time

from pypeit import msgs
from pypeit import utils
from pypeit.core import framematch
from pypeit.core import flux
from pypeit.par import PypeItPar
from pypeit.par.util import make_pypeit_file
from pypeit.bitmask import BitMask
from pypeit.spectrographs.util import load_spectrograph

from pypeit import debugger

# Initially tried to subclass this from astropy.table.Table, but that
# proved too difficult.
class PypeItMetaData:
    """
    Provides a table and interface to the relevant fits file metadata
    used during the reduction.

    The content of the fits table is dictated by the header keywords
    specified for the provided spectrograph.  It is expected that this
    table can be used to set the frame type of each file.

    The metadata is validated using checks specified by the provided
    spectrograph class.

    For the data table, one should typically provided either the file
    list from which to grab the data from the fits headers or the data
    directly.  If neither are provided the table is instantiated without
    any data.

    .. todo::
        - Implement additional keywords for `bkg_pairs` that provide
          some nominal pairings that the user can edit.

    Args:
        spectrograph (:obj:`str`,
            :class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph used to collect the data save to each file.
            The class is used to provide the header keyword data to
            include in the table and specify any validation checks.
        par (:obj:`pypeit.par.pypeitpar.PypeItPar`, optional):
            PypeIt parameters used to set the code behavior.  If not
            provided, the default parameters specific to the provided
            spectrograph are used.
        file_list (:obj:`list`, optional):
            The list of files to include in the table.
        data (table-like, optional):
            The data to incude in the table.  The type can be anything
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
            Function will fault if :func:`fits.getheader` fails to read
            any of the headers in the provided file list.  Set to False
            to instead report a warning and continue.

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
    def __init__(self, spectrograph, par=None, file_list=None, data=None, usrdata=None, strict=True):
        if data is None and file_list is None:
            msgs.warn('Both data and file_list are None in the instantiation of PypeItMetaData.'
                      '  The table will be empty!')
        self.spectrograph = load_spectrograph(spectrograph)
        self.par = self.spectrograph.default_pypeit_par() if par is None else par
        if not isinstance(self.par, PypeItPar):
            raise TypeError('Input parameter set must be of type PypeItPar.')
        self.type_bitmask = framematch.FrameTypeBitMask()
        self.table = table.Table(data if file_list is None 
                                 else self._new_build(file_list, strict=strict))
                        #else self._build(file_list, strict=strict))
        if usrdata is not None:
            self.merge(usrdata)
        # Instrument-specific validation of the header metadata. This
        # alters self.table in place!
        self.spectrograph.validate_metadata(self.table)

        # Initialize internal attributes
        self.configs = None
        self.calib_bitmask = None
    
    def _add_bkg_pairs(self, data, bkg_pairs):
        """
        Add the background-pair columns to the table.
        
        This include two columns called `comb_id` and `bkg_id` that
        identify object and background frame pairs.  The string
        indicates how these these columns should be added::

            - `empty`: The columns are added but their values are all
              originally set to -1.  **This is currently the only
              option.**

        Args:
            data (:obj:`dict`):
                The current data table.  *This is edited in place.*
            bkg_pairs (:obj:`str`):
                The method to use for including these pairs.

        Raises:
            PypeItError:
                Raised if the method selected to add the pairs is not
                defined.
        """
        numfiles = len(data['filename'])
        if bkg_pairs == 'empty':
            data['comb_id'] = [-1]*numfiles
            data['bkg_id'] = [-1]*numfiles
            return
        msgs.error('{0} not a defined method for the background pair columns.'.format(bkg_pairs))

    def _new_build(self, file_list, strict=True):
        mandatory_meta = define_mandatory_meta()
        meta_keys = list(mandatory_meta.keys())

        # Now additional spectrograph specific items
        all_meta = mandatory_meta.copy()
        other_meta = define_other_meta()
        config_keys = self.spectrograph.configuration_keys()
        for key in config_keys:
            if key not in meta_keys:
                meta_keys += [key]
                all_meta[key] = other_meta[key]

        # Build
        data = {k:[] for k in meta_keys}

        for ifile in file_list:
            # Read the fits headers
            headarr = self.spectrograph.get_headarr(ifile, strict=strict)
            # Add the directory, file name, and instrument to the table
            d,f = os.path.split(ifile)
            data['directory'].append(d)
            data['filename'].append(f)
            data['instrume'].append(self.spectrograph.spectrograph)
            # Grab Meta
            for meta_key in data.keys():
                # Skip external ones
                if meta_key in ['filename', 'directory', 'instrume']:
                    continue
                # Grab it
                value = self.spectrograph.get_meta(ifile, meta_key, headarr=headarr, required=True)
                if all_meta[meta_key]['dtype'] == str:
                    data[meta_key].append(str(value).strip())
                elif all_meta[meta_key]['dtype'] == int:
                    data[meta_key].append(int(value))
                elif all_meta[meta_key]['dtype'] == float:
                    try:
                        data[meta_key].append(float(value))
                    except TypeError:
                        debugger.set_trace()
                elif all_meta[meta_key]['dtype'] == tuple:
                    assert isinstance(value, tuple)
                    data[meta_key].append(value)
                else:
                    debugger.set_trace()

        # Additional bits and pieces
        self._add_bkg_pairs(data, 'empty')
        return data

    def _build(self, file_list, strict=True, bkg_pairs='empty'):
        """
        Returns a dictionary with the data to be included in the table.
        """
        # Get the header keywords specific to the provided spectrograph.
        # head_keys is a nested dictionary listing the header keywords
        # for each extension in the file.  The top-level dictionary key
        # is just the 0-indexed number of the extension
        head_keys = self.spectrograph.header_keys()

        # The table is declared based on this input dictionary: The
        # directory, filename, instrument are always included
        data = {k:[] for k in PypeItMetaData.default_keys()}

        # Add columns to the output table for each keyword.  The
        # keywords from all extensions must be unique.
        ext = {}
        for i in head_keys.keys():
            for k in head_keys[i].keys():
                if k in data.keys():
                    raise ValueError('Keywords are not unique across all extensions!')
                ext[k] = i
                data[k] = []

        # TODO: Stopgap for required keys in fitstbl used by other parts of
        # the code.  Need to decide how to handle these.
        required_columns = ['time', 'date', 'target']
        required_for_ABBA = ['ra', 'dec'] # and for standards?
        if any([ c not in data.keys() for c in required_columns]):
            msgs.warn('Columns are missing.')

        # Number of files to read
        numfiles = len(file_list)

        # Loop on files
        for ifile in file_list:

            # Read the fits headers
            headarr = self.spectrograph.get_headarr(ifile, strict=strict)

            # Check that the header is valid
            # TODO: The check_headers function needs to be implemented
            # for each instrument.  spectrograph.check_headers() should
            # raise an exception with an appropriate message.
            try:
                # TODO: Move this into spectrograph.validate_fitstbl()
                self.spectrograph.check_headers(headarr)
            except Exception as e:
                msgs.warn('Reading of headers from file:' + msgs.newline() + ifile
                          + msgs.newline() + 'failed with the following exception'
                          + msgs.newline() + e.__repr__() + msgs.newline() +
                          'Please check that the file was taken with the provided instrument:'
                          + msgs.newline() + '{0}'.format(self.spectrograph.spectrograph)
                          + msgs.newline() + 'Then either change the instrument or remove/skip '
                          + 'the file.' + msgs.newline()+ 'Continuing by ignoring this file...')
                numfiles -= 1
                continue

            # Add the directory, file name, and instrument to the table
            d,f = os.path.split(ifile)
            data['directory'].append(d)
            data['filename'].append(f)
            data['instrume'].append(self.spectrograph.spectrograph)

            # Now get the rest of the keywords
            for key in data.keys():
                if key in PypeItMetaData.default_keys():
                    continue

                # Try to read the header item
                try:
                    value = headarr[ext[key]][head_keys[ext[key]][key]]
                except KeyError:
                    # Keyword not found in header
                    msgs.warn("{0} keyword not in header. Setting to None".format(key))
                    value = 'None'
                except TypeError:
                    import pdb; pdb.set_trace()

                # Convert the time to days
                # TODO: Done here or as a method in Spectrograph?
                if key == 'time' and value != 'None':
                    # HERE IS A KULDGE
                    if 'date' in data.keys():
                        idate = headarr[ext['date']][head_keys[ext['date']]['date']]
                    else:
                        idate = None
                    value = self.convert_time(value, date=idate)

                # Set the value
                vtype = type(value)
                if np.issubdtype(vtype, np.str_):
                    value = value.strip()
                if np.issubdtype(vtype, np.integer) or np.issubdtype(vtype, np.floating) \
                        or np.issubdtype(vtype, np.str_) or np.issubdtype(vtype, np.bool_):
                    data[key].append(value)
                else:
                    msgs.bug('Unexpected type, {1}, for key {0}'.format(key,
                             vtype).replace('<type ','').replace('>',''))

            msgs.info('Successfully loaded headers for file:' + msgs.newline() + ifile)

        # Report
        msgs.info("Headers loaded for {0:d} files successfully".format(numfiles))
        if numfiles != len(file_list):
            msgs.warn("Headers were not loaded for {0:d} files".format(len(file_list) - numfiles))
        if numfiles == 0:
            msgs.error("The headers could not be read from the input data files."
                       + msgs.newline() + "Please check that the settings file matches the data.")

#        # Add the index column
#        data['index'] = np.arange(numfiles).tolist()

        # Add the background_index column
        #if bkg_pairs is not None:
        self._add_bkg_pairs(data, bkg_pairs)

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
                                                                    self.spectrograph.spectrograph),
                                        '              length={0}\n'.format(len(self))])

    def _repr_html_(self):
        return self.table._base_repr_(html=True, max_width=-1,
                            descr_vals=['PypeItMetaData: spectrograph={0}, length={1}\n'.format(
                                                    self.spectrograph.spectrograph, len(self))])

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
        # Check the input
        if not isinstance(usrdata, table.Table):
            raise TypeError('Must provide an astropy.io.table.Table instance.')
        if 'filename' not in usrdata.keys():
            raise KeyError('The user-provided table must have \'filename\' column!')

        # Make sure the data are correctly ordered
        srt = [np.where(f == self.table['filename'])[0][0] for f in usrdata['filename']]

        # Convert types if possible
        existing_keys = list(set(self.table.keys()) & set(usrdata.keys()))
        if len(existing_keys) > 0 and match_type:
            for key in existing_keys:
                if self.table[key].dtype.type != usrdata[key].dtype.type:
                    try:
                        usrdata[key] = usrdata[key].astype(self.table[key].dtype)
                    except:
                        pass
        
        # Include the user data in the table
        for key in usrdata.keys():
            self.table[key] = usrdata[key][srt]

#    @staticmethod
#    def get_utc(headarr):
#        """
#        Find and return the UTC for a file based on the headers read from
#        all extensions.
#    
#        The value returned is the first UT or UTC keyword found any in any
#        header object.
#    
#        Args:
#            headarr (list):
#                List of :obj:`astropy.io.fits.Header` objects to search
#                through for a UTC of the observation.
#        Returns:
#            object: The value of the header keyword.
#        """
#        for h in headarr:
#            if h == 'None':
#                continue
#            if 'UTC' in h.keys():
#                return h['UTC']
#            elif 'UT' in h.keys():
#                return h['UT']
#        return None

    def convert_time(self, in_time, date=None):
        """
        Convert the time read from a file header to MJD for all spectrographs.
    
        Args:
            in_time (str):
                The time read from the file header
            date (str, optional):
                The date read from the file header

        Returns:
            float: The MJD time
        """
        # Convert seconds to hours
        if self.spectrograph.timeunit == 's':
            msgs.error("SHOULD NOT GET HERE")
            return float(in_time)/3600.0
    
        # Convert minutes to hours
        if self.spectrograph.timeunit == 'm':
            msgs.error("SHOULD NOT GET HERE")
            return float(in_time)/60.0

        # Convert from an astropy.Time format
        if self.spectrograph.timeunit in time.Time.FORMATS.keys():
            if date is not None:
                if 'T' not in date:
                    in_time = date+'T'+in_time
                else:
                    in_time = date
            ival = float(in_time) if self.spectrograph.timeunit == 'mjd' else in_time
            tval = time.Time(ival, scale='tt', format=self.spectrograph.timeunit)
            # Put MJD
            return tval.mjd
        
        msgs.error('Bad time unit')

    def get_configuration(self, indx, cfg_keys=None):
        """
        Return the configuration dictionary for a given frame.

        This is not the same as the backwards compatible "setup"
        dictionary.

        Args:
            indx (:obj:`int`):
                The index of the table row to use to construct the
                configuration.
            cfg_keys (:obj:`list`, optional):
                The list of metadata keys to use to construct the
                configuration.  If None, the `configuration_keys` of
                :attr:`spectrograph` is used.

        Returns:
            dict: A dictionary with the metadata values from the
            selected row.
        """
        _cfg_keys = self.spectrograph.configuration_keys() if cfg_keys is None else cfg_keys
        return {k:self.table[k][indx] for k in _cfg_keys}

    def master_key(self, row, det=1):
        """
        Construct the master key for the file in the provided row.

        The master key is the combination of the configuration, the
        calibration group, and the detector.  The configuration ID is
        the same as included in the configuration column (A, B, C, etc),
        the calibration group is the same as the calibration bit number,
        and the detector number is provided as an argument and converted
        to a zero-filled string with two digits (the maximum number of
        detectors is 99).

        Using the calibration bit in the keyword allows MasterFrames to
        be used with multiple calibration groups.

        Args:
            row (:obj:`int`):
                The 0-indexed row used to construct the key.
            det (:obj:`int`, optional):
                The 1-indexed detector.  Default is 1.

        Returns:
            str: Master key with configuration, calibration group(s),
            and detector.

        Raises:
            PypeItError:
                Raised if the 'setup' or 'calibbit' columns
                haven't been defined.
        """
        if 'setup' not in self.keys() or 'calibbit' not in self.keys():
            msgs.error('Cannot provide master key string without setup and calibbit; '
                       'run set_configurations and set_calibration_groups.')
        return '{0}_{1}_{2}'.format(self['setup'][row], self['calibbit'][row], str(det).zfill(2))

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
        try:
            return time.Time(self['time'][row], format='mjd')
        except:
            msgs.warn('There is no time column in the metadata table!' + msgs.newline() +
                      'The time and heliocentric corrections will be wrong!' + msgs.newline() +
                      'This is a bad idea. Continuing with a dummy time value')
            return '2010-01-01'

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
        return '{0}_{1}_{2}{3}'.format(self['target'][row].replace(" ", ""),
                                       self.spectrograph.camera,
                                       datetime.datetime.strftime(dtime, '%Y%b%dT'),
                                       tiso.value.split("T")[1].replace(':',''))

    def get_setup(self, row, det=None, config_only=False):
        """
        Construct the setup dictionary.

        .. todo::
            - This is for backwards compatibility, but we should
              consider reformatting it.  And it may be something to put
              in the relevant spectrograph class.

        Args:
            row (:obj:`int`):
                The 0-indexed row used to construct the setup.
            det (:obj:`int`, optional):
                The 1-indexed detector to include.  If None, all
                detectors are included.
            config_only (:obj:`bool`, optional):
                Just return the dictionary with the configuration, don't
                include the top-level designation of the configuration
                itself.

        Returns:
            dict: The pypeit setup dictionary with the default format.

        Raises:
            PypeItError:
                Raised if the 'setup' isn't been defined.
        """
        if 'setup' not in self.keys():
            msgs.error('Cannot provide instrument setup without \'setup\' column; '
                       'run set_configurations.')
        dispname = 'none' if 'dispname' not in self.keys() else self['dispname'][row]
        dispangle = 'none' if 'dispangle' not in self.keys() else self['dispangle'][row]
        dichroic = 'none' if 'dichroic' not in self.keys() else self['dichroic'][row]
        decker = 'none' if 'decker' not in self.keys() else self['decker'][row]
        slitwid = 'none' if 'slitwid' not in self.keys() else self['slitwid'][row]
        slitlen = 'none' if 'slitlen' not in self.keys() else self['slitlen'][row]
        binning = '1,1' if 'binning' not in self.keys() else self['binning'][row]

        setup = {self['setup'][row]:
                    {'--':
                        {'disperser': {'name': dispname, 'angle':dispangle},
                         'dichroic': dichroic,
                         'slit': {'decker': decker, 'slitwid':slitwid, 'slitlen':slitlen}}}}
        _det = np.arange(self.spectrograph.ndet)+1 if det is None else [det]
        for d in _det:
            setup[self['setup'][row]][str(d).zfill(2)] \
                    = {'binning': binning, 'det': d,
                       'namp': self.spectrograph.detector[d-1]['numamplifiers']}
        return setup[self['setup'][row]] if config_only else setup

    def get_configuration_names(self, ignore=None, return_index=False):
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
            return_index (:obj:`bool, optional):
                Return row indices with the first occurence of these
                configurations.

        Returns:
            numpy.array: The list of unique setup names.  A second
            returned object provides the indices of the first occurrence
            of these setups, if requested.

        Raises:
            PypeItError:
                Raised if the 'setup' isn't been defined.
        """
        if 'setup' not in self.keys():
            msgs.error('Cannot get setup names; run set_configurations.')

        # Unique configurations
        setups, indx = np.unique(self['setup'], return_index=True)

        if ignore is not None:
            # Remove the selected configurations to ignore
            rm = np.invert(np.isin(setups, ignore))
            setups = setups[rm]
            indx = indx[rm]

        return setups, indx if return_index else setups

    def get_all_setups(self, ignore=None):
        """
        Construct the full dictionary with all the instrument setup
        detail.

        This is mostly a convenience function for the writing routines.

        .. todo::
            - This is for backwards compatibility, but we should
              consider reformatting/removing it.

        Args:
            ignore (:obj:`list`, optional):
                Ignore configurations in the provided list.
        """
        setups, indx = self.get_configuration_names(ignore=ignore, return_index=True)
        return {setup: self.get_setup(i, config_only=True) for setup,i in zip(setups,indx)}

    def unique_configurations(self, ignore_frames=None, force=False):
        """
        Return the unique instrument configurations.

        If run before the 'setup' column is initialized, this function
        determines the unique instrument configurations by finding
        unique combinations of the items in the metadata table listed by
        the spectrograph `configuration_keys` method.

        If run after the 'setup' column has been set, this simply
        constructs the configuration dictionary using the unique
        configurations in that column.

        This is used to set the internal :attr:`configs`.

        Args:
            ignore_frames (:obj:`list`, optional):
                A list of strings with the frame types (e.g., ['bias',
                'dark']) to ignore when determining the unique
                configurations.  This is only possible if the frame
                types have been previously defined and can be selected
                using :func:`find_frames`.  Default is to use all frame
                types.
            force (:obj:`bool`, optional):
                Force the configurations to be redetermined.  Otherwise
                the configurations are only determined if
                :attr:`configs` has not yet been define.

        Returns:
            dict: A nested dictionary, one dictionary per configuration
            with the associated metadata for each.

        Raises:
            PypeItError:
                Raised if `ignore_frames` is provided but the frame
                types have not been defined yet.
        """
        if self.configs is not None and not force:
            return self.configs

        if 'setup' in self.keys():
            msgs.info('Setup column already set.  Finding unique configurations.')
            uniq, indx = np.unique(self['setup'], return_index=True)
            ignore = uniq == 'None'
            if np.sum(ignore) > 0:
                msgs.warn('Ignoring {0} frames with configuration set to None.'.format(
                            np.sum(ignore)))
            self.configs = {}
            for i in range(len(uniq)):
                if ignore[i]:
                    continue
                self.configs[uniq[i]] = self.get_configuration(indx[i])
            msgs.info('Found {0} unique configurations.'.format(len(self.configs)))
            return self.configs

        msgs.info('Using metadata to determine unique configurations.')

        # If the frame types have been set, ignore anything listed in
        # the ignore_frames
        indx = np.arange(len(self))
        if ignore_frames is not None:
            if 'frametype' not in self.keys():
                msgs.error('To ignore frames, types must have been defined; run get_frame_types.')
            msgs.info('Unique configurations ignore frames with type: {0}'.format(ignore_frames))
            use = np.ones(len(self), dtype=bool)
            for ftype in ignore_frames:
                use &= np.invert(self.find_frames(ftype))
            indx = indx[use]
        if len(indx) == 0:
            msgs.error('No frames to use to define configurations!')

        # Get the list of keys to use
        cfg_keys = self.spectrograph.configuration_keys()

        # Configuration identifiers are iterations through the
        # upper-case letters: A, B, C, etc.
        cfg_iter = string.ascii_uppercase
        cfg_indx = 0

        # TODO: Placeholder: Allow an empty set of configuration keys
        # meaning that the instrument setup has only one configuration.
        if len(cfg_keys) == 0:
            self.configs = {}
            self.configs[cfg_iter[cfg_indx]] = {}
            msgs.info('All files assumed to be from a single configuration.')
            return self.configs

        # Use the first file to set the first unique configuration
        self.configs = {}
        self.configs[cfg_iter[cfg_indx]] = self.get_configuration(indx[0], cfg_keys=cfg_keys)
        cfg_indx += 1

        # Check if any of the other files show a different
        # configuration.  The check is for *exact* equality, meaning
        # *any* difference in the values for the keywords listed in
        # `cfg_keys` will lead to a new configuration.
        # TODO: Add a tolerance for floating point values?
        for i in indx[1:]:
            j = 0
            for c in self.configs.values():
                if np.all([c[k] == self.table[k][i] for k in cfg_keys]):
                    break
                j += 1
            unique = j == len(self.configs)
            if unique:
                if cfg_indx == len(cfg_iter):
                    msgs.error('Cannot assign more than {0} configurations!'.format(len(cfg_iter)))
                self.configs[cfg_iter[cfg_indx]] = self.get_configuration(i, cfg_keys=cfg_keys)
                cfg_indx += 1

        msgs.info('Found {0} unique configurations.'.format(len(self.configs)))
        return self.configs

    def set_configurations(self, configs=None, force=False):
        """
        Assign each frame to a configuration (setup) and include it in
        the metadata table.

        The internal table is edited *in place*.  If the 'setup' column
        already exists, the configurations are **not** reset unless you
        call the function with `force=True`.

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

        Raises:
            PypeItError:
                Raised if none of the keywords in the provided
                configuration match with the metadata keywords.
        """
        # Configurations have already been set
        if 'setup' in self.keys() and not force:
            return

        _configs = self.unique_configurations() if configs is None else configs
        for k, cfg in _configs.items():
            if len(set(cfg.keys()) - set(self.keys())) > 0:
                msgs.error('Configuration {0} defined using unavailable keywords!'.format(k))

        self.table['setup'] = 'None'
        nrows = len(self)
        for i in range(nrows):
            for d, cfg in _configs.items():
                if np.all([cfg[d] == self.table[d][i] for d in cfg.keys()]):
                    self.table['setup'][i] = d

    def _set_calib_group_bits(self):
        grp = np.empty(len(self), dtype=object)
        ngroups = 0
        for i in range(len(self)):
            if self['calib'][i] == 'None':
                continue
            grp[i] = eval('['+self['calib'][i]+']') if self['calib'][i].find(',') > -1 \
                        else [int(self['calib'][i])]
            ngroups = np.amax(np.append(np.array(grp[i])+1,ngroups))
        
        self.calib_bitmask = BitMask(np.arange(ngroups))
        self['calibbit'] = 0
        for i in range(len(self)):
            if self['calib'][i] == 'None':
                continue
            self['calibbit'][i] = self.calib_bitmask.turn_on(self['calibbit'][i], grp[i])

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
        return None if self.calib_bitmask is None else self.calib_bitmask.nbits
                
    def set_calibration_groups(self, global_frames=None, force=False):
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
            force (:obj:`bool`, optional):
                Force the calibration groups to be reconstructed if
                the 'calib' column already exists.

        Raises:
            PypeItError:
                Raised if 'setup' column is not defined, or if
                `global_frames` is provided but the frame types have not
                been defined yet.
        """
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
        configs = np.unique(self['setup'].data).tolist()
        try:
            configs.remove('None')      # Ignore frames with undefined configurations
        except:
            pass
        n_cfg = len(configs)

        # TODO: Science frames can only have one calibration group

        # Assign everything from the same configuration to the same
        # calibration group
        self.table['calib'] = 'None'
        for i in range(n_cfg):
            self['calib'][(self['setup'] == configs[i]) & (self['framebit'] > 0)] = str(i)
        
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
                :class:`pypeit.core.framematch.FrameTypeBitMask`.
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
        """
        if 'framebit' not in self.keys():
            raise ValueError('Frame types are not set.  First run get_frame_types.')
        if ftype is 'None':
            return self['framebit'] == 0
        indx = self.type_bitmask.flagged(self['framebit'], ftype)
        # Calib ID?
        if calib_ID is not None:
            in_grp = self.find_calib_group(calib_ID)
            #
            indx &= in_grp
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
            indx: ndarray of rows in fitstbl
        """
        indx = self.find_frames(ftype, calib_ID=calib_ID, index=True)
        return [os.path.join(d,f) for d,f in zip(self['directory'][indx], self['filename'][indx])], indx

    def frame_paths(self, indx):
        """
        Return the full paths to one or more frames.

        Args:
            indx (:obj:`int`, array-like):
                One or more 0-indexed rows in the table with the frames
                to return.  Can be an array of indices or a boolean
                array of the correct length.
        Returns:
            str, list: The full paths of one or more frames.
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
                types are overwitten by the provided types.
        """
        if not append:
            self['framebit'][indx] = 0
        self['framebit'][indx] = self.type_bitmask.turn_on(self['framebit'][indx], flag=frame_type)
        self['frametype'][indx] = self.type_bitmask.type_names(self['framebit'][indx])

    def get_frame_types(self, flag_unknown=False, user=None, useIDname=False, merge=True):
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
                therefore *must* match the number of files in the
                provided `fitstbl`.  For frames that have multiple
                types, the types should be provided as a string with
                comma-separated types.
            useIDname (:obj:`bool`, optional):
                Use ID name in the Header to image type
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

        # TODO: This needs to be moved into each Spectrograph
        if useIDname and 'idname' not in self.keys():
            raise ValueError('idname is not set in table; cannot use it for file typing.')

        # Start
        msgs.info("Typing files")
        type_bits = np.zeros(len(self), dtype=self.type_bitmask.minimum_dtype())
    
        # Use the user-defined frame types from the input dictionary
        if user is not None:
            if len(user.keys()) != len(self):
                raise ValueError('The user-provided dictionary does not match table length.')
            msgs.info('Using user-provided frame types.')
            for ifile,ftypes in user.items():
                indx = self['filename'] == ifile
                type_bits[indx] = self.type_bitmask.turn_on(type_bits[indx], flag=ftypes.split(','))
            return self.set_frame_types(type_bits, merge=merge)
    
        # Loop over the frame types
        for i, ftype in enumerate(self.type_bitmask.keys()):
    
            # Initialize: Flag frames with the correct ID name or start by
            # flagging all as true
            indx = self['idname'] == self.spectrograph.idname(ftype) if useIDname \
                        else np.ones(len(self), dtype=bool)
    
            # Include a combination of instrument-specific checks using
            # combinations of the full set of metadata
            exprng = self.par['scienceframe']['exprng'] if ftype == 'science' \
                        else self.par['calibrations']['{0}frame'.format(ftype)]['exprng']
            # TODO: Use & or | ?  Using idname above gets overwritten by
            # this if the frames to meet the other checks in this call.
            indx &= self.spectrograph.check_frame_type(ftype, self.table, exprng=exprng)
    
            # Turn on the relevant bits
            type_bits[indx] = self.type_bitmask.turn_on(type_bits[indx], flag=ftype)
    
        # Find the nearest standard star to each science frame
        # TODO: Should this be 'standard' or 'science' or both?
        if 'ra' not in self.keys() or 'dec' not in self.keys():
            msgs.warn('Cannot associate standard with science frames without sky coordinates.')
        else:
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
                foundstd = flux.find_standard_file(ra, dec, check=True)
                b = self.type_bitmask.turn_off(b, flag='science' if foundstd else 'standard')
    
        # Find the files without any types
        indx = np.invert(self.type_bitmask.flagged(type_bits))
        if np.any(indx):
            msgs.info("Couldn't identify the following files:")
            for f in self['filename'][indx]:
                msgs.info(f)
            if not flag_unknown:
                msgs.error("Check these files before continuing")
    
        # Now identify the dark frames
        # TODO: This should not be here.  Move to instrument specific
        # selection as above
        indx = self.type_bitmask.flagged(type_bits, flag='bias') \
                        & (self['exptime'].data.astype(float) > self.spectrograph.minexp)
        type_bits[indx] = self.type_bitmask.turn_on(type_bits[indx], 'dark')
    
        # Finish up (note that this is called above if user is not None!)
        msgs.info("Typing completed!")
        return self.set_frame_types(type_bits, merge=merge)

    def write_setups(self, ofile, overwrite=True, ignore=None):
        """
        Write the *.setups file.

        The *.setups file lists all the unique instrument configurations
        (setups).

        .. todo::
            - This is for backwards compatibility, but we should
              consider reformatting/removing it.

        Args:
            ofile (:obj:`str`):
                Name for the output sorted file.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file with the same name.
            ignore (:obj:`list`, optional):
                Ignore configurations in the provided list.
        """
        if os.path.isfile(ofile) and not overwrite:
            msgs.error('{0} already exists.  Use ovewrite=True to overwrite.'.format(ofile))
        
        # Construct file
        cfg = self.get_all_setups(ignore=ignore)
        ff = open(ofile, 'w')
        ff.write(yaml.dump(utils.yamlify(cfg)))
        ff.close()

    def write_sorted(self, ofile, overwrite=True, ignore=None):
        """
        Write the *.sorted file.

        The *.sorted file lists all the unique instrument configurations
        (setups) and the frames associated with each configuration.  The
        output data table is identical to the pypeit file output.

        .. todo::
            - This is for backwards compatibility, but we should
              consider reformatting/removing it.

        Args:
            ofile (:obj:`str`):
                Name for the output sorted file.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file with the same name.
            ignore (:obj:`list`, optional):
                Ignore configurations in the provided list.

        Raises:
            PypeItError:
                Raised if the 'setup' isn't been defined.
        """
        if 'setup' not in self.keys():
            msgs.error('Cannot write sorted instrument configuration table without \'setup\' '
                       'column; run set_configurations.')

        if os.path.isfile(ofile) and not overwrite:
            msgs.error('{0} already exists.  Use ovewrite=True to overwrite.'.format(ofile))

        # Columns for output
        output_cols = np.array(self.spectrograph.metadata_keys())
        output_cols = output_cols[np.isin(output_cols, self.keys())].tolist()

        # Unique configurations
        setups, indx = self.get_configuration_names(ignore=ignore, return_index=True)
        
        # Construct file
        ff = open(ofile, 'w')
        for setup,i in zip(setups,indx):
            # Get the configuration dictionary
            cfg = self.get_setup(i, config_only=True)
            # Get the subtable of frames taken in this configuration
            subtbl = self.table[output_cols][self['setup'] == setup]
#            subtbl.sort(['frametype','filename'])
            # Write the file
            ff.write('##########################################################\n')
            ff.write('Setup {:s}\n'.format(setup))
            ff.write(yaml.dump(utils.yamlify(cfg)))
            ff.write('#---------------------------------------------------------\n')
            subtbl.write(ff, format='ascii.fixed_width')
        ff.write('##end\n')
        ff.close()

    def write_calib(self, ofile, overwrite=True, ignore=None):
        """
        Write the *.calib file.

        The *.calib file provides the unique instrument configurations
        (setups) and the association of each frame from that
        configuration with a given calibration group.

        .. todo::
            - This is for backwards compatibility, but we should
              consider reformatting/removing it.
            - This is complicated by allowing some frame types to have
              no association with an instrument configuration

        Args:
            ofile (:obj:`str`):
                Name for the output sorted file.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file with the same name.
            ignore (:obj:`list`, optional):
                Ignore calibration groups in the provided list.

        Raises:
            PypeItError:
                Raised if the 'setup' or 'calibbit' columns haven't been
                defined.
        """
        if 'setup' not in self.keys() or 'calibbit' not in self.keys():
            msgs.error('Cannot write calibration groups without \'setup\' and \'calibbit\' '
                       'columns; run set_configurations and set_calibration_groups.')

        if os.path.isfile(ofile) and not overwrite:
            msgs.error('{0} already exists.  Use ovewrite=True to overwrite.'.format(ofile))

        # Construct the setups dictionary
        cfg = self.get_all_setups(ignore=['None'])

        # Iterate through the calibration groups and add frames to the
        # appropriate configuration
        for i in range(self.n_calib_groups):
            # Skip this group
            if ignore is not None and i in ignore:
                continue

            # Find the frames in this group
            in_group = self.find_calib_group(i)

            # Find the unique configurations in this group, ignoring any
            # undefined ('None') configurations
            setup = np.unique(self['setup'][in_group]).tolist()
            if 'None' in setup:
                setup.remove('None')

            # Make sure that each calibration group should only contain
            # frames from a single configuration
            if len(setup) != 1:
                msgs.error('Each calibration group must be from one and only one instrument '
                           'configuration with a valid letter identifier; i.e., the '
                           'configuration cannot be None.')

            # Find the frames of each type in this group
            cfg[setup[0]][i] = {}
            for key in self.type_bitmask.keys():
                ftype_in_group = self.find_frames(key) & in_group
                cfg[setup[0]][i][key] = [ os.path.join(d,f) 
                                                for d,f in zip(self['directory'][ftype_in_group],
                                                               self['filename'][ftype_in_group])]
        # Write it
        ff = open(ofile, 'w')
        ff.write(yaml.dump(utils.yamlify(cfg)))
        ff.close()

    def write_pypeit(self, ofile, split=False, overwrite=False, ignore=None, cfg_lines=None,
                     write_bkg_pairs=False):
        """
        Write a *.pypeit file in data-table format.

        The *.pypeit file is the main configuration file for PypeIt,
        configuring the control-flow and algorithmic parameters and
        listing the data files to read.  This function writes the
        columns selected by the
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.metadata_keys`,
        which can be specific to each instrument.

        Users can write the file, edit it, and then re-read it into
        PypeIt to change the designated frame type, instrument setup,
        calibration group, or source/background groupings.

        A pypeit file can contain multiple configurations or be split
        into one pypeit file per instrument configuration.

        Args:
            ofile (:obj:`str`):
                Name for the output pypeit file.
            split (:obj:`bool`, optional):

                Split the data by instrument configuration leading to
                one PypeIt file per setup.  The name of the output
                directory is the root of the provided output file name
                and the setup character.  The pypeit file in each
                directory is based on the same, but with the '.pypeit'
                extension.  By default, all instrument setups are kept
                in the same PypeIt file.

            overwrite (:obj:`bool`, optional):
                Overwrite any existing file(s).
            ignore (:obj:`list`, optional):
                Ignore configurations in the provided list.
            cfg_lines (:obj:`list`, optional):
                The list of configuration lines to include in the file.
                If None are provided, the vanilla configuration is
                included.
            write_bkg_pairs (:obj:`bool`, optional):
                When constructing the
                :class:`pypeit.metadata.PypeItMetaData` object, include
                two columns called `comb_id` and `bkg_id` that identify
                object and background frame pairs.  The string indicates
                how these these columns should be added::
                    - `empty`: The columns are added but their values
                      are all originally set to -1.  **This is
                      currently the only option.**

        Raises:
            PypeItError:
                Raised if the 'setup' isn't defined and split is True.
        """
        if 'setup' not in self.keys() and split:
            msgs.error('Cannot write pypeit file split by setup; run set_configurations.')

        # TODO: I don't think sortroot is ever used.

        # Columns for output
        columns = self.spectrograph.metadata_keys()

        # comb, bkg columns
        if write_bkg_pairs:
            for key in ['comb_id', 'bkg_id']:
                if key not in columns:
                    columns += [key]

        # Take only those present
        output_cols = np.array(columns)
        output_cols = output_cols[np.isin(output_cols, self.keys())].tolist()

        if not split:
            # Write all configurations to a single file

            # Get the setup lines
            cfg = self.get_all_setups(ignore=ignore)
            setup_lines = yaml.dump(utils.yamlify(cfg)).split('\n')[:-1]
            # Get the paths
            paths = np.unique(self['directory']).tolist()
            # Get the data lines
            subtbl = self.table[output_cols]
            subtbl.sort(['frametype','filename'])
            with io.StringIO() as ff:
                subtbl.write(ff, format='ascii.fixed_width')
                data_lines = ff.getvalue().split('\n')[:-1]
            # Write the file
            make_pypeit_file(ofile, self.spectrograph.spectrograph, [], cfg_lines=cfg_lines,
                             setup_lines=setup_lines, sorted_files=data_lines, paths=paths)
            return

        # Write each configuration to its own file

        # Unique configurations
        setups, indx = self.get_configuration_names(ignore=ignore, return_index=True)

        for setup,i in zip(setups, indx):
            # Create the output directory
            root = '{0}_{1}'.format(self.spectrograph.spectrograph, setup)
            odir = os.path.join(os.path.split(ofile)[0], root)
            if not os.path.isdir(odir):
                os.makedirs(odir)
            # Create the output file name
            _ofile = os.path.join(odir, '{0}.pypeit'.format(root))
            # Get the setup lines
            cfg = self.get_setup(i, config_only=True)
            setup_lines = ['Setup {0}'.format(setup)]
            setup_lines += yaml.dump(utils.yamlify(cfg)).split('\n')[:-1]
            # Get the paths
            in_cfg = self['setup'] == setup
            paths = np.unique(self['directory'][in_cfg]).tolist()
            # Get the data lines
            subtbl = self.table[output_cols][self['setup'] == setup]
            subtbl.sort(['frametype','filename'])
            with io.StringIO() as ff:
                subtbl.write(ff, format='ascii.fixed_width')
                data_lines = ff.getvalue().split('\n')[:-1]
            # Write the file
            make_pypeit_file(_ofile, self.spectrograph.spectrograph, [], cfg_lines=cfg_lines,
                             setup_lines=setup_lines, sorted_files=data_lines, paths=paths)


    def write(self, ofile, columns=None, format=None, overwrite=False):
        """
        Write the metadata for the files to reduce.
    
        The table is written with the filename and frametype columns
        first.  All remaining columns, or a subset of them selected by
        `columns`, follow these first two.
   
        For a pypit file, use 'ascii.fixed_width' to print the table.
        
        Args:
            ofile (:obj:`str`, file-like):
                Output file name or file stream.  Passed directly to the
                `astropy.table.Table.write` function.
            columns (:obj:`list`, optional):
                A list of columns to include in the output file.  If None,
                all columns in `fitstbl` are written.
            format (:obj:`str`, optional):
                Format for the file output.  See
                :func:`astropy.table.Table.write`.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file; otherwise raise an
                exception.

        Raises:
            ValueError:
                Raised if the columns to include are not unique.
            FileExistsError:
                Raised if overwrite is False and the file exists.
        """
        if os.path.isfile(ofile) and not overwrite:
            raise FileExistsError('File {0} already exists.'.format(ofile) 
                                    + '  Change file name or set overwrite=True.')

        msgs.info('Writing fits file metadata to {0}.'.format(ofile))
    
        # Set the columns to include and check that they are unique
        _columns = list(self.keys()) if columns is None else columns
        if len(np.unique(_columns)) != len(_columns):
            # TODO: A warning may suffice...
            raise ValueError('Column names must be unique!')

        # Force the filename and frametype columns to go first
        col_order = [ 'filename', 'frametype' ]
        col_order += list(set(_columns) - set(col_order))

        # Remove any columns that don't exist
        popme = []
        for kk,c in enumerate(col_order):
            if c not in self.keys():
                msgs.warn('{0} is not a valid column!  Removing from output.'.format(c))
                popme.append(kk)
        popme.reverse()
        for index in popme:
            col_order.pop(index)

        # Write the output
        self.table[col_order].write(ofile, format=format, overwrite=overwrite)

    def find_calib_group(self, grp):
        """
        Find all the frames associated with the provided calibration group.
        """
        return self.calib_bitmask.flagged(self['calibbit'].data, grp)

    def find_frame_calib_groups(self, row):
        """
        Find the calibration groups associated with a specific frame.
        """
        return self.calib_bitmask.flagged_bits(self['calibbit'][row])

    def calib_to_science(self):
        """
        Construct the science ID based on the calibration group.
        """
        self['sci_ID'] = 0
        self['failures'] = False                    # TODO: placeholder

        is_science = self.find_frames('science')
        is_not_science = np.invert(is_science)
        sciid = np.arange(np.sum(is_science))       # Science image IDs
        self['sci_ID'][is_science] = 1 << sciid
        j = 0
        for i in range(len(self)):
            if is_not_science[i]:
                continue
            
            calibs = self.find_calib_group(int(self['calib'][i])) & is_not_science
            self['sci_ID'][calibs] |= (1 << sciid[j])
            j += 1
    
    def match_to_science(self, calib_par, calwin, fluxcalib_par, setup=False, verbose=True):
        """

        For a given set of identified data, match calibration frames to
        science frames.
        
        The matching is performed using criteria provided by
        `self.spectrograph.get_match_criteria()`; see `these functions`.

        For calwin to work, must have time column.

        Args:
            calib_par (:class:`pypeit.par.pypitpar.CalibrationsPar`):
                The calibration parameter set, which houses all the
                :class:`pypeit.par.pypeitpar.FrameGroupPar` instances
                that list how many of each calibration frame are needed.
            calwin (:obj:`float`):
                The minimum time difference between a calibration and
                science exposure allowed.
            fluxcalib_par (:class:`pypeit.par.pypitpar.FluxCalibPar`):
                Needed for worrying about standard stars
            setup (:obj:`bool`, optional):
                The operation is being executed only for an initial
                setup.  **Need a good description of the behavior this
                causes.**
            verbose (:obj:`bool`, optional):
                Provide verbose output.

        """
        if verbose:
            msgs.info("Matching calibrations to Science frames")
            if 'time' not in self.keys():
                end = ', and cannot match calibrations based on time difference.' \
                        if calwin > 0.0 else '.'
                msgs.warn('Table does not include time column!  Cannot sort calibration frames '
                          'by observation time with respect to each science frame{0}'.format(end))
        
        match_dict = self.spectrograph.get_match_criteria()

        # New columns
        self['failures'] = False
        self['sci_ID'] = 0

        # Loop on science frames
        indx = np.arange(len(self))[self.find_frames('science')]
        for ss, sci_idx in enumerate(np.where(self.find_frames('science'))[0]):
            target = self['target'][sci_idx] if 'target' in self.keys() \
                            else self['filename'][sci_idx]
            if verbose:
                msgs.info('='*50)
                msgs.info('Matching calibrations to frame: {0}'.format(target))

            # Science ID (trivial but key to the bit-wise that follows)
            self['sci_ID'][sci_idx] = 2**ss

            # Find matching (and nearby) calibration frames
            for ftag in self.type_bitmask.keys():
                # Skip science frames
                if ftag == 'science':
                    continue

                # bias/dark check to make sure we need to find matching frames
                if ftag == 'dark' and calib_par['biasframe']['useframe'] != 'dark':
                    if verbose:
                        msgs.info('  Dark frames not required.  Not matching...')
                    continue
                if ftag == 'bias' and calib_par['biasframe']['useframe'] != 'bias' \
                            and not calib_par['badpix']:
                    if verbose:
                        msgs.info('  Bias frames not required.  Not matching...')
                    continue

                # Find the files of this type
                gd_match = self.find_frames(ftag)


                # How many matching frames are required?  This is instrument specific
                numfr = (1 if ftag == 'arc' else 0) if setup \
                            else calib_par['{0}frame'.format(ftag)]['number']

                # If not required and not doing setup, continue
                if numfr == 0 and not setup:
                    if verbose:
                        msgs.info('   No {0} frames are required.  Not matching...'.format(ftag))
                    continue

                # Now go ahead and match the frames
                if 'match' not in match_dict[ftag].keys() and (not setup):
                    msgs.error('Match criteria needed for {0}!'.format(ftag))
                elif 'match' not in match_dict[ftag].keys() and verbose:
                    msgs.info('No matching criteria for {0} frames'.format(ftag))
                else:
                    chkk = match_dict[ftag]['match'].keys()
                    for ch in chkk:
                        tmtch = match_dict[ftag]['match'][ch]
                        gd_match &= framematch.match_logic(ch, tmtch, self.table, sci_idx)

                # Find the time difference between the calibrations and science frames
                if calwin > 0.0 and 'time' in self.keys():
                    tdiff = np.abs(self['time']-self['time'][sci_idx])
                    gd_match &= tdiff <= calwin

                # Now find which of the remaining n are the appropriate calibration frames
                nmatch = np.sum(gd_match)
                if verbose:
                    msgs.info('  Found {0} {1} frame for {2} ({3} required)'.format(nmatch, ftag,
                                                                                    target, numfr))

                # Have we identified enough of these calibration frames to continue?
                if nmatch < np.abs(numfr):
                    code = framematch.match_warnings(calib_par, ftag, nmatch, numfr, target,
                                                     fluxpar=fluxcalib_par)
                    if code == 'break':
                        self['failures'][sci_idx] = True
                        self['sci_ID'][sci_idx] = -1  # This might break things but who knows..
                        break
                else:
                    # Select the closest calibration frames to the science frame
                    wa = np.where(gd_match)[0]
                    if 'time' in self.keys():
                        tdiff = np.abs(self['time'][wa]-self['time'][sci_idx])
                        wa = wa[np.argsort(tdiff)]
                    # Identify the relevant calibration frames with this
                    # science frame
                    if setup or (numfr < 0):
                        self['sci_ID'][wa] |= 2**ss
                    else:
                        self['sci_ID'][wa[:numfr]] |= 2**ss
        if verbose:
            msgs.info('Science frames successfully matched to calibration frames.')

    def match_ABBA(self, max_targ_sep=30, max_nod_sep=2, merge=True):
        """
        Find ABBA exposure series based on the target and nodding
        separation.

        The :attr:`table` must have `ra` and `dec` columns to execute
        this function.

        Args:
            max_targ_sep (:obj:`int`, optional):
                The maximum separation (arcsec) between frame sky
                coordinates allowed for identifying frames taken of the
                same object.  Note that the default (30 arcsec, 1e-2
                deg) is arbitrary.
            max_nod_sep (:obj:`int`, optional):
                The maximum separation (arcsec) between the 1st and 4th
                frame sky coordinates in the ABBA sequence that is
                allowed when identifying the sequence.  Note that the
                default (2 arcsec) is arbitrary.
            merge (:obj:`bool`, optional):
                Merge the frame typing into the exiting table in a
                `AB_frame` column.  Otherwise, a single column table is
                returned with this column.

        Returns:
            :class:`astropy.table.Table`: A single column Table is
            returned if `merge` is False.  Otherwise, nothing is
            returned.

        Raises:
            KeyError:
                Raised if the internal table does not have `ra` and
                `dec` columns.
        """

        if 'ra' not in self.keys() or 'dec' not in self.keys():
            raise KeyError('Table must have ra and dec columns to find ABBA sequences.')

        # Make coords for all observations, including calibration frames
        coords = coordinates.SkyCoord(ra=self['ra'], dec=self['dec'], unit='deg')

        # Indices of science frames
        sci = self.find_frames('science')

        # Create empty dictionary of targets
        targets = {}

        # The identifying key for each target; allows for target to be
        # undefined
        tkey = 'target' if 'target' in self.keys() else 'filename'

        # Use a boolean array to identify each frame as matched to a
        # specific target.  All calibration frames are ignored, and the
        # science frames are all initialized as not being matched.
        free = np.zeros(len(self), dtype=bool)
        free[sci] = True

        while np.any(free):
            # Find the first science frame that hasn't been matched
            indx = np.where(free)[0][0]
            # Find other science frames that match the coordinates of
            # this science frame within the specified tolerance
            match = sci & (coords[indx].separation(coords).arcsec < max_targ_sep)
            # Add this target or filename to the dictionary and identify
            # the list of table rows matched to it.
            targets[self[tkey][indx]] = np.where(match)[0]
            # Set the relevant frames as no longer free to be matched
            free[match] = False

        # Group the frames
        t = table.Table({'AB_frame': framematch.group_AB_frames(self['filename'], targets, coords,
                                                                max_nod_sep=max_nod_sep)})

        if merge:
            # Merge the column into the existing table
            self['AB_frame'] = t['AB_frame']
            return
        # Return a new Table with the single column
        return t

def dummy_fitstbl(nfile=10, spectrograph='shane_kast_blue', directory='', notype=False):
    """
    Generate a dummy fitstbl for testing

    Parameters
    ----------
    nfile : int, optional
      Number of files to mimic
    spectrograph : str, optional
      Name of spectrograph to mimic
    notype : bool (optional)
      If True, do not add image type info to the fitstbl

    Returns
    -------
    fitstbl : PypeItMetaData

    """
    fitsdict = dict({'directory': [], 'filename': [], 'utc': []})
    fitsdict['index'] = np.arange(nfile)
    fitsdict['utc'] = ['2015-01-23']*nfile
    fitsdict['directory'] = [directory]*nfile
    fitsdict['filename'] = ['b{:03d}.fits.gz'.format(i) for i in range(nfile)]
    # TODO: The below will fail at 60
    fitsdict['date'] = ['2015-01-23T00:{:02d}:11.04'.format(i) for i in range(nfile)]
    fitsdict['time'] = [(1432085758+i*60)/3600. for i in range(nfile)]
    fitsdict['target'] = ['Dummy']*nfile
    fitsdict['ra'] = ['00:00:00']*nfile
    fitsdict['dec'] = ['+00:00:00']*nfile
    fitsdict['exptime'] = [300.] * nfile
    fitsdict['naxis0'] = [2048] * nfile
    fitsdict['naxis1'] = [2048] * nfile
    fitsdict['dispname'] = ['600/4310'] * nfile
    fitsdict['dichroic'] = ['560'] * nfile
    fitsdict['dispangle'] = ['none'] * nfile
    fitsdict["binning"] = ['1x1']*nfile
    fitsdict["airmass"] = [1.0]*nfile

    if spectrograph == 'shane_kast_blue':
        fitsdict['numamplifiers'] = [1] * nfile
        fitsdict['naxis0'] = [2112] * nfile
        fitsdict['naxis1'] = [2048] * nfile
        fitsdict['slitwid'] = [1.] * nfile
        fitsdict['slitlen'] = ['none'] * nfile
        # Lamps
        for i in range(1,17):
            fitsdict['lampstat{:02d}'.format(i)] = ['off'] * nfile
        fitsdict['exptime'][0] = 0        # Bias
        fitsdict['lampstat06'][1] = 'on'  # Arc
        fitsdict['exptime'][1] = 30       # Arc
        fitsdict['lampstat01'][2] = 'on'  # Trace, pixel, slit flat
        fitsdict['lampstat01'][3] = 'on'  # Trace, pixel, slit flat
        fitsdict['exptime'][2] = 30     # flat
        fitsdict['exptime'][3] = 30     # flat
        fitsdict['ra'][4] = '05:06:36.6'  # Standard
        fitsdict['dec'][4] = '52:52:01.0'
        fitsdict['airmass'][4] = 1.2
        fitsdict['ra'][5] = '07:06:23.45' # Random object
        fitsdict['dec'][5] = '+30:20:50.5'
        fitsdict['decker'] = ['0.5 arcsec'] * nfile

    # arrays
    for k in fitsdict.keys():
        fitsdict[k] = np.array(fitsdict[k])

    fitstbl = PypeItMetaData(load_spectrograph(spectrograph), data=fitsdict)
    fitstbl['instrume'] = spectrograph
    type_bits = np.zeros(len(fitstbl), dtype=fitstbl.type_bitmask.minimum_dtype())

    # Image typing
    if not notype:
        if spectrograph == 'shane_kast_blue':
            #fitstbl['sci_ID'] = 1  # This links all the files to the science object
            type_bits[0] = fitstbl.type_bitmask.turn_on(type_bits[0], flag='bias')
            type_bits[1] = fitstbl.type_bitmask.turn_on(type_bits[1], flag='arc')
            type_bits[2:4] = fitstbl.type_bitmask.turn_on(type_bits[2:4], flag=['pixelflat', 'trace'])
            type_bits[4] = fitstbl.type_bitmask.turn_on(type_bits[4], flag='standard')
            type_bits[5:] = fitstbl.type_bitmask.turn_on(type_bits[5:], flag='science')
            fitstbl.set_frame_types(type_bits)
            # Calibration groups
            cfgs = fitstbl.unique_configurations(ignore_frames=['bias', 'dark'])
            fitstbl.set_configurations(cfgs)
            fitstbl.set_calibration_groups(global_frames=['bias', 'dark'])

    return fitstbl


def define_mandatory_meta():
    mandatory_meta = {}
    # Filename
    mandatory_meta['directory'] = dict(dtype=str, comment='Path to raw data file')
    mandatory_meta['filename'] = dict(dtype=str, comment='Basename of raw data file')

    # Instrument related
    mandatory_meta['instrume'] = dict(dtype=str, comment='Spectrograph name')
    mandatory_meta['dispname'] = dict(dtype=str, comment='Disperser name')
    mandatory_meta['decker'] = dict(dtype=str, comment='Slit/mask/decker name')
    mandatory_meta['binning'] = dict(dtype=tuple, comment='(spatial,spectral) binning')

    # Target
    mandatory_meta['target'] = dict(dtype=str, comment='Name of the target')
    mandatory_meta['ra'] = dict(dtype=str, comment='Colon separated (J2000) RA')
    mandatory_meta['dec'] = dict(dtype=str, comment='Colon separated (J2000) DEC')

    # Obs
    mandatory_meta['time'] = dict(dtype=str, comment='Date+time readable by astropy.time.Time')
    mandatory_meta['airmass'] = dict(dtype=float, comment='Airmass')
    mandatory_meta['exptime'] = dict(dtype=float, comment='Exposure time')

    # Return
    return mandatory_meta

def define_other_meta():
    other_meta = {}

    other_meta['dichroic'] = dict(dtype=str, comment='Beam splitter')

    return other_meta
