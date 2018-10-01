"""
Provides a class that handles the fits metadata required by PypeIt.
"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from __future__ import unicode_literals

import os

import numpy as np

from astropy import table, coordinates, time

from pypeit import msgs
from pypeit.core import framematch
from pypeit.core import flux
from pypeit.par import PypeItPar
from pypeit.spectrographs.util import load_spectrograph

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
        This should get abstracted to be independent of the
        spectrograph, with all the spectrograph dependent keywords be
        passed into the function.  The validation of the table would
        then happen in PypeItSetup.

    Args:
        spectrograph (:obj:`str`,
            :class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph used to collect the data save to each file.
            The class is used to provide the header keyword data to
            include in the table and specify any validation checks.

        file_list (:obj:`list`, optional):
            The list of files to include in the table.
        data (table-like, optional):
            The data to incude in the table.  The type can be anything
            allowed by the instantiation of
            :class:`astropy.table.Table`.
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
        par
        bitmask (:class:`pypeit.core.framematch.FrameTypeBitMask`):
            The bitmask used to set the frame type of each fits file.
        table (:class:`astropy.table.Table`):
            The table with the relevant metadata for each fits file to
            use in the data reduction.
    """
    def __init__(self, spectrograph, par=None, file_list=None, data=None, strict=True):

        self.spectrograph = load_spectrograph(spectrograph)
        self.par = self.spectrograph.default_pypeit_par() if par is None else par
        if not isinstance(self.par, PypeItPar):
            raise TypeError('Input parameter set must be of type PypeItPar.')
        self.bitmask = framematch.FrameTypeBitMask()
        self.table = table.Table(data if file_list is None 
                                 else self._build(file_list, strict=strict))
        # Instrument-specific validation of the header metadata. This
        # alters self.table in place!
        self.spectrograph.validate_metadata(self.table)
    
    def _build(self, file_list, strict=True):
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
        added_by_fsort = ['frametype', 'framebit']
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

                # Convert the time to hours
                # TODO: Done here or as a method in Spectrograph?
                if key == 'time' and value != 'None':
                    value = self.convert_time(value)
    
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

    def convert_time(self, in_time):
        """
        Convert the time read from a file header to hours for all
        spectrographs.
    
        Args:
            in_time (str):
                The time read from the file header

        Returns:
            float: The time in hours.
        """
        # Convert seconds to hours
        if self.spectrograph.timeunit == 's':
            return float(in_time)/3600.0
    
        # Convert minutes to hours
        if self.spectrograph.timeunit == 'm':
            return float(in_time)/60.0

        # Convert from an astropy.Time format
        if self.spectrograph.timeunit in time.Time.FORMATS.keys():
            ival = float(in_time) if self.spectrograph.timeunit == 'mjd' else in_time
            tval = time.Time(ival, scale='tt', format=self.spectrograph.timeunit)
            # Put MJD in hours
            return tval.mjd * 24.0
        
        msgs.error('Bad time unit')

    def find_frames(self, ftype, sci_ID=None, index=False):
        """
        Find the rows with the associated frame type.

        If the index is provided, the frames must also be matched the
        relevant science frame.

        Args:
            ftype (str):
                The frame type identifier.  See the keys for
                :class:`pypeit.core.framematch.FrameTypeBitMask`.
            sci_ID (:obj:`int`, optional):
                Index of the science frame that it must match.  If None,
                any row of the specified frame type is included.

        Returns:
            numpy.ndarray: Boolean array with the rows that contain the
            appropriate frames matched to the science frame, if
            provided.
        """
        if 'framebit' not in self.keys():
            raise ValueError('Frame types are not set.  First run get_frame_types.')
        if ftype is 'unknown':
            return self['framebit'] == 0
        indx = self.bitmask.flagged(self['framebit'], ftype)
        if sci_ID is not None:
            indx &= (self['sci_ID'] & sci_ID > 0)
        return np.where(indx)[0] if index else indx

    def find_frame_files(self, ftype, sci_ID=None):
        """
        Return the list of files with a given frame type.

        The frames must also match the science frame index, if it is
        provided.

        Args:
            ftype (str):
                The frame type identifier.  See the keys for
                :class:`pypeit.core.framematch.FrameTypeBitMask`.
            sci_ID (:obj:`int`, optional):
                Index of the science frame that it must match.  If None,
                any row of the specified frame type is included.

        Returns:
            list: List of file paths that match the frame type and
            science frame ID, if the latter is provided.
        """
        indx = self.find_frames(ftype, sci_ID=sci_ID)
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
        t = table.Table({'frametype':self.bitmask.type_names(type_bits), 'framebit':type_bits})
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
        self['framebit'][indx] = self.bitmask.turn_on(self['framebit'][indx], flag=frame_type)
        self['frametype'][indx] = self.bitmask.type_names(self['framebit'][indx])

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
            del self.table['frametype']
            del self.table['framebit']
        if useIDname and 'idname' not in self.keys():
            raise ValueError('idname is not set in table; cannot use it for file typing.')

        # Start
        msgs.info("Typing files")
        type_bits = np.zeros(len(self), dtype=self.bitmask.minimum_dtype())
    
        # Use the user-defined frame types from the input dictionary
        if user is not None:
            if len(user.keys()) != len(self):
                raise ValueError('The user-provided dictionary does not match table length.')
            msgs.info('Using user-provided frame types.')
            for ifile,ftypes in user.items():
                indx = self['filename'] == ifile
                type_bits[indx] = self.bitmask.turn_on(type_bits[indx], flag=ftypes.split(','))
            return self.set_frame_types(type_bits, merge=merge)
    
        # Loop over the frame types
        for i, ftype in enumerate(self.bitmask.keys()):
    
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
            type_bits[indx] = self.bitmask.turn_on(type_bits[indx], flag=ftype)
    
        # Find the nearest standard star to each science frame
        # TODO: Should this be 'standard' or 'science' or both?
        if 'ra' not in self.keys() or 'dec' not in self.keys():
            msgs.warn('Cannot associate standard with science frames without sky coordinates.')
        else:
            indx = self.bitmask.flagged(type_bits, flag='standard')
            for b, f, ra, dec in zip(type_bits[indx], self['filename'][indx], self['ra'][indx],
                                     self['dec'][indx]):
                if ra == 'None' or dec == 'None':
                    msgs.warn('RA and DEC must not be None for file:' + msgs.newline() + f)
                    msgs.warn('The above file could be a twilight flat frame that was'
                              + msgs.newline() + 'missed by the automatic identification.')
                    b = self.bitmask.turn_off(b, flag='standard')
                    continue
    
                # If an object exists within 20 arcmins of a listed standard,
                # then it is probably a standard star
                foundstd = flux.find_standard_file(ra, dec, check=True)
                b = self.bitmask.turn_off(b, flag='science' if foundstd else 'standard')
    
        # Find the files without any types
        indx = np.invert(self.bitmask.flagged(type_bits))
        if np.any(indx):
            msgs.info("Couldn't identify the following files:")
            for f in self['filename'][indx]:
                msgs.info(f)
            if not flag_unknown:
                msgs.error("Check these files before continuing")
    
        # Now identify the dark frames
        # TODO: This should not be here.  Move to instrument specific
        # selection as above
        indx = self.bitmask.flagged(type_bits, flag='bias') \
                        & (self['exptime'].data.astype(float) > self.spectrograph.minexp)
        type_bits[indx] = self.bitmask.turn_on(type_bits[indx], 'dark')
    
        # Finish up (note that this is called above if user is not None!)
        msgs.info("Typing completed!")
        return self.set_frame_types(type_bits, merge=merge)

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
        for c in col_order:
            if c not in self.keys():
                msgs.warn('{0} is not a valid column!  Removing from output.'.format(c))
                col_order.remove(c)
    
        # Write the output
        self.table[col_order].write(ofile, format=format, overwrite=overwrite)
    
    def match_to_science(self, calib_par, calwin, setup=False, verbose=True):
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
        # TODO: Failures is never set to True!  Iterate through science
        # frames and make sure they have an associated arc?
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
            for ftag in self.bitmask.keys():
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
                    code = framematch.match_warnings(calib_par, ftag, nmatch, numfr, target)
                    if code == 'break':
                        self['failure'][sci_idx] = True
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
    type_bits = np.zeros(len(fitstbl), dtype=fitstbl.bitmask.minimum_dtype())

    # Image typing
    if not notype:
        if spectrograph == 'shane_kast_blue':
            fitstbl['sci_ID'] = 1  # This links all the files to the science object
            type_bits[0] = fitstbl.bitmask.turn_on(type_bits[0], flag='bias')
            type_bits[1] = fitstbl.bitmask.turn_on(type_bits[1], flag='arc')
            type_bits[2:4] = fitstbl.bitmask.turn_on(type_bits[2:4], flag=['pixelflat', 'trace'])
            type_bits[4] = fitstbl.bitmask.turn_on(type_bits[4], flag='standard')
            type_bits[5:] = fitstbl.bitmask.turn_on(type_bits[5:], flag='science')
            fitstbl.set_frame_types(type_bits)

    return fitstbl


