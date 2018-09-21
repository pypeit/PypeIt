"""
Provides a class that handles the fits metadata required by PypeIt.
"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from __future__ import unicode_literals

import os
import re
import sys
import shutil

import numpy as np

import astropy.table
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units

from pypeit import msgs
from pypeit import utils
from pypeit.core.flux import find_standard_file
from pypeit import debugger

from pypeit.bitmask import BitMask
from pypeit.spectrographs.util import load_spectrograph

class FrameTypeBitMask(BitMask):
    """
    Define a bitmask to set the frame types.

    Frame types can be arc, bias, dark, pinhole, pixelflat, science,
    standard, or trace.
    """
    def __init__(self):
        frame_types = {         'arc': 'Arc lamp observation used for wavelength calibration',
                               'bias': 'Bias readout for detector bias subtraction',
                               'dark': 'Shuttered exposure to measure dark current',
                            'pinhole': 'Pinhole observation used for tracing slit centers',
                          'pixelflat': 'Flat-field exposure used for pixel-to-pixel response',
                            'science': 'On-sky observation of a primary target',
                           'standard': 'On-sky observation of a flux calibrator',
                              'trace': 'High-count exposure used to trace slit positions'
                      }
        super(FrameTypeBitMask, self).__init__(list(frame_types.keys()),
                                               descr=list(frame_types.values()))

    def type_names(self, type_bits, join=True):
        """
        Use the type bits to get the type names for each frame.

        .. todo::
            - This should probably be a general function in
              :class:`pypeit.bitmask.BitMask`
    
        Args:
            type_bits (int):
                The bit mask for each frame.
            bitmask (:class:`pypeit.bitmask.BitMask`, optional):
                The bit mask used to pull out the bit names.  Uses
                :class:`FrameTypeBitMask` by default.
            join (:obj:`bool`, optional):
                Instead of providing a list of type names for items with
                multiple bits tripped, joint the list into a single,
                comma-separated string.
    
        Returns:
            list: List of the frame types for each frame.  Each frame can
            have multiple types, meaning the 2nd axis is not necessarily the
            same length for all frames.
        """
        out = []
        for b in type_bits:
            n = self.flagged_bits(b)
            if len(n) == 0:
                n = ['None']
            out += [','.join(n)] if join else [n]
        return out
    
# Initially tried to subclass this from astropy.table.Table, but that
# proved too difficult.
class FitsMetaData:
    """
    Replacement for fitstbl

    Args:

    Attributes:
        spectrograph
        bitmask

    Raises:

    ---

    Create a table of relevant fits file metadata used during the
    reduction.

    The content of the fits table is dictated by the header keywords
    specified for the provided spectrograph.  It is expected that this
    table can be used to set the frame type of each file.

    The metadata is validated using checks specified by the provided
    spectrograph class.

    .. todo::
        This should get abstracted to be independent of the
        spectrograph, with all the spectrograph dependent keywords be
        passed into the function.  The validation of the table would
        then happen in PypeItSetup.

    Args:
        file_list (list):
            The list of files to include in the table.
        spectrograph
            (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph used to collect the data save to each file.
            The class is used to provide the header keyword data to
            include in the table and specify any validation checks.
            
    Returns:
        :class:`astropy.table.Table`: A table with the relevant metadata
        for the provided fits files.
    """
    def __init__(self, spectrograph, file_list=None, data=None, strict=True):
        self.spectrograph = spectrograph
        self.bitmask = FrameTypeBitMask()
        self.table = astropy.table.Table(data if file_list is None 
                                            else self._build_data(file_list, strict=strict))
    
    def _build_data(self, file_list, strict=True):
        """
        Returns a dictionary with the data to be included in the table.
        """
        # Get the header keywords specific to the provided spectrograph.
        # head_keys is a nested dictionary listing the header keywords
        # for each extension in the file.  The top-level dictionary key
        # is just the 0-indexed number of the extension
        head_keys = self.spectrograph.header_keys()

        # The table is declared based on this input dictionary: The
        # directory, filename, instrument, utc are always included
        data = {k:[] for k in FitsMetaData.default_keys()}

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
        for i in range(numfiles):

            # Read the fits headers
            headarr = self.spectrograph.get_headarr(file_list[i], strict=strict)

            # Check that the header is valid
            # TODO: The check_headers function needs to be implemented
            # for each instrument.  spectrograph.check_headers() should
            # raise an exception with an appropriate message.
            try:
                # TODO: Move this into spectrograph.validate_fitstbl()
                self.spectrograph.check_headers(headarr)
            except Exception as e:
                msgs.warn('Reading of headers from file:' + msgs.newline() + file_list[i]
                          + msgs.newline() + 'failed with the following exception'
                          + msgs.newline() + e.__repr__() + msgs.newline() +
                          'Please check that the file was taken with the provided instrument:'
                          + msgs.newline() + '{0}'.format(self.spectrograph.spectrograph)
                          + msgs.newline() + 'Then either change the instrument or remove/skip '
                          + 'the file.' + msgs.newline()+ 'Continuing by ignoring this file...')
                numfiles -= 1
                continue

            # Add the directory, file name, and instrument to the table
            d,f = os.path.split(file_list[i])
            data['directory'].append(d)
            data['filename'].append(f)
            data['instrume'].append(self.spectrograph.spectrograph)

            # Add the time of the observation
            utc = self.get_utc(headarr)
            data['utc'].append('None' if utc is None else utc)
            if utc is None:
                msgs.warn('UTC is not listed as a header keyword in file:' + msgs.newline()
                          + file_list[i])

            # TODO: Read binning-dependent detector properties here? (maybe read speed too)
    
            # Now get the rest of the keywords
            for k in data.keys():
                if k in FitsMetaData.default_keys():
                    continue
    
                # Try to read the header item
                try:
                    value = headarr[ext[k]][head_keys[ext[k]][k]]
                except KeyError:
                    # Keyword not found in header
                    msgs.warn("{:s} keyword not in header. Setting to None".format(k))
                    value = 'None'
    
                # Convert the time to hours
                # TODO: Done here or as a method in Spectrograph?
                if k == 'time' and value != 'None':
                    value = self.convert_time(value)
    
                # Set the value
                vtype = type(value)
                if np.issubdtype(vtype, str):
                    value = value.strip()
                if np.issubdtype(vtype, np.integer) or np.issubdtype(vtype, np.floating) \
                        or np.issubdtype(vtype, str) or np.issubdtype(vtype, np.bool_):
                    data[k].append(value)
                else:
                    msgs.bug('Unexpected type, {1:s}, for key {0:s}'.format(k,
                             vtype).replace('<type ','').replace('>',''))
    
            msgs.info('Successfully loaded headers for file:' + msgs.newline() + file_list[i])

        # Report
        msgs.info("Headers loaded for {0:d} files successfully".format(numfiles))
        if numfiles != len(file_list):
            msgs.warn("Headers were not loaded for {0:d} files".format(len(file_list) - numfiles))
        if numfiles == 0:
            msgs.error("The headers could not be read from the input data files." + msgs.newline() +
                    "Please check that the settings file matches the data.")
        
        return data

    # TODO:  In this implementation, slicing the FitsMetaData object
    # will return an astropy.table.Table, not a FitsMetaData object.
    def __getitem__(self, item):
        return self.table.__getitem__(item)

    def __setitem__(self, item):
        return self.table.__setitem__(item)

    def __len__(self):
        return self.table.__len__()

    def __repr__(self):
        return self.table._base_repr_(html=False,
                            descr_vals=['FitsMetaData:\n',
                                        '              spectrograph={0}\n'.format(
                                                                    self.spectrograph.spectrograph),
                                        '              length={0}\n'.format(len(self))])

    def _repr_html_(self):
        return self.table._base_repr_(html=True, max_width=-1,
                            descr_vals=['FitsMetaData: spectrograph={0}, length={1}\n'.format(
                                                    self.spectrograph.spectrograph, len(self))])

    @staticmethod
    def default_keys():
        return [ 'directory', 'filename', 'instrume', 'utc' ]

    @staticmethod
    def get_utc(headarr):
        """
        Find and return the UTC for a file based on the headers read from
        all extensions.
    
        The value returned is the first UT or UTC keyword found any in any
        header object.
    
        Args:
            headarr (list):
                List of :obj:`astropy.io.fits.Header` objects to search
                through for a UTC of the observation.
        Returns:
            object: The value of the header keyword.
        """
        for h in headarr:
            if h == 'None':
                continue
            if 'UTC' in h.keys():
                return h['UTC']
            elif 'UT' in h.keys():
                return h['UT']
        return None

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
        if self.spectrograph.timeunit in Time.FORMATS.keys():
            ival = float(in_time) if self.spectrograph.timeunit == 'mjd' else in_time
            tval = Time(ival, scale='tt', format=self.spectrograph.timeunit)
            # Put MJD in hours
            return tval.mjd * 24.0
        
        msgs.error('Bad time unit')

    def find_frames(self, ftype, sci_ID=None):
        """
        Find the rows with the associated frame type.

        The frames must also match the science frame index, if it is
        provided.

        Args:
            ftype (str):
                The frame type identifier.  See the keys for
                :class:`FrameTypeBitMask`.

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
        indx = self.bitmask.flagged(self['framebit'], ftype)
        return indx if sci_ID is None else indx & (self['sci_ID'] == sci_ID)

    def find_frame_files(self, ftype, sci_ID=None):
        """
        Return the list of files with a given frame type.

        The frames must also match the science frame index, if it is
        provided.

        Args:
            ftype (str):
                The frame type identifier.  See the keys for
                :class:`FrameTypeBitMask`.
            sci_ID (:obj:`int`, optional):
                Index of the science frame that it must match.  If None,
                any row of the specified frame type is included.

        Returns:
            list: List of file paths that match the frame type and
            science frame ID, if the latter is provided.
        """
        indx = self.find_frames(ftype, sci_ID=sci_ID)
        return [ os.path.join(d,f) for d,f in zip(self['directory'][indx], self['filename'][indx])]

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
            type name and bits.  Nothing is returned if merge is True.
        """
        t = astropy.table.Table({'frametype': get_type_names(type_bits, bitmask=bm), 
                                 'framebit': type_bits})
        if merge:
            self['frametype'] = t['frametype']
            self['framebit'] = t['framebit']
            return
        return t

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
            type names and the type bits.  See :class:`FrameTypeBitMask`
            for the allowed frame types.
        """
        # Checks
        if 'frametype' in self.keys() or 'framebit' in self.keys():
            msgs.warn('Removing existing frametype and framebit columns.')
            del self['frametype']
            del self['framebit']
        if useIDname and 'idname' not in self.keys():
            raise ValueError('idname is not set in table; cannot use it for file typing.')

        # Start
        msgs.info("Typing files")
        type_bits = np.zeros(self.__len__(), dtype=self.bitmask.minimum_dtype())
    
        # Use the user-defined frame types from the input dictionary
        if user is not None:
            if len(user.keys()) != self.__len__():
                raise ValueError('The user-provided dictionary does not match table length.')
            msgs.info('Using user-provided frame types.')
            for ifile,ftypes in ftdict.items():
                indx = self['filename'] == ifile
                type_bits[indx] = self.bitmask.turn_on(type_bits[indx], flag=ftypes.split(','))
            return self.set_frame_types(type_bits, merge=merge)
    
        # Loop over the frame types
        for i, ftype in enumerate(self.bitmask.keys()):
    
            # Initialize: Flag frames with the correct ID name or start by
            # flagging all as true
            indx = self['idname'] == self.spectrograph.idname(ftype) if useIDname \
                        else np.ones(self.__len__(), dtype=bool)
    
            # Include a combination of instrument-specific checks using
            # combinations of the full set of metadata
            indx &= self.spectrograph.check_ftype(ftype, self)
    
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
                foundstd = find_standard_file(ra, dec, check=True)
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
        # TODO: Move this to Spectrograph.check_ftype?
        indx = self.bitmask.flagged(type_bits, flag='bias') \
                        & (self['exptime'].data.astype(float) > self.spectrograph.minexp)
        type_bits[indx] = self.bitmask.turn_on(type_bits[indx], 'dark')
    
        # Finish up (note that this is called above if user is not None!)
        msgs.info("Typing completed!")
        return self.set_frame_types(type_bits, merge=merge)

    def write(self, ofile, columns=None, format=None):
        """
        Write the metadata for the files to reduce.
    
        The table is written with the filename and frametype columns
        first.  All remaining columns, or a subset of them selected by
        `columns`, follow these first two.
    
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
        
        Raises:
            ValueError:
                Raised if the columns to include are not unique.
        """
        msgs.info('Writing fits file metadata to {0}.'.format(ofile))
    
        # Set the columns to include and check that they are unique
        _columns = list(self.keys()) if columns is None else columns
        if len(np.unique(_columns)) != len(_columns):
            # TODO: A warning may suffice...
            raise ValueError('Column names must be unique!')
    
        # Force the filename and frametype columns to go first
        col_order = [ 'filename', 'frametype' ]
        col_order.append(list(set(_columns) - set(col_order)))
    
        # Remove any columns that don't exist
        for c in col_order:
            if c not in self.keys():
                msgs.warn('{0} is not a valid column!  Removing from output.'.format(c))
                col_order.remove(c)
    
        # Write the output
        self[col_order].write(ofile, format=format)
#        'ascii.fixed_width')
    

def chk_all_conditions(fitstbl, cond_dict):
    """ Loop on the conditions for this given file type

    Parameters
    ----------
    fitstbl : Table
    cond_dict : dict

    Returns
    -------
    gd_chk : ndarray (bool)
      True = Passes all checks
    """
    gd_chk = np.ones(len(fitstbl), dtype=bool)
    # Loop on the items to check
    chkk = cond_dict.keys()
    for ch in chkk:
        if ch[0:9] == 'condition':
            # Deal with a conditional argument
            conds = re.split("(\||\&)", cond_dict[ch])
            ntmp = chk_condition(fitstbl, conds[0])
            # And more
            for cn in range((len(conds)-1)//2):
                if conds[2*cn+1] == "|":
                    ntmp = ntmp | chk_condition(fitstbl, conds[2*cn+2])
                elif conds[2*cn+1] == "&":
                    ntmp = ntmp & chk_condition(fitstbl, conds[2*cn+2])
            gd_chk = gd_chk & ntmp
        else:
            if fitstbl[ch].dtype.char in ['S','U']:  # Numpy string array
                # Strip numpy string array of all whitespace
                gd_chk = gd_chk & (np.char.strip(fitstbl[ch]) == cond_dict[ch])
            else:
                gd_chk = gd_chk & (fitstbl[ch] == cond_dict[ch])
    # Return
    return gd_chk


def chk_condition(fitstbl, cond):
    """
    Code to perform condition.  A bit messy so a separate definition
    was generated.
    Create an exposure class for every science frame

    Parameters
    ----------
    fitsdict : Table
      Contains relevant information from fits header files
    cond : str
      A user-specified condition that is used to identify filetypes.
      This string is the fourth argument of the frame conditions that
      is specified in the settings file. For example, in the line:
      'bias check condition1 exptime=0'
      cond = 'exptime=0'

    Returns
    -------
    ntmp: bool array
      A boolean array of all frames that satisfy the input condition
    """
    if "<=" in cond:
        tcond = cond.split("<=")
        ntmp = fitstbl[tcond[0]] <= float(tcond[1])
    elif ">=" in cond:
        tcond = cond.split(">=")
        ntmp = fitstbl[tcond[0]] >= float(tcond[1])
    elif "!=" in cond:
        tcond = cond.split("!=")
        if 'int' in fitstbl[tcond[0]].dtype.name:
            ntmp = fitstbl[tcond[0]] != int(tcond[1])
        elif 'float' in fitstbl[tcond[0]].dtype.name:
            ntmp = fitstbl[tcond[0]] != float(tcond[1])
        else:
            ntmp = fitstbl[tcond[0]] != tcond[1]
    elif "<" in cond:
        tcond = cond.split("<")
        ntmp = fitstbl[tcond[0]] < float(tcond[1])
    elif ">" in cond:
        tcond = cond.split(">")
        ntmp = fitstbl[tcond[0]] > float(tcond[1])
    elif "=" in cond:
        tcond = cond.split("=")
        if 'int' in fitstbl[tcond[0]].dtype.name:
            ntmp = fitstbl[tcond[0]] == int(tcond[1])
        elif 'float' in fitstbl[tcond[0]].dtype.name:
            ntmp = fitstbl[tcond[0]] == float(tcond[1])
        else:
            ntmp = fitstbl[tcond[0]] == tcond[1]
    else:
        ntmp = None
    return ntmp



# TODO: Does this function work?  while statement seems weird
def match_ABBA(fitstbl, max_targ_sep=30, max_nod_sep=2):
    """

    Parameters
    ----------
    fitstbl : Table
        Contains relevant information from fits header files
    Returns
    -------
    fitstbl : Table (w/new column)

    """
    # Make coords for all observations, including calibration frames
    coords = SkyCoord(ra=fitstbl['ra'], dec=fitstbl['dec'], unit='deg')

    # Indices of files classified as science frames
    sci_idx = np.where(fitstbl['science'])[0]

    # Create a mask, currently all True
    mask = np.ones(np.sum(len(fitstbl)), dtype=bool)

    # Create empty dictionary of targets
    targets = {}

    # Set science frames to False in the mask
    mask[sci_idx] = False

    while np.any(~mask):
        idx = np.where(~mask)[0][0]
        # Find indices corresponding to any others files that match the coordinates of current science frame in consideration
        match = coords[idx].separation(coords[sci_idx]).arcsec < max_targ_sep # 30 arcsec separation is arbitrary (from 1e-2 deg) for target separation
        # Add this science target to the dictionary
        targ_name = fitstbl[idx]['target']  # Grab the target name from fitstbl
        # Create that target in dict and add corresponding indices that point to the relevant science frames
        targets[targ_name] = sci_idx[match]
        # Turn mask 'off' (i.e., to True) now that they have been grouped to this target
        mask[targets[targ_name]] = True

    # Create an empty list of filenames that will store each frame's buddy A/B frame
    AB_frame = [''] * len(fitstbl)

    for key, value in targets.items():

        files = fitstbl['filename'][value]

        # Check here that there are more than 1 files and that # of files is even
        if len(files) == 1:
            msgs.warn('Cannot perform NIR A-B reduction on targets with 1 file')
        elif len(files) % 2 != 0:
            msgs.warn('Expected an even number of files associated with target ' + key)
        #### Check for increasing time? Files are read in numerical sequential order -- should be in order of increasing time anyway..

        # Assume that the files are initially in ABBA order and proceed
        ABBA_coords = coords[value]

        # Break files into ABBA groups (includes 'remainder' if there are only 2 files)
        ABBA_groups = [ABBA_coords[i:i + 4] for i in range(0, len(ABBA_coords), 4)]
        value_groups = [value[i:i + 4] for i in range(0, len(ABBA_coords), 4)]

        for group in range(len(ABBA_groups)):
            # Warn user that if there are any groups of only 2 files, assuming they are in order of A and B
            if len(ABBA_groups[group]) == 2:
                msgs.info('Assuming these two frames are A and B frame')
            # Check that we have a 4-file, ABBA sequence
            elif len(ABBA_groups[group]) == 4:
                # Check that frames 1, 4 of an ABBA sequence are at the same nod position (A) based on their RA, DEC
                if ABBA_coords[0].separation(ABBA_coords[-1]).arcsec > max_nod_sep: # 5e-4 deg --> 1.8 arcsec --> ~2 arcsec set to be max nod separation
                    msgs.info('Are frames ' + str((group * 4) + 1) + ' and ' + str(
                        (group * 4) + 4) + ' for target ' + key + ' both A frames?')
                # Check that frames 2, 3 of an ABBA sequence are both at the same nod position (B)
                if ABBA_coords[1].separation(ABBA_coords[2]).arcsec > max_nod_sep:
                    msgs.info('Are frames ' + str((group * 4) + 2) + ' and ' + str(
                        (group * 4) + 3) + ' for target ' + key + ' both B frames?')
            else:
                msgs.error('Check number of frames for this target -- files are not grouped in ABBA or AB')

            # Create a copy of the array value_groups[group] (which gives the indices corresponding to the 4/2 ABBA/AB files in consideration)
            AB_idx_flip = np.copy(value_groups[group])
            # Flip such that, for example, (1, 2, 3, 4) --> (2, 1, 4, 3)
            AB_idx_flip[::2], AB_idx_flip[1::2] = value_groups[group][1::2], value_groups[group][::2]

            # Fill in AB_frame list
            for i in range(len(value_groups[group])):
                AB_frame[value_groups[group][i]] = fitstbl['filename'][AB_idx_flip[i]]

    fitstbl['AB_frame'] = AB_frame

    return fitstbl


def match_logic(ch, tmtch, fitstbl, idx):
    """ Perform logic on matching with fitsdict
    Parameters
    ----------
    ch : str
      Header card alias, eg. exptime
    tmtch : str
      Defines the logic
      any
      ''
      >, <, >=, <=, =, !=
      If tmtch begins with a "|", the match compares to the science frame
      else the value is added to the science frame
    fitstbl : Table
    idx : int
      Science index

    Returns
    -------
    w : ndarray, bool
      True/False for the rows in fitstbl satisfying the condition
    """
    if tmtch == "any":   # Anything goes
        w = np.ones_like(fitstbl, dtype=bool)
    elif tmtch == '':  # Header value must match that of science
        w = fitstbl[ch] == fitstbl[ch][idx]
    elif tmtch[0] in ['=','<','>','|']: # Numerics
        mtch = np.float64(fitstbl[ch][idx]) + float(
            ''.join(c for c in tmtch if c not in ['=', '<', '>', '|']))
        operand = ''.join(c for c in tmtch if c in ['=', '<', '>'])
        if operand == '=':
            operand += '='
        #
        if tmtch[0] != '|':
            w = eval('fitstbl[ch].data.astype(np.float64) {:s} {:f}'.format(operand, mtch))
        else:
            w = eval('np.abs(fitstbl[ch].data.astype(np.float64) - np.float64(fitstbl[ch][idx])) {:s} {:f}'.format(operand, mtch))
    elif tmtch[0:2] == '%,':  # Splitting a header keyword
        splcom = tmtch.split(',')
        debugger.set_trace()
        spltxt, argtxt, valtxt = splcom[1], np.int(splcom[2]), splcom[3]
        tspl = []
        for sp in fitstbl[ch]:
            tmpspl = str(re.escape(spltxt)).replace("\\|", "|")
            tmpspl = re.split(tmpspl, sp)
            if len(tmpspl) < argtxt+1:
                tspl.append("-9999999")
            else:
                tspl.append(tmpspl[argtxt])
        tspl = np.array(tspl)
        #                        debugger.set_trace()
        tmpspl = str(re.escape(spltxt)).replace("\\|", "|")
        tmpspl = re.split(tmpspl, fitstbl[ch][idx])
        msgs.warn("HAS NOT BEEN DEVELOPED SINCE THE SetupClass refactor;  no test case..")
        debugger.set_trace()  # HAS NOT BEEN DEVELOPED SINCE THE SetupClass refactor;  no test case..
        if len(tmpspl) < argtxt + 1:
            return None
        else:
            scispl = tmpspl[argtxt]
        if valtxt == "''":
            w = np.where(tspl == scispl)[0]
        elif valtxt[0] == '=':
            mtch = np.float64(scispl) + np.float64(valtxt[1:])
            w = np.where(tspl.astype(np.float64) == mtch)[0]
        elif valtxt[0] == '<':
            if valtxt[1] == '=':
                mtch = np.float64(scispl) + np.float64(valtxt[2:])
                w = np.where(tspl.astype(np.float64) <= mtch)[0]
            else:
                mtch = np.float64(scispl) + np.float64(valtxt[1:])
                w = np.where(tspl.astype(np.float64) < mtch)[0]
        elif valtxt[0] == '>':
            if valtxt[1] == '=':
                mtch = np.float64(scispl) + np.float64(valtxt[2:])
                w = np.where(tspl.astype(np.float64) >= mtch)[0]
            else:
                mtch = np.float64(scispl) + np.float64(valtxt[1:])
                w = np.where(tspl.astype(np.float64) > mtch)[0]
    # Return
    return w


def match_to_science(calib_par, match_dict, fitstbl, calwin, setup=False, verbose=True,
                     match_nods=False):
    """
    For a given set of identified data, match calibration frames to science frames

    Parameters
    ----------
    fitstbl : Table
      Contains relevant information from fits header files
    settings_spect : dict
    setup: bool, optional
      Running in setup mode?
    flux_calibrate : bool, optional
      Do checks related to flux calibration
    wave_calib : str
    calwin : float

    Returns
    -------
    fitstbl : Table
      Updated with failures and sci_ID columns
    """
    msgs.info("Matching calibrations to Science frames")

    # New columns
    fitstbl['failures'] = False
    fitstbl['sci_ID'] = 0

    # Frame type masking
    bm = FrameTypeBitMask()

    # Loop on science frames
    for ss, sci_idx in enumerate(np.where(fitstbl['science'])[0]):
        msgs.info("=================================================")
        msgs.info("Matching calibrations to {:s}: {:s}".format(
                fitstbl['target'][sci_idx], fitstbl['filename'][sci_idx]))

        # Science ID (trivial but key to the bit-wise that follows)
        fitstbl['sci_ID'][sci_idx] = 2**ss

        # Find matching (and nearby) calibration frames
        for ftag in bm.keys():
            # Skip science frames
            if ftag == 'science':
                continue

            # bias/dark check to make sure we need to find matching frames
            if ftag == 'dark' and calib_par['biasframe']['useframe'] != 'dark':
                msgs.info("  Dark frames not required.  Not matching..")
                continue
            if ftag == 'bias' and calib_par['biasframe']['useframe'] != 'bias' \
                        and not calib_par['badpix']:
                msgs.info("  Bias frames not required.  Not matching..")
                continue

            # Find the files of this type
            gd_match = bm.flagged(fitstbl['framebit'], flag=ftag)

            # How many matching frames are required?  This is instrument specific
            numfr = (1 if ftag == 'arc' else 0) if setup \
                        else calib_par['{0}frame'.format(ftag)]['number']

            # If not required and not doing setup, continue
            if numfr == 0 and not setup:
                msgs.info("   No {0:s} frames are required.  Not matching..".format(ftag))
                continue

            # Now go ahead and match the frames
            if 'match' not in match_dict[ftag].keys() and (not setup):
                msgs.error("Need match criteria for {0:s}!!".format(ftag))
            elif 'match' not in match_dict[ftag].keys():
                msgs.info("No matching criteria for {0:s} frames with this instrument".format(ftag))
            else:
                chkk = match_dict[ftag]['match'].keys()
                for ch in chkk:
                    tmtch = match_dict[ftag]['match'][ch]
                    gd_match &= match_logic(ch, tmtch, fitstbl, sci_idx)

            # Find the time difference between the calibrations and science frames
            if calwin > 0.0 and 'time' in fitstbl.keys():
                tdiff = np.abs(fitstbl['time']-fitstbl['time'][sci_idx])
                gd_match &= tdiff <= calwin

            # Now find which of the remaining n are the appropriate calibration frames
            nmatch = np.sum(gd_match)
            if verbose and 'target' in fitstbl.keys():
                msgs.info('  Found {0:d} {1:s} frame for {2:s} ({3:d} required)'.format(
                                nmatch, ftag, fitstbl['target'][sci_idx], numfr))

            # Have we identified enough of these calibration frames to continue?
            if nmatch < np.abs(numfr):
                code = match_warnings(calib_par, ftag, nmatch, numfr,
                                      fitstbl['target'][sci_idx] if 'target' in fitstbl.keys()
                                            else 'unknown')
                if code == 'break':
                    fitstbl['failure'][sci_idx] = True
                    fitstbl['sci_ID'][sci_idx] = -1  # This might break things but who knows..
                    break
            else:
                wa = np.where(gd_match)[0]
                # Select the closest calibration frames to the science frame
                tdiff = np.abs(fitstbl['time'][wa]-fitstbl['time'][sci_idx])
                wa = wa[np.argsort(tdiff)]
                #if ftag == 'bias':
                #    debugger.set_trace()
                if setup or (numfr < 0):
                    fitstbl['sci_ID'][wa] |= 2**ss  # Flip the switch (if need be)
                else:
                    fitstbl['sci_ID'][wa[:numfr]] |= 2**ss  # Flip the switch (if need be)

    msgs.info("Science frames successfully matched to calibration frames")

    # Return with nods matched if requested
    return match_ABBA(fitstbl) if match_nods else fitstbl

#    # How to do this if statement only if '--custom' is on?
#    if spectrograph.spectrograph == 'keck_nirspec':
#        fitstbl = match_ABBA(fitstbl)


def insufficient_frame_error(frametype):
    msgs.error('Insufficient {0} frames found. Include more frames, '.format(frametype)
                + 'reduce the required amount by setting'
                + msgs.newline() + '[calibrations]'
                + msgs.newline() + '    [[{0}frame]]'.format(frametype)
                + msgs.newline() + '        number = XX'
                + msgs.newline() + 'in the pypeit file, or specify a specific'
                + 'pixelflat file by setting'
                + msgs.newline() + '[calibrations]'
                + msgs.newline() + '    [[{0}frame]]'.format(frametype)
                + msgs.newline() + '        useframe = XX'
                + msgs.newline() + 'in the pypeit file')


def match_warnings(calib_par, ftag, nmatch, numfr, target, setup=False):
    """
    Provide match warnings

    Parameters
    ----------
    ftag : str
      frametype, e.g. bias
    nmatch : int
    numfr : int
    target : str
      Name of the target
    settings_argflag : dict

    Returns
    -------
    code : str
      'None' = no further action required
    """
    code = 'None'
    msgs.warn("  Only {0:d}/{1:d} {2:s} frames for {3:s}".format(nmatch, numfr, ftag, target))

    # TODO: Why does number of pixelflat, trace, and standard not matter
    # if you're not flat-fielding the data?  Particularly for trace...
    flatfield = calib_par['flatfield']['method'] is not None

    # Errors for insufficient BIAS frames
    if calib_par['biasframe']['useframe'].lower() == ftag:
        insufficient_frame_error(ftag)

    # Errors for insufficient PIXELFLAT frames
    if ftag == 'pixelflat' and flatfield and calib_par['flatfield']['frame'] == 'pixelflat':
        if calib_par['masters'] == 'force':
            msgs.warn('Fewer {0} frames than expected for {1}'.format(ftag, target)
                      +', but will use MasterFrames.')
        else:
            insufficient_frame_error(ftag)

    # Errors for insufficient PINHOLE frames
    if ftag == 'pinhole':
        insufficient_frame_error(ftag)

    # Errors for insufficient TRACE frames
    if ftag == 'trace' and flatfield:
        if calib_par['masters'] == 'force':
            msgs.warn('Fewer {0} frames than expected for {1}'.format(ftag, target)
                      +', but will use MasterFrames.')
        else:
            insufficient_frame_error(ftag)

    # Errors for insufficient standard frames
    if ftag == 'standard' and flatfield:
        if calib_par['masters'] == 'force':
            msgs.warn('No {0} frames for {1}'.format(ftag, target)
                      + ', but will use MasterFrames.')
        else:
            insufficient_frame_error(ftag)

    # Errors for insufficient ARC frames
    if ftag == 'arc' and (calib_par['wavelengths']['reference'] not in ['pixel', 'sky']):
        if setup:
            msgs.warn('No {0} frame for {1}. '.format(ftag, target)
                      + 'Removing it from list of science frames.  Add an arc and rerun if '
                      + 'you wish to reduce this with PYPIT!!')
            return 'break'
        elif calib_par['masters'] == 'force':
            msgs.warn('No {0} frames for {1}'.format(ftag, target)
                      + ', but will use MasterFrames.')
        else:
            insufficient_frame_error(ftag)

    return code


def make_dirs(spectrograph, caldir, scidir, qadir, overwrite=False):
    """
    Make the directories for the pypeit output.

    .. todo::
        I think this should just fault if the directories exist and
        `overwrite` is False.

    Args:
        spectrograph (str):
            The name of the spectrograph that provided the data to be
            reduced.
        caldir (str):
            The directory to use for saving the master calibration
            frames.
        scidir (str):
            The directory to use for the main reduction output files.
        qadir (str):
            The directory to use for the quality assessment output.
        overwrite(:obj:`bool`, optional):
            Flag to overwrite any existing files/directories.
    """

    # First, get the current working directory
    currDIR = os.getcwd()
    newdir = "{0:s}/{1:s}".format(currDIR, scidir)
    msgs.info('Creating Science directory:' + msgs.newline() + newdir)
    if os.path.exists(newdir):
        msgs.warn('Directory already exists.  Files may be ovewritten.')
    else:
        os.mkdir(newdir)

    # Create a directory where all of the master calibration frames are stored.
    newdir = "{:s}/{:s}_{:s}".format(currDIR, caldir, spectrograph)
    msgs.info('Creating Master Calibrations directory:' + msgs.newline() + newdir)
    if os.path.exists(newdir):
        msgs.warn('Directory already exists.  Files may be ovewritten.')
    else:
        os.mkdir(newdir)

    # Create a directory where all of the QA is stored
    newdir = "{0:s}/{1:s}".format(currDIR, qadir)
    msgs.info('Creating QA directory:' + msgs.newline() + newdir)
    if os.path.exists(newdir):
        msgs.warn('Directory already exists.  Files may be ovewritten.')
    else:
        os.mkdir(newdir)
        os.mkdir(newdir+'/PNGs')


def dummy_fitstbl(nfile=10, spectrograph='shane_kast_blue', directory='./', notype=False):
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
    fitstbl : Table

    """
    fitsdict = dict({'directory': [], 'filename': [], 'utc': []})
    fitsdict['utc'] = ['2015-01-23']*nfile
    fitsdict['directory'] = [directory]*nfile
    fitsdict['filename'] = ['b{:03d}.fits.gz'.format(i) for i in range(nfile)]
    fitsdict['date'] = ['2015-01-23T00:{:02d}:11.04'.format(i) for i in range(nfile)]  # Will fail at 60
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
    #
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
    elif spectrograph == 'none':
        pass
    # arrays
    for k in fitsdict.keys():
        fitsdict[k] = np.array(fitsdict[k])
    # Table me
    fitstbl = Table(fitsdict)
    fitstbl['instrume'] = spectrograph
    # Image typing
    if not notype:
        for ftype in ftype_list:
            fitstbl[ftype] = np.zeros(len(fitstbl), dtype=bool)
        if spectrograph == 'shane_kast_blue':
            fitstbl['sci_ID'] = 1  # This links all the files to the science object
            fitstbl['bias'][0] = True
            fitstbl['arc'][1] = True
            fitstbl['trace'][2:4] = True
            fitstbl['pixelflat'][2:4] = True
            fitstbl['standard'][4] = True
            fitstbl['science'][5:] = True
    # Return
    return FitsMetaData(load_spectrograph(spectrograph), data=fitstbl)


