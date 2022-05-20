""" Class for I/O of PypeIt file"""
import os
import numpy as np
import yaml

import configobj

from astropy.table import Table

from pypeit import msgs

from IPython import embed

class PypeItFile:
    """
    Class to load, process, and write PypeIt files

    Args:
        config (:obj:`dict` or :obj:`list`):
            Configuration lines or dict
        data_table (:class:`astropy.table.Table`):
            Data block
        setup (:obj:`dict`, optional):
            dict defining the Setup
            Defaults to A and otherwise a null description
    """
    def __init__(self, config:dict, data_table:Table,
                 setup:dict=None):
        # Load up
        self.data = data_table
        self.setup = setup

        # Load up config
        self.config = config
        
    @classmethod
    def from_file(cls, pypeit_file:str): 
        """
        Parse the user-provided .pypeit reduction file.

        Args:
            ifile (:obj:`str`):
                Name of pypeit file

        #Returns:
        #    :obj:`tuple`:  Provides (1) a list of configuration lines, (2) a list of
        #    datafiles to read, (3) a list of frametypes for each file, (4) an
        #    `astropy.table.Table`_ with the user supplied metadata for the data
        #    files, (5) a list of setup lines, and (6) a instrument configuration
        #    (setup) dictionary.
        Returns:
            pypeItFile (:class:`PypeItFile`):
        """
        # Read in the pypeit reduction file
        msgs.info('Loading the reduction file')
        lines = _read_pypeit_file_lines(pypeit_file)

        # Used to select the configuration lines: Anything that isn't part
        # of the data or setup blocks is assumed to be part of the
        # configuration
        is_config = np.ones(len(lines), dtype=bool)

        # Parse data block
        s, e = _find_block(lines, 'data')
        if s >= 0 and e < 0:
            msgs.error("Missing 'data end' in {0}".format(pypeit_file))
        if s < 0:
            msgs.error("You haven't specified any data!")
        usrtbl = cls._read_data_file_table(lines[s:e])
        is_config[s-1:e+1] = False

        # Parse the setup block
        s, e = _find_block(lines, 'setup')
        if s >= 0 and e < 0:
            msgs.error(f"Missing 'setup end' in {pypeit_file}")
        elif s < 0:
            msgs.error(f"Missing 'setup read' in {pypeit_file}")

        # Proceed
        setups, sdict = cls._parse_setup_lines(lines[s:e])
        is_config[s-1:e+1] = False

        # Build the configobj 
        confObj = configobj.ConfigObj(list(lines[is_config]))
        confdict = confObj.dict()

        # vet
        for key in ['filename', 'frametype']:
            if key not in usrtbl.keys():
                msgs.error("Add {:s} to your PypeIt file before using run_pypeit".format(key))
        if 'rdx' not in confdict.keys() or 'spectrograph' not in confdict['rdx'].keys():
            msgs.error(f"Missing spectrograph in the Parameter block of your PypeIt file.  Add it!")
        msgs.info('PypeIt file loaded successfully')

        # Instantiate
        slf = cls(confdict, usrtbl, setup=sdict)

        return slf

    @staticmethod
    def _parse_setup_lines(lines):
        """
        Return a list of the setup names and corresponding dict

        Args:
            lines (`numpy.ndarray`_): Setup lines as an array

        Returns:
            tuple: list, dict

        """
        setups = []

        # Kludge for backwards compatability
        line_list = lines.tolist()
        for ss, line in enumerate(line_list):
            if 'Setup' in line and ':' not in line:
                line_list[ss] = line+':'

        # Slurp
        ystr = '\n'.join(line_list)
        sdict = yaml.safe_load(ystr)
        for key in sdict:
            if 'Setup' in key:
                tsetup = key.split()[1].strip()
                setups.append(tsetup)

        # Check
        if len(setups) > 1:
            msgs.error("Setup block contains more than one Setup!")
        elif len(setups) != 1:
            msgs.error("Add setup info to your PypeIt file in the setup block!")

        return setups, sdict

    @staticmethod
    def _read_data_file_table(lines):
        """
        Read the file table format.
        
        Args:
            lines (:obj:`list`):
                List of lines *within the data* block read from the pypeit
                file.
        
        Returns:
            list, dict, Table:  Returns the list of data file names, a
            dictionary with the frame types of each file where the key of
            the dictionary is the file name, and a Table with the data
            provided in the pypeit file.  Note that the files listed in the
            first object contain the full path, whereas the file names in
            the frame type dictionary and the data table do not include the
            full path to the file.

            PypeItError:
                Raised if file does not have a 'filename' or 'frametype' column.
        """

        # Allow for multiple paths
        paths = []
        for l in lines:
            space_ind = l.index(" ")
            if l[:space_ind].strip() != 'path':
                break
            paths += [ l[space_ind+1:] ]

        npaths = len(paths)
        header = [ l.strip() for l in lines[npaths].split('|') ][1:-1]

        # Minimum columns required
        if 'filename' not in header:
            msgs.error('Table format failure: No \'filename\' column.')
        if 'frametype' not in header:
            msgs.error('Table format failure: No \'frametype\' column.')

        # Build the table
        nfiles = len(lines) - npaths - 1
        tbl = np.empty((nfiles, len(header)), dtype=object)

        for i in range(nfiles):
            row = np.array([ l.strip() for l in lines[i+npaths+1].split('|') ])[1:-1]
            if len(row) != tbl.shape[1]:
                raise ValueError('Data and header lines have mismatched columns!')
            tbl[i,:] = row
        data = {}
        for i,key in enumerate(header):
            data[key] = tbl[:,i]
        tbl = Table(data)

        ## Build full paths to file and set frame types
        #frametype = {}
        #data_files = []
        #for i in range(nfiles):
        #    frametype[tbl['filename'][i]] = tbl['frametype'][i]
        #    for p in paths:
        #        filename = os.path.join(p, tbl['filename'][i])
        #        if os.path.isfile(filename):
        #            break
        #    data_files.append(filename)
        #    #if not os.path.isfile(filename) and file_check:
        #    #    msgs.error('File does not exist: {0}'.format(filename))

        return tbl



def _find_block(lines, block):
    """
    Find a specific block of lines

    These must be wrapped within "block read" and "block end", e.g.

    setup read
    Setup A: 
      blah blah
    setup end

    Args:
        lines (:obj:`list`):
            List of file lines
        block (:obj:`str`):
            Name of block to parse

    Returns:
        int, int: Starting,ending line of the block;  
            -1 if not present

    """
    start = -1
    end = -1
    for i, l in enumerate(lines):
        entries = l.split()
        if start < 0 and entries[0] == block and entries[1] == 'read':
            start = i+1
            continue
        if entries[0] == block and entries[1] == 'end':
            end = i
            continue
        if start >= 0 and end >= 0:
            break
    return start, end


#    def _parse_data_file_name(inp, current_path=None, sort=True):
#        """
#        Expand the data file name as necessary and
#        then search for all data files.
#
#        Args:
#            inp (:obj:`str`):
#                Primary search string used by `glob.glob`_ to find the data files.
#            current_path (:obj:`str`, optional):
#                If provided, this is used as the root path for ``inp`` search
#                string.
#            sort (:obj:`bool`, optional):
#                Sort the list of data file names.
#
#        Returns:
#            :obj:`list`: Glob list of files in the generated a path
#        """
#        out = os.path.expanduser(inp) if inp[0] == '~' else inp
#        if current_path is not None:
#            out = os.path.join(current_path, out)
#        files = glob.glob(out)
#        return sorted(files) if sort else files
#        
#
#    def _read_data_file_names(lines, file_check=True, sort=True):
#        """
#        Read the raw data file format
#
#        Args:
#            lines (:obj:`list`):
#                The list of string lines to parse.
#            file_check (:obj:`bool`, optional):
#                Check that the parsed file names are valid files on disk.
#            sort (:obj:`bool`, optional):
#                Sort the list of file names on output
#
#        Returns:
#            :obj:`list`: List of data file names
#        """
#        # Pass through all the lines and:
#        #   - Determine if a path is set
#        #   - Gather the files to skip, can include wildcards
#        #   - Gather the files to read, can include wildcards
#        current_path = None
#        skip_inp = []
#        read_inp = []
#        for l in lines:
#            _l = l.split(' ')
#
#            if _l[0] == 'skip':
#                space_ind = l.index(" ")
#                path = l[space_ind + 1:]
#                skip_inp += _parse_data_file_name(path, current_path=current_path, sort=sort)
#                continue
#
#            if _l[0] == 'path':
#                space_ind = l.index(" ")
#                current_path = l[space_ind + 1:]
#                continue
#
#            read_inp += _parse_data_file_name(l, current_path=current_path, sort=sort)
#
#        # Remove any repeated lines
#        if len(skip_inp) > 0 and len(skip_inp) != len(set(skip_inp)):
#            msgs.warn('There are duplicated files to skip.')
#            skip_inp = list(set(skip_inp))
#        if len(read_inp) > 0 and len(read_inp) != len(set(read_inp)):
#            msgs.warn('There are duplicated files to read.')
#            read_inp = list(set(read_inp))
#
#        # Remove any files to skip
#        for _skip in skip_inp:
#            if _skip in read_inp:
#                read_inp.remove(_skip)
#
#        # Check that the files exist
#        if file_check:
#            for f in read_inp:
#                if not os.path.isfile(f):
#                    raise FileNotFoundError('{0} does not exist!'.format(f))
#
#        return read_inp

#    def _determine_data_format(lines):
#        """
#        Determine the format of the data block in the .pypeit file.
#
#        The test used in this function is pretty basic.  A table format is
#        assumed if the first character in *any* line is `|`.
#
#        Args:
#            lines (:obj:`list`):
#                The list of lines read from the data block of the pypeit
#                file.
#        
#        Returns:
#            str: The syntax of the data files to read::
#
#                'raw': A (list of) file roots to be read or found using
#                `glob`.
#
#                'table': ASCII output of an astropy.table.Table
#
#        """
#        for l in lines:
#            if l[0] == '|':
#                return 'table'
#        return 'raw'

def _read_pypeit_file_lines(ifile):
    """
    General parser for a pypeit file.

    - Checks that the file exists.
    - Reads all the lines in the file
    - Removes comments, empty lines, and replaces special characters.
    
    Applies to settings, setup, and user-level reduction files.

    Args:
        ifile (str): Name of the file to parse.

    Returns:
        :obj:`numpy.ndarray`: Returns a list of the valid lines in the
        files.
    """
    # Check the files
    if not os.path.isfile(ifile):
        msgs.error('The filename does not exist -' + msgs.newline() + ifile)

    # Read the input lines and replace special characters
    with open(ifile, 'r') as f:
        lines = np.array([l.replace('\t', ' ').replace('\n', ' ').strip() \
                                for l in f.readlines()])
    # Remove empty or fully commented lines
    lines = lines[np.array([ len(l) > 0 and l[0] != '#' for l in lines ])]
    # Remove appended comments and return
    return np.array([ l.split('#')[0] for l in lines ])
