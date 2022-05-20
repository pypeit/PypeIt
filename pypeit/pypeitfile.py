""" Class for I/O of PypeIt file"""
import os
import numpy as np
import yaml
import time
import io

import configobj

from astropy.table import Table

from pypeit import utils
from pypeit import msgs, __version__

from IPython import embed

class PypeItFile:
    """
    Class to load, process, and write PypeIt files

    Args:
        config (:obj:`dict`):
            Configuration dict
        data_table (:class:`astropy.table.Table`):
            Data block
        setup (:obj:`dict`, optional):
            dict defining the Setup
            Defaults to A and otherwise a null description
    """
    def __init__(self, config:dict, file_paths:list,
                 data_table:Table,
                 setup:dict=None):
        # Load up
        self.data = data_table
        self.file_paths = file_paths
        self.setup = setup

        # Load up config
        self.config = config

        # Vet
        self.vet()
        
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
        paths, usrtbl = cls._read_data_file_table(lines[s:e])
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
        msgs.info('PypeIt file loaded successfully.')

        # Instantiate
        slf = cls(confdict, paths, usrtbl, setup=sdict)

        return slf

    def vet(self):
        """ Check for required bits and pieces of the PypeIt file
        besides the input objects themselves
        """
        # Data table
        for key in ['filename', 'frametype']:
            if key not in self.data.keys():
                msgs.error("Add {:s} to your PypeIt file before using run_pypeit".format(key))

        # Config
        # Confirm we can generate a ConfiObj
        confObj = configobj.ConfigObj(self.config)
        # Confirm spectrograph is present
        if 'rdx' not in self.config.keys() or 'spectrograph' not in self.config['rdx'].keys():
            msgs.error(f"Missing spectrograph in the Parameter block of your PypeIt file.  Add it!")

        # Done
        msgs.info('PypeIt file vetted.')
        

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
            Table:  A Table with the data provided in 
            the pypeit file.  
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
        #if 'filename' not in header:
        #    msgs.error('Table format failure: No \'filename\' column.')
        #if 'frametype' not in header:
        #    msgs.error('Table format failure: No \'frametype\' column.')

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

        return paths, tbl
    def write(self, pypeit_file):
        """
        Generate a PypeIt file

        Args:
            pypeit_file (str): Name of PypeIt file to be generated
        """

        # Here we go
        with open(pypeit_file, 'w') as f:
            f.write('# Auto-generated PypeIt file using PypeIt version: {}\n'.format(__version__))
            #f.write('# {0}\n'.format(time.strftime("%a %d %b %Y %H:%M:%S",time.localtime())))
            f.write('# {0}\n'.format(time.strftime("%Y-%m-%d",time.localtime())))
            f.write("\n")

            # Parameter block
            f.write("# User-defined execution parameters\n")
            confObj = configobj.ConfigObj(self.config)
            f.write('\n'.join(confObj.write()))
            f.write('\n')
            f.write('\n')

            # Setup block
            if self.setup is not None:
                setup_lines = yaml.dump(utils.yamlify(
                    self.setup)).split('\n')[:-1]
            else: # Default
                setup_lines = ['Setup A:']

            f.write("# Setup\n")
            f.write("setup read\n")
            f.write('\n'.join(setup_lines)+'\n')
            f.write("setup end\n")
            f.write("\n")
            
            # Data block
            f.write("# Read in the data\n")
            f.write("data read\n")
            # paths and Setupfiles
            for path in self.file_paths:
                f.write(' path '+path+'\n')
            with io.StringIO() as ff:
                self.data.write(ff, format='ascii.fixed_width')
                data_lines = ff.getvalue().split('\n')[:-1]
            f.write('\n'.join(data_lines))
            f.write('\n')
            f.write("data end\n")
            f.write("\n")

        msgs.info('PypeIt file written to: {0}'.format(pypeit_file))



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


def _read_pypeit_file_lines(ifile):
    """
    General parser for a pypeit file.
    Used for many of our input files, including the main PypeIt file.

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
