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

class InputFile:
    """
    Generic class to load, process, and write PypeIt input files

    In practice, we use one of the children of this class

    Args:
        config (:obj:`dict` or :obj:`list`):
            Configuration dict or list of config lines
            Converted to a ConfigObj
        file_paths (list):
            List of file paths for the data files
        data_table (:class:`astropy.table.Table`):
            Data block
        setup (:obj:`dict`):
            dict defining the Setup
            The first key contains the name
    """
    data_block = None  # Defines naming of data block
    setup_required = False
    flavor = 'generic'

    def __init__(self, 
                 config=None, 
                 file_paths:list=None,
                 data_table:Table=None,
                 setup:dict=None):
        # Load up
        self.data = data_table
        self.file_paths = file_paths
        self.setup = setup

        # Load up ConfigObj
        if config is not None:
            self.config = configobj.ConfigObj(config)
        else:
            self.config = None

        # Vet
        self.vet()
        
    @classmethod
    def from_file(cls, input_file:str): 
        """
        Parse the user-provided input file.

        Args:
            input_file (:obj:`str`):
                Name of input file

        Returns:
            pypeItFile (:class:`PypeItFile`):
        """
        # Read in the pypeit reduction file
        msgs.info('Loading the reduction file')
        lines = read_pypeit_file_lines(input_file)

        # Used to select the configuration lines: Anything that isn't part
        # of the data or setup blocks is assumed to be part of the
        # configuration
        is_config = np.ones(len(lines), dtype=bool)

        # Parse data block
        s, e = find_block(lines, cls.data_block) #'data')
        if s >= 0 and e < 0:
            msgs.error(
                f"Missing '{cls.data_block} end' in {input_file}")
        if s < 0:
            msgs.error("You haven't specified any data!")
        paths, usrtbl = cls._read_data_file_table(lines[s:e])
        is_config[s-1:e+1] = False

        # Parse the setup block
        setup_found = False
        s, e = find_block(lines, 'setup')
        if s >= 0 and e < 0 and cls.setup_required:
            msgs.error(f"Missing 'setup end' in {input_file}")
        elif s < 0 and cls.setup_required:
            msgs.error(f"Missing 'setup read' in {input_file}")
        elif s >= 0 and e > 0:
            setup_found = True

        # Proceed
        if setup_found:
            setups, sdict = cls._parse_setup_lines(lines[s:e])
            is_config[s-1:e+1] = False
        else:
            sdict = None

        # vet
        msgs.info('PypeIt input file loaded successfully.')

        # Instantiate
        slf = cls(config=list(lines[is_config]), 
                  file_paths=paths, 
                  data_table=usrtbl, 
                  setup=sdict)

        return slf

    def vet(self):
        """ Check for required bits and pieces of the PypeIt file
        besides the input objects themselves
        """
        pass


    @property
    def setup_name(self):
        """Return the setup name

        Returns:
            str: Setup name.  Usually a single, capitalized Letter
        """
        keys = list(self.setup.keys())
        return keys[0].split(' ')[-1]

    @property
    def cfg_lines(self):
        """Return a list containing the configuration lines

        Returns:
            list: List of the configuration lines held in self.config
        """
        return self.config.write()

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

#    @staticmethod
#    def _read_data_file_table(lines):
#        """
#        Read the file table format.
#        
#        Args:
#            lines (:obj:`list`):
#                List of lines *within the data* block read from the pypeit
#                file.
#        
#        Returns:
#            Table:  A Table with the data provided in 
#            the pypeit file.  
#        """
#        return None, None

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


        return paths, tbl

    def write(self, pypeit_input_file):
        """
        Write a PypeIt input file to disk

        Args:
            pypeit_file (str): Name of PypeIt file to be generated
        """

        # Here we go
        with open(pypeit_input_file, 'w') as f:
            f.write(f'# Auto-generated {self.flavor} input file using PypeIt version: {__version__}\n')
            #f.write('# {0}\n'.format(time.strftime("%a %d %b %Y %H:%M:%S",time.localtime())))
            f.write('# {0}\n'.format(time.strftime("%Y-%m-%d",time.localtime())))
            f.write("\n")

            # Parameter block
            if self.config is not None:
                f.write("# User-defined execution parameters\n")
                f.write('\n'.join(self.cfg_lines))
                f.write('\n')
                f.write('\n')

            # Setup block
            if self.setup is not None:
                setup_lines = yaml.dump(utils.yamlify(
                    self.setup)).split('\n')[:-1]
            elif self.setup_required: # Default
                setup_lines = ['Setup A:']
            else:
                setup_lines = None

            if setup_lines is not None:
                f.write("# Setup\n")
                f.write("setup read\n")
                f.write('\n'.join(setup_lines)+'\n')
                f.write("setup end\n")
                f.write("\n")
            
            # Data block
            f.write("# File block \n")
            f.write(f"{self.data_block} read\n")
            # paths and Setupfiles
            if self.file_paths is not None:
                for path in self.file_paths:
                    f.write(' path '+path+'\n')
            with io.StringIO() as ff:
                self.data.write(ff, format='ascii.fixed_width')
                data_lines = ff.getvalue().split('\n')[:-1]
            f.write('\n'.join(data_lines))
            f.write('\n')
            f.write(f"{self.data_block} end\n")
            f.write("\n")

        msgs.info(f'{self.flavor} input file written to: {pypeit_input_file}')


def find_block(lines, block):
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


def read_pypeit_file_lines(ifile):
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

class PypeItFile(InputFile):
    data_block = 'data'  # Defines naming of data block
    flavor = 'PypeIt'  # Defines naming of data block
    setup_required = True

    def vet(self):
        """ Check for required bits and pieces of the PypeIt file
        besides the input objects themselves
        """

        # Data table
        for key in ['filename', 'frametype']:
            if key not in self.data.keys():
                msgs.error("Add {:s} to your PypeIt file before using run_pypeit".format(key))

        # Confirm spectrograph is present
        if 'rdx' not in self.config.keys() or 'spectrograph' not in self.config['rdx'].keys():
            msgs.error(f"Missing spectrograph in the Parameter block of your PypeIt file.  Add it!")

        # Setup
        setup_keys = list(self.setup)
        assert 'Setup' in setup_keys[0]

        # Done
        msgs.info('PypeIt file successfully vetted.')

        
    @property
    def data_files(self):
        """Generate a list of the data files with 
        the full path.  The files must exist and be 
        within one of the paths for this to succeed.

        Raises:
            ValueError: _description_

        Returns:
            list: List of full path to each data file
        """
        ## Build full paths to file and set frame types
        data_files = []
        for row in self.data:
            #frametype[tbl['filename'][i]] = tbl['frametype'][i]
            for p in self.file_paths:
                filename = os.path.join(
                    p, row['filename'])
                if os.path.isfile(filename):
                    break
            # Check we got a good hit
            if not os.path.isfile(filename): 
                msgs.error(f"{row['filename']} does not exist in one of the provided paths.  Remove from your PypeIt file")
            data_files.append(filename)
        # Return
        return data_files

    @property
    def frametypes(self):
        """Return a dict of the frametypes
        with key, item the filename, frametype 

        Returns:
            dict: 
        """
        frametypes = {}
        for row in self.data:
            frametypes[row['filename']] = row['frametype']
        #
        return frametypes

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


        return paths, tbl

class FluxFile(InputFile):
    data_block = 'flux'  # Defines naming of data block
    flavor = 'Flux'  # Defines naming of data block
    setup_required = False

    @property
    def data_files(self):
        """Generate a list of the data files with 
        the full path.  The files must exist and be 
        within one of the paths for this to succeed.

        Raises:
            ValueError: _description_

        Returns:
            list: List of full path to each data file
        """
        return list(self.data['filename'].data)

'''
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
            the input file.  
        """

        # Allow for multiple paths
        paths = None

        spec1dfiles = []
        sensfiles_in = []
        for ctr, line in enumerate(lines):
            prs = line.split(' ')
            spec1dfiles.append(prs[0])
            if len(prs) > 1:
                sensfiles_in.append(prs[1])
        # 
        if len(sensfiles_in) == 1:
            sensfiles_in = sensfiles_in*len(spec1dfiles)
        elif len(sensfiles_in) == 0:
            # This allows for auto-sensitivity files (e.g. DEIMOS)
            sensfiles_in = None

        # Build the Table
        tbl = Table()
        tbl['filename'] = spec1dfiles
        tbl['sensfile'] = sensfiles_in

        return paths, tbl
'''