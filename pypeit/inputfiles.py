""" Class for I/O of PypeIt input files"""

import os
import numpy as np
import yaml
import time
import io
import warnings

import configobj

from astropy.table import Table
from astropy.io import ascii

from pypeit import utils
from pypeit import msgs, __version__

from IPython import embed

def read_pypeit_file_lines(ifile:str):
    """
    General parser for a PypeIt input file.
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

class InputFile:
    """
    Generic class to load, process, and write PypeIt input files

    In practice, we use one of the children of this clases,
    e.g. PypeItFile

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
    flavor = 'Generic'

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
        s, e = cls.find_block(lines, cls.data_block) #'data')
        if s >= 0 and e < 0:
            msgs.error(
                f"Missing '{cls.data_block} end' in {input_file}")
        if s < 0:
            msgs.error("You haven't specified any data!")
        paths, usrtbl = cls._read_data_file_table(lines[s:e])
        is_config[s-1:e+1] = False

        # Parse the setup block
        setup_found = False
        s, e = cls.find_block(lines, 'setup')
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
        msgs.info(f'{cls.flavor} input file loaded successfully.')

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
            or None if self.config is None
        """
        if self.config is None:
            return None
        else:
            return self.config.write()

    @property
    def filenames(self):
        """ List of path + filename's
        Wrapper to path_and_files().
        See that function for a full description.

        Returns:
            list: List of the full paths to each data file
            or None if `filename` is not part of the data table
            or there is no data table!
        """
        # Return
        return self.path_and_files('filename')

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

        Because we allow (even encourage!) the users to modify entries by hand, 
        we have a custom way to parse what is largely a standard fixed_width 
        ASCII table
        
        Args:
            lines (:obj:`list`):
                List of lines *within the data* block read from the input file.
        
        Returns:
            tuple: list, Table.  A list of the paths provided (can be empty)
            and a Table with the data provided in the input file.  
        """

        # Allow for multiple paths
        paths = []
        for l in lines:
            space_ind = l.index(" ")
            if l[:space_ind].strip() != 'path':
                break
            paths += [ l[space_ind+1:] ]

        npaths = len(paths)

        # Read the table
        #embed(header='263')
        #tbl = ascii.read(lines[npaths:].tolist(), 
        #                 header_start=0, 
        #                 data_start=1, 
        #                 format='fixed_width')
        ## Recast each as "object" in case the user has mucked with the Table
        ##  e.g. a mix of floats and None
        #for key in tbl.keys():
        #    tbl[key] = tbl[key].data.astype(object)
        ##embed(header='249 of inputfiles')

        # Build the table
        #  Because we allow (even encourage!) the users to modify entries by hand, 
        #   we have a custom way to parse what is largely a standard fixed_width table
        nfiles = len(lines) - npaths - 1
        header = [ l.strip() for l in lines[npaths].split('|') ][1:-1]
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

        # Return
        return paths, tbl

    @staticmethod
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

    def path_and_files(self, key:str, skip_blank=False):
        """Generate a list of the filenames with 
        the full path.  The files must exist and be 
        within one of the paths for this to succeed.

        An error is raised if the path+file does not exist

        Args:
            key (str): Column of self.data with the filenames of interest
            skip_blank (bool, optional): If True, ignore any
            entry that is '', 'none' or 'None'. Defaults to False.

        Returns:
            list: List of the full paths to each data file
            or None if `filename` is not part of the data table
            or there is no data table!

        """
        if self.data is None or key not in self.data.keys():
            return None

        ## Build full paths to file and set frame types
        data_files = []
        for row in self.data:

            # Skip Empty entries?
            if skip_blank and row[key].strip() in ['', 'none', 'None']:
                continue

            # Searching..
            if len(self.file_paths) > 0:
                for p in self.file_paths:
                    filename = os.path.join(
                        p, row[key])
                    if os.path.isfile(filename):
                        break
            else:
                filename = row[key]

            # Check we got a good hit
            if not os.path.isfile(filename): 
                msgs.error(f"{row[key]} does not exist in one of the provided paths.  Modify your input {self.flavor} file")
            data_files.append(filename)

        # Return
        return data_files

    def write(self, pypeit_input_file):
        """
        Write a PypeIt input file to disk

        Args:
            pypeit_input_file (str): Name of PypeIt file to be generated
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



class PypeItFile(InputFile):
    """Child class for the PypeIt file
    """
    data_block = 'data'  # Defines naming of data block
    flavor = 'PypeIt'  # Defines naming of file
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

class FluxFile(InputFile):
    """Child class for the Fluxing input file
    """
    data_block = 'flux'  # Defines naming of data block
    flavor = 'Flux'  # Defines naming of file
    setup_required = False

    @property
    def sensfiles(self):
        """Generate a list of the sensitivity files with 
        the full path.  The files must exist and be 
        within one of the paths (or the current
        folder with not other paths specified) for this to succeed.

        Raises:
            ValueError: _description_

        Returns:
            list: List of full path to each data file
            or None if `filename` is not part of the data table
            or there is no data table!
        """
        # Grab em
        sens_files = self.path_and_files('sensfile', skip_blank=True)
        # Pad out
        if len(sens_files) == 1 and len(self.filenames) > 1:
            sens_files = sens_files*len(self.filenames)
            
        # Return
        return sens_files

    @classmethod
    def read_old_fluxfile(cls, ifile):
        """
        Read an old style ``PypeIt`` flux file, akin to a standard ``PypeIt`` file.

        The top is a config block that sets ParSet parameters.  

        Args:
            ifile (:obj:`str`):
                Name of the flux file

        Returns:
            pypeit.inputfiles.FluxFile:
        """
        # Warn
        warnings.warn("The old file type is deprecated and this code may disappear", DeprecationWarning)

        # Read in the pypeit reduction file
        msgs.info('Loading the fluxcalib file')
        lines = read_pypeit_file_lines(ifile)
        is_config = np.ones(len(lines), dtype=bool)

        # Parse the fluxing block
        spec1dfiles = []
        sensfiles_in = []
        s, e = InputFile.find_block(lines, 'flux')
        if s >= 0 and e < 0:
            msgs.error("Missing 'flux end' in {0}".format(ifile))
        elif (s < 0) or (s==e):
            msgs.error("Missing flux block in {0}. Check the input format for the .flux file".format(ifile))
        else:
            for ctr, line in enumerate(lines[s:e]):
                prs = line.split(' ')
                spec1dfiles.append(prs[0])
                if len(prs) > 1:
                    sensfiles_in.append(prs[1])
            is_config[s-1:e+1] = False

        # data table
        data = Table()
        data['filename'] = spec1dfiles
        # Sensfiles are a fussy matter
        if len(sensfiles_in) == 0:
            sensfiles = ['']*len(spec1dfiles)
        elif len(sensfiles_in) == 1:
            sensfiles = sensfiles_in*len(spec1dfiles)
        else:
            sensfiles = sensfiles_in
        data['sensfile'] = sensfiles

        # Instantiate
        slf = cls(config=list(lines[is_config]), 
                  data_table=data)

        return slf