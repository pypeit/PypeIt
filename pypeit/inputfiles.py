""" Class for I/O of PypeIt input files

.. include:: ../include/links.rst
"""

import os
import glob
import numpy as np
import yaml
from datetime import datetime
import io
import warnings

import configobj

from astropy.table import Table, column
from astropy.io import ascii

from pypeit import utils
from pypeit import msgs, __version__

from IPython import embed


class InputFile:
    """
    Generic class to load, process, and write PypeIt input files

    In practice, we use one of the children of this class,
    e.g. PypeItFile

    Args:
        config (:obj:`dict` or :obj:`list`):
            Configuration dict or list of config lines
            Converted to a ConfigObj
        file_paths (list):
            List of file paths for the data files
        data_table (`astropy.table.Table`_):
            Data block
        setup (:obj:`dict`):
            dict defining the Setup
            The first key contains the name
    """
    flavor = 'Generic' # Type of InputFile

    # Data block items
    required_columns = []
    data_block = None  # Defines naming of data block
    datablock_required = False # Denotes whether the data block is required
    
    setup_required = False # Denotes whether the setup block is required

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

    @staticmethod
    def readlines(ifile:str):
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
        
    @classmethod
    def from_file(cls, input_file:str): 
        """
        Parse the user-provided input file.

        Args:
            input_file (:obj:`str`):
                Name of input file

        Returns:
            :class:`InputFile`: An instance of the InputFile class
        """
        # Read in the pypeit reduction file
        msgs.info('Loading the reduction file')
        lines = cls.readlines(input_file)

        # Used to select the configuration lines: Anything that isn't part
        # of the data or setup blocks is assumed to be part of the
        # configuration
        is_config = np.ones(len(lines), dtype=bool)

        # Parse data block
        s, e = cls.find_block(lines, cls.data_block) 
        if s >= 0 and e < 0:
            msgs.error(
                f"Missing '{cls.data_block} end' in {input_file}")
        if s < 0 and e>0:
            msgs.error("You have not specified the start of the data block!")
        # Read it, if it exists
        if s>0 and e>0:
            paths, usrtbl = cls._read_data_file_table(lines[s:e])
            is_config[s-1:e+1] = False
        else:
            if cls.datablock_required:
                msgs.error("You have not specified the data block!")
            paths, usrtbl = [], None

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
        return cls(config=list(lines[is_config]), 
                  file_paths=paths, 
                  data_table=usrtbl, 
                  setup=sdict)

    def vet(self):
        """ Check for required bits and pieces of the Input file
        besides the input objects themselves
        """
        # Data table
        for key in self.required_columns:
            if key not in self.data.keys():
                msgs.error(f'Add {key} to the Data block of your {self.flavor} file before running.')


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
           If no configuration is available (:attr:`config` is 
           `None`), `None` is returned.

        Returns:
            :obj:`list`: List of the configuration lines prepared for
            writing to a file (and other usages).

        """
        return None if self.config is None else self.config.write()

    @property
    def filenames(self):
        """ List of path + filename's
        Wrapper to :func:`~pypeit.inputfiles.InputFile.path_and_files`.
        See that function for a full description.

        Returns:
            list or None: List of the full paths to each data file
            or None if `filename` is not part of the data table
            or there is no data table!
        """
        # Return
        return self.path_and_files('filename')

    @staticmethod
    def _parse_setup_lines(lines):
        """
        Return a list of the setup names and corresponding Setup dict

        Args:
            lines (`numpy.ndarray`_): Setup lines as an array

        Returns:
            tuple: list, dict.  List of setup name, setup dict

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
            # Strip allows for preceding spaces before path
            prs = l.strip().split(" ")
            if prs[0] != 'path':
                break
            paths += [ prs[1] ]

        npaths = len(paths)

        # Read the table
        tbl = ascii.read(lines[npaths:].tolist(), 
                         header_start=0, 
                         data_start=1, 
                         delimiter='|', 
                         format='basic')

        # Backwards compatability (i.e. the leading |)
        if list(tbl.keys())[0] == 'col0':
            message = 'Your PypeIt file has leading | characters in the data table, which is the old '\
                'format.  Please update your file to remove leading and trailing | characters '\
                'because their inclusion will be deprecated.'
            warnings.warn(message, DeprecationWarning)
            tbl.remove_column('col0')
            tbl.remove_column('_1')

        ## Recast each as "object" in case the user has mucked with the Table
        ##  e.g. a mix of floats and None
        ##  Also handle Masked columns -- fill with ''
        for key in tbl.keys():
            # Object
            tbl[key] = tbl[key].data.astype(object)
            if isinstance(tbl[key], column.MaskedColumn):
                # Fill with empty string
                tbl[key].fill_value = ''
                tbl[key] = tbl[key].filled()

        # Build the table -- Old code
        #  Because we allow (even encourage!) the users to modify entries by hand, 
        #   we have a custom way to parse what is largely a standard fixed_width table
        #nfiles = len(lines) - npaths - 1
        #header = [ l.strip() for l in lines[npaths].split('|') ][1:-1]
        #tbl = np.empty((nfiles, len(header)), dtype=object)
#
#        for i in range(nfiles):
#            row = np.array([ l.strip() for l in lines[i+npaths+1].split('|') ])[1:-1]
#            if len(row) != tbl.shape[1]:
#                raise ValueError('Data and header lines have mismatched columns!')
#            tbl[i,:] = row
#        data = {}
#        for i,key in enumerate(header):
#            data[key] = tbl[:,i]
#        tbl = Table(data)

        # Return
        return paths, tbl

    @staticmethod
    def find_block(lines, block):
        """
        Find a specific block of lines

        These must be wrapped within ``block read`` and ``block end``, e.g.::

            setup read
            Setup A: 
            ...
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
        the full path from the column of the data Table
        specified by `key`.  The files must exist and be 
        within one of the paths for this to succeed.

        Args:
            key (str): Column of self.data with the filenames of interest
            skip_blank (bool, optional): If True, ignore any
                entry that is '', 'none' or 'None'. Defaults to False.

        Returns:
            list: List of the full paths to each data file
            or None if `filename` is not part of the data table
            or there is no data table!

        Raises:
            PypeItError:
                Raised if the path+file does not exist

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
                    filename = os.path.join(p, row[key])
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

    def write(self, input_file):
        """
        Write an Input file to disk

        Args:
            input_file (str): Name of PypeIt file to be generated
        """

        # Here we go
        with open(input_file, 'w') as f:
            f.write(f'# Auto-generated {self.flavor} input file using PypeIt version: {__version__}\n')
            #f.write('# {0}\n'.format(time.strftime("%Y-%m-%d",time.localtime())))
            f.write('# UTC {0}\n'.format(datetime.utcnow().isoformat(timespec='milliseconds')))
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
            if self.data is not None:
                f.write("# Data block \n")
                f.write(f"{self.data_block} read\n")
                # paths and Setupfiles
                if self.file_paths is not None:
                    for path in self.file_paths:
                        f.write(' path '+path+'\n')
                with io.StringIO() as ff:
                    self.data.write(ff, format='ascii.fixed_width',
                                    bookend=False)
                    data_lines = ff.getvalue().split('\n')[:-1]
                f.write('\n'.join(data_lines))
                f.write('\n')
                f.write(f"{self.data_block} end\n")
                f.write("\n")

        msgs.info(f'{self.flavor} input file written to: {input_file}')



class PypeItFile(InputFile):
    """Child class for the PypeIt file
    """
    flavor = 'PypeIt'  # Defines naming of file

    required_columns = ['filename', 'frametype']
    data_block = 'data'  # Defines naming of data block
    datablock_required = True

    setup_required = True

    def vet(self):
        """ Check for required bits and pieces of the PypeIt file
        besides the input objects themselves
        """
        super().vet()

        # Confirm spectrograph is present
        if 'rdx' not in self.config.keys() or 'spectrograph' not in self.config['rdx'].keys():
            msgs.error(f"Missing spectrograph in the Parameter block of your PypeIt file.  Add it!")

        # Setup
        setup_keys = list(self.setup)
        if 'Setup' not in setup_keys[0]:
            msgs.error("Setup does not appear in your setup block! Add it")

        # Done
        msgs.info('PypeIt file successfully vetted.')

    @property
    def frametypes(self):
        """Return a dict of the frametypes
        with key, item the filename, frametype 

        Returns:
            dict: 
        """
        return {row['filename']:row['frametype'] for row in self.data}

class SensFile(InputFile):
    """Child class for the Sensitivity input file
    """
    data_block = 'sens'  # Defines naming of data block
    flavor = 'Sens'  # Defines naming of file
    datablock_required = False
    setup_required = False

class FluxFile(InputFile):
    """Child class for the Fluxing input file
    """
    data_block = 'flux'  # Defines naming of data block
    flavor = 'Flux'  # Defines naming of file
    setup_required = False
    datablock_required = True
    required_columns = ['filename']

    def vet(self):
        """ Check for required parts of the Fluxing input
        file and handle the various options for sensfile
        """
        super().vet()

        # Add a dummy sensfile column?
        #  This is allowed if using an archived sensitivity function
        #  And the checking has to be done in the script as the specgtrograph must be known..
        if 'sensfile' not in self.data.keys():
            msgs.warn("sensfile column not provided.  Fluxing will crash if an archived sensitivity function does not exist")
            self.data['sensfile'] = ''

    @property
    def sensfiles(self):
        """Generate a list of the sensitivity files with 
        the full path.  The files must exist and be 
        within one of the paths (or the current
        folder with not other paths specified) for this to succeed.

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

class Coadd1DFile(InputFile):
    """Child class for coaddition in 1D
    """
    data_block = 'coadd1d'  # Defines naming of data block
    flavor = 'Coadd1D'  # Defines naming of file
    setup_required = False
    datablock_required = True
    required_columns = ['filename', 'obj_id']

    @property
    def objids(self):
        # Generate list, scrubbing empty entries
        oids = [item for item in self.data['obj_id'] if item.strip() not in ['', 'none', 'None']]

        # Inflate as needed
        if len(oids) == 1 and len(oids) < len(self.data):
            oids = oids*len(self.data)
        # Return
        return oids


class Coadd2DFile(InputFile):
    """Child class for coaddition in 2D
    """
    data_block = 'spec2d'  # Defines naming of data block
    flavor = 'Coadd2D'  # Defines naming of file
    setup_required = False
    datablock_required = True
    required_columns = ['filename'] 

    def vet(self):
        """ Check for required bits and pieces of the .coadd2d file
        besides the input objects themselves
        """
        super().vet()

        # Confirm spectrograph is present
        if 'rdx' not in self.config.keys() or 'spectrograph' not in self.config['rdx'].keys():
            msgs.error(f"Missing spectrograph in the Parameter block of your .coadd2d file.  Add it!")

        # Done
        msgs.info('.coadd2d file successfully vetted.')


class Coadd3DFile(InputFile):
    """Child class for coadding spec2d files into datacubes
    """
    data_block = 'spec2d'  # Defines naming of data block
    flavor = 'Coadd3d'  # Defines naming of file
    setup_required = False
    datablock_required = True
    required_columns = ['filename'] 

    def vet(self):
        """ Check for required bits and pieces of the .coadd2d file
        besides the input objects themselves
        """
        super().vet()

        # Confirm spectrograph is present
        if 'rdx' not in self.config.keys() or 'spectrograph' not in self.config['rdx'].keys():
            msgs.error(f"Missing spectrograph in the Parameter block of your .coadd2d file.  Add it!")

        # Done
        msgs.info('.cube file successfully vetted.')

    @property
    def options(self):
        """
        Parse the options associated with a cube block.
        Here is a description of the available options:

        scale_corr     : The name of an alternative spec2d file that is used for the relative spectral scale correction.
                        This parameter can also be set for all frames with the default command::

                            [reduce]
                                [[cube]]
                                    scale_corr = spec2d_alternative.fits

        Returns
        -------
        opts: dict
            Dictionary containing cube options.
        """
        # Define the list of allowed parameters
        opts = dict(scale_corr=None)

        # Get the scale correction files
        scale_corr = self.path_and_files('scale_corr', skip_blank=True)
        if scale_corr is None:
            opts['scale_corr'] = [None]*len(self.filenames)
        elif len(scale_corr) == 1 and len(self.filenames) > 1:
            scale_corr = scale_corr*len(self.filenames)
        elif len(scale_corr) != 0:
            opts['scale_corr'] = scale_corr
        # Return all options
        return opts


class TelluricFile(InputFile):
    """Child class for telluric corrections
    """
    data_block = None  # Defines naming of data block
    flavor = 'Telluric'  # Defines naming of file
    setup_required = False
    datablock_required = False


class FlexureFile(InputFile):
    """Child class for flexure corrections
    """
    data_block = 'flexure'  # Defines naming of data block
    flavor = 'Flexure'  # Defines naming of file
    setup_required = False
    datablock_required = True
    required_columns = ['filename']

    def vet(self):
        """ Check for required bits and pieces of the .flex file
        besides the input objects themselves
        """
        super().vet()

        # Confirm spectrograph is present
        if 'rdx' not in self.config.keys() or 'spectrograph' not in self.config['rdx'].keys():
            msgs.error(f"Missing spectrograph in the Parameter block of your .flex file.  Add it!")

        # Done
        msgs.info('.flex file successfully vetted.')

class Collate1DFile(InputFile):
    """Child class for collate 1D script
    """
    data_block = 'spec1d'  # Defines naming of data block
    flavor = 'Collate1D'  # Defines naming of file
    setup_required = False
    datablock_required = True
    required_columns = ['filename']

    @property
    def filenames(self):
        """ List of path + filename's

        Allows for wildcads

        Returns:
            list or None: List of the full paths to each data file
            or None if `filename` is not part of the data table
            or there is no data table!
        """
        all_files = []
        paths = [os.getcwd()] if len(self.file_paths) == 0 else self.file_paths
        # Paths?
        for p in paths:
            for row in self.data['filename']:
                all_files += glob.glob(os.path.join(p, row))

        # Return
        return all_files