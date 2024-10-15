""" Class for I/O of PypeIt input files

.. include:: ../include/links.rst
"""
from pathlib import Path
import os
import glob
import numpy as np
import yaml
from datetime import datetime
import io
import warnings
from collections.abc import Sequence
import configobj

# TODO: datetime.UTC is not defined in python 3.10.  Remove this when we decide
# to no longer support it.
try:
    __UTC__ = datetime.UTC
except AttributeError as e:
    from datetime import timezone
    __UTC__ = timezone.utc

from astropy.table import Table, column
from astropy.io import ascii

from pypeit import utils
from pypeit.io import files_from_extension
from pypeit import msgs, __version__
from pypeit.spectrographs.util import load_spectrograph
from pypeit.par.pypeitpar import PypeItPar

from IPython import embed


class InputFile:
    """
    Generic class to load, process, and write PypeIt input files

    In practice, we use one of the children of this class,
    e.g. PypeItFile

    This class has limited support for preserving comments in the input files.
    We currently support comments in the configuration section, and commented
    out data lines in the data sections. Other comments are not preserved.

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
        vet (bool,Optional):
            Whether or not to vet the file after iniitialization. Defaults to True.
            The vet() method can be called after initialization if needed.
        preserve_comments (bool, Optional):
            Whether or not to preserve comments in the input file
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
                 setup:dict=None,
                 vet:bool=True,
                 preserve_comments:bool=False):
        # Load up
        self.data = data_table
        self.file_paths = file_paths
        self.setup = setup
        self.preserve_comments=preserve_comments
        self.config = None if config is None else configobj.ConfigObj(config)

        # Vet
        if vet:
            self.vet()

    @staticmethod
    def remove_comments_and_blanks(lines : np.ndarray) -> np.ndarray:
        # Remove comment and blank lines
        lines = lines[np.array([ len(l.strip()) > 0 and l.strip()[0] != '#' for l in lines ])]

        # Remove appended comments and return
        return np.array([ l.split('#')[0].strip() for l in lines ])
        
    @staticmethod
    def readlines(ifile:str):
        """
        General parser for a PypeIt input file.
        Used for many of our input files, including the main PypeIt file.

        - Checks that the file exists.
        - Reads all the lines in the file
        - Replaces special characters.
        
        Applies to settings, setup, and user-level reduction files.

        Args:
            ifile (str): Name of the file to parse.

        Returns:
            `numpy.ndarray`_: Returns a list of the valid lines in the files.
        """
        # Check the files
        if not os.path.isfile(ifile):
            msgs.error('The filename does not exist -' + msgs.newline() + ifile)

        # Read the input lines and replace special characters
        with open(ifile, 'r') as f:
            return np.array([l.replace('\t', ' ').rstrip()  for l in f.readlines()])
        
    @classmethod
    def from_file(cls, input_file:str, vet:bool=True, preserve_comments:bool=False):
        """
        Parse the user-provided input file.

        Args:
            input_file (:obj:`str`):
                Name of input file
            vet (bool,Optional):
                Whether or not to vet the file after iniitialization. Defaults to True.
                The vet() method can be called after initialization if needed.
            preserve_comments (bool,Optional):
                Whether or not to preserve comments and blank lines 
                in the congiguration section and the data table. Defaults to False.

        Returns:
            :class:`InputFile`: An instance of the InputFile class
        """
        # Read in the pypeit reduction file
        msgs.info('Loading the reduction file')
        lines = cls.readlines(input_file)
       
        if not preserve_comments:
            lines = InputFile.remove_comments_and_blanks(lines)

        # Used to select the configuration lines: Anything that isn't part
        # of the data or setup blocks is assumed to be part of the
        # configuration
        is_config = np.ones(len(lines), dtype=bool)

        # Parse data block
        data_start, data_end = cls.find_block(lines, cls.data_block) 
        if data_start >= 0 and data_end < 0:
            msgs.error(
                f"Missing '{cls.data_block} end' in {input_file}")
        if data_start < 0 and data_end>0:
            msgs.error("You have not specified the start of the data block!")
        # Read it, if it exists
        if data_start>0 and data_end>0:
            paths, usrtbl = cls._read_data_file_table(lines[data_start:data_end], preserve_comments)
            is_config[data_start-1:data_end+1] = False
            data_block_found = True
        else:
            if cls.datablock_required:
                msgs.error("You have not specified the data block!")            
            paths, usrtbl = [], None
            data_block_found = False

        # Parse the setup block
        setup_found = False
        setup_start, setup_end = cls.find_block(lines, 'setup')
        if setup_start >= 0 and setup_end < 0 and cls.setup_required:
            msgs.error(f"Missing 'setup end' in {input_file}")
        elif setup_start < 0 and cls.setup_required:
            msgs.error(f"Missing 'setup read' in {input_file}")
        elif setup_start >= 0 and setup_end > 0:
            setup_found = True

        # Proceed
        if setup_found:
            setup_lines = lines[setup_start:setup_end]

            # We don't currently support preserving comments in the
            # setup block as doing so causess parsing problems and
            # backwards compatibility problems with some of the
            # files in the dev-suite. So we remove those 
            # from the setup_lines if they weren't removed before
            if preserve_comments:
                setup_lines = InputFile.remove_comments_and_blanks(setup_lines)

            setups, sdict = cls._parse_setup_lines(setup_lines)
            is_config[setup_start-1:setup_end+1] = False
        else:
            sdict = None

        # Lines between setup and data blocks are currently considered
        # part of the configuration. As are lines after the data block.
        # If we're preserving comments, ConfigObj may pick up those comments,
        # and then write them to a different place when resaving the file.
        # So we ignore those lines.
        if preserve_comments:
            if data_block_found and setup_found:
                # Ignore lines between datablock and setup block.
                if data_end+1 < setup_start-1:
                    # Data block is before setup block.
                    # This never seems to happen but it's technically supported
                    is_config[data_end+1:setup_start-1] = False
                elif setup_end+1 < data_start+1:
                    # Setup block is before data block
                    # Clear the lines between them and after the datablock
                    is_config[setup_end+1:data_start-1] = False
                    is_config[data_end+1:] = False
                # Else the two blocks are adjacent and there's no lines to preserve
            elif data_block_found:
                is_config[data_end+1:] = False

        # vet
        msgs.info(f'{cls.flavor} input file loaded successfully.')

        # Instantiate
        return cls(config=list(lines[is_config]), 
                  file_paths=paths, 
                  data_table=usrtbl, 
                  setup=sdict,
                  vet=vet,
                  preserve_comments=preserve_comments)

    def vet(self):
        """ Check for required bits and pieces of the Input file
        besides the input objects themselves
        """
        # Data table
        if self.data is None:
            if self.datablock_required:
                msgs.error("You have not specified the data block!")
        else:
            for key in self.required_columns:
                if key not in self.data.keys():
                    msgs.error(f'Add {key} to the Data block of your {self.flavor} file before running.')

        if self.setup_required and self.setup is None:
            msgs.error("Add setup info to your PypeIt file in the setup block!")

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
        return self.path_and_files('filename',include_commented_out=self.preserve_comments)

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

        return setups, sdict

    @staticmethod
    def _read_data_file_table(lines, preserve_comments):
        """
        Read the file table format.

        Because we allow (even encourage!) the users to modify entries by hand, 
        we have a custom way to parse what is largely a standard fixed_width 
        ASCII table
        
        Args:
            lines (:obj:`list`):
                List of lines *within the data* block read from the input file.
            preserve_comments (bool):
                Whether or not to preserve comments in the input file.

        Returns:
            tuple: A :obj:`list` with the paths provided (can be empty) and an
            `astropy.table.Table`_ with the data provided in the input file.  
        """
        # Allow for multiple paths
        paths = []
        for i, l in enumerate(lines):
            # Ignore comments, blank lines in the paths section 
            # of the data block
            #  strip allows for preceding spaces before path
            
            l = l.strip().split("#")[0]
            
            if len(l) == 0:
                # Skip empty line
                continue

            # Split at only the first space to allow paths that contain spaces
            prs = [l[:l.find(' ')], l[l.find(' ')+1:]]
            if prs[0] != 'path':
                break
            paths += [ prs[1] ]

        # We currently preserve commented out lines of the table, but not other comments in the
        # table. To prevent them from breaking ascii.read when preserve_comments is true
        # we remove lines without the correct number of delimiters
        delimiter_num = lines[i].count("|")

        table_lines = [l for l in lines[i:] if l.count("|") == delimiter_num]

        # Read the table
        tbl = ascii.read(table_lines, 
                         header_start=0, 
                         data_start=1, 
                         delimiter='|', 
                         comment=None if preserve_comments else "#",
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
        for key in tbl.keys():
            # Object
            tbl[key] = tbl[key].data.astype(object)
            # RJC -- Not sure why we need to fill masked columns with empty data. Let's just retain the masked values.
            # if isinstance(tbl[key], column.MaskedColumn):
                # Fill with empty string
                # tbl[key].fill_value = np.ma.masked
                # tbl[key] = tbl[key].filled()

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
            # Ignore comments/empty lines/leading whitespace
            line = l.split("#")[0].strip()
            if len(line) == 0:
                continue
            entries = line.split()
            if start < 0 and entries[0] == block and entries[1] == 'read':
                start = i+1
                continue
            if entries[0] == block and entries[1] == 'end':
                end = i
                continue
            if start >= 0 and end >= 0:
                break
        return start, end

    def path_and_files(self, key:str, skip_blank=False, include_commented_out=False, check_exists=True):
        """Generate a list of the filenames with 
        the full path from the column of the data `astropy.table.Table`_
        specified by `key`.  The files must exist and be 
        within one of the paths for this to succeed.

        Args:
            key (str): Column of self.data with the filenames of interest
            skip_blank (bool, optional): If True, ignore any
                entry that is '', 'none' or 'None'. Defaults to False.
            check_exists (bool, optional):If False, PypeIt will not
                check if 'key' exists as a file. Defaults to True.
            include_commented_out (bool,Optional): If False, commented out files will not be included. If True, they
                are included, without the "#" character.

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
            rowkey = '' if np.ma.is_masked(row[key]) else row[key]

            # Skip Empty entries?
            if skip_blank and rowkey.strip() in ['', 'none', 'None']:
                continue

            # Skip commented out entries
            if rowkey.strip().startswith("#"):
                if not include_commented_out:
                    continue
                # Strip the comment character and any whitespace following it
                # from the filename
                name = rowkey.strip("# ")
            else:
                name = rowkey

            # Searching..
            if len(self.file_paths) > 0:
                for p in self.file_paths:
                    filename = os.path.join(p,name)
                    if os.path.isfile(filename):
                        break
            else:
                filename = rowkey

            # Check we got a good hit
            if check_exists and not os.path.isfile(filename):
                msgs.error(f"{name} does not exist in one of the provided paths.  Modify your input {self.flavor} file")
            data_files.append(filename)

        # Return
        return data_files

    def write(self, input_file, version_override=None, date_override=None):
        """
        Write an Input file to disk

        Args:
            input_file (:obj:`str`):
                Name of PypeIt file to be generated
            version_override (:obj:`str`, optional):
                Override the current version and use this one instead.  **For
                documentation purposes only!**
            date_override (:obj:`str`, optional):
                Override the current date and use this one instead.  **For
                documentation purposes only!**
        """
        _version = __version__ if version_override is None else version_override
        _date = datetime.now(__UTC__).isoformat(timespec='milliseconds') \
                    if date_override is None else date_override

        # Here we go
        with open(input_file, 'wb') as bf:

            # We use ConfigObj to write the original config lines with comments.
            # But it wants a binary stream.
            # So use a TextIOWrapper to make it easier to write to it.
            with io.TextIOWrapper(bf, encoding='utf-8') as f:

                # The ConfigObj will have the original comments so new ones
                # aren't needed
                if not self.preserve_comments:
                    f.write(f'# Auto-generated {self.flavor} input file using PypeIt version: {_version}\n')
                    f.write(f'# UTC {_date}\n')
                    f.write("\n")

                # Parameter block
                if self.config is not None:
                    if not self.preserve_comments:
                        f.write("# User-defined execution parameters\n")
                        f.write('\n'.join(self.cfg_lines))
                        f.write('\n')
                        f.write('\n')
                    else:
                        self.config.write(bf)

                        # Flush the binary stream just to make sure future writes to the text
                        # stream will be after it in the output file.
                        bf.flush()

                # Setup block
                if self.setup is not None:
                    setup_lines = yaml.dump(utils.yamlify(
                        self.setup)).split('\n')[:-1]
                elif self.setup_required: # Default
                    setup_lines = ['Setup A:']
                else:
                    setup_lines = None

                if setup_lines is not None:
                    if not self.preserve_comments:
                        # This comment is part of the configuration lines,
                        # and so is preserved by ConfigObj
                        f.write("# Setup\n")
                    f.write("setup read\n")
                    f.write('\n'.join(setup_lines)+'\n')
                    f.write("setup end\n")
                    f.write("\n")
                
                # Data block                
                if self.data is not None or self.datablock_required:
                    f.write("# Data block \n")
                    f.write(f"{self.data_block} read\n")
                    # paths and Setupfiles
                    if self.file_paths is not None:
                        for path in self.file_paths:
                            f.write(' path '+path+'\n')
                    if self.data is not None:
                        with io.StringIO() as ff:
                            self.data.write(ff, format='ascii.fixed_width',
                                            bookend=False)
                            data_lines = ff.getvalue().split('\n')[:-1]
                        f.write('\n'.join(data_lines))
                        f.write('\n')
                    f.write(f"{self.data_block} end\n")
                    f.write("\n")

        msgs.info(f'{self.flavor} input file written to: {input_file}')

    def get_spectrograph(self):
        """
        Use the configuration lines to instantiate the relevant
        :class:`~pypeit.spectrographs.spectrograph.Spectrograph` subclass.

        Returns:
            :class:`~pypeit.spectrographs.spectrograph.Spectrograph`:
            Spectrograph subclass instance using the relevant configuration
            parameter.

        Raises:
            :class:`~pypeit.pypmsgs.PypeItError`:
                Raised if the relevant configuration parameter is not available.
        """
        if 'rdx' not in self.config.keys() or 'spectrograph' not in self.config['rdx'].keys():
            msgs.error('Cannot define spectrograph.  Configuration file missing \n'
                       '    [rdx]\n    spectrograph=\n entry.')
        return load_spectrograph(self.config['rdx']['spectrograph'])

    def get_pypeitpar(self, config_specific_file=None):
        """
        Use the configuration lines and a configuration-specific example file to
        build the full parameter set.

        Args:
            config_specific_file (:obj:`str`, `Path`_, optional):
                The file to use to generate the default, configuration-specific
                parameters.  If None and instance contains filenames, use the
                first file.  If None and instance provides no filenames,
                configuration-specific parameters are not set.

        Returns:
            :obj:`tuple`: A tuple with the spectrograph instance, the
            parameters, and the file name used to generate the
            configuration-specific parameters.  That latter will be None if the
            no example file was available.
        """
        spec = self.get_spectrograph()

        if config_specific_file is None:
            _files = self.filenames
            if _files is not None:
                config_specific_file = _files[0]

        spec_par = spec.default_pypeit_par() if config_specific_file is None \
                    else spec.config_specific_par(config_specific_file)
        par = PypeItPar.from_cfg_lines(cfg_lines=spec_par.to_config(),
                                       merge_with=(self.cfg_lines,))
        return spec, par, config_specific_file


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

    def get_pypeitpar(self):
        """
        Override the base class function to use files with specific frametypes
        for the config-specific parameters.

        Returns:
            :obj:`tuple`: A tuple with the spectrograph instance, the
            parameters, and the file name used to generate the
            configuration-specific parameters.  That latter will be None if the
            no example file was available.
        """
        if 'frametype' not in self.data.keys():
            msgs.error('PypeIt file must provide the frametype column.')

        # NOTE: self.filenames is a property function that generates the full
        # set of file names each time they are requested.  However, this should
        # only be done once in the code below because as soon as a relevant file
        # is found the loops are discontinued using `break`.

        config_specific_file = None
        for idx, row in enumerate(self.data):
            if 'science' in row['frametype'] or 'standard' in row['frametype']:
                config_specific_file = self.filenames[idx]
                break

        # If no science/standard frames available, search for an arc/trace
        # instead.
        if config_specific_file is None:
            for idx, row in enumerate(self.data):
                if 'arc' in row['frametype'] or 'trace' in row['frametype']:
                    config_specific_file = self.filenames[idx]
                    break

        if config_specific_file is not None:
            self.get_spectrograph()._check_extensions(config_specific_file)
        return super().get_pypeitpar(config_specific_file=config_specific_file)


class SensFile(InputFile):
    """Child class for the Sensitivity input file
    """
    data_block = 'sens'  # Defines naming of data block
    flavor = 'Sens'  # Defines naming of file
    datablock_required = False
    setup_required = False


class ExtractFile(InputFile):
    """Child class for the Extraction input file
    """
    data_block = 'extract'  # Defines naming of data block
    flavor = 'Extract'  # Defines naming of file
    setup_required = False
    datablock_required = False


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

    # TODO is the correct way to treat optional table entries?

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

        if 'sensfile' not in self.data.keys():
            return None

        # Grab em
        sens_files = self.path_and_files('sensfile', skip_blank=True)
        # Pad out
        if len(sens_files) == 1 and len(self.filenames) > 1:
            sens_files = sens_files * len(self.filenames)
            # Return
        return sens_files


    @property
    def setup_id(self):

        if 'setup_id' not in self.data.keys():
            return None

        # Generate list, scrubbing empty entries
        sid = [str(item) for item in self.data['setup_id'] if str(item).strip() not in ['', 'none', 'None']]

        # Inflate as needed
        if len(sid) == 1 and len(sid) < len(self.data):
            sid = sid * len(self.data)
        # Return
        return sid





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
        msgs.info('.coadd3d file successfully vetted.')

    @property
    def options(self):
        """
        Parse the options associated with a cube block.
        Here is a description of the available options:

        - ``sensfunc``: The name of an a sensitivity function file that is used
            for the flux calibration.  The file provided here should be generated by
            (or of the same format as the output of) the command :ref:`pypeit_sensfunc`.
            This parameter can also be set for all frames
            with the default command:

            .. code-block:: ini

                [reduce]
                    [[cube]]
                        sensfile = sens.fits

        - ``scale_corr``: The name of an alternative spec2d file that is used for
          the relative spectral scale correction.  This parameter can also be set
          for all frames with the default command:

          .. code-block:: ini

                [reduce]
                    [[cube]]
                        scale_corr = spec2d_alternative.fits

        - ``grating_corr``: The name of a Flat calibrations file that is used
            for the grating tilt correction.  This parameter can also be set for all frames
            with the default command:

            .. code-block:: ini

                [reduce]
                    [[cube]]
                        grating_corr = Flat_A_0_DET01.fits

        - ``skysub_frame``: The name of an alternative spec2d file that is used
          for the sky subtraction.  This parameter can also be set for all frames
          with the default command:

          .. code-block:: ini

                [reduce]
                    [[cube]]
                        skysub_frame = spec2d_alternative.fits

        - ``ra_offset``: The RA offset to apply to the WCS of the cube.

        - ``dec_offset``: The DEC offset to apply to the WCS of the cube.


        Returns
        -------
        opts: dict
            Dictionary containing cube options.
        """
        # Define the list of allowed parameters
        opts = dict(sensfile=None, scale_corr=None, grating_corr=None, skysub_frame=None,
                    ra_offset=None, dec_offset=None)

        # Get the sensfunc files
        sensfile = self.path_and_files('sensfile', skip_blank=False, check_exists=False)
        if sensfile is None:
            opts['sensfile'] = None
        elif len(sensfile) == 1 and len(self.filenames) > 1:
            opts['sensfile'] = sensfile*len(self.filenames)
        elif len(sensfile) != 0:
            opts['sensfile'] = sensfile

        # Get the scale correction files
        scale_corr = self.path_and_files('scale_corr', skip_blank=False, check_exists=False)
        if scale_corr is None:
            opts['scale_corr'] = [None]*len(self.filenames)
        elif len(scale_corr) == 1 and len(self.filenames) > 1:
            opts['scale_corr'] = scale_corr*len(self.filenames)
        elif len(scale_corr) != 0:
            opts['scale_corr'] = scale_corr

        # Get the grating correction files
        grating_corr = self.path_and_files('grating_corr', skip_blank=False, check_exists=False)
        if grating_corr is None:
            opts['grating_corr'] = [None]*len(self.filenames)
        elif len(grating_corr) == 1 and len(self.filenames) > 1:
            msgs.error("You cannot specify a single grating correction file for multiple input files.")
        elif len(grating_corr) != 0:
            opts['grating_corr'] = grating_corr

        # Get the skysub files
        skysub_frame = self.path_and_files('skysub_frame', skip_blank=False, check_exists=False)
        if skysub_frame is None:
            opts['skysub_frame'] = ["default"]*len(self.filenames)
        elif len(skysub_frame) == 1 and len(self.filenames) > 1:
            opts['skysub_frame'] = skysub_frame*len(self.filenames)
        elif len(skysub_frame) != 0:
            opts['skysub_frame'] = skysub_frame

        # Load coordinate offsets for each file. This is "Delta RA cos(dec)" and "Delta Dec"
        # Get the RA offset of each file
        off_ra, off_dec = None, None
        if 'ra_offset' in self.data.keys():
            off_ra = self.data['ra_offset'].tolist()
            if len(off_ra) == 1 and len(self.filenames) > 1:
                # Convert from arcsec to degrees
                opts['ra_offset'] = [off_ra[0]/3600.0 for _ in range(len(self.filenames))]
            elif len(off_ra) != 0:
                # Convert from arcsec to degrees
                opts['ra_offset'] = [ora/3600.0 for ora in off_ra]
        # Get the DEC offset of each file
        if 'dec_offset' in self.data.keys():
            off_dec = self.data['dec_offset'].tolist()
            if len(off_dec) == 1 and len(self.filenames) > 1:
                # Convert from arcsec to degrees
                opts['dec_offset'] = [off_dec[0]/3600.0 for _ in range(len(self.filenames))]
            elif len(off_dec) != 0:
                # Convert from arcsec to degrees
                opts['dec_offset'] = [odec/3600.0 for odec in off_dec]
        # Check that both have been set or both are not set
        if (off_ra is not None and off_dec is None) or (off_ra is None and off_dec is not None):
            msgs.error("You must specify both or neither of the following arguments: ra_offset, dec_offset")

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

class RawFiles(InputFile):
    """Child class for a list of raw files
    """
    data_block = 'raw'  # Defines naming of data block
    flavor = 'Raw'      
    setup_required = False
    datablock_required = True
    required_columns = ['filename'] 

    def vet(self):
        """ Check for required bits and pieces of the .coadd2d file
        besides the input objects themselves
        """
        super().vet()

        # Done
        msgs.info('.rawfiles file successfully vetted.')


# NOTE: I originally had this in pypeit/io.py, but I think it was causing a
# circular import.  Moving it here solved the issue.
def grab_rawfiles(file_of_files:str=None, list_of_files:list=None, raw_paths:list=None, 
                  extension:str='.fits'):
    """
    Parse a set of raw files from the input.

    Although all arguments are optional, one of ``file_of_files``,
    ``list_of_files``, or ``raw_paths`` must be true.  Precedence is given in
    that order; i.e., if ``file_of_files`` is provided, all other arguments are
    ignored.

    Args:
        file_of_files (str, optional): 
            File with list of raw files.  Format must follow the
            :ref:`input-files-data-block` of a PypeIt file, and the only
            required column is the filename.
        list_of_files (list, optional): 
            List of raw files (str).  Ignored if ``file_of_files`` is provided.
            If ``raw_paths`` is None, the path is assumed to be the current
            working directory.
        raw_paths (list, optional): 
            One or more paths with the raw files.  Ignored if ``file_of_files``
            is provided.  If ``list_of_files`` is None, all files with the
            provided extension are assumed to be raw files.
        extension (str, optional): 
            File extension to search on.  Ignored if ``file_of_files`` or
            ``list_of_files`` is provided.

    Returns:
        list: List of raw data filenames with full path
    """
    if file_of_files is not None:
        # PypeIt formatted list of files
        return RawFiles.from_file(file_of_files).filenames

    _raw_paths = [Path().absolute()] if raw_paths is None \
                    else [Path(p).absolute() for p in raw_paths]

    if list_of_files is not None:
        # An actual list
        return [str(p / f) for p in _raw_paths for f in list_of_files if (p / f).exists()]

    # Find all files that have the correct extension
    return files_from_extension(_raw_paths, extension=extension)

