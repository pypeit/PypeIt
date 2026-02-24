""" Class for I/O of PypeIt input files

.. include:: ../include/links.rst
"""
import datetime
import io
from pathlib import Path

from astropy.table import Table, column
from astropy.io import ascii
import configobj
from IPython import embed
import numpy as np
import yaml

from pypeit import log
from pypeit import PypeItError
from pypeit import utils
from pypeit import __version__
from pypeit.io import files_from_extension
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph


class InputFile:
    """
    Generic class to load, process, and write PypeIt input files

    In practice, we use one of the children of this class,
    e.g. PypeItFile

    This class has limited support for preserving comments in the input files.
    We currently support comments in the configuration section, and commented
    out data lines in the data sections. Other comments are not preserved.

    Parameters
    ----------
    config : :obj:`dict`, :obj:`list`, optional
        Configuration dictionary or list of config lines.  This is converted to
        a `ConfigObj`_ instance.
    file_paths : list, optional
        List of string or `Path`_ objects with the file paths for the data files
    data_table : `astropy.table.Table`_, optional
        Data block
    setup : :obj:`dict`, optional
        Dictionary with the instrument setup/configuration.  The first key
        contains the setup name.
    vet : bool, optional
        If True, vet the file after iniitialization.  The :func:`vet` method can
        be called after initialization, if needed.
    preserve_comments : bool, optional
        If True, (try to) preserve the comments in the input file.

    Attributes
    ----------
    data : `astropy.table.Table`_
        Table with the data in the data block
    file_paths : list
        List of `Path` objects with the directory paths for data files
    setup : dict
        Dictionary with the instrument setup/configuration.
    preserve_comments : bool
        Flag to preserve comments in the input file.
    config : `ConfigObj`_
        Parameters parsed from the configuration block of the input file.
    """

    flavor = 'Generic'
    """
    Defines the type of the input file.
    """

    required_columns = []
    """
    Sets the required columns in the data table.
    """

    data_block = None
    """
    Defines the name of data block.
    """

    datablock_required = False
    """
    Sets whether or not the data block is required.
    """
    
    setup_required = False
    """
    Sets whether or not the setup block is required.
    """

    def __init__(
        self, 
        config=None, 
        file_paths:list=None,
        data_table:Table=None,
        setup:dict=None,
        vet:bool=True,
        preserve_comments:bool=False
    ):
        # Load up
        self.data = data_table
        self.file_paths = (
            None if file_paths is None else [Path(fp).absolute() for fp in file_paths]
        )
        self.setup = setup
        self.preserve_comments = preserve_comments
        self.config = None if config is None else configobj.ConfigObj(config)

        # Vet
        if vet:
            self.vet()

    @staticmethod
    def remove_comments_and_blanks(lines : np.ndarray) -> np.ndarray:
        # Remove comment and blank lines
        lines = lines[np.array([ len(l.strip()) > 0 and l.strip()[0] != '#' for l in lines ])]

        # Remove appended comments and return
        return np.array([l.split('#')[0].strip() for l in lines])
        
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
        _ifile = Path(ifile).absolute()
        if not _ifile.is_file():
            raise FileNotFoundError(f'Input file {ifile} does not exist!')

        # Read the input lines and replace special characters
        with open(_ifile, 'r') as f:
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
        log.info('Loading the reduction file')
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
            raise PypeItError(
                f"Missing '{cls.data_block} end' in {input_file}")
        if data_start < 0 and data_end>0:
            raise PypeItError("You have not specified the start of the data block!")
        # Read it, if it exists
        if data_start>0 and data_end>0:
            paths, usrtbl = cls._read_data_file_table(lines[data_start:data_end], preserve_comments)
            is_config[data_start-1:data_end+1] = False
            data_block_found = True
        else:
            if cls.datablock_required:
                raise PypeItError("You have not specified the data block!")            
            paths, usrtbl = [], None
            data_block_found = False

        # Parse the setup block
        setup_found = False
        setup_start, setup_end = cls.find_block(lines, 'setup')
        if setup_start >= 0 and setup_end < 0 and cls.setup_required:
            raise PypeItError(f"Missing 'setup end' in {input_file}")
        elif setup_start < 0 and cls.setup_required:
            raise PypeItError(f"Missing 'setup read' in {input_file}")
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
        log.info(f'{cls.flavor} input file loaded successfully.')

        # Instantiate
        return cls(
            config=list(lines[is_config]), 
            file_paths=paths, 
            data_table=usrtbl, 
            setup=sdict,
            vet=vet,
            preserve_comments=preserve_comments
        )

    def vet(self):
        """ Check for required bits and pieces of the Input file
        besides the input objects themselves
        """
        # Data table
        if self.data is None:
            if self.datablock_required:
                raise PypeItError("You have not specified the data block!")
        else:
            for key in self.required_columns:
                if key not in self.data.keys():
                    raise PypeItError(f'Add {key} to the Data block of your {self.flavor} file before running.')

        if self.setup_required and self.setup is None:
            raise PypeItError("Add setup info to your PypeIt file in the setup block!")

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
    def filenames(self) -> list[str]:
        """ List of path + filename's
        Wrapper to :func:`~pypeit.inputfiles.InputFile.path_and_files`.
        See that function for a full description.

        Returns:
            list or None: List of the full paths to each data file
            or None if `filename` is not part of the data table
            or there is no data table!
        """
        # Return
        return self.path_and_files('filename', include_commented_out=self.preserve_comments)

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
            raise PypeItError("Setup block contains more than one Setup!")

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
            log.warning(
                'Your PypeIt file has leading | characters in the data table, which is the old '
                'format.  Please update your file to remove leading and trailing | characters '
                'because their inclusion will be deprecated.'
            )
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

    def path_and_files(
            self, key:str, skip_blank=False, include_commented_out=False, check_exists=True
    ) -> list[str]:
        """
        Generate a list of the filenames with the full path from the column of
        the data `astropy.table.Table`_ specified by ``key``.  The files must
        exist and be within one of the paths for this to succeed.

        Parameters
        ----------
        key : str
            Column in :attr:`data` with the filenames of interest
        skip_blank : bool, optional
            If True, ignore any entry that is '', 'none' or 'None'.
        check_exists : bool, optional
            If False, PypeIt will not check if the files exist.
        include_commented_out : bool, optional
            If False, commented out files will not be included. If True, they
            are included, without the "#" character.

        Returns
        -------
        list

            List of strings, where each string provides the full path to each
            data file.  None is returned if :attr:`data` is not defined or if
            ``key`` is not one of its columns.
              
        Raises
        ------
        FileNotFoundError:
            Raised if ``check_exists`` is True and any of the files do not exist.
        """
        if self.data is None or key not in self.data.keys():
            return None

        ## Build full paths to file and set frame types
        data_files = []
        for row in self.data:

            # Skip Empty entries?
            if skip_blank and row[key].strip() in ['', 'none', 'None']:
                continue

            # Skip commented out entries
            if row[key].strip().startswith("#"):
                if not include_commented_out:
                    continue
                # Strip the comment character and any whitespace following it
                # from the filename
                name = row[key].strip("# ")
            else:
                name = row[key]

            # Searching..
            if self.file_paths is not None and len(self.file_paths) > 0:
                for p in self.file_paths:
                    filename = p / name
                    if filename.is_file():
                        break
            else:
                filename = Path(name).absolute()

            # Check we got a good hit
            if check_exists and not filename.is_file():
                raise FileNotFoundError(
                    f'{name} does not exist in one of the provided paths.  Modify your input '
                    f'{self.flavor} file.'
                )
            data_files.append(str(filename))

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
        _date = datetime.datetime.now(datetime.UTC).isoformat(timespec='milliseconds') \
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
                            f.write(f' path {path}\n')
                    if self.data is not None:
                        with io.StringIO() as ff:
                            self.data.write(ff, format='ascii.fixed_width',
                                            bookend=False)
                            data_lines = ff.getvalue().split('\n')[:-1]
                        f.write('\n'.join(data_lines))
                        f.write('\n')
                    f.write(f"{self.data_block} end\n")
                    f.write("\n")

        log.info(f'{self.flavor} input file written to: {input_file}')

    # TODO: Should we add the spectrograph instance as an attribute?
    def get_spectrograph(self, spectrograph_name:str=None, pypeit_fits:bool=False):
        """
        Use the configuration lines to instantiate the relevant
        :class:`~pypeit.spectrographs.spectrograph.Spectrograph` subclass.

        Parameters
        ----------
        spectrograph_name : str, optional
            The name of the spectrograph to load.  If None, the spectrograph
            name is pulled from the configuration.  If the configuration is None
            or the ``['rdx']['spectrograph']`` configuration entry is not
            available, an exception is raised.
        pypeit_fits : :obj:`bool`, optional
            The spectrograph loader is being called from a post-processing
            script where the expected input files are PypeIt-written FITS files
            only.  This has the effect of overriding the :attr:`allowed_extensions`
            attribute to be ``[".fits"]``.

        Returns
        -------
        :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
            Spectrograph subclass instance using the relevant configuration
            parameter.

        Raises
        ------
        :class:`~pypeit.exceptions.PypeItError`
            Raised if the configuration or the relevant parameter is not available.
        """
        if spectrograph_name is None:
            if (
                self.config is None
                or 'rdx' not in self.config.keys()
                or 'spectrograph' not in self.config['rdx'].keys()
            ):
                raise PypeItError(
                    'Cannot define spectrograph.  Configuration block is missing or does not '
                    'include the spectrograph definition (in the rdx parameter block).'
                )
            _spec = self.config['rdx']['spectrograph']
        else:
            _spec = spectrograph_name

        return load_spectrograph(_spec, pypeit_fits=pypeit_fits)

    def get_pypeitpar(
        self, config_specific_file=None, spectrograph_name:str=None, pypeit_fits:bool=False
    ):
        """
        Use the configuration lines and a configuration-specific example file to
        build the full parameter set.

        Parameters
        ----------
        config_specific_file : :obj:`str`, :obj:`list`, `Path`_, `astropy.io.fits.Header`_, `astropy.table.Table`_, optional
            The object used to generate the default, configuration-specific
            parameters using
            :func:`~pypeit.spectrographs.spectrograph.Spectrograph.config_specific_par`.
            If None and :attr:`filenames` are available, use the first file.  If
            None and :attr:`filenames` are *not* available, use
            :func:`~pypeit.spectrographs.spectrograph.Spectrograph.default_pypeit_par`
            to generate the parameters instead.
        spectrograph_name : str, optional
            The name of the spectrograph to load.  If None, the spectrograph
            name is pulled from the configuration.  If the configuration is None
            or the ``['rdx']['spectrograph']`` configuration entry is not
            available, an exception is raised.
        pypeit_fits : :obj:`bool`, optional
            The spectrograph loader is being called from a post-processing
            script where the expected input files are PypeIt-written FITS files
            only.  This has the effect of overriding the :attr:`allowed_extensions`
            attribute to be ``[".fits"]``.

        Returns
        -------
        spec : :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
            Spectrograph subclass instance.
        par : :class:`~pypeit.par.pypeitpar.PypeItPar`
            Parameter set generated by merging the configuration-specific parameters
            with the configuration parameters.
        config_specific_file : :obj:`str`, :obj:`list`, `Path`_, `astropy.io.fits.Header`_, `astropy.table.Table`_
            The object used to generate the configuration-specific parameters,
            which can be None.
        """
        spec = self.get_spectrograph(spectrograph_name=spectrograph_name, pypeit_fits=pypeit_fits)

        if config_specific_file is None:
            _files = self.filenames
            if _files is not None:
                config_specific_file = _files[0]

        # Get the configuration-specific parameters based on the file
        spec_par = spec.default_pypeit_par() if config_specific_file is None \
                    else spec.config_specific_par(config_specific_file)
        par = pypeitpar.PypeItPar.from_cfg_lines(
            cfg_lines=spec_par.to_config(), merge_with=(self.cfg_lines,)
        )
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
            raise PypeItError(f"Missing spectrograph in the Parameter block of your PypeIt file.  Add it!")

        # Setup
        setup_keys = list(self.setup)
        if 'Setup' not in setup_keys[0]:
            raise PypeItError("Setup does not appear in your setup block! Add it")

        # Done
        log.info('PypeIt file successfully vetted.')

    @property
    def frametypes(self):
        """Return a dict of the frametypes
        with key, item the filename, frametype 

        Returns:
            dict: 
        """
        return {row['filename']:row['frametype'] for row in self.data}

    def get_pypeitpar(
        self, config_specific_file=None, spectrograph_name:str=None, pypeit_fits:bool=False
    ):
        """
        Use the configuration lines and a configuration-specific example file to
        build the full parameter set.

        This overrides the base class function to use files with specific
        frametypes and the contents of the PypeIt File (including and
        modifications away from values in the FITS headers) for the
        config-specific parameters.

        .. warning::

            The calling sequence is purposely kept consistent with the base
            class, however only the ``config_specific_file`` is used.  An
            exception will be raised if ``spectrograph_name`` is not None or
            ``pypeit_fits`` is True.

        Parameters
        ----------
        config_specific_file : :obj:`str`, `Path`_, optional
            The file used to generate the default, configuration-specific
            parameters.  Note that this *must* be a file name, unlike the type
            options available to the base class.  The behavior of this parameter
            is defined as follows:
                
                - If None and filenames are present, the file to use is
                  determined as follows:

                    - If the frametypes are not defined or there are no science,
                      standard, arc, or trace frames, the first data file is
                      used.

                    - If the frametypes are defined, the first science or
                      standard frame is used.

                    - If the frametypes are defined and there are no science or
                      standard frames, the first arc or trace frame is used.

                - If not None, the filenames are expected to be present and the
                  provided string must be a file that is listed by the data
                  table.  If it is not, an exception is raised.

        spectrograph_name : str, optional
            Provided for consistency with the base class.  This must always be
            ``None``.  The name of the spectrograph in pypeit files is always
            defined by the "rdx" parameter set.
        pypeit_fits : :obj:`bool`, optional
            Provided for consistency with the base class.  The files provided by
            the pypeit file are, by definition, raw files, meaning this
            parameter must always be False.

        Returns
        -------
        spec : :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
            Spectrograph subclass instance.
        par : :class:`~pypeit.par.pypeitpar.PypeItPar`
            Parameter set used to controle the data processing.  If a valid
            ``config_specific_file`` is defined (see above), the paraemters are
            determined using
            :func:`~pypeit.spectrographs.spectrograph.Spectrograph.config_specific_par`.
            Otherwise, the parameters are determined by
            :func:`~pypeit.spectrographs.spectrograph.Spectrograph.default_pypeit_par`
        config_specific_file : `Path`_
            The path to the file used to generate the configuration-specific
            parameters, which can be None.
        """
        # Check the input
        if spectrograph_name is not None:
            raise PypeItError(
                f'For instances of {self.__class__.__name__}, the spectrograph must always be '
                'defined by the parameter block.  Calls to `get_pypeitpar` cannot independently '
                'define the spectrograph name.'
            )
        if pypeit_fits:
            raise PypeItError(
                f'For instances of {self.__class__.__name__}, the provided files are, by '
                'definition, raw data files.  Calls to `get_pypeitpar` cannot indicate that the '
                'files are output FITS files produced by PypeIt.'
            )

        # NOTE: self.filenames is a property function that generates the full
        # set of file names each time they are requested.  Generate the set here
        # so that it only needs to be done once.
        try:
            filenames = self.filenames
        except Exception as e:
            log.warning(
                'Unable to define filenames.  Proceeding by assuming filenames are not defined.  '
                f'Original exception raised is: {e}'
            )
            filenames = None

        if config_specific_file is None:
            if filenames is None:
                _config_specific_file = None
            elif 'frametype' not in self.data.keys():
                _config_specific_file = filenames[0]
            else:
                # Search for the first science/standard frame
                _config_specific_file = None
                csf_indx = None
                for idx, row in enumerate(self.data):
                    if 'science' in row['frametype'] or 'standard' in row['frametype']:
                        _config_specific_file = filenames[idx]
                        csf_indx = np.array([idx])
                        break

                # If no science/standard frames available, search for an arc/trace
                # instead.
                if _config_specific_file is None:
                    for idx, row in enumerate(self.data):
                        if 'arc' in row['frametype'] or 'trace' in row['frametype']:
                            _config_specific_file = filenames[idx]
                            csf_indx = np.array([idx])
                            break
        else:
            if filenames is None:
                raise PypeItError(
                    f'When providing a config_specific_file for a {self.__class__.__name__}, a '
                    'data table with a set of raw files must be included.'
                )
            _config_specific_file = Path(config_specific_file).absolute()
            csf_indx = np.where(self.data['filename'] == _config_specific_file.name)[0]
            if len(csf_indx) == 0:
                raise PypeItError(
                    f'The provided file ({config_specific_file}) is not one of the files '
                    'included in the data table.'
                )
            if len(csf_indx) > 1:
                raise PypeItError(
                    f'The provided file ({config_specific_file}) matches to more than one file '
                    'in the data table.'
                )

        # Load the spectrograph
        spec = self.get_spectrograph()

        if _config_specific_file is None:
            spec_par = spec.default_pypeit_par() 
        else:
            # Check file extensions
            spec._check_extensions(_config_specific_file)

            # Send the Row of the metadata table corresponding to the file
            # NOTE: csf_indx needs to be a single element array/list such that
            # `self.data[csf_indx]` returns a Table object.  If it is an
            # integer, `self.data[csf_indx]` returns a Row object, which does
            # not have a `copy()` method.
            data_row = self.data[csf_indx].copy()
            # Use the full path to the ``config_specific_file`` for insurance
            data_row['filename'] = str(_config_specific_file)
            spec_par = spec.config_specific_par(data_row)

        par = pypeitpar.PypeItPar.from_cfg_lines(
            cfg_lines=spec_par.to_config(), merge_with=(self.cfg_lines,)
        )
        return spec, par, _config_specific_file        


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
            log.warning("sensfile column not provided.  Fluxing will crash if an archived sensitivity function does not exist")
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
            raise PypeItError(f"Missing spectrograph in the Parameter block of your .coadd2d file.  Add it!")

        # Done
        log.info('.coadd2d file successfully vetted.')


class Coadd3DFile(InputFile):
    """Child class for coadding spec2d files into datacubes
    """
    data_block = 'spec2d'  # Defines naming of data block
    flavor = 'Coadd3d'  # Defines naming of file
    setup_required = False
    datablock_required = True
    required_columns = ['filename'] 

    def vet(self):
        """ Check for required bits and pieces of the .coadd3d file
        besides the input objects themselves
        """
        super().vet()

        # Confirm spectrograph is present
        if 'rdx' not in self.config.keys() or 'spectrograph' not in self.config['rdx'].keys():
            raise PypeItError(f"Missing spectrograph in the Parameter block of your .coadd2d file.  Add it!")

        # Done
        log.info('.coadd3d file successfully vetted.')

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
            raise PypeItError("You cannot specify a single grating correction file for multiple input files.")
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
            raise PypeItError("You must specify both or neither of the following arguments: ra_offset, dec_offset")

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
            raise PypeItError(f"Missing spectrograph in the Parameter block of your .flex file.  Add it!")

        # Done
        log.info('.flex file successfully vetted.')


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

        Allows for wildcards

        Returns:
            list or None: List of the full paths to each data file
            or None if `filename` is not part of the data table
            or there is no data table!
        """
        all_files = []
        paths = (
            [Path().absolute()]
            if self.file_paths is None or len(self.file_paths) == 0
            else self.file_paths
        )
        # Paths?
        for p in paths:
            for row in self.data['filename']:
                all_files += sorted(p.glob(row))

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
        log.info('.rawfiles file successfully vetted.')


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

    # Find all files that have the correct extension.  Force the returned list
    # to contain strings, not Path objects.
    return [str(f) for f in files_from_extension(_raw_paths, extension=extension)]
