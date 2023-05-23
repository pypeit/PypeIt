"""
Class for organizing PypeIt setup

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
from pathlib import Path
import time
import os

from IPython import embed

from pypeit import msgs
from pypeit.metadata import PypeItMetaData
from pypeit import inputfiles
from pypeit.par import PypeItPar
from pypeit import io
from pypeit.spectrographs.util import load_spectrograph


class PypeItSetup:
    """
    Prepare for a pypeit run.

    The main deliverables are the set of parameters used for PypeIt's algorithms
    (:attr:`par`), a table with metadata pulled from the files to be reduced
    (:attr:`fitstbl`), and a class used to perform instrument-specific
    operations (:attr:`spectrograph`).

    Args:
        file_list (list):
            A list of strings with the full path to each file to be
            reduced.
        frametype (:obj:`dict`, optional):
            A dictionary that associates the name of the file (just the fits
            file name without the full path) to a specific frame type (e.g.,
            arc, bias, etc.).  The file name and type are expected to be the key
            and value of the dictionary, respectively.  If None, this is
            determined by the :func:`get_frame_types` method.
        usrdata (:obj:`astropy.table.Table`, optional):
            A user provided set of data used to supplement or overwrite
            metadata read from the file headers.  The table must have a
            `filename` column that is used to match to the metadata
            table generated within PypeIt.
        setups (:obj:`list`, optional):
            A list of setups that each file can be associated with.  If
            None, all files are expected to be for a single setup.
        cfg_lines (:obj:`list`, optional):
            A list of strings that provide a set of user-defined
            parameters for executing pypeit.  These are the lines of a
            configuration file.  See the documentation for the
            `configobj`_ package.  One of the user-level inputs should
            be the spectrograph that provided the data to be reduced.
            One can get the list of spectrographs currently served by
            running::
                
                from pypeit.spectrographs.util import valid_spectrographs
                print(valid_spectrographs())

            To use all the default parameters when reducing data from a
            given spectrograph, you can set `cfg_lines = None`, but you
            then *must* provide the `spectrograph_name`.
        spectrograph_name (:obj:`str`, optional):
            If not providing a list of configuration lines
            (`cfg_lines`), this sets the spectrograph.  The spectrograph
            defined in `cfg_lines` takes precedent over anything
            provided by this argument.
        pypeit_file (:obj:`str`, optional):
            The name of the pypeit file used to instantiate the
            reduction. This can be None, and will lead to default
            names for output files (TODO: Give list). Setting
            :ref:`pypeit_file` here *only sets the name of the file*.
            To instantiate a :class:`~pypeit.pypeitsetup.PypeItSetup`
            object directly from a pypeit file (i.e. by reading the
            file), use the :func:`from_pypeit_file` method; i.e.::
                
                setup = PypeItSetup.from_pypeit_file('myfile.pypeit')

    Attributes:
        file_list (list):
            See description of class argument.
        frametype (dict):
            See description of class argument.
        setups (list):
            See description of class argument.
        pypeit_file (str):
            See description of class argument.
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            An instance of the `Spectrograph` class used throughout the
            reduction procedures.
        par (:class:`pypeit.par.pypeitpar.PypeItPar`):
            An instance of the `PypeitPar` class that provides the
            parameters to all the algorthms that pypeit uses to reduce
            the data.
        fitstbl (:class:`pypeit.metadata.PypeItMetaData`):
            A `Table` that provides the salient metadata for the fits
            files to be reduced.
    """
    def __init__(self, file_list, frametype=None, usrdata=None, setups=None, cfg_lines=None,
                 spectrograph_name=None, pypeit_file=None):

        # The provided list of files cannot be None
        if file_list is None or len(file_list) == 0:
            msgs.error('Must provide a list of files to be reduced!')

        # Save input
        self.file_list = file_list
        self.frametype = frametype
        self.usrdata = usrdata
        self.setups = setups
        self.user_cfg = cfg_lines

        # Determine the spectrograph name
        _spectrograph_name = spectrograph_name if cfg_lines is None \
                    else PypeItPar.from_cfg_lines(merge_with=(cfg_lines,))['rdx']['spectrograph']

        # Cannot proceed without spectrograph name
        if _spectrograph_name is None:
            msgs.error('Must provide spectrograph name directly or using configuration lines.')
       
        # Instantiate the spectrograph
        self.spectrograph = load_spectrograph(_spectrograph_name)

        # Set the pypeit file name, if none is provided
        self.pypeit_file = self.spectrograph.name + '.pypeit' \
                                if pypeit_file is None else pypeit_file

        # Get the spectrograph specific configuration to be merged with
        # the user modifications.
        spectrograph_cfg_lines = self.spectrograph.default_pypeit_par().to_config()

        # Instantiate the pypeit parameters.  The user input
        # configuration (cfg_lines) can be None.
        self.par = PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines, 
                                            merge_with=(cfg_lines,))

        # Prepare internals for execution
        self.fitstbl = None

    @classmethod
    def from_pypeit_file(cls, filename):
        """
        Instantiate the :class:`PypeitSetup` object using a pypeit file.

        Args:
            filename (str):
                Name of the pypeit file to read.  PypeIt files have a
                specific set of valid formats. A description can be
                found :ref:`pypeit_file`.
        
        Returns:
            :class:`PypeitSetup`: The instance of the class.
        """
        # Load up a pypeItFile object
        pypeItFile = inputfiles.PypeItFile.from_file(filename)
        # Instantiate
        return cls(pypeItFile.filenames, 
                   frametype=pypeItFile.frametypes, 
                   usrdata=pypeItFile.data, 
                   setups=[pypeItFile.setup_name],
                   cfg_lines=pypeItFile.cfg_lines, 
                   pypeit_file=filename)

    # TODO: Make the default here match the default used by
    # io.files_from_extension?
    @classmethod
    def from_file_root(cls, root, spectrograph, extension='.fits'):
        """
        Instantiate the :class:`PypeItSetup` object by providing a file
        root.
        
        Args:
            root (:obj:`str`):
                String used to find the raw files; see
                :func:`~pypeit.io.files_from_extension`.
            spectrograph (:obj:`str`):
                The PypeIt name of the spectrograph used to take the
                observations.  This should be one of the available
                options in
                :func:`~pypeit.spectrographs.available_spectrographs`.
            extension (:obj:`str`, optional):
                The extension common to all the fits files to reduce; see
                :func:`~pypeit.io.files_from_extension`.
        
        Returns:
            :class:`PypeitSetup`: The instance of the class.
        """
        return cls.from_rawfiles(io.files_from_extension(root, extension=extension), spectrograph)

    @classmethod
    def from_rawfiles(cls, data_files:list, spectrograph:str, frametype=None):
        """ Instantiate the :class:`PypeItSetup` object by providing a list of raw files.

        Args:
            data_files (list): 
                List of raw files to be reduced.
            spectrograph (str): 
                The PypeIt name of the spectrograph used 
            frametype (:obj:`dict`, optional):
                A dictionary that associates the name of the file (just the fits
                file name without the full path) to a specific frame type (e.g.,
                arc, bias, etc.).  The file name and type are expected to be the
                key and value of the dictionary, respectively.  If None, this is
                determined by the :func:`get_frame_types` method.

        Returns:
            :class:`PypeItSetup`: The instance of the class.
        """

        # Configure me
        cfg_lines = ['[rdx]']
        cfg_lines += ['    spectrograph = {0}'.format(spectrograph)]

        # Instantiate
        return cls(data_files, cfg_lines=cfg_lines, frametype=frametype)

    @property
    def nfiles(self):
        """The number of files to reduce."""
        if self.fitstbl is None:
            msgs.warn('No fits files have been read!')
        return 0 if self.fitstbl is None else len(self.fitstbl)

    def __repr__(self):
        return '<{:s}: nfiles={:d}>'.format(self.__class__.__name__, self.nfiles)

    def build_fitstbl(self, strict=True):
        """
        Construct the table with metadata for the frames to reduce.

        Largely a wrapper for :func:`pypeit.core.load.create_fitstbl`.

        Args:
            strict (:obj:`bool`, optional):
                Function will fault if :func:`fits.getheader` fails to
                read the headers of any of the files in
                :attr:`file_list`.  Set to False to only report a
                warning and continue.

        Returns:
            `astropy.table.Table`_: Table with the metadata for each fits file
            to reduce.  Note this is different from :attr:`fitstbl`, which is a
            :class:`~pypeit.metadata.PypeItMetaData` object.
        """
        # Build and sort the table
        self.fitstbl = PypeItMetaData(self.spectrograph, par=self.par, 
                                      files=self.file_list,
                                      usrdata=self.usrdata, strict=strict)
        # Sort by the time
        if 'time' in self.fitstbl.keys():
            self.fitstbl.sort('time')

        # Return the table
        return self.fitstbl.table

    def get_frame_types(self, flag_unknown=False):
        """
        Include the frame types in the metadata table.

        This is mainly a wrapper for
        :func:`PypeItMetaData.get_frame_types`.

        .. warning::

            Because this merges the frame types with the existing
            :attr:`fitstbl` this should only be run once.

        Args:
            flag_unknown (:obj:`bool`, optional):
                Allow for frames to have unknown types instead of crashing.
                This should be True for initial setup and False otherwise.
                Passed to :func:`~pypeit.metadata.PypeItMetaData.get_frame_types`.
        """
        # Use PypeItMetaData methods to get the frame types
        self.fitstbl.get_frame_types(flag_unknown=flag_unknown, user=self.frametype)

    def run(self, setup_only=False, clean_config=True, groupings=True):
        """
        Perform the main setup operations.
        
        The code flow is as follows::
            - Build the metadata table from an input file_list (if it hasn't
              been already)
            - Remove frames (if requested using ``clean_config``) that cannot be
              associated with a configuration because of incomplete metadata or
              a configuration that PypeIt cannot reduce 
            - Type the files (bias, arc, etc.)
            - Match calibration files to the science files (if ``groupings`` is
              True)

        It is expected that a user will run this twice if they're being careful.
        Once with `setup_only=True` to confirm the images are properly typed and
        grouped together for calibration.  A second time to do the actual setup
        that groups frames by instrument configuration, preparing the user to
        proceed with the reductions using the ``run_pypeit`` script.

        Args:
            setup_only (:obj:`bool`, optional):
                Only this setup will be performed. ``PypeIt`` is
                expected to execute in a way that ends after this
                class is fully instantiated such that the user can
                inspect the results before proceeding. This has the
                effect of providing more output describing the
                success of the setup and how to proceed, and provides
                warnings (instead of errors) for issues that may
                cause the reduction itself to fail.
                If True, this also allows for bad headers.
            clean_config (:obj:`bool`, optional):
                Remove files with metadata that indicate an
                instrument configuration that ``PypeIt`` cannot
                reduce. See
                :func:`~pypeit.spectrographs.spectrograph.Spectrograph.valid_configuration_values`.
            groupings (:obj:`bool`, optional):
                Group frames into instrument configurations and calibration
                sets, and add the default combination-group columns.

        Returns:
            :obj:`tuple`: Returns, respectively, the
            :class:`~pypeit.par.pypeitpar.PypeItPar` object with the reduction
            parameters, the
            :class:`~pypeit.spectrographs.spectrograph.Spectrograph` object with
            the spectrograph instance, and the 
            :func:`~pypeit.metadata.PypeItMetaData` object with the an file
            metadata.
        """
        # Build the minimal metadata table if it doesn't exist already
        if self.fitstbl is None:
            self.build_fitstbl(strict=not setup_only)

        # Remove frames that have invalid values for
        # configuration-defining metadata
        if clean_config:
            self.fitstbl.clean_configurations()
            if len(self.fitstbl) == 0:
                msgs.error('Cleaning the configurations removed all the files!  Rerun '
                           'pypeit_setup with the --keep_bad_frames option.')

        # Determine the type of each frame.
        self.get_frame_types(flag_unknown=setup_only)

        if groupings:
            # Determine the configurations and assign each frame to the
            # specified configuration
            self.fitstbl.set_configurations()

            # Assign frames to calibration groups
            self.fitstbl.set_calibration_groups()

            # Set default comb_id
            self.fitstbl.set_combination_groups()

        return self.par, self.spectrograph, self.fitstbl

    def remove_table_rows(self, rows, regroup=False):
        """
        Remove rows from :attr:`fitstbl`.

        This is a wrapper for
        :func:`~pypeit.metadata.PypeItMetaData.remove_rows` that propagates the
        files removed from the fits table to the other attributes of the class:
        :attr:`file_list`, :attr:`frametype`, and :attr:`usrdata`.

        This *directly* modifies the attributes of the instance.

        Args:
            rows (:obj:`int`, array-like):
                One or more rows that should be *removed* from the datatable.
                This is passed directly to `astropy.table.Table.remove_rows`_;
                see astropy documentation to confirm allowed types.
            regroup (:obj:`bool`, optional):
                If True, reset the setup/configuration, calibration, and
                combination groups.
        """
        # Get the names of the files to be removed
        removed_files = self.fitstbl.frame_paths(rows)
        # Remove 'em
        self.fitstbl.remove_rows(rows, regroup=regroup)
        # Remove the files from the file list
        self.file_list = [f for f in self.file_list 
                            if Path(f).resolve().name in self.fitstbl['filename']]
        # Remove the files from the frametype
        if self.frametype is not None:
            self.frametype = {k : v for k,v in self.frametype.items()
                                if Path(k).resolve().name in self.fitstbl['filename']}
        # Remove the files from the user data
        if self.usrdata is not None:
            keep = [i for i in range(len(self.usrdata)) 
                        if self.usrdata['filename'][i] in self.fitstbl['filename']]
            self.usrdata = self.usrdata[keep]

    def generate_ql_calib_pypeit_files(self, output_path:str, 
                                       det:str=None, 
                                       configs:str='all',
                                       bkg_redux:bool=False,
                                       overwrite:bool=False):
        """
        Generate the PypeIt files for the calibrations for quicklook purposes.

        Args:
            output_path (str): 
                Output path for the PypeIt files
            det (str, optional):
                Detector/mosaic.
            configs (str, optional):
                Which configurations to generate.
            bkg_redux (bool, optional):
                Setup for A-B subtraction.
            overwrite (bool, optional):
                Overwrite existing files.

        Returns:
            list: List of calib PypeIt files
        """
        # Grab setups
        setups = self.fitstbl.get_configuration_names()

        # Restrict on detector -- May remove this
        self.user_cfg = ['[rdx]', f'spectrograph = {self.spectrograph.name}']
        if det is not None:
            self.user_cfg += [f'detnum = {det}']
        self.user_cfg += ['quicklook = True']

        # Set basic image processing flags (some of these may be the default,
        # but they're set here just to make sure)
        self.user_cfg += ['[baseprocess]',
                          'use_specillum = False',
                          'use_pattern = False']
        if not any(self.fitstbl.find_frames('bias')):
            # Turn-off use of bias by default
            self.user_cfg += ['use_biasimage = False']
        if not any(self.fitstbl.find_frames('dark')):
            # Turn-off use of bias by default
            self.user_cfg += ['use_darkimage = False']

        if not any(self.fitstbl.find_frames('pixelflat')):
            self.user_cfg += ['use_pixelflat = False']
            self.user_cfg += ['use_illumflat = False']
        elif not any(self.fitstbl.find_frames('illumflat')):
            self.user_cfg += ['use_illumflat = False']

        # Turn off the fine correction for the slit illumination
        self.user_cfg += ['[calibrations]', '[[flatfield]]', 'slit_illum_finecorr = False']

        # Write the PypeIt files
        # TODO: Exclude science/standard files from file?
        pypeit_files = self.fitstbl.write_pypeit(
            output_path=output_path, 
            cfg_lines=self.user_cfg, 
            write_bkg_pairs=bkg_redux, 
            configs=configs)

        # Rename name the pypeit files so that they're specific to the
        # calibrations
        calib_pypeit_files = []
        for pypeit_file, setup in zip(pypeit_files, setups):
            calib_pypeit_file = pypeit_file.replace(f'_{setup}.pypeit', f'_calib_{setup}.pypeit')
            if not os.path.isfile(calib_pypeit_file) or overwrite:
                # The file doesn't exist or we want to overwrite it, so move the
                # written pypeit file to its new name
                os.rename(pypeit_file, calib_pypeit_file)
            else:
                # The calibration file already exists or we don't want to
                # overwrite it, so remove the pypeit file that was written.
                os.remove(pypeit_file)
            # At this point, the 'calib_pypeit_file' should always exist.
            calib_pypeit_files.append(calib_pypeit_file)

        return calib_pypeit_files

