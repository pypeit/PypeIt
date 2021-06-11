"""
Class for organizing PypeIt setup

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import time
import os
import glob
import inspect
import datetime

from IPython import embed

import numpy as np

from astropy.table import hstack, Table

from pypeit import msgs
from pypeit.metadata import PypeItMetaData

from pypeit.par import PypeItPar
from pypeit.par.util import parse_pypeit_file, make_pypeit_file
from pypeit.spectrographs.util import load_spectrograph


# TODO: Instantiation should just automatically trigger the run
# method...
class PypeItSetup:
    """
    Prepare for a pypeit run.

    .. todo::
        - This is now mostly a wrapper for PypeItMetaData.  Should we
          remove this class, or merge PypeItSetup and PypeItMetaData.

    The main deliverables are the set of parameters used for pypeit's
    algorithms (:attr:`par`), an :obj:`astropy.table.Table` with the
    details of the files to be reduced (:attr:`fitstbl`), and a
    dictionary with the list of instrument setups.

    Args:
        file_list (list):
            A list of strings with the full path to each file to be
            reduced.
        frametype (:obj:`dict`, optional):
            A dictionary that associates the name of the file (just the
            fits file name without the full path) to a specific frame
            type (e.g., arc, bias, etc.).  If None, this is determined
            by the :func:`get_frame_types` method.
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
            An instance of the `Spectograph` class used throughout the
            reduction procedures.
        par (:class:`pypeit.par.pypeitpar.PypitPar`):
            An instance of the `PypitPar` class that provides the
            parameters to all the algorthms that pypeit uses to reduce
            the data.
        fitstbl (:class:`pypeit.metadata.PypeItMetaData`):
            A `Table` that provides the salient metadata for the fits
            files to be reduced.
        setup_dict (dict):
            The dictionary with the list of instrument setups.
        steps (list):
            The steps run to provide the pypeit setup.
    """
    def __init__(self, file_list, path=None, frametype=None, usrdata=None, setups=None,
                 cfg_lines=None, spectrograph_name=None, pypeit_file=None):

        # The provided list of files cannot be None
        if file_list is None or len(file_list) == 0:
            msgs.error('Must provide a list of files to be reduced!')

        # Save input
        self.file_list = file_list
        # TODO: It seems like path is never used...
        self.path = os.getcwd() if path is None else path
        self.frametype = frametype
        self.usrdata = usrdata
        self.setups = setups
        self.pypeit_file = pypeit_file
        self.user_cfg = cfg_lines

        # Determine the spectrograph name
        _spectrograph_name = spectrograph_name if cfg_lines is None \
                    else PypeItPar.from_cfg_lines(merge_with=cfg_lines)['rdx']['spectrograph']

        # Cannot proceed without spectrograph name
        if _spectrograph_name is None:
            msgs.error('Must provide spectrograph name directly or using configuration lines.')
       
        # Instantiate the spectrograph
        self.spectrograph = load_spectrograph(_spectrograph_name)#, ifile=file_list[0])

        # Get the spectrograph specific configuration to be merged with
        # the user modifications.
        spectrograph_cfg_lines = self.spectrograph.default_pypeit_par().to_config()

        # Instantiate the pypeit parameters.  The user input
        # configuration (cfg_lines) can be None.
        self.par = PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines, merge_with=cfg_lines)

        # Prepare internals for execution
        self.fitstbl = None
        self.setup_dict = None
        self.steps = []

    @classmethod
    def from_pypeit_file(cls, filename):
        """
        Instantiate the :class:`PypeitSetup` object using a pypeit file.

        Args:
            filename (str):
                Name of the pypeit file to read.  Pypit files have a
                specific set of valid formats. A description can be
                found :ref:`pypeit_file`.
        
        Returns:
            :class:`PypeitSetup`: The instance of the class.
        """
        cfg_lines, data_files, frametype, usrdata, setups = parse_pypeit_file(filename)
        return cls(data_files, frametype=frametype, usrdata=usrdata, setups=setups,
                   cfg_lines=cfg_lines, pypeit_file=filename)

    @classmethod
    def from_file_root(cls, root, spectrograph, extension='.fits', output_path=None):
        """
        Instantiate the :class:`PypeItSetup` object by providing a file
        root.
        
        Args:
            root (:obj:`str`):
                The root path to all the files for PypeIt to reduce.
                This should be everything up to the wild-card before the
                file extension to use to find the relevant files.  The
                root itself can have wild cards to read through multiple
                directories.
            spectrograph (:obj:`str`):
                The PypeIt name of the spectrograph used to take the
                observations.  This should be one of the available
                options in
                :func:`pypeit.spectrographs.available_spectrographs`.
            extension (:obj:`str`, optional):
                The extension common to all the fits files to reduce.
                Default is '.fits', meaning anything with `root*.fits*`
                will be be included.
            output_path (:obj:`str`, optional):
                Path to use for the output. If None, no output is written.
                Otherwise, the method first writes a vanilla ``PypeIt`` file
                and then uses :func:`from_pypeit_file` to instantiate the
                object. If the path doesn't yet exist, it is created.
        
        Returns:
            :class:`PypitSetup`: The instance of the class.
        """
        if output_path is None:
            # Find all the files
            dfname = os.path.join(root, '*{0}*'.format(extension)) \
                        if os.path.isdir(root) else '{0}*{1}*'.format(root, extension)
            # No pypeit file is written and we just construct the setup object
            # directly
            return cls(glob.glob(dfname), spectrograph_name=spectrograph)

        # Set the output directory
        if not os.path.isdir(output_path):
            os.makedirs(output_path)
        # Set the output file name
        date = str(datetime.date.today().strftime('%Y-%m-%d'))
        pypeit_file = os.path.join(output_path, '{0}_{1}.pypeit'.format(spectrograph, date))
        msgs.info('A vanilla pypeit file will be written to: {0}'.format(pypeit_file))
        
        # Generate the pypeit file
        cls.vanilla_pypeit_file(pypeit_file, root, spectrograph, extension=extension)

        # Now setup PypeIt using that file
        return cls.from_pypeit_file(pypeit_file)

    @staticmethod
    def vanilla_pypeit_file(pypeit_file, root, spectrograph, extension='.fits'):
        """
        Write a vanilla PypeIt file.

        Args:
            pypeit_file (str):
              Name of PypeIt file to be generated
            root (str):
            spectrograph (str):
              Name of spectrograph
            extension (str, optional):
              File extension

        Returns:

        """
        # Generate
        dfname = os.path.join(root, '*{0}*'.format(extension)) \
                    if os.path.isdir(root) else '{0}*{1}*'.format(root, extension)
        # configuration lines
        cfg_lines = ['[rdx]']
        cfg_lines += ['    spectrograph = {0}'.format(spectrograph)]
#        cfg_lines += ['    sortroot = {0}'.format(root)]
        make_pypeit_file(pypeit_file, spectrograph, [dfname], cfg_lines=cfg_lines, setup_mode=True)

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
            :obj:`astropy.table.Table`: Table with the metadata for each
            fits file to reduce.  Note this is different from
            :attr:`fitstbl` which is a :obj:`PypeItMetaData` object
        """
        # Build and sort the table
        self.fitstbl = PypeItMetaData(self.spectrograph, par=self.par, files=self.file_list,
                                      usrdata=self.usrdata, strict=strict)
        # Sort by the time
        if 'time' in self.fitstbl.keys():
            self.fitstbl.sort('time')

        # Add this to the completed steps
        self.steps.append(inspect.stack()[0][3])

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
                Allow for frames to have unknown types instead of
                crashing.  This should be True for initial setup and
                False otherwise.  Passed to get_frame_types()

        """
        # Use PypeItMetaData methods to get the frame types
        self.fitstbl.get_frame_types(flag_unknown=flag_unknown, user=self.frametype)
        # Include finished processing step
        self.steps.append(inspect.stack()[0][3])

    def run(self, setup_only=False, calibration_check=False, sort_dir=None, write_bkg_pairs=False,
            clean_config=True, groupings=True, obslog=False, write_files=True):
        """
        Once instantiated, this is the main method used to construct the
        object.
        
        The code flow is as follows::
            - Build the fitstbl from an input file_list (optional)
            - Type the files (bias, arc, etc.)
            - Match calibration files to the science files
            - Generate the setup_dict
                - Write group info to disk
                - Write calib info to disk (if main run)

        It is expected that a user will run this three times if they're
        being careful.  Once with `setup_only=True` to confirm the
        images are properly typed and grouped together for calibration.
        A second time with `calibration_check=True` to confirm the
        appropriate calibrations frames are available.  And a third time
        to do the actual setup before proceeding with the reductions.

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
            calibration_check (obj:`bool`, optional):
                Only check that the calibration frames are
                appropriately setup and exist on disk. ``PypeIt`` is
                expected to execute in a way that ends after this
                class is fully instantiated such that the user can
                inspect the results before proceeding.
            sort_dir (:obj:`str`, optional):
                The directory to put the '.sorted' file.
            write_bkg_pairs (:obj:`bool`, optional):
                Include columns with the (unassigned) background
                image pairs. This is a convenience so that users can
                more easily add/edit the background pair ID numbers.
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
            :class:`~pypeit.par.pypeitpar.PypeItPar` object with the
            reduction parameters, the
            :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
            object with the spectrograph instance, and an
            `astropy.table.Table`_ with the frame metadata. If
            ``setup_only`` is True, these are all returned as None
            values.
        """
        # Build the minimal metadata table if it doesn't exist already
        if self.fitstbl is None:
            self.build_fitstbl(strict=not setup_only)

        # Remove frames that have invalid values for
        # configuration-defining metadata
        if clean_config:
            self.fitstbl.clean_configurations()

        # Determine the type of each frame.
        self.get_frame_types(flag_unknown=setup_only or calibration_check)

        if groupings:
            # Determine the configurations and assign each frame to the
            # specified configuration
            self.fitstbl.set_configurations(self.fitstbl.unique_configurations())

            # Assign frames to calibration groups
            self.fitstbl.set_calibration_groups()

            # Set default comb_id
            self.fitstbl.set_combination_groups()

#        # TODO: Are we planning to do this?
#        # Assign science IDs based on the calibrations groups (to be
#        # deprecated)
#        self.fitstbl['failures'] = False                    # TODO: placeholder

        # Write the output files
        if write_files:
            self.write(output_path=sort_dir, setup_only=setup_only,
                       calibration_check=calibration_check, write_bkg_pairs=write_bkg_pairs,
                       obslog=obslog)

        return (None, None, None) if setup_only else (self.par, self.spectrograph, self.fitstbl)

    def write(self, output_path=None, setup_only=False, calibration_check=False, 
              write_bkg_pairs=False, obslog=False):
        """
        Write the set of pypeit setup files.

        Args:
            output_path (:obj:`str`, optional):
                Directory for the output files.
            setup_only (:obj:`bool`, optional):
                Only this setup will be performed. ``PypeIt`` is
                expected to execute in a way that ends after this
                class is fully instantiated such that the user can
                inspect the results before proceeding. This has the
                effect of providing more output describing the
                success of the setup and how to proceed, and provides
                warnings (instead of errors) for issues that may
                cause the reduction itself to fail.
            calibration_check (obj:`bool`, optional):
                Only check that the calibration frames are
                appropriately setup and exist on disk. ``PypeIt`` is
                expected to execute in a way that ends after this
                class is fully instantiated such that the user can
                inspect the results before proceeding.
            write_bkg_pairs (:obj:`bool`, optional):
                Include columns with the (unassigned) background
                image pairs. This is a convenience so that users can
                more easily add/edit the background pair ID numbers.
            clean_config (:obj:`bool`, optional):
                Remove files with metadata that indicate an
                instrument configuration that ``PypeIt`` cannot
                reduce. See
                :func:`~pypeit.spectrographs.spectrograph.Spectrograph.valid_configuration_values`.
            groupings (:obj:`bool`, optional):
                Group frames into instrument configurations and calibration
                sets, and add the default combination-group columns.
        """
        pypeit_file = self.spectrograph.name + '.pypeit' \
                            if self.pypeit_file is None or len(self.pypeit_file) == 0 \
                            else self.pypeit_file
        _output_path = os.path.split(os.path.abspath(pypeit_file))[0] if output_path is None \
                            else os.path.abspath(output_path)

        if setup_only:
            # Collate all matching files and write .sorted Table (on pypeit_setup only)
            sorted_file = pypeit_file.replace('.pypeit', '.sorted')
            sorted_file = os.path.join(_output_path, os.path.split(sorted_file)[1])
            self.fitstbl.write_sorted(sorted_file, write_bkg_pairs=write_bkg_pairs)
            msgs.info("Wrote sorted file data to {:s}".format(sorted_file))
        else:
            # Write the calib file
            # This is currently needed for QA
            calib_file = pypeit_file.replace('.pypeit', '.calib')
            self.fitstbl.write_calib(os.path.join(_output_path, os.path.split(calib_file)[1]))

        if obslog:
            log_file = pypeit_file.replace('.pypeit', '.obslog')
            header = ['Auto-generated PypeIt Observing Log',
                      '{0}'.format(time.strftime("%a %d %b %Y %H:%M:%S", time.localtime()))]
            self.fitstbl.write(output=log_file, columns='pypeit', sort_col='mjd',
                               overwrite=True, header=header)

        # Finish (depends on PypeIt run mode)
        # TODO: Do we need this?
        if calibration_check:
            msgs.info("Inspect the .calib file: {:s}".format(calib_file))
            msgs.info("*********************************************************")
            msgs.info("Calibration check complete and successful!")
            msgs.info("*********************************************************")

        if setup_only:
#            for idx in np.where(self.fitstbl['failures'])[0]:
#                msgs.warn("No Arc found: Skipping object {:s} with file {:s}".format(
#                            self.fitstbl['target'][idx],self.fitstbl['filename'][idx]))
            msgs.info("Setup is complete.")
            msgs.info("Inspect the .sorted file")


