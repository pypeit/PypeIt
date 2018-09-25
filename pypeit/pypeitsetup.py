#  Class for organizing PYPIT setup
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np

#from importlib import reload

from astropy.table import hstack, Table

from pypeit import msgs
from pypeit.core import pypsetup
from pypeit.metadata import PypeItMetaData

from pypeit.par import PypeItPar
from pypeit.par.util import parse_pypeit_file
from pypeit.spectrographs.util import load_spectrograph

from pypeit import debugger

class PypeItSetup(object):
    """
    Prepare for a pypeit run.

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
            by the :func:`type_data` method.
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
            given spectrogaph, you can set `cfg_lines = None`, but you
            then *must* provide the `spectrograph_name`.
        spectrograph_name (:obj:`str`, optional):
            If not providing a list of configuration lines
            (`cfg_lines`), this sets the spectrograph.  The spectrograph
            defined in `cfg_lines` takes precedent over anything
            provided by this argument.
        pypeit_file (:obj:`str`, optional):
            The name of the pypeit file used to instantiate the
            reduction.  This can be None, and will lead to default names
            for output files (TODO: Give list).  Setting `pypeit_file`
            here *only sets the name of the file*.  To instantiate a
            `PypitSetup` object directly from a pypeit file (i.e. by
            reading the file), use the :func:`from_pypeit_file` method;
            i.e.::
                
                setup = PypitSetup.from_pypeit_file('myfile.pypeit')

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
        fitstbl (:class:`astropy.table.Table`):
            A `Table` that provides the salient metadata for the fits
            files to be reduced.
        filetypeflags(:class:`astropy.table.Table`):
            A `Table` that flags the frame types of each fits file.
            TODO: Is it necessary to keep this?
        setup_dict (dict):
            The dictionary with the list of instrument setups.
        steps (list):
            The steps run to provide the pypeit setup.

    .. _configobj: http://configobj.readthedocs.io/en/latest/
    """
    def __init__(self, file_list, frametype=None, setups=None, cfg_lines=None,
                 spectrograph_name=None, pypeit_file=None):

        # The provided list of files cannot be None
        if file_list is None or len(file_list) == 0:
            msgs.error('Must provide a list of files to be reduced!')

        # Save input
        self.file_list = file_list
        self.frametype = frametype
        self.setups = setups
        self.pypeit_file = pypeit_file

        # Determine the spectrograph name
        _spectrograph_name = spectrograph_name if cfg_lines is None \
                    else PypeItPar.from_cfg_lines(merge_with=cfg_lines)['rdx']['spectrograph']

        # Cannot proceed without spectrograph name
        if _spectrograph_name is None:
            msgs.error('Must provide spectrograph name directly or using configuration lines.')
       
        # Instantiate the spectrograph
        self.spectrograph = load_spectrograph(_spectrograph_name)

        # Get the spectrograph specific configuration to be merged with
        # the user modifications.
        spectrograph_cfg_lines = self.spectrograph.default_pypeit_par().to_config()

        # Instantiate the pypeit parameters.  The user input
        # configuration (cfg_lines) can be None.
        self.par = PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines, merge_with=cfg_lines)

        # Prepare internals for execution
        self.fitstbl = None
        self.filetypeflags = None
        self.setup_dict = None
        self.steps = []

    @classmethod
    def from_pypeit_file(cls, filename):
        """
        Instantiate the :class:`PypitSetup` object using a pypeit file.

        Args:
            filename (str):
                Name of the pypeit file to read.  Pypit files have a
                specific set of valid formats. A description can be
                found `here`_ (include doc link).
        
        Returns:
            :class:`PypitSetup`: The instance of the class.
        """
        cfg_lines, data_files, frametype, setups = parse_pypeit_file(filename)
        return cls(data_files, frametype=frametype, setups=setups, cfg_lines=cfg_lines,
                   pypeit_file=filename)

    @property
    def nfiles(self):
        """The number of files to reduce."""
        return 0 if self.fitstbl is None else len(self.fitstbl)

    def __repr__(self):
        txt = '<{:s}: nfiles={:d}>'.format(self.__class__.__name__,
                                           self.nfiles)
        return txt

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
            fits file to reduce.
        """
        # Build and sort the table
        self.fitstbl = PypeItMetaData(self.spectrograph, par=self.par, file_list=self.file_list,
                                      strict=strict)
        # Sort by the time
        if 'time' in self.fitstbl.keys():
            self.fitstbl.sort('time')

        # Add this to the completed steps
        self.steps.append(inspect.stack()[0][3])

        # Return the table
        return self.fitstbl.table

    def build_group_dict(self, pypeit_file=None):
        """
        Builds a group dict and writes to disk
          This may be Deprecated (if the .sorted files are deemed too unintersting)

        Returns
        -------
        group_dict : dict
          Dict describing the various setups
        """
        #
        all_sci_idx = np.where(self.fitstbl.find_frames('science'))[0]
        all_sci_ID = self.fitstbl['sci_ID'][all_sci_idx]
        self.group_dict = pypsetup.build_group_dict(self.fitstbl, self.setupIDs, all_sci_idx,
                                                    all_sci_ID)

        # TODO: Move this to a method that writes the sorted file
        # Write .sorted file
        if len(self.group_dict) > 0:
            group_file = self.spectrograph.spectrograph + '.sorted' \
                                if pypeit_file is None or len(pypeit_file) == 0 \
                                else pypeit_file.replace('.pypeit', '.sorted')
            pypsetup.write_sorted(group_file, self.fitstbl, self.group_dict, self.setup_dict)
            msgs.info("Wrote group dict to {:s}".format(group_file))
        else:
            msgs.warn("No group dict entries and therefore no .sorted file")

        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.group_dict


    def build_setup_dict(self, setup_only=False):
        """
        Generate the setup_dict
          Mainly a Wrapper to new_instr_setup

        Returns
        -------
        setup_dict :

        """
        # Run with masters?
        if self.par['calibrations']['masters'] == 'force':
            print('forcing!')
            # Generate a dummy setup_dict
            self.setup_dict = pypsetup.dummy_setup_dict(self.fitstbl.find_frame_files('science'),
                                                        self.par['calibrations']['setup'])
            # Step
            self.steps.append(inspect.stack()[0][3])
            # Return
            return self.setup_dict

        # Run through the setups to fill setup_dict
        self.setupIDs = []
        all_sci_ID = self.fitstbl['sci_ID'][self.fitstbl.find_frames('science')]
        for sc in all_sci_ID:
            for kk in range(len(self.spectrograph.detector)):
                cname = None if self.par['calibrations']['setup'] is None \
                                    else self.par['calibrations']['setup'][0]
                # Amplifiers
                namp = self.spectrograph.detector[kk]['numamplifiers']
                # Run
                setupID, self.setup_dict = pypsetup.instr_setup(sc, kk+1, self.fitstbl,
                                                                setup_dict=self.setup_dict,
                                                                skip_cset=setup_only,
                                                                config_name=cname)
                # Only save the first detector for run setup
                if kk == 0:
                    self.setupIDs.append(setupID)
        # Step
        self.steps.append(inspect.stack()[0][3])
        return self.setup_dict

    def match_ABBA(self):
        """
          Matches science frames to their partner A/B frame
          Mainly a wrapper to arsort.match_ABBA()

        Returns
        -------
        self.fitstbl -- Updated with 'AB_frame' column

        """
        self.fitstbl.match_ABBA()
        self.steps.append(inspect.stack()[0][3])
        return self.fitstbl

    def match_to_science(self, setup_only=False):
        """
          Matches calibration frames to the Science
          Mainly a wrapper to arsort.match_to_science()

        Returns
        -------
        self.fitstbl -- Updated with 'sci_ID' and 'failures' columns

        """
        self.fitstbl.match_to_science(self.par['calibrations'], self.par['rdx']['calwin'],
                                      setup=setup_only)
        self.steps.append(inspect.stack()[0][3])
        return self.fitstbl

    def get_frame_types(self, flag_unknown=False, use_header_id=False):
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
                False otherwise.

        Returns
        -------
        self.filetypeflags
        """
        # Use PypeItMetaData methods to get the frame types
        self.filetypeflags = self.fitstbl.get_frame_types(flag_unknown=flag_unknown,
                                                          user=self.frametype,
                                                          useIDname=use_header_id)
        # Include finished processing step
        self.steps.append(inspect.stack()[0][3])

    def load_fitstbl(self, fits_file):
        """
          Load the fitstbl from disk (a binary FITS table)

        Parameters
        ----------
        fits_file : str

        Returns
        -------
        self.fitstbl

        """
        self.fitstbl = PypeItMetaData(self.spectrograph, data=Table.read(fits_file))
#        # Need to convert bytestrings back to unicode
#        try:
#            self.fitstbl.convert_bytestring_to_unicode()
#        except:
#            pass
        msgs.info("Loaded fitstbl from {:s}".format(fits_file))
        return self.fitstbl

    def write_metadata(self, setup_only=False, sort_dir=None):
        """
        Write :attr:`fitstbl` to a file.

        The default file is determined as follows:
            - if both :attr:`pypeit_file` and `sort_dir` are None, the
              output file is the name of the spectrograph with a '.lst'
              extension.
            - if `setup_only` is true and the pypeit file is valid, the
              output file is the pypeit file name with the ".pypeit"
              extension replaced by ".lst".
            - if `setup_only` is false and `sort_dir` is not None, the
              output file is simply `sort_dir` with '.lst' appended.
            - if none of these conditions are met a ValueError is
              raised.

        .. todo::
            - Is the use of sort_dir right here?

        Args:
            setup_only (:obj:`bool`, optional):
                Code is being run for setup only, not to reduce the data.
            sort_dir (:obj:`str`, optional):
                The full root of the name for the metadata table
                ('.lst') file.

        Raises:
            ValueError:
                Raised if the output file cannot be determined from the
                input and attributes.
        """
        if self.pypeit_file is None and sort_dir is None:
            ofile = self.spectrograph.spectrograph + '.lst'
        elif setup_only and self.pypeit_file is not None:
            ofile = self.pypeit_file.replace('.pypeit', '.lst')
        elif not setup_only and sort_dir is not None:
            ofile = sort_dir + '.lst'
        else:
            raise ValueError('Could not determine name for output file.')

        self.fitstbl.write(ofile, columns=self.spectrograph.metadata_keys(),
                           format='ascii.fixed_width', overwrite=True)

    def run(self, setup_only=False, calibration_check=False, use_header_id=False, sort_dir=None):
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
                Only this setup will be performed.  Pypit is expected to
                execute in a way that ends after this class is fully
                instantiated such that the user can inspect the results
                before proceeding.  This has the effect of providing
                more output describing the success of the setup and how
                to proceed, and provides warnings (instead of errors)
                for issues that may cause the reduction itself to fail.
            calibration_check (obj:`bool`, optional):
                Only check that the calibration frames are appropriately
                setup and exist on disk.  Pypit is expected to execute
                in a way that ends after this class is fully
                instantiated such that the user can inspect the results
                before proceeding. 
            use_header_id (:obj:`bool`, optional):
                Allow setup to use the frame types drawn from single,
                instrument-specific header keywords set to `idname` in
                the metadata table (:attr:`fitstbl`).
            sort_dir (:obj:`str`, optional):
                The directory to put the '.sorted' file.

        Returns:
            :class:`pypeit.par.pypeitpar.PypitPar`,
            :class:`pypeit.spectrographs.spectrograph.Spectrograph`,
            :class:`astropy.table.Table`, dict: Returns the attributes
            :attr:`par`, :attr:`spectrograph`, :attr:`fitstbl`, and
            :attr:`setup_dict`.  If running with `setup_only` or
            `calibrations_check`, these are all returned as `None`
            values.
        """
        # Kludge
        pypeit_file = '' if self.pypeit_file is None else self.pypeit_file

        # Build fitstbl
        if self.fitstbl is None:
            self.build_fitstbl(strict=not setup_only)

        # File typing
        if 'frametype' not in self.fitstbl.keys():
            # Add the frame type if it isn't already in the table
            self.get_frame_types(flag_unknown=setup_only or calibration_check,
                                 use_header_id=use_header_id)

        # Match calibs to science
        self.match_to_science(setup_only=setup_only)

        if self.par['scienceimage'] is not None and self.par['scienceimage']['nodding']:
            self.match_ABBA()

        # Write metadata
        self.write_metadata(setup_only=setup_only, sort_dir=sort_dir)

        # Setup dict
        self.build_setup_dict(setup_only=setup_only)

        if setup_only:
            # Collate all matching files and write .sorted Table (on pypeit_setup only)
            self.build_group_dict(pypeit_file=pypeit_file)

            # Write the setup file
            setup_file = self.spectrograph.spectrograph + '.setups' \
                                if pypeit_file is None or len(pypeit_file) == 0 \
                                else pypeit_file.replace('.pypeit', '.setups')
            pypsetup.write_setup(self.setup_dict, setup_file)
        else:
            # Write the calib file
            calib_file = self.spectrograph.spectrograph + '.calib' \
                                if pypeit_file is None or len(pypeit_file) == 0 \
                                else pypeit_file.replace('.pypeit', '.calib')
            pypsetup.write_calib(calib_file, self.setup_dict)

        # Finish (depends on PypeIt run mode)
        if calibration_check:
            msgs.info("Inspect the .calib file: {:s}".format(calib_file))
            msgs.info("*********************************************************")
            msgs.info("Calibration check complete and successful!")
            msgs.info("Set 'run calcheck False' to continue with data reduction")
            msgs.info("*********************************************************")
            return None, None, None, None

        if setup_only:
            for idx in np.where(self.fitstbl['failures'])[0]:
                msgs.warn("No Arc found: Skipping object {:s} with file {:s}".format(
                    self.fitstbl['target'][idx],self.fitstbl['filename'][idx]))
            msgs.info("Setup is complete.")
            msgs.info("Inspect the .setups file")
            return None, None, None, None

        return self.par, self.spectrograph, self.fitstbl, self.setup_dict



