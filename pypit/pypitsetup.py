#  Class for organizing PYPIT setup
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np

#from importlib import reload

from astropy.table import hstack, Table

from pypit import msgs
from pypit import arload
from pypit import arparse
from pypit.core import arsort
from pypit.core import arsetup

from pypit.par import PypitPar
from pypit.par.util import parse_pypit_file
from pypit.spectrographs.util import load_spectrograph

from pypit import ardebug as debugger

class PypitSetup(object):
    """
    Prepare for a pypit run.

    The main deliverables are the set of parameters used for pypit's
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
            parameters for executing pypit.  These are the lines of a
            configuration file.  See the documentation for the
            `configobj`_ package.  One of the user-level inputs should
            be the spectrograph that provided the data to be reduced.
            One can get the list of spectrographs currently served by
            running::
                
                from pypit.spectrographs.util import valid_spectrographs
                print(valid_spectrographs())

            To use all the default parameters when reducing data from a
            given spectrogaph, you can set `cfg_lines = None`, but you
            then *must* provide the `spectrograph_name`.
        spectrograph_name (:obj:`str`, optional):
            If not providing a list of configuration lines
            (`cfg_lines`), this sets the spectrograph.  The spectrograph
            defined in `cfg_lines` takes precedent over anything
            provided by this argument.
        pypit_file (:obj:`str`, optional):
            The name of the pypit file used to instantiate the
            reduction.  This can be None, and will lead to default names
            for output files (TODO: Give list).  Setting `pypit_file`
            here *only sets the name of the file*.  To instantiate a
            `PypitSetup` object directly from a pypit file (i.e. by
            reading the file), use the :func:`from_pypit_file` method;
            i.e.::
                
                setup = PypitSetup.from_pypit_file('myfile.pypit')

    Attributes:
        file_list (list):
            See description of class argument.
        frametype (dict):
            See description of class argument.
        setups (list):
            See description of class argument.
        pypit_file (str):
            See description of class argument.
        spectrograph (:class:`pypit.spectrographs.spectrograph.Spectrograph`):
            An instance of the `Spectograph` class used throughout the
            reduction procedures.
        par (:class:`pypit.par.pypitpar.PypitPar`):
            An instance of the `PypitPar` class that provides the
            parameters to all the algorthms that pypit uses to reduce
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
            The steps run to provide the pypit setup.

    .. _configobj: http://configobj.readthedocs.io/en/latest/
    """
    def __init__(self, file_list, frametype=None, setups=None, cfg_lines=None,
                 spectrograph_name=None, pypit_file=None, par=None):

        # The provided list of files cannot be None
        if file_list is None or len(file_list) == 0:
            msgs.error('Must provide a list of files to be reduced!')

        # Save input
        self.file_list = file_list
        self.frametype = frametype
        self.setups = setups
        self.pypit_file = pypit_file

        # Determine the spectrograph name
        _spectrograph_name = spectrograph_name if cfg_lines is None \
                    else PypitPar.from_cfg_lines(merge_with=cfg_lines)['rdx']['spectrograph'] 

        # Cannot proceed without spectrograph name
        if _spectrograph_name is None:
            msgs.error('Must provide spectrograph name directly or using configuration lines.')
       
        # Instantiate the spectrograph
        self.spectrograph = load_spectrograph(_spectrograph_name)

        # Get the spectrograph specific configuration to be merged with
        # the user modifications.
        spectrograph_cfg_lines = self.spectrograph.default_pypit_par().to_config()

        # Instantiate the pypit parameters.  The user input
        # configuration (cfg_lines) can be None.
        if par is None:
            self.par = PypitPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines, merge_with=cfg_lines)
        else:
            self.par = par

        # Prepare internals for execution
        self.fitstbl = None
        self.filetypeflags = None
        self.setup_dict = {}
        self.steps = []

    @classmethod
    def from_pypit_file(cls, filename):
        """
        Instantiate the :class:`PypitSetup` object using a pypit file.

        Args:
            filename (str):
                Name of the pypit file to read.  Pypit files have a
                specific set of valid formats. A description can be
                found `here`_ (include doc link).
        
        Returns:
            :class:`PypitSetup`: The instance of the class.
        """
        cfg_lines, data_files, frametype, setups = parse_pypit_file(filename)
        return cls(data_files, frametype=frametype, setups=setups, cfg_lines=cfg_lines,
                   pypit_file=filename)

    @property
    def nfiles(self):
        """The number of files to reduce."""
        return 0 if self.fitstbl is None else len(self.fitstbl)

    def build_fitstbl(self, strict=True):
        """

        Parameters
        ----------
        file_list : list
          List of file names for generating fitstbl

        Returns
        -------
        fitstbl : Table

        """
        # Build and sort the table
        self.fitstbl = arload.load_headers(self.file_list, self.spectrograph, strict=strict)
        self.fitstbl.sort('time')
        # Step
        self.steps.append(inspect.stack()[0][3])
        return self.fitstbl

    def build_group_dict(self, pypit_file=None):
        """
        Builds a group dict and writes to disk
          This may be Deprecated (if the .sorted files are deemed too unintersting)

        Returns
        -------
        group_dict : dict
          Dict describing the various setups
        """
        #
        all_sci_idx = np.where(self.fitstbl['science'])[0]
        all_sci_ID = self.fitstbl['sci_ID'][self.fitstbl['science']]
        self.group_dict = arsetup.build_group_dict(self.fitstbl, self.setupIDs, all_sci_idx,
                                                   all_sci_ID)

        # TODO: Move this to a method that writes the sorted file
        # Write .sorted file
        if len(self.group_dict) > 0:
            group_file = 'tmp.sorted' if pypit_file is None or len(pypit_file) == 0 \
                                else pypit_file.replace('.pypit', '.sorted')
            arsetup.write_sorted(group_file, self.fitstbl, self.group_dict, self.setup_dict)
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
            print(self.par['calibrations']['masters'])
            # TODO: This is now checked when validating the parameter
            # set.  See CalibrationsPar.validate()
#            # Check that setup was input
#            if len(self.spectrograph.calib_par['setup']) == 0:
#                msgs.error("When forcing use of master frames, you need to specify the You need to specify the following parameter in your PYPIT file:" 
#                           + msgs.newline() + "reduce masters setup")
            # Generate a dummy setup_dict
            self.setup_dict = arsetup.dummy_setup_dict(self.fitstbl,
                                                       self.par['calibrations']['setup'])
            # Step
            self.steps.append(inspect.stack()[0][3])
            # Return
            return self.setup_dict

        # Run through the setups to fill setup_dict
        self.setupIDs = []
        all_sci_ID = self.fitstbl['sci_ID'].data[self.fitstbl['science']]
        for sc in all_sci_ID:
            for kk in range(len(self.spectrograph.detector)):
                cname = None if self.par['calibrations']['setup'] is None \
                                    else self.par['calibrations']['setup'][0]
                # Amplifiers
                namp = self.spectrograph.detector[kk]["numamplifiers"]
                # Run
                det = kk+1
                setupID = arsetup.instr_setup(sc, det, self.fitstbl, self.setup_dict, namp,
                                              skip_cset=setup_only, config_name=cname)
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
        self.fitstbl = arsort.match_ABBA(self.fitstbl)

        # Step
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
        self.fitstbl = arsort.match_to_science(self.par['calibrations'],
                                               self.spectrograph.get_match_criteria(),
                                               self.fitstbl, self.par['rdx']['calwin'],
                                               setup=setup_only,
                                               match_nods=self.par['skysubtract'] is not None \
                                                            and self.par['skysubtract']['nodding'])
        # Step
        self.steps.append(inspect.stack()[0][3])
        return self.fitstbl

    # TODO: This appends the data to fitstbl meaning that it should not
    # be run multiple times.  Make it a "private" function?
    def type_data(self, flag_unknown=False, use_header_frametype=False):
        """
          Perform image typing on the full set of input files
          Mainly a wrapper to arsort.type_data()

        The table (filetypeflags) returned is horizontally stacked
          onto the fitstbl.

        Parameters
        ----------
        flag_unknown: bool, optional
          Mark a frame as UNKNOWN instead of crashing out
          Required when doing initial setup

        Returns
        -------
        self.filetypeflags

        """
        # Allow for input file types from the PYPIT file
        self.filetypeflags = arsort.type_data(self.spectrograph, self.fitstbl,
                                              ftdict=self.frametype, flag_unknown=flag_unknown,
                                              useIDname=use_header_frametype)

        # hstack me -- Might over-write self.fitstbl here
        msgs.info("Adding file type information to the fitstbl")
        self.fitstbl = hstack([self.fitstbl, self.filetypeflags])

        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.filetypeflags

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
        self.fitstbl = Table.read(fits_file)
        # Need to convert bytestrings back to unicode
        try:
            self.fitstbl.convert_bytestring_to_unicode()
        except:
            pass
        msgs.info("Loaded fitstbl from {:s}".format(fits_file))
        return self.fitstbl

    def write_fitstbl(self, outfile=None, overwrite=True):
        """
        Write fitstbl to FITS

        Parameters
        ----------
        outfile : str
        overwrite : bool (optional)
        """
        if outfile is None:
            outfile = self.pypit_file.replace('.pypit', '.fits')
        self.fitstbl.write(outfile, overwrite=overwrite)

    def run(self, setup_only=False, calibration_check=False, use_header_frametype=False,
            sort_dir=None):
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
            setup_only (bool):
                Only this setup will be performed.  Pypit is expected to
                execute in a way that ends after this class is fully
                instantiated such that the user can inspect the results
                before proceeding.  This has the effect of providing
                more output describing the success of the setup and how
                to proceed, and provides warnings (instead of errors)
                for issues that may cause the reduction itself to fail.
            calibration_check (bool):
                Only check that the calibration frames are appropriately
                setup and exist on disk.  Pypit is expected to execute
                in a way that ends after this class is fully
                instantiated such that the user can inspect the results
                before proceeding. 
            use_header_frametype (bool):
                Allow setup to use the frame types drawn from the file
                headers using the instrument specific keywords.
            sort_dir (str):
                The directory to put the '.sorted' file.

        Returns:
            :class:`pypit.par.pypitpar.PypitPar`,
            :class:`pypit.spectrographs.spectrograph.Spectrograph`,
            :class:`astropy.table.Table`, dict: Returns the attributes
            :attr:`par`, :attr:`spectrograph`, :attr:`fitstbl`, and
            :attr:`setup_dict`.  If running with `setup_only` or
            `calibrations_check`, these are all returned as `None`
            values.
        """
        # Kludge
        pypit_file = '' if self.pypit_file is None else self.pypit_file

        # Build fitstbl
        if self.fitstbl is None:
            self.build_fitstbl(strict=not setup_only)

        # File typing
        self.type_data(flag_unknown=setup_only or calibration_check,
                       use_header_frametype=use_header_frametype)

        # Write?
        if sort_dir is not None:
            print('WRITING: {0}'.format(sort_dir))
            arsort.write_lst(self.fitstbl, self.spectrograph.header_keys(), pypit_file,
                             setup=setup_only, sort_dir=sort_dir)

        # Match calibs to science
        self.match_to_science(setup_only=setup_only)

        # Setup dict
        self.build_setup_dict(setup_only=setup_only)

        if setup_only:
            # Collate all matching files and write .sorted Table (on pypit_setup only)
            self.build_group_dict(pypit_file=pypit_file)

            # Write the setup file
            setup_file = 'tmp.setups' if pypit_file is None or len(pypit_file) == 0 \
                                else pypit_file.replace('.pypit', '.setups')
            arsetup.write_setup(self.setup_dict, setup_file=setup_file)
        else:
            # Write the calib file
            calib_file = 'tmp.calib' if pypit_file is None or len(pypit_file) == 0 \
                                else pypit_file.replace('.pypit', '.calib')
            arsetup.write_calib(calib_file, self.setup_dict)

        # Finish (depends on PYPIT run mode)
        if calibration_check:
            msgs.info("Inspect the .calib file: {:s}".format(calib_file))
            msgs.info("*********************************************************")
            msgs.info("Calibration check complete and successful!")
            msgs.info("Set 'run calcheck False' to continue with data reduction")
            msgs.info("*********************************************************")
            # Instrument specific (might push into a separate file)
            # TODO: Move to spectrograph class
#            if self.spectrograph.spectrograph in ['keck_lris_blue']:
#                if self.spectrograph.calib_par['flatfield']['useframe'] in ['pixelflat']:
#                    msgs.warn("We recommend a slitless flat for your instrument.")
            return None, None, None, None

        if setup_only:
            for idx in np.where(self.fitstbl['failures'])[0]:
                msgs.warn("No Arc found: Skipping object {:s} with file {:s}".format(
                    self.fitstbl['target'][idx],self.fitstbl['filename'][idx]))
            msgs.info("Setup is complete.")
            msgs.info("Inspect the .setups file")
            return None, None, None, None

        return self.par, self.spectrograph, self.fitstbl, self.setup_dict

    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: nfiles={:d}>'.format(self.__class__.__name__,
                                           self.nfiles)
        return txt



