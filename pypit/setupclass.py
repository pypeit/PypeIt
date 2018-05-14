#  Class for organizing PYPIT setup
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np

from importlib import reload

from astropy.table import hstack, Table

from linetools import utils as ltu

from pypit import msgs
from pypit import ardebug as debugger
from pypit import arload
from pypit import arparse
from pypit import arsort
from pypit import arsetup

# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)


class SetupClass(object):
    """Class to handle setup

    Parameters
    ----------

    Attributes
    ----------
    """
    def __init__(self, settings_argflag, settings_spect, fitstbl=None):

        # Required parameters
        self.settings_argflag = settings_argflag
        self.settings_spect = settings_spect

        # Other parameters
        self.fitstbl = fitstbl

        # Outputs
        self.setup_dict = {}

        # Attributes
        self.ftypes = arsort.ftype_list

    def build_fitstbl(self, file_list):
        self.fitstbl, updates = arload.load_headers(file_list, self.settings_spect,
                                                    self.settings_argflag)
        self.fitstbl.sort('time')
        return self.fitstbl

    def build_group_dict(self):
        reload(arsetup)
        #
        all_sci_idx = self.fitstbl['sci_idx'].data[self.fitstbl['science']]
        self.group_dict = arsetup.new_build_group_dict(self.fitstbl, self.setupIDs, all_sci_idx)

        # Write .sorted file
        if len(self.group_dict) > 0:
            if len(self.settings_argflag['run']['redname']) == 0: # Stop gap
                group_file = 'tmp.sorted'
            else:
                group_file = self.settings_argflag['run']['redname'].replace('.pypit', '.sorted')
            arsetup.write_sorted(group_file, self.fitstbl, self.group_dict, self.setup_dict)
            msgs.info("Wrote group dict to {:s}".format(group_file))
        else:
            msgs.warn("No group dict entries and therefore no .sorted file")

    def build_setup_dict(self):
        #reload(arsetup)

        # Setup?
        if self.settings_argflag['run']['setup']:
            skip_cset = True
        else:
            skip_cset = False

        # Run with masters?
        if self.settings_argflag['reduce']['masters']['force']:
            # Check that setup was input
            if len(self.settings_argflag['reduce']['masters']['setup']) == 0:
                msgs.error("You need to specify the following parameter in your PYPIT file:"+msgs.newline()+"reduce masters setup")
            # Generate a dummy setup_dict
            self.setup_dict = arsetup.new_dummy_setup_dict(self.fitstbl)
            # Return
            return self.setup_dict

        # Run through the setups to fill setup_dict
        self.setupIDs = []
        all_sci_idx = self.fitstbl['sci_idx'].data[self.fitstbl['science']]
        for sc in all_sci_idx:
            for kk in range(self.settings_spect['mosaic']['ndet']):
                det = kk+1
                try:
                    cname = self.settings_argflag['setup']['name']
                except KeyError:
                    cname = None
                # Amplifiers
                dnum = arparse.get_dnum(det)
                namp = self.settings_spect[dnum]["numamplifiers"]
                # Run
                setupID = arsetup.new_instr_setup(sc, det, self.fitstbl, self.setup_dict, namp,
                                                  skip_cset=skip_cset, config_name=cname)
                # Only save the first detector for run setup
                if kk == 0:
                    self.setupIDs.append(setupID)

    def match_to_science(self):
        self.fitstbl = arsort.match_to_science(self.fitstbl,
                                         self.settings_spect,
                                         self.settings_argflag)
        return self.fitstbl

    def type_data(self):
        if len(arparse.ftdict) > 0:  # This is ugly!
            ftdict = arparse.ftdict
        else:
            ftdict = None
        self.filetypes = arsort.type_data(self.fitstbl, self.settings_spect,
                                     self.settings_argflag,
                                     ftdict=ftdict)

        # hstack me -- Might over-write self.fitstbl here
        msgs.info("Adding file type information to the fitstbl")
        self.fitstbl = hstack([self.fitstbl, self.filetypes])

        # Return
        return self.filetypes

    def load_fitstbl(self, fits_file):
        self.fitstbl = Table.read(fits_file)
        msgs.info("Loaded fitstbl from {:s}".format(fits_file))

    def write_fitstbl(self, outfile, overwrite=True):
        """
        Write to FITS

        Parameters
        ----------
        outfile : str
        overwrite : bool (optional)
        """
        self.fitstbl.write(outfile, overwrite=overwrite)

    def run(self, file_list=None):
        """ Main driver for file typing and sorting

          Code flow:

        Parameters
        ----------

        Returns
        -------
        fitstbl : Table
        setup_dict : Table
        """
        # Build fitstbl
        if self.fitstbl is None:
            _ = self.build_fitstbl(file_list)

        # Write?

        # File typing
        _ = self.type_data()

        # Match calibs to science
        _ = self.match_to_science()

        # Make dirs
        #arsort.make_dirs()

        # Setup dict
        _ = self.build_setup_dict()

        # .sorted Table (on pypit_setup only)
        if self.settings_argflag['run']['setup']:  # Collate all matching files
            _ = self.build_group_dict()

        # Write setup -- only if not present
        #setup_file, nexist = arsetup.get_setup_file(spectrograph=self.settings.argflag['run']['spectrograph'])

        # Write calib file (not in setup mode) or setup file (in setup mode)
        if not self.settings_argflag['run']['setup']:
            if len(self.settings_argflag['run']['redname']) == 0: # Stop gap
                calib_file = 'tmp.calib'
            else:
                calib_file = self.settings_argflag['run']['redname'].replace('.pypit', '.calib')
            arsetup.write_calib(calib_file, self.setup_dict)
        else:
            if len(self.settings_argflag['run']['redname']) == 0: # Stop gap
                setup_file = 'tmp.setups'
            else:
                setup_file = self.settings_argflag['run']['redname'].replace('.pypit', '.setups')
            arsetup.write_setup(self.setup_dict, setup_file=setup_file)

        # Finish (depends on PYPIT run mode)
        if self.settings_argflag['run']['calcheck']:
            msgs.info("Inspect the .calib file: {:s}".format(calib_file))
            msgs.info("*********************************************************")
            msgs.info("Calibration check complete and successful!")
            msgs.info("Set 'run calcheck False' to continue with data reduction")
            msgs.info("*********************************************************")
            # Instrument specific (might push into a separate file)
            if self.settings_argflag['run']['spectrograph'] in ['keck_lris_blue']:
                if self.settings.argflag['reduce']['flatfield']['useframe'] in ['pixelflat']:
                    msgs.warn("We recommend a slitless flat for your instrument.")
            return 'calcheck', None, None
        elif self.settings_argflag['run']['setup']:
            for idx in np.where(self.fitstbl['failures'])[0]:
                msgs.warn("No Arc found: Skipping object {:s} with file {:s}".format(
                    self.fitstbl['target'][idx],self.fitstbl['filename'][idx]))
            msgs.info("Setup is complete.")
            msgs.info("Inspect the .setups file")
            return 'setup', None, None
        else:
            return 'run', self.fitstbl, self.setup_dict

    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: >'.format(self.__class__.__name__)
        return txt



