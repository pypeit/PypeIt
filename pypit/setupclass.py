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
    def __init__(self, settings, fitstbl=None):

        # Required parameters
        self.settings = settings

        # Other parameters
        self.fitstbl = fitstbl

        # Outputs
        self.setup_dict = {}

        # Attributes
        self.ftypes = arsort.ftype_list

    def build_fitstbl(self, file_list):
        self.fitstbl, updates = arload.load_headers(file_list, self.settings.spect,
                                                    self.settings.argflag)
        self.fitstbl.sort('time')
        return self.fitstbl

    def build_group_dict(self):
        reload(arsetup)
        #
        all_sci_idx = self.fitstbl['sci_idx'].data[self.fitstbl['science']]
        self.group_dict = arsetup.new_build_group_dict(self.fitstbl, self.setupIDs, all_sci_idx)

        # Write .sorted file
        if len(self.group_dict) > 0:
            if len(self.settings.argflag['run']['redname']) == 0: # Stop gap
                group_file = 'tmp.sorted'
            else:
                group_file = self.settings.argflag['run']['redname'].replace('.pypit', '.sorted')
            arsetup.write_sorted(group_file, self.fitstbl, self.group_dict, self.setup_dict)
            msgs.info("Wrote group dict to {:s}".format(group_file))
        else:
            msgs.warn("No group dict entries and therefore no .sorted file")

    def build_setup_dict(self):
        #reload(arsetup)

        # Setup?
        if self.settings.argflag['run']['setup']:
            skip_cset = True
        else:
            skip_cset = False

        # Run with masters?
        if self.settings.argflag['reduce']['masters']['force']:
            # Check that setup was input
            if len(self.settings.argflag['reduce']['masters']['setup']) == 0:
                msgs.error("You need to specify the following parameter in your PYPIT file:"+msgs.newline()+"reduce masters setup")
            # Generate a dummy setup_dict
            self.setup_dict = arsetup.new_dummy_setup_dict(self.fitstbl)
            # Return
            return self.setup_dict

        # Run through the setups to fill setup_dict
        self.setupIDs = []
        all_sci_idx = self.fitstbl['sci_idx'].data[self.fitstbl['science']]
        for sc in all_sci_idx:
            for kk in range(self.settings.spect['mosaic']['ndet']):
                det = kk+1
                try:
                    cname = self.settings.argflag['setup']['name']
                except KeyError:
                    cname = None
                # Amplifiers
                dnum = self.settings.get_dnum(det)
                namp = self.settings.spect[dnum]["numamplifiers"]
                # Run
                setupID = arsetup.new_instr_setup(sc, det, self.fitstbl, self.setup_dict, namp,
                                                  skip_cset=skip_cset, config_name=cname)
                # Only save the first detector for run setup
                if kk == 0:
                    self.setupIDs.append(setupID)

    def match_to_science(self):
        self.fitstbl = arsort.match_to_science(self.fitstbl,
                                         self.settings.spect,
                                         self.settings.argflag)
        return self.fitstbl

    def type_data(self):
        if len(self.settings.ftdict) > 0:
            ftdict = self.settings.ftdict
        else:
            ftdict = None
        self.filetypes = arsort.type_data(self.fitstbl, self.settings.spect,
                                     self.settings.argflag,
                                     ftdict=ftdict)

        # hstack me -- Might over-write self.fitstbl here
        msgs.info("Adding file type information to the fitstbl")
        self.fitstbl = hstack([self.fitstbl, self.filetypes])

        # Return
        return self.filetypes

    def load_fitstbl(self, fits_file):
        self.fitstbl = Table.read(fits_file)
        msgs.info("Loaded fitstbl from {:s}".format(fits_file))

    def write(self, outfile, overwrite=True):
        """
        Write to FITS

        Parameters
        ----------
        outfile : str
        overwrite : bool (optional)
        """
        self.fitstbl.write(outfile, overwrite=overwrite)

    def run(self, file_list):
        """ Main driver for file typing and sorting

          Code flow:

        Parameters
        ----------

        Returns
        -------
        fitstbl : Table
        setup_dict : Table
        """
        # Build
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
        if self.settings.argflag['run']['setup']:  # Collate all matching files
            _ = self.build_group_dict()

        # Return
        return self.fitstbl, self.setup_dict

    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: >'.format(self.__class__.__name__)
        return txt



