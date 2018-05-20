# Module for guiding Slit/Order tracing
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import yaml

from importlib import reload



from pypit import msgs
from pypit import ardebug as debugger
from pypit.core import arflux
from pypit import arload
from pypit import arsave
from pypit import arutils
from pypit import masterframe


# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

frametype = 'sensfunc'

# Default settings, if any

# TODO - Remove this kludge
def kludge_settings(instr):
    from pypit import arparse as settings
    settings.dummy_settings()
    settings.argflag['run']['spectrograph'] = instr
    settings.argflag['run']['directory']['master'] = 'MF'
    #settings.argflag['reduce']['masters']['setup'] = 'C_01_aa'
    #
    # Load default spectrograph settings
    spect = settings.get_spect_class(('ARMLSD', instr, 'pypit'))  # '.'.join(redname.split('.')[:-1])))
    lines = spect.load_file(base=True)  # Base spectrograph settings
    spect.set_paramlist(lines)
    lines = spect.load_file()  # Instrument specific
    spect.set_paramlist(lines)
    # Kludge deeper
    for key in ['run', 'reduce']:
        settings.spect[key] = settings.argflag[key].copy()
    return settings.spect

class FluxSpec(masterframe.MasterFrame):
    """Class to guide fluxing

    Parameters
    ----------

    Attributes
    ----------
    """
    def __init__(self, std_spec1d_file=None, sci_spec1d_file=None, spectrograph=None, sens_file=None, setup=None):

        # Parameters
        self.std_spec1d_file = std_spec1d_file
        self.sci_spec1d_file = sci_spec1d_file
        self.sens_file = sens_file
        self.setup = setup


        # Attributes
        self.frametype = frametype

        # Main outputs
        self.sensfunc = None  # dict

        # Settings (for now)
        if spectrograph is not None:
            self.settings = kludge_settings(spectrograph)

        # Key Internals
        self.std = None         # Standard star spectrum (SpecObj object)
        self.std_specobjs = []
        self.std_header = None
        self.std_idx = None     # Nested indices for the std_specobjs list that corresponds to the star!
        self.sci_specobjs = []
        self.sci_header = None

        # Attributes
        self.steps = []

        # Load
        if self.std_spec1d_file is not None:
            self.std_specobjs, self.std_header = arload.load_specobj(self.std_spec1d_file)
            msgs.info("Loaded {:d} spectra from the spec1d standard star file: {:s}".format(len(self.std_specobjs),
                                                                                            self.std_spec1d_file))
        if self.std_spec1d_file is not None:
            self.sci_specobjs, self.sci_header = arload.load_specobj(self.sci_spec1d_file)
            msgs.info("Loaded spectra from the spec1d science file: {:s}".format(self.sci_spec1d_file))
        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, self.setup, self.settings)

    def find_standard(self):
        # Find brightest object in the exposures
        # Search over all slits (over all detectors), and all objects
        #
        self.std_idx = arflux.find_standard(self.std_specobjs)
        #
        self.std = self.std_specobjs[self.std_idx]
        # Step
        self.steps.append(inspect.stack()[0][3])
        return self.std

    def generate_sensfunc(self):
        if self.std is None:
            msgs.warn("You need to identify the Star first (with find_standard)") #load in the Standard spectra first")
            return None
        # Run
        self.sensfunc = arflux.generate_sensfunc(self.std, self.std_header['RA'],
                                            self.std_header['DEC'], self.std_header['AIRMASS'],
                                            self.std_header['EXPTIME'], self.settings)
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.sensfunc

    def flux_science(self):
        reload(arflux)
        for sci_obj in self.sci_specobjs:
            arflux.new_apply_sensfunc(sci_obj, self.sensfunc,
                                  self.sci_header['AIRMASS'],
                                  self.sci_header['EXPTIME'],
                                  self.settings)
        # Step
        self.steps.append(inspect.stack()[0][3])

    def _set_std_obj(self):
        pass

    def run(self):
        #
        if self.std_specobjs is None:
            msgs.warn("You need to load in the Standard spectra first!")
            return None

        # Find the star automatically?
        _ = self.find_standard()

        # Sensitivity
        _ = self.generate_sensfunc()

        # Apply to science
        if len(self.sci_specobjs) > 0:
            self.flux_sci()

    def save_master(self):
        self.sensfunc['steps'] = self.steps
        # yamlify
        ysens = arutils.yamlify(self.sensfunc)
        with open(self.ms_name, 'w') as yamlf:
            yamlf.write(yaml.dump(ysens))
        #
        msgs.info("Wrote sensfunc to MasterFrame: {:s}".format(self.ms_name))
        # Step
        self.steps.append(inspect.stack()[0][3])

    def write_science(self, outfile):
        if len(self.sci_specobjs) == 0:
            msgs.warn("No science spectra to write to disk!")
        #
        helio_dict = dict(refframe=self.sci_header['VEL-TYPE'], vel_correction=self.sci_header['VEL'])
        arsave.new_save_1d_spectra_fits(self.sci_specobjs, self.sci_header, self.settings, outfile,
                             helio_dict=helio_dict, clobber=True)
        # Step
        self.steps.append(inspect.stack()[0][3])


    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: >'.format(self.__class__.__name__)
        return txt



