# Module for guiding Slit/Order tracing
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import yaml

#from importlib import reload

try:
    basestring
except NameError:  # For Python 3
    basestring = str

from pypit import msgs
from pypit import ardebug as debugger
from pypit.core import arflux
from pypit import arload
from pypit import armasters
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
    std_spec1d_file : str
      Filename of the spec1d file containing the standard star spectrum
      One or more of these are used to generate the sensitivity function
    std_specobjs : list
      List of SpecObj objects for the standard star spectrum/spectra
      May be input instead of std_spec1d_file to generate the sensitivity function
    sci_spec1d_file : str
      Filename of a spec1d file to be fluxed
    spectrograph : str
      Name of the spectrograph, e.g. shane_kast_blue
      Used only to set settings for calls to the Class outside of PYPIT
    sens_file : str
      Filename of a sensitivity function file to be input
    setup : str
      Setup name (for MasterFrame)
    settings : dict-like
      Settings to guide the fluxing
      Key ones are ['mosaic']['longitude', 'latitude', 'elevation']


    Attributes
    ----------
    sensfunc : dict
      Sensitivity function

    steps : list
      List of steps performed
    frametype : str
      Set to 'sensfunc'

    std : SpecObj
      The chosen one for generating the sensitivity function
    std_header : dict-like
      Used for the RA, DEC, AIRMASS, EXPTIME of the standard star spectrum
    std_idx : int
      Index that std is within the std_specbojs list
    sci_specobjs : list
      List of SpecObj objects to be fluxed (or that were fluxed)
    sci_header : dict-like
      Used for the airmass, exptime of the science spectra
    """
    def __init__(self, std_spec1d_file=None, sci_spec1d_file=None, spectrograph=None,
                 sens_file=None, setup=None, settings=None, std_specobjs=None):

        # Parameters
        self.std_spec1d_file = std_spec1d_file
        self.sci_spec1d_file = sci_spec1d_file
        self.sens_file = sens_file
        self.setup = setup
        if std_specobjs is not None:
            # Need to unwrap these (sometimes)..
            self.std_specobjs = arutils.unravel_specobjs(std_specobjs)

        # Main outputs
        self.sensfunc = None  # dict

        # Attributes
        self.frametype = frametype
        self.steps = []

        # Settings (for now)
        if spectrograph is not None:
            if settings is None:
                self.settings = kludge_settings(spectrograph)
        else:
            self.settings = settings

        # Key Internals
        self.std = None         # Standard star spectrum (SpecObj object)
        self.std_header = None
        self.std_idx = None     # Nested indices for the std_specobjs list that corresponds to the star!
        self.sci_specobjs = []
        self.sci_header = None

        # Load
        if self.std_spec1d_file is not None:
            self.std_specobjs, self.std_header = arload.load_specobj(self.std_spec1d_file)
            msgs.info("Loaded {:d} spectra from the spec1d standard star file: {:s}".format(len(self.std_specobjs),
                                                                                            self.std_spec1d_file))
        if self.sci_spec1d_file is not None:
            self.sci_specobjs, self.sci_header = arload.load_specobj(self.sci_spec1d_file)
            msgs.info("Loaded spectra from the spec1d science file: {:s}".format(self.sci_spec1d_file))
        if self.sens_file is not None:
            self.sensfunc, _, _ = armasters._load(self.sens_file, frametype=self.frametype)

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, self.setup, self.settings)

    def find_standard(self):
        """
        Identify the standard star from the list of all spectra in the specobjs
          Wrapper to arflux.find_standard which simply takes the brightest

        Returns
        -------
        self.std : SpecObj
          Corresponds to the chosen spectrum
        """
        # Find brightest object in the exposures
        # Searches over all slits (over all detectors), and all objects
        self.std_idx = arflux.find_standard(self.std_specobjs)
        # Set internal
        self.std = self.std_specobjs[self.std_idx]
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.std

    def generate_sensfunc(self):
        """
        Generate the senstivity function

        Wrapper to arflux.generate_sensfunc
          Requires self.std has been set

        Returns
        -------
        self.sensfunc : dict

        """
        if self.std is None:
            msgs.warn("You need to identify the Star first (with find_standard)") #load in the Standard spectra first")
            return None
        if self.std_header is None:
            msgs.warn("You need to set the std_header with a dict-like object holding RA, DEC, AIRMASS, EXPTIME")
            return None
        # Run
        self.sensfunc = arflux.generate_sensfunc(self.std, self.std_header['RA'],
                                            self.std_header['DEC'], self.std_header['AIRMASS'],
                                            self.std_header['EXPTIME'], self.settings)
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.sensfunc

    def _flux_specobjs(self, specobjs, airmass, exptime):
        """
        Flux an input list of SpecObj objects
          Can be packed in detector, slit, etc. or as a simple list

        Wrapper to arflux.apply_sensfunc()

        Parameters
        ----------
        specobjs : list
        airmass : float
        exptime : float

        Returns
        -------

        """
        # Note the unravel here
        for sci_obj in arutils.unravel_specobjs(specobjs):
            if sci_obj is not None:
                # Do it
                arflux.new_apply_sensfunc(sci_obj, self.sensfunc, airmass, exptime, self.settings)


    def flux_science(self):
        """
        Flux the internal list of sci_specobjs

        Wrapper to arflux.apply_sensfunc()

        Returns
        -------

        """
        # Settings kludge if the Class is being used external to PYPIT
        if 'mosaic' not in self.settings.keys():
            self.settings['mosaic'] = {}
            self.settings['mosaic']['longitude'] = self.sci_header['LON-OBS']
            self.settings['mosaic']['latitude'] = self.sci_header['LAT-OBS']
            self.settings['mosaic']['elevation'] = self.sci_header['ALT-OBS']
        # Flux
        for sci_obj in self.sci_specobjs:
            arflux.new_apply_sensfunc(sci_obj, self.sensfunc,
                                  self.sci_header['AIRMASS'],
                                  self.sci_header['EXPTIME'],
                                  self.settings)
        # Step
        self.steps.append(inspect.stack()[0][3])

    def _set_std_obj(self, obj_id):
        """
        Method which allows the user to identify the standard star
          with an input

        Parameters
        ----------
        obj_id : str or int
          If str, object id of the spectrum
          If int, index in the internal std_specobjs list

        Returns
        -------
        self.std : SpecObj

        """
        if self.std_specobjs is None:
            msgs.warn("You need to load in the Standard spectra first!")
            return None
        #
        if isinstance(obj_id, basestring):
            names = [spobj.idx for spobj in self.std_specobjs]
            self.std_idx = names.index(obj_id)
        elif isinstance(obj_id, int): # Extension number
            self.std_idx = obj_id-1
        self.std = self.std_specobjs[self.std_idx]
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.std

    def master(self, row_fitstbl, clobber=False):
        """
        Load or generate+save the MasterFrame sensitivity function
        The solution is applied to the list of science spectra loaded (if any)

        Parameters
        ----------
        row_fitstbl : Row
          Row of the fitstbl corresponding to the Standard star
          Used to parse RA, DEC, AIRMASS, EXPTIME
        clobber : bool, optional

        Returns
        -------
        self.sensfunc

        """
        self.sensfunc, _, _ = self.load_master_frame()
        # Sensitivity Function
        if (self.sensfunc is None) or clobber:
            if self.std_specobjs is None:
                msgs.warn("You need to load in the Standard spectra first!")
                return None

            # Find the star automatically?
            _ = self.find_standard()

            # Kludge a header
            if self.std_header is None:
                self.std_header={}
                for key in ['ra','dec','airmass','exptime']:
                    self.std_header[key.upper()] = row_fitstbl[key]

            # Sensitivity
            _ = self.generate_sensfunc()

            # Save to master
            self.save_master()

        # Apply to science
        if len(self.sci_specobjs) > 0:
            self.flux_science()
        # Return
        return self.sensfunc

    def save_master(self, outfile=None):
        """
        Over-load the save_master() method in MasterFrame to write a YAML file

        Parameters
        ----------
        outfile : str, optional
          Use this input instead of the 'proper' (or unattainable) MasterFrame name

        Returns
        -------

        """
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Allow one to over-ride output name
        if outfile is None:
            outfile = self.ms_name
        # Add steps
        self.sensfunc['steps'] = self.steps
        # yamlify
        ysens = arutils.yamlify(self.sensfunc)
        with open(outfile, 'w') as yamlf:
            yamlf.write(yaml.dump(ysens))
        #
        msgs.info("Wrote sensfunc to MasterFrame: {:s}".format(outfile))

    def show_sensfunc(self):
        """
        Plot the sensitivity function
        """
        if self.sensfunc is None:
            msgs.warn("You need to generate the sensfunc first!")
            return None
        # Generate from model
        wave = np.linspace(self.sensfunc['wave_min'], self.sensfunc['wave_max'], 1000)
        mag_func = arutils.func_val(self.sensfunc['c'], wave, self.sensfunc['func'])
        sens = 10.0**(0.4*mag_func)
        # Plot
        debugger.plot1d(wave, sens, xlbl='Wavelength', ylbl='Sensitivity Function')


    def write_science(self, outfile):
        """
        Write the flux-calibrated science spectra

        Parameters
        ----------
        outfile : str

        Returns
        -------

        """
        if len(self.sci_specobjs) == 0:
            msgs.warn("No science spectra to write to disk!")
        #
        if 'VEL-TYPE' in self.sci_header.keys():
            helio_dict = dict(refframe=self.sci_header['VEL-TYPE'], vel_correction=self.sci_header['VEL'])
        else:
            helio_dict = None
        if 'LON-OBS' in self.sci_header.keys():
            obs_dict = dict(longitude=self.sci_header['LON-OBS'],
                            latitude=self.sci_header['LAT-OBS'],
                            elevation=self.sci_header['ALT-OBS'],
                            )
        arsave.new_save_1d_spectra_fits(self.sci_specobjs, self.sci_header, outfile,
                             helio_dict=helio_dict, clobber=True, obs_dict=obs_dict)
        # Step
        self.steps.append(inspect.stack()[0][3])


    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: '.format(self.__class__.__name__)
        if len(self.steps) > 0:
            txt+= ' steps: ['
            for step in self.steps:
                txt += '{:s}, '.format(step)
            txt = txt[:-2]+']'  # Trim the trailing comma
        txt += '>'
        return txt



