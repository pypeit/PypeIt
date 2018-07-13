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

from astropy import units

from pypit import msgs
from pypit.core import arflux
from pypit import arload
from pypit.core import armasters
from pypit.core import arsave
from pypit import arutils
from pypit import masterframe

from pypit.spectrographs.spectrograph import Spectrograph
from pypit.spectrographs.util import load_spectrograph
from pypit.par.pypitpar import TelescopePar

from pypit import ardebug as debugger

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

    multi_det : tuple, optional
      List of detector numbers to splice together for multi-detector instruments (e.g. DEIMOS)
      They are assumed to be in order of increasing wavelength
      And that there is *no* overlap in wavelength across detectors (might be ok if there is)
    std : SpecObj
      The chosen one for generating the sensitivity function
    std_header : dict-like
      Used for the RA, DEC, AIRMASS, EXPTIME of the standard star spectrum
    std_idx : int or list
      Index that std is within the std_specbojs list
    sci_specobjs : list
      List of SpecObj objects to be fluxed (or that were fluxed)
    sci_header : dict-like
      Used for the airmass, exptime of the science spectra
    """

    # Frametype is a class attribute
    frametype = 'sensfunc'

    def __init__(self, std_spec1d_file=None, sci_spec1d_file=None, sens_file=None,
                 std_specobjs=None, std_header=None, spectrograph=None, multi_det=None,
                 setup=None, root_path=None, mode=None):

        # Load standard files
        std_spectro = None
        self.std_spec1d_file = std_spec1d_file
        # Need to unwrap these (sometimes)..
        self.std_specobjs = std_specobjs if std_specobjs is None \
                                            else arutils.unravel_specobjs(std_specobjs)
        self.std_header = std_header
        if self.std_spec1d_file is not None:
            self.std_specobjs, self.std_header = arload.load_specobj(self.std_spec1d_file)
            msgs.info('Loaded {0} spectra from the spec1d standard star file: {1}'.format(
                                len(self.std_specobjs), self.std_spec1d_file))
            std_spectro = self.std_header['INSTRUME']

        # Load the science files
        sci_spectro = None
        self.sci_spec1d_file = sci_spec1d_file
        self.sci_specobjs = []
        self.sci_header = None
        if self.sci_spec1d_file is not None:
            self.sci_specobjs, self.sci_header = arload.load_specobj(self.sci_spec1d_file)
            msgs.info('Loaded {0} spectra from the spec1d science file: {1}'.format(
                                len(self.sci_specobjs), self.sci_spec1d_file))
            sci_spectro = self.sci_header['INSTRUME']

        # Compare instruments if they exist
        if std_spectro is not None and sci_spectro is not None and std_spectro != sci_spectro:
            msgs.error('Standard spectra are not the same instrument as science!!')

        # Instantiate the spectograph
        _spectrograph = spectrograph
        spectrograph_name = None
        if _spectrograph is None:
            _spectrograph = std_spectro
            spectrograph_name = _spectrograph
            if _spectrograph is not None:
                msgs.info("Spectrograph set to {0} from standard file".format(_spectrograph))
        if _spectrograph is None:
            _spectrograph = sci_spectro
            spectrograph_name = _spectrograph
            if _spectrograph is not None:
                msgs.info("Spectrograph set to {0} from science file".format(_spectrograph))
        if isinstance(_spectrograph, basestring):
            spectrograph_name = _spectrograph
            msgs.info("Spectrograph set to {0}, from argument string".format(_spectrograph))
        elif isinstance(_spectrograph, Spectrograph):
            spectrograph_name = _spectrograph.spectrograph
            msgs.info("Spectrograph set to {0}, from argument object".format(_spectrograph))
    
        # MasterFrame
        directory_path = None
        if root_path is not None:
            directory_path = root_path
            if spectrograph_name is not None:
                directory_path += '_'+spectrograph_name
        masterframe.MasterFrame.__init__(self, self.frametype, setup,
                                         directory_path=directory_path, mode=mode)
        # Get the extinction data
        self.extinction_data = None
        if _spectrograph is not None:
            telescope = load_spectrograph(spectrograph=_spectrograph).telescope \
                                    if isinstance(_spectrograph, basestring) \
                                    else _spectrograph.telescope
            self.extinction_data \
                    = arflux.load_extinction_data(telescope['longitude'], telescope['latitude'])
        elif self.sci_header is not None and 'LON-OBS' in self.sci_header.keys():
            self.extinction_data \
                    = arflux.load_extinction_data(self.sci_header['LON-OBS'],
                                                  self.sci_header['LAT-OBS'])
       
        # Once the spectrograph is instantiated, can also set the
        # extinction data
        # Parameters
        self.sens_file = sens_file
        self.multi_det = multi_det

        # Main outputs
        self.sensfunc = None if self.sens_file is None \
                            else armasters._load(self.sens_file, frametype=self.frametype)[0] 

        # Attributes
        self.steps = []

        # Key Internals
        self.std = None         # Standard star spectrum (SpecObj object)
        self.std_idx = None     # Nested indices for the std_specobjs list that corresponds
                                # to the star!


    def find_standard(self):
        """
        Identify the standard star from the list of all spectra in the specobjs

          Wrapper to arflux.find_standard which simply takes the brightest

        Returns
        -------
        self.std : SpecObj
          Corresponds to the chosen spectrum
        """
        if self.multi_det is not None:
            sv_stds = []
            # Find the standard in each detector
            for det in self.multi_det:
                stds = [sobj for sobj in self.std_specobjs if sobj.det == det]
                idx = arflux.find_standard(stds)
                sv_stds.append(stds[idx])
                msgs.info("Using standard {} for det={}".format(stds[idx], det))
            # Now splice
            msgs.info("Splicing the standards -- The name will be for the first detector")
            std_splice = sv_stds[0].copy()
            # Append
            for ostd in sv_stds[1:]:
                std_splice.boxcar['wave'] = np.append(std_splice.boxcar['wave'].value,
                                                      ostd.boxcar['wave'].value) * units.AA
                for key in ['counts', 'var']:
                    std_splice.boxcar[key] = np.append(std_splice.boxcar[key], ostd.boxcar[key])
            self.std = std_splice
        else:
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
        # Check internals
        if self.std is None:
            msgs.warn('First identify the star first (with find_standard).')
            return None
        if self.std_header is None:
            msgs.warn('First set std_header with a dict-like object holding RA, DEC, '
                      'AIRMASS, EXPTIME.')
            return None

        # Get extinction correction
        extinction_corr = arflux.extinction_correction(self.std.boxcar['wave'],
                                                       self.std_header['AIRMASS'],
                                                       self.extinction_data)
        self.sensfunc = arflux.generate_sensfunc(self.std, self.std_header['RA'],
                                                 self.std_header['DEC'],
                                                 self.std_header['EXPTIME'], extinction_corr)

        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.sensfunc

    def flux_specobjs(self, specobjs, airmass, exptime):
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
                arflux.apply_sensfunc(sci_obj, self.sensfunc, airmass, exptime,
                                      self.extinction_data)

    def flux_science(self):
        """
        Flux the internal list of sci_specobjs

        Wrapper to arflux.apply_sensfunc()

        Returns
        -------

        """
        for sci_obj in self.sci_specobjs:
            arflux.apply_sensfunc(sci_obj, self.sensfunc, self.sci_header['AIRMASS'],
                                  self.sci_header['EXPTIME'], self.extinction_data)
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

    def master(self, row_fitstbl, clobber=False, save=True):
        """
        Load or generate+save the MasterFrame sensitivity function
        The solution is applied to the list of science spectra loaded (if any)

        Parameters
        ----------
        row_fitstbl : Row
          Row of the fitstbl corresponding to the Standard star
          Used to parse RA, DEC, AIRMASS, EXPTIME
        clobber : bool, optional
        save : bool, optional
          Save to master frame?

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
            if save:
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
            helio_dict = dict(refframe=self.sci_header['VEL-TYPE'],
                              vel_correction=self.sci_header['VEL'])
        else:
            helio_dict = None
        telescope=None
        if 'LON-OBS' in self.sci_header.keys():
            telescope = TelescopePar(longitude=self.sci_header['LON-OBS'],
                                     latitude=self.sci_header['LAT-OBS'],
                                     elevation=self.sci_header['ALT-OBS'])
        arsave.save_1d_spectra_fits(self.sci_specobjs, self.sci_header, outfile,
                                    helio_dict=helio_dict, telescope=telescope, clobber=True)
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

