# Module for guiding Slit/Order tracing
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import linetools
import os
import json

#from importlib import reload

from astropy import units
from astropy.io import fits

from pypeit import msgs
from pypeit.core import flux
from pypeit.core import load
from pypeit.core import save
from pypeit import utils
from pypeit import masterframe
from pypeit import specobjs

from pypeit.spectrographs.util import load_spectrograph
from pypeit.par.pypeitpar import TelescopePar

from pypeit import debugger

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
      Used only to set settings for calls to the Class outside of PyepIt
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
    sci_specobjs : SpecObjs
      List of SpecObj objects to be fluxed (or that were fluxed)
    sci_header : dict-like
      Used for the airmass, exptime of the science spectra
    spectrograph : Spectrograph
      Used for extinction correction
    """

    # Frametype is a class attribute
    frametype = 'sensfunc'

    def __init__(self, std_spec1d_file=None, sci_spec1d_file=None, sens_file=None,
                 std_specobjs=None, std_header=None, spectrograph=None, multi_det=None,
                 telluric=False, setup=None, master_dir=None, mode=None):

        # Load standard files
        std_spectro = None
        self.std_spec1d_file = std_spec1d_file
        # Need to unwrap these (sometimes)..
        self.std_specobjs = std_specobjs
        self.std_header = std_header
        if self.std_spec1d_file is not None:
            self.std_specobjs, self.std_header = load.load_specobj(self.std_spec1d_file)
            msgs.info('Loaded {0} spectra from the spec1d standard star file: {1}'.format(
                                len(self.std_specobjs), self.std_spec1d_file))
            std_spectro = self.std_header['INSTRUME']

        try:
            self.std_ra = self.std_header['RA']
        except:
            self.std_ra = None
        try:
            self.std_dec = self.std_header['DEC']
        except:
            self.std_dec = None
        try:
            self.std_file = self.std_header['FILENAME']
        except:
            self.std_file = None

        # Load the science files
        sci_spectro = None
        self.sci_spec1d_file = sci_spec1d_file
        self.sci_specobjs = []
        self.sci_header = None
        if self.sci_spec1d_file is not None:
            self.sci_specobjs, self.sci_header = load.load_specobj(self.sci_spec1d_file)
            msgs.info('Loaded {0} spectra from the spec1d science file: {1}'.format(
                                len(self.sci_specobjs), self.sci_spec1d_file))
            sci_spectro = self.sci_header['INSTRUME']

        # Compare instruments if they exist
        if std_spectro is not None and sci_spectro is not None and std_spectro != sci_spectro:
            msgs.error('Standard spectra are not the same instrument as science!!')

        # Instantiate the spectrograph
        _spectrograph = spectrograph
        if _spectrograph is None:
            _spectrograph = std_spectro
            if _spectrograph is not None:
                msgs.info("Spectrograph set to {0} from standard file".format(_spectrograph))
        if _spectrograph is None:
            _spectrograph = sci_spectro
            if _spectrograph is not None:
                msgs.info("Spectrograph set to {0} from science file".format(_spectrograph))
        self.spectrograph = load_spectrograph(_spectrograph)

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, setup,
                                         master_dir=master_dir, mode=mode)
        # Get the extinction data
        self.extinction_data = None
        if self.spectrograph is not None:
            self.extinction_data \
                    = flux.load_extinction_data(self.spectrograph.telescope['longitude'],
                                                self.spectrograph.telescope['latitude'])
        elif self.sci_header is not None and 'LON-OBS' in self.sci_header.keys():
            self.extinction_data \
                    = flux.load_extinction_data(self.sci_header['LON-OBS'],
                                                self.sci_header['LAT-OBS'])
       
        # Once the spectrograph is instantiated, can also set the
        # extinction data
        # Parameters
        self.sens_file = sens_file
        self.multi_det = multi_det

        # Set telluric option
        self.telluric = telluric

        # Main outputs
        self.sens_dict = None if self.sens_file is None \
                            else self.load_master(self.sens_file)

        # Attributes
        self.steps = []

        # Key Internals
        self.std = None         # Standard star spectrum (SpecObj object)
        self.std_idx = None     # Nested indices for the std_specobjs list that corresponds
                                # to the star!


    def find_standard(self):
        """
        Identify the standard star from the list of all spectra in the specobjs

          Wrapper to flux.find_standard which simply takes the brightest

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
                if len(stds) == 0:
                    debugger.set_trace()
                idx = flux.find_standard(stds)
                sv_stds.append(stds[idx])
                msgs.info("Using standard {} for det={}".format(stds[idx], det))
            # Now splice
            msgs.info("Splicing the standards -- The name will be for the first detector")
            std_splice = sv_stds[0].copy()
            # Append
            for ostd in sv_stds[1:]:
                std_splice.boxcar['WAVE'] = np.append(std_splice.boxcar['WAVE'].value,
                                                      ostd.boxcar['WAVE'].value) * units.AA
                for key in ['COUNTS', 'COUNTS_IVAR']:
                    std_splice.boxcar[key] = np.append(std_splice.boxcar[key], ostd.boxcar[key])
            self.std = std_splice
        else:
            # Find brightest object in the exposures
            # Searches over all slits (over all detectors), and all objects
            self.std_idx = flux.find_standard(self.std_specobjs)
            # Set internal
            self.std = self.std_specobjs[self.std_idx]
            # Step
            self.steps.append(inspect.stack()[0][3])
            # Return
            return self.std

    def generate_sensfunc(self):
        """
        Generate the senstivity function

        Wrapper to flux.generate_sensfunc
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
        #extinction_corr = flux.extinction_correction(self.std.boxcar['WAVE'],
                                                       #self.std_header['AIRMASS'],
                                                       #self.extinction_data)
        #self.sensfunc = flux.generate_sensfunc(self.std, self.std_header['RA'],
                                                 #self.std_header['DEC'],
                                                 #self.std_header['EXPTIME'], extinction_corr)
        # extinction_corr = flux.extinction_correction(self.std.boxcar['WAVE'],
        #                                              self.std_header['AIRMASS'],
        #                                              self.extinction_data)

        self.sens_dict = flux.generate_sensfunc(self.std.boxcar['WAVE'],
                                               self.std.boxcar['COUNTS'],
                                               self.std.boxcar['COUNTS_IVAR'],
                                               self.std_header['AIRMASS'],
                                               self.std_header['EXPTIME'],
                                               self.spectrograph,
                                               telluric=self.telluric,
                                               ra=self.std_ra,
                                               dec=self.std_dec,
                                               std_file = self.std_file)
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.sens_dict

    def flux_specobjs(self, specobjs, airmass, exptime):
        """
        Flux an input list of SpecObj objects
          Can be packed in detector, slit, etc. or as a simple list

        Wrapper to flux.apply_sensfunc()

        Parameters
        ----------
        specobjs : list
        airmass : float
        exptime : float

        Returns
        -------

        """
        # Note the unravel here
        for sci_obj in utils.unravel_specobjs(specobjs):
            if sci_obj is not None:
                # Do it
                flux.apply_sensfunc(sci_obj, self.sens_dict, airmass, exptime,
                                    self.spectrograph)

    def flux_science(self):
        """
        Flux the internal list of sci_specobjs

        Wrapper to flux.apply_sensfunc()

        Returns
        -------

        """
        for sci_obj in self.sci_specobjs:
            flux.apply_sensfunc(sci_obj, self.sens_dict, self.sci_header['AIRMASS'],
                                  self.sci_header['EXPTIME'], self.spectrograph)
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
        if isinstance(obj_id, str):
            names = [spobj.idx for spobj in self.std_specobjs]
            self.std_idx = names.index(obj_id)
        elif isinstance(obj_id, int): # Extension number
            self.std_idx = obj_id-1
        self.std = self.std_specobjs[self.std_idx]
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.std

    def get_sens_dict(self, row_fitstbl, clobber=False, save=True):
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
        self.sens_dict = self.master()
        # Sensitivity Function
        if (self.sens_dict is None) or clobber:
            if self.std_specobjs is None:
                msgs.warn("You need to load in the Standard spectra first!")
                return None

            # Find the star automatically?
            _ = self.find_standard()

            # Kludge a header
            if self.std_header is None:
                self.std_header={}
                for key in ['ra','dec','airmass','exptime','filename']:
                    self.std_header[key.upper()] = row_fitstbl[key]

            # Sensitivity
            _ = self.generate_sensfunc()

            # Save to master
            if save:
                self.save_master(self.sens_dict)

        # Apply to science
        if len(self.sci_specobjs) > 0:
            self.flux_science()
        # Return
        return self.sens_dict

    def load_master(self, filename, force=False):


        # Does the master file exist?
        if not os.path.isfile(filename):
            msgs.warn("No Master frame found of type {:s}: {:s}".format(self.frametype, filename))
            if force:
                msgs.error("Crashing out because reduce-masters-force=True:" + msgs.newline() + filename)
            return None
        else:
            msgs.info("Loading a pre-existing master calibration frame of type: {:}".format(self.frametype) + " from filename: {:}".format(filename))

            hdu = fits.open(filename)
            head = hdu[0].header
            tbl = hdu['SENSFUNC'].data
            sens_dict = {}
            sens_dict['wave'] = tbl['WAVE']
            sens_dict['sensfunc'] = tbl['SENSFUNC']
            for key in ['wave_min','wave_max','exptime','airmass','std_file','std_ra','std_dec','std_name','cal_file']:
                try:
                    sens_dict[key] = head[key.upper()]
                except:
                    pass
            return sens_dict

    def save_master(self, sens_dict, outfile=None):
        """
        Over-load the save_master() method in MasterFrame to write a JSON file

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
        self.sens_dict['steps'] = self.steps
        # Do it
        prihdu = fits.PrimaryHDU()
        hdus = [prihdu]
        # Add critical keys from sens_dict to header
        for key in ['wave_min', 'wave_max', 'exptime', 'airmass', 'std_file', 'std_ra',
                    'std_dec', 'std_name', 'cal_file']:
            try:
                prihdu.header[key.upper()] = sens_dict[key].value
            except AttributeError:
                prihdu.header[key.upper()] = sens_dict[key]
            except KeyError:
                pass  # Will not require all of these

        cols = []
        cols += [fits.Column(array=sens_dict['wave'], name=str('WAVE'), format=sens_dict['wave'].dtype)]
        cols += [
            fits.Column(array=sens_dict['sensfunc'], name=str('SENSFUNC'), format=sens_dict['sensfunc'].dtype)]
        # Finish
        coldefs = fits.ColDefs(cols)
        tbhdu = fits.BinTableHDU.from_columns(coldefs)
        tbhdu.name = 'SENSFUNC'
        hdus += [tbhdu]
        # Finish
        hdulist = fits.HDUList(hdus)
        hdulist.writeto(outfile, overwrite=True)

        # Finish
        msgs.info("Wrote sensfunc to MasterFrame: {:s}".format(outfile))


    def save_master_old(self, outfile=None):
        """
        Over-load the save_master() method in MasterFrame to write a JSON file

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
        self.sens_dict['steps'] = self.steps
        # # yamlify
        # ysens = utils.yamlify(self.sensfunc)
        # with open(outfile, 'w') as yamlf:
        #     yamlf.write(yaml.dump(ysens))
        #  #
        # msgs.info("Wrote sensfunc to MasterFrame: {:s}".format(outfile))

        # jsonify
        self.sensfunc['mag_set'] = self.sensfunc['mag_set'].to_dict()
        jsensfunc = self.sensfunc.copy()
        jsensfunc['bspline'] = jsensfunc['bspline'].to_dict()
        jsensfunc = linetools.utils.jsonify(jsensfunc)
        with open(outfile, 'w') as jsonf:
            linetools.utils.savejson(outfile, jsensfunc, overwrite=True, indent=4, easy_to_read=True)

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
        mag_func = utils.func_val(self.sensfunc['c'], wave, self.sensfunc['func'])
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
        # KLUDGE ME
        if isinstance(self.sci_specobjs, list):
            specObjs = specobjs.SpecObjs(self.sci_specobjs)
        elif isinstance(self.sci_specobjs, specobjs.SpecObjs):
            specObjs = self.sci_specobjs
        else:
            msgs.error("BAD INPUT")
        save.save_1d_spectra_fits(specObjs, self.sci_header, outfile,
                                  helio_dict=helio_dict,
                                  telescope=telescope, overwrite=True)
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

