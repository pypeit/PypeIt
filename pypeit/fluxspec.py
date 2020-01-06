# Module for guiding Slit/Order tracing
import inspect
import numpy as np
import linetools
import os
import json
import matplotlib.pyplot as plt

#from importlib import reload

from astropy import units
from astropy.io import fits

from pypeit import msgs
from pypeit.core import flux_calib
from pypeit.core import load
from pypeit.core import save
from pypeit import specobjs
from pypeit import debugger

from IPython import embed

class FluxSpec(object):
    """Class to guide fluxing

    Parameters
    ----------
    spectrograph : Spectrograph or str
      Name of the spectrograph, e.g. shane_kast_blue
      Used only to set settings for calls to the Class outside of PypeIt
      This includes extinction data..
    par: FluxCalib parset
    sens_file : str, optional
      Filename of a sensitivity function file to be input

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

    # JFH TODO In my opinion the things that the class operates on should be arguments and the parameters
    # should be in the parset. I really don't like to see parsets used in this way where everything is passed in
    # under the hood. It makes it challenging to understand what the code is doing.
    def __init__(self, spectrograph, par, sens_file=None, debug=False):

        # Init
        self.spectrograph = spectrograph
        self.par = par

        # Get the extinction data
        self.extinction_data = flux_calib.load_extinction_data(
                self.spectrograph.telescope['longitude'], self.spectrograph.telescope['latitude'])

        # Parameters
        if sens_file is None:
            self.sens_file = par['sensfunc']
        else:
            self.sens_file = sens_file
        self.multi_det = par['multi_det']

        # Set telluric option
        self.telluric = par['telluric']

        # Main outputs
        self.sens_dict = None if self.sens_file is None \
                            else self.load_sens_dict(self.sens_file)
        # Attributes
        self.steps = []

        # Key Internals
        self.std = None         # Standard star spectrum (SpecObj object)
        self.std_idx = None     # Nested indices for the std_specobjs list that corresponds
                                # to the star!
        # standard/telluric star information
        self.star_type = par['star_type']
        self.star_mag = par['star_mag']
        # telluric mask keywords
        self.BALM_MASK_WID = par['balm_mask_wid']
        # sensfunc fitting parameters
        self.poly_norder = par['poly_norder']
        self.polycorrect = par['polycorrect']
        self.debug = debug

    def load_objs(self, spec1d_file, std=True):
        """
        Load specobjs and heade from an input spec1d_file

        Args:
            spec1d_file: str
            std: bool, optional
              If True, load up standard star file and header
              If False, load up science file and header

        Returns:
            Loads up self.std_specobjs or self.sci_specobjs

        """
        sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file)
        header = fits.open(spec1d_file)[0].header
        #
        if std:
            self.std_specobjs, self.std_header = sobjs, header
            msgs.info('Loaded {0} spectra from the spec1d standard star file: {1}'.format(
                len(self.std_specobjs), spec1d_file))
            self.std_ra = self.std_header['RA']
            self.std_dec = self.std_header['DEC']
            self.std_file = self.std_header['FILENAME']
        else:
            self.sci_specobjs, self.sci_header = sobjs, header
            msgs.info('Loaded {0} spectra from the spec1d science file: {1}'.format(
                len(self.sci_specobjs), spec1d_file))
        # Check instrument
        spectro = header['PYP_SPEC']
        assert spectro == self.spectrograph.spectrograph
        return sobjs, header

    def find_standard(self):
        """
        Identify the standard star from the list of all spectra in the specobjs

          Wrapper to flux.find_standard which simply takes the brightest

        Returns
        -------
        self.std : SpecObj
          Corresponds to the chosen spectrum
        """
        if self.par['std_obj_id'] is not None:
            _ = self._set_std_obj()
            return
        if self.multi_det is not None:
            sv_stds = []
            # Find the standard in each detector
            for det in self.multi_det:
                stds = [sobj for sobj in self.std_specobjs if sobj.det == det]
                if len(stds) == 0:
                    debugger.set_trace()
                idx = flux_calib.find_standard(stds)
                sv_stds.append(stds[idx])
                msgs.info("Using standard {} for det={}".format(stds[idx], det))

            # Now splice
            msgs.info("Splicing the standards -- The name will be for the first detector")
            std_splice = sv_stds[0].copy()
            # Append
            for ostd in sv_stds[1:]:
                try:
                    std_splice.optimal['WAVE_GRID'] = np.append(std_splice.optimal['WAVE_GRID'].value,
                                                           ostd.optimal['WAVE_GRID'].value) * units.AA
                except KeyError:
                    std_splice.optimal['WAVE'] = np.append(std_splice.optimal['WAVE'].value,
                                                           ostd.optimal['WAVE'].value) * units.AA
                for key in ['COUNTS', 'COUNTS_IVAR']:
                    std_splice.optimal[key] = np.append(std_splice.optimal[key], ostd.optimal[key])
            self.std = std_splice
        elif self.spectrograph.pypeline == 'Echelle':
            # Find brightest object in each order
            std_brightest = self.std_specobjs[flux_calib.find_standard(self.std_specobjs)]
            self.std_idx = np.isclose(std_brightest.ECH_FRACPOS, self.std_specobjs.ECH_FRACPOS)
            # Set internal
            self.std = self.std_specobjs[self.std_idx]
            # Step
            self.steps.append(inspect.stack()[0][3])
            # Return
            return self.std
        else:
            # Find brightest object in the exposures
            # Searches over all slits (over all detectors), and all objects
            self.std_idx = flux_calib.find_standard(self.std_specobjs)
            # Set internal
            self.std = self.std_specobjs[self.std_idx]
            # Step
            self.steps.append(inspect.stack()[0][3])
            # Return
            return self.std

    def generate_sensfunc(self):
        """
         Dummy method for generating sensfunc. Overloaded by class specific sensfunc generator.
         Returns:
         """
        return None

    def flux_science(self, sci_file):
        """
        Dummy method for fluxing. Overloaded by class specific fluxing
        Returns
        """
        return

    def _set_std_obj(self, obj_id=None):
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
        # From parameters?
        if obj_id is None:
            obj_id = self.par['std_obj_id']
        #
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

    def save_sens_dict(self, sens_dict, outfile):
        """
        Save the sens_dict. Wrapper for save.save_sens_dict
        Args:
            sens_dict:
            outfile:

        Returns:

        """
        save.save_sens_dict(sens_dict, outfile)

    def load_sens_dict(self, filename):
        self.sens_dict = load.load_sens_dict(filename)
        return self.sens_dict


    # TODO Need to improve this QA, it is really not informative
    # ToDo: either make it a dummy function or make it works for both multislit and echelle.
    def show_sensfunc(self):
        """
        Plot the sensitivity function
        """
        if self.sens_dict is None:
            msgs.warn("You need to generate the sens_dict first!")
            return None
        # Generate from model
        #wave = np.linspace(self.sens_dict['wave_min'], self.sens_dict['wave_max'], 1000)
        #mag_func = utils.func_val(self.sens_dict['c'], wave, self.sens_dict['func'])
        #sens = 10.0**(0.4*mag_func)
        # Plot
        if self.spectrograph.pypeline == 'Echelle':
            for iord in range(len(self.std)):
                sord = str(iord)
                debugger.plot1d(self.sens_dict[sord]['wave'],
                                self.sens_dict[sord]['sensfunc'],
                                xlbl='Wavelength', ylbl='Sensitivity Function')
        else:
            debugger.plot1d(self.sens_dict['0']['wave'], self.sens_dict['0']['sensfunc'], xlbl='Wavelength', ylbl='Sensitivity Function')

    def write_science(self, outfile):
        """
        Write the flux-calibrated science spectra

        Args:
            outfile (str):
        """
        if len(self.sci_specobjs) == 0:
            msgs.warn("No science spectra to write to disk!")
            return

        self.sci_specobjs.write_to_fits(outfile, header=self.sci_header,
                                        spectrograph=self.spectrograph, overwrite=True)
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

class MultiSlit(FluxSpec):
    """
    Child of FluxSpec for Multislit and Longslit reductions

    """
    def __init__(self, spectrograph, par, **kwargs):
        super(MultiSlit, self).__init__(spectrograph, par, **kwargs)

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

        self.sens_dict = {}
        this_wave = self.std.OPT_WAVE
        sens_dict_long = flux_calib.generate_sensfunc(this_wave,
                                                      self.std.OPT_COUNTS,
                                                      self.std.OPT_COUNTS_IVAR,
                                                      self.std_header['AIRMASS'],
                                                      self.std_header['EXPTIME'],
                                                      self.spectrograph.telescope['longitude'],
                                                      self.spectrograph.telescope['latitude'],
                                                      BALM_MASK_WID=self.par['balm_mask_wid'],
                                                      star_type=self.star_type,
                                                      star_mag=self.star_mag,
                                                      telluric=self.telluric,
                                                      ra=self.std_ra,
                                                      dec=self.std_dec,
                                                      std_file = self.std_file,
                                                      poly_norder=self.poly_norder,
                                                      polycorrect=self.polycorrect,
                                                      debug=self.debug)
        self.sens_dict['0'] = sens_dict_long
        self.sens_dict['nslits'] = 1

        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.sens_dict

    def flux_science(self, sci_file):
        """
        Flux the internal list of sci_specobjs

        Wrapper to flux.apply_sensfunc()

        Returns
        -------

        """
        # Load into self.sci_specobjs
        self.load_objs(sci_file, std=False)

        # Run
        for sci_obj in self.sci_specobjs:
            sci_obj.apply_flux_calib(self.sens_dict['0'],
                                     self.sci_header['EXPTIME'],
                                     telluric_correct=self.par['telluric_correct'],
                                     extinct_correct=self.par['extinct_correct'],
                                     longitude=self.spectrograph.telescope['longitude'],
                                     latitude=self.spectrograph.telescope['latitude'],
                                     airmass=self.sci_header['AIRMASS'])
        self.steps.append(inspect.stack()[0][3])


class Echelle(FluxSpec):
    """
    Child of FluxSpec for Echelle reductions

    """
    def __init__(self, spectrograph, par, **kwargs):
        super(Echelle, self).__init__(spectrograph, par, **kwargs)

        # Echelle key
        # TODO add these to the parameters to the parset or try to get rid of these parameters in flux_calib.py
        self.resolution = 3000. #par['resolution']
        self.nresln = 10.0 #par['nresln']

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

        norder = len(self.std)

        self.sens_dict = {}
        for iord in range(norder):
            #std_specobjs, std_header = load.load_specobjs(self.par['std_file'], order=iord)
            #embed(header='447 of fluxspec')
            #std_idx = flux.find_standard(std_specobjs)
            std = self.std[iord] #std_specobjs[std_idx]

            # THIS IS A CRAZY KLUDGE....
            if 'OPT_WAVE_GRID' in std._data.keys():
                # Deal with this
                embed(header='449 of fluxspec')
            #try:
            #    wavemask = std.optimal['WAVE_GRID'] > 0.0 #*units.AA
            #except KeyError:
            #    wavemask = std.optimal['WAVE'] > 1000.0 * units.AA
            #    this_wave = std.optimal['WAVE'][wavemask]
            #else:
            #    this_wave = std.optimal['WAVE_GRID'][wavemask]
            wavemask = std.OPT_WAVE > 1000.0#*units.AA)
            this_wave = std.OPT_WAVE[wavemask]

            counts, ivar = std.OPT_COUNTS[wavemask], std.OPT_COUNTS_IVAR[wavemask]
            sens_dict_iord = flux_calib.generate_sensfunc(this_wave, counts, ivar,
                                                          float(self.std_header['AIRMASS']),
                                                          self.std_header['EXPTIME'],
                                                          self.spectrograph.telescope['longitude'],
                                                          self.spectrograph.telescope['latitude'],
                                                          star_type=self.star_type,
                                                          star_mag=self.star_mag,
                                                          telluric=self.telluric, ra=self.std_ra, dec=self.std_dec,
                                                          resolution=self.resolution,
                                                          BALM_MASK_WID=self.BALM_MASK_WID, std_file=self.std_file,
                                                          poly_norder=self.poly_norder,
                                                          polycorrect=self.polycorrect, debug=self.debug)
            sens_dict_iord['ech_orderindx'] = iord
            sens_dict_iord['ech_order'] = std.ECH_ORDER
            self.sens_dict[str(iord)] = sens_dict_iord
        ## add some keys to be saved into primary header in masterframe
        for key in ['wave_max', 'exptime', 'airmass', 'std_file', 'std_ra', 'std_dec',
                    'std_name', 'cal_file']:
            try:
                self.sens_dict[key] = sens_dict_iord[key]
            except:
                pass
        self.sens_dict['meta'] = {}
        self.sens_dict['meta']['nslits'] = norder
        self.sens_dict['wave_min'] = self.sens_dict['0']['wave_min']

        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.sens_dict

    def flux_science(self, sci_file):
        """
        Flux the internal list of sci_specobjs

        Wrapper to flux.apply_sensfunc()

        Args:
            sci_file (str):

        """
        # Load
        self.load_objs(sci_file, std=False)
        # Loop on Echelle order
        norder = self.sens_dict['meta']['nslits']
        for iord in range(norder):
            sens_dict_iord = self.sens_dict[str(iord)]
            ech_order = sens_dict_iord['ech_order']
            for sci_obj in self.sci_specobjs:
                if sci_obj.ECH_ORDER == ech_order:
                    #flux_calib.apply_standard_sens(sci_obj, sens_dict_iord, float(self.sci_header['AIRMASS']),
                    #                               self.sci_header['EXPTIME'], extinct_correct=self.par['extinct_correct'],
                    #                               longitude=self.spectrograph.telescope['longitude'],
                    #                               latitude=self.spectrograph.telescope['latitude'])
                    sci_obj.apply_flux_calib(sens_dict_iord,
                                         self.sci_header['EXPTIME'],
                                         #telluric_correct=self.par['telluric_correct'],
                                         extinct_correct=self.par['extinct_correct'],
                                         longitude=self.spectrograph.telescope['longitude'],
                                         latitude=self.spectrograph.telescope['latitude'],
                                         airmass=float(self.sci_header['AIRMASS']))

        self.steps.append(inspect.stack()[0][3])


def instantiate_me(spectrograph, par, **kwargs):
    """
    Instantiate the FluxSpec subclass appropriate for the provided spectrograph.

    The class must be subclassed from FluxSpec.  See :class:`FluxSpec` for
    the description of the valid keyword arguments.

    Args:
      spectrograph : Spectrograph or str
        Name of the spectrograph, e.g. shane_kast_blue
        Used only to set settings for calls to the Class outside of PypeIt
        This includes extinction data..
      par: FluxCalib parset
      sens_file : str, optional
        Filename of a sensitivity function file to be input
    Returns:
        :class:`PypeIt`: One of the classes with :class:`PypeIt` as its
        base.
    """
    indx = [ c.__name__ == spectrograph.pypeline for c in FluxSpec.__subclasses__() ]
    if not np.any(indx):
        msgs.error('Pipeline {0} is not defined!'.format(spectrograph.pypeline))
    return FluxSpec.__subclasses__()[np.where(indx)[0][0]](spectrograph, par, **kwargs)
