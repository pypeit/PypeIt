# Module for echelle fluxing
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import linetools
import os
import json
import matplotlib.pyplot as plt

from astropy import units
from astropy.io import fits
from astropy.table import Table

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

class EchFluxSpec(masterframe.MasterFrame):
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
                 std_specobjs=None, std_header=None, spectrograph=None,
                 telluric=False, setup=None, master_dir=None, mode=None,
                 star_type=None, star_mag = None, BALM_MASK_WID = 5.0, nresln = None, debug=False):

        # Load standard files
        std_spectro = None
        self.std_spec1d_file = std_spec1d_file
        # Need to unwrap these (sometimes)..
        self.std_specobjs = std_specobjs
        self.std_header = std_header
        if self.std_spec1d_file is not None:
            self.std_specobjs, self.std_header = load.ech_load_specobj(self.std_spec1d_file)
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
            self.sci_specobjs, self.sci_header = load.ech_load_specobj(self.sci_spec1d_file)
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
        # Echelle key
        self.norder = None
        # ToDo: Need to parse the following in
        self.star_type = star_type
        self.star_mag = star_mag
        self.BALM_MASK_WID = BALM_MASK_WID
        self.nresln = nresln
        self.debug = debug


    def load_master(self, filename, force=False):

        # Does the master file exist?
        if not os.path.isfile(filename):
            # msgs.warn("No Master frame found of type {:s}: {:s}".format(self.frametype, filename))
            msgs.warn("No Master frame found of {:s}".format(filename))
            if force:
                msgs.error("Crashing out because reduce-masters-force=True:" + msgs.newline() + filename)
            return None
        else:
            # msgs.info("Loading a pre-existing master calibration frame of type: {:}".format(self.frametype) + " from filename: {:}".format(filename))
            msgs.info("Loading a pre-existing master calibration frame of SENSFUNC from filename: {:}".format(filename))

            hdu = fits.open(filename)
            norder = hdu[0].header['NORDER']
            sens_dicts = {}
            for iord in range(norder):
                head = hdu[iord + 1].header
                tbl = hdu['SENSFUNC-ORDER{0:04}'.format(iord)].data
                sens_dict = {}
                sens_dict['wave'] = tbl['WAVE']
                sens_dict['sensfunc'] = tbl['SENSFUNC']
                for key in ['wave_min', 'wave_max', 'exptime', 'airmass', 'std_file', 'std_ra', 'std_dec',
                            'std_name', 'cal_file', 'ech_orderindx']:
                    try:
                        sens_dict[key] = head[key.upper()]
                    except:
                        pass
                sens_dicts[str(iord)] = sens_dict
            sens_dicts['norder'] = norder
            return sens_dicts

    def save_master(self, sens_dicts, outfile=None):
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

        norder = self.sens_dict['norder']

        # Do it
        prihdu = fits.PrimaryHDU()
        hdus = [prihdu]

        for iord in range(norder):
            sens_dict_iord = self.sens_dict[str(iord)]
            cols = []
            cols += [fits.Column(array=sens_dict_iord['wave'], name=str('WAVE'), format=sens_dict_iord['wave'].dtype)]
            cols += [
                fits.Column(array=sens_dict_iord['sensfunc'], name=str('SENSFUNC'), format=sens_dict_iord['sensfunc'].dtype)]
            # Finish
            coldefs = fits.ColDefs(cols)
            tbhdu = fits.BinTableHDU.from_columns(coldefs)
            tbhdu.name = 'SENSFUNC-ORDER{0:04}'.format(iord)
            # Add critical keys from sens_dict to header
            for key in ['wave_min', 'wave_max', 'exptime', 'airmass', 'std_file', 'std_ra',
                        'std_dec', 'std_name', 'cal_file', 'ech_orderindx']:
                try:
                    tbhdu.header[key.upper()] = sens_dict_iord[key].value
                except AttributeError:
                    tbhdu.header[key.upper()] = sens_dict_iord[key]
                except KeyError:
                    pass  # Will not require all of these
            hdus += [tbhdu]

        # Add critical keys from sens_dict to primary header
        for key in ['exptime', 'airmass', 'std_file', 'std_ra',
                    'std_dec', 'std_name', 'cal_file']:
            try:
                prihdu.header[key.upper()] = sens_dict_iord[key].value
            except AttributeError:
                prihdu.header[key.upper()] = sens_dict_iord[key]
            except KeyError:
                pass  # Will not require all of these
        prihdu.header['NORDER'] = norder

        # Finish
        hdulist = fits.HDUList(hdus)
        hdulist.writeto(outfile, overwrite=True)

        # Finish
        msgs.info("Wrote sensfunc to MasterFrame: {:s}".format(outfile))

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
        #if self.std is None:
        #    msgs.warn('First identify the star first (with find_standard).')
        #    return None
        if self.std_header is None:
            msgs.warn('First set std_header with a dict-like object holding RA, DEC, '
                      'AIRMASS, EXPTIME.')
            return None
        if self.norder is None:
            ext_final = fits.getheader(self.std_spec1d_file, -1)
            self.norder = ext_final['ORDER'] + 1

        self.sens_dict = {}
        for iord in range(self.norder):
            std_specobjs, std_header = load.ech_load_specobj(self.std_spec1d_file, order=iord)
            std_idx = flux.find_standard(std_specobjs)
            std = std_specobjs[std_idx]
            wavemask = std.boxcar['WAVE'] > 1000.0 * units.AA
            wave, counts, ivar = std.boxcar['WAVE'][wavemask], std.boxcar['COUNTS'][wavemask], \
                                 std.boxcar['COUNTS_IVAR'][wavemask]
            sens_dict_iord = flux.generate_sensfunc(wave, counts, ivar, std_header['AIRMASS'], std_header['EXPTIME'],
                                               self.spectrograph, star_type=self.star_type, star_mag=self.star_mag,
                                               telluric=self.telluric, ra=self.std_ra, dec=self.std_dec,
                                               BALM_MASK_WID=self.BALM_MASK_WID,
                                               nresln=self.nresln, std_file=self.std_file, debug=self.debug)
            sens_dict_iord['ech_orderindx'] = iord
            self.sens_dict[str(iord)] = sens_dict_iord
        self.sens_dict['norder'] = self.norder

        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.sens_dict

    def flux_science(self):
        """
        Flux the internal list of sci_specobjs

        Wrapper to flux.apply_sensfunc()

        Returns
        -------

        """
        norder = self.norder
        for iord in range(norder):
            sens_dict = self.sens_dict[str(iord)]
            for sci_obj in self.sci_specobjs:
                if sci_obj.ech_orderindx == iord:
                    flux.apply_sensfunc(sci_obj, self.sens_dict, self.sci_header['AIRMASS'],
                                        self.sci_header['EXPTIME'], self.spectrograph)
        self.steps.append(inspect.stack()[0][3])

    def show_sensfunc(self):
        """
        Plot the sensitivity function
        """
        if self.sens_dict is None:
            msgs.warn("You need to generate the sensfunc first!")
            return None
        # Generate from model
        norder = self.sens_dict['norder']
        for iord in range(norder):
            sens_dict_iord = self.sens_dict[str(iord)]
            plt.plot(sens_dict_iord['wave'],sens_dict_iord['sensfunc'])
        plt.xlabel('Wavelength [ang]')
        plt.ylabel('Sensfunc')
        plt.ylim([0.,100.0])
        plt.show()

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

