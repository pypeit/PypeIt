""" Module for the SpecObjs and SpecObj classes
"""
import copy
from IPython import embed

import numpy as np

from scipy import interpolate

from astropy import units
from astropy.table import Table

from linetools.spectra import xspectrum1d

from pypeit import msgs
from pypeit.core import parse
from pypeit.core import flux_calib
from pypeit import utils

naming_model = {}
for skey in ['SPAT', 'SLIT', 'DET', 'SCI','OBJ', 'ORDER']:
    naming_model[skey.lower()] = skey

# Data model -- Put here to be able to reach it without instantiating the class
#  These are outward facing items, i.e. items that the user will receive and use.
#  These are upper case to distinguish them from internal attributes
data_model = {
    'TRACE_SPAT': dict(otype=np.ndarray, atype=float, desc='Object trace along the spec (spatial pixel)'),
    'FWHM': dict(otype=float, desc='Spatial FWHM of the object (pixels)'),
    'FWHMFIT': dict(otype=np.ndarray, desc='Spatial FWHM across the detector (pixels)'),
    'OPT_WAVE': dict(otype=np.ndarray, atype=float, desc='Optimal Wavelengths (Angstroms)'),
    'OPT_FLAM': dict(otype=np.ndarray, atype=float, desc='Optimal flux (erg/s/cm^2/Ang)'),
    'OPT_FLAM_SIG': dict(otype=np.ndarray, atype=float, desc='Optimal flux uncertainty (erg/s/cm^2/Ang)'),
    'OPT_FLAM_IVAR': dict(otype=np.ndarray, atype=float, desc='Optimal flux inverse variance (erg/s/cm^2/Ang)^-2'),
    'OPT_COUNTS': dict(otype=np.ndarray, atype=float, desc='Optimal flux (counts)'),
    'OPT_COUNTS_IVAR': dict(otype=np.ndarray, atype=float,
                            desc='Inverse variance of optimally extracted flux using modelivar image (counts^2)'),
    'OPT_COUNTS_SIG': dict(otype=np.ndarray, atype=float,
                             desc='Optimally extracted noise from IVAR (counts)'),
    'OPT_COUNTS_NIVAR': dict(otype=np.ndarray, atype=float,
                             desc='Optimally extracted noise variance, sky+read noise only (counts^2)'),
    'OPT_MASK': dict(otype=np.ndarray, atype=bool, desc='Mask for optimally extracted flux'),
    'OPT_COUNTS_SKY': dict(otype=np.ndarray, atype=float, desc='Optimally extracted sky (counts)'),
    'OPT_COUNTS_RN': dict(otype=np.ndarray, atype=float, desc='Optimally extracted RN squared (counts)'),
    'OPT_FRAC_USE': dict(otype=np.ndarray, atype=float,
                         desc='Fraction of pixels in the object profile subimage used for this extraction'),
    'OPT_CHI2': dict(otype=np.ndarray, atype=float,
                     desc='Reduced chi2 of the model fit for this spectral pixel'),
    'BOX_WAVE': dict(otype=np.ndarray, atype=float, desc='Boxcar Wavelengths (Angstroms)'),
    'BOX_FLAM': dict(otype=np.ndarray, atype=float, desc='Boxcar flux (erg/s/cm^2/Ang)'),
    'BOX_FLAM_SIG': dict(otype=np.ndarray, atype=float, desc='Boxcar flux uncertainty (erg/s/cm^2/Ang)'),
    'BOX_FLAM_IVAR': dict(otype=np.ndarray, atype=float, desc='Boxcar flux inverse variance (erg/s/cm^2/Ang)^-2'),
    'BOX_COUNTS': dict(otype=np.ndarray, atype=float, desc='Boxcar flux (counts)'),
    'BOX_COUNTS_IVAR': dict(otype=np.ndarray, atype=float,
                            desc='Inverse variance of optimally extracted flux using modelivar image (counts^2)'),
    'BOX_COUNTS_SIG': dict(otype=np.ndarray, atype=float,
                           desc='Boxcar extracted noise from IVAR (counts)'),
    'BOX_COUNTS_NIVAR': dict(otype=np.ndarray, atype=float,
                             desc='Boxcar extracted noise variance, sky+read noise only (counts^2)'),
    'BOX_MASK': dict(otype=np.ndarray, atype=bool, desc='Mask for optimally extracted flux'),
    'BOX_COUNTS_SKY': dict(otype=np.ndarray, atype=float, desc='Boxcar extracted sky (counts)'),
    'BOX_COUNTS_RN': dict(otype=np.ndarray, atype=float, desc='Boxcar extracted RN squared (counts)'),
    'BOX_FRAC_USE': dict(otype=np.ndarray, atype=float,
                         desc='Fraction of pixels in the object profile subimage used for this extraction'),
    'BOX_CHI2': dict(otype=np.ndarray, atype=float,
                     desc='Reduced chi2 of the model fit for this spectral pixel'),
    'BOX_RADIUS': dict(otype=float, desc='Size of boxcar radius (pixels)'),
    #
    'FLEX_SHIFT': dict(otype=float, desc='Shift of the spectrum to correct for flexure (pixels)'),
    'VEL_TYPE': dict(otype=str, desc='Type of heliocentric correction (if any)'),
    'VEL_CORR': dict(otype=float, desc='Relativistic velocity correction for wavelengths'),
    #
    'DET': dict(otype=(int,np.integer), desc='Detector number'),
    'PYPELINE': dict(otype=str, desc='Name of the PypeIt pipeline mode'),
    'OBJTYPE': dict(otype=str, desc='PypeIt type of object (standard, science)'),
    'SPAT_PIXPOS': dict(otype=(float,np.floating), desc='Spatial location of the trace on detector (pixel)'),
    'SPAT_FRACPOS': dict(otype=(float,np.floating), desc='Fractional location of the object on the slit'),
    #
    'SLITID': dict(otype=(int,np.integer), desc='Slit ID. Increasing from left to right on detector. Zero based.'),
    'OBJID': dict(otype=(int, np.integer), desc='Object ID for multislit data. Each object is given an index for the slit '
                                                  'it appears increasing from from left to right. These are one based.'),
    'NAME': dict(otype=str, desc='Name of the object following the naming model'),
    #
    'ECH_OBJID': dict(otype=(int, np.integer),
                      desc='Object ID for echelle data. Each object is given an index in the order '
                           'it appears increasing from from left to right. These are one based.'),
    'ECH_ORDERINDX': dict(otype=(int, np.integer), desc='Order indx, analogous to SLITID for echelle. Zero based.'),
    'ECH_FRACPOS': dict(otype=(float,np.floating), desc='Synced echelle fractional location of the object on the slit'),
    'ECH_ORDER': dict(otype=(int, np.integer), desc='Physical echelle order'),
    'ECH_NAME': dict(otype=str, desc='Name of the object for echelle data. Same as NAME above but order numbers are '
                                     'omitted giving a unique name per object.')
}


class SpecObj(object):
    """Class to handle object spectra from a single exposure
    One generates one of these Objects for each spectrum in the exposure. They are instantiated by the object
    finding routine, and then all spectral extraction information for the object are assigned as attributes

    Args:
        pypeline (str): Name of the PypeIt pypeline method
            Allowed options are:  MultiSlit, Echelle
        det (int): Detector number
        copy_dict (dict, optional): Used to set the entire internal dict of the object.
            Only used in the copy() method so far.
        objtype (str, optional)
           Type of object ('unknown', 'standard', 'science')
        slitid (int, optional):
           Identifier for the slit (max=9999).
           Multislit only
        specobj_dict (dict, optional):
           Uswed in the objfind() method of extract.py to Instantiate
        orderindx (int, optional):
           Running index for the order
        ech_order (int, optional):
           Physical order number

    Attributes:
        slitcen (float): Center of slit in fraction of total (trimmed) detector size at ypos
        objid (int): Identifier for the object (max=999)
        flex_shift (float): Flexure correction in pixels

    Extraction dict's
        'WAVE' : wave_opt  # Optimally extracted wavelengths
        'COUNTS' : flux_opt  # Optimally extracted flux
        'COUNTS_IVAR' : mivar_opt  # Inverse variance of optimally extracted flux using modelivar image
        'COUNTS_NIVAR' : nivar_opt  # Optimally extracted noise variance (sky + read noise) only
        'MASK' : mask_opt  # Mask for optimally extracted flux
        'COUNTS_SKY' : sky_opt  # Optimally extracted sky
        'COUNTS_RN' : rn_opt  # Square root of optimally extracted read noise squared
        'FRAC_USE' : frac_use  # Fraction of pixels in the object profile subimage used for this extraction
        'CHI2' : chi2  # Reduced chi2 of the model fit for this spectral pixel
    """
    @classmethod
    def from_table(cls, table, copy_dict=None):
        if table.meta['PYPELINE'] == 'MultiSlit':
            # Instantiate
            slf = cls(table.meta['PYPELINE'], table.meta['DET'],
                      slitid=table.meta['SLITID'], copy_dict=copy_dict)
        else:
            slf = cls(table.meta['PYPELINE'], table.meta['DET'],
                      copy_dict=copy_dict, ech_order=table.meta['ECH_ORDER'], orderindx=table.meta['ECH_ORDERINDX'])
        # Pop a few that land in standard FITS header
        # Loop me -- Do this to deal with checking the data model
        for key in table.keys():
            setattr(slf, key, table[key].data)
        for key in table.meta.keys():
            # Skip ones that can appear in FITS header
            if key in ['EXTNAME']:
                continue
            #
            setattr(slf, key, table.meta[key])
        return slf

    # TODO: JFH I really don't like this copy_dict implementation and I don't know why you added it. This should simply be
    # done via the copy method as it was before.
    def __init__(self, pypeline, det, objtype='unknown',
                 copy_dict=None,
                 slitid=None,
                 ech_order=None,
                 orderindx=None,
                 specobj_dict=None):

        self._data = Table()

        # For copying the object
        if copy_dict is not None:
            if '_SpecObj_initialised' in copy_dict:
                copy_dict.pop('_SpecObj_initialised')
            self.__dict__ = copy_dict
        else:
            # set any attributes here - before initialisation
            # these remain as normal attributes
            # We may wish to eliminate *all* of these

            # Object finding
            self.smash_peakflux = None
            self.smash_nsig = None
            self.maskwidth = None
            self.hand_extract_flag = False

            # Object profile
            self.prof_nsigma = None
            self.sign = 1.0
            self.min_spat = None
            self.max_spat = None

            # Trace
            self.trace_spec = None  # Only for debuggin, internal plotting

            # Echelle
            #self.ech_orderindx = None #': dict(otype=(int,np.int64), desc='Order index.  Mainly for internal PypeIt usage'),
            #self.ech_objid = None # 'ECH_OBJID': dict(otype=(int,np.int64), desc='Echelle Object ID'),
            self.ech_frac_was_fit = None #
            self.ech_snr = None #

        # after initialisation, setting attributes is the same as setting an item
        self.__initialised = True

        # TODO: JFH This is not very elegant and error prone. Perhaps we should loop over the copy dict keys
        # and populate whatever keys are actually there. For instance, if the user does not pass in a copy dict,
        # and forgets to pass in ech_order, the code is going to crash becuase ech_order cannot be None.

        # Initialize a few, if we aren't copying
        if copy_dict is None:
            self.DET = det
            if specobj_dict is not None:
                self.PYPELINE = specobj_dict['pypeline']
                self.OBJTYPE = specobj_dict['objtype']
                if self.PYPELINE == 'MultiSlit':
                    self.SLITID = specobj_dict['slitid']
                elif self.PYPELINE == 'Echelle':
                    self.ECH_ORDER = specobj_dict['order']
                    self.ECH_ORDERINDX = specobj_dict['orderindx']
            else:
                self.PYPELINE = pypeline
                self.OBJTYPE = objtype
                # pypeline specific
                if self.PYPELINE == 'MultiSlit':
                    self.SLITID = slitid
                elif self.PYPELINE == 'Echelle':
                    self.ECH_ORDER = ech_order
                    self.ECH_ORDERINDX = orderindx
                else:
                    msgs.error("Uh oh")

            self.FLEX_SHIFT = 0.

        # Name
        self.set_name()

    @property
    def slit_order(self):
        if self.PYPELINE == 'Echelle':
            return self.ECH_ORDER
        elif self.PYPELINE == 'MultiSlit':
            return self.SLITID
        else:
            msgs.error("Uh oh")


    @property
    def slit_orderindx(self):
        if self.PYPELINE == 'Echelle':
            return self.ECH_ORDERINDX
        elif self.PYPELINE == 'MultiSlit':
            return self.SLITID
        else:
            msgs.error("Uh oh")

    def __getattr__(self, item):
        """Maps values to attributes.
        Only called if there *isn't* an attribute with this name
        """
        try:
            return self.__getitem__(item)
        except KeyError:
            raise AttributeError(item)

    def __setattr__(self, item, value):

        if not '_SpecObj__initialised' in self.__dict__:  # this test allows attributes to be set in the __init__ method
            return dict.__setattr__(self, item, value)
        elif item in self.__dict__:       # any normal attributes are handled normally
            dict.__setattr__(self, item, value)
        else:
            self.__setitem__(item, value)

    def __setitem__(self, item, value):
        if item not in data_model.keys():
            raise IOError("Cannot set {} attribute.  It is not in the data model".format(item))
        if not isinstance(value, data_model[item]['otype']):
            print("Wrong data type for attribute: {}".format(item))
            print("Allowed type(s) are: {}".format(data_model[item]['otype']))
            raise IOError("Try again")
        if isinstance(value, np.ndarray):
            self._data[item] = value
        else:
            self._data.meta[item] = value

    def __getitem__(self, item):
        if item in self._data.keys():
            return self._data[item].data
        elif item in self._data.meta.keys():
            return self._data.meta[item]
        else:
            raise KeyError

    # TODO: JFH Can we change this functio name to data_model, and have a separate function called
    # keys which returns sobjs._dict__.keys() which is the list of attributes of the actual object
    # excluding the data model. Or even better, use keys for both, but append the two lists with one set uppercase and
    # the other set lowercase
    def keys(self):
        """
        Simple method to return the keys of _data

        Returns:
            list
        """
        return self._data.keys()

    # TODO JFH Please describe the naming model somewhere in this module.
    def set_name(self):
        """
        Generate a unique index for this spectrum based on the
        slit/order, its position and for multi-slit the detector.

        Sets name

        Returns:
            str

        """
        if 'Echelle' in self.PYPELINE:
            # ObjID
            name = naming_model['obj']
            ech_name = naming_model['obj']
            if 'ECH_FRACPOS' not in self._data.meta.keys():
                name += '----'
            else:
                # JFH TODO Why not just write it out with the decimal place. That is clearer than this??
                name += '{:04d}'.format(int(np.rint(1000*self.ECH_FRACPOS)))
                ech_name += '{:04d}'.format(int(np.rint(1000*self.ECH_FRACPOS)))
            sdet = parse.get_dnum(self.DET, prefix=False)
            name += '-{:s}{:s}'.format(naming_model['det'], sdet)
            ech_name += '-{:s}{:s}'.format(naming_model['det'], sdet)
            # Order number
            name += '-'+naming_model['order']
            name += '{:04d}'.format(self.ECH_ORDER)
            self.ECH_NAME = ech_name
            self.NAME = name
        elif 'MultiSlit' in self.PYPELINE:
            # Spat
            name = naming_model['spat']
            if 'SPAT_PIXPOS' not in self._data.meta.keys():
                name += '----'
            else:
                name += '{:04d}'.format(int(np.rint(self.SPAT_PIXPOS)))
            # Slit
            name += '-'+naming_model['slit']
            name += '{:04d}'.format(self.SLITID)
            sdet = parse.get_dnum(self.DET, prefix=False)
            name += '-{:s}{:s}'.format(naming_model['det'], sdet)
            self.NAME = name
        else:
            msgs.error("Bad PYPELINE")


    def copy(self):
        """
        Generate a copy of this object

        Returns:
            SpecObj

        """
        #sobj_copy = SpecObj(self.PYPELINE, self.DET,
        #                    copy_dict=self.__dict__.copy())
        # JFH Without doing a deepcopy here, this does not make a true copy. It is somehow using pointers, and so changing the
        # copy changes the original object which wreaks havoc. That is why it was deepcopy before (I think).
        sobj_copy = SpecObj(self.PYPELINE, self.DET,
                            copy_dict=copy.deepcopy(self.__dict__))
        # Return
        return sobj_copy

    def flexure_interp(self, sky_wave, fdict):
        """
        Apply interpolation with the flexure dict

        Args:
            sky_wave (np.ndarray): Wavelengths of the extracted sky
            fdict (dict): Holds the various flexure items

        Returns:
            np.ndarray:  New sky spectrum (mainly for QA)

        """
        # Simple interpolation to apply
        npix = len(sky_wave)
        x = np.linspace(0., 1., npix)
        # Apply
        for attr in ['BOX', 'OPT']:
            if attr+'_WAVE' in self._data.keys():
                msgs.info("Applying flexure correction to {0:s} extraction for object:".format(attr) +
                          msgs.newline() + "{0:s}".format(str(self.NAME)))
                f = interpolate.interp1d(x, sky_wave, bounds_error=False, fill_value="extrapolate")
                self[attr+'_WAVE'] = f(x + fdict['shift'] / (npix - 1)) * units.AA
        # Shift sky spec too
        cut_sky = fdict['sky_spec']
        x = np.linspace(0., 1., cut_sky.npix)
        f = interpolate.interp1d(x, cut_sky.wavelength.value, bounds_error=False, fill_value="extrapolate")
        twave = f(x + fdict['shift'] / (cut_sky.npix - 1)) * units.AA
        new_sky = xspectrum1d.XSpectrum1D.from_tuple((twave, cut_sky.flux))
        # Save
        self.FLEX_SHIFT = fdict['shift']
        # Return
        return new_sky

    # TODO This should be a wrapper calling a core algorithm.
    def apply_flux_calib(self, wave_sens, sensfunc, exptime, telluric=None, extinct_correct=False,
                         airmass=None, longitude=None, latitude=None):
        """
        Apply a sensitivity function to our spectrum

        FLAM, FLAM_SIG, and FLAM_IVAR are generated

        Args:
            sens_dict (dict):
                Sens Function dict
            exptime (float):
            telluric_correct:
            extinct_correct:
            airmass (float, optional):
            longitude (float, optional):
                longitude in degree for observatory
            latitude:
                latitude in degree for observatory
                Used for extinction correction

        """
        # Loop on extraction modes
        for attr in ['BOX', 'OPT']:
            if attr+'_WAVE' not in self._data.keys():
                continue
            msgs.info("Fluxing {:s} extraction for:".format(attr) + msgs.newline() + "{}".format(self))

            wave = self[attr+'_WAVE']
            # Interpolate the sensitivity function onto the wavelength grid of the data

            # TODO Telluric corrections via this method are deprecated
            # Did the user request a telluric correction?
            if telluric is not None:
                # This assumes there is a separate telluric key in this dict.
                msgs.info('Applying telluric correction')
                sensfunc = sensfunc * (telluric > 1e-10) / (telluric + (telluric < 1e-10))

            sensfunc_obs = np.zeros_like(wave)
            wave_mask = wave > 1.0  # filter out masked regions or bad wavelengths
            try:
                sensfunc_obs[wave_mask] = interpolate.interp1d(wave_sens, sensfunc, bounds_error=True)(wave[wave_mask])
            except ValueError:
                msgs.error("Your data extends beyond the bounds of your sensfunc. " + msgs.newline() +
                           "Adjust the par['sensfunc']['extrap_blu'] and/or par['sensfunc']['extrap_red'] to extrapolate "
                           "further and recreate your sensfunc.")

            if extinct_correct:
                if longitude is None or latitude is None:
                    msgs.error('You must specify longitude and latitude if we are extinction correcting')
                # Apply Extinction if optical bands
                msgs.info("Applying extinction correction")
                msgs.warn("Extinction correction applyed only if the spectra covers <10000Ang.")
                extinct = flux_calib.load_extinction_data(longitude, latitude)
                ext_corr = flux_calib.extinction_correction(wave * units.AA, airmass, extinct)
                senstot = sensfunc_obs * ext_corr
            else:
                senstot = sensfunc_obs.copy()

            flam = self[attr+'_COUNTS'] * senstot / exptime
            flam_sig = (senstot / exptime) / (np.sqrt(self[attr+'_COUNTS_IVAR']))
            flam_var = self[attr+'_COUNTS_IVAR'] / (senstot / exptime) ** 2

            # Mask bad pixels
            msgs.info(" Masking bad pixels")
            msk = np.zeros_like(senstot).astype(bool)
            msk[senstot <= 0.] = True
            msk[self[attr+'_COUNTS_IVAR'] <= 0.] = True
            flam[msk] = 0.
            flam_sig[msk] = 0.
            flam_var[msk] = 0.
            # TODO JFH We need to update the mask here. I think we need a mask for the counts and a mask for the flam,
            # since they can in principle be different. We are masking bad sensfunc locations.

            # Finish
            self[attr+'_FLAM'] = flam
            self[attr+'_FLAM_SIG'] = flam_sig
            self[attr+'_FLAM_IVAR'] = flam_var


    def apply_helio(self, vel_corr, refframe):
        """
        Apply a heliocentric correction

        Wavelength arrays are modified in place

        Args:
            vel_corr (float):
            refframe (str):

        """
        # Apply
        for attr in ['BOX', 'OPT']:
            if attr+'_WAVE' in self._data.keys():
                msgs.info('Applying {0} correction to '.format(refframe)
                          + '{0} extraction for object:'.format(attr)
                          + msgs.newline() + "{0}".format(str(self.NAME)))
                self[attr+'_WAVE'] *= vel_corr
                # Record
                self['VEL_TYPE'] = refframe
                self['VEL_CORR'] = vel_corr

    def to_arrays(self, extraction='OPT', fluxed=True):
        """

        Args:
            extraction (str): Extraction method to convert
            fluxed:

        Returns:
            tuple: wave, flux, ivar, mask arrays

        """
        swave = extraction+'_WAVE'
        smask = extraction+'_MASK'
        if swave not in self._data.keys():
            msgs.error("This object has not been extracted with extract={}.".format(extraction))
        # Fluxed?
        if fluxed:
            sflux = extraction+'_FLAM'
            sivar = extraction+'_FLAM_IVAR'
        else:
            sflux = extraction+'_COUNTS'
            sivar = extraction+'_COUNTS_IVAR'
        # Return
        return self[swave], self[sflux], self[sivar], self[smask]

    def to_xspec1d(self, **kwargs):
        """
        Push the data in SpecObj into an XSpectrum1D object


        Returns:
            linetools.spectra.xspectrum1d.XSpectrum1D:  Spectrum object

        """
        wave, flux, ivar, _ = self.to_arrays(**kwargs)
        sig = np.sqrt(utils.inverse(ivar))
        # Create
        xspec = xspectrum1d.XSpectrum1D.from_tuple((wave, flux, sig))
        # Return
        return xspec

    # TODO JFH: This method does not work
    def show(self, extraction='optimal'):
        """
        Show the spectrum by converting it to a XSpectrum1D object

        Args:
            extraction (str): Extraction option 'optimal' or 'boxcar'

        Returns:

        """
        extract = getattr(self, extraction)
        # Generate an XSpec
        xspec = self.to_xspec1d(extraction=extraction)
        if xspec is None:
            return
        xspec.plot(xspec=True)


    def __repr__(self):
        txt = '<{:s}: {:s}'.format(self.__class__.__name__, self.NAME)
        txt += '>'
        return txt
