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
    #
    'BOX_WAVE': dict(otype=np.ndarray, atype=float, desc='Optimal Wavelengths (Angstroms)'),
    'BOX_COUNTS': dict(otype=np.ndarray, atype=float, desc='Optimal flux (counts)'),
    'BOX_COUNTS_IVAR': dict(otype=np.ndarray, atype=float,
                            desc='Inverse variance of optimally extracted flux using modelivar image (counts^2)'),
    'BOX_COUNTS_SIG': dict(otype=np.ndarray, atype=float,
                           desc='Optimally extracted noise from IVAR (counts)'),
    'BOX_COUNTS_NIVAR': dict(otype=np.ndarray, atype=float,
                             desc='Optimally extracted noise variance, sky+read noise only (counts^2)'),
    'BOX_MASK': dict(otype=np.ndarray, atype=bool, desc='Mask for optimally extracted flux'),
    'BOX_COUNTS_SKY': dict(otype=np.ndarray, atype=float, desc='Optimally extracted sky (counts)'),
    'BOX_COUNTS_RN': dict(otype=np.ndarray, atype=float, desc='Optimally extracted RN squared (counts)'),
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
    'DET': dict(otype=(int,np.int64), desc='Detector number'),
    'PYPELINE': dict(otype=str, desc='Name of the PypeIt pipeline mode'),
    'OBJTYPE': dict(otype=str, desc='PypeIt type of object (standard, science)'),
    'SPAT_PIXPOS': dict(otype=(float,np.float32), desc='Spatial location of the trace on detector (pixel)'),
    #
    'SLITID': dict(otype=(int,np.int64), desc='Slit ID'),
    #
    'ECH_OBJID': dict(otype=(int,np.int64), desc='Echelle Object ID'),
    'ECH_ORDERINDX': dict(otype=(int,np.int64), desc='Order index.  Mainly for internal PypeIt usage'),
    'ECH_ORDER': dict(otype=(int,np.int64), desc='Physical echelle order'),
}


class SpecObj(object):
    """Class to handle object spectra from a single exposure
    One generates one of these Objects for each spectrum in the exposure. They are instantiated by the object
    finding routine, and then all spectral extraction information for the object are assigned as attributes

    Args:
        det (int): Detector number
        objtype (str, optional)
           Type of object ('unknown', 'standard', 'science')
        slitid (int, optional):
           Identifier for the slit (max=9999)

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
    def from_table(cls, table, indict=None):
        if table.meta['PYPELINE'] == 'MultiSlit':
            # Instantiate
            slf = cls(table.meta['PYPELINE'], table.meta['DET'], None,
                      slitid=table.meta['SLITID'], indict=indict)
        else:
            embed(header='112')
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
        # Name
        slf.set_name()
        # Return
        return slf

    def __init__(self, pypeline, det, slit_spat,
                 slitid=-1,
                 indict=None, objtype='unknown', orderindx=None):

        self._data = Table()

        # For copying the object
        if indict is not None:
            if '_SpecObj_initialised' in indict:
                indict.pop('_SpecObj_initialised')
            self.__dict__ = indict
        else:
            # set any attributes here - before initialisation
            # these remain as normal attributes
            # We may wish to eliminate *all* of these

            self.objid = 999
            self.name = None

            # Slit
            self.slit_spat = slit_spat  # (left, right)

            # Object finding
            self.spat_fracpos = None
            self.smash_peakflux = None
            self.smash_nsig = None
            self.maskwidth = None

            # Object profile
            self.prof_nsigma = None
            self.sign = 1.0
            self.min_spat = None
            self.max_spat = None

            # Trace
            self.trace_spec = None  # Only for debuggin, internal plotting

        # after initialisation, setting attributes is the same as setting an item
        self.__initialised = True

        # Initialize a few, if we aren't copying
        if indict is None:
            self.FLEX_SHIFT = 0.
            self.OBJTYPE = objtype
            self.DET = det
            self.PYPELINE = pypeline

            # pypeline specific
            if self.PYPELINE == 'MultiSlit':
                self.SLITID = slitid

        # Name
        self.set_name()

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

    '''
    @staticmethod
    def sobjs_key():
        """
        This function returns the dictionary that defines the mapping between specobjs attributes and the fits header
        cards

        Returns:
            dict:

        """
        sobjs_key_dict = dict(det='DET',
                              objid='OBJID',
                              slitid='SLITID',
                              ech_objid='ECHOBJID',
                              ech_orderindx='ECHOINDX',
                              ech_order='ECHORDER',
                              pypeline='PYPELINE')

        return sobjs_key_dict
    '''

    def set_name(self):
        """
        Generate a unique index for this spectrum based on the
        slit/order, its position and for multi-slit the detector.

        Sets self.name internally

        Returns:
            str: :attr:`self.name`

        """
        if 'Echelle' in self.PYPELINE:
            # ObjID
            self.name = naming_model['obj']
            if self.ech_objid is None:
                self.name += '----'
            else:
                self.name += '{:04d}'.format(self.ech_objid)
            self.name += '-'+naming_model['order']
            # Order
            if self.ech_orderindx is None:
                self.name += '----'
            else:
                self.name += '{:04d}'.format(self.ech_order)
        else:
            # Spat
            self.name = naming_model['spat']
            if 'SPAT_PIXPOS' not in self._data.meta.keys():
                self.name += '----'
            else:
                self.name += '{:04d}'.format(int(np.rint(self.SPAT_PIXPOS)))
            # Slit
            self.name += '-'+naming_model['slit']
            self.name += '{:04d}'.format(self.SLITID)
        # Detector
        sdet = parse.get_dnum(self.DET, prefix=False)
        self.name += '-{:s}{:s}'.format(naming_model['det'], sdet)
        # Return
        return self.name

    def copy(self):
        """
        Generate a copy of this object

        Returns:
            SpecObj

        """
        sobj_copy = SpecObj(self.PYPELINE, self.DET, self.slit_spat,
                            indict=self.__dict__.copy())
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
                          msgs.newline() + "{0:s}".format(str(self.name)))
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
                          + msgs.newline() + "{0}".format(str(self.name)))
                self[attr+'_WAVE'] *= vel_corr
                # Record
                self['VEL_TYPE'] = refframe
                self['VEL_CORR'] = vel_corr

    def to_xspec1d(self, extraction='OPT', fluxed=True):
        """
        Push the data in SpecObj into an XSpectrum1D object

        Args:
            extraction (str): Extraction method to convert

        Returns:
            linetools.spectra.xspectrum1d.XSpectrum1D:  Spectrum object

        """
        swave = extraction+'_WAVE'
        if swave not in self._data.keys():
            msgs.error("This object has not been extracted with extract={}.".format(extraction))
        # Fluxed?
        if fluxed:
            sflux = extraction+'_FLAM'
            ssig = extraction+'_FLAM_SIG'
        else:
            sflux = extraction+'_COUNTS'
            ssig = extraction+'_COUNTS_SIG'
        # Create
        xspec = xspectrum1d.XSpectrum1D.from_tuple((self[swave], self[sflux], self[ssig]))
        # Return
        return xspec

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



