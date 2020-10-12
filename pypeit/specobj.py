"""
Module for the SpecObj classes

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import copy
import inspect
from IPython import embed

import numpy as np

from scipy import interpolate

from astropy import units

from linetools.spectra import xspectrum1d

from pypeit import msgs
from pypeit.core import flexure
from pypeit.core import parse
from pypeit.core import flux_calib
from pypeit import utils
from pypeit import datamodel
from pypeit.images import detector_container

naming_model = {}
for skey in ['SPAT', 'SLIT', 'DET', 'SCI','OBJ', 'ORDER']:
    naming_model[skey.lower()] = skey

def det_hdu_prefix(det):
    return 'DET{:02d}-'.format(det)

class SpecObj(datamodel.DataContainer):
    """Class to handle object spectra from a single exposure
    One generates one of these Objects for each spectrum in the exposure. They are instantiated by the object
    finding routine, and then all spectral extraction information for the object are assigned as attributes

    Args:
        pypeline (str): Name of the PypeIt pypeline method
            Allowed options are:  MultiSlit, Echelle, IFU
        DET (int): Detector number
        copy_dict (dict, optional): Used to set the entire internal dict of the object.
            Only used in the copy() method so far.
        objtype (str, optional)
           Type of object ('unknown', 'standard', 'science')
        slitid (int, optional):
           Identifier for the slit (max=9999).
           Multislit and IFU
        specobj_dict (dict, optional):
           Uswed in the objfind() method of extract.py to Instantiate
        orderindx (int, optional):
           Running index for the order
        ech_order (int, optional):
           Physical order number

    Attributes:
        See datamodel and _init_internals()
    """
    version = '1.1.2'
    hdu_prefix = None

    datamodel = {'TRACE_SPAT': dict(otype=np.ndarray, atype=float,
                                    descr='Object trace along the spec (spatial pixel)'),
                 'FWHM': dict(otype=float, descr='Spatial FWHM of the object (pixels)'),
                 'FWHMFIT': dict(otype=np.ndarray,
                                 descr='Spatial FWHM across the detector (pixels)'),
                 'OPT_WAVE': dict(otype=np.ndarray, atype=float,
                                  descr='Optimal Wavelengths (Angstroms)'),
                 'OPT_FLAM': dict(otype=np.ndarray, atype=float,
                                  descr='Optimal flux (erg/s/cm^2/Ang)'),
                 'OPT_FLAM_SIG': dict(otype=np.ndarray, atype=float,
                                      descr='Optimal flux uncertainty (erg/s/cm^2/Ang)'),
                 'OPT_FLAM_IVAR': dict(otype=np.ndarray, atype=float,
                                       descr='Optimal flux inverse variance (erg/s/cm^2/Ang)^-2'),
                 'OPT_COUNTS': dict(otype=np.ndarray, atype=float, descr='Optimal flux (counts)'),
                 'OPT_COUNTS_IVAR': dict(otype=np.ndarray, atype=float,
                                         descr='Inverse variance of optimally extracted flux '
                                               'using modelivar image (counts^2)'),
                 'OPT_COUNTS_SIG': dict(otype=np.ndarray, atype=float,
                                        descr='Optimally extracted noise from IVAR (counts)'),
                 'OPT_COUNTS_NIVAR': dict(otype=np.ndarray, atype=float,
                                          descr='Optimally extracted noise variance, sky+read '
                                                'noise only (counts^2)'),
                 'OPT_MASK': dict(otype=np.ndarray, atype=np.bool_,
                                  descr='Mask for optimally extracted flux'),
                 'OPT_COUNTS_SKY': dict(otype=np.ndarray, atype=float,
                                        descr='Optimally extracted sky (counts)'),
                 'OPT_COUNTS_RN': dict(otype=np.ndarray, atype=float,
                                       descr='Optimally extracted RN squared (counts)'),
                 'OPT_FRAC_USE': dict(otype=np.ndarray, atype=float,
                                      descr='Fraction of pixels in the object profile subimage '
                                            'used for this extraction'),
                 'OPT_CHI2': dict(otype=np.ndarray, atype=float,
                                  descr='Reduced chi2 of the model fit for this spectral pixel'),
                 # TODO -- Confirm BOX_NPIX should be a float and not int!
                 'BOX_NPIX': dict(otype=np.ndarray, atype=float,
                                  descr='Number of pixels used for the boxcar extraction; can be '
                                        'fractional'),
                 'BOX_WAVE': dict(otype=np.ndarray, atype=float,
                                  descr='Boxcar Wavelengths (Angstroms)'),
                 'BOX_FLAM': dict(otype=np.ndarray, atype=float,
                                  descr='Boxcar flux (erg/s/cm^2/Ang)'),
                 'BOX_FLAM_SIG': dict(otype=np.ndarray, atype=float,
                                      descr='Boxcar flux uncertainty (erg/s/cm^2/Ang)'),
                 'BOX_FLAM_IVAR': dict(otype=np.ndarray, atype=float,
                                       descr='Boxcar flux inverse variance (erg/s/cm^2/Ang)^-2'),
                 'BOX_COUNTS': dict(otype=np.ndarray, atype=float, descr='Boxcar flux (counts)'),
                 'BOX_COUNTS_IVAR': dict(otype=np.ndarray, atype=float,
                                         descr='Inverse variance of optimally extracted flux '
                                               'using modelivar image (counts^2)'),
                 'BOX_COUNTS_SIG': dict(otype=np.ndarray, atype=float,
                                        descr='Boxcar extracted noise from IVAR (counts)'),
                 'BOX_COUNTS_NIVAR': dict(otype=np.ndarray, atype=float,
                                          descr='Boxcar extracted noise variance, sky+read noise '
                                                'only (counts^2)'),
                 'BOX_MASK': dict(otype=np.ndarray, atype=np.bool_,
                                  descr='Mask for optimally extracted flux'),
                 'BOX_COUNTS_SKY': dict(otype=np.ndarray, atype=float,
                                        descr='Boxcar extracted sky (counts)'),
                 'BOX_COUNTS_RN': dict(otype=np.ndarray, atype=float,
                                       descr='Boxcar extracted RN squared (counts)'),
                 'BOX_FRAC_USE': dict(otype=np.ndarray, atype=float,
                                      descr='Fraction of pixels in the object profile subimage '
                                            'used for this extraction'),
                 'BOX_CHI2': dict(otype=np.ndarray, atype=float,
                                  descr='Reduced chi2 of the model fit for this spectral pixel'),
                 'BOX_RADIUS': dict(otype=float, descr='Size of boxcar radius (pixels)'),
                 #
                 'FLEX_SHIFT_GLOBAL': dict(otype=float, descr='Global shift of the spectrum to correct for spectral'
                                                              'flexure (pixels). This is based on the sky spectrum at'
                                                              'the center of the slit'),
                 'FLEX_SHIFT_LOCAL': dict(otype=float, descr='Local shift of the spectrum to correct for spectral'
                                                             'flexure (pixels). This should be a small correction to'
                                                             'the global value, and is based on the sky spectrum'
                                                             'extracted near the object'),
                 'FLEX_SHIFT_TOTAL': dict(otype=float, descr='Total shift of the spectrum to correct for spectral'
                                                             'flexure (pixels). This is the sum of the global and'
                                                             'local FLEX_SHIFT'),
                 'VEL_TYPE': dict(otype=str, descr='Type of heliocentric correction (if any)'),
                 'VEL_CORR': dict(otype=float,
                                  descr='Relativistic velocity correction for wavelengths'),
                 # Detector
                 # TODO: Why are both det and detector attributes, isn't det in detector?
                 'DET': dict(otype=(int, np.integer), descr='Detector number'),
                 'DETECTOR': dict(otype=detector_container.DetectorContainer,
                                  descr='Detector DataContainer'),
                 #
                 'PYPELINE': dict(otype=str, descr='Name of the PypeIt pipeline mode'),
                 'OBJTYPE': dict(otype=str, descr='PypeIt type of object (standard, science)'),
                 'SPAT_PIXPOS': dict(otype=(float, np.floating),
                                     descr='Spatial location of the trace on detector (pixel)'),
                 'SPAT_FRACPOS': dict(otype=(float, np.floating),
                                      descr='Fractional location of the object on the slit'),
                 # Slit and Object
                 'SLITID': dict(otype=(int, np.integer),
                                descr='PypeIt slit ID. Increasing from left to right on '
                                      'detector. Zero based.'),
                 'OBJID': dict(otype=(int, np.integer),
                               descr='Object ID for multislit data. Each object is given an index '
                                     'for the slit it appears increasing from from left to right. '
                                     'These are one based.'),
                 'NAME': dict(otype=str, descr='Name of the object following the naming model'),
                 'RA': dict(otype=float, descr='Right Ascension (J2000) decimal degree'),
                 'DEC': dict(otype=float, descr='Declination (J2000) decimal degree'),
                 'MASKDEF_ID': dict(otype=(int, np.integer), descr='Slitmask definition ID'),
                 'MASKDEF_OBJNAME': dict(otype=str, descr='Name of the object from the slitmask definition'),
                 #
                 'ECH_OBJID': dict(otype=(int, np.integer),
                                   descr='Object ID for echelle data. Each object is given an '
                                         'index in the order it appears increasing from from left '
                                         'to right. These are one based.'),
                 'ECH_ORDERINDX': dict(otype=(int, np.integer),
                                       descr='Order indx, analogous to SLITID for echelle. '
                                             'Zero based.'),
                 'ECH_FRACPOS': dict(otype=(float, np.floating),
                                     descr='Synced echelle fractional location of the object on '
                                           'the slit'),
                 'ECH_ORDER': dict(otype=(int, np.integer), descr='Physical echelle order'),
                 'ECH_NAME': dict(otype=str,
                                  descr='Name of the object for echelle data. Same as NAME above '
                                        'but order numbers are omitted giving a unique name per '
                                        'object.')}

    def __init__(self, PYPELINE, DET, OBJTYPE='unknown',
                 SLITID=None, ECH_ORDER=None, ECH_ORDERINDX=None):

        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=_d)

        self.FLEX_SHIFT_GLOBAL = 0.
        self.FLEX_SHIFT_LOCAL = 0.
        self.FLEX_SHIFT_TOTAL = 0.

        # Name
        self.set_name()

    def _init_internals(self):
        # Object finding
        self.smash_peakflux = None
        self.smash_nsig = None
        self.maskwidth = None

        # Hand
        self.hand_extract_flag = False
        self.hand_extract_spec = None
        self.hand_extract_spat = None
        self.hand_extract_det = None
        self.hand_extract_fwhm = None

        # Object profile
        self.prof_nsigma = None
        self.sign = 1.0
        self.min_spat = None
        self.max_spat = None

        # Trace
        self.trace_spec = None  # Only for debuggin, internal plotting

        # Echelle
        self.ech_frac_was_fit = None #
        self.ech_snr = None #

    def _bundle(self, **kwargs):
        """
        Over-ride DataContainer._bundle() to deal with DETECTOR

        Args:
            kwargs:
                Passed to DataContainer._bundle()

        Returns:
            list:

        """
        _d = super(SpecObj, self)._bundle(**kwargs)
        # Move DetectorContainer into its own HDU
        if _d[0]['DETECTOR'] is not None:
            _d.append(dict(detector=_d[0].pop('DETECTOR')))
        # Return
        return _d


    def to_hdu(self, hdr=None, add_primary=False, primary_hdr=None,
               limit_hdus=None, force_to_bintbl=True):
        """
        Over-ride :func:`pypeit.datamodel.DataContainer.to_hdu` to force to
        a BinTableHDU

        See that func for Args and Returns
        """
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = dict([(k,values[k]) for k in args[1:]])
        # Force
        _d['force_to_bintbl'] = True
        # Do it
        return super(SpecObj, self).to_hdu(**_d)

    @property
    def slit_order(self):
        if self.PYPELINE == 'Echelle':
            return self.ECH_ORDER
        elif self.PYPELINE == 'MultiSlit':
            return self.SLITID
        elif self.PYPELINE == 'IFU':
            return self.SLITID
        else:
            msgs.error("Bad PYPELINE")


    @property
    def slit_orderindx(self):
        if self.PYPELINE == 'Echelle':
            return self.ECH_ORDERINDX
        elif self.PYPELINE == 'MultiSlit':
            return self.SLITID
        elif self.PYPELINE == 'IFU':
            return self.SLITID
        else:
            msgs.error("Bad PYPELINE")

    def set_name(self):
        """
        Generate a unique index for this spectrum based on the
        slit/order, its position and for multi-slit the detector.

        Multi-slit

            Each object is named by its:
             - spatial position (pixel number) on the reduced image [SPAT]
             - the slit number based on SPAT center of the slit or SlitMask ID [SLIT]
             - the detector number [DET]

            For example::

                SPAT0176-SLIT0185-DET01

        Echelle

        Returns:
            str:

        """
        if 'Echelle' in self.PYPELINE:
            # ObjID
            name = naming_model['obj']
            ech_name = naming_model['obj']
            if self['ECH_FRACPOS'] is None:
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
            if self['SPAT_PIXPOS'] is None:
                name += '----'
            else:
                name += '{:04d}'.format(int(np.rint(self.SPAT_PIXPOS)))
            # Slit
            name += '-'+naming_model['slit']
            name += '{:04d}'.format(self.SLITID)
            sdet = parse.get_dnum(self.DET, prefix=False)
            name += '-{:s}{:s}'.format(naming_model['det'], sdet)
            self.NAME = name
        elif 'IFU' in self.PYPELINE:
            # Spat
            name = naming_model['spat']
            if self['SPAT_PIXPOS'] is None:
                name += '----'
            else:
                name += '{:04d}'.format(int(np.rint(self.SPAT_PIXPOS)))
            # Slit
            name += '-' + naming_model['slit']
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
            :class:`SpecObj`:

        """
        # Return
        return copy.deepcopy(self)

    def apply_spectral_flexure(self, shift, sky_spec):
        """
        Apply interpolation with the flexure dict

        Args:
            shift (float):
                additive spectral flexure in pixels
            sky_spec (`xspectrum1d.XSpectrum1D`_):
                Sky Spectrum

        Returns:
            xspectrum1d.XSpectrum1D:  New sky spectrum (mainly for QA)
        """
        # Simple interpolation to apply
        # Apply
        for attr in ['BOX', 'OPT']:
            if self[attr+'_WAVE'] is not None:
                msgs.info("Applying flexure correction to {0:s} extraction for object:".format(attr) +
                          msgs.newline() + "{0:s}".format(str(self.NAME)))
                self[attr+'_WAVE'] = flexure.flexure_interp(shift, self[attr+'_WAVE']).copy()
        # Shift sky spec too
        twave = flexure.flexure_interp(shift, sky_spec.wavelength.value) * units.AA
        new_sky = xspectrum1d.XSpectrum1D.from_tuple((twave, sky_spec.flux))
        # Save - since flexure may have been applied/calculated twice, this needs to be additive
        self.update_flex_shift(shift, flex_type='local')
        # Return
        return new_sky

    def update_flex_shift(self, shift, flex_type='local'):
        """Store the total spectral flexure shift in pixels

        Args:
            shift (float):
                additive spectral flexure in pixels
        """
        if flex_type == 'global':
            self.FLEX_SHIFT_GLOBAL = shift
        elif flex_type == 'local':
            self.FLEX_SHIFT_LOCAL = shift
        else:
            msgs.error("Spectral flexure type must be 'global' or 'local' only")
        # Now update the total flexure
        self.FLEX_SHIFT_TOTAL += shift

    # TODO This should be a wrapper calling a core algorithm.
    def apply_flux_calib(self, wave_sens, sensfunc, exptime, telluric=None, extinct_correct=False,
                         airmass=None, longitude=None, latitude=None, extrap_sens=False):
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
            extrap_sens (bool, optional):
                Extrapolate the sensitivity function (instead of crashing out)

        """
        # Loop on extraction modes
        for attr in ['BOX', 'OPT']:
            if self[attr+'_WAVE'] is None:
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
                if extrap_sens:
                    sensfunc_obs[wave_mask] = interpolate.interp1d(wave_sens, sensfunc, bounds_error=False)(wave[wave_mask])
                    msgs.warn("Your data extends beyond the bounds of your sensfunc. You should be adjusting the par['sensfunc']['extrap_blu'] and/or par['sensfunc']['extrap_red'] to extrapolate further and recreate your sensfunc. But we are extrapolating per your direction. Good luck!")
                else:
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
            if self[attr+'_WAVE'] is not None:
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
        if self[swave] is None:
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
        Push the data in :class:`SpecObj` into an XSpectrum1D object


        Returns:
            linetools.spectra.xspectrum1d.XSpectrum1D:  Spectrum object

        """
        wave, flux, ivar, _ = self.to_arrays(**kwargs)
        sig = np.sqrt(utils.inverse(ivar))
        # Create
        return xspectrum1d.XSpectrum1D.from_tuple((wave, flux, sig))

