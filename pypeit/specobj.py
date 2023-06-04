"""
Module for the SpecObj classes

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import copy
import inspect
from IPython import embed

import numpy as np

from astropy import units

from linetools.spectra import xspectrum1d

from pypeit import msgs
from pypeit.core import flexure
from pypeit.core import flux_calib
from pypeit import utils
from pypeit import datamodel
from pypeit.images.detector_container import DetectorContainer
from pypeit.images.mosaic import Mosaic


class SpecObj(datamodel.DataContainer):
    """
    Class to handle single spectra from a single exposure.

    One generates one of these objects for each spectrum in the exposure. They
    are instantiated by the object finding routine, and then all spectral
    extraction information for the object are assigned as attributes

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_specobj.rst

    Args:
        PYPELINE (:obj:`str`):
            Name of the ``PypeIt`` pipeline method.  Allowed options are
            MultiSlit, Echelle, or IFU.
        DET (:obj:`str`):
            The name of the detector or mosaic from which the spectrum was
            extracted.  For example, DET01.
        OBJTYPE (:obj:`str`, optional):
            Object type.  For example: 'unknown', 'standard', 'science'.
        SLITID (:obj:`int`, optional):
            For multislit and IFU reductions, this is an identifier for the slit
            (max=9999).
        ECH_ORDER (:obj:`int`, optional):
            Physical order number.
        ECH_ORDERINDX (:obj:`int`, optional):
            Running index for the order.
    """

    version = '1.1.8'
    """
    Current datamodel version number.
    """

    datamodel = {'TRACE_SPAT': dict(otype=np.ndarray, atype=float,
                                    descr='Object trace along the spec (spatial pixel)'),
                 'FWHM': dict(otype=float, descr='Spatial FWHM of the object (pixels)'),
                 'FWHMFIT': dict(otype=np.ndarray,
                                 descr='Spatial FWHM across the detector (pixels)'),
                 'smash_peakflux': dict(otype=float,
                                        descr='Peak value of the spectral direction collapsed spatial profile'),
                 'smash_snr': dict(otype=float,
                                        descr='Peak S/N ratio of the spectral direction collapsed patial profile'),
                 'OPT_WAVE': dict(otype=np.ndarray, atype=float,
                                  descr='Optimal Wavelengths in vacuum (Angstroms)'),
                 'OPT_FLAM': dict(otype=np.ndarray, atype=float,
                                  descr='Optimal flux (1e-17 erg/s/cm^2/Ang)'),
                 'OPT_FLAM_SIG': dict(otype=np.ndarray, atype=float,
                                      descr='Optimal flux uncertainty (1e-17 erg/s/cm^2/Ang)'),
                 'OPT_FLAM_IVAR': dict(otype=np.ndarray, atype=float,
                                       descr='Optimal flux inverse variance (1e-17 erg/s/cm^2/Ang)^-2'),
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
                                  descr='Mask for optimally extracted flux. True=good'),
                 'OPT_COUNTS_SKY': dict(otype=np.ndarray, atype=float,
                                        descr='Optimally extracted sky (counts)'),
                 'OPT_COUNTS_SIG_DET': dict(otype=np.ndarray, atype=float,
                                            descr='Optimally extracted detector noise (counts)'),
                 'OPT_FRAC_USE': dict(otype=np.ndarray, atype=float,
                                      descr='Fraction of pixels in the object profile subimage '
                                            'used for this extraction'),
                 'OPT_CHI2': dict(otype=np.ndarray, atype=float,
                                  descr='Reduced chi2 of the model fit for this spectral pixel'),
                 'BOX_NPIX': dict(otype=np.ndarray, atype=float,
                                  descr='Number of pixels used for the boxcar extraction; can be '
                                        'fractional'),
                 'BOX_WAVE': dict(otype=np.ndarray, atype=float,
                                  descr='Boxcar Wavelengths in vacuum (Angstroms)'),
                 'BOX_FLAM': dict(otype=np.ndarray, atype=float,
                                  descr='Boxcar flux (erg/s/cm^2/Ang)'),
                 'BOX_FLAM_SIG': dict(otype=np.ndarray, atype=float,
                                      descr='Boxcar flux uncertainty (1e-17 erg/s/cm^2/Ang)'),
                 'BOX_FLAM_IVAR': dict(otype=np.ndarray, atype=float,
                                       descr='Boxcar flux inverse variance (1e-17 erg/s/cm^2/Ang)^-2'),
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
                                  descr='Mask for boxcar extracted flux. True=good'),
                 'BOX_COUNTS_SKY': dict(otype=np.ndarray, atype=float,
                                        descr='Boxcar extracted sky (counts)'),
                 'BOX_COUNTS_SIG_DET': dict(otype=np.ndarray, atype=float,
                                            descr='Boxcar extracted detector noise (counts)'),
                 'BOX_FRAC_USE': dict(otype=np.ndarray, atype=float,
                                      descr='Fraction of pixels in the object profile subimage '
                                            'used for this extraction'),
                 'BOX_CHI2': dict(otype=np.ndarray, atype=float,
                                  descr='Reduced chi2 of the model fit for this spectral pixel'),
                 'BOX_RADIUS': dict(otype=float, descr='Size of boxcar radius (pixels)'),
                 'S2N': dict(otype=float, descr='Median signal to noise ratio of the extracted spectrum'
                                                '(OPT if available, otherwise BOX)'),
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
                 # TODO: Change this to DETNAME?
                 # NOTE: DET (or DETNAME) is needed in the case when DETECTOR is None.
                 'DET': dict(otype=str,
                             descr='A string identifier for the reduced detector or mosaic.'),
                 'DETECTOR': dict(otype=(DetectorContainer, Mosaic),
                                  descr='Object with the detector or mosaic metadata'),
                 'PYPELINE': dict(otype=str, descr='Name of the PypeIt pipeline mode'),
                 # TODO: It's unclear if OBJTYPE has to be one among a set of
                 # specific keywords.
                 'OBJTYPE': dict(otype=str, descr='Object type (e.g., standard, science)'),
                 'SPAT_PIXPOS': dict(otype=(float, np.floating),
                                     descr='Spatial location of the trace on detector (pixel) at half-way'),
                 'SPAT_FRACPOS': dict(otype=(float, np.floating),
                                      descr='Fractional location of the object on the slit'),
                 'trace_spec': dict(otype=np.ndarray, atype=(int,np.integer),
                                      descr='Array of pixels along the spectral direction'),
                 'maskwidth': dict(otype=(float, np.floating),
                                      descr='Size (in units of fwhm) of the region used for local sky subtraction'),
                 # Slit and Object
                 'WAVE_RMS': dict(otype=(float, np.floating),
                                     descr='RMS (pix) for the wavelength solution for this slit.'),
                 'SLITID': dict(otype=(int, np.integer),
                                descr='PypeIt slit ID (aka SPAT_ID).'),
                 'OBJID': dict(otype=(int, np.integer),
                               descr='Object ID for multislit data. Each object is given an index '
                                     'for the slit it appears increasing from from left to right. '
                                     'These are one based.'),
                 'NAME': dict(otype=str, descr='Name of the object following the naming model'),
                 'RA': dict(otype=float, descr='Right Ascension (J2000) decimal degree'),
                 'DEC': dict(otype=float, descr='Declination (J2000) decimal degree'),
                 'MASKDEF_ID': dict(otype=(int, np.integer), descr='Slitmask definition ID'),
                 'MASKDEF_OBJNAME': dict(otype=str, descr='Name of the object from the slitmask definition'),
                 'MASKDEF_OBJMAG': dict(otype=float, descr='Magnitude of the object from the slitmask definition'),
                 'MASKDEF_OBJMAG_BAND': dict(otype=str, descr='Magnitude band of the object from the slitmask '
                                                              'definition'),
                 'MASKDEF_EXTRACT': dict(otype=bool, descr='Boolean indicating if this is a forced extraction '
                                                           'at the expected location from slitmask design. '),
                 'hand_extract_flag': dict(otype=bool, descr='Boolean indicating if this is a forced extraction '
                                                             'at the location provided by the user. '),
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
    """
    Defines the current datmodel.
    """

    internals = [# Object finding
                 'smash_peakflux',
                 'smash_snr',
                 # Hand
                 'hand_extract_flag',
                 'hand_extract_spec',
                 'hand_extract_spat',
                 'hand_extract_det',
                 'hand_extract_fwhm',
                 # Object profile
                 'prof_nsigma',
                 'sign',
                 'min_spat',
                 'max_spat',
                 # Echelle
                 'ech_frac_was_fit',
                 'ech_snr'
                ]

    def __init__(self, PYPELINE, DET, OBJTYPE='unknown',
                 SLITID=None, ECH_ORDER=None, ECH_ORDERINDX=None):

        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = dict([(k,values[k]) for k in args[1:]])
        super().__init__(d=_d)

        # Initialize internal values that are not None
        self.hand_extract_flag = False
        self.sign = 1.

        self.FLEX_SHIFT_GLOBAL = 0.
        self.FLEX_SHIFT_LOCAL = 0.
        self.FLEX_SHIFT_TOTAL = 0.

        # Name
        self.set_name()

    @classmethod
    def from_arrays(cls, PYPELINE:str, wave:np.ndarray, counts:np.ndarray, ivar:np.ndarray,
                    mode='OPT', DET='DET01', SLITID=0, **kwargs):
        # Instantiate
        slf = cls(PYPELINE, DET, SLITID=SLITID)
        # Add in arrays
        for item, attr in zip([wave, counts, ivar], ['_WAVE', '_COUNTS', '_COUNTS_IVAR']):
            setattr(slf, mode+attr, item.astype(float))
        # Mask. Watch out for places where ivar is infinite due to a divide by 0
        slf[mode+'_MASK'] = (slf[mode+'_COUNTS_IVAR'] > 0.) & np.isfinite(slf[mode+'_COUNTS_IVAR'])
        return slf

    def _validate(self):
        """
        Validate the object.
        """
        pypelines = ['MultiSlit', 'IFU', 'Echelle']
        if self.PYPELINE not in pypelines:
            msgs.error(f'{self.PYPELINE} is not a known pipeline procedure.  Options are: '
                       f"{', '.join(pypelines)}")

    def _bundle(self, **kwargs):
        """
        Override base class to handle inclusion of
        :class:`~pypeit.images.detector_container.DetectorContainer` or
        :class:`~pypeit.images.mosaic.Mosaic` objects for each spectrum.

        Args:
            kwargs (:obj:`dict`):
                Passed directly to the base class method.

        Returns:
            :obj:`list`: List of dictionaries used to construct fits extensions.
        """
        # Use the base class for most of the data model
        _d = super()._bundle(**kwargs)

        # Separate out the detector into its own HDU
        if _d[0]['DETECTOR'] is not None:
            _d.append(dict(detector=_d[0].pop('DETECTOR')))

        return _d

    def to_hdu(self, **kwargs):
        """
        Override the base class function to force the main datamodel attributes
        to be written to an `astropy.io.fits.BinTableHDU`_ object.  This is
        identical to the base class method except ``force_to_bintbl`` is always
        set to True.
        """
        if 'force_to_bintbl' in kwargs and not kwargs['force_to_bintbl']:
            msgs.warn(f'Writing a {self.__class__.__name__} always requires force_to_bintbl=True')
            del kwargs['force_to_bintbl']
        return super().to_hdu(force_to_bintbl=True, **kwargs)

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

    @property
    def mnx_wave(self):
        """Return min, max wavelength of the spectrum
        Uses OPT_WAVE if present and then BOX_WAVE

        Returns:
            tuple: min, max (float)
        """
        mnx = (0., 0.)
        for pref in ['OPT', 'BOX']:
            if self[pref+'_WAVE'] is not None:
                mnx = self[pref+'_WAVE'].min(), self[pref+'_WAVE'].max() 
            if mnx[0] != 0.:
                break
        return mnx

    def med_s2n(self):
        """Return median S/N of the spectrum
        Uses OPT_COUNTS if present and then BOX_COUNTS

        Returns:
            float
        """
        SN = 0.
        for pref in ['OPT', 'BOX']:
            if self[pref+'_COUNTS'] is not None and np.any(self[pref+'_MASK']):
                SN = np.median(self[pref+'_COUNTS'][self[pref+'_MASK']] *
                               np.sqrt(self[pref+'_COUNTS_IVAR'][self[pref+'_MASK']]))
            if SN != 0.:
                break
        return SN

    def set_name(self):
        """
        Construct the ``PypeIt`` name for this object.

        The ``PypeIt`` name depends on the type of data being processed:

            - For multislit and IFU data, the name is
              ``SPATnnnn-SLITmmmm-{DET}``, where ``nnnn`` is the nearest integer
              pixel in the spatial direction (at the spectral midpoint) where
              the object was extracted, ``mmmm`` is the slit identification
              number, and ``{DET}`` is the string identifier for the detector or
              mosaic.

            - For echelle data, the name is ``OBJnnnn-{DET}-ORDERoooo``, where
              ``nnnn`` is 1000 times the fractional position along the spatial
              direction rounded to the nearest integer, ``{DET}`` is the string
              identifier for the detector or mosaic, and ``oooo`` is the echelle
              order.

        For echelle data, this modifies both :attr:`ECH_NAME` and :attr:`NAME`;
        otherwise, only the latter is set.
        """
        naming_model = {}
        for skey in ['SPAT', 'SLIT', 'SCI', 'OBJ', 'ORDER']:
            naming_model[skey.lower()] = skey

        if self.PYPELINE == 'Echelle':
            # ObjID
            name = naming_model['obj']
            ech_name = naming_model['obj']
            if self['ECH_FRACPOS'] is None:
                name += '----'
            else:
                # JFH TODO Why not just write it out with the decimal place. That is clearer than this??
                name += '{:04d}'.format(int(np.rint(1000*self.ECH_FRACPOS)))
                ech_name += '{:04d}'.format(int(np.rint(1000*self.ECH_FRACPOS)))
            name += f'-{self.DET}'
            ech_name += f'-{self.DET}'
            # Order number
            name += '-'+naming_model['order']
            name += '{:04d}'.format(self.ECH_ORDER)
            self.ECH_NAME = ech_name
            self.NAME = name
        elif self.PYPELINE in ['MultiSlit', 'IFU']:
            # Spat
            name = naming_model['spat']
            if self['SPAT_PIXPOS'] is None:
                name += '----'
            else:
                name += '{:04d}'.format(int(np.rint(self.SPAT_PIXPOS)))
            # Slit
            name += '-'+naming_model['slit']
            name += '{:04d}'.format(self.SLITID)
            name += f'-{self.DET}'
            self.NAME = name
        else:
            msgs.error(f'{self.PYPELINE} is not an understood pipeline.')

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
            sky_spec (`linetools.spectra.xspectrum1d.XSpectrum1D`_):
                Sky Spectrum

        Returns:
            `linetools.spectra.xspectrum1d.XSpectrum1D`_: New sky
            spectrum (mainly for QA)
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
    def apply_flux_calib(self, wave_zp, zeropoint, exptime, tellmodel=None, extinct_correct=False,
                         airmass=None, longitude=None, latitude=None, extinctfilepar=None, extrap_sens=False):
        """
        Apply a sensitivity function to our spectrum

        FLAM, FLAM_SIG, and FLAM_IVAR are generated

        Args:
            wave_zp (float array)
                Zeropoint wavelength array
            zeropoint (float array):
                zeropoint array
            exptime (float):
                Exposure time
            tellmodel:
                Telluric correction
            extinct_correct:
                If True, extinction correct
            airmass (float, optional):
                Airmass
            longitude (float, optional):
                longitude in degree for observatory
            latitude:
                latitude in degree for observatory
                Used for extinction correction
            extinctfilepar (str):
                [sensfunc][UVIS][extinct_file] parameter
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
            sens_factor = flux_calib.get_sensfunc_factor(
                wave, wave_zp, zeropoint, exptime, tellmodel=tellmodel, extinct_correct=extinct_correct, airmass=airmass,
                longitude=longitude, latitude=latitude, extinctfilepar=extinctfilepar, extrap_sens=extrap_sens)
            flam = self[attr+'_COUNTS']*sens_factor
            flam_sig = sens_factor/np.sqrt(self[attr+'_COUNTS_IVAR'])
            flam_ivar = self[attr+'_COUNTS_IVAR']/sens_factor**2

            # Mask bad pixels
            msgs.info(" Masking bad pixels")
            msk = np.zeros_like(sens_factor).astype(bool)
            msk[sens_factor <= 0.] = True
            msk[self[attr+'_COUNTS_IVAR'] <= 0.] = True
            flam[msk] = 0.
            flam_sig[msk] = 0.
            flam_ivar[msk] = 0.
            # TODO JFH We need to update the mask here. I think we need a mask for the counts and a mask for the flam,
            # since they can in principle be different. We are masking bad sensfunc locations.

            # Finish
            self[attr+'_FLAM'] = flam
            self[attr+'_FLAM_SIG'] = flam_sig
            self[attr+'_FLAM_IVAR'] = flam_ivar


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
        Convert spectrum into np.ndarray arrays

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

    def to_xspec1d(self, masked=False, **kwargs):
        """
        Create an `XSpectrum1D <linetools.spectra.xspectrum1d.XSpectrum1D>`_
        using this spectrum.

        Args:
            masked (:obj:`bool`, optional):
                If True, only unmasked data are included.
            kwargs (:obj:`dict`, optional):
                Passed directly to :func:`to_arrays`.

        Returns:
            `linetools.spectra.xspectrum1d.XSpectrum1D`_: Spectrum object
        """
        wave, flux, ivar, gpm = self.to_arrays(**kwargs)
        sig = np.sqrt(utils.inverse(ivar))
        if masked:
            wave = wave[gpm]
            flux = flux[gpm]
            sig = sig[gpm]
        # Create
        return xspectrum1d.XSpectrum1D.from_tuple((wave, flux, sig))

    def ready_for_extraction(self):
        """ Simple method to check all the items are filled
        and ready for skysub and extraction.

        Returns:
            bool: True if all checks have passed
        """
        required = ['TRACE_SPAT', 'SPAT_PIXPOS', 'SPAT_FRACPOS',
            'trace_spec', 'OBJID', 'FWHM', 'maskwidth', 'NAME',
            'smash_peakflux', 'smash_snr',
            'SLITID', 'DET', 'PYPELINE', 'OBJTYPE']
        if 'Echelle' in self.PYPELINE:
            required += ['ECH_NAME']

        passed = True
        for key in required:
            if self[key] is None:
                msgs.warn("Item {} is missing from SpecObj. Failing vette".format(key))
                msgs.warn('{}'.format(self))
                passed = False
        #
        return passed

    def __repr__(self):
        """ Over-ride print representation

        Returns:
            str: Basics of the Data Container
        """
        repr = '<{:s}: '.format(self.__class__.__name__)
        # Image
        rdict = {}
        for attr in self.datamodel.keys():
            if hasattr(self, attr) and getattr(self, attr) is not None:
                # Special ones
                if attr in ['DET', 'SLITID', 'SPAT_PIXPOS', 'NAME', 'RA', 
                            'DEC', 'MASKDEF_ID', 'MASKDEF_OBJNAME', 'MASKDEF_EXTRACT',
                            'MASKDEF_OBJMAG', 'MASKDEF_OBJMAG_BAND']:
                    rdict[attr] = getattr(self,attr)
                else:
                    rdict[attr] = True
            else:
                rdict[attr] = False
        #repr += ' items={}'.format(rdict)
        repr += ' items={'
        for key in rdict.keys():
            if rdict[key] is not False:
                repr += '{}: {}\n'.format(key, rdict[key])
        return repr + '>'

    def has_opt_ext(self, fluxed=False):
        """
        Check that all the values of the optimal extraction exist

        Args:
            fluxed (:obj:`bool`, optional):
                Check that the flux-calibrated data exist.

        Returns:
            :obj:`bool`: True if all OPT values are available
        """
        flx = 'FLAM' if fluxed else 'COUNTS'
        return self.check_populated(['OPT_WAVE', f'OPT_{flx}', f'OPT_{flx}_IVAR', 'OPT_MASK'])

    def get_opt_ext(self, fluxed=False):
        """
        Return the optimal extraction values

        Args:
            fluxed (:obj:`bool`, optional):
                Return the flux-calibrated data.

        Returns:
            :obj:`tuple`: OPT_WAVE, OPT_COUNTS, OPT_COUNTS_IVAR, OPT_MASK
            attributes of SpecObj
        """
        if fluxed:
            # TODO: This does not check if the fluxed data exists!
            return self.OPT_WAVE, self.OPT_FLAM, self.OPT_FLAM_IVAR, self.OPT_MASK
        return self.OPT_WAVE, self.OPT_COUNTS, self.OPT_COUNTS_IVAR, self.OPT_MASK

    def has_box_ext(self, fluxed=False):
        """
        Check that all the values of the boxcar extraction exist

        Args:
            fluxed (:obj:`bool`, optional):
                Check that the flux-calibrated data exist.

        Returns:
            :obj:`bool`: True if all BOX values are available
        """
        flx = 'FLAM' if fluxed else 'COUNTS'
        return self.check_populated(['BOX_WAVE', f'BOX_{flx}', f'BOX_{flx}_IVAR', 'BOX_MASK'])

    def get_box_ext(self, fluxed=False):
        """
        Return the boxcar extraction values

        Args:
            fluxed (:obj:`bool`, optional):
                Return the flux-calibrated data.

        Returns:
            :obj:`tuple`: BOX_WAVE, BOX_COUNTS, BOX_COUNTS_IVAR, BOX_MASK
            attributes of SpecObj
        """
        if fluxed:
            # TODO: This does not check if the fluxed data exists!
            return self.BOX_WAVE, self.BOX_FLAM, self.BOX_FLAM_IVAR, self.BOX_MASK
        return self.BOX_WAVE, self.BOX_COUNTS, self.BOX_COUNTS_IVAR, self.BOX_MASK

    def best_ext_match(self, extract=None, fluxed=True):
        """
        Determine the extraction and calibration type that best matches a user
        request.

        Precedence is given to the requested extraction and calibration types.
        Beyond that, optimal extraction takes precedence over boxcar extraction,
        and flux-calibrated data take precedence over uncalibrated counts.

        Args:
            extract (:obj:`str`, optional):
                The extraction used to produce the spectrum.  Must be either
                None, ``'BOX'`` (for a boxcar extraction), or ``'OPT'`` for
                optimal extraction.  If None, the optimal extraction will be
                returned, if it exists, otherwise the boxcar extraction will be
                returned.
            fluxed (:obj:`bool`, optional):
                If True, return the flux calibrated spectrum, if it exists.  If
                the flux calibration hasn't been performed or ``fluxed=False``,
                the spectrum is returned in counts.

        Returns:
            :obj:`tuple`: The adjusted extraction type (``'BOX'`` or ``'OPT'``)
            and the adjusted calibration type (True for flux-calibrated, False
            for uncalibrated counts).
        """
        # If not set, prefer the optimal extraction over the boxcar one.
        _extract = 'OPT' if extract is None else extract
        if _extract not in ['OPT', 'BOX']:
            msgs.error(f'Extraction type ({_extract}) not understood; must be OPT or BOX.')
        if _extract == 'OPT':
            if self.has_opt_ext(fluxed=fluxed):
                return 'OPT', fluxed
            # If we make it here, expect that fluxed was True.  Try flipping it.
            if self.has_opt_ext(fluxed=False):
                return 'OPT', False

        # If we make it here, either extract was BOX from the start, or none of
        # the optimal extraction options were available
        if self.has_box_ext(fluxed=fluxed):
            return 'BOX', fluxed
        # If we make it here, assume fluxed was True.  Try flipping it.
        if self.has_box_ext(fluxed=False):
            return 'BOX', False
        # If we make it here, we've got a problem!
        msgs.error('Unable to find a relevant set of data!')

