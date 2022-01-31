"""
Module for LBT/LUCI specific methods.

.. include:: ../include/links.rst
"""

from IPython import embed

import numpy as np

from astropy.coordinates import SkyCoord

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container


class LBTLUCISpectrograph(spectrograph.Spectrograph):
    """
    Class to handle LBT/LUCI specific code

    The current implementation of the LUCI spectrograph sets the plate scale
    appropriate for the N1.8 camera. Other cameras, such as the N30 camera
    used in adative optics mode is not implemented.

    The provided default reduction parameters have been tailored and tested
    for LUCI1 and LUCI2 using the G200 low resolution grating and the zJspec
    and HKspec filters.

    """
    ndet = 1
    telescope = telescopes.LBTTelescopePar()


    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)

        # Target meta
        # CORE
        self.meta['ra'] = dict(ext=0, card='OBJRA')
        self.meta['dec'] = dict(ext=0, card='OBJDEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')

        # Instrument meta
        # CORE
        self.meta['dispname'] = dict(ext=0, card='GRATNAME')
        self.meta['decker'] = dict(ext=0, card='MASKID')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')
        # additional
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')
        self.meta['filter1'] = dict(ext=0, card='FILTER2')
        # populating dispangle with the grating wavelength
        self.meta['dispangle'] = dict(ext=0, card='GRATWLEN', rtol=0.1)

        # Observation meta
        # CORE
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        # additioanl
        self.meta['readmode'] = dict(ext=0, card='READMODE')
        self.meta['savemode'] = dict(ext=0, card='SAVEMODE')
        self.meta['dit'] = dict(ext=0, card='DIT')
        self.meta['ndit'] = dict(ext=0, card='NDIT')

        # Define compound meta
        self.meta['idname'] = dict(card=None, compound=True)

    # TODO: Deal with isot time here.
    def compound_meta(self, headarr, meta_key):
        """
        Methods to generate metadata requiring interpretation of the header
        data, instead of simply reading the value of a header card.

        Args:
            headarr (:obj:`list`):
                List of `astropy.io.fits.Header`_ objects.
            meta_key (:obj:`str`):
                Metadata keyword to construct.

        Returns:
            object: Metadata value read from the header(s).
        """
        # Populate the idname based on the header information of LUCI
        # This is an implicit way of pre-typing without adding too many
        # variables to the self.meta.

        if meta_key == 'idname':
            targetname = (headarr[0].get('OBJECT'))
            dispname = (headarr[0].get('GRATNAME'))
            calib_unit = (headarr[0].get('CALIB'))
            filter1 = (headarr[0].get('FILTER1'))
            filter2 = (headarr[0].get('FILTER2'))
            lamp1 = (headarr[0].get('LAMP1'))
            lamp2 = (headarr[0].get('LAMP2'))
            lamp3 = (headarr[0].get('LAMP3'))
            lamp4 = (headarr[0].get('LAMP4'))
            lamp5 = (headarr[0].get('LAMP5'))
            lamp6 = (headarr[0].get('LAMP6'))

            # object frame -> will be typed as science
            # This currently includes sky flats, science and standard images
            # We will guess standards using the beginning of their names.
            if ((dispname != 'Mirror') and
                (calib_unit == False) and
                (lamp1 == False) and
                (lamp2 == False) and
                (lamp3 == False) and
                (lamp4 == False) and
                (lamp5 == False) and
                (lamp6 == False)):
                if (targetname[:3] == 'HIP' or
                    targetname[:2] == 'HD' or
                    targetname[:5] == 'Feige'):
                    return 'standard'
                else:
                    return 'object'
            # flat frames -> will be typed as pixelflat, trace
            elif ((calib_unit == True) and
                  ((lamp4 == True) or
                   (lamp5 == True) or
                   (lamp6 == True))):
                return 'flat'
            # arcs -> will be typed as arc, tilt
            elif ((dispname != 'Mirror') and
                  (calib_unit == True) and
                  ((lamp1 == True) or
                   (lamp2 == True) or
                   (lamp3 == True))):
                return 'arc'
            # pixelflat off -> will be typed as bias
            elif ((dispname != 'Mirror') and
                (calib_unit == True) and
                (lamp1 == False) and
                (lamp2 == False) and
                (lamp3 == False) and
                (lamp4 == False) and
                (lamp5 == False) and
                (lamp6 == False) and
                (filter1 != 'blind') and
                (filter2 != 'blind')):
                return 'flat_off'
            # darks -> will not be typed currently
            elif ((filter1 == 'blind') or
                  (filter2 == 'blind')):
                return 'dark'

        msgs.error("Not ready for this compound meta")

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        This list is used by :class:`~pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            :obj:`list`: List of keywords of data pulled from file headers
            and used to constuct the :class:`~pypeit.metadata.PypeItMetaData`
            object.
        """
        return ['decker', 'dispname', 'dispangle', 'filter1']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard ``PypeIt`` file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        pypeit_keys = super().pypeit_file_keys()
        pypeit_keys += ['calib', 'comb_id', 'bkg_id', 'idname']
        return pypeit_keys

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.

        Args:
            ftype (:obj:`str`):
                Type of frame to check. Must be a valid frame type; see
                frame-type :ref:`frame_type_defs`.
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for one or more frames to check.
            exprng (:obj:`list`, optional):
                Range in the allowed exposure time for a frame of type
                ``ftype``. See
                :func:`pypeit.core.framematch.check_frame_exptime`.

        Returns:
            `numpy.ndarray`_: Boolean array with the flags selecting the
            exposures in ``fitstbl`` that are ``ftype`` type frames.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        # ATTENTION: Standards have to be added manually for LUCI because
        # there is not unique flag that allows to distinguish between targets
        # and standards
        if ftype in ['science']:
            return good_exp & (fitstbl['idname'] == 'object')
        if ftype in ['standard']:
            return good_exp & (fitstbl['idname'] == 'standard')
        if ftype == 'bias':
            # for NIR data we type off lamp flats as biases
            return good_exp & (fitstbl['idname'] == 'flat_off')
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & (fitstbl['idname'] == 'flat')
        if ftype in ['dark']:
            # NOT Typing dark frames
            # return np.zeros(len(fitstbl), dtype=bool)
            # for testing dark typing uncommen the following line and comment
            # out the previous line
            return good_exp & (fitstbl['idname'] == 'dark')
        if ftype in ['arc', 'tilt']:
            return (good_exp & ((fitstbl['idname'] == 'object') |
                    (fitstbl['idname'] == 'arc')))

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def tweak_standard(self, wave_in, counts_in, counts_ivar_in, gpm_in,
                       meta_table, debug=False):
        """

        This routine is for performing instrument/disperser specific tweaks
        to standard stars so that sensitivity function fits will be well
        behaved.

        With regard to LUCI 1/2 the very edges of the spectra for both zJspec
        and HKspec configurations will be masked removing the extremely steep
        drop in sensitivity, predominantly at the blue edge of the spectrum.

        This function is copied and adapted from keck_mosfire.py

        Parameters
        ----------
        wave_in: (float np.ndarray) shape = (nspec,)
            Input standard star wavelenghts
        counts_in: (float np.ndarray) shape = (nspec,)
            Input standard star counts
        counts_ivar_in: (float np.ndarray) shape = (nspec,)
            Input inverse variance of standard star counts
        gpm_in: (bool np.ndarray) shape = (nspec,)
            Input good pixel mask for standard
        meta_table: (astropy.table)
            Table containing meta data that is slupred from the specobjs object. See unpack_object routine in specobjs.py
            for the contents of this table.

        Returns
        -------
        wave_out: (float np.ndarray) shape = (nspec,)
            Output standard star wavelenghts
        counts_out: (float np.ndarray) shape = (nspec,)
            Output standard star counts
        counts_ivar_out: (float np.ndarray) shape = (nspec,)
            Output inverse variance of standard star counts
        gpm_out: (bool np.ndarray) shape = (nspec,)
            Output good pixel mask for standard

        """

        wvmin = np.min(wave_in)
        wvmax = np.max(wave_in)

        if wvmin < 10000 and wvmax > 13000:
            coverage = 'zJspec'

        elif wvmin < 15500 and wvmax > 22000:
            coverage = 'HKspec'

        else:
            coverage = None
            msgs.warn('Standard tweaks not available for provided wavelength '
                      'range')
            msgs.warn('Continuing without tweaking standard spectrum.')

        if 'G200 LoRes' in meta_table['DISPNAME'] and coverage is not None:

            # Setting the wavelength edges depending on the wavelength
            # and the spectrograph LUCI1/LUCI2

            if 'lbt_luci2' in meta_table['PYP_SPEC'] and coverage == 'zJspec':

                wave_blue = 9800  # blue wavelength below which there is a strong
                # drop off in flux
                wave_red = 13400  # red wavelength above which there is a strong
                # drop off in flux combined with telluric absorption

            elif 'lbt_luci2' in meta_table['PYP_SPEC'] and coverage == 'HKspec':

                wave_blue = 13000  # blue wavelength below which there is a
                # strong
                # drop off in flux
                wave_red = 25000 # red wavelength above which there is a strong
                # drop off in flux combined with telluric absorption

            if 'lbt_luci1' in meta_table['PYP_SPEC'] and coverage == 'zJspec':

                wave_blue = 9400  # blue wavelength below which there is a
                # strong
                # drop off in flux
                wave_red = 13400  # red wavelength above which there is a strong
                # drop off in flux combined with telluric absorption

            elif 'lbt_luci1' in meta_table['PYP_SPEC'] and coverage == 'HKspec':

                wave_blue = 15200  # blue wavelength below which there is a
                # strong drop off in flux
                wave_red = 23100 # red wavelength above which there is a strong
                # drop off in flux combined with telluric absorption


            msgs.warn('Tweaking standard spectrum coverage.')
            msgs.warn('Standard spectrum reduced wavelenghts: {} to {'
                         '}'.format(wave_blue, wave_red))
            msgs.warn('For the telluric correction (pypeit_tellfilt), '
                      'make sure to mask the excluded regions.')




            exclusion_region = (wave_in < wave_blue) | (wave_in > wave_red)

            wave = wave_in.copy()
            counts = counts_in.copy()
            gpm = gpm_in.copy()
            counts_ivar = counts_ivar_in.copy()
            # By setting the wavelengths to zero, we guarantee that the sensitvity function will only be computed
            # over the valid wavelength region. While we could mask, this would still produce a wave_min and wave_max
            # for the zeropoint that includes the bad regions, and the polynomial fits will extrapolate crazily there
            wave[exclusion_region] = 0.0
            counts[exclusion_region] = 0.0
            counts_ivar[exclusion_region] = 0.0
            gpm[exclusion_region] = False
            if debug:
               from matplotlib import pyplot as plt
               counts_sigma = np.sqrt(utils.inverse(counts_ivar_in))
               plt.plot(wave_in, counts, color='red', alpha=0.7, label='apodized flux')
               plt.plot(wave_in, counts_in, color='black', alpha=0.7, label='flux')
               plt.plot(wave_in, counts_sigma, color='blue', alpha=0.7, label='flux')
               plt.axvline(wave_blue, color='blue')
               plt.axvline(wave_red, color='red')
               plt.legend()
               plt.show()
            return wave, counts, counts_ivar, gpm
        else:
            return wave_in, counts_in, counts_ivar_in, gpm_in

# Detector information from official LBT LUCI website
# https://sites.google.com/a/lbto.org/luci/instrument-characteristics/detector
 # Parameter 	 LUCI1 	 LUCI2
 # Gain (e-/ADU) 	 2.21 	 2.15
 # Dark Current (e-/sec/pix) 	  0.03 	 0.006
 # Read Noise (e-) 	 LIR 	  9.6 	9.2
 # MER 	  5.1 	 4.5
 # Full Well Capacity (e-) 	  TBD 	 122,000
 # Linearity (% at 40k ADU) 	  11.6% 	 11.1%
 # Crosstalk at Saturation 	  TBD 	 ≤0.2%
 # Persistence (% after 5min) 	  TBD 	 0.03%
 # Minimum DIT (sec) 	 LIR 	  2.503 	 2.503
 # MER 	  6.31 	 6.31



class LBTLUCI1Spectrograph(LBTLUCISpectrograph):
    """
    Child to handle LBT/LUCI1 specific code
    """
    name = 'lbt_luci1'
    camera = 'LUCI1'
    header_name = 'LUCI1'
    supported = True

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """

        # Get detector parameters from hdu header if available

        if hdu is not None:
            readmode = self.get_meta_value(self.get_headarr(hdu), 'readmode')

            if readmode == 'LIR':
                ronoise = np.atleast_1d(9.6)
            elif readmode == 'MER':
                ronoise = np.atleast_1d(5.1)
            else:
                msgs.error("Read mode not recognized (options: LIR, MER)")
                raise ValueError()

        # Detector 1
        detector_dict = dict(
            binning='1,1',
            det=1,
            dataext=0,
            specaxis=1,
            specflip=False,
            spatflip=False,
            platescale=0.2496,
            # Dark current nominally is < 360 electrons per hours
            # but the dark subtraction will effectively bring this to 0
            darkcurr=0.0,
            # Saturation is 55000, but will be set to dummy value for
            # now as integrated exposures over multiple detector integrations
            # will provide higher counts.
            saturation=1e+8,
            # NIR detectors are non-linear even in lower percentages
            # of the full well, thus for precision measurements one
            # should take into account a general non-linearity
            # correction.
            nonlinear=0.80,
            mincounts=-1e10,
            # In fact there are 32 amplifiers, which gain and ronoise
            # are extremely similar to each other, thus it will be
            # mimicked as 1
            numamplifiers=1,
            gain=np.atleast_1d(2.21),
            # The readout noise for LUCI are different for
            # different readout modes. We will be adopting the read out noise
            # for the LIR mode as it is higher the MER readout noise will be
            # commented out.
            ronoise=ronoise,  # variable populated from readmode meta above
            datasec=np.atleast_1d('[5:2044,5:2044]'),
            # For Luci the first 4 pixels on each side can
            # technically be used for as a biassec. This is not
            # included here.
            oscansec=np.atleast_1d('[5:2044,1:4]'),
            )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths'][
            'rms_threshold'] = 0.20  # 0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['fwhm'] = 4.0
        par['calibrations']['wavelengths']['n_final'] = 4
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'

        # Reduced from 0.93 to 0.9
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.9

        # Identification of slit edges
        par['calibrations']['slitedges']['edge_thresh'] = 25.0
        par['calibrations']['slitedges']['minimum_slit_length'] = 30.
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Extraction
        # Model full slit currently turned on
        par['reduce']['extraction']['model_full_slit'] = True
        # Tailored profile nsigma parameter for the standard, trying 100 (30
        # was standard
        par['reduce']['extraction']['std_prof_nsigma'] = 100.
        # Do not perform global sky subtraction for standard stars
        par['reduce']['skysub']['global_sky_std'] = False
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0

        # Processing steps
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Flexure
        par['flexure']['spec_method'] = 'skip'

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] = 'nothing'
        # par['scienceframe']['process']['satpix'] = 'reject'

        return par


    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the ``PypeIt`` parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            scifile (:obj:`str`):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """

        par = super().config_specific_par(scifile, inp_par=inp_par)

        # Set the wavelength identification
        dispname = self.get_meta_value(scifile, 'dispname')
        filter = self.get_meta_value(scifile, 'filter1')
        cenwave = self.get_meta_value(scifile, 'dispangle')

        if dispname == 'G200 LoRes' and filter == 'zJspec' and cenwave > \
                1.165 and cenwave < 1.175:
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = \
                'lbt_luci1_g200_zJ.fits'

        elif dispname == 'G200 LoRes' and filter == 'HKspec' and cenwave > \
                1.925 and cenwave < 1.935:
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = \
                'lbt_luci1_g200_HK.fits'
            # par['calibrations']['wavelengths']['ech_fix_format'] = False

        return par

class LBTLUCI2Spectrograph(LBTLUCISpectrograph):
    """
    Child to handle LBT/LUCI2 specific code
    """
    name = 'lbt_luci2'
    camera = 'LUCI2'
    header_name = 'LUCI2'
    supported = True

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """

        if hdu is not None:
            readmode = self.get_meta_value(self.get_headarr(hdu), 'readmode')

            if readmode == 'LIR':
                ronoise = np.atleast_1d(9.2)
            elif readmode == 'MER':
                ronoise = np.atleast_1d(4.5)
            else:
                msgs.error("Read mode not recognized (options: LIR, MER)")
                raise ValueError()


        # Detector 1
        detector_dict = dict(
            binning='1,1',
            det=1,
            dataext=0,
            specaxis=1,
            specflip=False,
            spatflip=False,
            platescale=0.2496,
            darkcurr=0.0,
            # Saturation is 55000, but will be set to dummy value for
            # now as integrated exposures over multiple detector integrations
            # will provide higher counts.
            saturation=1e+8,
            nonlinear=0.80,
            mincounts=-1e10,
            numamplifiers=1,
            gain=np.atleast_1d(2.15),
            ronoise=ronoise,  # variable populated from readmode meta above
            datasec= np.atleast_1d('[5:2044,5:2044]'),
            oscansec= np.atleast_1d('[5:2044,1:4]'),
            )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths'][
            'rms_threshold'] = 0.20  # 0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['fwhm'] = 4.0
        par['calibrations']['wavelengths']['n_final'] = 4
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'

        # Identification of slit edges
        par['calibrations']['slitedges']['edge_thresh'] = 25.0
        par['calibrations']['slitedges']['minimum_slit_length'] = 30.
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Reduced from 0.93 to 0.9
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.9

        # Extraction
        # Model full slit currently turned on
        par['reduce']['extraction']['model_full_slit'] = True
        # Tailored profile nsigma parameter for the standard
        par['reduce']['extraction']['std_prof_nsigma'] = 100.
        # Do not perform global sky subtraction for standard stars
        par['reduce']['skysub']['global_sky_std'] = False
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0

        # Processing steps
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Flexure
        par['flexure']['spec_method'] = 'skip'

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] = 'nothing'
        # par['scienceframe']['process']['satpix'] = 'reject'

        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the ``PypeIt`` parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            scifile (:obj:`str`):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """

        par = super().config_specific_par(scifile, inp_par=inp_par)

        # Set the wavelength identification
        dispname = self.get_meta_value(scifile, 'dispname')
        filter = self.get_meta_value(scifile, 'filter1')
        cenwave = self.get_meta_value(scifile, 'dispangle')

        if dispname == 'G200 LoRes' and filter == 'zJspec':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = \
                'lbt_luci2_g200_zJ.fits'

        elif dispname == 'G200 LoRes' and filter == 'HKspec':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = \
                'lbt_luci2_g200_HK.fits'
            # par['calibrations']['wavelengths']['ech_fix_format'] = False

        return par