"""
Module for P200/DBSP specific methods.

.. include:: ../include/links.rst
"""
from typing import List, Optional

import numpy as np

from astropy.io import fits
from astropy.coordinates import Angle
from astropy import units as u
from astropy.time import Time

from pypeit import msgs
from pypeit import io
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container


def flip_fits_slice(s: str) -> str:
    return '[' + ','.join(s.strip('[]').split(',')[::-1]) + ']'


class P200DBSPSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle P200/DBSP specific code
    """
    ndet = 1
    telescope = telescopes.P200TelescopePar()
    url = 'https://sites.astro.caltech.edu/palomar/observer/200inchResources/dbspoverview.html'

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA', required_ftypes=['science', 'standard'])
        self.meta['dec'] = dict(ext=0, card='DEC', required_ftypes=['science', 'standard'])
        self.meta['target'] = dict(ext=0, card='OBJECT')

        self.meta['dispname'] = dict(ext=0, card='GRATING')
        self.meta['decker'] = dict(ext=0, card='APERTURE')
        self.meta['binning'] = dict(card=None, compound=True)

        self.meta['mjd'] = dict(card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS', required_ftypes=['science', 'standard'])

        # Extras for config and frametyping
        self.meta['dichroic'] = dict(ext=0, card='DICHROIC')
        self.meta['dispangle'] = dict(card=None, rtol=1e-2, compound=True)
        self.meta['slitwid'] = dict(ext=0, card='APERTURE')
        self.meta['idname'] = dict(ext=0, card='IMGTYPE')
        self.meta['instrument'] = dict(ext=0, card='FPA')
        # Lamps
        self.meta['lampstat01'] = dict(ext=0, card='LAMPS')

    def compound_meta(self, headarr: List[fits.Header], meta_key: str):
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
        if meta_key == 'mjd':
            return Time(headarr[0]['UTSHUT']).mjd
        elif meta_key == 'dispangle':
            try:
                return Angle(headarr[0]['ANGLE'].lower()).deg
            except Exception as e:
                msgs.warn("Could not read dispangle from header:" + msgs.newline() + str(headarr[0]['ANGLE']))
                raise e
        else:
            return None

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
        return ['dispname', 'binning', 'dispangle', 'dichroic']

    def raw_header_cards(self):
        """
        Return additional raw header cards to be propagated in
        downstream output files for configuration identification.

        The list of raw data FITS keywords should be those used to populate
        the :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.configuration_keys`
        or are used in :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.config_specific_par`
        for a particular spectrograph, if different from the name of the
        PypeIt metadata keyword.

        This list is used by :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.subheader_for_spec`
        to include additional FITS keywords in downstream output files.

        Returns:
            :obj:`list`: List of keywords from the raw data files that should
            be propagated in output files.
        """
        return ['GRATING', 'ANGLE', 'APERTURE']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys() + ['slitwid']
    
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
        if ftype in ['science', 'standard']:
            return good_exp & (fitstbl['lampstat01'] == '0000000') & (fitstbl['idname'] == 'object')
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'bias')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (fitstbl['idname'] == 'flat')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['lampstat01'] != '0000000') & (fitstbl['idname'] == 'cal')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


class P200DBSPBlueSpectrograph(P200DBSPSpectrograph):
    """
    Child to handle P200/DBSP blue specific code
    """
    name = 'p200_dbsp_blue'
    camera = 'DBSPb'
    header_name = 'DBSP_BLUE'
    supported = True
    comment = 'Blue camera'
    
    def compound_meta(self, headarr: List[fits.Header], meta_key: str):
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
        # Handle dispangle and mjd from superclass method
        retval = super().compound_meta(headarr, meta_key)

        # If superclass could not handle the meta key
        if retval is not None:
            return retval
        if meta_key == 'binning':
            binspatial, binspec = headarr[0]['CCDSUM'].split(' ')
            return parse.binning2string(binspec, binspatial)
        msgs.error("Not ready for this compound meta")

    def get_detector_par(self, det: int, hdu: Optional[fits.HDUList] = None):
        """
        Return metadata for the selected detector.

        .. warning::

            Many of the necessary detector parameters are read from the file
            header, meaning the ``hdu`` argument is effectively **required** for
            P200/DBSPb.  The optional use of ``hdu`` is only viable for
            automatically generated documentation.

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
        if hdu is None:
            binning = '1,1'
            datasec = None
            oscansec = None
        else:
            # TODO: Could this be detector dependent??
            binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
            datasec = np.atleast_1d(flip_fits_slice(hdu[0].header['TSEC1']))
            oscansec = np.atleast_1d(flip_fits_slice(hdu[0].header['BSEC1']))

        # Detector 1
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 0,
            specflip        = True,
            spatflip        = False, # check
            platescale      = 0.389,
            darkcurr        = 0.0,
            saturation      = 65000.,
            nonlinear       = 62./65.,
            mincounts       = -1e10, # cross-check
            numamplifiers   = 1,
            gain            = np.atleast_1d(0.72),
            ronoise         = np.atleast_1d(2.5),
            datasec         = datasec,
            oscansec        = oscansec
            )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        par['scienceframe']['process']['combine'] = 'median'
        par['calibrations']['standardframe']['process']['combine'] = 'median'

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.55


        par['scienceframe']['process']['use_overscan'] = True
        # Make a bad pixel mask
        par['calibrations']['bpm_usebias'] = True
        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        # Change the wavelength calibration method
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['lamps'] = ['FeI', 'ArI', 'ArII']

        #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        #par['calibrations']['wavelengths']['n_first'] = 3
        #par['calibrations']['wavelengths']['n_final'] = 5
        #par['calibrations']['wavelengths']['sigdetect'] = 10.0
        # Do not flux calibrate
        par['fluxcalib'] = None
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, 120]
        par['calibrations']['standardframe']['exprng'] = [None, 120]
        par['scienceframe']['exprng'] = [90, None]

        par['sensfunc']['UVIS']['nresln'] = 5

        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
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
        # Start with instrument wide
        par = super().config_specific_par(scifile, inp_par=inp_par)

        grating = self.get_meta_value(scifile, 'dispname')
        dichroic = self.get_meta_value(scifile, 'dichroic')

        angle = Angle(self.get_meta_value(scifile, 'dispangle'), unit=u.deg).rad
        slitwidth = self.get_meta_value(scifile, 'slitwid') * u.arcsec
        lines_mm = float(grating.split('/')[0]) / u.mm

        theta_m = 38.5 * 2*np.pi / 360. - angle
        order = 2. if lines_mm == 158. * u.mm else 1.
        platescale = 0.389 * u.arcsec / u.pix
        pix_size = 15 * u.um

        disp = np.cos(theta_m)/(lines_mm * 9 * u.imperial.inch) * 1e7 * u.AA / u.mm
        cen_wv = np.abs(1/lines_mm * (np.sin(theta_m) - np.sin(angle)) / order)
        dlam = slitwidth / platescale * pix_size / u.pix * disp

        resolving_power = cen_wv / dlam

        par['sensfunc']['UVIS']['resolution'] = resolving_power.decompose().value

        cen_wv_AA = cen_wv.to(u.AA).value

        high_res_reids = {
            '1200/5000': {
                'D55': {
                    4700: 'p200_dbsp_blue_1200_5000_d55_4700.fits'
                },
                'D68': {
                    6000: 'p200_dbsp_blue_1200_5000_d68_6000.fits'
                }
            }
        }

        if grating.split("/")[0] == "1200":
            # high resolution grating!
            # here we need to select the reid_arxiv that most closely matches the central wavelength
            # and emit a warning if the difference is too great / the wavelength overlap is too small
            try:
                reids = high_res_reids[grating][dichroic]
                cen_wvs = np.array(list(reids))
                best_wv = cen_wvs[np.argmin(np.abs(cen_wvs - cen_wv_AA))]

                # blue wavelength coverage with a 1200 lines/mm grating is about 1550 A
                diff = np.abs(best_wv - cen_wv_AA)
                if diff > 775:
                    msgs.warn("Closest matching archived wavelength solutions"
                        f"differs in central wavelength by {diff:4.0f} A. The"
                        "wavelength solution may be unreliable. If wavelength"
                        "calibration fails, try using the holy grail method by"
                        "adding the following to your PypeIt file:\n"
                        "[calibrations]\n"
                        "\t[[wavelengths]]\n"
                        "\t\tmethod = holy-grail")
                par['calibrations']['wavelengths']['reid_arxiv'] = reids[best_wv]
            except KeyError:
                msgs.warn("Your grating " + grating + " doesn't have a template spectrum for the blue arm of DBSP.")
        else:
            if grating == '600/4000' and dichroic == 'D55':
                par['calibrations']['wavelengths']['reid_arxiv'] = 'p200_dbsp_blue_600_4000_d55.fits'
            elif grating == '600/4000' and dichroic == 'D68':
                par['calibrations']['wavelengths']['reid_arxiv'] = 'p200_dbsp_blue_600_4000_d68.fits'
            elif grating == '300/3990' and dichroic == 'D55':
                par['calibrations']['wavelengths']['reid_arxiv'] = 'p200_dbsp_blue_300_3990_d55.fits'
            else:
                msgs.warn("Your grating " + grating + " doesn't have a template spectrum for the blue arm of DBSP.")
        
        return par


class P200DBSPRedSpectrograph(P200DBSPSpectrograph):
    """
    Child to handle P200/DBSPr red specific code
    """
    name = 'p200_dbsp_red'
    camera = 'DBSPr'
    header_name = 'DBSP_RED2'
    supported = True
    comment = 'Red camera'
    
    def compound_meta(self, headarr: List[fits.Header], meta_key: str):
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
        # Handle dispangle and mjd from superclass method
        retval = super().compound_meta(headarr, meta_key)
        
        # If superclass could not handle the meta key
        if retval is not None:
            return retval
        if meta_key == 'binning':
            binspec, binspatial = headarr[0]['CCDSUM'].split(' ')
            return parse.binning2string(binspec, binspatial)
        else:
            msgs.error("Not ready for this compound meta")

    def get_detector_par(self, det: int, hdu: Optional[fits.HDUList] = None):
        """
        Return metadata for the selected detector.

        .. warning::

            Many of the necessary detector parameters are read from the file
            header, meaning the ``hdu`` argument is effectively **required** for
            P200/DBSPr.  The optional use of ``hdu`` is only viable for
            automatically generated documentation.

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
        if hdu is None:
            binning = '1,1'
            datasec = None
            oscansec = None
        else:
            # TODO: Could this be detector dependent??
            binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
            datasec = np.atleast_1d(flip_fits_slice(hdu[0].header['TSEC1']))
            oscansec = np.atleast_1d(flip_fits_slice(hdu[0].header['BSEC1']))

        # Detector 1
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = False,
            spatflip        = False, # check
            platescale      = 0.293,
            darkcurr        = 0.0,
            saturation      = 45000.,
            nonlinear       = 40./45.,
            mincounts       = -1e10, # check
            numamplifiers   = 1,
            gain            = np.atleast_1d(2.8),
            ronoise         = np.atleast_1d(8.5),
            datasec         = datasec,
            oscansec        = oscansec
        )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        par['scienceframe']['process']['combine'] = 'median'
        par['calibrations']['standardframe']['process']['combine'] = 'median'

        par['scienceframe']['process']['use_overscan'] = True
        par['scienceframe']['process']['sigclip'] = 4.0 # Tweaked downward from 4.5. 
        par['scienceframe']['process']['objlim'] = 1.5 # Tweaked downward from 3.0. Same value as Keck KCWI and DEIMOS
        # Make a bad pixel mask
        par['calibrations']['bpm_usebias'] = True
        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        # Change the wavelength calibration method
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['lamps'] = ['ArI', 'ArII', 'NeI', 'HeI']
        # par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        #par['calibrations']['wavelengths']['n_first'] = 3
        #par['calibrations']['wavelengths']['n_final'] = 5
        #par['calibrations']['wavelengths']['sigdetect'] = 10.0
        # Do not flux calibrate
        par['fluxcalib'] = None
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, 120]
        par['calibrations']['standardframe']['exprng'] = [None, 120]
        par['scienceframe']['exprng'] = [90, None]

        par['sensfunc']['algorithm'] = 'UVIS'
        par['sensfunc']['UVIS']['polycorrect'] = False
        par['sensfunc']['IR']['telgridfile'] = 'TelFit_Lick_3100_11100_R10000.fits'
        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
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
        # Start with instrument wide
        par = super().config_specific_par(scifile, inp_par=inp_par)

        grating = self.get_meta_value(scifile, 'dispname')
        dichroic = self.get_meta_value(scifile, 'dichroic')

        angle = Angle(self.get_meta_value(scifile, 'dispangle'), unit=u.deg).rad
        slitwidth = self.get_meta_value(scifile, 'slitwid') * u.arcsec
        lines_mm = float(grating.split('/')[0]) / u.mm

        theta_m = 35.0 * 2*np.pi / 360. - angle
        order = 1.
        platescale = 0.293 * u.arcsec / u.pix
        pix_size = 15 * u.um

        disp = np.cos(theta_m)/(lines_mm * 12 * u.imperial.inch) * 1e7 * u.AA / u.mm
        cen_wv = np.abs(1/lines_mm * (np.sin(theta_m) - np.sin(angle)) / order)
        dlam = slitwidth / platescale * pix_size / u.pix * disp

        resolving_power = cen_wv / dlam

        par['sensfunc']['UVIS']['resolution'] = resolving_power.decompose().value

        cen_wv_AA = cen_wv.to(u.AA).value

        high_res_reids = {
            '1200/7100': {
                'D68': {
                    7600: 'p200_dbsp_red_1200_7100_d68.fits',
                    8200: 'p200_dbsp_red_1200_7100_d68.fits'
                }
            },
            '1200/9400': {
                'D55': {
                    8800: 'p200_dbsp_red_1200_9400_d55_8800.fits'
                }
            }
        }

        if grating.split("/")[0] == "1200":
            # high resolution grating!
            # here we need to select the reid_arxiv that most closely matches the central wavelength
            # and emit a warning if the difference is too great / the wavelength overlap is too small
            try:
                reids = high_res_reids[grating][dichroic]
                cen_wvs = np.array(list(reids))
                best_wv = cen_wvs[np.argmin(np.abs(cen_wvs - cen_wv_AA))]

                # red wavelength coverage with a 1200 lines/mm grating is about 1600 A
                diff = np.abs(best_wv - cen_wv_AA)
                if diff > 800:
                    msgs.warn("Closest matching archived wavelength solutions"
                        f"differs in central wavelength by {diff:4.0f} A. The"
                        "wavelength solution may be unreliable. If wavelength"
                        "calibration fails, try using the holy grail method by"
                        "adding the following to your PypeIt file:\n"
                        "[calibrations]\n"
                        "\t[[wavelengths]]\n"
                        "\t\tmethod = holy-grail")
                par['calibrations']['wavelengths']['reid_arxiv'] = reids[best_wv]
            except KeyError:
                msgs.warn("Your grating " + grating + " doesn't have a template spectrum for the red arm of DBSP.")
        else:
            if grating == '316/7500' and dichroic == 'D55':
                par['calibrations']['wavelengths']['reid_arxiv'] = 'p200_dbsp_red_316_7500_d55.fits'
            elif grating == '600/10000' and dichroic == 'D55':
                par['calibrations']['wavelengths']['reid_arxiv'] = 'p200_dbsp_red_600_10000_d55.fits'
            else:
                msgs.warn("Your grating " + grating + " doesn't have a template spectrum for the red arm of DBSP.")

        return par

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Override parent bpm function with BPM specific to P200 DBSPr.

        Parameters
        ----------
        det : int
            Detector number
        msbias : numpy.ndarray
            Processed bias frame used when constructing the bpm (see :func:`bpm_frombias`)

        Returns
        -------
        bpix : ndarray
          0 = ok; 1 = Mask
        """
        msgs.info("Custom bad pixel mask for DBSPr")
        bpm_img = self.empty_bpm(filename, det, shape=shape)

        # Fill in bad pixels if a processed bias frame is provided
        if msbias is not None:
            return self.bpm_frombias(msbias, bpm_img)

        # Red CCD detector defect is present in data taken 2020-05-22
        # and absent in data taken 2020-04-21
        DEFECT_DATE = Time('2020-05-21')
        # TODO: Model the growth of the detector defect with time.
        # TODO: Get more precise date range for detector.
        with io.fits_open(filename) as hdul:
            if Time(hdul[0].header['UTSHUT']) > DEFECT_DATE:
                spec_binning = int(self.get_meta_value([hdul[0].header], 'binning').split(',')[0])
                bpm_img[464 // spec_binning : 723 // spec_binning, :] = 1

        return bpm_img


