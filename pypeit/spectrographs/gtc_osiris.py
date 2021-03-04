"""
Module for GTC OSIRIS specific methods.

.. include:: ../include/links.rst
"""
import glob
from pkg_resources import resource_filename

import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container

def flip_fits_slice(s: str) -> str:
    return '[' + ','.join(s.strip('[]').split(',')[::-1]) + ']'


class GTCOSIRISSpectrograph(spectrograph.Spectrograph):
    """
    Child of Spectrograph to handle GTC/OSIRIS specific code
    """
    ndet = 2
    name = 'gtc_osiris'
    telescope = telescopes.GTCTelescopePar()
    camera = 'OSIRIS'
    supported = False
#    comment = 'Grisms R2500V, R2500R'

    def get_detector_par(self, hdu, det):
        """
        Return metadata for the selected detector.

        Detector data from `here
        <http://www.gtc.iac.es/instruments/osiris/>`__.

        Args:
            hdu (`astropy.io.fits.HDUList`_):
                The open fits file with the raw image of interest.
            det (:obj:`int`):
                1-indexed detector number.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """

        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
        datasec1 = hdu[1].header['datasec']
        datasec2 = hdu[2].header['datasec']


        # Detector 1
        detector_dict1 = dict(
            binning         = binning,
            det             = 1,
            dataext         = 1,  # Not sure this is used
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.127,  # arcsec per pixel
            darkcurr        = 0.0,
            saturation      = 65535., # ADU
            nonlinear       = 0.95,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d([0.95]),
            ronoise         = np.atleast_1d([4.5]),
            datasec         = np.atleast_1d('[1:4102,52:1880]')
#            datasec         = np.atleast_1d(flip_fits_slice(datasec1)),
            )
        # Detector 2
        detector_dict2 = dict(
            binning         = binning,
            det             = 2,
            dataext         = 2,  # Not sure this is used
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.127,
            darkcurr        = 0.0,
            saturation      = 65535., # ADU
            nonlinear       = 0.95,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d([0.95]),
            ronoise         = np.atleast_1d([4.5]),
            datasec         = np.atleast_1d('[1:4102,52:1880]')
#            datasec         = np.atleast_1d(flip_fits_slice(datasec2)),
            )

        detectors = [detector_dict1, detector_dict2]
        # Return
        return detector_container.DetectorContainer(**detectors[det-1])


    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        # Wavelength calibration methods
        par['calibrations']['wavelengths']['method'] = 'holy-grail'
        #par['calibrations']['wavelengths']['method'] = 'full_template'
        # par['calibrations']['wavelengths']['method'] = 'identify'
        par['calibrations']['wavelengths']['lamps'] = ['GTC_OSIRIS']
        # par['calibrations']['wavelengths']['rms_threshold'] = 0.5
        # par['calibrations']['wavelengths']['sigdetect'] = 5.
        # par['calibrations']['wavelengths']['fwhm']= 5.0
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, None]  # Long arc exposures on this telescope
        par['calibrations']['arcframe']['process']['sigclip']=2.
        par['calibrations']['arcframe']['process']['objlim']=2.
        par['calibrations']['arcframe']['process']['cr_sigrej']=10.
        par['calibrations']['arcframe']['process']['replace']='mean'
        par['calibrations']['arcframe']['process']['clip']=False
        par['calibrations']['arcframe']['process']['combine']='weightmean'
        par['calibrations']['standardframe']['exprng'] = [None, 120]
        par['scienceframe']['exprng'] = [90, None]

        # No ovescan region!
        turn_off = dict(use_overscan=False)
        par.reset_all_processimages_par(**turn_off)
        par['scienceframe']['process']['use_overscan'] = False

       # Extraction
        # par['reduce']['skysub']['bspline_spacing'] = 0.8
        # par['reduce']['skysub']['no_poly'] = True
        # par['reduce']['skysub']['bspline_spacing'] = 0.6
        # par['reduce']['skysub']['joint_fit'] = True
        # par['reduce']['skysub']['global_sky_std']  = False
        #
        # par['reduce']['extraction']['sn_gauss'] = 4.0
        # par['reduce']['findobj']['sig_thresh'] = 3.0
        # par['reduce']['skysub']['sky_sigrej'] = 5.0
        # par['reduce']['findobj']['find_trim_edge'] = [5,5]

        # cosmic ray rejection parameters for science frames
        # par['scienceframe']['process']['sigclip'] = 3.0
        # par['scienceframe']['process']['objlim'] = 2.0



        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='object')
        self.meta['idname'] = dict(ext=0, card='obsmode')
        self.meta['decker'] = dict(ext=0, card='SLITW')
        self.meta['binning'] = dict(card=None, compound=True)  # Uses CCDSUM
        self.meta['detector']=dict(ext=0,card='detector')
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='GRISM')
#        self.meta['dichroic'] = dict(ext=0, card='FILTER1')
        self.meta['datasec'] = dict(ext=1, card='DATASEC')

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
        if meta_key == 'binning':
            binspatial, binspec = parse.parse_binning(headarr[0]['CCDSUM'])
            binning = parse.binning2string(binspec, binspatial)
            return binning

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
        return ['dispname', 'decker', 'binning']

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
        if ftype in ['science','standard']:
            return good_exp & (fitstbl['target'] != 'ArcLamp_Xe') \
            & (fitstbl['target'] != 'ArcLamp_HgAr') \
            & (fitstbl['target'] != 'ArcLamp_Ne') \
            & (fitstbl['target'] != 'SpectralFlat') \
                    & (fitstbl['target'] != 'BIAS')
        if ftype in ['arc', 'tilt']:
            return good_exp & ((fitstbl['target'] == 'ArcLamp_Xe') \
            | (fitstbl['target'] == 'ArcLamp_HgAr') \
            | (fitstbl['target'] == 'ArcLamp_Ne'))
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (fitstbl['target'] == 'SpectralFlat')
        if ftype == 'bias':
            return good_exp & (fitstbl['target'] == 'BIAS')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


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
        # Start with instrument wide
        par = super().config_specific_par(scifile, inp_par=inp_par)

        # Wavelength calibrations
        if self.get_meta_value(scifile, 'dispname') == 'R2500V':
            par['calibrations']['wavelengths']['wv_cen'] = 5185.
            par['calibrations']['wavelengths']['disp'] = 0.8

        elif self.get_meta_value(scifile, 'dispname') == 'R2500R':
            par['calibrations']['wavelengths']['wv_cen'] = 6560.
            par['calibrations']['wavelengths']['disp'] = 1.04
        else:
            msgs.warn('gtc_osiris.py: wv_cen and disp missing for this grism')

        # Return
        return par
