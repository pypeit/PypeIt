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


class GTCOSIRISSpectrograph(spectrograph.Spectrograph):
    """
    Child of Spectrograph to handle GTC/OSIRIS specific code
    """
    ndet = 2
    name = 'gtc_osiris'
    telescope = telescopes.GTCTelescopePar()
    camera = 'OSIRIS'
    supported = True
    comment = 'See :doc:`gtc_osiris`'

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
            mincounts       = 0,
            numamplifiers   = 1,
            gain            = np.atleast_1d([0.95]),
            ronoise         = np.atleast_1d([4.5]),
            datasec         = np.atleast_1d('[1:4102,52:1880]')
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
            mincounts       = 0,
            numamplifiers   = 1,
            gain            = np.atleast_1d([0.95]),
            ronoise         = np.atleast_1d([4.5]),
            datasec         = np.atleast_1d('[1:4102,52:1880]')
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
        par['calibrations']['slitedges']['bound_detector'] = True

        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        # Wavelength calibration methods
#        par['calibrations']['wavelengths']['method'] = 'holy-grail'
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI,NeI,ArI']

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, None]  # Long arc exposures
        par['calibrations']['arcframe']['process']['clip']=False
        par['calibrations']['arcframe']['process']['combine']='weightmean' #Multiple arcs with different lamps, so can't median combine
        par['calibrations']['standardframe']['exprng'] = [None, 120]
        par['scienceframe']['exprng'] = [90, None]

        # No ovescan region
        turn_off = dict(use_overscan=False)
        par.reset_all_processimages_par(**turn_off)
        par['scienceframe']['process']['use_overscan'] = False





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
        self.meta['decker'] = dict(ext=0, card='MASKNAME')
        self.meta['binning'] = dict(card=None, compound=True)  # Uses CCDSUM
        self.meta['detector']=dict(ext=0,card='detector')
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        self.meta['dispname'] = dict(ext=0, card='GRISM')
        self.meta['datasec'] = dict(ext=1, card='DATASEC')
        self.meta['dichroic'] = dict(ext=0, card='FILTER1')

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

    def config_independent_frames(self):
        """
        Define frame types that are independent of the fully defined
        instrument configuration.

        Bias and dark frames are considered independent of a configuration.
        Standards are assigned to the correct configuration frame group by
        grism (i.e. ignoring that they are taken with a wider slit).
        See :func:`~pypeit.metadata.PypeItMetaData.set_configurations`.

        Returns:
            :obj:`dict`: Dictionary where the keys are the frame types that
            are configuration independent and the values are the metadata
            keywords that can be used to assign the frames to a configuration
            group.
        """
        return {'standard': 'dispname','bias': None, 'dark': None}


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

        if self.get_meta_value(scifile, 'idname') == 'OsirisMOS':
            par['reduce']['findobj']['find_trim_edge'] = [1,1]
            par['calibrations']['slitedges']['sync_predict'] = 'pca'
            par['calibrations']['slitedges']['det_buffer'] = 1

        # Wavelength calibrations
        if self.get_meta_value(scifile, 'dispname') == 'R300B':
            par['calibrations']['wavelengths']['wv_cen'] = 4405.
            par['calibrations']['wavelengths']['disp'] = 4.96
            par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R300B.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R300R':
            par['calibrations']['wavelengths']['wv_cen'] = 6635.
            par['calibrations']['wavelengths']['disp'] = 7.74
            par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R300R.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R500B':
            par['calibrations']['wavelengths']['wv_cen'] = 4745.
            par['calibrations']['wavelengths']['disp'] = 3.54
            par['calibrations']['wavelengths']['lamps'] = ['HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R500B.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R500R':
            par['calibrations']['wavelengths']['wv_cen'] = 7165.
            par['calibrations']['wavelengths']['disp'] = 4.88
            par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R500R.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R1000B':
            par['calibrations']['wavelengths']['wv_cen'] = 5455.
            par['calibrations']['wavelengths']['disp'] = 2.12
            par['calibrations']['wavelengths']['lamps'] = ['ArI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R1000B.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R1000R':
            par['calibrations']['wavelengths']['wv_cen'] = 7430.
            par['calibrations']['wavelengths']['disp'] = 2.62
            par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R1000R.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R2000B':
            par['calibrations']['wavelengths']['wv_cen'] = 4755.
            par['calibrations']['wavelengths']['disp'] = 0.86
            par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R2000B.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R2500U':
            par['calibrations']['wavelengths']['wv_cen'] = 3975.
            par['calibrations']['wavelengths']['disp'] = 0.62
            par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R2500U.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R2500V':
            par['calibrations']['wavelengths']['wv_cen'] = 5185.
            par['calibrations']['wavelengths']['disp'] = 0.85
            par['calibrations']['wavelengths']['lamps'] = ['HgI','NeI','XeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R2500V.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R2500R':
            par['calibrations']['wavelengths']['wv_cen'] = 6560.
            par['calibrations']['wavelengths']['disp'] = 1.04
            par['calibrations']['wavelengths']['lamps'] = ['ArI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R2500R.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R2500I':
            par['calibrations']['wavelengths']['wv_cen'] = 8650.
            par['calibrations']['wavelengths']['disp'] = 1.36
            par['calibrations']['wavelengths']['lamps'] = ['ArI,XeI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R2500I.fits'
        else:
            msgs.warn('gtc_osiris.py: template arc missing for this grism! Trying holy-grail...')
            par['calibrations']['wavelengths']['method'] = 'holy-grail'

        # Return
        return par
