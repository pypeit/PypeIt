"""
Module for NOT ALFOSC spectrograph

.. include:: ../include/links.rst
"""
from IPython import embed

import numpy as np

from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container


class NOTALFOSCSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle NOT ALFOSC spectrograph
    """
    ndet = 1
    name = 'not_alfosc'
    telescope = telescopes.NOTTelescopePar()
    camera = 'ALFOSC'
    url = 'https://www.not.iac.es/instruments/alfosc/'
    header_name = 'ALFOSC_FASU'
    supported = True
    comment = 'For use with the standard horizontal slits only. Grisms 3, 4, 5, 7, 8, 10, 11, 17, 18, 19, 20'

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Detector data from `here
        <http://www.not.iac.es/instruments/detectors/CCD14/>`__.

        .. warning::

            Many of the necessary detector parameters are read from the file
            header, meaning the ``hdu`` argument is effectively **required** for
            NOT/ALFOSC.  The optional use of ``hdu`` is only viable for
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
        # http://www.not.iac.es/instruments/detectors/CCD14/

        if hdu is None:
            binning = '1,1'
            gain = None
            ronoise = None
        else:
            binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
            gain = np.atleast_1d(hdu[1].header['GAIN'])  # e-/ADU
            ronoise = np.atleast_1d(hdu[1].header['RDNOISE'])  # e-

        # Detector 1
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 1,
            specaxis        = 0,
            specflip        = True,
            spatflip        = False,
            xgap            = 0.,
            ygap            = 0.,
            ysize           = 1.,
            platescale      = 0.2138,
            mincounts       = -1e10,
            darkcurr        = 1.3,      # e-/pix/hr
            saturation      = 700000.,  # ADU
            nonlinear       = 0.86,
            datasec         = np.atleast_1d('[:,{}:{}]'.format(1, 2148)),  # Unbinned
            oscansec        = None,
            numamplifiers   = 1,
            gain            = gain,     # e-/ADU
            ronoise         = ronoise   # e-
        )

#        # Parse datasec, oscancsec from the header
#        head1 = hdu[1].header
#        detector_dict['gain'] = np.atleast_1d(head1['GAIN'])  # e-/ADU
#        detector_dict['ronoise'] = np.atleast_1d(head1['RDNOISE'])  # e-

        # Return
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
        par['calibrations']['slitedges']['bound_detector'] = True
        # Flats are sometimes quite ugly due to dust on the slit which leads to the erroneous detection of multiple slits. So set a higher edge_thresh and minimum_slit_gap.
        par['calibrations']['slitedges']['edge_thresh'] = 30
        par['calibrations']['slitedges']['minimum_slit_gap'] = 15

        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        # Wavelength calibration methods
        #par['calibrations']['wavelengths']['method'] = 'holy-grail'
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['lamps'] = ['HeI', 'NeI', 'ArI']
        par['calibrations']['wavelengths']['sigdetect'] = 10.0
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, None]  # Long arc exposures on this telescope
        par['calibrations']['standardframe']['exprng'] = [None, 120]
        par['scienceframe']['exprng'] = [10, None]

        # Multiple arcs with different lamps, so can't median combine nor clip, also need to remove continuum
        par['calibrations']['arcframe']['process']['clip'] = False
        par['calibrations']['arcframe']['process']['combine'] = 'mean'
        par['calibrations']['arcframe']['process']['subtract_continuum'] = True
        par['calibrations']['tiltframe']['process']['clip'] = False
        par['calibrations']['tiltframe']['process']['combine'] = 'mean'
        par['calibrations']['tiltframe']['process']['subtract_continuum'] = True


        # No overscan region!
        turn_off = dict(use_overscan=False)
        par.reset_all_processimages_par(**turn_off)

        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        #self.meta['ra'] = dict(ext=0, card='OBJRA')
        self.meta['ra'] = dict(card=None, compound=True)
        self.meta['dec'] = dict(ext=0, card='OBJDEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(ext=0, card='ALAPRTNM')
        self.meta['binning'] = dict(card=None, compound=True)

        self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='ALGRNM')
        self.meta['idname'] = dict(ext=0, card='IMAGETYP')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')

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
            # PypeIt frame
            binspatial = headarr[0]['DETXBIN']
            binspec = headarr[0]['DETYBIN']
            return parse.binning2string(binspec, binspatial)
        elif meta_key == 'mjd':
            time = headarr[0]['DATE-AVG']
            ttime = Time(time, format='isot')
            return ttime.mjd
        elif meta_key == 'ra':
            objra = headarr[0]['OBJRA'] # Given in hours, not deg
            return objra*15.
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
        return ['dispname', 'decker', 'binning']

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
        return ['ALGRNM', 'ALAPRTNM', 'DETXBIN', 'DETYBIN']

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
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'OBJECT')
        if ftype == 'standard':
            return good_exp & ((fitstbl['idname'] == 'STD') | (fitstbl['target'] == 'STD') | (fitstbl['target'] == 'STD,SLIT'))
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'BIAS')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (fitstbl['idname'] == 'FLAT,LAMP')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc','tilt']:
            return good_exp & (fitstbl['idname'] == 'WAVE,LAMP')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

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

        # Wavelength calibrations
        if self.get_meta_value(scifile, 'dispname') == 'Grism_#3':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'not_alfosc_grism3.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'Grism_#4':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'not_alfosc_grism4.fits'
            par['calibrations']['wavelengths']['lamps'] = ['HeI','NeI']
        elif self.get_meta_value(scifile, 'dispname') == 'Grism_#5':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'not_alfosc_grism5.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'Grism_#7':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'not_alfosc_grism7.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'Grism_#8':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'not_alfosc_grism8.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'Grism_#10':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'not_alfosc_grism10.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'Grism_#11':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'not_alfosc_grism11.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'Grism_#17':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'not_alfosc_grism17.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'Grism_#18':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'not_alfosc_grism18.fits'
            par['calibrations']['wavelengths']['lamps'] = ['HeI','NeI','ArI','ArII']
        elif self.get_meta_value(scifile, 'dispname') == 'Grism_#19':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'not_alfosc_grism19.fits'
            par['calibrations']['wavelengths']['lamps'] = ['HeI','NeI','ArI','ArII']
        elif self.get_meta_value(scifile, 'dispname') == 'Grism_#20':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'not_alfosc_grism20.fits'
        else:
            msgs.warn('not_alfosc.py: YOU NEED TO ADD IN THE WAVELENGTH SOLUTION FOR THIS GRISM')

        # Return
        return par




class NOTALFOSCSpectrographVert(NOTALFOSCSpectrograph):
    """
    Child to handle Vertical slits for NOT ALFOSC spectrograph
    """
    name = 'not_alfosc_vert'
    comment = 'Grisms 3, 4, 5, 7, 8, 10, 11, 17, 18, 19, 20. For vertical slits only'

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Detector data from `here
        <http://www.not.iac.es/instruments/detectors/CCD14/>`__.

        .. warning::

            Many of the necessary detector parameters are read from the file
            header, meaning the ``hdu`` argument is effectively **required** for
            NOT/ALFOSC.  The optional use of ``hdu`` is only viable for
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
        # http://www.not.iac.es/instruments/detectors/CCD14/

        if hdu is None:
            binning = '1,1'
            gain = None
            ronoise = None
        else:
            binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
            gain = np.atleast_1d(hdu[1].header['GAIN'])  # e-/ADU
            ronoise = np.atleast_1d(hdu[1].header['RDNOISE'])  # e-

        # Detector 1
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 1,
            specaxis        = 1, #Vertical slits have horizontal spectral dispersion
            specflip        = False,
            spatflip        = False,
            xgap            = 0.,
            ygap            = 0.,
            ysize           = 1.,
            platescale      = 0.2138,
            mincounts       = -1e10,
            darkcurr        = 1.3,      # e-/pix/hr
            saturation      = 700000.,  # ADU
            nonlinear       = 0.86,
            datasec         = np.atleast_1d('[{}:{},:]'.format(1, 2102)),  # Unbinned
            oscansec        = None,
            numamplifiers   = 1,
            gain            = gain,     # e-/ADU
            ronoise         = ronoise   # e-
        )

        # Return
        return detector_container.DetectorContainer(**detector_dict)
