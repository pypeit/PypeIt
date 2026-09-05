"""
Module for INT/IDS specific methods.

The Intermediate Dispersion Spectrograph (IDS) is a long-slit
spectrograph mounted at the Cassegrain focus of the 2.5m Isaac Newton
Telescope (INT) at the Roque de los Muchachos Observatory, La Palma.
It offers two 4k x 2k CCDs (the blue-sensitive EEV10 and the
red-sensitive, default RED+2) and a suite of 16 gratings.

.. include:: ../include/links.rst
"""
from pathlib import Path

import numpy as np

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table

from pypeit import log
from pypeit import PypeItError
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.core import parse
from pypeit.par import parset
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container


class INTIDSSpectrograph(spectrograph.Spectrograph):
    """
    Parent to handle INT/IDS specific code.

    This base class holds everything shared by the two IDS cameras;
    use the detector-specific children
    :class:`INTIDSEEV10Spectrograph` (``int_ids_eev10``) and
    :class:`INTIDSREDPLUS2Spectrograph` (``int_ids_redplus2``).
    """
    ndet = 1
    header_name = 'IDS'
    telescope = telescopes.INTTelescopePar()
    url = 'https://www.ing.iac.es/astronomy/instruments/ids/'
    # UltraDAS writes raw ING frames with a .fit extension; also accept
    # the standard .fits variants
    allowed_extensions = ['.fit', '.fits', '.fits.gz']

    # Detector-specific attributes set by each child class:
    # value of the DETECTOR header card
    detector_name = None
    # unbinned spatial plate scale [arcsec/pixel]
    detector_platescale = None
    # is the raw spectral axis reversed (wavelength decreasing)?
    detector_specflip = None
    # dark current [e-/pixel/hour]
    detector_darkcurr = None

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        .. warning::

            Several detector parameters (gain, read noise, and the data
            and overscan sections) are read from the file header,
            meaning the ``hdu`` argument is effectively **required** for
            INT/IDS.  The optional use of ``hdu`` is only viable for
            automatically generated documentation.

        Parameters
        ----------
        det : :obj:`int`
            1-indexed detector number.
        hdu : `astropy.io.fits.HDUList`_, optional
            The open fits file with the raw image of interest.  If not
            provided, frame-dependent parameters are set to a default.

        Returns
        -------
        :class:`~pypeit.images.detector_container.DetectorContainer`
            Object with the detector metadata.
        """
        if hdu is None:
            binning = '1,1'
            gain = None
            ronoise = None
            datasec = None
            oscansec = None
        else:
            binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
            # Gain [e-/ADU] and read noise [e-] are recorded in the
            # image extension header
            gain = np.atleast_1d(hdu[1].header['GAIN'])
            ronoise = np.atleast_1d(hdu[1].header['READNOIS'])
            # TRIMSEC/BIASSEC are given in the FITS (x,y) convention and
            # in the read-out (windowed and binned) pixel coordinates;
            # flip them to the numpy order expected by PypeIt
            datasec = np.atleast_1d(
                self.fits_to_pypeit_section(hdu[1].header['TRIMSEC']))
            oscansec = np.atleast_1d(
                self.fits_to_pypeit_section(hdu[1].header['BIASSEC']))
            # Sanity check that the data match this camera-specific
            # child class (also enforced by check_spectrograph)
            detector = hdu[0].header['DETECTOR'].strip()
            if detector != self.detector_name:
                log.warning(f'DETECTOR card is {detector}, but the '
                            f'{self.name} class expects '
                            f'{self.detector_name}.')

        # Single-amplifier readout (AMPNAME='LH' in the test data)
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 1,
            specaxis        = 0,        # Spectral axis is NAXIS2
            specflip        = self.detector_specflip,
            spatflip        = False,
            xgap            = 0.,
            ygap            = 0.,
            ysize           = 1.,
            platescale      = self.detector_platescale,
            darkcurr        = self.detector_darkcurr,
            saturation      = 65535.,
            # Both IDS detectors are linear to better than ~1% over
            # their full dynamic range per the ING linearity reports
            nonlinear       = 0.99,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = gain,
            ronoise         = ronoise,
            datasec         = datasec,
            oscansec        = oscansec
        )

        return detector_container.DetectorContainer(**detector_dict)

    def check_spectrograph(self, filename):
        """
        Check that the selected spectrograph is the correct one for the
        input data, i.e. that the ``DETECTOR`` header card matches the
        camera handled by this class.

        Parameters
        ----------
        filename : :obj:`str`
            File to use when determining if the input spectrograph is
            the correct one.  Raises :class:`~pypeit.PypeItError` if
            the file was taken with the other IDS camera.
        """
        detector = self.get_meta_value(filename, 'detector').strip()
        if detector != self.detector_name:
            # Suggest the sibling class for the other camera
            other = ('int_ids_redplus2' if self.name == 'int_ids_eev10'
                     else 'int_ids_eev10')
            raise PypeItError(f'{filename} was taken with the '
                              f'{detector} detector, not '
                              f'{self.detector_name}. You may want to '
                              f'use {other} instead of {self.name}.')

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        # RA/DEC are sexagesimal strings (hours and degrees) in the
        # header; convert them to decimal degrees
        self.meta['ra'] = dict(ext=0, card=None, compound=True)
        self.meta['dec'] = dict(ext=0, card=None, compound=True)
        self.meta['target'] = dict(ext=0, card='OBJECT')
        # The dekker slide position serves as the aperture identifier
        self.meta['decker'] = dict(ext=0, card='DEKPOS')
        self.meta['binning'] = dict(ext=0, card=None, compound=True)
        # MJD at the start of the observation is provided directly
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')

        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='GRATNAME')
        self.meta['dispangle'] = dict(ext=0, card='GRATANGL', rtol=1e-3)
        self.meta['cenwave'] = dict(ext=0, card='CENWAVE', rtol=1e-2)
        self.meta['idname'] = dict(ext=0, card='IMAGETYP')
        # Slit width projected onto the sky [arcsec]
        self.meta['slitwid'] = dict(ext=0, card='SLTWDSKY')
        self.meta['detector'] = dict(ext=0, card='DETECTOR')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')

        # Lamps: AGARCLMP lists the arc lamps in use (e.g.
        # 'CuAr CuAr CuNe CuN'); 'W' indicates the tungsten flat lamp
        self.meta['lampstat01'] = dict(ext=0, card=None, compound=True)

    def compound_meta(self, headarr, meta_key):
        """
        Methods to generate metadata requiring interpretation of the
        header data, instead of simply reading the value of a header
        card.

        Parameters
        ----------
        headarr : :obj:`list`
            List of `astropy.io.fits.Header`_ objects.
        meta_key : :obj:`str`
            Metadata keyword to construct.

        Returns
        -------
        object
            Metadata value read from the header(s).
        """
        if meta_key in ['ra', 'dec']:
            # Convert the sexagesimal strings (RA in hours, DEC in
            # degrees) to decimal degrees
            coord = SkyCoord(headarr[0]['RA'], headarr[0]['DEC'],
                             unit=('hourangle', 'deg'))
            return coord.ra.deg if meta_key == 'ra' else coord.dec.deg
        if meta_key == 'binning':
            # CCDSUM is space-separated with the spatial (x) binning
            # first and the spectral (y) binning second
            binspatial, binspec = headarr[0]['CCDSUM'].split()
            return parse.binning2string(binspec, binspatial)
        if meta_key == 'lampstat01':
            # AGARCLMP is blank when no comparison lamp is in use
            lamps = headarr[0].get('AGARCLMP', '').strip()
            return lamps if len(lamps) > 0 else 'off'

        raise PypeItError(f'Not ready for compound meta, {meta_key}, '
                          'for INT/IDS.')

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        This list is used by :class:`~pypeit.metadata.PypeItMetaData`
        to identify the unique configurations among the list of frames
        read for a given reduction.

        Returns
        -------
        :obj:`list`
            List of keywords of data pulled from file headers and used
            to construct the :class:`~pypeit.metadata.PypeItMetaData`
            object.
        """
        return ['dispname', 'cenwave', 'detector', 'binning']

    def raw_header_cards(self):
        """
        Return additional raw header cards to be propagated in
        downstream output files for configuration identification.

        Returns
        -------
        :obj:`list`
            List of keywords from the raw data files that should be
            propagated in output files.
        """
        return ['GRATNAME', 'CENWAVE', 'DETECTOR', 'CCDSUM']

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns
        -------
        :class:`~pypeit.par.pypeitpar.PypeItPar`
            Parameters required by all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Wavelength calibration: the test data use CuAr+CuNe lamps.
        # No archived template exists yet, so start with holy-grail.
        par['calibrations']['wavelengths']['lamps'] = ['NeI', 'ArI',
                                                       'ArII', 'CuI']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'
        # ~1.5 arcsec slit at ~0.4 arcsec/pixel gives a ~4 pixel FWHM
        par['calibrations']['wavelengths']['fwhm'] = 4.0
        par['calibrations']['wavelengths']['sigdetect'] = 5.

        # Processing: overscan and bias frames are available; no dark
        # frames in the typical calibration plan
        turn_off = dict(use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Single long slit: do not attempt left/right edge syncing from
        # multiple slits and bound the slit at the detector edge
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['bound_detector'] = True

        # Frame-typing exposure-time ranges.  Standards are typically
        # short exposures of bright stars; science exposures are longer.
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['arcframe']['exprng'] = [None, 600]
        par['calibrations']['pixelflatframe']['exprng'] = [None, 600]
        par['calibrations']['traceframe']['exprng'] = [None, 600]
        par['calibrations']['standardframe']['exprng'] = [None, 120]
        par['scienceframe']['exprng'] = [120, None]

        return par

    def config_specific_par(
            self,
            inp:str|list|Path|fits.Header|Table,
            inp_par:parset.ParSet|None=None
        ) -> parset.ParSet:
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Parameters
        ----------
        inp : :obj:`str`, :obj:`list`, `Path`_, `astropy.io.fits.Header`_, `astropy.table.Table`_
            Input filename, an `astropy.io.fits.Header`_ object, or a
            list of `astropy.io.fits.Header`_ objects.  Or a row from
            the metadata table.
        inp_par : :class:`~pypeit.par.parset.ParSet`, optional
            Parameter set used for the full run of PypeIt.  If None,
            use :func:`default_pypeit_par`.

        Returns
        -------
        :class:`~pypeit.par.parset.ParSet`
            The PypeIt parameter set adjusted for configuration
            specific parameter values.
        """
        # Start with the instrument-wide parameters
        par = super().config_specific_par(inp, inp_par=inp_par)

        # Select the archived wavelength template by grating.  Only
        # R1200B has a template so far (built from an IRAF solution of
        # a CuAr+CuNe arc; see the PypeIt-development-suite
        # pypeitdev/int_ids directory); other gratings fall back on
        # the holy-grail default.
        grating = self.get_meta_value(inp, 'dispname')
        if grating == 'R1200B':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] \
                    = 'int_ids_R1200B.fits'

        return par

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt
        file.

        Returns
        -------
        :obj:`list`
            The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print
            to the :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys() + ['slitwid', 'lampstat01']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.

        Frame typing relies primarily on the ``IMAGETYP`` header card
        (mapped to ``idname``): 'zero' (bias), 'flat', 'arc', 'dark',
        and 'object'.  Science and standard frames are both typed as
        'object' and are separated by their exposure times.

        Parameters
        ----------
        ftype : :obj:`str`
            Type of frame to check.  Must be a valid frame type; see
            frame-type :ref:`frame_type_defs`.
        fitstbl : `astropy.table.Table`_
            The table with the metadata for one or more frames to
            check.
        exprng : :obj:`list`, optional
            Range in the allowed exposure time for a frame of type
            ``ftype``.  See
            :func:`pypeit.core.framematch.check_frame_exptime`.

        Returns
        -------
        `numpy.ndarray`_
            Boolean array with the flags selecting the exposures in
            ``fitstbl`` that are ``ftype`` type frames.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'],
                                                  exprng)
        if ftype == 'bias':
            return fitstbl['idname'] == 'zero'
        if ftype == 'dark':
            return good_exp & (fitstbl['idname'] == 'dark')
        if ftype in ['science', 'standard']:
            # Both are 'object' frames with no comparison lamp;
            # exposure-time ranges (set in default_pypeit_par)
            # distinguish them
            return (good_exp
                    & (fitstbl['idname'] == 'object')
                    & (fitstbl['lampstat01'] == 'off'))
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'arc')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            # Tungsten ('W') lamp flats
            return (good_exp
                    & (fitstbl['idname'] == 'flat')
                    & (fitstbl['lampstat01'] == 'W'))

        log.debug(f'Cannot determine if frames are of type {ftype}.')
        return np.zeros(len(fitstbl), dtype=bool)

    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces that
        are key for image processing.

        This simply calls the base-class method, but with
        ``sec_includes_binning=True`` because the ``TRIMSEC`` and
        ``BIASSEC`` cards in the IDS headers are written in the
        read-out (i.e., windowed and binned) pixel coordinates.

        Parameters
        ----------
        raw_file : :obj:`str`
            File to read.
        det : :obj:`int`
            1-indexed detector to read.

        Returns
        -------
        detector_par : :class:`~pypeit.images.detector_container.DetectorContainer`
            Detector metadata parameters.
        raw_img : `numpy.ndarray`_
            Raw image for this detector.
        hdu : `astropy.io.fits.HDUList`_
            Opened fits file.
        exptime : :obj:`float`
            Exposure time read from the file header.
        rawdatasec_img : `numpy.ndarray`_
            Data (science) section of the detector as provided by
            setting the (1-indexed) number of the amplifier used to
            read each detector pixel.  Pixels unassociated with any
            amplifier are set to 0.
        oscansec_img : `numpy.ndarray`_
            Overscan section of the detector as provided by setting
            the (1-indexed) number of the amplifier used to read each
            detector pixel.  Pixels unassociated with any amplifier
            are set to 0.
        """
        return super().get_rawimage(raw_file, det,
                                    sec_includes_binning=True)


class INTIDSEEV10Spectrograph(INTIDSSpectrograph):
    """
    Child to handle the blue-sensitive EEV10 detector of INT/IDS.

    This is the camera used by the test data (grating R1200B); its
    detector layout, plate scale, and spectral flip have been verified
    against real frames.
    """
    name = 'int_ids_eev10'
    camera = 'EEV10'
    supported = False
    comment = 'Blue-sensitive EEV10 detector (0.40 arcsec/pixel)'

    # DETECTOR header card value
    detector_name = 'EEV10'
    # Unbinned spatial plate scale [arcsec/pixel] from the IDS web pages
    detector_platescale = 0.40
    # Wavelength decreases with row number (verified against an IRAF
    # solution for the R1200B grating)
    detector_specflip = True
    # Dark current [e-/pixel/hour] from the ING EEV10A detector page
    detector_darkcurr = 4.0


class INTIDSREDPLUS2Spectrograph(INTIDSSpectrograph):
    """
    Child to handle the red-sensitive RED+2 detector of INT/IDS.

    .. warning::

        This class is a placeholder: no RED+2 data have been reduced
        yet, so the exact ``DETECTOR`` header string, amplifier layout,
        gain/read-noise cards, and spectral flip all still need to be
        verified against real RED+2 frames.
    """
    name = 'int_ids_redplus2'
    camera = 'RED+2'
    supported = False
    comment = 'Red-sensitive RED+2 detector (0.44 arcsec/pixel); ' \
              'placeholder, not yet verified against real data'

    # TODO: verify the DETECTOR header card value with real RED+2 data
    detector_name = 'REDPLUS2'
    # Unbinned spatial plate scale [arcsec/pixel] from the IDS web pages
    detector_platescale = 0.44
    # TODO: verify the dispersion direction with real RED+2 data
    detector_specflip = True
    # TODO: verify the dark current with real RED+2 data
    detector_darkcurr = 0.0
