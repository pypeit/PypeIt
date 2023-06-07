"""
Module for MDM/OSMOS specific methods.

.. include:: ../include/links.rst
"""
import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container


class MDMOSMOSMDM4KSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle MDM OSMOS MDM4K instrument+detector
    """
    ndet = 1
    name = 'mdm_osmos_mdm4k'
    telescope = telescopes.KPNOTelescopePar()
    camera = 'MDM4K'
    url = 'https://www.astronomy.ohio-state.edu/martini.10/osmos/'
    header_name = 'OSMOS'
    supported = True
    comment = 'MDM OSMOS spectrometer'

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
        # Detector 1
        detector_dict = dict(
            binning         = '1,1' if hdu is None 
                                    else self.get_meta_value(self.get_headarr(hdu), 'binning'),
            det=1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = True,
            spatflip        = False,
            xgap            = 0.,
            ygap            = 0.,
            ysize           = 1.,
            platescale      = 0.273,
            mincounts       = -1e10,
            darkcurr        = 0.0,
            saturation      = 65535.,
            nonlinear       = 0.86,
            numamplifiers   = 4,
            gain            = np.atleast_1d([2.2, 2.2, 2.2, 2.2]),
            ronoise         = np.atleast_1d([5.0, 5.0, 5.0, 5.0]),
            datasec         = np.atleast_1d(['[9:509,33:2064]', '[509:,33:2064]',
                '[9:509, 2065:4092', '[509:, 2065:4092']),
            oscansec        = np.atleast_1d(['[9:509, 1:32]', '[509:, 1:32]',
                '[9:509, 4098:]', '[509:, 4098:]']),
        )
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

        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        # Wavelength calibration methods
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['lamps'] = ['ArI', 'XeI']
        par['calibrations']['wavelengths']['reid_arxiv'] = 'mdm_osmos_mdm4k.fits'
        par['calibrations']['wavelengths']['sigdetect'] = 10.0
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, None]  # Long arc exposures on this telescope
        par['calibrations']['standardframe']['exprng'] = [None, 120]
        par['scienceframe']['exprng'] = [90, None]

        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(ext=0, card='SLITID')
        self.meta['binning'] = dict(card=None, compound=True)

        self.meta['mjd'] = dict(ext=0, card='MJD')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='SECZ')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='DISPID')
        self.meta['idname'] = dict(ext=0, card='IMAGETYP')
        # Lamps
        self.meta['lampstat01'] = dict(ext=0, card='LAMPS')
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
            binspatial = headarr[0]['CCDXBIN']
            binspec = headarr[0]['CCDYBIN']
            return parse.binning2string(binspec, binspatial)
        else:
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
        return ['DISPID', 'SLITID', 'CCDXBIN', 'CCDYBIN']

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
            return good_exp & (fitstbl['idname'] == 'OBJECT')
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'zero')
        if ftype in ['pixelflat', 'trace']:
            return good_exp & (fitstbl['lampstat01'] == 'Flat') & (fitstbl['idname'] == 'FLAT')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc','tilt']:
            return good_exp & np.array([ilamp in ['Ar','Xe'] for ilamp in fitstbl['lampstat01']]) & (fitstbl['idname'] == 'COMP')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

