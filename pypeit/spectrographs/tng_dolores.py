"""
Module for TNG/Dolores

.. include:: ../include/links.rst
"""

import numpy as np

from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container


class TNGDoloresSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Shane/Kast specific code
    """
    ndet = 1
    name = 'tng_dolores'
    telescope = telescopes.TNGTelescopePar()
    camera = 'DOLORES'
    url = 'https://oapd.inaf.it/mos/'
    comment = 'DOLORES (LRS) spectrograph; LR-R'


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
            binning='1,1',
            det             = 1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            xgap            = 0.,
            ygap            = 0.,
            ysize           = 1.,
            platescale      = 0.252,
            darkcurr        = 0.0,
            saturation      = 65500.,
            nonlinear       = 0.99,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(0.97),
            ronoise         = np.atleast_1d(9.0),
            datasec         = np.atleast_1d('[1:2045,51:]'),
            oscansec        = np.atleast_1d('[2054:,51:]'),
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

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 0.1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['scienceframe']['exprng'] = [1, None]
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['arcframe']['process']['clip'] = False
        par['calibrations']['arcframe']['process']['combine'] = 'mean'
        par['calibrations']['arcframe']['process']['subtract_continuum'] = True
        par['calibrations']['tiltframe']['process']['clip'] = False
        par['calibrations']['tiltframe']['process']['combine'] = 'mean'
        par['calibrations']['tiltframe']['process']['subtract_continuum'] = True

        return par

    
    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            scifile (str):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`pypeit.par.parset.ParSet`: The PypeIt paramter set
            adjusted for configuration specific parameter values.
        """
        
        par = self.default_pypeit_par() if inp_par is None else inp_par

        if self.get_meta_value(scifile, 'dispname') == 'LR-B':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'tng_dolores_LR-B_arx_v2.fits'
            # Add CdI
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['lamps'] = ['NeI', 'HgI', 'HeI']
        else:
            par['calibrations']['wavelengths']['method'] = 'holy-grail'
            msgs.warn('Check wavelength calibration file.')

        # Return
        return par
    

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card=None, compound=True)
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='OBJCAT')
        self.meta['decker'] = dict(ext=0, card='SLT_ID')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')

        self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='GRM_ID')
        #self.meta['dispangle'] = dict(card=None, compound=True, rtol=1e-5)
        self.meta['idname'] = dict(ext=0, card='OBS-TYPE')
        # Lamps
        self.meta['lampstat01'] = dict(ext=0, card='LMP_ID')

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
        if meta_key == 'mjd':
            time = headarr[0]['DATE-OBS']
            ttime = Time(time, format='isot')
            return ttime.mjd
        elif meta_key == 'ra':
            radeg = headarr[0]['RA-RAD']*180.0/np.pi  # Convert radians to decimal degrees
            return radeg
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
        return ['dispname', 'decker']

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
        return ['GRM_ID', 'SLT_ID']

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
            return good_exp & (fitstbl['idname'] == 'OBJECT') & (fitstbl['lampstat01'] == 'Parking') \
                        & (fitstbl['dispname'] != 'OPEN')
        if ftype == 'bias':
            return good_exp & (fitstbl['dispname'] == 'OPEN')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (fitstbl['idname'] == 'CALIB') & (fitstbl['lampstat01'] == 'Halogen') \
                        & (fitstbl['dispname'] != 'OPEN')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'CALIB') & ( (fitstbl['lampstat01'] == 'Ne+Hg') | (fitstbl['lampstat01'] == 'Helium') ) \
                        & (fitstbl['dispname'] != 'OPEN')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


