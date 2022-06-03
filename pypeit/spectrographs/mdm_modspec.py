"""
Module for MDM/Modspec specific methods.

.. include:: ../include/links.rst
"""
import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container

class MDMModspecEchelleSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle MDM Modspec Echelle instrument+detector
    """
    ndet = 1
    name = 'mdm_modspec_echelle'
    telescope = telescopes.KPNOTelescopePar()
    camera = 'Echelle'
    header_name = 'ModSpec'
    supported = True
    comment = 'MDM Modspec spectrometer'

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
        gain = np.atleast_1d(1.3)      # Hardcoded in the header
        ronoise = np.atleast_1d(7.90)    # Hardcoded in the header
        len1 = hdu[0].header['NAXIS1']
        len2 = hdu[0].header['NAXIS2']
        datasec = np.atleast_1d([
            '[{0:d}:{1:d},{2:d}:{3:d}'.format(1+5, len1-5, 1, len2)])
        oscansec = np.atleast_1d([
            '[{0:d}:{1:d},{2:d}:{3:d}'.format(1, 1+5, 1, len2),
            '[{0:d}:{1:d},{2:d}:{3:d}'.format(len1-5, len1, 1, len2),
        ])
        if hdu is None:
            binning = '1,1'                 # Most common use mode
        else:
            binning = "{0},{1}".format(
                hdu[0].header['CCDBIN1'], hdu[0].header['CCDBIN2'])

        # Detector
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 0,        # Native spectrum is along the x-axis
            specflip        = True,     # DeVeny CCD has blue at the right
            spatflip        = False,
            platescale      = 0.28,     # Arcsec / pixel
            darkcurr        = 0.0,      # Electrons per hour
            saturation      = 65535.,   # 16-bit ADC
            nonlinear       = 0.97,     # Linear to ~97% of saturation
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = gain,     # See above
            ronoise         = ronoise,  # See above
            # Data & Overscan Sections -- Edge tracing can handle slit edges
            datasec         = datasec,  # See above
            oscansec        = oscansec  # See above
            )
        return detector_container.DetectorContainer(**detector_dict)
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
            all of ``PypeIt`` methods.
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

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = 
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
            return good_exp & (fitstbl['imagetyp'] == 'OBJECT')
        if ftype == 'bias':
            return good_exp & (fitstbl['imagetyp'] == 'zero')
        if ftype in ['arc','tilt']:
            return good_exp & np.array([ilamp in ['Ar','Xe'] for ilamp in fitstbl['lampstat01']]) & (fitstbl['idname'] == 'COMP')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

# I have added a comment!