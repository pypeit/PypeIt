"""
Module for Gemini flamingos

.. include: ../include/links.rst
"""
import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.images import detector_container
from pypeit.spectrographs import spectrograph
from pkg_resources import resource_filename


class GeminiFLAMINGOSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Gemini/Flamingos specific code
    """
    ndet = 1
    # NOTE: This is a base class for the two main spectrograph classes below,
    # so this one does not have a name.
    name = None
    telescope = telescopes.GeminiSTelescopePar()

    def __init__(self):
        super().__init__()

    def init_meta(self):
        """
        Generate the metadata dictionary.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(ext=0, card='MASKNAME')
        self.meta['dichroic'] = dict(ext=0, card='FILTER')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')

        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='GRISM')
        self.meta['idname'] = dict(ext=0, card='OBSTYPE')


class GeminiFLAMINGOS2Spectrograph(GeminiFLAMINGOSSpectrograph):
    """
    Child to handle Gemini/Flamingos2 Echelle data
    """
    name = 'gemini_flamingos2'
    camera = 'FLAMINGOS'

    def __init__(self):
        super().__init__()

    # TODO: Can these be made static or class methods?
    def get_detector_par(self, hdu, det):
        """
        Return a DectectorContainer for the current image

        Args:
            hdu (`astropy.io.fits.HDUList`_):
                HDUList of the image of interest.
                Ought to be the raw file, or else..
            det (int):
                1-indexed detector number.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
        """
        # Detector 1
        detector_dict = dict(
            binning         = '1,1',
            det             = 1,
            dataext         = 1,
            specaxis        = 0,
            specflip        = True,
            spatflip        = False,
            platescale      = 0.1787,
            darkcurr        = 0.5,
            saturation      = 700000., #155400.,
            nonlinear       = 1.0,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(4.44),
            ronoise         = np.atleast_1d(5.0), #8 CDS read
            datasec         = np.atleast_1d('[:,:]'),
            oscansec        = np.atleast_1d('[:,:]'),
            )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Set default parameters for the reductions.
        """
        par = super().default_pypeit_par()

        # Image processing steps
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Wavelengths
        # 1D wavelength solution with arc lines
        par['calibrations']['wavelengths']['rms_threshold'] = 0.5
        par['calibrations']['wavelengths']['sigdetect']=5
        par['calibrations']['wavelengths']['fwhm'] = 5
        par['calibrations']['wavelengths']['n_first']=2
        par['calibrations']['wavelengths']['n_final']=4
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        par['calibrations']['wavelengths']['match_toler']=5.0

        # Set slits and tilts parameters
        par['calibrations']['tilts']['tracethresh'] = 5
        par['calibrations']['tilts']['spat_order'] = 4
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['edge_thresh'] = 200.
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.4
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 30]
        par['calibrations']['tiltframe']['exprng'] = [50, None]
        par['calibrations']['arcframe']['exprng'] = [50, None]
        par['calibrations']['darkframe']['exprng'] = [20, None]
        par['scienceframe']['exprng'] = [20, None]

        # Scienceimage parameters
        par['reduce']['findobj']['sig_thresh'] = 5.0
        par['reduce']['skysub']['sky_sigrej'] = 5.0
        par['reduce']['findobj']['find_trim_edge'] = [10,10]
        # Do not correct for flexure
        par['flexure']['spec_method'] = 'skip'

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 8
        # TODO: replace the telluric grid file for Gemini-S site.
        par['sensfunc']['IR']['telgridfile'] \
                = resource_filename('pypeit',
                                    '/data/telluric/TelFit_LasCampanas_3100_26100_R20000.fits')

        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the ``PypeIt`` parameters to hard-wired values used for
        specific instrument configurations.

        .. todo::
            Document the changes made!

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
        par = self.__class__.default_pypeit_par() if inp_par is None else inp_par
        # TODO: Should we allow the user to override these?

        if self.get_meta_value(scifile, 'dispname') == 'JH_G5801':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'Flamingos2_JH_JH.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'HK_G5802':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'Flamingos2_HK_HK.fits'
        return par

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.

        Args:
            ftype (:obj:`str`):
                Type of frame to check.  Must be a valid frame type.

        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['pinhole', 'bias']:
            # No pinhole or bias frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['pixelflat', 'trace']:
            return good_exp & (fitstbl['idname'] == 'FLAT')
        if ftype == 'standard':
            return good_exp & (fitstbl['idname'] == 'OBJECT')
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'OBJECT')
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'OBJECT')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


class GeminiFLAMINGOS1Spectrograph(GeminiFLAMINGOSSpectrograph):
    """
    TODO: Place holder, NOT works yet.
    """

    name = 'gemini_flamingos1'
    camera = None
    # TODO: Assign the camera name when the spectrograph is supported!
    #camera = 'FLAMINGOS'

    def __init__(self):
        super().__init__()

    def get_detector_par(self, hdu, det):
        """
        Return a DectectorContainer for the current image

        Args:
            hdu (`astropy.io.fits.HDUList`):
                HDUList of the image of interest.
                Ought to be the raw file, or else..
            det (int):
                1-indexed detector number.

        Returns:
            :class:`pypeit.images.detector_container.DetectorContainer`:

        """
        # Detector 1
        detector_dict = dict(
            binning='1,1',
            det             = 1,
            dataext         = 1,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.15,
            darkcurr        = 0.01,
            saturation      = 320000., #155400.,
            nonlinear       = 0.875,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(3.8),
            ronoise         = np.atleast_1d(6.0), # SUTR readout
            datasec= np.atleast_1d('[5:2044, 900:1250]'),
            oscansec= np.atleast_1d('[:5, 900:1250]'),
            )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Set default parameters.
        """
        par = super().default_pypeit_par()

        # Image processing steps
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Wavelengths
        # 1D wavelength solution with arc lines
        par['calibrations']['wavelengths']['rms_threshold'] = 1.0
        par['calibrations']['wavelengths']['sigdetect']=3
        par['calibrations']['wavelengths']['fwhm'] = 20
        par['calibrations']['wavelengths']['n_first']=2
        par['calibrations']['wavelengths']['n_final']=4
        par['calibrations']['wavelengths']['lamps'] = ['ArI', 'ArII', 'ThAr', 'NeI']
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['reid_arxiv'] = 'magellan_fire_long.fits'
        par['calibrations']['wavelengths']['match_toler']=5.0

        # Set slits and tilts parameters
        par['calibrations']['tilts']['tracethresh'] = 5
        par['calibrations']['slitedges']['trace_thresh'] = 5.
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Scienceimage parameters
        par['reduce']['findobj']['sig_thresh'] = 5.0
        # TODO: I think this parameter was removed
        par['reduce']['findobj']['find_trim_edge'] = [50,50]

        # Do not correct for flexure
        par['flexure']['spec_method'] = 'skip'

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 60]
        par['calibrations']['arcframe']['exprng'] = [1, 50]
        par['calibrations']['darkframe']['exprng'] = [20, None]
        par['scienceframe']['exprng'] = [20, None]
        return par

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['pinhole', 'bias']:
            # No pinhole or bias frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['pixelflat', 'trace']:
            return good_exp & (fitstbl['idname'] == 'PixFlat')
        if ftype == 'standard':
            return good_exp & (fitstbl['idname'] == 'Telluric')
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'Science')
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'Arc')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

