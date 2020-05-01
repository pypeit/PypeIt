""" Module for WHT/ISIS specific codes
"""
import numpy as np

from astropy.io import fits
from pkg_resources import resource_filename

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container

from pypeit import debugger


class WHTISISSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle WHT/ISIS specific code
    """
    ndet = 1

    def __init__(self):
        # Get it started
        super(WHTISISSpectrograph, self).__init__()
        self.spectrograph = 'wht_isis'
        self.telescope = telescopes.WHTTelescopePar()

    def configuration_keys(self):
        """
        Return the metadata keys that defines a unique instrument
        configuration.

        This list is used by :class:`pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            list: List of keywords of data pulled from meta
        """
        return ['dispname', 'decker', 'binning', 'dispangle', 'dichroic']

    def init_meta(self):
        """
        Generate the meta data dict
        Note that the children can add to this

        Returns:
            self.meta: dict (generated in place)

        """
        meta = {}
        # Required (core)
        meta['ra'] = dict(ext=0, card='RA')
        meta['dec'] = dict(ext=0, card='DEC')
        meta['target'] = dict(ext=0, card='OBJECT')
        meta['decker'] = dict(card=None, compound=True)
        meta['binning'] = dict(card=None, compound=True)

        meta['mjd'] = dict(ext=0, card='MJD-OBS')
        meta['exptime'] = dict(ext=0, card='EXPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')
        meta['decker'] = dict(ext=0, card='ISISLITU')
        # Extras for config and frametyping
        meta['dispname'] = dict(ext=0, card='ISIGRAT')
        meta['dichroic'] = dict(ext=0, card='ISIDICHR')
        meta['dispangle'] = dict(ext=0, card='CENWAVE', rtol=1e-3)
        meta['slitwid'] = dict(ext=0, card='ISISLITW')
        meta['idname'] = dict(ext=0, card='IMAGETYP')
        # Lamps
        meta['lampstat01'] = dict(ext=0, card='CAGLAMPS')

        # Ingest
        self.meta = meta

    def compound_meta(self, headarr, meta_key):
        if meta_key == 'binning':
            binspatial = headarr[0]['CCDXBIN']
            binspec = headarr[0]['CCDYBIN']
            return parse.binning2string(binspec, binspatial)
        else:
            msgs.error("Not ready for this compound meta")

    def pypeit_file_keys(self):
        pypeit_keys = super(WHTISISSpectrograph, self).pypeit_file_keys()
        pypeit_keys += ['slitwid']
        return pypeit_keys


class WHTISISBlueSpectrograph(WHTISISSpectrograph):
    """
    Child to handle WHT/ISIS blue specific code
    """
    def __init__(self):
        # Get it started
        super(WHTISISBlueSpectrograph, self).__init__()
        self.spectrograph = 'wht_isis_blue'
        self.camera = 'ISISb'

    def get_detector_par(self, hdu, det):
        """
        Return a DectectorContainer for the current image

        Args:
            hdu (`astropy.io.fits.HDUList`):
                HDUList of the image of interest.
                Ought to be the raw file, or else..
            det (int):

        Returns:
            :class:`pypeit.images.detector_container.DetectorContainer`:

        """
        # Binning
        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')  # Could this be detector dependent??

        # Detector 1
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 1,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.20,
            darkcurr        = 0.0,
            saturation      = 65535.,
            nonlinear       = 0.76,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(1.2),
            ronoise         = np.atleast_1d(5.0),
            datasec         = np.atleast_1d('[:,2:4030]'),
            )
        return detector_container.DetectorContainer(**detector_dict)


    def default_pypeit_par(self):
        """
        Set default parameters for Keck LRISb reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'wht_isis_blue'

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'


        # JFH Is this correct?
        # Processing steps
        turn_off = dict(use_overscan=False)
        par.reset_all_processimages_par(**turn_off)

        # Turn off the overscan
        #for ftype in par['calibrations'].keys():
        #    try:
        #        par['calibrations'][ftype]['process']['overscan'] = 'none'
        #    except (TypeError, KeyError):
        #        pass
        par['scienceframe']['process']['use_overscan'] = False
        # Make a bad pixel mask
        par['calibrations']['bpm_usebias'] = True
        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['pixelflatframe']['process']['sig_lohi'] = [10.,10.]
        # Change the wavelength calibration method
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['lamps'] = ['NeI', 'ArI', 'ArII', 'CuI']
        #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 3
        par['calibrations']['wavelengths']['n_final'] = 5
        par['calibrations']['wavelengths']['sigdetect'] = 10.0
        par['calibrations']['wavelengths']['wv_cen'] = 4859.0
        par['calibrations']['wavelengths']['disp'] = 0.2
        # Do not flux calibrate
        par['fluxcalib'] = None
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, 120]
        par['calibrations']['standardframe']['exprng'] = [None, 120]
        par['scienceframe']['exprng'] = [90, None]

        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        .. todo::
            Document the changes made!

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

        # Wavelength calibrations
        if self.get_meta_value(scifile, 'dispname') == 'R1200B':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'wht_isis_blue_1200_4800.fits'

        # Return
        return par

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['science', 'standard']:
            return good_exp & (fitstbl['lampstat01'] == 'Off') & (fitstbl['idname'] == 'object')
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'zero')
        if ftype in ['pixelflat', 'trace']:
            return good_exp & (fitstbl['lampstat01'] == 'W') & (fitstbl['idname'] == 'flat')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['lampstat01'] == 'CuNe+CuAr') & (fitstbl['idname'] == 'arc')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


class WHTISISRedSpectrograph(WHTISISSpectrograph):
    """
    Child to handle WHT/ISISr red specific code
    """
    def __init__(self):
        # Get it started
        super(WHTISISRedSpectrograph, self).__init__()
        self.spectrograph = 'wht_isis_red'
        self.camera = 'ISISr'

    def get_detector_par(self, hdu, det):
        """
        Return a DectectorContainer for the current image

        Args:
            hdu (`astropy.io.fits.HDUList`):
                HDUList of the image of interest.
                Ought to be the raw file, or else..
            det (int):

        Returns:
            :class:`pypeit.images.detector_container.DetectorContainer`:

        """
        # Binning
        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')  # Could this be detector dependent??

        # Detector 1
        detector_dict = dict(
            binning=binning,
            det=1,
            dataext=1,
            specaxis=0,
            specflip=False,
            spatflip=False,
            platescale=0.22,
            darkcurr=0.0,
            saturation=65535.,
            nonlinear=0.76,
            mincounts=-1e10,
            numamplifiers=1,
            gain=np.atleast_1d(0.98),
            ronoise=np.atleast_1d(4.0),
            datasec=np.atleast_1d('[:,:]'),
        )
        return detector_container.DetectorContainer(**detector_dict)

    def default_pypeit_par(self):
        """
        Set default parameters for WHT ISISr reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'wht_isis_red'

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Turn off the overscan
        for ftype in par['calibrations'].keys():
            try:
                par['calibrations'][ftype]['process']['overscan'] = 'none'
            except (TypeError, KeyError):
                pass
        par['scienceframe']['process']['use_overscan'] = False
        # Make a bad pixel mask
        par['calibrations']['bpm_usebias'] = True
        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['pixelflatframe']['process']['sig_lohi'] = [10.,10.]
        # Change the wavelength calibration method
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['lamps'] = ['NeI', 'ArI', 'ArII', 'CuI']
        #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['sigdetect'] = 10.0
        par['calibrations']['wavelengths']['wv_cen'] = 6000.0
        par['calibrations']['wavelengths']['disp'] = 0.2
        # Do not flux calibrate
        par['fluxcalib'] = None
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, 120]
        par['calibrations']['standardframe']['exprng'] = [None, 120]
        par['scienceframe']['exprng'] = [90, None]

        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        .. todo::
            Document the changes made!

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

        # Wavelength calibrations
        if self.get_meta_value(scifile, 'dispname') == 'R1200R':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'wht_isis_red_1200_6000.fits'

        # Return
        return par

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['science', 'standard']:
            return good_exp & (fitstbl['lampstat01'] == 'Off') & (fitstbl['idname'] == 'object')
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'zero')
        if ftype in ['pixelflat', 'trace']:
            return good_exp & (fitstbl['lampstat01'] == 'W') & (fitstbl['idname'] == 'flat')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['lampstat01'] == 'CuNe+CuAr') & (fitstbl['idname'] == 'arc')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)
