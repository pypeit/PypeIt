""" Module for P200/DBSP specific codes
"""
import numpy as np

from astropy.io import fits
from astropy.coordinates import Angle
from astropy import units as u
from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container

from typing import List
from pkg_resources import resource_filename

def flip_fits_slice(s: str) -> str:
    return '[' + ','.join(s.strip('[]').split(',')[::-1]) + ']'

class P200DBSPSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle P200/DBSP specific code
    """
    ndet = 1

    def __init__(self):
        # Get it started
        super(P200DBSPSpectrograph, self).__init__()
        self.spectrograph = 'p200_dbsp'
        self.telescope = telescopes.P200TelescopePar()

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
        return ['dispname', 'binning', 'dispangle', 'dichroic']

    def init_meta(self):
        """
        Generate the meta data dict
        Note that the children can add to this

        Returns:
            self.meta: dict (generated in place)

        """
        meta = {}
        # Required (core)
        meta['ra'] = dict(ext=0, card='RA', required_ftypes=['science', 'standard'])
        meta['dec'] = dict(ext=0, card='DEC', required_ftypes=['science', 'standard'])
        meta['target'] = dict(ext=0, card='OBJECT')

        meta['dispname'] = dict(ext=0, card='GRATING')
        meta['decker'] = dict(ext=0, card='APERTURE')
        meta['binning'] = dict(card=None, compound=True)

        meta['mjd'] = dict(card=None, compound=True)
        meta['exptime'] = dict(ext=0, card='EXPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS', required_ftypes=['science', 'standard'])

        # Extras for config and frametyping
        meta['dichroic'] = dict(ext=0, card='DICHROIC')
        meta['dispangle'] = dict(card=None, rtol=1e-2, compound=True)
        meta['slitwid'] = dict(ext=0, card='APERTURE')
        meta['idname'] = dict(ext=0, card='IMGTYPE')
        # Lamps
        meta['lampstat01'] = dict(ext=0, card='LAMPS')

        # Ingest
        self.meta = meta

    def compound_meta(self, headarr: List[fits.Header], meta_key: str):
        """
        Methods to generate meta in a more complex manner than simply
        reading from the header.

        mjd is converted from UTSHUT header
        dispangle is parsed from ANGLE header

        Args:
            headarr: List[fits.Header]
              List of headers
            meta_key: str

        Returns:
            value:

        """
        if meta_key == 'mjd':
            return Time(headarr[0]['UTSHUT']).mjd
        elif meta_key == 'dispangle':
            try:
                return Angle(headarr[0]['ANGLE']).deg
            except Exception as e:
                msgs.warn("Could not read dispangle from header:" + msgs.newline() + str(headarr[0]['ANGLE']))
                raise e
        else:
            return None

    def pypeit_file_keys(self):
        pypeit_keys = super(P200DBSPSpectrograph, self).pypeit_file_keys()
        pypeit_keys += ['slitwid']
        return pypeit_keys
    
    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
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
    def __init__(self):
        # Get it started
        super(P200DBSPBlueSpectrograph, self).__init__()
        self.spectrograph = 'p200_dbsp_blue'
        self.camera = 'DBSPb'
    
    def compound_meta(self, headarr: List[fits.Header], meta_key: str):
        """
        Methods to generate meta in a more complex manner than simply
        reading from the header. Super method handles mjd and dispangle

        binning is parsed from CCDSUM header

        Args:
            headarr: List[fits.Header]
              List of headers
            meta_key: str

        Returns:
            value:

        """
        # Handle dispangle and mjd from superclass method
        retval = super(P200DBSPBlueSpectrograph, self).compound_meta(headarr, meta_key)

        # If superclass could not handle the meta key
        if retval is not None:
            return retval
        if meta_key == 'binning':
            binspatial, binspec = headarr[0]['CCDSUM'].split(' ')
            return parse.binning2string(binspec, binspatial)
        msgs.error("Not ready for this compound meta")

    def get_detector_par(self, hdu: fits.HDUList, det: int):
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
            ronoise         = np.atleast_1d(2.5)
            )
        
        header = hdu[0].header
        datasec = header['TSEC1']
        oscansec = header['BSEC1']

        detector_dict['datasec'] = np.atleast_1d(flip_fits_slice(datasec))
        detector_dict['oscansec'] = np.atleast_1d(flip_fits_slice(oscansec))

        return detector_container.DetectorContainer(**detector_dict)


    def default_pypeit_par(self):
        """
        Set default parameters for P200 DBSPb reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'p200_dbsp_blue'

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.55


        par['scienceframe']['process']['use_overscan'] = True
        # Make a bad pixel mask
        par['calibrations']['bpm_usebias'] = True
        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['pixelflatframe']['process']['sig_lohi'] = [10.,10.]
        # Change the wavelength calibration method
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['lamps'] = ['FeI', 'FeII', 'ArI', 'ArII']

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

        disp = self.get_meta_value(scifile, 'dispname')
        if disp == '600/4000':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'p200_dbsp_blue_600_4000_d55.fits'
        else:
            msgs.error("Your grating " + disp + ' needs a template spectrum for the blue arm of DBSP.')

        angle = Angle(self.get_meta_value(scifile, 'dispangle'), unit=u.deg).rad
        slitwidth = self.get_meta_value(scifile, 'slitwid') * u.arcsec
        lines_mm = float(self.get_meta_value(scifile, 'dispname').split('/')[0]) / u.mm

        theta_m = 38.5 * 2*np.pi / 360. - angle
        order = 2. if lines_mm == 158. * u.mm else 1.
        platescale = 0.389 * u.arcsec / u.pix
        pix_size = 15 * u.um

        disp = np.cos(theta_m)/(lines_mm * 9 * u.imperial.inch) * 1e7 * u.AA / u.mm
        cen_wv = np.abs(1/lines_mm * (np.sin(theta_m) - np.sin(angle)) / order)
        dlam = slitwidth / platescale * pix_size / u.pix * disp

        resolving_power = cen_wv / dlam

        par['sensfunc']['UVIS']['resolution'] = resolving_power.decompose().value
        
        return par


class P200DBSPRedSpectrograph(P200DBSPSpectrograph):
    """
    Child to handle P200/DBSPr red specific code
    """
    def __init__(self):
        # Get it started
        super(P200DBSPRedSpectrograph, self).__init__()
        self.spectrograph = 'p200_dbsp_red'
        self.camera = 'DBSPr'
    
    def compound_meta(self, headarr: List[fits.Header], meta_key: str):
        """
        Methods to generate meta in a more complex manner than simply
        reading from the header. Super method handles mjd and dispangle

        binning is parsed from CCDSUM header

        Args:
            headarr: List[fits.Header]
              List of headers
            meta_key: str

        Returns:
            value:

        """
        # Handle dispangle and mjd from superclass method
        retval = super(P200DBSPRedSpectrograph, self).compound_meta(headarr, meta_key)
        
        # If superclass could not handle the meta key
        if retval is not None:
            return retval
        if meta_key == 'binning':
            binspec, binspatial = headarr[0]['CCDSUM'].split(' ')
            return parse.binning2string(binspec, binspatial)
        else:
            msgs.error("Not ready for this compound meta")

    def get_detector_par(self, hdu: fits.HDUList, det: int):
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
            ronoise         = np.atleast_1d(8.5)
        )

        header = hdu[0].header

        datasec = header['TSEC1']
        oscansec = header['BSEC1']

        detector_dict['datasec'] = np.atleast_1d(flip_fits_slice(datasec))
        detector_dict['oscansec'] = np.atleast_1d(flip_fits_slice(oscansec))

        return detector_container.DetectorContainer(**detector_dict)

    def default_pypeit_par(self):
        """
        Set default parameters for P200 DBSPr reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'p200_dbsp_red'

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'


        par['scienceframe']['process']['use_overscan'] = True
        par['scienceframe']['process']['sigclip'] = 4.0 # Tweaked downward from 4.5. 
        par['scienceframe']['process']['objlim'] = 1.5 # Tweaked downward from 3.0. Same value as Keck KCWI and DEIMOS
        # Make a bad pixel mask
        par['calibrations']['bpm_usebias'] = True
        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['pixelflatframe']['process']['sig_lohi'] = [10.,10.]
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
        par['sensfunc']['IR']['telgridfile'] = resource_filename('pypeit', '/data/telluric/TelFit_Lick_3100_11100_R10000.fits')
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

        disp = self.get_meta_value(scifile, 'dispname')
        if disp == '316/7500':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'p200_dbsp_red_316_7500_d55.fits'
        else:
            msgs.error("Your grating " + disp + ' needs a template spectrum for the red arm of DBSP.')


        angle = Angle(self.get_meta_value(scifile, 'dispangle'), unit=u.deg).rad
        slitwidth = self.get_meta_value(scifile, 'slitwid') * u.arcsec
        lines_mm = float(self.get_meta_value(scifile, 'dispname').split('/')[0]) / u.mm

        theta_m = 35.0 * 2*np.pi / 360. - angle
        order = 1.
        platescale = 0.293 * u.arcsec / u.pix
        pix_size = 15 * u.um

        disp = np.cos(theta_m)/(lines_mm * 12 * u.imperial.inch) * 1e7 * u.AA / u.mm
        cen_wv = np.abs(1/lines_mm * (np.sin(theta_m) - np.sin(angle)) / order)
        dlam = slitwidth / platescale * pix_size / u.pix * disp

        resolving_power = cen_wv / dlam

        par['sensfunc']['UVIS']['resolution'] = resolving_power.decompose().value

        return par
