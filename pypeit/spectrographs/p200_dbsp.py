""" Module for P200/DBSP specific codes
"""
import numpy as np

from astropy.io import fits
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz, ICRS
from astropy import units as u
from astropy.time import Time
from pkg_resources import resource_filename

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container

from pypeit import debugger
from typing import List

loc = EarthLocation.of_site('Palomar')

def get_zenith_ra_dec(time) -> SkyCoord:
    time = Time(time)
    altaz = AltAz(alt=Angle(90, unit=u.deg), az=Angle(0, unit=u.deg), obstime=time, location=loc)
    return altaz.transform_to(ICRS)

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
        # VERY HACKY!!!!
        meta['ra'] = dict(card=None, compound=True)
        meta['dec'] = dict(card=None, compound=True)
        #meta['ra'] = dict(ext=0, card='RA')
        #meta['dec'] = dict(ext=0, card='DEC')
        meta['target'] = dict(ext=0, card='OBJECT')

        meta['dispname'] = dict(ext=0, card='GRATING')
        meta['decker'] = dict(ext=0, card='APERTURE')
        meta['binning'] = dict(card=None, compound=True)

        meta['mjd'] = dict(card=None, compound=True)
        meta['exptime'] = dict(ext=0, card='EXPTIME')
        #meta['airmass'] = dict(ext=0, card='AIRMASS')
        # VERY HACKY!!!!
        meta['airmass'] = dict(card=None, compound=True)

        # Extras for config and frametyping
        meta['dichroic'] = dict(ext=0, card='DICHROIC')
        meta['dispangle'] = dict(card=None, rtol=1e-2, compound=True)
        meta['slitwid'] = dict(ext=0, card='APERTURE')
        meta['idname'] = dict(ext=0, card='IMGTYPE')
        # Lamps
        meta['lampstat01'] = dict(ext=0, card='LAMPS')

        # Ingest
        self.meta = meta

    def compound_meta(self, headarr: List[fits.Header], meta_key):
        if meta_key == 'mjd':
            return Time(headarr[0]['UTSHUT']).mjd
        elif meta_key == 'dispangle':
            try:
                return Angle(headarr[0]['ANGLE']).deg
            except Exception as e:
                print(headarr[0]['ANGLE'])
                raise e
        elif meta_key == 'ra':
            ra = headarr[0].get('RA', default=None)
            frametype = headarr[0]['IMGTYPE']
            if ra is None and frametype in ['bias', 'flat', 'cal']:
                ra = get_zenith_ra_dec(headarr[0]['UTSHUT']).ra.to_string(unit='hour', sep=':')
            return ra
        elif meta_key == 'dec':
            dec = headarr[0].get('DEC', default=None)
            frametype = headarr[0]['IMGTYPE']
            if dec is None and frametype in ['bias', 'flat', 'cal']:
                dec = get_zenith_ra_dec(headarr[0]['UTSHUT']).dec.to_string(unit='deg', sep=':')
            return dec
        elif meta_key == 'airmass':
            am = headarr[0].get('AIRMASS', default=None)
            frametype = headarr[0]['IMGTYPE']
            if am is None and frametype in ['bias', 'flat', 'cal']:
                am = '1.000'
            elif headarr[0].get('RA', default=None) and headarr[0].get('DEC', default=None):
                ra = headarr[0]['RA']
                dec = headarr[0]['DEC']
                altaz = SkyCoord(ra, dec, unit=(u.hour, u.deg)).transform_to(AltAz(obstime=headarr[0]['UTSHUT'], location=loc))
                am = altaz.secz
            return am
        else:
            return None
            #msgs.error("Not ready for this compound meta")

    def pypeit_file_keys(self):
        pypeit_keys = super(P200DBSPSpectrograph, self).pypeit_file_keys()
        pypeit_keys += ['slitwid']
        return pypeit_keys


class P200DBSPBlueSpectrograph(P200DBSPSpectrograph):
    """
    Child to handle P200/DBSP blue specific code
    """
    def __init__(self):
        # Get it started
        super(P200DBSPBlueSpectrograph, self).__init__()
        self.spectrograph = 'p200_dbsp_blue'
        self.camera = 'DBSPb'
    
    def compound_meta(self, headarr, meta_key):
        retval = super(P200DBSPBlueSpectrograph, self).compound_meta(headarr, meta_key)
        if retval is not None:
            return retval
        if meta_key == 'binning':
            binspatial, binspec = headarr[0]['CCDSUM'].split(' ')
            return parse.binning2string(binspec, binspatial)
        else:
            msgs.error("Not ready for this compound meta")

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
            dataext         = 0, # check
            specaxis        = 0, # this should be wrong
            specflip        = True, # check # Apparently this needs to be true?????? don't like this.....
            spatflip        = False, # check
            platescale      = 0.389,
            darkcurr        = 0.0,
            saturation      = 65000.,
            nonlinear       = 62./65.,
            mincounts       = -1e10, # cross-check
            numamplifiers   = 1,
            gain            = np.atleast_1d(0.72),
            ronoise         = np.atleast_1d(2.5),
            datasec         = np.atleast_1d('[1:2835,1:410]'),       # TODO: from datamodel, should just be able to point this to DSEC1/BSEC1
            oscansec        = np.atleast_1d('[1:2835,411:460]')      # but that just doesn't work, so we will hardcode.
            )
        return detector_container.DetectorContainer(**detector_dict)


    def default_pypeit_par(self):
        """
        Set default parameters for P200 DBSPb reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'p200_dbsp_blue'

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'


        # JFH Is this correct?
        # Processing steps
        # turn_off = dict(use_illumflat=False)
        # par.reset_all_processimages_par(**turn_off)

        # Turn off the overscan
        #for ftype in par['calibrations'].keys():
        #    try:
        #        par['calibrations'][ftype]['process']['overscan'] = 'none'
        #    except (TypeError, KeyError):
        #        pass
        par['scienceframe']['process']['use_overscan'] = True
        # Make a bad pixel mask
        par['calibrations']['bpm_usebias'] = True
        # Set pixel flat combination method
        par['calibrations']['master_dir'] = 'Masters_Blue'
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['pixelflatframe']['process']['sig_lohi'] = [10.,10.]
        # Change the wavelength calibration method
        par['calibrations']['wavelengths']['method'] = 'holy-grail' # ????
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['reid_arxiv'] = 'p200_dbsp_blue.fits'
        par['calibrations']['wavelengths']['lamps'] = ['FeI', 'FeII', 'ArI', 'ArII']

        # Need to set calibrations -> wavelength -> lamps by inspecting the header
        # par['calibrations']['wavelengths']['lamps'] = ['FeI', 'FeII', 'ArI', 'ArII', 'HgI', 'NeI', 'HeI']
        #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        #par['calibrations']['wavelengths']['n_first'] = 3
        #par['calibrations']['wavelengths']['n_final'] = 5
        #par['calibrations']['wavelengths']['sigdetect'] = 10.0
        #par['calibrations']['wavelengths']['wv_cen'] = 4859.0 # set under config specific par
        #par['calibrations']['wavelengths']['disp'] = 0.2 # also set under config specific par
        # Do not flux calibrate
        # par['fluxcalib'] = None
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

        angle_str = self.get_meta_value(scifile, 'dispangle')
        angle = Angle(angle_str, 'degree').rad
        lines_mm = float(self.get_meta_value(scifile, 'dispname').split('/')[0]) # dispname is 600/4000 e.g.

        if lines_mm == 158.:
            order = 2.
        else:
            order = 1.

        theta_m = 38.5 * 2*np.pi / 360. - angle
        par['calibrations']['wavelengths']['wv_cen'] = np.abs(1e7/lines_mm * (np.sin(theta_m) - np.sin(angle)) / order)
        par['calibrations']['wavelengths']['disp'] = np.cos(theta_m)/(lines_mm*1e-7*228.6*order)

        lampstr = self.get_meta_value(scifile, 'lampstat01')
        # Lamp status is D-FeAr-Hg-Ar-Ne-He-InCand
        lamps = [[], ['FeI', 'FeII', 'ArI', 'ArII'], ['HgI'], ['ArI', 'ArII'], ['NeI'], ['HeI'], []]
        lamp_arr = []
        for lamp_idx in range(len(lampstr)):
            lamp_char = lampstr[lamp_idx]
            if lamp_char == "1":
                lamp_arr.extend(lamps[lamp_idx])

        par['calibrations']['wavelengths']['lamps'] = lamp_arr
        par['calibrations']['wavelengths']['lamps'] = ['FeI', 'FeII', 'ArI', 'ArII']

        # Return
        return par

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


class P200DBSPRedSpectrograph(P200DBSPSpectrograph):
    """
    Child to handle P200/DBSPr red specific code
    """
    def __init__(self):
        # Get it started
        super(P200DBSPRedSpectrograph, self).__init__()
        self.spectrograph = 'p200_dbsp_red'
        self.camera = 'DBSPr'
    
    def compound_meta(self, headarr, meta_key):
        retval = super(P200DBSPRedSpectrograph, self).compound_meta(headarr, meta_key)
        if retval is not None:
            return retval
        if meta_key == 'binning':
            binspec, binspatial = headarr[0]['CCDSUM'].split(' ')
            return parse.binning2string(binspec, binspatial)
        else:
            msgs.error("Not ready for this compound meta")

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
            dataext         = 0,
            specaxis        = 1, # aaaaaaaa
            specflip        = False, # check
            spatflip        = False, # check
            platescale      = 0.293,
            darkcurr        = 0.0,
            saturation      = 45000.,
            nonlinear       = 40./45.,
            mincounts       = -1e10, # check
            numamplifiers   = 1,
            gain            = np.atleast_1d(2.8),
            ronoise         = np.atleast_1d(8.5),
            datasec         = np.atleast_1d('[1:440,1:4121]'),       # TODO: from datamodel, should just be able to point this to DSEC1/BSEC1
            oscansec        = np.atleast_1d('[1:440,4122:4141]')     # but that just doesn't work, so we will hardcode.
        )
        return detector_container.DetectorContainer(**detector_dict)

    def default_pypeit_par(self):
        """
        Set default parameters for P200 DBSPr reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'p200_dbsp_red'

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        #par['calibrations']['slitedges']['edge_thresh'] = 'nearest'

        # Turn off the overscan
        #for ftype in par['calibrations'].keys():
        #    try:
        #        par['calibrations'][ftype]['process']['overscan'] = 'none'
        #    except (TypeError, KeyError):
        #        pass
        par['scienceframe']['process']['use_overscan'] = True
        par['scienceframe']['process']['sigclip'] = 4.0 # Tweaked downward from 4.5. 
        par['scienceframe']['process']['objlim'] = 1.5 # Tweaked downward from 3.0. Same value as Keck KCWI and DEIMOS
        # Make a bad pixel mask
        par['calibrations']['bpm_usebias'] = True
        # Set pixel flat combination method
        par['calibrations']['master_dir'] = 'Masters_Red'
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['pixelflatframe']['process']['sig_lohi'] = [10.,10.]
        # Change the wavelength calibration method
        par['calibrations']['wavelengths']['method'] = 'holy-grail' # ????
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['reid_arxiv'] = 'p200_dbsp_red.fits'
        par['calibrations']['wavelengths']['lamps'] = ['ArI', 'ArII', 'NeI', 'HeI']
        #par['calibrations']['wavelengths']['lamps'] = ['FeI', 'FeII', 'ArI', 'ArII', 'HgI', 'NeI', 'HeI']
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

        angle_str = self.get_meta_value(scifile, 'dispangle')
        angle = Angle(angle_str, 'degree').rad
        lines_mm = float(self.get_meta_value(scifile, 'dispname').split('/')[0]) # dispname is 600/4000 e.g.

        order = 1.

        theta_m = 35. * 2*np.pi / 360. - angle
        par['calibrations']['wavelengths']['wv_cen'] = (1e7 / lines_mm) * (np.sin(theta_m) - np.sin(angle)) / order
        par['calibrations']['wavelengths']['disp'] = np.cos(theta_m)/(lines_mm * 1e-7 * 304.8 * order)

        par['calibrations']['wavelengths']['lamps'] = ['ArI', 'ArII', 'NeI', 'HeI']

        # Return
        return par

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
