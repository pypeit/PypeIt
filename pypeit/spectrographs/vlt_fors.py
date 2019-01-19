""" Module for VLT FORS (1 and 2)
"""
from __future__ import absolute_import, division, print_function

import glob

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.core import pixels

from pkg_resources import resource_filename

from pypeit import debugger

class VLTFORSSpectrograph(spectrograph.Spectrograph):
    """
    Child of Spectrograph to handle VLT/FORS specific code
    Parent for FORS1 and FORS2
    """
    def __init__(self):
        # Get it started
        super(VLTFORSSpectrograph, self).__init__()
        self.spectrograph = 'vlt_fors_base'
        self.telescope = telescopes.VLTTelescopePar()

    @property
    def pypeline(self):
        return 'MultiSlit'

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for FORS Reductions
        """
        par = pypeitpar.PypeItPar()

        # Always correct for flexure, starting with default parameters
        par['flexure']['method'] = 'boxcar'

        # Adjustments to slit and tilts for NIR
        par['calibrations']['slits']['sigdetect'] = 50.
        par['calibrations']['slits']['polyorder'] = 3
        par['calibrations']['slits']['maxshift'] = 0.5

        # Tilt parameters
        par['calibrations']['tilts']['tracethresh'] = 25.0
        par['calibrations']['tilts']['spat_order'] = 3
        par['calibrations']['tilts']['spec_order'] = 4

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['HeI', 'ArI']  # Grating dependent
        par['calibrations']['wavelengths']['rms_threshold'] = 0.25
        par['calibrations']['wavelengths']['sigdetect'] = 10.0
        par['calibrations']['wavelengths']['fwhm'] = 4.0  # Good for 2x binning
        par['calibrations']['wavelengths']['n_final'] = 4
        # Reidentification parameters
        #par['calibrations']['wavelengths']['method'] = 'reidentify'
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_xshooter_nir.json'

        # Flats
        par['calibrations']['flatfield']['illumflatten'] = False
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.90
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.10

        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()

        return par

    def init_meta(self):
        """
        Generate the meta data dict
        Note that the children can add to this

        Returns:
            self.meta: dict (generated in place)

        """
        meta = {}
        # Required (core)
        meta['ra'] = dict(card=None, compound=True, required_ftypes=['science', 'standard'])  # Need to convert to : separated
        meta['dec'] = dict(card=None, compound=True, required_ftypes=['science', 'standard'])
        meta['target'] = dict(ext=0, card='OBJECT')
        meta['binning'] = dict(card=None, compound=True)

        meta['mjd'] = dict(ext=0, card='MJD-OBS')
        meta['exptime'] = dict(ext=0, card='EXPTIME')
        meta['airmass'] = dict(ext=0, card='HIERARCH ESO TEL AIRM START', required_ftypes=['science', 'standard'])
        #
        meta['decker'] = dict(ext=0, card='HIERARCH ESO INS SLIT NAME', required_ftypes=['science', 'standard'])
        # Extras for config and frametyping
        meta['dispname'] = dict(ext=0, card='HIERARCH ESO INS GRIS1 NAME', required_ftypes=['science', 'standard'])
        meta['dispangle'] = dict(ext=0, card='HIERARCH ESO INS GRIS1 WLEN', rtol=2.0, required_ftypes=['science', 'standard'])
        meta['idname'] = dict(ext=0, card='HIERARCH ESO DPR CATG')
        meta['detector'] = dict(ext=0, card='EXTNAME')

        # Ingest
        self.meta = meta

    def compound_meta(self, headarr, meta_key):
        if meta_key == 'binning':
            binspatial = headarr[0]['HIERARCH ESO DET WIN1 BINX']
            binspec = headarr[0]['HIERARCH ESO DET WIN1 BINY']
            binning = parse.binning2string(binspatial, binspec)
            return binning
        elif meta_key in ['ra', 'dec']:
            try:  # Calibs do not have RA values
                coord = SkyCoord(ra=headarr[0]['RA'], dec=headarr[0]['DEC'], unit='deg')
            except:
                return None
            if meta_key == 'ra':
                return coord.ra.to_string(unit=units.hour,sep=':',pad=True,precision=2)
            else:
                return coord.dec.to_string(sep=':',pad=True,alwayssign=True,precision=1)
        else:
            msgs.error("Not ready for this compound meta")

    def configuration_keys(self):
        return []

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        # TODO: Allow for 'sky' frame type, for now include sky in
        # 'science' category
        if ftype == 'science':
            return good_exp & ((fitstbl['idname'] == 'SCIENCE')
                                | (fitstbl['target'] == 'STD,TELLURIC')
                                | (fitstbl['target'] == 'STD,SKY'))
        if ftype == 'standard':
            return good_exp & (fitstbl['target'] == 'STD,FLUX')
        if ftype == 'bias':
            return good_exp & (fitstbl['target'] == 'BIAS')
        if ftype == 'dark':
            return good_exp & (fitstbl['target'] == 'DARK')
        if ftype == 'pixelflat' or ftype == 'trace':
            # Flats and trace frames are typed together
            return good_exp & ((fitstbl['target'] == 'LAMP,DFLAT')
                               | (fitstbl['target'] == 'LAMP,QFLAT')
                               | (fitstbl['target'] == 'FLAT,LAMP')
                               | (fitstbl['target'] == 'LAMP,FLAT'))
        if ftype == 'pinhole':
            # Don't type pinhole
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype == 'arc':
            return good_exp & ((fitstbl['target'] == 'LAMP,WAVE')
                               | (fitstbl['target'] == 'WAVE,LAMP'))

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


class VLTFORS2Spectrograph(VLTFORSSpectrograph):
    """
    Child to handle VLT/FORS2 specific code
    """
    def __init__(self):
        # Get it started
        super(VLTFORS2Spectrograph, self).__init__()
        self.spectrograph = 'vlt_fors2'
        self.camera = 'vlt_fors2'
        self.numhead = 1

    def set_detector(self, chip):
        detectors = [
            # Detector 1 (Thor)  -- http://www.eso.org/sci/php/optdet/instruments/fors2/index.html
            pypeitpar.DetectorPar(
                dataext         = 0,
                specaxis        = 1,
                specflip        = False,
                xgap            = 0.,
                ygap            = 0.,
                ysize           = 1.,
                platescale      = 0.126,  # average between order 11 & 30, see manual
                darkcurr        = 0.0,
                saturation      = 2.0e5,  # I think saturation may never be a problem here since there are many DITs
                nonlinear       = 0.80,
                numamplifiers   = 1,
                gain            = 0.70,
                ronoise         = 2.9, # High gain
                datasec         = '[:,10:]',  # For 1x binning, I think
                oscansec        = '[:,0:10]',
                suffix          = '_Thor'),
            # Detector 2 (Belenos)
            pypeitpar.DetectorPar(
                dataext         = 0,
                specaxis        = 1,
                specflip        = False,
                xgap            = 0.,
                ygap            = 0.,
                ysize           = 1.,
                platescale      = 0.126,  # average between order 11 & 30, see manual
                darkcurr        = 0.0,
                saturation      = 2.0e5,  # I think saturation may never be a problem here since there are many DITs
                nonlinear       = 0.80,
                numamplifiers   = 1,
                gain            = 0.70,
                ronoise         = 3.0,  # High gain
                datasec         = '[20:,0:2048]',
                oscansec        = '[4:20,4:2044]',
                suffix          = '_Belenos'
                )]
        if chip == 'CHIP1':
            self.detector = [detectors[0]]
        elif chip == 'CHIP2':
            debugger.set_trace()  # NEED TO SET DATASEC
            self.detector = [detectors[1]]


    def default_pypeit_par(self):
        """
        Set default parameters for Keck LRISr reductions.
        """
        par = VLTFORSSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = self.spectrograph

        return par

    def config_specific_par(self, par, scifile):
        detector = self.get_meta_value(scifile, 'detector')
        self.set_detector(detector)
        # Wavelengths
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']

        if self.get_meta_value(scifile, 'dispname') == 'GRIS_300I':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_fors2_300I.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'

        return par

    def configuration_keys(self):
        #return ['dispname', 'dispangle', 'decker', 'detector']
        return ['dispname', 'dispangle', 'decker', 'detector']

    def compound_meta(self, headarr, meta_key):
        if meta_key == 'binning':
            binspatial = headarr[0]['HIERARCH ESO DET WIN1 BINX']
            binspec = headarr[0]['HIERARCH ESO DET WIN1 BINY']
            binning = parse.binning2string(binspatial, binspec)
            return binning
        elif meta_key in ['ra', 'dec']:
            try:  # Calibs do not have RA values
                coord = SkyCoord(ra=headarr[0]['RA'], dec=headarr[0]['DEC'], unit='deg')
            except:
                return None
            if meta_key == 'ra':
                return coord.ra.to_string(unit=units.hour,sep=':',pad=True,precision=2)
            else:
                return coord.dec.to_string(sep=':',pad=True,alwayssign=True,precision=1)
        else:
            msgs.error("Not ready for this compound meta")


    def bpm(self, shape=None, filename=None, det=None, **null_kwargs):
        """
        Override parent bpm function with BPM specific to X-ShooterNIR.

        .. todo::
            Allow for binning changes.

        Parameters
        ----------
        det : int, REQUIRED
        **null_kwargs:
            Captured and never used

        Returns
        -------
        bpix : ndarray
          0 = ok; 1 = Mask

        """
        self.empty_bpm(shape=shape, filename=filename, det=det)
        return self.bpm_img



