""" Module for VLT FORS (1 and 2)
"""
import glob

import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container

from pkg_resources import resource_filename

from pypeit import debugger

class VLTFORSSpectrograph(spectrograph.Spectrograph):
    """
    Child of Spectrograph to handle VLT/FORS specific code
    Parent for FORS1 and FORS2
    """
    ndet = 1  # Because each detector is written to a separate FITS file

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
        par['flexure']['spec_method'] = 'boxcar'

        # Median overscan
        #   IF YOU CHANGE THIS, YOU WILL NEED TO DEAL WITH THE OVERSCAN GOING ALONG ROWS
        for key in par['calibrations'].keys():
            if 'frame' in key:
                par['calibrations'][key]['process']['overscan_method'] = 'median'

        # Adjustments to slit and tilts for NIR
        par['calibrations']['slitedges']['edge_thresh'] = 50.
        par['calibrations']['slitedges']['fit_order'] = 3
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5

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

        # Flats
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.90
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.10

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
        meta['ra'] = dict(ext=0, card='RA', required_ftypes=['science', 'standard'])  # Need to convert to : separated
        meta['dec'] = dict(ext=0, card='DEC', required_ftypes=['science', 'standard'])
        meta['target'] = dict(ext=0, card='OBJECT')
        meta['binning'] = dict(card=None, compound=True)

        meta['mjd'] = dict(ext=0, card='MJD-OBS')
        meta['exptime'] = dict(ext=0, card='EXPTIME')
        meta['airmass'] = dict(ext=0, card='HIERARCH ESO TEL AIRM START', required_ftypes=['science', 'standard'])
        #
        meta['decker'] = dict(card=None, compound=True, required_ftypes=['science', 'standard'])
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
            binning = parse.binning2string(binspec, binspatial)
            return binning
        elif meta_key == 'decker':
            try:  # Science
                decker = headarr[0]['HIERARCH ESO INS SLIT NAME']
            except KeyError:  # Standard!
                try:
                    decker = headarr[0]['HIERARCH ESO SEQ SPEC TARG']
                except KeyError:
                    return None
            return decker
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
            return good_exp & ((fitstbl['target'] == 'STD,FLUX')
                               | (fitstbl['target'] == 'STD'))
        if ftype == 'bias':
            return good_exp & (fitstbl['target'] == 'BIAS')
        if ftype == 'dark':
            return good_exp & (fitstbl['target'] == 'DARK')
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & ((fitstbl['target'] == 'LAMP,DFLAT')
                               | (fitstbl['target'] == 'LAMP,QFLAT')
                               | (fitstbl['target'] == 'FLAT,LAMP')
                               | (fitstbl['target'] == 'LAMP,FLAT'))
        if ftype == 'pinhole':
            # Don't type pinhole
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
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
        self.camera = 'FORS2'
        self.numhead = 1

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

        # Detector 1 (Thor)  -- http://www.eso.org/sci/php/optdet/instruments/fors2/index.html
        detector_dict1 = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.126,  # average between order 11 & 30, see manual
            darkcurr        = 0.0,
            saturation      = 2.0e5,  # I think saturation may never be a problem here since there are many DITs
            nonlinear       = 0.80,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(0.70),
            ronoise         = np.atleast_1d(2.9), # High gain
            datasec         = np.atleast_1d('[11:2059,:]'),  # For 1x binning, I think
            oscansec        = np.atleast_1d('[2062:,:]'),
            #suffix          = '_Thor',
        )
        # Detector 2 (Belenos)
        detector_dict2 = dict(
            binning         = binning,
            det             = 1,  # ESO writes these to separate FITS images!!
            dataext         = 0,
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.126,  # average between order 11 & 30, see manual
            darkcurr        = 0.0,
            saturation      = 2.0e5,  # I think saturation may never be a problem here since there are many DITs
            nonlinear       = 0.80,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(0.70),
            ronoise         = np.atleast_1d(3.0),  # High gain
            datasec         = np.atleast_1d('[20:,0:2048]'),
            oscansec        = np.atleast_1d('[4:20,4:2044]'),
            #suffix          = '_Belenos'
        )

        # Finish
        chip = self.get_meta_value(self.get_headarr(hdu), 'detector')
        if chip == 'CHIP1':
            return detector_container.DetectorContainer(**detector_dict1)
        elif chip == 'CHIP2':
            return detector_container.DetectorContainer(**detector_dict2)


    def default_pypeit_par(self):
        """
        Set default parameters for Keck LRISr reductions.
        """
        par = VLTFORSSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = self.spectrograph

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
        # TODO: Should we allow the user to override these?

        #detector = self.get_meta_value(scifile, 'detector')
        #self.set_detector(detector)
        # Wavelengths
        #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']

        if self.get_meta_value(scifile, 'dispname') == 'GRIS_300I':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_fors2_300I.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'
        elif self.get_meta_value(scifile, 'dispname') == 'GRIS_300V':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_fors2_300V.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'

        if 'lSlit' in self.get_meta_value(scifile, 'decker'):
            par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        return par

    def configuration_keys(self):
        #return ['dispname', 'dispangle', 'decker', 'detector']
        return ['dispname', 'dispangle', 'decker', 'detector']


