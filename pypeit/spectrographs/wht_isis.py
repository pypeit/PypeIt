""" Module for Shane/Kast specific codes
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

from pypeit import debugger

#class WhtIsisSpectrograph(spectrograph.Spectrograph):
#    """
#    Child to handle Shane/Kast specific code
#    """
#
#    def __init__(self):
#        super(WhtIsisSpectrograph, self).__init__()
#        self.spectrograph = 'wht_isis_base'
#        self.telescope = telescopes.WHTTelescopePar()
#
#    def metadata_keys(self):
#        return super(KeckLRISSpectrograph, self).metadata_keys() \
#                    + ['binning', 'dichroic', 'dispangle']

class WHTISISBlueSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle WHT/ISIS blue specific code
    """
    def __init__(self):
        # Get it started
        super(WHTISISBlueSpectrograph, self).__init__()
        self.spectrograph = 'wht_isis_blue'
        self.telescope = telescopes.WHTTelescopePar()
        self.camera = 'ISISb'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 1,
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.225,
                            darkcurr        = 0.0,
                            saturation      = 65535.,
                            nonlinear       = 0.76,
                            numamplifiers   = 1,
                            gain            = 1.2,
                            ronoise         = 5.0,
                            datasec         = '[:,2:4030]',
                            oscansec        = None,
                            suffix          = '_blue'
                            )]
        self.numhead = 2
        # Uses default timeunit
        # Uses default primary_hdrext
        # self.sky_file = ?

    def default_pypeit_par(self):
        """
        Set default parameters for Keck LRISb reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'wht_isis_blue'

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Turn off the overscan
        for ftype in par['calibrations'].keys():
            try:
                par['calibrations'][ftype]['process']['overscan'] = 'none'
            except (TypeError, KeyError):
                pass
        par['scienceframe']['process']['overscan'] = 'none'
        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['pixelflatframe']['process']['sig_lohi'] = [10.,10.]
        # Change the wavelength calibration method
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['lamps'] = ['NeI', 'ArI', 'ArII', 'CuI']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 3
        par['calibrations']['wavelengths']['n_final'] = 5
        par['calibrations']['wavelengths']['sigdetect'] = 10.0
        par['calibrations']['wavelengths']['wv_cen'] = 4859.0
        par['calibrations']['wavelengths']['disp'] = 0.2
        # Do not flux calibrate
        par['fluxcalib'] = None
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
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

    def pypeit_file_keys(self):
        pypeit_keys = super(WHTISISBlueSpectrograph, self).pypeit_file_keys()
        pypeit_keys += ['slitwid']
        return pypeit_keys

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

    def bpm(self, filename=None, det=None, shape=None, msbias=None, **null_kwargs):
        """ Generate a BPM

        Parameters
        ----------
        shape : tuple, REQUIRED
        filename : str, REQUIRED for binning
        det : int, REQUIRED
        **null_kwargs:
           Captured and never used

        Returns
        -------
        badpix : ndarray

        """
        # Get the empty bpm: force is always True
        #import pdb
        #pdb.set_trace()
        self.bpm_img = self.empty_bpm(filename, det=det, shape=shape)

        # Only defined for det=2
        if msbias is not None:
            msgs.info("Generating a BPM for det={0:d} on ISISb".format(det))
            medval = np.median(msbias.image)
            madval = 1.4826 * np.median(np.abs(medval - msbias.image))
            ww = np.where(np.abs(msbias.image-medval) > 10.0*madval)
            self.bpm_img[ww] = 1

        return self.bpm_img
