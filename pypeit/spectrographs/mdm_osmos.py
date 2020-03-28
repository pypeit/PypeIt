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
from pypeit.images import detector_container

from IPython import embed


class MDMOSMOSMDM4KSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle MDM OSMOS MDM4K instrument+detector
    """
    ndet = 1

    def __init__(self):
        # Get it started
        super(MDMOSMOSMDM4KSpectrograph, self).__init__()
        self.spectrograph = 'mdm_osmos_mdm4k'
        self.telescope = telescopes.KPNOTelescopePar()
        self.camera = 'MDM4K'

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

        # Detector 1
        detector_dict = dict(
            binning         =self.get_meta_value(self.get_headarr(hdu), 'binning'),
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

    def default_pypeit_par(self):
        """
        Set default parameters for reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'mdm_osmos_mdm4k'

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['pixelflatframe']['process']['sig_lohi'] = [10.,10.]
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
        meta['decker'] = dict(ext=0, card='SLITID')
        meta['binning'] = dict(card=None, compound=True)

        meta['mjd'] = dict(ext=0, card='MJD')
        meta['exptime'] = dict(ext=0, card='EXPTIME')
        meta['airmass'] = dict(ext=0, card='SECZ')
        # Extras for config and frametyping
        meta['dispname'] = dict(ext=0, card='DISPID')
        meta['idname'] = dict(ext=0, card='IMAGETYP')
        # Lamps
        meta['lampstat01'] = dict(ext=0, card='LAMPS')

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
        return ['dispname', 'decker', 'binning']

    def pypeit_file_keys(self):
        pypeit_keys = super(MDMOSMOSMDM4KSpectrograph, self).pypeit_file_keys()
        pypeit_keys += ['slitwid']
        return pypeit_keys

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
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

