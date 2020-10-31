""" Module for MMT/Blue Channel specific codes
"""
import glob
import numpy as np
from astropy.io import fits
from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container

from pkg_resources import resource_filename

class MMTBlueChannelSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle MMT/Blue Channel specific code
    """
    ndet = 1

    def __init__(self):
        # Get it started
        super(MMTBlueChannelSpectrograph, self).__init__()
        self.spectrograph = 'mmt_bluechannel'
        self.telescope = telescopes.MMTTelescopePar()
        self.camera = 'Blue Channel'

    def get_detector_par(self, hdu, det):
        header = hdu[0].header

        # Binning
        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')

        # Detector 1
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            xgap            = 0.,
            ygap            = 0.,
            ysize           = 1.,
            platescale      = 0.3,
            darkcurr        = header['DARKCUR'],
            saturation      = 65535.,
            nonlinear       = 0.95,  # need to look up and update
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(header['GAIN']),
            ronoise         = np.atleast_1d(header['RDNOISE']),
            datasec         = np.atleast_1d(header['DATASEC']),
            oscansec        = np.atleast_1d(header['BIASSEC'])
        )

        return detector_container.DetectorContainer(**detector_dict)

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
        meta['decker'] = dict(ext=0, card='APERTURE')
        meta['dichroic'] = dict(ext=0, card='INSFILTE')
        meta['binning'] = dict(ext=0, card=None, compound=True)
        meta['mjd'] = dict(ext=0, card=None, compound=True)
        meta['exptime'] = dict(ext=0, card='EXPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')

        # Extras for config and frametyping
        meta['dispname'] = dict(ext=0, card='DISPERSE')
        meta['idname'] = dict(ext=0, card='IMAGETYP')

        # used for arc and continuum lamps
        meta['lampstat01'] = dict(ext=0, card=None, compound=True)

        # Ingest
        self.meta = meta


    def compound_meta(self, headarr, meta_key):
        """

        Args:
            headarr: list
            meta_key: str

        Returns:
            value

        """
        if meta_key == 'binning':
            """
            Binning in blue channel headers is space-separated rather than comma-separated.
            """
            binspec, binspatial = headarr[0]['CCDSUM'].split()
            binning = parse.binning2string(binspec, binspatial)
            return binning
        elif meta_key == 'mjd':
            """
            Need to combine 'DATE-OBS' and 'UT' headers and then use astropy to make an mjd.
            """
            date = headarr[0]['DATE-OBS']
            ut = headarr[0]['UT']
            ttime = Time(f"{date}T{ut}", format='isot')
            return ttime.mjd
        elif meta_key == 'lampstat01':
            """
            If the comparison mirror is in, there will be a 'COMPLAMP' header entry containing the lamps
            that are turned on. However, if the comparison mirror is out, then this header entry doesn't exist.
            So need to test for it and set to 'Off' if it's not there.
            """
            if 'COMPLAMP' in headarr[0]:
                return headarr[0]['COMPLAMP']
            else:
                return 'off'
        else:
            msgs.error(f"Not ready for compound meta, {meta_key}, for MMT Blue Channel.")


    def default_pypeit_par(self):
        """
        Set default parameters for MMT/Blue Channel reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'mmt_bluechannel'
        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.5
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['fwhm']= 5.0
        par['calibrations']['wavelengths']['lamps'] = ['ArI', 'ArII', 'HeI', 'NeI', 'ThAr']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'

        # Tilt and slit parameters
        par['calibrations']['tilts']['tracethresh'] =  10.0
        par['calibrations']['tilts']['spat_order'] = 6
        par['calibrations']['tilts']['spec_order'] = 6
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Processing steps
        turn_off = dict(use_biasimage=False, use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0
        ## Do not perform global sky subtraction for standard stars
        par['reduce']['skysub']['global_sky_std']  = False

        # Flexure
        par['flexure']['spec_method'] = 'boxcar'

        # cosmic ray rejection parameters for science frames
        par['scienceframe']['process']['sigclip'] = 5.0
        par['scienceframe']['process']['objlim'] = 2.0

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 600]
        par['calibrations']['arcframe']['exprng'] = [10, None]
        par['calibrations']['darkframe']['exprng'] = [300, None]
        par['scienceframe']['exprng'] = [10, None]

        # Sensitivity function parameters
        par['sensfunc']['polyorder'] = 7
        par['sensfunc']['IR']['telgridfile'] = resource_filename('pypeit', '/data/telluric/TelFit_MaunaKea_3100_26100_R20000.fits')

        return par

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a BPM

        Parameters
        ----------
        det : int, REQUIRED
        **null_kwargs:
           Captured and never used

        Returns
        -------
        badpix : ndarray

        """
        # Get the empty bpm: force is always True
        bpm_img = self.empty_bpm(filename, det, shape=shape)

        # Fill in bad pixels if a master bias frame is provided
        if msbias is not None:
            return self.bpm_frombias(msbias, det, bpm_img)

        if det == 1:
            msgs.info("Using hard-coded BPM for  Blue Channel")

            bpm_img[:, -1] = 1

        else:
            msgs.error(f"Invalid detector number, {det}, for MMT Blue Channel (only one detector).")

        return bpm_img

    def configuration_keys(self):
        return ['dispname']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype == 'science':
            return good_exp & (fitstbl['lampstat01'] == 'off') & (fitstbl['idname'] == 'object')
        if ftype == 'standard':
            return good_exp & (fitstbl['lampstat01'] == 'off') & (fitstbl['idname'] == 'object')
        if ftype in ['arc', 'tilt']:
            # should flesh this out to include all valid arc lamp combos
            return good_exp & (fitstbl['lampstat01'] != 'off') & (fitstbl['idname'] == 'comp') & (fitstbl['decker'] != '5.0x180')
        if ftype in ['trace']:
            # i think the bright lamp, BC, is the only one ever used for this
            return good_exp & (fitstbl['lampstat01'] == 'BC')
        if ftype in ['pixelflat', 'illumflat']:
            # imagetyp should always be set to flat for these
            return good_exp & (fitstbl['idname'] == 'flat')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)
