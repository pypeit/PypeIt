""" Module for TNG/Dolores
"""
import numpy as np

from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container

from pypeit import debugger

class TNGDoloresSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Shane/Kast specific code
    """
    ndet = 1

    def __init__(self):
        super(TNGDoloresSpectrograph, self).__init__()
        self.spectrograph = 'tng_dolores'
        self.telescope = telescopes.TNGTelescopePar()
        self.camera = 'DOLORES'
        self.timeunit = 'isot'


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
            binning='1,1',
            det             = 1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            xgap            = 0.,
            ygap            = 0.,
            ysize           = 1.,
            platescale      = 0.252,
            darkcurr        = 0.0,
            saturation      = 65535.,
            nonlinear       = 0.76,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(0.97),
            ronoise         = np.atleast_1d(9.0),
            datasec         = np.atleast_1d('[51:,1:2045]'),
            oscansec        = np.atleast_1d('[51:,2054:]'),
            )
        return detector_container.DetectorContainer(**detector_dict)

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for TNG Dolores reductions.
        """
        par = pypeitpar.PypeItPar()
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 0.1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['scienceframe']['exprng'] = [1, None]
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
        meta['target'] = dict(ext=0, card='OBJCAT')
        meta['decker'] = dict(ext=0, card='SLMSKNAM')
        meta['binning'] = dict(ext=0, card=None, default='1,1')

        meta['mjd'] = dict(ext=0, card=None, compound=True)
        meta['exptime'] = dict(ext=0, card='EXPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        meta['dispname'] = dict(ext=0, card='GRM_ID')
        #meta['dispangle'] = dict(card=None, compound=True, rtol=1e-5)
        meta['idname'] = dict(ext=0, card='IMAGETYP')
        # Lamps
        meta['lampstat01'] = dict(ext=0, card='LMP_ID')

        # Ingest
        self.meta = meta

    def compound_meta(self, headarr, meta_key):
        if meta_key == 'mjd':
            time = headarr[0]['DATE-OBS']
            ttime = Time(time, format='isot')
            return ttime.mjd
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
        return ['dispname', 'decker']


    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['science', 'standard']:
            return good_exp & (fitstbl['idname'] == 'OBJECT') & (fitstbl['lamps'] == 'Parking') \
                        & (fitstbl['dispname'] != 'OPEN')
        if ftype == 'bias':
            return good_exp & (fitstbl['dispname'] == 'OPEN')
        if ftype in ['pixelflat', 'trace']:
            return good_exp & (fitstbl['idname'] == 'CALIB') & (fitstbl['lamps'] == 'Halogen') \
                        & (fitstbl['dispname'] != 'OPEN')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'arc') & (fitstbl['lamps'] == 'Ne+Hg') \
                        & (fitstbl['dispname'] != 'OPEN')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

