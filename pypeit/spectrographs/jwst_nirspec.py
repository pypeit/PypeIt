"""
Module for JWST NIRSpec specific methods.

.. include:: ../include/links.rst
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
from IPython import embed

class JWSTNIRSpecSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle JWST NIRSpec specific code
    """
    ndet = 2
    name = 'jwst_nirspec'
    header_name = 'jwst_nirspec'
    telescope = telescopes.JWSTTelescopePar()
    camera = 'NIRSPEC'
    url = 'https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph'
    supported = True

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """

        # Detector 1, i.e. NRS1 from
        # https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-instrumentation/nirspec-detectors/nirspec-detector-performance
        detector_dict1 = dict(
            binning='1,1',
            det=1,
            dataext=0,
            specaxis=1,
            specflip=False,
            spatflip=False,
            xgap=0.,
            ygap=0.,
            ysize=1.,
            platescale=0.1,
            darkcurr=0.0092,
            saturation=55100.,
            nonlinear=0.95,  # need to look up and update
            mincounts=-1e10,
            numamplifiers=1,
            gain=np.atleast_1d(0.996),
            ronoise=np.atleast_1d(5.17),
            datasec=None,
            oscansec=None
        )

        # Detector 2
        detector_dict2 = detector_dict1.copy()
        detector_dict2.update(dict(
            det=2,
            dataext=0,
            darkcurr=0.0057,
            saturation=60400.,
            gain=np.atleast_1d(1.137),
            ronoise=np.atleast_1d(6.60),
        ))
        detector_dicts = [detector_dict1, detector_dict2]
        return detector_container.DetectorContainer(**detector_dicts[det-1])

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(ext=0, card='APERTURE')
        self.meta['dichroic'] = dict(ext=0, card='INSFILTE')
        self.meta['binning'] = dict(ext=0, card=None, compound=True)
        self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')

        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='DISPERSE')
        self.meta['idname'] = dict(ext=0, card='IMAGETYP')

        # used for arc and continuum lamps
        self.meta['lampstat01'] = dict(ext=0, card=None, compound=True)
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')


    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()


        # Reduce
        par['reduce']['trim_edge'] = [0,0]

        # Object finding
        par['reduce']['findobj']['find_trim_edge'] = [0,0]
        par['reduce']['findobj']['maxnumber_sci'] = 2
        par['reduce']['findobj']['snr_thresh'] = 10.0
        par['reduce']['findobj']['trace_npoly'] = 5
        par['reduce']['findobj']['snr_thresh'] = 10.0
        par['reduce']['findobj']['find_fwhm'] = 2.0


        # Sky-subtraction
        par['reduce']['skysub']['bspline_spacing'] = 5.0 # JWST sky is smooth
        par['reduce']['skysub']['max_mask_frac'] = 0.95
        par['reduce']['skysub']['mask_by_boxcar'] = True
        par['reduce']['skysub']['sky_sigrej'] = 4.0

        # Extraction
        par['reduce']['extraction']['model_full_slit'] = True
        par['reduce']['extraction']['sn_gauss'] = 6.0
        par['reduce']['extraction']['boxcar_radius'] = 0.2 # extent in calwebb is 0.55" source and on NIRSpec website
        par['reduce']['extraction']['use_2dmodel_mask'] = False # Don't use 2d mask in local skysub

        # Cosmic ray rejection parameters for science frames
        par['scienceframe']['process']['sigclip'] = 5.0
        par['scienceframe']['process']['objlim'] = 2.0
        par['scienceframe']['process']['mask_cr'] = False # Turn off for now since we coadd.

        # Skip reference frame correction for now.
        par['calibrations']['wavelengths']['refframe'] = 'observed'



        return par




