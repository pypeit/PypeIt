"""
Module for MMT/Blue Channel specific methods.

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
    Child to handle MMT/Blue Channel specific code
    """
    ndet = 1
    name = 'jwst_nirspec'
    header_name = 'jwst_nirspec'
    telescope = telescopes.JWSTTelescopePar()
    camera = 'NIRSPEC'
    supported = True

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        .. warning::

            Many of the necessary detector parameters are read from the file
            header, meaning the ``hdu`` argument is effectively **required** for
            MMT/BlueChannel.  The optional use of ``hdu`` is only viable for
            automatically generated documentation.

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
        detector = detector_container.DetectorContainer(**detector_dicts[det-1])
        return detector


    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
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
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.5
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['fwhm'] = 5.0
        # HeNeAr is by far most commonly used, though ThAr is used for some situations.
        par['calibrations']['wavelengths']['lamps'] = ['ArI', 'ArII', 'HeI', 'NeI']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'

        # Processing steps
        turn_off = dict(use_biasimage=False, use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0
        ## Do not perform global sky subtraction for standard stars
        par['reduce']['skysub']['global_sky_std'] = False

        # cosmic ray rejection parameters for science frames
        par['scienceframe']['process']['sigclip'] = 5.0
        par['scienceframe']['process']['objlim'] = 2.0

        # Set the default exposure time ranges for the frame typing

        # Appropriate exposure times for Blue Channel can vary a lot depending
        # on grating and wavelength. E.g. 300 and 500 line gratings need very
        # short exposures for flats to avoid saturation, but the 1200 and 832
        # can use much longer exposures due to the higher resolution and the
        # continuum lamp not being very bright in the blue/near-UV.
        par['calibrations']['pixelflatframe']['exprng'] = [None, 100]
        par['calibrations']['traceframe']['exprng'] = [None, 100]
        par['calibrations']['standardframe']['exprng'] = [None, 600]
        par['calibrations']['arcframe']['exprng'] = [10, None]
        par['calibrations']['darkframe']['exprng'] = [300, None]

        # less than 30 sec implies conditions are bright enough for scattered
        # light to be significant which affects the illumination of the slit.
        par['calibrations']['illumflatframe']['exprng'] = [30, None]

        # Need to specify this for long-slit data
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['bound_detector'] = True

        # Sensitivity function parameters
        par['sensfunc']['polyorder'] = 7

        return par




