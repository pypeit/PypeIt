"""
Module for Shane/Kast specific methods.

.. include:: ../include/links.rst
"""
import os

from IPython import embed

import numpy as np

from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container
from pypeit import data


class AATUHRFSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle AAT/UHRF specific code
    """
    ndet = 1
    telescope = telescopes.AATTelescopePar()
    url = 'https://aat.anu.edu.au/science/instruments/decomissioned/uhrf/overview'
    ql_supported = False
    name = 'aat_uhrf'
    camera = 'UHRF'
    supported = True
    header_name = 'uhrf'
    allowed_extensions = [".FTS"]

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
        # Retrieve the binning
        binning = '1,1' if hdu is None else self.compound_meta(self.get_headarr(hdu), "binning")
        dsec = 1 + 1024//int(binning.split(',')[0])
        # Detector 1
        detector_dict = dict(
            binning=binning,
            det=1,
            dataext=0,
            specaxis=0,
            specflip=False,
            spatflip=False,
            platescale=0.05,  # Not sure about this value
            saturation=65535.,
            mincounts=-1e10,
            nonlinear=0.76,
            numamplifiers=1,
            gain=np.asarray([1.0]),  # Not sure about this value
            ronoise=np.asarray([0.0]),  # Determine the read noise from the overscan region
            xgap=0.,
            ygap=0.,
            ysize=1.,
            darkcurr=0.0,  # e-/pixel/hour
            # These are rows, columns on the raw frame, 1-indexed
            datasec=np.asarray(['[:, 1:{:d}]'.format(dsec)]),
            oscansec=np.asarray(['[:, {:d}:]'.format(dsec+1)])
        )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        # Bound the detector with slit edges if no edges are found
        par['calibrations']['slitedges']['bound_detector'] = True

        # Never correct for flexure - the sky is subdominant compared to the object and basically never detected.
        par['flexure']['spec_method'] = 'skip'

        # Sky subtraction parameters - this instrument has no sky lines, but we still use the sky subtraction
        # routine to subtract scattered light.
        par['reduce']['skysub']['no_poly'] = True
        par['reduce']['skysub']['bspline_spacing'] = 3.0
        par['reduce']['skysub']['user_regions'] = ':10,75:'  # This is about right for most setups tested so far
        par['scienceframe']['process']['sigclip'] = 10.0

        # Set some parameters for the calibrations
        # par['calibrations']['wavelengths']['reid_arxiv'] = 'None'
        par['calibrations']['wavelengths']['lamps'] = ['ThAr']
        par['calibrations']['wavelengths']['n_final'] = 3
        par['calibrations']['tilts']['spat_order'] = 4
        par['calibrations']['tilts']['spec_order'] = 1

        # Set the default exposure time ranges for the frame typing
        # Trace frames should be the same as arc frames - it will force a bound detector and this
        # allows the scattered light to be subtracted. A pixel-to-pixel sensitivity correction is
        # not needed for this instrument, since it's a small slicer that projects the target onto
        # multiple pixels. This instrument observes bright objects only, so sky subtraction is not
        # important, but the sky subtraction routine is used to subtract scattered light, instead.
        par['calibrations']['arcframe']['exprng'] = [None, 60.0]
        par['calibrations']['tiltframe']['exprng'] = [None, 60.0]
        par['calibrations']['traceframe']['exprng'] = [None, 60.0]
        par['scienceframe']['exprng'] = [61, None]

        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='MEANRA')
        self.meta['dec'] = dict(ext=0, card='MEANDEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        # dispname is arm specific (blue/red)
        self.meta['decker'] = dict(ext=0, card='WINDOW')
        self.meta['dispname'] = dict(ext=0, card='WINDOW')
        self.meta['binning'] = dict(ext=0, card=None, compound=True)
        self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='TOTALEXP')
        self.meta['airmass'] = dict(ext=0, card=None, compound=True)
        # Additional ones, generally for configuration determination or time
        # self.meta['dichroic'] = dict(ext=0, card='BSPLIT_N')
        # self.meta['instrument'] = dict(ext=0, card='VERSION')

    def compound_meta(self, headarr, meta_key):
        """
        Methods to generate metadata requiring interpretation of the header
        data, instead of simply reading the value of a header card.

        Args:
            headarr (:obj:`list`):
                List of `astropy.io.fits.Header`_ objects.
            meta_key (:obj:`str`):
                Metadata keyword to construct.

        Returns:
            object: Metadata value read from the header(s).
        """
        if meta_key == 'mjd':
            date = headarr[0]['UTDATE'].replace(":","-")
            time = headarr[0]['UTSTART']
            ttime = Time(f'{date}T{time}', format='isot')
            return ttime.mjd
        elif meta_key == 'binning':
            binspat = int(np.ceil(1024/headarr[0]['NAXIS1']))
            binspec = int(np.ceil(1024/headarr[0]['NAXIS2']))
            return f'{binspat},{binspec}'
        elif meta_key == 'airmass':
            # Calculate the zenith distance
            zendist = 0.5*(headarr[0]['ZDSTART']+headarr[0]['ZDEND'])
            # Return the airmass based on the zenith distance
            return 1./np.cos(np.deg2rad(zendist))
        msgs.error("Not ready for this compound meta")

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        This list is used by :class:`~pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            :obj:`list`: List of keywords of data pulled from file headers
            and used to constuct the :class:`~pypeit.metadata.PypeItMetaData`
            object.
        """
        # decker is not included because arcs are often taken with a 0.5" slit
        return ['dispname']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.

        Args:
            ftype (:obj:`str`):
                Type of frame to check. Must be a valid frame type; see
                frame-type :ref:`frame_type_defs`.
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for one or more frames to check.
            exprng (:obj:`list`, optional):
                Range in the allowed exposure time for a frame of type
                ``ftype``. See
                :func:`pypeit.core.framematch.check_frame_exptime`.

        Returns:
            `numpy.ndarray`_: Boolean array with the flags selecting the
            exposures in ``fitstbl`` that are ``ftype`` type frames.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['science']:
            return good_exp
        if ftype in ['standard']:
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype == 'bias':
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            # Flats and trace frames are typed together
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            scifile (:obj:`str`):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        par = super().config_specific_par(scifile, inp_par=inp_par)

        if par['calibrations']['wavelengths']['reid_arxiv'] is None:
            msgs.warn("Wavelength setup not supported!" + msgs.newline() + msgs.newline() +
                       "Please perform your own wavelength calibration, and provide the path+filename using:" + msgs.newline() +
                       msgs.pypeitpar_text(['calibrations', 'wavelengths', 'reid_arxiv = <insert path+fileanme>']))
        # Return
        return par
