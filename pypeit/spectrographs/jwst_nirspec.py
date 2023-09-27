"""
Module for JWST NIRSpec specific methods.

.. include:: ../include/links.rst
"""
import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit import utils
from pypeit.core import framematch
from pypeit import io
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
            dataext=1,
            specaxis=1,
            specflip=False,
            spatflip=False,
            xgap=180.,
            ygap=0.,
            ysize=1.,
            platescale=0.1,
            darkcurr=33.12,  # e-/pixel/hour  (=0.0092 e-/pixel/s)
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
            dataext=1,
            darkcurr=20.52,  # e-/pixel/hour,  (=0.0057 e-/pixel/s)
            saturation=60400.,
            gain=np.atleast_1d(1.137),
            ronoise=np.atleast_1d(6.60),
        ))
        detector_dicts = [detector_dict1, detector_dict2]
        return detector_container.DetectorContainer(**detector_dicts[det-1])

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
        par['reduce']['extraction']['sn_gauss'] = 5.0
        par['reduce']['extraction']['boxcar_radius'] = 0.2 # extent in calwebb is 0.55" source and on NIRSpec website
        par['reduce']['extraction']['use_2dmodel_mask'] = False # Don't use 2d mask in local skysub

        # Cosmic ray rejection parameters for science frames
        par['scienceframe']['process']['sigclip'] = 5.0
        par['scienceframe']['process']['objlim'] = 2.0
        par['scienceframe']['process']['mask_cr'] = False # Turn off for now since we coadd.

        # Skip reference frame correction for now.
        par['calibrations']['wavelengths']['refframe'] = 'observed'

        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='TARG_RA')
        self.meta['dec'] = dict(ext=0, card='TARG_DEC')
        self.meta['target'] = dict(ext=0, card='TARGPROP')
        self.meta['mode'] = dict(ext=0, card='EXP_TYPE')
        self.meta['decker'] = dict(ext=0, card='APERNAME')

        self.meta['binning'] = dict(ext=0, card=None, default='1,1')
        self.meta['mjd'] = dict(ext=0, card='EXPMID')
        self.meta['exptime'] = dict(ext=0, card='EFFEXPTM')
        self.meta['airmass'] = dict(ext=0, card=None, compound=True)

        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='GRATING')
        self.meta['filter1'] = dict(ext=0, card='FILTER')
        self.meta['idname'] = dict(ext=0, card=None, compound=True)
        self.meta['dithpat'] = dict(ext=0, card=None, compound=True)
        self.meta['dithpos'] = dict(ext=0, card='YOFFSET')

        # used for arc and continuum lamps
        self.meta['lampstat01'] = dict(ext=0, card=None, compound=True)
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')



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

        if meta_key == 'dithpat':
            exp_type = headarr[0].get('EXP_TYPE')
            if exp_type == 'NRS_MSASPEC':
                return headarr[0].get('NOD_TYPE')
            elif exp_type == 'NRS_FIXEDSLIT':
                return headarr[0].get('PATTTYPE')


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
        return ['dispname', 'filter1', 'decker']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard ``PypeIt`` file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        pypeit_keys = super().pypeit_file_keys()
        pypeit_keys.remove('airmass')
        pypeit_keys.remove('binning')
        return pypeit_keys


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

        if ftype == 'science':
            return np.ones(len(fitstbl), dtype=bool)
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)



    @property
    def allowed_mosaics(self):
        """
        Return the list of allowed detector mosaics.

        JWST/NIRSpec only allows for mosaicing the NRS1 and NRS2 detectors.

        Returns:
            :obj:`list`: List of tuples, where each tuple provides the 1-indexed
            detector numbers that can be combined into a mosaic and processed by
            ``PypeIt``.
        """
        return [(1,2)]

    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces
        that are key for image processing.

        Based on readmhdufits.pro

        Parameters
        ----------
        raw_file : :obj:`str`
            File to read
        det : :obj:`int`
            1-indexed detector to read

        Returns
        -------
        detector_par : :class:`pypeit.images.detector_container.DetectorContainer`
            Detector metadata parameters.
        raw_img : `numpy.ndarray`_
            Raw image for this detector.
        hdu : `astropy.io.fits.HDUList`_
            Opened fits file
        exptime : :obj:`float`
            Exposure time read from the file header
        rawdatasec_img : `numpy.ndarray`_
            Data (Science) section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.
        oscansec_img : `numpy.ndarray`_
            Overscan section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.
        """
        fil = utils.find_single_file(f'{raw_file}*', required=True)

        # Read
        msgs.info(f'Reading JWST/NIRSpec file: {fil}')
        hdu = io.fits_open(fil)
        head0 = hdu[0].header

        detector = self.get_detector_par(det if det is not None else 1, hdu=hdu)
        raw_img = hdu[detector['dataext']].data.astype(float)

        # Need the exposure time
        exptime = hdu[self.meta['exptime']['ext']].header[self.meta['exptime']['card']]

        # Return
        return detector, raw_img.T, hdu, exptime, None, None
