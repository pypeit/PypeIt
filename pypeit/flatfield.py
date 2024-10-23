"""
Implements the flat-field class.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
from pathlib import Path
from copy import deepcopy
import inspect
import numpy as np

from scipy import interpolate, ndimage

from astropy.io import fits

from matplotlib import pyplot as plt
from matplotlib import gridspec

from IPython import embed

from pypeit import msgs
from pypeit.pypmsgs import PypeItDataModelError
from pypeit import utils
from pypeit import bspline

from pypeit import datamodel
from pypeit import calibframe
from pypeit import edgetrace
from pypeit import io
from pypeit.display import display
from pypeit.images import buildimage
from pypeit.core import qa
from pypeit.core import flat
from pypeit.core import tracewave
from pypeit.core import basis
from pypeit.core import fitting
from pypeit.core import parse
from pypeit.core.mosaic import build_image_mosaic
from pypeit.spectrographs.util import load_spectrograph
from pypeit import slittrace
from pypeit import dataPaths
from pypeit import cache


class FlatImages(calibframe.CalibFrame):
    """
    Container for the processed flat-field calibrations.

    All of the items in the datamodel are required for instantiation, although
    they can be None (but shouldn't be).

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_flatimages.rst

    """
    version = '1.1.2'

    # Calibration frame attributes
    calib_type = 'Flat'
    calib_file_format = 'fits'

    # Datamodel already includes PYP_SPEC, so no need to combine it with the
    # CalibFrame base datamodel.
    datamodel = {'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 'pixelflat_raw': dict(otype=np.ndarray, atype=np.floating,
                                       descr='Processed, combined pixel flats'),
                 'pixelflat_norm': dict(otype=np.ndarray, atype=np.floating,
                                        descr='Normalized pixel flat'),
                 'pixelflat_model': dict(otype=np.ndarray, atype=np.floating, descr='Model flat'),
                 'pixelflat_spat_bsplines': dict(otype=np.ndarray, atype=bspline.bspline,
                                                 descr='B-spline models for pixel flat; see '
                                                       ':class:`~pypeit.bspline.bspline.bspline`'),
                 'pixelflat_finecorr': dict(otype=np.ndarray, atype=fitting.PypeItFit,
                                       descr='PypeIt 2D polynomial fits to the fine correction of '
                                             'the spatial illumination profile'),
                 'pixelflat_bpm': dict(otype=np.ndarray, atype=np.integer,
                                       descr='Mirrors SlitTraceSet mask for flat-specific flags'),
                 'pixelflat_spec_illum': dict(otype=np.ndarray, atype=np.floating,
                                              descr='Relative spectral illumination'),
                 'pixelflat_waveimg': dict(otype=np.ndarray, atype=np.floating,
                                           descr='Waveimage for pixel flat'),
                 'illumflat_raw': dict(otype=np.ndarray, atype=np.floating,
                                       descr='Processed, combined illum flats'),
                 'illumflat_spat_bsplines': dict(otype=np.ndarray, atype=bspline.bspline,
                                                 descr='B-spline models for illum flat; see '
                                                       ':class:`~pypeit.bspline.bspline.bspline`'),
                 'illumflat_finecorr': dict(otype=np.ndarray, atype=fitting.PypeItFit,
                                       descr='PypeIt 2D polynomial fits to the fine correction of '
                                             'the spatial illumination profile'),
                 'illumflat_bpm': dict(otype=np.ndarray, atype=np.integer,
                                       descr='Mirrors SlitTraceSet mask for flat-specific flags'),
                 'spat_id': dict(otype=np.ndarray, atype=np.integer, descr='Slit spat_id')}

    def __init__(self, pixelflat_raw=None, pixelflat_norm=None, pixelflat_bpm=None,
                 pixelflat_model=None, pixelflat_spat_bsplines=None, pixelflat_finecorr=None,
                 pixelflat_spec_illum=None, pixelflat_waveimg=None, illumflat_raw=None,
                 illumflat_spat_bsplines=None, illumflat_bpm=None, illumflat_finecorr=None,
                 PYP_SPEC=None, spat_id=None):
        # Parse
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=d)

    def _validate(self):
        """
        Validate the instantiation of the flat-field calibrations.
        """
        if self.pixelflat_spat_bsplines is not None and len(self.pixelflat_spat_bsplines) > 0:
            if len(self.spat_id) != len(self.pixelflat_spat_bsplines):
                msgs.error("Pixelflat Bsplines are out of sync with the slit IDs")
        if self.illumflat_spat_bsplines is not None and len(self.illumflat_spat_bsplines) > 0:
            if len(self.spat_id) != len(self.illumflat_spat_bsplines):
                msgs.error("Illumflat Bsplines are out of sync with the slit IDs")

    def is_synced(self, slits):
        """
        Confirm the slits in WaveTilts are aligned to that in SlitTraceSet

        Barfs if not

        Args:
            slits (:class:`~pypeit.slittrace.SlitTraceSet`):

        """
        if not np.array_equal(self.spat_id, slits.spat_id):
            msgs.error('Your flat solutions are out of sync with your slits.  Remove Calibrations'
                       'and restart from scratch.')

    def _bundle(self):
        """
        Over-write default _bundle() method to write one
        HDU per image.  Any extras are in the HDU header of
        the primary image.

        Returns:
            :obj:`list`: A list of dictionaries, each list element is
            written to its own fits extension. See the description
            above.
        """

        # NOTE: Everything in the datamodel is currently (21 Mar 2023) a numpy
        # array, *except* PYP_SPEC.  We need to avoid adding an element to `d`
        # below that *only* contains the PYP_SPEC because it leads to an unnamed
        # empty extension where the only added value is that PYP_SPEC is in the
        # header.  I deal with this by skipping PYP_SPEC in the list of keys and
        # adding it to the dictionaries of the simple objects; i.e., it's not
        # added to entries in `d` that are themselves DataContainer objects
        # (because that causes havoc).

        d = []
        for key in self.keys():
            # Skip None
            if self[key] is None or key == 'PYP_SPEC':
                continue
            if self.datamodel[key]['otype'] == np.ndarray and 'bsplines' not in key \
                    and 'finecorr' not in key:
                d += [{key: {'PYP_SPEC':self.PYP_SPEC, key : self[key]}}]
            elif 'bsplines' in key:
                flattype = 'pixelflat' if 'pixelflat' in key else 'illumflat'
                d += [{f'{flattype}_spat_id-{self.spat_id[ss]}_bspline': self[key][ss]}
                      for ss in range(len(self[key]))]
            elif 'finecorr' in key:
                flattype = 'pixelflat' if 'pixelflat' in key else 'illumflat'
                d += [{f'{flattype}_spat_id-{self.spat_id[ss]}_finecorr': self[key][ss]}
                      for ss in range(len(self[key]))]
            else:
                if len(d) > 0:
                    d[0][key] = self[key]
                else:
                    d += [{key: self[key]}]
        return d

    # NOTE: Previously we had code to override the default to_hdu function,
    # forcing the HDUs to be binary tables.  I don't know why we had done that,
    # but it's not necessary and it was causing havoc with me trying to avoid
    # the empty PYP_SPEC extension.

    @classmethod
    def _parse(cls, hdu, ext=None, transpose_table_arrays=False, hdu_prefix=None, **kwargs):
        """
        Override base-class function to deal with the many idiosyncracies of the
        datamodel.

        See :func:`~pypeit.datamodel.DataContainer._parse` for the argument
        descriptions, and returned items.
        """
        # Grab everything but the bsplines. The bsplines are not parsed
        # because the tailored extension names do not match any of the
        # datamodel keys.
        d, version_passed, type_passed, parsed_hdus = super()._parse(hdu)

        # Find bsplines, if they exist
        nspat = len(d['spat_id'])
        hdunames = [h.name for h in hdu]
        for flattype in ['pixelflat', 'illumflat']:
            # Parse the bspline hdus
            ext_bspl = ['{0}_SPAT_ID-{1}_BSPLINE'.format(flattype.upper(), d['spat_id'][i])
                        for i in range(nspat)]
            indx = np.isin(ext_bspl, hdunames)
            if np.any(indx) and not np.all(indx):
                msgs.error('Expected {0} {1} bspline extensions, but only found {2}.'.format(
                           nspat, flattype, np.sum(indx)))
            if np.all(indx):
                key = '{0}_spat_bsplines'.format(flattype)
                try:
                    d[key] = np.array([bspline.bspline.from_hdu(hdu[k]) for k in ext_bspl])
                except Exception as e:
                    msgs.warn('Error in bspline extension read:\n {0}: {1}'.format(
                                e.__class__.__name__, str(e)))
                    # Assume this is because the type failed
                    type_passed = False
                else:
                    version_passed &= np.all([d[key][i].version == bspline.bspline.version 
                                              for i in range(nspat)])
                    parsed_hdus += ext_bspl
            # Parse the finecorr fits
            ext_fcor = ['{0}_SPAT_ID-{1}_FINECORR'.format(flattype.upper(), d['spat_id'][i])
                        for i in range(nspat)]
            indx = np.isin(ext_fcor, hdunames)
            if np.any(indx) and not np.all(indx):
                msgs.error('Expected {0} {1} finecorr extensions, but only found {2}.'.format(
                           nspat, flattype, np.sum(indx)))
            if np.all(indx):
                key = '{0}_finecorr'.format(flattype)
                try:
                    allfit = []
                    for k in ext_fcor:
                        if hdu[k].data.size == 0:
                            allfit.append(fitting.PypeItFit(None))
                        else:
                            allfit.append(fitting.PypeItFit.from_hdu(hdu[k]))
                    d[key] = np.array(allfit)
                except Exception as e:
                    msgs.warn('Error in finecorr extension read:\n {0}: {1}'.format(
                                e.__class__.__name__, str(e)))
                    # Assume this is because the type failed
                    type_passed = False
                else:
                    version_passed &= np.all([d[key][i].version == fitting.PypeItFit.version
                                              for i in range(nspat)])
                    parsed_hdus += ext_fcor
        return d, version_passed, type_passed, parsed_hdus

    @property
    def shape(self):
        """
        Shape of the image arrays.
        """
        if self.pixelflat_raw is not None:
            return self.pixelflat_raw.shape
        if self.illumflat_raw is not None:
            return self.illumflat_raw.shape
        msgs.error("Shape of FlatImages could not be determined")

    def get_procflat(self, frametype='pixel'):
        """
        Get the processed flat data.

        Args:
            frametype (:obj:`str`, optional):
                The type of flat to return.  Must be either 'illum' for the
                illumination flat or 'pixel' for the pixel flat.

        Returns:
            `numpy.ndarray`_: The selected flat.  Can be None if the flat has
            not been instantiated/processed.
        """
        return self.illumflat_raw if frametype == 'illum' else self.pixelflat_raw

    def get_bpmflats(self, frametype='pixel'):
        """
        Get the processed bad-pixel mask.

        Args:
            frametype (:obj:`str`, optional):
                The type of mask to return.  Must be either 'illum' for the
                illumination flat mask or 'pixel' for the pixel flat mask.

        Returns:
            `numpy.ndarray`_: The selected mask.  If neither the illumination
            flat or pixel flat mask exist, the returned array is fully unmasked
            (all values are False).
        """
        # Check if both BPMs are none
        if self.pixelflat_bpm is None and self.illumflat_bpm is None:
            msgs.warn("FlatImages contains no BPM - trying to generate one")
            return np.zeros(self.shape, dtype=int)
        # Now return the requested case, checking for None
        if frametype == 'illum':
            if self.illumflat_bpm is not None:
                return self.illumflat_bpm
            msgs.warn("illumflat has no BPM - using the pixelflat BPM")
            return self.pixelflat_bpm
        if self.pixelflat_bpm is not None:
            return self.pixelflat_bpm
        msgs.warn("pixelflat has no BPM - using the illumflat BPM")
        return self.illumflat_bpm

    def get_spat_bsplines(self, frametype='illum', finecorr=False):
        """
        Grab a list of bspline fits

        Args:
            frametype (:obj:`str`, optional):
                The type of mask to return.  Must be either 'illum' for the
                illumination flat mask or 'pixel' for the pixel flat mask.
            finecorr (:obj:`bool`, optional):
                If True, return the fine correction bsplines; otherwise, return
                the zeroth order correction.

        Returns:
            :obj:`list`: The selected list of spatial bsplines.  Can be None if
            the requested data (or the fall-back) do not exist.
        """
        # Decide if the finecorrection splines are needed, or the zeroth order correction
        if finecorr:
            fctxt = 'fine correction to the '
            pixel_bsplines = self.pixelflat_finecorr
            illum_bsplines = self.illumflat_finecorr
            # Do a quick check if no data exist
            if pixel_bsplines is not None and pixel_bsplines[0].xval is None:
                pixel_bsplines = None
            if illum_bsplines is not None and illum_bsplines[0].xval is None:
                illum_bsplines = None
        else:
            fctxt = ''
            pixel_bsplines = self.pixelflat_spat_bsplines
            illum_bsplines = self.illumflat_spat_bsplines
        # Ensure that at least one has been generated
        if pixel_bsplines is None and illum_bsplines is None:
            msgs.warn(f'FlatImages contains no {fctxt}spatial bspline fit.')
            return None
        # Now return the requested case, checking for None
        if frametype == 'illum':
            if illum_bsplines is not None:
                return illum_bsplines
            msgs.warn(f'illumflat has no {fctxt}spatial bspline fit - using the pixelflat.')
            return pixel_bsplines
        if pixel_bsplines is not None:
            return pixel_bsplines
        msgs.warn(f'pixelflat has no {fctxt}spatial bspline fit - using the illumflat.')
        return illum_bsplines

    def fit2illumflat(self, slits, frametype='illum', finecorr=False, initial=False,
                      spat_flexure=None):
        """
        Construct the model flat using the spatial bsplines.

        Args:
            slits (:class:`~pypeit.slittrace.SlitTraceSet`):
                Definition of the slit edges
            frametype (str, optional):
                The frame type should be 'illum' to return the illumflat
                version, or 'pixel' to return the pixelflat version.
            finecorr (bool, optional):
                Return the fine correction bsplines (finecorr=True), or the
                zeroth order correction (finecorr=False)
            initial (bool, optional):
                If True, the initial slit edges will be used
            spat_flexure (float, optional):
                Spatial flexure in pixels

        Returns:
            `numpy.ndarray`_: An image of the spatial illumination profile for all slits.
        """
        # Check spatial flexure type
        if spat_flexure is not None and not isinstance(spat_flexure, float):
            msgs.error('Spatial flexure must be None or float.')
        # Initialise the returned array
        illumflat = np.ones(self.shape, dtype=float)
        # Load spatial bsplines
        spat_bsplines = self.get_spat_bsplines(frametype=frametype, finecorr=finecorr)
        # Check that the bsplines exist
        if spat_bsplines is None:
            if finecorr:
                return np.ones(self.shape, dtype=float)
            msgs.error('Cannot continue without spatial bsplines.')

        # Loop
        for slit_idx in range(slits.nslits):
            # Skip masked
            if slits.mask[slit_idx] != 0:
                continue
            # Skip those without a bspline
            # DO it
            _slitid_img = slits.slit_img(slitidx=slit_idx, initial=initial, flexure=spat_flexure)
            onslit = _slitid_img == slits.spat_id[slit_idx]
            spat_coo = slits.spatial_coordinate_image(slitidx=slit_idx,
                                                      initial=initial,
                                                      slitid_img=_slitid_img,
                                                      flexure_shift=spat_flexure)
            if finecorr:
                spec_coo = np.where(onslit)[0] / (slits.nspec - 1)
                illumflat[onslit] = spat_bsplines[slit_idx].eval(spat_coo[onslit], spec_coo)
            else:
                illumflat[onslit] = spat_bsplines[slit_idx].value(spat_coo[onslit])[0]
        # TODO -- Update the internal one?  Or remove it altogether??
        return illumflat

    def show(self, frametype='all', slits=None, wcs_match=True, chk_version=True):
        """
        Simple wrapper to :func:`show_flats`.

        Args:
            frametype (str, optional):
                String used to select the flats to be displayed.  The frame type
                should be 'illum' to show the illumflat version, 'pixel' to show
                the pixelflat version, or 'all' to show both.
            slits (:class:`~pypeit.slittrace.SlitTraceSet`):
                Definition of the slit edges
            wcs_match (:obj:`bool`, optional):
                (Attempt to) Match the WCS coordinates of the output images in
                the `ginga`_ viewer.
            chk_version (:obj:`bool`, optional):
                When reading in existing files written by PypeIt, perform strict
                version checking to ensure a valid file.  If False, the code
                will try to keep going, but this may lead to faults and quiet
                failures.  User beware!
        """
        illumflat_pixel, illumflat_illum = None, None
        pixelflat_finecorr, illumflat_finecorr = None, None
        pixelflat_totalillum, illumflat_totalillum = None, None

        if slits is None and self.calib_dir is not None or self.calib_key is not None:
            # If the slits are not defined, and the relevant attributes are set,
            # try to read the associated SlitTraceSet
            slits_file = slittrace.SlitTraceSet.construct_file_name(self.calib_key,
                                                                    calib_dir=self.calib_dir)
            try:
                slits = slittrace.SlitTraceSet.from_file(slits_file, chk_version=chk_version)
            except (FileNotFoundError, PypeItDataModelError):
                msgs.warn('Could not load slits to include when showing flat-field images.  File '
                          'was either not provided directly, or it could not be read based on its '
                          f'expected name: {slits_file}.')

        if slits is not None:
            slits.mask_flats(self)
            illumflat_pixel = self.fit2illumflat(slits, frametype='pixel', finecorr=False)
            pixelflat_finecorr = self.fit2illumflat(slits, frametype='pixel', finecorr=True)
            if self.illumflat_spat_bsplines is not None:
                illumflat_illum = self.fit2illumflat(slits, frametype='illum', finecorr=False)
            if self.illumflat_finecorr is not None:
                illumflat_finecorr = self.fit2illumflat(slits, frametype='illum', finecorr=True)

        # Construct a total illumination flat if the fine correction has been computed
        if pixelflat_finecorr is not None:
            pixelflat_totalillum = illumflat_pixel*pixelflat_finecorr
        if illumflat_finecorr is not None:
            illumflat_totalillum = illumflat_illum*illumflat_finecorr

        # Decide which frames should be displayed
        if frametype == 'pixel':
            image_list = zip([self.pixelflat_norm, illumflat_pixel, pixelflat_finecorr, pixelflat_totalillum,
                              self.pixelflat_raw, self.pixelflat_model, self.pixelflat_spec_illum],
                             ['pixelflat_norm', 'pixelflat_spat_illum', 'pixelflat_finecorr', 'pixelflat_totalillum',
                              'pixelflat_raw', 'pixelflat_model', 'pixelflat_spec_illum'],
                             [(0.9, 1.1), (0.9, 1.1), (0.95, 1.05), (0.9, 1.1),
                              None, None, (0.8, 1.2)])
        elif frametype == 'illum':
            image_list = zip([illumflat_illum, illumflat_finecorr, illumflat_totalillum, self.illumflat_raw],
                             ['illumflat_spat_illum', 'illumflat_finecorr', 'illumflat_totalillum', 'illumflat_raw'],
                             [(0.9, 1.1), (0.95, 1.05), (0.9, 1.1), None])
        else:
            # Show everything that's available (anything that is None will not be displayed)
            image_list = zip([self.pixelflat_norm, illumflat_pixel, pixelflat_finecorr, pixelflat_totalillum,
                              self.pixelflat_raw, self.pixelflat_model, self.pixelflat_spec_illum,
                              illumflat_illum, illumflat_finecorr, illumflat_totalillum, self.illumflat_raw],
                             ['pixelflat_norm', 'pixelflat_spat_illum', 'pixelflat_finecorr', 'pixelflat_totalillum',
                              'pixelflat_raw', 'pixelflat_model', 'pixelflat_spec_illum',
                              'illumflat_spat_illum', 'illumflat_finecorr', 'illumflat_totalillum', 'illumflat_raw'],
                             [(0.9, 1.1), (0.9, 1.1), (0.95, 1.05), (0.9, 1.1),
                              None, None, (0.8, 1.2),
                              (0.9, 1.1), (0.95, 1.05), (0.9, 1.1), None])
        # Display frames
        show_flats(image_list, wcs_match=wcs_match, slits=slits, waveimg=self.pixelflat_waveimg)


class FlatField:
    """
    Builds pixel-level flat-field and the illumination flat-field.

    For the primary methods, see :func:`run`.

    Args:
        rawflatimg (:class:`~pypeit.images.pypeitimage.PypeItImage`):
            Processed, combined set of pixelflat images
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The `Spectrograph` instance that sets the instrument used to
            take the observations.
        flatpar (:class:`~pypeit.par.pypeitpar.FlatFieldPar`):
            User-level parameters for constructing the flat-field
            corrections.  If None, the default parameters are used.
        slits (:class:`~pypeit.slittrace.SlitTraceSet`):
            The current slit traces.
        wavetilts (:class:`~pypeit.wavetilts.WaveTilts`, optional):
            The current fit to the wavelength tilts. I can be None,
            for example, if slitless is True.
        wv_calib (:class:`~pypeit.wavecalib.WaveCalib`, optional):
            Wavelength calibration object. It can be None, for example, if
            slitless is True.
        slitless (bool, optional):
            True if the input rawflatimg is a slitless flat. Default is False.
        spat_illum_only (bool, optional):
            Only perform the spatial illumination calculation, and ignore
            the 2D bspline fit. This should only be set to true if you
            want the spatial illumination profile only. If you want to
            simultaneously generate a pixel flat and a spatial
            illumination profile from the same input, this should be
            False (which is the default).
        qa_path (str, optional):
            Path to QA directory

    Attributes:
        rawflatimg (:class:`~pypeit.images.pypeitimage.PypeItImage`):
        mspixelflat (`numpy.ndarray`_):
            Normalized flat
        msillumflat (`numpy.ndarray`_):
            Illumination flat
        flat_model (`numpy.ndarray`_):
            Model of the flat
        list_of_spat_bsplines (list):
        spec_illum (`numpy.ndarray`_):
            Image of the relative spectral illumination for a multislit spectrograph

    """
    def __init__(self, rawflatimg, spectrograph, flatpar, slits, wavetilts=None, wv_calib=None,
                 slitless=False, spat_illum_only=False, qa_path=None, calib_key=None):

        # Defaults
        self.spectrograph = spectrograph
        # TODO: We're passing this key only so we can create the output QA
        # filename.  Can we find a different way?
        self.calib_key = calib_key
        self.qa_path = qa_path
        # FieldFlattening parameters
        self.flatpar = flatpar

        # Input data
        self.slits = slits
        self.wavetilts = wavetilts
        self.wv_calib = wv_calib

        # Worth a check
        if self.wavetilts is not None and not slitless:
            self.wavetilts.is_synced(self.slits)

        # Attributes unique to this Object
        self.rawflatimg = rawflatimg      # Un-normalized pixel flat as a PypeItImage
        self.mspixelflat = None     # Normalized pixel flat
        self.msillumflat = None     # Illumination flat
        self.flat_model = None      # Model flat
        self.list_of_spat_bsplines = None
        self.list_of_finecorr_fits = None
        self.spat_illum_only = spat_illum_only
        self.spec_illum = None      # Relative spectral illumination image
        self.waveimg = None
        self.slitless = slitless    # is this a slitless flat?

        # get waveimg here if available
        if self.wavetilts is None or self.wv_calib is None:
            msgs.warn("Wavelength calib or tilts are not available.  Wavelength image not generated.")
        else:
            self.build_waveimg()   # this set self.waveimg

        # Completed steps
        self.steps = []

    @property
    def nslits(self):
        """
        Return the number of slits.  Pulled directly from :attr:`slits`, if it exists.
        """
        return 0 if self.slits is None else self.slits.nslits

    # TODO: Need to add functionality to use a different frame for the
    # ilumination flat, e.g. a sky flat
    def run(self, doqa=False, debug=False, show=False):
        """
        Generate normalized pixel and illumination flats.

        This is a simple wrapper for the main flat-field methods:

            - Full 2D model, illumination flat, and pixel flat images are
              constructed by :func:`fit`.

            - The results can be shown in a ginga window using :func:`show`.

        The method is a simple wrapper for :func:`fit` and :func:`show`.

        Args:
            doqa (:obj:`bool`, optional):
                Save the QA?
            debug (:obj:`bool`, optional):
                Run in debug mode.
            show (:obj:`bool`, optional):
                Show the results in the ginga viewer.

        Returns:
            :class:`FlatImages`: Container with the results of the flat-field
            analysis.
        """

        # check if self.wavetilts is available. It can be None if the flat is slitless, but it's needed otherwise
        if self.wavetilts is None and not self.slitless:
            msgs.warn("Wavelength tilts are not available.  Cannot generate this flat image.")
            return None

        # Fit it
        # NOTE: Tilts do not change and self.slits is updated internally.
        if not self.flatpar['fit_2d_det_response']:
            # This spectrograph does not have a structure correction
            # implemented. Ignore detector structure.
            self.fit(spat_illum_only=self.spat_illum_only, doqa=doqa, debug=debug)
        elif self.waveimg is not None:
            # Iterate on the pixelflat if required by the spectrograph
            # User has requested a structure correction.
            # Note: This will only be performed if it is coded for each individual spectrograph.
            # Make a copy of the original flat
            rawflat_orig = self.rawflatimg.image.copy()
            # TODO: Should this be *any* flag, or just BPM?
            gpm = self.rawflatimg.select_flag(flag='BPM', invert=True)
            # Just get the spatial and spectral profiles for now
            self.fit(spat_illum_only=self.spat_illum_only, doqa=doqa, debug=debug)
            # If we're only doing the spatial illumination profile, the detector structure
            # has already been divided out by the pixel flat. No need to calculate structure
            if not self.spat_illum_only:
                niter = 1  # Just do one iteration... two is too long, and doesn't significantly improve the fine spatial illumination correction.
                det_resp_model = 1  # Initialise detector structure to a value of 1 (i.e. no detector structure)
                onslits = self.slits.slit_img(pad=-self.flatpar['slit_trim'], initial=False) != -1
                for ff in range(niter):
                    # If we're only doing the spatial illumination profile, the detector structure
                    # has already been divided out by the pixel flat.
                    if self.spat_illum_only:
                        break
                    msgs.info("Iteration {0:d}/{1:d} of 2D detector response extraction".format(ff+1, niter))
                    # Extract a detector response image
                    det_resp = self.extract_structure(rawflat_orig)
                    # Trim the slits to avoid edge effects
                    gpmask = (self.waveimg != 0.0) & gpm & onslits
                    # Model the 2D detector response in an instrument specific way
                    det_resp_model = self.spectrograph.fit_2d_det_response(det_resp, gpmask)
                    # Apply this model
                    self.rawflatimg.image = rawflat_orig * utils.inverse(det_resp_model)
                    # Perform a 2D fit with the cleaned image
                    self.fit(spat_illum_only=self.spat_illum_only, doqa=doqa, debug=debug)
                # Save the QA, if requested
                if doqa:
                    # TODO :: Probably need to pass in det when more spectrographs implement a structure correction...
                    outfile = qa.set_qa_filename("DetectorStructure_" + self.calib_key, 'detector_structure',
                                                 det="DET01", out_dir=self.qa_path)
                    detector_structure_qa(det_resp, det_resp_model, outfile=outfile)
                # Include the structure in the flat model and the pixelflat
                self.mspixelflat *= det_resp_model
                # Reset the rawimg
                self.rawflatimg.image = rawflat_orig

        # Show the flatfield images if requested
        if show:
            self.show(wcs_match=True)

        # Build the mask
        bpmflats = self.build_mask()

        # Return
        if self.spat_illum_only:
            # Illumination correction only
            return FlatImages(illumflat_raw=self.rawflatimg.image,
                              illumflat_spat_bsplines=np.asarray(self.list_of_spat_bsplines),
                              illumflat_finecorr=np.asarray(self.list_of_finecorr_fits),
                              illumflat_bpm=bpmflats, PYP_SPEC=self.spectrograph.name,
                              spat_id=self.slits.spat_id)

        # Pixel and illumination correction only
        return FlatImages(pixelflat_raw=self.rawflatimg.image,
                          pixelflat_norm=self.mspixelflat,
                          pixelflat_model=self.flat_model,
                          pixelflat_spat_bsplines=np.asarray(self.list_of_spat_bsplines),
                          pixelflat_finecorr=np.asarray(self.list_of_finecorr_fits),
                          pixelflat_bpm=bpmflats, pixelflat_spec_illum=self.spec_illum,
                          pixelflat_waveimg=self.waveimg,
                          PYP_SPEC=self.spectrograph.name, spat_id=self.slits.spat_id)

    def build_mask(self):
        """
        Generate bad pixel mask.

        Returns:
            `numpy.ndarray`_ : bad pixel mask
        """
        bpmflats = np.zeros_like(self.slits.mask, dtype=self.slits.bitmask.minimum_dtype())
        for flag in ['SKIPFLATCALIB', 'BADFLATCALIB']:
            bpm = self.slits.bitmask.flagged(self.slits.mask, flag)
            if np.any(bpm):
                bpmflats[bpm] = self.slits.bitmask.turn_on(bpmflats[bpm], flag)
        return bpmflats

    def build_waveimg(self):
        """
        Generate an image of the wavelength of each pixel.
        """
        msgs.info("Generating wavelength image")
        if self.wavetilts is None or self.wv_calib is None:
            msgs.error("Wavelength calib or tilts are not available.  Cannot generate wavelength image.")
        else:
            flex = self.wavetilts.spat_flexure
            slitmask = self.slits.slit_img(initial=True, flexure=flex)
            tilts = self.wavetilts.fit2tiltimg(slitmask, flexure=flex)
            # Save to class attribute for inclusion in the Flat calibration frame
            self.waveimg = self.wv_calib.build_waveimg(tilts, self.slits, spat_flexure=flex)

    def show(self, wcs_match=True):
        """
        Show all of the flat field products in ginga.

        Args:
            wcs_match (:obj:`bool`, optional):
                Match the WCS of the flat-field images
        """
        # Prepare the images to show, their names and their cuts
        image_list = zip([self.mspixelflat, self.msillumflat, self.rawflatimg.image, self.flat_model],
                         ['pixelflat', 'spat_illum', 'raw', 'model', 'spec_illum'],
                         [(0.9, 1.1), (0.9, 1.1), None, None, (0.8, 1.2)])
        show_flats(image_list, wcs_match=wcs_match, slits=self.slits, waveimg=self.waveimg)

    def fit(self, spat_illum_only=False, doqa=True, debug=False):
        """
        Construct a model of the flat-field image.

        For this method to work, :attr:`rawflatimg` must have been
        previously constructed.

        The method loops through all slits provided by the :attr:`slits`
        object, except those that have been masked (i.e., slits with
        ``self.slits.mask == True`` are skipped).  For each slit:

            - Collapse the flat-field data spatially using the
              wavelength coordinates provided by the fit to the arc-line
              traces (:class:`~pypeit.wavetilts.WaveTilts`), and fit the
              result with a bspline.  This provides the
              spatially-averaged spectral response of the instrument.
              The data used in the fit is trimmed toward the slit
              spatial center via the ``slit_trim`` parameter in
              :attr:`flatpar`.
            - Use the bspline fit to construct and normalize out the
              spectral response.
            - Collapse the normalized flat-field data spatially using a
              coordinate system defined by the left slit edge.  The data
              included in the spatial (illumination) profile calculation
              is expanded beyond the nominal slit edges using the
              ``slit_illum_pad`` parameter in :attr:`flatpar`.  The raw,
              collapsed data is then median filtered (see ``spat_samp``
              in :attr:`flatpar`) and Gaussian filtered; see
              :func:`pypeit.core.flat.illum_filter`.  This creates an
              empirical, highly smoothed representation of the
              illumination profile that is fit with a bspline using
              the :func:`spatial_fit` method.  The
              construction of the empirical illumination profile (i.e.,
              before the bspline fitting) can be done iteratively, where
              each iteration sigma-clips outliers; see the
              ``illum_iter`` and ``illum_rej`` parameters in
              :attr:`flatpar` and
              :func:`pypeit.core.flat.construct_illum_profile`.
            - If requested, the 1D illumination profile is used to
              "tweak" the slit edges by offsetting them to a threshold
              of the illumination peak to either side of the slit center
              (see ``tweak_slits_thresh`` in :attr:`flatpar`), up to a
              maximum allowed shift from the existing slit edge (see
              ``tweak_slits_maxfrac`` in :attr:`flatpar`).  See
              :func:`~pypeit.core.flat.tweak_slit_edges`.  If tweaked, the
              :func:`spatial_fit` is repeated to place it on the tweaked
              slits reference frame.
            - Use the bspline fit to construct the 2D illumination image
              (:attr:`msillumflat`) and normalize out the spatial
              response.
            - Fit the residuals of the flat-field data that has been
              independently normalized for its spectral and spatial
              response with a 2D bspline-polynomial fit.  The order of
              the polynomial has been optimized via experimentation; it
              can be changed but you should use extreme caution when
              doing so (see ``twod_fit_npoly``).  The multiplication of
              the 2D spectral response, 2D spatial response, and joint
              2D fit to the high-order residuals define the final flat
              model (:attr:`flat_model`).
            - Finally, the pixel-to-pixel response of the instrument is
              defined as the ratio of the raw flat data to the
              best-fitting flat-field model (:attr:`mspixelflat`)

        This method is the primary method that builds the
        :class:`FlatField` instance, constructing :attr:`mspixelflat`,
        :attr:`msillumflat`, and :attr:`flat_model`.  All of these
        attributes are altered internally.  If the slit edges are to be
        tweaked using the 1D illumination profile (``tweak_slits`` in
        :attr:`flatpar`), the tweaked slit edge arrays in the internal
        :class:`~pypeit.slittrace.SlitTraceSet` object, :attr:`slits`,
        are also altered.

        Used parameters from :attr:`flatpar`
        (:class:`~pypeit.par.pypeitpar.FlatFieldPar`) are
        ``spec_samp_fine``, ``spec_samp_coarse``, ``spat_samp``,
        ``tweak_slits``, ``tweak_slits_thresh``,
        ``tweak_slits_maxfrac``, ``rej_sticky``, ``slit_trim``,
        ``slit_illum_pad``, ``illum_iter``, ``illum_rej``, and
        ``twod_fit_npoly``, ``saturated_slits``.

        **Revision History**:

            - 11-Mar-2005  First version written by Scott Burles.
            - 2005-2018    Improved by J. F. Hennawi and J. X. Prochaska
            - 3-Sep-2018 Ported to python by J. F. Hennawi and significantly improved

        Args:
            spat_illum_only (:obj:`bool`, optional):
                If true, only the spatial illumination profile will be calculated.
                The 2D bspline fit will not be performed. This is primarily used
                to build an illumflat.
            doqa (:obj:`bool`, optional):
                Save the QA?
            debug (:obj:`bool`, optional):
                Show plots useful for debugging. This will block
                further execution of the code until the plot windows
                are closed.

        """
        # TODO: break up this function!  Can it be partitioned into a series of "core" methods?
        # TODO: JFH I wrote all this code and will have to maintain it and I don't want to see it broken up.
        # TODO: JXP This definitely needs breaking up..

        # Initialise with a series of bad splines (for when slits go wrong)
        if self.list_of_spat_bsplines is None:
            self.list_of_spat_bsplines = [bspline.bspline(None) for all in self.slits.spat_id]
        if self.list_of_finecorr_fits is None:
            self.list_of_finecorr_fits = [fitting.PypeItFit(None) for all in self.slits.spat_id]

        # Set parameters (for convenience;
        spec_samp_fine = self.flatpar['spec_samp_fine']
        spec_samp_coarse = self.flatpar['spec_samp_coarse']
        tweak_method = self.flatpar['tweak_method']
        tweak_slits = self.flatpar['tweak_slits']
        tweak_slits_thresh = self.flatpar['tweak_slits_thresh']
        tweak_slits_maxfrac = self.flatpar['tweak_slits_maxfrac']
        # If sticky, points rejected at each stage (spec, spat, 2d) are
        # propagated to the next stage
        sticky = self.flatpar['rej_sticky']
        trim = self.flatpar['slit_trim']
        pad = self.flatpar['slit_illum_pad']
        # Iteratively construct the illumination profile by rejecting outliers
        npoly = self.flatpar['twod_fit_npoly']
        saturated_slits = self.flatpar['saturated_slits']

        # Setup images
        nspec, nspat = self.rawflatimg.image.shape
        rawflat = self.rawflatimg.image
        # Good pixel mask
        # TODO: Should this be *any* flag, or just BPM?
        gpm = self.rawflatimg.select_flag(flag='BPM', invert=True)

        # Flat-field modeling is done in the log of the counts
        flat_log = np.log(np.fmax(rawflat, 1.0))
        gpm_log = (rawflat > 1.0) & gpm
        # set errors to just be 0.5 in the log
        ivar_log = gpm_log.astype(float)/0.5**2

        # Get the non-linear count level
        if self.rawflatimg.is_mosaic:
            # if this is a mosaic we take the maximum value among all the detectors
            nonlinear_counts = np.max([rawdets.nonlinear_counts() for rawdets in self.rawflatimg.detector.detectors])
        else:
            nonlinear_counts = self.rawflatimg.detector.nonlinear_counts()

        # TODO -- JFH -- CONFIRM THIS SHOULD BE ON INIT
        # It does need to be *all* of the slits
        median_slit_widths = np.median(self.slits.right_init - self.slits.left_init, axis=0)

        if tweak_slits:
            # NOTE: This copies the input slit edges to a set that can be tweaked.
            self.slits.init_tweaked()

        # TODO: This needs to include a padding check
        # Construct three versions of the slit ID image, all of unmasked slits!
        #   - an image that uses the padding defined by self.slits
        slitid_img_init = self.slits.slit_img(initial=True)
        #   - an image that uses the extra padding defined by
        #     self.flatpar. This was always 5 pixels in the previous
        #     version.
        padded_slitid_img = self.slits.slit_img(initial=True, pad=pad)
        #   - and an image that trims the width of the slit using the
        #     parameter in self.flatpar. This was always 3 pixels in
        #     the previous version.
        # TODO: Fix this for when trim is a tuple
        trimmed_slitid_img = self.slits.slit_img(pad=-trim, initial=True)

        # Prep for results
        self.mspixelflat = np.ones_like(rawflat)
        self.msillumflat = np.ones_like(rawflat)
        self.flat_model = np.zeros_like(rawflat)

        # Allocate work arrays only once
        spec_model = np.ones_like(rawflat)
        norm_spec = np.ones_like(rawflat)
        norm_spec_spat = np.ones_like(rawflat)
        twod_model = np.ones_like(rawflat)
        twod_gpm_out = np.ones_like(rawflat, dtype=bool)

        # #################################################
        # Model each slit independently
        for slit_idx, slit_spat in enumerate(self.slits.spat_id):
            # Is this a good slit??
            if self.slits.bitmask.flagged(self.slits.mask[slit_idx], flag=['SHORTSLIT', 'USERIGNORE', 'BADTILTCALIB']):
                msgs.info('Skipping bad slit: {}'.format(slit_spat))
                self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(self.slits.mask[slit_idx], 'BADFLATCALIB')
                continue
            elif self.slits.bitmask.flagged(self.slits.mask[slit_idx], flag=['BOXSLIT']):
                msgs.info('Skipping alignment slit: {}'.format(slit_spat))
                continue
            elif self.slits.bitmask.flagged(self.slits.mask[slit_idx], flag=['BADWVCALIB']) and \
                    (self.flatpar['pixelflat_min_wave'] is not None or self.flatpar['pixelflat_max_wave'] is not None):
                msgs.info('Skipping slit with bad wavecalib: {}'.format(slit_spat))
                self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(self.slits.mask[slit_idx], 'BADFLATCALIB')
                continue

            msgs.info('Modeling the flat-field response for slit spat_id={}: {}/{}'.format(
                        slit_spat, slit_idx+1, self.slits.nslits))

            # Find the pixels on the initial slit
            onslit_init = slitid_img_init == slit_spat

            # Check for saturation of the flat. If there are not enough
            # pixels do not attempt a fit, and continue to the next
            # slit.
            # TODO: set the threshold to a parameter?
            good_frac = np.sum(onslit_init & (rawflat < nonlinear_counts))/np.sum(onslit_init)
            if good_frac < 0.5:
                common_message = 'To change the behavior, use the \'saturated_slits\' parameter ' \
                                 'in the \'flatfield\' parameter group; see here:\n\n' \
                                 'https://pypeit.readthedocs.io/en/latest/pypeit_par.html \n\n' \
                                 'You could also choose to use a different flat-field image ' \
                                 'for this calibration group.'
                if saturated_slits == 'crash':
                    msgs.error('Only {:4.2f}'.format(100*good_frac)
                               + '% of the pixels on slit {0} are not saturated.  '.format(slit_spat)
                               + 'Selected behavior was to crash if this occurred.  '
                               + common_message)
                elif saturated_slits == 'mask':
                    self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(self.slits.mask[slit_idx], 'BADFLATCALIB')
                    msgs.warn('Only {:4.2f}'.format(100*good_frac)
                                                + '% of the pixels on slit {0} are not saturated.  '.format(slit_spat)
                              + 'Selected behavior was to mask this slit and continue with the '
                              + 'remainder of the reduction, meaning no science data will be '
                              + 'extracted from this slit.  ' + common_message)
                elif saturated_slits == 'continue':
                    self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(self.slits.mask[slit_idx], 'SKIPFLATCALIB')
                    msgs.warn('Only {:4.2f}'.format(100*good_frac)
                              + '% of the pixels on slit {0} are not saturated.  '.format(slit_spat)
                              + 'Selected behavior was to simply continue, meaning no '
                              + 'field-flatting correction will be applied to this slit but '
                              + 'pypeit will attempt to extract any objects found on this slit.  '
                              + common_message)
                else:
                    # Should never get here
                    raise NotImplementedError('Unknown behavior for saturated slits: {0}'.format(
                                              saturated_slits))
                continue

            # Demand at least 10 pixels per row (on average) per degree
            # of the polynomial.
            # NOTE: This is not used until the 2D fit. Defined here to
            # be close to the definition of ``onslit``.
            if npoly is None:
                # Approximate number of pixels sampling each spatial pixel
                # for this (original) slit.
                npercol = np.fmax(np.floor(np.sum(onslit_init)/nspec),1.0)
                npoly  = np.clip(7, 1, int(np.ceil(npercol/10.)))
            
            # TODO: Always calculate the optimized `npoly` and warn the
            #  user if npoly is provided but higher than the nominal
            #  calculation?

            # Create an image with the spatial coordinates relative to the left edge of this slit
            spat_coo_init = self.slits.spatial_coordinate_image(slitidx=slit_idx, full=True, initial=True)

            # Find pixels on the padded and trimmed slit coordinates
            onslit_padded = padded_slitid_img == slit_spat
            onslit_trimmed = trimmed_slitid_img == slit_spat

            # ----------------------------------------------------------
            # Collapse the slit spatially and fit the spectral function
            # TODO: Put this stuff in a self.spectral_fit method?

            # Create the tilts image for this slit
            if self.slitless:
                tilts = np.tile(np.arange(rawflat.shape[0]) / rawflat.shape[0], (rawflat.shape[1], 1)).T
            else:
                # TODO -- JFH Confirm the sign of this shift is correct!
                _flexure = 0. if self.wavetilts.spat_flexure is None else self.wavetilts.spat_flexure
                tilts = tracewave.fit2tilts(rawflat.shape, self.wavetilts['coeffs'][:,:,slit_idx],
                                            self.wavetilts['func2d'], spat_shift=-1*_flexure)
            # Convert the tilt image to an image with the spectral pixel index
            spec_coo = tilts * (nspec-1)

            # Only include the trimmed set of pixels in the flat-field
            # fit along the spectral direction.
            spec_gpm = onslit_trimmed & gpm_log  # & (rawflat < nonlinear_counts)
            spec_nfit = np.sum(spec_gpm)
            spec_ntot = np.sum(onslit_init)
            msgs.info('Spectral fit of flatfield for {0}/{1} '.format(spec_nfit, spec_ntot)
                      + ' pixels in the slit.')
            # Set this to a parameter?
            if spec_nfit/spec_ntot < 0.5:
                # TODO: Shouldn't this raise an exception or continue to the next slit instead?
                msgs.warn('Spectral fit includes only {:.1f}'.format(100*spec_nfit/spec_ntot)
                          + '% of the pixels on this slit.' + msgs.newline()
                          + '          Either the slit has many bad pixels or the number of '
                            'trimmed pixels is too large.')

            # Sort the pixels by their spectral coordinate.
            # TODO: Include ivar and sorted gpm in outputs?
            spec_gpm, spec_srt, spec_coo_data, spec_flat_data \
                    = flat.sorted_flat_data(flat_log, spec_coo, gpm=spec_gpm)
            # NOTE: By default np.argsort sorts the data over the last
            # axis. Just to avoid the possibility (however unlikely) of
            # spec_coo[spec_gpm] returning an array, all the arrays are
            # explicitly flattened.
            spec_ivar_data = ivar_log[spec_gpm].ravel()[spec_srt]
            spec_gpm_data = gpm_log[spec_gpm].ravel()[spec_srt]

            # Rejection threshold for spectral fit in log(image)
            # TODO: Make this a parameter?
            logrej = 0.5

            # Fit the spectral direction of the blaze.
            # TODO: Figure out how to deal with the fits going crazy at
            #  the edges of the chip in spec direction
            # TODO: Can we add defaults to bspline_profile so that we
            #  don't have to instantiate invvar and profile_basis
            spec_bspl, spec_gpm_fit, spec_flat_fit, _, exit_status \
                    = fitting.bspline_profile(spec_coo_data, spec_flat_data, spec_ivar_data,
                                            np.ones_like(spec_coo_data), ingpm=spec_gpm_data,
                                            nord=4, upper=logrej, lower=logrej,
                                            kwargs_bspline={'bkspace': spec_samp_fine},
                                            kwargs_reject={'groupbadpix': True, 'maxrej': 5})

            if exit_status > 1:
                # TODO -- MAKE A FUNCTION
                msgs.warn('Flat-field spectral response bspline fit failed!  Not flat-fielding '
                          'slit {0} and continuing!'.format(slit_spat))
                self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(self.slits.mask[slit_idx], 'BADFLATCALIB')
                continue

            # Debugging/checking spectral fit
            if debug:
                fitting.bspline_qa(spec_coo_data, spec_flat_data, spec_bspl, spec_gpm_fit,
                                 spec_flat_fit, xlabel='Spectral Pixel', ylabel='log(flat counts)',
                                 title='Spectral Fit for slit={:d}'.format(slit_spat))

            if sticky:
                # Add rejected pixels to gpm
                gpm[spec_gpm] = (spec_gpm_fit & spec_gpm_data)[np.argsort(spec_srt, kind='stable')]

            # Construct the model of the flat-field spectral shape
            # including padding on either side of the slit.
            spec_model[...] = 1.
            spec_model[onslit_padded] = np.exp(spec_bspl.value(spec_coo[onslit_padded])[0])
            # ----------------------------------------------------------

            # ----------------------------------------------------------
            # To fit the spatial response, first normalize out the
            # spectral response, and then collapse the slit spectrally.

            # Normalize out the spectral shape of the flat
            norm_spec[...] = 1.
            norm_spec[onslit_padded] = rawflat[onslit_padded] \
                                            / np.fmax(spec_model[onslit_padded],1.0)

            # Find pixels fot fit in the spatial direction:
            #   - Fit pixels in the padded slit that haven't been masked
            #     by the BPM
            spat_gpm = onslit_padded & gpm #& (rawflat < nonlinear_counts)
            #   - Fit pixels with non-zero flux and less than 70% above
            #     the average spectral profile.
            spat_gpm &= (norm_spec > 0.0) & (norm_spec < 1.7)
            #   - Determine maximum counts in median filtered flat
            #     spectrum model.
            spec_interp = interpolate.interp1d(spec_coo_data, spec_flat_fit, kind='linear',
                                               assume_sorted=True, bounds_error=False,
                                               fill_value=-np.inf)
            spec_sm = utils.fast_running_median(np.exp(spec_interp(np.arange(nspec))),
                                                np.fmax(np.ceil(0.10*nspec).astype(int),10))
            #   - Only fit pixels with at least values > 10% of this maximum and no less than 1.
            spat_gpm &= (spec_model > 0.1*np.amax(spec_sm)) & (spec_model > 1.0)

            # Report
            spat_nfit = np.sum(spat_gpm)
            spat_ntot = np.sum(onslit_padded)
            msgs.info('Spatial fit of flatfield for {0}/{1} '.format(spat_nfit, spat_ntot)
                      + ' pixels in the slit.')
            if spat_nfit/spat_ntot < 0.5:
                # TODO: Shouldn't this raise an exception or continue to the next slit instead?
                msgs.warn('Spatial fit includes only {:.1f}'.format(100*spat_nfit/spat_ntot)
                          + '% of the pixels on this slit.' + msgs.newline()
                          + '          Either the slit has many bad pixels, the model of the '
                          'spectral shape is poor, or the illumination profile is very irregular.')

            # First fit -- With initial slits
            if not np.any(spat_gpm):
                msgs.warn('Flat-field failed during normalization!  Not flat-fielding '
                          'slit {0} and continuing!'.format(slit_spat))
                self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(
                    self.slits.mask[slit_idx], 'BADFLATCALIB')
                continue

            exit_status, spat_coo_data,  spat_flat_data, spat_bspl, spat_gpm_fit, \
                spat_flat_fit, spat_flat_data_raw \
                        = self.spatial_fit(norm_spec, spat_coo_init, median_slit_widths[slit_idx],
                                           spat_gpm, gpm, debug=debug)

            if tweak_slits:
                # TODO: Should the tweak be based on the bspline fit?
                # TODO: Will this break if
                left_thresh, left_shift, self.slits.left_tweak[:,slit_idx], right_thresh, \
                    right_shift, self.slits.right_tweak[:,slit_idx] \
                        = self.tweak_slit_edges(self.slits.left_init[:,slit_idx],
                                                self.slits.right_init[:,slit_idx],
                                                spat_coo_data, spat_flat_data,
                                                method=tweak_method,
                                                thresh=tweak_slits_thresh,
                                                maxfrac=tweak_slits_maxfrac, debug=debug)
                # TODO: Because the padding doesn't consider adjacent
                #  slits, calling slit_img for individual slits can be
                #  different from the result when you construct the
                #  image for all slits. Fix this...

                # Update the onslit mask
                _slitid_img = self.slits.slit_img(slitidx=slit_idx, initial=False)
                onslit_tweak = _slitid_img == slit_spat
                # Note, we need to get the full image with the coordinates similar to spat_coo_init, otherwise, the
                # tweaked locations are biased.
                spat_coo_tweak = self.slits.spatial_coordinate_image(slitidx=slit_idx, full=True, slitid_img=_slitid_img)

                # Construct the empirical illumination profile
                # TODO This is extremely inefficient, because we only need to re-fit the illumflat, but
                #  spatial_fit does both the reconstruction of the illumination function and the bspline fitting.
                #  Only the b-spline fitting needs be reddone with the new tweaked spatial coordinates, so that would
                #  save a ton of runtime. It is not a trivial change because the coords are sorted, etc.
                exit_status, spat_coo_data, spat_flat_data, spat_bspl, spat_gpm_fit, \
                    spat_flat_fit, spat_flat_data_raw = self.spatial_fit(
                    norm_spec, spat_coo_tweak, median_slit_widths[slit_idx], spat_gpm, gpm, debug=False)

                spat_coo_final = spat_coo_tweak
            else:
                _slitid_img = slitid_img_init
                spat_coo_final = spat_coo_init
                onslit_tweak = onslit_init

            # Add an approximate pixel axis at the top
            if debug:
                # TODO: Move this into a qa plot that gets saved
                ax = fitting.bspline_qa(spat_coo_data, spat_flat_data, spat_bspl, spat_gpm_fit,
                                      spat_flat_fit, show=False)
                ax.scatter(spat_coo_data, spat_flat_data_raw, marker='.', s=1, zorder=0, color='k',
                           label='raw data')
                # Force the center of the slit to be at the center of the plot for the hline
                ax.set_xlim(-0.1,1.1)
                ax.axvline(0.0, color='lightgreen', linestyle=':', linewidth=2.0,
                           label='original left edge', zorder=8)
                ax.axvline(1.0, color='red', linestyle=':', linewidth=2.0,
                           label='original right edge', zorder=8)
                if tweak_slits and left_shift > 0:
                    label = 'threshold = {:5.2f}'.format(tweak_slits_thresh) \
                                + ' % of max of left illumprofile'
                    ax.axhline(left_thresh, xmax=0.5, color='lightgreen', linewidth=3.0,
                               label=label, zorder=10)
                    ax.axvline(left_shift, color='lightgreen', linestyle='--', linewidth=3.0,
                               label='tweaked left edge', zorder=11)
                if tweak_slits and right_shift > 0:
                    label = 'threshold = {:5.2f}'.format(tweak_slits_thresh) \
                                + ' % of max of right illumprofile'
                    ax.axhline(right_thresh, xmin=0.5, color='red', linewidth=3.0, label=label,
                               zorder=10)
                    ax.axvline(1-right_shift, color='red', linestyle='--', linewidth=3.0,
                               label='tweaked right edge', zorder=20)
                ax.legend()
                ax.set_xlabel('Normalized Slit Position')
                ax.set_ylabel('Normflat Spatial Profile')
                ax.set_title('Illumination Function Fit for slit={:d}'.format(slit_spat))
                plt.show()

            # Perform a fine correction to the spatial illumination profile
            spat_illum_fine = 1  # Default value if the fine correction is not performed
            if exit_status <= 1 and self.flatpar['slit_illum_finecorr']:
                spat_model = np.ones_like(spec_model)
                spat_model[onslit_padded] = spat_bspl.value(spat_coo_final[onslit_padded])[0]
                specspat_illum = np.fmax(spec_model, 1.0) * spat_model
                norm_spatspec = rawflat / specspat_illum
                spat_illum_fine = self.spatial_fit_finecorr(norm_spatspec, onslit_tweak, slit_idx, slit_spat, gpm, doqa=doqa)[onslit_tweak]

            # ----------------------------------------------------------
            # Construct the illumination profile with the tweaked edges
            # of the slit
            if exit_status <= 1:
                # TODO -- JFH -- Check this is ok for flexure!!
                self.msillumflat[onslit_tweak] = spat_illum_fine * spat_bspl.value(spat_coo_final[onslit_tweak])[0]
                self.list_of_spat_bsplines[slit_idx] = spat_bspl
                # No need to proceed further if we just need the illumination profile
                if spat_illum_only:
                    continue
            else:
                # Save the nada
                msgs.warn('Slit illumination profile bspline fit failed!  Spatial profile not '
                          'included in flat-field model for slit {0}!'.format(slit_spat))
                self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(self.slits.mask[slit_idx], 'BADFLATCALIB')
                continue

            # ----------------------------------------------------------
            # Fit the 2D residuals of the 1D spectral and spatial fits.
            msgs.info('Performing 2D illumination + scattered light flat field fit')

            # Construct the spectrally and spatially normalized flat
            norm_spec_spat[...] = 1.
            norm_spec_spat[onslit_tweak] = rawflat[onslit_tweak] / np.fmax(spec_model[onslit_tweak], 1.0) \
                                                    / np.fmax(self.msillumflat[onslit_tweak], 0.01)

            # Sort the pixels by their spectral coordinate. The mask
            # uses the nominal padding defined by the slits object.
            twod_gpm, twod_srt, twod_spec_coo_data, twod_flat_data \
                    = flat.sorted_flat_data(norm_spec_spat, spec_coo, gpm=onslit_tweak)
            # Also apply the sorting to the spatial coordinates
            twod_spat_coo_data = spat_coo_final[twod_gpm].ravel()[twod_srt]
            # TODO: Reset back to origin gpm if sticky is true?
            twod_gpm_data = gpm[twod_gpm].ravel()[twod_srt]
            # Only fit data with less than 30% variations
            # TODO: Make 30% a parameter?
            twod_gpm_data &= np.absolute(twod_flat_data - 1) < 0.3
            # Here we ignore the formal photon counting errors and
            # simply assume that a typical error per pixel. This guess
            # is somewhat aribtrary. We then set the rejection
            # threshold with sigrej_twod
            # TODO: Make twod_sig and twod_sigrej parameters?
            twod_sig = 0.01
            twod_ivar_data = twod_gpm_data.astype(float)/(twod_sig**2)
            twod_sigrej = 4.0

            poly_basis = basis.fpoly(2.0*twod_spat_coo_data - 1.0, npoly)

            # Perform the full 2d fit
            twod_bspl, twod_gpm_fit, twod_flat_fit, _, exit_status \
                    = fitting.bspline_profile(twod_spec_coo_data, twod_flat_data, twod_ivar_data,
                                            poly_basis, ingpm=twod_gpm_data, nord=4,
                                            upper=twod_sigrej, lower=twod_sigrej,
                                            kwargs_bspline={'bkspace': spec_samp_coarse},
                                            kwargs_reject={'groupbadpix': True, 'maxrej': 10})
            if debug:
                # TODO: Make a plot that shows the residuals in the 2D
                # image
                resid = twod_flat_data - twod_flat_fit
                goodpix = twod_gpm_fit & twod_gpm_data
                badpix = np.invert(twod_gpm_fit) & twod_gpm_data

                plt.clf()
                ax = plt.gca()
                ax.plot(twod_spec_coo_data[goodpix], resid[goodpix], color='k', marker='o',
                        markersize=0.2, mfc='k', fillstyle='full', linestyle='None',
                        label='good points')
                ax.plot(twod_spec_coo_data[badpix], resid[badpix], color='red', marker='+',
                        markersize=0.5, mfc='red', fillstyle='full', linestyle='None',
                        label='masked')
                ax.axhline(twod_sigrej*twod_sig, color='lawngreen', linestyle='--',
                           label='rejection thresholds', zorder=10, linewidth=2.0)
                ax.axhline(-twod_sigrej*twod_sig, color='lawngreen', linestyle='--', zorder=10,
                           linewidth=2.0)
#                ax.set_ylim(-0.05, 0.05)
                ax.legend()
                ax.set_xlabel('Spectral Pixel')
                ax.set_ylabel('Residuals from pixelflat 2-d fit')
                ax.set_title('Spectral Residuals for slit={:d}'.format(slit_spat))
                plt.show()

                plt.clf()
                ax = plt.gca()
                ax.plot(twod_spat_coo_data[goodpix], resid[goodpix], color='k', marker='o',
                        markersize=0.2, mfc='k', fillstyle='full', linestyle='None',
                        label='good points')
                ax.plot(twod_spat_coo_data[badpix], resid[badpix], color='red', marker='+',
                        markersize=0.5, mfc='red', fillstyle='full', linestyle='None',
                        label='masked')
                ax.axhline(twod_sigrej*twod_sig, color='lawngreen', linestyle='--',
                           label='rejection thresholds', zorder=10, linewidth=2.0)
                ax.axhline(-twod_sigrej*twod_sig, color='lawngreen', linestyle='--', zorder=10,
                           linewidth=2.0)
#                ax.set_ylim((-0.05, 0.05))
#                ax.set_xlim(-0.02, 1.02)
                ax.legend()
                ax.set_xlabel('Normalized Slit Position')
                ax.set_ylabel('Residuals from pixelflat 2-d fit')
                ax.set_title('Spatial Residuals for slit={:d}'.format(slit_spat))
                plt.show()

            # Save the 2D residual model
            twod_model[...] = 1.
            if exit_status > 1:
                msgs.warn('Two-dimensional fit to flat-field data failed!  No higher order '
                          'flat-field corrections included in model of slit {0}!'.format(slit_spat))
                self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(self.slits.mask[slit_idx], 'BADFLATCALIB')
            else:
                twod_model[twod_gpm] = twod_flat_fit[np.argsort(twod_srt, kind='stable')]
                twod_gpm_out[twod_gpm] = twod_gpm_fit[np.argsort(twod_srt, kind='stable')]


            # Construct the full flat-field model
            # TODO: Why is the 0.05 here for the illumflat compared to the 0.01 above?
            self.flat_model[onslit_tweak] = twod_model[onslit_tweak] \
                                        * np.fmax(self.msillumflat[onslit_tweak], 0.05) \
                                        * np.fmax(spec_model[onslit_tweak], 1.0)

            # Check for infinities and NaNs in the flat-field model
            winfnan = np.where(np.logical_not(np.isfinite(self.flat_model[onslit_tweak])))
            if winfnan[0].size != 0:
                msgs.warn('There are {0:d} pixels with non-finite values in the flat-field model '
                          'for slit {1:d}!'.format(winfnan[0].size, slit_spat) + msgs.newline() +
                          'These model pixel values will be set to the raw pixel value.')
                self.flat_model[np.where(onslit_tweak)[0][winfnan]] = rawflat[np.where(onslit_tweak)[0][winfnan]]
            # Check for unrealistically high or low values of the model
            whilo = np.where((self.flat_model[onslit_tweak] >= nonlinear_counts) |
                             (self.flat_model[onslit_tweak] <= 0.0))
            if whilo[0].size != 0:
                msgs.warn('There are {0:d} pixels with unrealistically high or low values in the flat-field model '
                          'for slit {1:d}!'.format(whilo[0].size, slit_spat) + msgs.newline() +
                          'These model pixel values will be set to the raw pixel value.')
                self.flat_model[np.where(onslit_tweak)[0][whilo]] = rawflat[np.where(onslit_tweak)[0][whilo]]

            # Construct the pixel flat
            #trimmed_slitid_img_anew = self.slits.slit_img(pad=-trim, slitidx=slit_idx)
            #onslit_trimmed_anew = trimmed_slitid_img_anew == slit_spat
            self.mspixelflat[onslit_tweak] = rawflat[onslit_tweak] * utils.inverse(self.flat_model[onslit_tweak])
            # TODO: Add some code here to treat the edges and places where fits
            #  go bad?

            # Minimum wavelength?
            if self.flatpar['pixelflat_min_wave'] is not None and self.waveimg is not None:
                bad_wv = self.waveimg[onslit_tweak] < self.flatpar['pixelflat_min_wave']
                self.mspixelflat[np.where(onslit_tweak)[0][bad_wv]] = 1.
            # Maximum wavelength?
            if self.flatpar['pixelflat_max_wave'] is not None and self.waveimg is not None:
                bad_wv = self.waveimg[onslit_tweak] > self.flatpar['pixelflat_max_wave']
                self.mspixelflat[np.where(onslit_tweak)[0][bad_wv]] = 1.

        # No need to continue if we're just doing the spatial illumination
        if spat_illum_only:
            return

        # Set the pixelflat to 1.0 wherever the flat was nonlinear
        self.mspixelflat[rawflat >= nonlinear_counts] = 1.0
        # Set the pixelflat to 1.0 within trim pixels of all the slit edges
        trimmed_slitid_img_new = self.slits.slit_img(pad=-trim, initial=False)
        tweaked_slitid_img = self.slits.slit_img(initial=False)
        self.mspixelflat[(trimmed_slitid_img_new < 0) & (tweaked_slitid_img > 0)] = 1.0

        # Do not apply pixelflat field corrections that are greater than
        # 100% to avoid creating edge effects, etc.
        self.mspixelflat = np.clip(self.mspixelflat, 0.5, 2.0)

        # Calculate the relative spectral illumination, if requested
        if self.flatpar['slit_illum_relative']:
            self.spec_illum = self.spectral_illumination(twod_gpm_out, debug=debug)

    def spatial_fit(self, norm_spec, spat_coo, median_slit_width, spat_gpm, gpm, debug=False):
        """
        Perform the spatial fit

        Args:
            norm_spec (`numpy.ndarray`_):
            spat_coo (`numpy.ndarray`_):
                Spatial coordinate array
            median_slit_width (:obj:`float`):
            spat_gpm (`numpy.ndarray`_):
            gpm (`numpy.ndarray`_):
            debug (bool, optional):

        Returns:
            tuple: 7 objects
                 - exit_status (int):
                 - spat_coo_data
                 - spat_flat_data
                 - spat_bspl (:class:`~pypeit.bspline.bspline.bspline`): Bspline model of the spatial fit.  Used for illumflat
                 - spat_gpm_fit
                 - spat_flat_fit
                 - spat_flat_data_raw
        """

        # Construct the empirical illumination profile
        _spat_gpm, spat_srt, spat_coo_data, spat_flat_data_raw, spat_flat_data \
            = flat.construct_illum_profile(norm_spec, spat_coo, median_slit_width,
                                           spat_gpm=spat_gpm,
                                           spat_samp=self.flatpar['spat_samp'],
                                           illum_iter=self.flatpar['illum_iter'],
                                           illum_rej=self.flatpar['illum_rej'],
                                           debug=debug)

        if self.flatpar['rej_sticky']:
            # Add rejected pixels to gpm
            gpm[spat_gpm] &= (spat_gpm & _spat_gpm)[spat_gpm]

        # Make sure that the normalized and filtered flat is finite!
        if np.any(np.invert(np.isfinite(spat_flat_data))):
            msgs.error('Inifinities in slit illumination function computation!')

        # Determine the breakpoint spacing from the sampling of the
        # spatial coordinates. Use breakpoints at a spacing of a
        # 1/10th of a pixel, but do not allow a bsp smaller than
        # the typical sampling. Use the bspline class to determine
        # the breakpoints:
        spat_bspl = bspline.bspline(spat_coo_data, nord=4,
                                    bkspace=np.fmax(1.0 / median_slit_width / 10.0,
                                                    1.2 * np.median(np.diff(spat_coo_data))))
        # TODO: Can we add defaults to bspline_profile so that we
        #  don't have to instantiate invvar and profile_basis
        spat_bspl, spat_gpm_fit, spat_flat_fit, _, exit_status \
            = fitting.bspline_profile(spat_coo_data, spat_flat_data,
                                    np.ones_like(spat_flat_data),
                                    np.ones_like(spat_flat_data), nord=4, upper=5.0,
                                    lower=5.0, fullbkpt=spat_bspl.breakpoints)
        # Return
        return exit_status, spat_coo_data, spat_flat_data, spat_bspl, spat_gpm_fit, \
               spat_flat_fit, spat_flat_data_raw

    def spatial_fit_finecorr(self, normed, onslit_tweak, slit_idx, slit_spat, gpm,
                             slit_trim=3, tolerance=0.1, doqa=False):
        """
        Generate a relative scaling image for a slicer IFU. All
        slits are scaled relative to a reference slit, specified in
        the spectrograph settings file.

        Parameters
        ----------
        normed : `numpy.ndarray`_
            Raw flat field image, normalized by the spectral and spatial illuminations.
        onslit_tweak : `numpy.ndarray`_
            mask indicating which pixels are on the slit (True = on slit)
        slit_idx : int
            Slit number (0-indexed)
        slit_spat : int
            Spatial ID of the slit
        gpm : `numpy.ndarray`_
            Good pixel mask
        slit_txt : str
            if pypeline is "Echelle", then slit_txt should be set to "order", otherwise, use "slit"
        slit_trim : int, optional
            Trim the slit edges by this number of pixels during the fitting. Note that the
            fit will be evaluated on the pixels indicated by onslit_tweak.
            A positive number trims the slit edges, a negative number pads the slit edges.
        tolerance : float, optional
            Tolerance for the relative scaling of the slits. A value of 0.1 means that the
            relative scaling of the slits must be within 10% of unity. Any data outside of
            this tolerance will be masked.
        doqa : :obj:`bool`, optional:
            Save the QA?

        Returns
        -------
        illumflat_finecorr: `numpy.ndarray`_
            An image (same shape as normed) containing the fine correction to the spatial illumination profile
        """
        # check id self.waveimg is available
        if self.waveimg is None:
            msgs.warn("Cannot perform the fine correction to the spatial illumination without the wavelength image.")
            return
        # TODO :: Include fit_order in the parset??
        fit_order = np.array([3, 6])
        slit_txt = self.slits.slitord_txt
        slit_ordid = self.slits.slitord_id[slit_idx]
        msgs.info(f"Performing a fine correction to the spatial illumination ({slit_txt} {slit_ordid})")
        # initialise
        illumflat_finecorr = np.ones_like(self.rawflatimg.image)
        # Trim the edges by a few pixels to avoid edge effects
        onslit_tweak_trim = self.slits.slit_img(pad=-slit_trim, slitidx=slit_idx, initial=False) == slit_spat
        # Setup
        slitimg = (slit_spat + 1) * onslit_tweak.astype(int) - 1  # Need to +1 and -1 so that slitimg=-1 when off the slit

        left, right, msk = self.slits.select_edges(flexure=self.wavetilts.spat_flexure if self.wavetilts is not None else 0.0)
        this_left = left[:, slit_idx]
        this_right = right[:, slit_idx]
        slitlen = int(np.median(this_right - this_left))

        # Generate the coordinates to evaluate the fit
        this_slit = np.where(onslit_tweak & self.rawflatimg.select_flag(invert=True) & (self.waveimg!=0.0))
        this_wave = self.waveimg[this_slit]
        xpos_img = self.slits.spatial_coordinate_image(slitidx=slit_idx,
                                                       slitid_img=slitimg,
                                                       flexure_shift=self.wavetilts.spat_flexure if self.wavetilts is not None else 0.0)
        # Generate the trimmed versions for fitting
        this_slit_trim = np.where(onslit_tweak_trim & self.rawflatimg.select_flag(invert=True))
        this_wave_trim = self.waveimg[this_slit_trim]
        wave_min, wave_max = this_wave_trim.min(), this_wave_trim.max()
        ypos_fit = (this_wave_trim - wave_min) / (wave_max - wave_min)
        xpos_fit = xpos_img[this_slit_trim]
        # Evaluation coordinates
        ypos = (this_wave - wave_min) / (wave_max - wave_min)  # Need to use the same wave_min and wave_max as the fitting coordinates
        xpos = xpos_img[this_slit]

        # Mask the edges and fit
        gpmfit = gpm[this_slit_trim]
        # Trim by 5% of the slit length, or at least slit_trim pixels
        xfrac = 0.05
        if xfrac * slitlen < slit_trim:
            xfrac = slit_trim/slitlen
        gpmfit[np.where((xpos_fit < xfrac) | (xpos_fit > 1-xfrac))] = False
        # If the data deviate too much from unity, mask them. We're only interested in a
        # relative correction that's less than ~10% from unity.
        gpmfit[np.where((normed[this_slit_trim] < 1-tolerance) | (normed[this_slit_trim] > 1 + tolerance))] = False
        # Perform the full fit
        fullfit = fitting.robust_fit(xpos_fit, normed[this_slit_trim], fit_order, x2=ypos_fit,
                                     in_gpm=gpmfit, function='legendre2d', upper=2, lower=2, maxdev=1.0,
                                     minx=0.0, maxx=1.0, minx2=0.0, maxx2=1.0)

        # Generate the fine correction image and store the result
        if fullfit.success == 1:
            self.list_of_finecorr_fits[slit_idx] = fullfit
            illumflat_finecorr[this_slit] = fullfit.eval(xpos, ypos)
        else:
            msgs.warn(f"Fine correction to the spatial illumination failed for {slit_txt} {slit_ordid}")
            return illumflat_finecorr

        # If corrections exceed the tolerance, then clip them to the level of the tolerance
        illumflat_finecorr = np.clip(illumflat_finecorr, 1-tolerance, 1+tolerance)

        # Prepare QA
        if doqa:
            prefix = "Spatillum_FineCorr_"
            if self.spat_illum_only:
                prefix += "illumflat_"
            outfile = qa.set_qa_filename(prefix+self.calib_key, 'spatillum_finecorr', slit=slit_spat,
                                         out_dir=self.qa_path)
            title = f"Fine correction to spatial illumination ({slit_txt} {slit_ordid})"
            normed[np.logical_not(onslit_tweak)] = 1  # For the QA, make everything off the slit equal to 1
            spatillum_finecorr_qa(normed, illumflat_finecorr, this_left, this_right, ypos_fit, this_slit_trim,
                                  outfile=outfile, title=title, half_slen=slitlen//2)
        return illumflat_finecorr

    def extract_structure(self, rawflat_orig, slit_trim=3):
        """
        Generate a relative scaling image for a slicer IFU. All
        slits are scaled relative to a reference slit, specified in
        the spectrograph settings file.

        Parameters
        ----------
        rawflat_orig : `numpy.ndarray`_
            The original raw image of the flatfield
        slit_trim : int, optional
            Trim the slit edges by this number of pixels during the fitting. Note that the
            fit will be evaluated on the pixels indicated by onslit_tweak.
            A positive number trims the slit edges, a negative number pads the slit edges.

        Returns
        -------
        ff_struct: `numpy.ndarray`_
            An image containing the detector structure (i.e. the raw flatfield image
            divided by the spectral and spatial illumination profile fits).
        """
        msgs.info("Extracting flatfield structure")

        # check if the waveimg is available
        if self.waveimg is None:
            msgs.error("Cannot perform the extraction of the flatfield structure without the wavelength image.")

        # Build the mask and make a temporary instance of FlatImages
        bpmflats = self.build_mask()
        # Initialise bad splines (for when the fit goes wrong)
        if self.list_of_spat_bsplines is None:
            self.list_of_spat_bsplines = [bspline.bspline(None) for all in self.slits.spat_id]
        if self.list_of_finecorr_fits is None:
            self.list_of_finecorr_fits = [fitting.PypeItFit(None) for all in self.slits.spat_id]
        # Generate a dummy FlatImages
        tmp_flats = FlatImages(illumflat_raw=self.rawflatimg.image,
                               illumflat_spat_bsplines=np.asarray(self.list_of_spat_bsplines),
                               illumflat_finecorr=np.asarray(self.list_of_finecorr_fits),
                               illumflat_bpm=bpmflats, PYP_SPEC=self.spectrograph.name,
                               spat_id=self.slits.spat_id)
        # Divide by the spatial profile
        spat_illum = tmp_flats.fit2illumflat(self.slits, frametype='illum', finecorr=False)
        spat_illum *= tmp_flats.fit2illumflat(self.slits, frametype='illum', finecorr=True)
        rawflat = rawflat_orig * utils.inverse(spat_illum)
        # Now fit the spectral profile
        # TODO: Should this be *any* flag, or just BPM?
        gpm = self.rawflatimg.select_flag(flag='BPM', invert=True)
        scale_model = illum_profile_spectral(rawflat, self.waveimg, self.slits,
                                             slit_illum_ref_idx=self.flatpar['slit_illum_ref_idx'],
                                             model=None, gpmask=gpm, skymask=None, trim=self.flatpar['slit_trim'],
                                             flexure=self.wavetilts.spat_flexure if self.wavetilts is not None else 0.0,
                                             smooth_npix=self.flatpar['slit_illum_smooth_npix'])
        # Trim the edges by a few pixels to avoid edge effects
        onslits_trim = gpm & (self.slits.slit_img(pad=-slit_trim, initial=False) != -1)
        onslits = (self.waveimg != 0.0) & gpm
        # Construct a wavelength array
        minwv = np.min(self.waveimg[onslits])
        maxwv = np.max(self.waveimg)
        wavebins = np.linspace(minwv, maxwv, self.slits.nspec)
        # Correct the raw flat for spatial illumination, then generate a spectrum
        rawflat_corr = rawflat * utils.inverse(scale_model)
        hist, edge = np.histogram(self.waveimg[onslits_trim], bins=wavebins, weights=rawflat_corr[onslits_trim])
        cntr, edge = np.histogram(self.waveimg[onslits_trim], bins=wavebins)
        cntr = cntr.astype(float)
        spec_ref = hist * utils.inverse(cntr)
        wave_ref = 0.5 * (wavebins[1:] + wavebins[:-1])
        # Create a 1D model of the spectrum and assign a flux to each detector pixel
        spec_model = np.ones_like(rawflat)
        spec_model[onslits] = interpolate.interp1d(wave_ref, spec_ref, kind='linear', bounds_error=False,
                                                   fill_value="extrapolate")(self.waveimg[onslits])
        # Apply relative scale
        spec_model *= scale_model
        # Divide model spectrum of pixelflat, and the small-scale pixel-to-pixel sensitivity variations,
        # to uncover the large scale detector structure
        ff_struct = rawflat * utils.inverse(spec_model) * utils.inverse(self.mspixelflat)
        return ff_struct

    def spectral_illumination(self, gpm=None, debug=False):
        """
        Generate a relative scaling image for a slicer IFU. All
        slits are scaled relative to a reference slit, specified in
        the spectrograph settings file.

        Parameters
        ----------
        gpm : `numpy.ndarray`_, None
            Good pixel mask
        debug : bool
            Debug the routine

        Returns
        -------
        scale_model: `numpy.ndarray`_
            An image containing the appropriate scaling
        """
        msgs.info("Deriving spectral illumination profile")
        # check if the waveimg is available
        if self.waveimg is None:
            msgs.warn("Cannot perform the spectral illumination without the wavelength image.")
            return None
        msgs.info('Performing a joint fit to the flat-field response')
        # Grab some parameters
        trim = self.flatpar['slit_trim']
        rawflat = self.rawflatimg.image / (self.msillumflat * self.mspixelflat)
        # Grab the GPM and the slit images
        if gpm is None:
            # TODO: Should this be *any* flag, or just BPM?
            gpm = self.rawflatimg.select_flag(flag='BPM', invert=True)

        # Obtain relative spectral illumination
        return illum_profile_spectral(rawflat, self.waveimg, self.slits,
                                      slit_illum_ref_idx=self.flatpar['slit_illum_ref_idx'],
                                      model=None, gpmask=gpm, skymask=None, trim=trim,
                                      flexure=self.wavetilts.spat_flexure  if self.wavetilts is not None else 0.0,
                                      smooth_npix=self.flatpar['slit_illum_smooth_npix'],
                                      debug=debug)

    def tweak_slit_edges(self, left, right, spat_coo, norm_flat, method='threshold', thresh=0.93,
                         maxfrac=0.1, debug=False):
        r"""
        Tweak the slit edges based on the normalized slit illumination profile.

        Args:
            left (`numpy.ndarray`_):
                Array with the left slit edge for a single slit. Shape is
                :math:`(N_{\rm spec},)`.
            right (`numpy.ndarray`_):
                Array with the right slit edge for a single slit. Shape
                is :math:`(N_{\rm spec},)`.
            spat_coo (`numpy.ndarray`_):
                Spatial pixel coordinates in fractions of the slit width
                at each spectral row for the provided normalized flat
                data. Coordinates are relative to the left edge (with the
                left edge at 0.). Shape is :math:`(N_{\rm flat},)`.
                Function assumes the coordinate array is sorted.
            norm_flat (`numpy.ndarray`_)
                Normalized flat data that provide the slit illumination
                profile. Shape is :math:`(N_{\rm flat},)`.
            method (:obj:`str`, optional):
                Method to use for tweaking the slit edges. Options are:

                    - ``'threshold'``: Use the threshold to set the slit edge
                      and then shift it to the left or right based on the
                      illumination profile.

                    - ``'gradient'``: Use the gradient of the illumination
                      profile to set the slit edge and then shift it to the left
                      or right based on the illumination profile.

            thresh (:obj:`float`, optional):
                Threshold of the normalized flat profile at which to
                place the two slit edges.
            maxfrac (:obj:`float`, optional):
                The maximum fraction of the slit width that the slit edge
                can be adjusted by this algorithm. If ``maxfrac = 0.1``,
                this means the maximum change in the slit width (either
                narrowing or broadening) is 20% (i.e., 10% for either
                edge).
            debug (:obj:`bool`, optional):
                Show flow interrupting plots that show illumination
                profile in the case of a failure and the placement of the
                tweaked edge for each side of the slit regardless.

        Returns:
            tuple: Returns six objects:

                - The threshold used to set the left edge
                - The fraction of the slit that the left edge is shifted to the
                  right
                - The adjusted left edge
                - The threshold used to set the right edge
                - The fraction of the slit that the right edge is shifted to the
                  left
                - The adjusted right edge

        """
        # TODO :: Since this is just a wrapper, and not really "core", maybe it should be moved to pypeit.flatfield?
        # Tweak the edges via the specified method
        if method == "threshold":
            return flat.tweak_slit_edges_threshold(left, right, spat_coo, norm_flat,
                                                   thresh=thresh, maxfrac=maxfrac, debug=debug)
        elif method == "gradient":
            return flat.tweak_slit_edges_gradient(left, right, spat_coo, norm_flat, maxfrac=maxfrac, debug=debug)
        else:
            msgs.error("Method for tweaking slit edges not recognized: {0}".format(method))


class SlitlessFlat:
    """
    Class to generate a slitless pixel flat-field calibration image.

    Args:
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
                The class holding the metadata for all the frames.
        slitless_rows (`numpy.ndarray`_):
                Boolean array selecting the rows in the fitstbl that
                correspond to the slitless frames.
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
                The spectrograph object.
        par (:class:`~pypeit.par.pypeitpar.CalibrationsPar`):
            Parameter set defining optional parameters of PypeIt's algorithms
            for Calibrations
        qa_path (`Path`_):
            Path for the QA diagnostics.

    """

    def __init__(self, fitstbl, slitless_rows, spectrograph, par, qa_path=None):

        self.fitstbl = fitstbl
        # Boolean array selecting the rows in the fitstbl that correspond to the slitless frames.
        self.slitless_rows = slitless_rows
        self.spectrograph = spectrograph
        self.par = par
        self.qa_path = qa_path

    def slitless_pixflat_fname(self):
        """
        Generate the name of the slitless pixel flat file.

        Returns:
            :obj:`str`: The name of the slitless pixel flat

        """
        if len(self.slitless_rows) == 0:
            msgs.error('No slitless_pixflat frames found. Cannot generate the slitless pixel flat file name.')

        # generate the slitless pixel flat file name
        spec_name = self.fitstbl.spectrograph.name
        date = self.fitstbl.construct_obstime(self.slitless_rows[0]).iso.split(' ')[0].replace('-', '') if \
            self.fitstbl[self.slitless_rows][0]['mjd'] is not None else '00000000'
        # setup info to add to the filename
        dispname = '' if 'dispname' not in self.spectrograph.configuration_keys() else \
            f"_{self.fitstbl[self.slitless_rows[0]]['dispname'].replace('/', '_').replace(' ', '_').replace('(', '').replace(')', '').replace(':', '_').replace('+', '_')}"
        dichroic = '' if 'dichroic' not in self.spectrograph.configuration_keys() else \
            f"_d{self.fitstbl[self.slitless_rows[0]]['dichroic']}"
        binning = self.fitstbl[self.slitless_rows[0]]['binning'].replace(',', 'x')
        # file name
        return f'pixelflat_{spec_name}{dispname}{dichroic}_{binning}_{date}.fits'

    def make_slitless_pixflat(self, msbias=None, msdark=None, calib_dir=None, write_qa=False, show=False):
        """
        Generate and save to disc a slitless pixel flat-field calibration images.
        The pixel flat file will have one extension per detector, even in the case of a mosaic.
        Contrary to the regular calibration flow, the slitless pixel flat is created for all detectors
        of the current spectrograph at once, and not only the one for the current detector.
        Since the slitless pixel flat images are saved to disc, this approach helps with the I/O
        This is a method is used in `~pypeit.calibrations.get_flats()`.

        Note: par['flatfield']['pixelflat_file'] is updated in this method.

        Args:
            msbias (:class:`~pypeit.images.buildimage.BiasImage`, optional):
                Bias image for bias subtraction; passed to
                :func:`~pypeit.images.buildimage.buildimage_fromlist()`
            msdark (:class:`~pypeit.images.buildimage.DarkImage`, optional):
                Dark-current image; passed to
                :func:`~pypeit.images.buildimage.buildimage_fromlist()`
            calib_dir (`Path`_):
                Path for the processed calibration files.
            write_qa (:obj:`bool`, optional):
                Write QA plots to disk?
            show (:obj:`bool`, optional):
                Show the diagnostic plots?

        Returns:
            :obj:`str`: The name of the slitless pixel flat file that was generated.

        """

        # First thing first, check if the user has provided slitless_pixflat frames
        if len(self.slitless_rows) == 0:
            # return unchanged self.par['flatfield']['pixelflat_file']
            return self.par['flatfield']['pixelflat_file']

        # all detectors of this spectrograph
        _detectors = np.array(self.spectrograph.select_detectors())

        # Check if a user-provided slitless pixelflat already exists for the current detectors
        if self.par['flatfield']['pixelflat_file'] is not None:
            _pixel_flat_file = dataPaths.pixelflat.get_file_path(self.par['flatfield']['pixelflat_file'],
                                                                 return_none=True)

            if _pixel_flat_file is not None:
                # get detector names
                detnames = np.array([self.spectrograph.get_det_name(_det) for _det in _detectors])
                # open the file
                with io.fits_open(_pixel_flat_file) as hdu:
                    # list of available detectors in the pixel flat file
                    file_detnames = [h.name.split('-')[0] for h in hdu]
                    # check if the current detnames are in the list
                    in_file = np.array([d in file_detnames for d in detnames])
                    # if all detectors are in the file, return
                    if np.all(in_file):
                        msgs.info(f"Both slitless_pixflat frames and user-defined file found. "
                                  f"The user-defined file will be used: {self.par['flatfield']['pixelflat_file']}")
                        # return unchanged self.par['flatfield']['pixelflat_file']
                        return self.par['flatfield']['pixelflat_file']
                    else:
                        # get the detectors that are not in the file
                        _detectors = _detectors[np.logical_not(in_file)]
                        detnames = detnames[np.logical_not(in_file)]
                        msgs.info(f'Both slitless_pixflat frames and user-defined file found, but the '
                                  f'following detectors are not in the file: {detnames}. Using the '
                                  f'slitless_pixflat frames to generate the missing detectors.')

        # make the slitless pixel flat
        pixflat_norm_list = []
        detname_list = []
        for _det in _detectors:
            # Parse the raw slitless pixelflat frames. Note that this is spectrograph dependent.
            # If the method does not exist in the specific spectrograph class, nothing will happen
            this_raw_idx = self.spectrograph.parse_raw_files(self.fitstbl[self.slitless_rows], det=_det,
                                                             ftype='slitless_pixflat')
            if len(this_raw_idx) == 0:
                msgs.warn(f'No raw slitless_pixflat frames found for {self.spectrograph.get_det_name(_det)}. '
                          f'Continuing...')
                continue
            this_raw_files = self.fitstbl.frame_paths(self.slitless_rows[this_raw_idx])
            msgs.info(f'Creating slitless pixel-flat calibration frame '
                      f'for {self.spectrograph.get_det_name(_det)} using files: ')
            for f in this_raw_files:
                msgs.prindent(f'{Path(f).name}')

            # Reset the BPM
            msbpm = self.spectrograph.bpm(this_raw_files[0], _det, msbias=msbias if self.par['bpm_usebias'] else None)

            # trace image
            traceimg = buildimage.buildimage_fromlist(self.spectrograph, _det, self.par['traceframe'],
                                                      [this_raw_files[0]], dark=msdark, bias=msbias, bpm=msbpm)
            # slit edges
            # we need to change some parameters for the slit edge tracing
            edges_par = deepcopy(self.par['slitedges'])
            # lower the threshold for edge detection
            edges_par['edge_thresh'] = 50.
            # this is used for longslit (i.e., no pca)
            edges_par['sync_predict'] = 'nearest'
            # remove spurious edges by setting a large minimum slit gap (20% of the detector size
            platescale = parse.parse_binning(traceimg.detector.binning)[1] * traceimg.detector['platescale']
            edges_par['minimum_slit_gap'] = 0.2 * traceimg.image.shape[1] * platescale
            # if no slits are found the bound_detector parameter add 2 traces at the detector edges
            edges_par['bound_detector'] = True
            # set the buffer to 0
            edges_par['det_buffer'] = 0
            _spectrograph = deepcopy(self.spectrograph)
            # need to treat this as a MultiSlit spectrograph (no echelle parameters used)
            _spectrograph.pypeline = 'MultiSlit'
            edges = edgetrace.EdgeTraceSet(traceimg, _spectrograph, edges_par, auto=True)
            slits = edges.get_slits()
            if show:
                edges.show(title='Slitless flat edge tracing')
            #
            # flat image
            slitless_pixel_flat = buildimage.buildimage_fromlist(self.spectrograph, _det, self.par['slitless_pixflatframe'],
                                                                 this_raw_files, dark=msdark, bias=msbias, bpm=msbpm)

            # increase saturation threshold (some hires slitless flats are very bright)
            slitless_pixel_flat.detector.saturation *= 1.5
            # Initialise the pixel flat
            flatpar = deepcopy(self.par['flatfield'])
            # do not tweak the slits
            flatpar['tweak_slits'] = False
            flatpar['slit_illum_finecorr'] = False
            pixelFlatField = FlatField(slitless_pixel_flat, self.spectrograph, flatpar, slits, wavetilts=None,
                                       wv_calib=None, slitless=True, qa_path=self.qa_path)

            # Generate
            pixelflatImages = pixelFlatField.run(doqa=write_qa, show=show)
            pixflat_norm_list.append(pixelflatImages.pixelflat_norm)
            detname_list.append(self.spectrograph.get_det_name(_det))

        if len(detname_list) > 0:
            # get the pixel flat file name
            if self.par['flatfield']['pixelflat_file'] is not None and _pixel_flat_file is not None:
                fname = self.par['flatfield']['pixelflat_file']
            else:
                fname = self.slitless_pixflat_fname()
                # file will be saved in the reduction directory, but also cached in the data/pixelflats folder
                # therefore we update self.par['flatfield']['pixelflat_file'] to the new file,
                # so that it can be used for the rest of the reduction and for the other files in the same run
                self.par['flatfield']['pixelflat_file'] = fname

            # Save the result
            write_pixflat_to_fits(pixflat_norm_list, detname_list, self.spectrograph.name,
                                  calib_dir.parent if calib_dir is not None else Path('.').absolute(),
                                  fname, to_cache=True)

        return self.par['flatfield']['pixelflat_file']


def spatillum_finecorr_qa(normed, finecorr, left, right, ypos, cut, outfile=None, title=None, half_slen=50):
    """
    Plot the QA for the fine correction fits to the spatial illumination profile

    Parameters
    ----------
    normed : `numpy.ndarray`_
        Image data with the coarse spatial illumination profile divided out (normalised in the spectral direction)
    finecorr : `numpy.ndarray`_
        Image containing the fine correction
    left : `numpy.ndarray`_
        Left slit edge
    right : `numpy.ndarray`_
        Right slit edge
    ypos : `numpy.ndarray`_
        Spectral coordinate (from 0 to 1, where 0=blue wavelength, 1=red wavelength)
    cut : tuple
        A 2-tuple, containing the (x, y) coordinates of the pixels on the slit.
    outfile : str, optional
        Output file name
    title : str, optional
        A title to be printed on the QA
    half_slen : int, optional
        The sampling size for the spatial profile. Should be about half the slit length.
        In this case, each output pixel shown contains about 2 detector pixels.
    """
    plt.rcdefaults()
    plt.rcParams['font.family'] = 'serif'

    msgs.info("Generating QA for spatial illumination fine correction")
    # Setup some plotting variables
    nseg = 10  # Number of segments to plot in QA - needs to be large enough so the fine correction is approximately linear in between adjacent segments
    colors = plt.cm.jet(np.linspace(0, 1, nseg))
    # Setup the spatial binning
    spatbins = np.linspace(0, 1, half_slen)
    spatmid = 0.5 * (spatbins[1:] + spatbins[:-1])

    xmn, xmx, ymn, ymx = np.min(cut[0]), 1+np.max(cut[0]), np.min(cut[1]), 1+np.max(cut[1])
    norm_cut = normed[xmn:xmx, ymn:ymx]
    fcor_cut = finecorr[xmn:xmx, ymn:ymx]
    vmin, vmax = max(0.95, np.min(fcor_cut)), min(1.05, np.max(fcor_cut))  # Show maximum corrections of ~5%

    # For display/visual purposes, apply a median filter to the data
    norm_cut = ndimage.median_filter(norm_cut, size=(normed.shape[0]//100, 5))

    # Plot
    fighght = 8.5
    cutrat = fighght*norm_cut.shape[1]/norm_cut.shape[0]
    plt.figure(figsize=(5 + 3.25*cutrat, fighght))
    plt.clf()
    # Single panel plot
    gs = gridspec.GridSpec(1, 5, height_ratios=[1], width_ratios=[4.0, cutrat, cutrat, cutrat, cutrat*0.25])
    ax_spec = plt.subplot(gs[0])
    # Setup the bin edges, and some plotting variables
    bins = np.linspace(0, 1, nseg + 1)
    minmod, maxmod, sep = 1.0, 1.0, 0.01
    for bb in range(nseg):
        # Histogram the data to be displayed
        wb = np.where((ypos > bins[bb]) & (ypos < bins[bb + 1]))
        cutQA = (cut[0][wb], cut[1][wb])
        xposQA = (cutQA[1] - left[cutQA[0]]) / (right[cutQA[0]] - left[cutQA[0]])
        cntr, _ = np.histogram(xposQA, bins=spatbins, weights=normed[cutQA])
        model, _ = np.histogram(xposQA, bins=spatbins, weights=finecorr[cutQA])
        nrm, _ = np.histogram(xposQA, bins=spatbins)
        cntr *= utils.inverse(nrm)
        model *= utils.inverse(nrm)
        # Make the model
        offs = bb * sep
        model += offs
        nonzero = model != offs
        if np.any(nonzero):
            minmod = minmod if minmod < np.min(model[nonzero]) else np.min(model[nonzero])
            maxmod = maxmod if maxmod > np.max(model[nonzero]) else np.max(model[nonzero])
        # Plot it!
        ax_spec.plot(spatmid, offs + cntr, linestyle='-', color=colors[bb])
        ax_spec.plot(spatmid, model, linestyle='-', color=colors[bb], alpha=0.5, linewidth=3)
    # Axes
    ax_spec.set_xlim(0.0, 1.0)
    ax_spec.set_ylim(minmod - sep, maxmod + sep)
    ax_spec.set_xlabel('Fraction of spatial slit length')
    ax_spec.minorticks_on()
    ax_spec.set_ylabel('Spatial illumination (fine correction), curves offset by {0:.3f}'.format(sep))
    if title is not None:
        ax_spec.text(0.04, 1.01, title, transform=ax_spec.transAxes,
                     ha='left', va='bottom', fontsize='medium')
    # Plot the image, model, and residual
    ax_normed = plt.subplot(gs[1])
    ax_normed.imshow(np.flipud(norm_cut), vmin=vmin, vmax=vmax)
    ax_normed.set_title("data", fontsize='small')
    ax_normed.axis('off')
    ax_fincor = plt.subplot(gs[2])
    ax_fincor.imshow(np.flipud(fcor_cut), vmin=vmin, vmax=vmax)
    ax_fincor.set_title("model", fontsize='small')
    ax_fincor.axis('off')
    ax_resid = plt.subplot(gs[3])
    # Express the deviations as a percentage
    im = ax_resid.imshow(np.flipud(norm_cut-fcor_cut)*100, vmin=(vmin-1)*100, vmax=(vmax-1)*100)
    ax_resid.set_title("diff", fontsize='small')
    ax_resid.axis('off')
    # Add a colorbar
    cax = plt.subplot(gs[4])
    cbar = plt.colorbar(im, cax=cax)#, fraction=0.046, pad=0.04)
    cbar.set_label('Percentage deviation', rotation=270, labelpad=10)
    # Finish
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
    plt.subplots_adjust(wspace=0.03, hspace=0, left=0.12, right=0.9, bottom=0.05, top=0.94)
    if outfile is None:
        plt.show()
    else:
        plt.savefig(outfile, dpi=400)
        msgs.info("Saved QA:"+msgs.newline()+outfile)

    plt.close()
    plt.rcdefaults()
    return


def detector_structure_qa(det_resp, det_resp_model, outfile=None, title="Detector Structure Correction"):
    """
    Plot the QA for the fine correction fits to the spatial illumination profile

    Parameters
    ----------
    det_resp : `numpy.ndarray`_
        Image data showing the detector structure, generated with extract_structure
    det_resp_model : `numpy.ndarray`_
        Image containing the structure correction model
    outfile : str, optional
        Output file name
    title : str, optional
        A title to be printed on the QA
    """
    plt.rcdefaults()
    plt.rcParams['font.family'] = 'serif'
    msgs.info("Generating QA for flat field structure correction")
    # Calculate the scale to be used in the plot
    # med = np.median(det_resp)
    # mad = 1.4826*np.median(np.abs(det_resp-med))
    # vmin, vmax = med-2*mad, med+2*mad
    dev = (1-np.min(det_resp_model))
    vmin, vmax = 1-2*dev, 1+2*dev

    # Plot
    fig_height = 3.0
    plt.figure(figsize=(3*fig_height, fig_height))
    plt.clf()
    # Prepare axes
    gs = gridspec.GridSpec(1, 4, height_ratios=[1], width_ratios=[1.0, 1.0, 1.0, 0.05])
    # Axes showing the observed detector response
    ax_data = plt.subplot(gs[0])
    ax_data.imshow(det_resp, origin='lower', vmin=vmin, vmax=vmax)
    ax_data.set_xlabel("data", fontsize='medium')
    ax_data.axes.xaxis.set_ticks([])
    ax_data.axes.yaxis.set_ticks([])
    # Axes showing the model fit to the detector response
    ax_modl = plt.subplot(gs[1])
    im = ax_modl.imshow(det_resp_model, origin='lower', vmin=vmin, vmax=vmax)
    ax_modl.set_title(title, fontsize='medium')
    ax_modl.set_xlabel("model", fontsize='medium')
    ax_modl.axes.xaxis.set_ticks([])
    ax_modl.axes.yaxis.set_ticks([])
    # Axes showing the residual of the detector response fit
    ax_resd = plt.subplot(gs[2])
    ax_resd.imshow(det_resp-det_resp_model, origin='lower', vmin=vmin-1, vmax=vmax-1)
    ax_resd.set_xlabel("1+data-model", fontsize='medium')
    ax_resd.axes.xaxis.set_ticks([])
    ax_resd.axes.yaxis.set_ticks([])
    # Add a colorbar
    cax = plt.subplot(gs[3])
    cbar = plt.colorbar(im, cax=cax)  # , fraction=0.046, pad=0.04)
    cbar.set_label('Deviation', rotation=270, labelpad=10)
    # Finish
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
    plt.subplots_adjust(wspace=0.03, hspace=0, left=0.05, right=0.9, bottom=0.1, top=0.9)
    if outfile is None:
        plt.show()
    else:
        plt.savefig(outfile, dpi=400)
        msgs.info("Saved QA:" + msgs.newline() + outfile)

    plt.close()
    plt.rcdefaults()
    return


def show_flats(image_list, wcs_match=True, slits=None, waveimg=None):
    """
    Show the flat-field images

    Parameters
    ----------
    image_list : list
        List of tuples containing the image data, image name and the cut levels
    wcs_match : bool, optional
        Match the WCS of the images
    slits : :class:`pypeit.slittrace.SlitTraceSet`, optional
        Slit traces to be overplotted on the images
    waveimg : `numpy.ndarray`_, optional
        Wavelength image to be overplotted on the images

    """
    display.connect_to_ginga(raise_err=True, allow_new=True)
    if slits is not None:
        left, right, mask = slits.select_edges()
        gpm = mask == 0
    # Loop me
    clear = True
    for img, name, cut in image_list:
        if img is None:
            continue
        # TODO: Add an option that shows the relevant stuff in a
        # matplotlib window.
        viewer, ch = display.show_image(img, chname=name, cuts=cut, wcs_match=wcs_match,
                                        waveimg=waveimg, clear=clear)
        if slits is not None:
            display.show_slits(viewer, ch, left[:, gpm], right[:, gpm],
                               slit_ids=slits.spat_id[gpm])
        # Turn off clear
        if clear:
            clear = False


# TODO :: This could possibly be moved to core.flat
def illum_profile_spectral(rawimg, waveimg, slits, slit_illum_ref_idx=0, smooth_npix=None, polydeg=None,
                           model=None, gpmask=None, skymask=None, trim=3, flexure=None, maxiter=5, debug=False):
    """
    Determine the relative spectral illumination of all slits.
    Currently only used for image slicer IFUs.

    Parameters
    ----------
    rawimg : `numpy.ndarray`_
        Image data that will be used to estimate the spectral relative sensitivity
    waveimg : `numpy.ndarray`_
        Wavelength image
    slits : :class:`~pypeit.slittrace.SlitTraceSet`
        Information stored about the slits
    slit_illum_ref_idx : int
        Index of slit that is used as the reference.
    smooth_npix : int, optional
        smoothing used for determining smoothly varying relative weights by sn_weights
    polydeg : int, optional
        Degree of polynomial to be used for determining relative spectral sensitivity. If None,
        coadd.smooth_weights will be used, with the smoothing length set to smooth_npix.
    model : `numpy.ndarray`_, None
        A model of the rawimg data. If None, rawimg will be used.
    gpmask : `numpy.ndarray`_, None
        Good pixel mask
    skymask : `numpy.ndarray`_, None
        Sky mask
    trim : int
        Number of pixels to trim from the edges of the slit
        when deriving the spectral illumination
    flexure : float, None
        Spatial flexure
    maxiter : :obj:`int`
        Maximum number of iterations to perform
    debug : :obj:`bool`
        Show the results of the relative spectral illumination correction

    Returns
    -------
    scale_model: `numpy.ndarray`_
        An image containing the appropriate scaling
    """
    msgs.info("Performing relative spectral sensitivity correction (reference slit = {0:d})".format(slit_illum_ref_idx))
    if polydeg is not None:
        msgs.info("Using polynomial of degree {0:d} for relative spectral sensitivity".format(polydeg))
    else:
        msgs.info("Using 'smooth_weights' algorithm for relative spectral sensitivity")
    # Setup some helpful parameters
    skymask_now = skymask if (skymask is not None) else np.ones_like(rawimg, dtype=bool)
    gpm = gpmask if (gpmask is not None) else np.ones_like(rawimg, dtype=bool)
    modelimg = model if (model is not None) else rawimg.copy()
    # Setup the slits
    slitid_img = slits.slit_img(pad=0, flexure=flexure)
    slitid_img_trim = slits.slit_img(pad=-trim, flexure=flexure)
    scaleImg = np.ones_like(rawimg)
    modelimg_copy = modelimg.copy()
    # Obtain the minimum and maximum wavelength of all slits
    mnmx_wv = np.zeros((slits.nslits, 2))
    for slit_idx, slit_spat in enumerate(slits.spat_id):
        onslit_init = (slitid_img == slit_spat)
        mnmx_wv[slit_idx, 0] = np.min(waveimg[onslit_init])
        mnmx_wv[slit_idx, 1] = np.max(waveimg[onslit_init])
    wavecen = np.mean(mnmx_wv, axis=1)
    # Sort the central wavelengths by those that are closest to the reference slit
    wvsrt = np.argsort(np.abs(wavecen - wavecen[slit_illum_ref_idx]), kind='stable')

    # Prepare wavelength array for all spectra
    dwav = np.max((mnmx_wv[:, 1] - mnmx_wv[:, 0])/slits.nspec)
    numsamp = int((np.max(mnmx_wv) - np.min(mnmx_wv)) / dwav)
    wavebins = np.linspace(np.min(mnmx_wv), np.max(mnmx_wv), numsamp)

    # Start by building a reference spectrum
    onslit_ref_trim = (slitid_img_trim == slits.spat_id[slit_illum_ref_idx]) & gpm & skymask_now
    hist, edge = np.histogram(waveimg[onslit_ref_trim], bins=wavebins, weights=modelimg_copy[onslit_ref_trim])
    cntr, edge = np.histogram(waveimg[onslit_ref_trim], bins=wavebins)
    cntr = cntr.astype(float)
    norm = utils.inverse(cntr)
    spec_ref = hist * norm
    wave_ref = 0.5 * (wavebins[1:] + wavebins[:-1])
    sn_smooth_npix = wave_ref.size // smooth_npix if (smooth_npix is not None) else wave_ref.size // 10

    # Iterate until convergence
    lo_prev, hi_prev = 1.0E-32, 1.0E32
    for rr in range(maxiter):
        # Perform two iterations:
        # (ii=0) dynamically build reference spectrum
        # (ii=1) used fixed reference spectrum to calculate the relative illumination.
        for ii in range(2):
            # Reset the relative scaling for this iteration
            relscl_model = np.ones_like(rawimg)
            # Build the relative illumination, by successively finding the slits closest in wavelength to the reference
            for ss in range(1, slits.spat_id.size):
                # Calculate the region of overlap
                onslit_b = (slitid_img_trim == slits.spat_id[wvsrt[ss]])
                onslit_b_init = (slitid_img == slits.spat_id[wvsrt[ss]])
                onslit_b_olap = onslit_b & gpm & (waveimg >= mnmx_wv[wvsrt[ss], 0]) & (waveimg <= mnmx_wv[wvsrt[ss], 1]) & skymask_now
                hist, edge = np.histogram(waveimg[onslit_b_olap], bins=wavebins, weights=modelimg_copy[onslit_b_olap])
                cntr, edge = np.histogram(waveimg[onslit_b_olap], bins=wavebins)
                cntr = cntr.astype(float)
                # Note, if ii=1 (i.e. the reference spectrum is fixed), we want to make sure that the
                # spec_ref will give a result of 1 for all wavelengths. Apply this correction first.
                if (ii == 1) and (slits.spat_id[wvsrt[ss]] == slit_illum_ref_idx):
                    # This must be the first element of the loop by construction, but throw an error just in case
                    if ss != 0:
                        msgs.error("CODING ERROR - An error has occurred in the relative spectral illumination." +
                                   msgs.newline() + "Please contact the developers.")
                    tmp_cntr = cntr * spec_ref
                    tmp_arr = hist * utils.inverse(tmp_cntr)
                    # Calculate a smooth version of the relative response
                    ref_relscale = flat.smooth_scale(tmp_arr, wave_ref=wave_ref, polydeg=polydeg, sn_smooth_npix=sn_smooth_npix)
                    # Update the reference spectrum
                    spec_ref /= ref_relscale
                # Normalise by the reference spectrum
                cntr *= spec_ref
                norm = utils.inverse(cntr)
                arr = hist * norm
                # Calculate a smooth version of the relative response
                relscale = flat.smooth_scale(arr, wave_ref=wave_ref, polydeg=polydeg, sn_smooth_npix=sn_smooth_npix)
                # Store the result
                relscl_model[onslit_b_init] = interpolate.interp1d(wave_ref, relscale, kind='linear', bounds_error=False,
                                                 fill_value="extrapolate")(waveimg[onslit_b_init])

                # Build a new reference spectrum to increase wavelength coverage of the reference spectrum (and improve S/N)
                if ii == 0:
                    onslit_ref_trim |= (onslit_b & gpm & skymask_now)
                    hist, edge = np.histogram(waveimg[onslit_ref_trim], bins=wavebins, weights=modelimg_copy[onslit_ref_trim]/relscl_model[onslit_ref_trim])
                    cntr, edge = np.histogram(waveimg[onslit_ref_trim], bins=wavebins)
                    cntr = cntr.astype(float)
                    norm = utils.inverse(cntr)
                    spec_ref = hist * norm
        minv, maxv = np.min(relscl_model), np.max(relscl_model)
        if 1/minv + maxv > lo_prev+hi_prev:
            # Adding noise, so break
            # NOTE : The best precision one might hope for is about:
            # 1.4826 * MAD(arr) / np.sqrt(sn_smooth_npix/ 10)  # /10 comes from the coadd.smooth_weights function
            break
        else:
            lo_prev, hi_prev = 1/minv, maxv
        msgs.info("Iteration {0:d} :: Minimum/Maximum scales = {1:.5f}, {2:.5f}".format(rr + 1, minv, maxv))
        # Store rescaling
        scaleImg *= relscl_model
        #rawimg_copy /= relscl_model
        modelimg_copy /= relscl_model
        if max(abs(1/minv), abs(maxv)) < 1.005:  # Relative accuracy of 0.5% is sufficient
            break
    if debug:
        embed()
        ricp = rawimg.copy()
        for ss in range(slits.spat_id.size):
            onslit_ref_trim = (slitid_img_trim == slits.spat_id[ss]) & gpm & skymask_now
            hist, edge = np.histogram(waveimg[onslit_ref_trim], bins=wavebins, weights=ricp[onslit_ref_trim]/scaleImg[onslit_ref_trim])
            histScl, edge = np.histogram(waveimg[onslit_ref_trim], bins=wavebins, weights=scaleImg[onslit_ref_trim])
            cntr, edge = np.histogram(waveimg[onslit_ref_trim], bins=wavebins)
            cntr = cntr.astype(float)
            norm = (cntr != 0) / (cntr + (cntr == 0))
            this_spec = hist * norm
            scale_ref = histScl * norm
            plt.subplot(211)
            plt.plot(wave_ref, this_spec)
            plt.xlim([3600, 4500])
            plt.subplot(212)
            plt.plot(wave_ref, scale_ref)
        plt.xlim([3600, 4500])
        plt.subplot(211)
        plt.plot(wave_ref, spec_ref, 'k--')
        plt.xlim([3600, 4500])
        plt.show()
        # Plot the relative scales of each slit
        scales_med, scales_avg = np.zeros(slits.spat_id.size), np.zeros(slits.spat_id.size)
        for ss in range(slits.spat_id.size):
            onslit_ref_trim = (slitid_img_trim == slits.spat_id[ss]) & gpm & skymask_now & (waveimg>3628) & (waveimg<4510)
            scales_med[ss] = np.median(ricp[onslit_ref_trim]/scaleImg[onslit_ref_trim])
            scales_avg[ss] = np.mean(ricp[onslit_ref_trim]/scaleImg[onslit_ref_trim])
        plt.plot(slits.spat_id, scales_med, 'bo-', label='Median')
        plt.plot(slits.spat_id, scales_avg, 'ro-', label='Mean')
        plt.legend()
        plt.show()
    return scaleImg


def merge(init_cls, merge_cls):
    """
    Merge ``merge_cls`` into ``init_cls``, and return a merged
    :class:`~pypeit.flatfield.FlatImages` class.

    If ``init_cls`` is None, the returned value is ``merge_cls``, and vice
    versa.  If an element exists in both init_cls and merge_cls, the merge_cls
    value is taken

    Parameters
    ----------
    init_cls : :class:`~pypeit.flatfield.FlatImages`
        Initial class (the elements of this class will be considered the
        default).  Can be None.
    merge_cls : :class:`~pypeit.flatfield.FlatImages`
        The non-zero elements will be merged into init_cls.  Can be None.

    Returns
    -------
    flats : :class:`~pypeit.flatfield.FlatImages`
        A new instance of the FlatImages class with merged properties.
    """
    # Check the class to be merged in is not None
    if merge_cls is None:
        return init_cls
    if init_cls is None:
        return merge_cls
    # Initialise variables
    # extract all elements that are prefixed with 'pixelflat_' or 'illumflat_'
    keys = [a for a in list(init_cls.__dict__.keys()) if '_' in a and a.split('_')[0] in ['illumflat', 'pixelflat']]
    dd = dict()
    for key in keys:
        dd[key] = None
    # Cherry pick the values from each class
    dd['PYP_SPEC'] = merge_cls.PYP_SPEC if init_cls.PYP_SPEC is None else init_cls.PYP_SPEC
    dd['spat_id'] = merge_cls.spat_id if init_cls.spat_id is None else init_cls.spat_id
    for key in keys:
        dd[key] = getattr(init_cls, key) if getattr(merge_cls, key) is None \
                    else getattr(merge_cls, key)
    # Construct the merged class
    return FlatImages(**dd)


def write_pixflat_to_fits(pixflat_norm_list, detname_list, spec_name, outdir, pixelflat_name, to_cache=True):
    """
    Write the pixel-to-pixel flat-field images to a FITS file.
    The FITS file will have an extension for each detector (never a mosaic).
    The `load_pixflat()` method read this file and transform it into a mosaic if needed.
    This image is generally used as a user-provided pixel flat-field image and ingested
    in the reduction using the `pixelflat_file` parameter in the PypeIt file.

    Args:
        pixflat_norm_list (:obj:`list`):
            List of 2D `numpy.ndarray`_ arrays containing the pixel-to-pixel flat-field images.
        detname_list (:obj:`list`):
            List of detector names.
        spec_name (:obj:`str`):
            Name of the spectrograph.
        outdir (:obj:`pathlib.Path`):
            Path to the output directory.
        pixelflat_name (:obj:`str`):
            Name of the output file to be written.
        to_cache (:obj:`bool`, optional):
            If True, the file will be written to the cache directory pypeit/data/pixflats.

    """

    msgs.info("Writing the pixel-to-pixel flat-field images to a FITS file.")

    # Check that the number of detectors matches the number of pixelflat_norm arrays
    if len(pixflat_norm_list) != len(detname_list):
        msgs.error("The number of detectors does not match the number of pixelflat_norm arrays. "
                   "The pixelflat file cannot be written.")

    # local output (reduction directory)
    pixelflat_file = outdir / pixelflat_name

    # Check if the file already exists
    old_hdus = []
    old_detnames = []
    old_hdr = None
    if pixelflat_file.exists():
        msgs.warn("The pixelflat file already exists. It will be overwritten/updated.")
        old_hdus = fits.open(pixelflat_file)
        old_detnames = [h.name.split('-')[0] for h in old_hdus]  # this has also 'PRIMARY'
        old_hdr = old_hdus[0].header

    # load spectrograph
    spec = load_spectrograph(spec_name)

    # Create the new HDUList
    _hdr = io.initialize_header(hdr=old_hdr)
    prihdu = fits.PrimaryHDU(header=_hdr)
    prihdu.header['CALIBTYP'] = (FlatImages.calib_type, 'PypeIt: Calibration frame type')
    new_hdus = [prihdu]

    extnum = 1
    for d in spec.select_detectors():
        detname = spec.get_det_name(d)
        extname = f'{detname}-PIXELFLAT_NORM'
        # update or add the detectors that we want to save
        if detname in detname_list:
            det_idx = detname_list.index(detname)
            pixflat_norm = pixflat_norm_list[det_idx]
            hdu = fits.ImageHDU(data=pixflat_norm, name=extname)
            prihdu.header[f'EXT{extnum:04d}'] = hdu.name
            new_hdus.append(hdu)
        # keep the old detectors that were not updated
        elif detname in old_detnames:
            old_det_idx = old_detnames.index(detname)
            hdu = old_hdus[old_det_idx]
            prihdu.header[f'EXT{extnum:04d}'] = hdu.name
            new_hdus.append(hdu)
        extnum += 1

    # Write the new HDUList
    new_hdulist = fits.HDUList(new_hdus)
    # Check if the directory exists
    if not pixelflat_file.parent.is_dir():
        pixelflat_file.parent.mkdir(parents=True)
    new_hdulist.writeto(pixelflat_file, overwrite=True)
    msgs.info(f'A slitless Pixel Flat file for detectors {detname_list} has been saved to {msgs.newline()}'
              f'{pixelflat_file}')

    # common msg
    add_msgs = f"add the following to your PypeIt Reduction File:{msgs.newline()}"  \
               f" [calibrations]{msgs.newline()}"  \
               f"   [[flatfield]]{msgs.newline()}"  \
               f"     pixelflat_file = {pixelflat_name}{msgs.newline()}{msgs.newline()}{msgs.newline()}"  \
               f"Please consider sharing your Pixel Flat file with the PypeIt Developers.{msgs.newline()}"  \


    if to_cache:
        # NOTE that the file saved in the cache is gzipped, while the one saved in the outdir is not
        # This prevents `dataPaths.pixelflat.get_file_path()` from returning the file saved in the outdir
        cache.write_file_to_cache(pixelflat_file, pixelflat_name+'.gz', f"pixelflats")
        msgs.info(f"The slitless Pixel Flat file has also been saved to the PypeIt cache directory {msgs.newline()}"
                  f"{str(dataPaths.pixelflat)} {msgs.newline()}"
                  f"It will be automatically used in this run. "
                  f"If you want to use this file in future runs, {add_msgs}")
    else:
        msgs.info(f"To use this file, move it to the PypeIt data directory {msgs.newline()}"
                  f"{str(dataPaths.pixelflat)} {msgs.newline()} and {add_msgs}")


def load_pixflat(pixel_flat_file, spectrograph, det, flatimages, calib_dir=None, chk_version=False):
    """
    Load a pixel flat from a file and add it to the flatimages object.
    The pixel flat file has one detector per extension, even in the case of a mosaic.
    Therefore, if this is a mosaic reduction, this script will construct a pixel flat
    mosaic. The Edges file needs to exist in the Calibration Folder, since the mosaic
    parameters are pulled from it.
    This is used in `~pypeit.calibrations.get_flats()`.

    Args:
        pixel_flat_file (:obj:`str`):
            Name of the pixel flat file.
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph object.
        det (:obj:`int`, :obj:`tuple`):
            The single detector or set of detectors in a mosaic to process.
        flatimages (:class:`~pypeit.flatfield.FlatImages`):
            The flat field images object.
        calib_dir (:obj:`str`, optional):
            The path to the calibration directory.
        chk_version (:obj:`bool`, optional):
            Check the version of the file.

    Returns:
        :class:`~pypeit.flatfield.FlatImages`: The flat images object with the pixel flat added.

    """
    # Check if the pixel flat file exists
    if pixel_flat_file is None:
        msgs.error('No pixel flat file defined. Cannot load the pixel flat!')

    # get the path
    _pixel_flat_file = dataPaths.pixelflat.get_file_path(pixel_flat_file, return_none=True)
    if _pixel_flat_file is None:
        msgs.error(f'Cannot load the pixel flat file, {pixel_flat_file}. It is not a direct path, '
                   f'a cached file, or a file that can be downloaded from a PypeIt repository.')

    # If this is a mosaic, we need to construct the pixel flat mosaic
    if isinstance(det, tuple):
        # We need to grab mosaic info from another existing calibration frame.
        # We use EdgeTraceSet image to get `tform` and `msc_ord`. Check if EdgeTraceSet file exists.
        edges_file = Path(edgetrace.EdgeTraceSet.construct_file_name(flatimages.calib_key,
                                                                     calib_dir=calib_dir)).absolute()
        if not edges_file.exists():
            msgs.error('Edges file not found in the Calibrations folder. '
                       'It is needed to grab the mosaic parameters to load and mosaic the input pixel flat!')

        # Load detector info from EdgeTraceSet file
        traceimg = edgetrace.EdgeTraceSet.from_file(edges_file, chk_version=chk_version).traceimg
        det_info = traceimg.detector
        # check that the mosaic parameters are defined
        if not np.all(np.in1d(['tform', 'msc_ord'], list(det_info.keys()))) or  \
                det_info.tform is None or det_info.msc_ord is None:
            msgs.error('Mosaic parameters are not defined in the Edges frame. Cannot load the pixel flat!')

        # read the file
        with io.fits_open(_pixel_flat_file) as hdu:
            # list of available detectors in the pixel flat file
            file_dets = [int(h.name.split('-')[0].split('DET')[1]) for h in hdu[1:]]
            # check if all detectors required for the mosaic are in the list
            if not np.all(np.in1d(list(det), file_dets)):
                msgs.error(f'Not all detectors in the mosaic are in the pixel flat file: '
                           f'{pixel_flat_file}. Cannot load the pixel flat!')

            # get the pixel flat images of only the detectors in the mosaic
            pixflat_images = np.concatenate([hdu[f'DET{d:02d}-PIXELFLAT_NORM'].data[None,:,:] for d in det])
            # construct the pixel flat mosaic
            pixflat_msc, _,_,_ = build_image_mosaic(pixflat_images, det_info.tform, order=det_info.msc_ord)
            # check that the mosaic has the correct shape
            if pixflat_msc.shape != traceimg.image.shape:
                msgs.error('The constructed pixel flat mosaic does not have the correct shape. '
                           'Cannot load this pixel flat as a mosaic!')
            msgs.info(f'Using pixelflat file: {pixel_flat_file} '
                      f'for {spectrograph.get_det_name(det)}.')
            nrm_image = FlatImages(pixelflat_norm=pixflat_msc)

    # If this is not a mosaic, we can simply read the pixel flat for the current detector
    else:
        # current detector name
        detname = spectrograph.get_det_name(det)
        # read the file
        with io.fits_open(_pixel_flat_file) as hdu:
            # list of available detectors in the pixel flat file
            file_detnames = [h.name.split('-')[0] for h in hdu] # this list has also the 'PRIMARY' extension
            # check if the current detector is in the list
            if detname in file_detnames:
                # get the index of the current detector
                idx = file_detnames.index(detname)
                # get the pixel flat image
                msgs.info(f'Using pixelflat file: {pixel_flat_file} for {detname}.')
                nrm_image = FlatImages(pixelflat_norm=hdu[idx].data)
            else:
                msgs.error(f'{detname} not found in the pixel flat file: '
                           f'{pixel_flat_file}. Cannot load the pixel flat!')
                nrm_image = None

    return merge(flatimages, nrm_image)

