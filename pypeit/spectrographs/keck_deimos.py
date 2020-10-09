""" Implements DEIMOS-specific functions, including reading in slitmask design files.
"""

import glob
import re
import os
import numpy as np
import warnings

from pkg_resources import resource_filename

from scipy import interpolate
from astropy.io import fits

from pkg_resources import resource_filename

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container

from pypeit.utils import index_of_x_eq_y

from pypeit.spectrographs.slitmask import SlitMask
from pypeit.spectrographs.opticalmodel import ReflectionGrating, OpticalModel, DetectorMap
from IPython import embed

class KeckDEIMOSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/DEIMOS specific code
    """
    ndet = 8

    def __init__(self):
        # Get it started
        super(KeckDEIMOSSpectrograph, self).__init__()
        self.spectrograph = 'keck_deimos'
        self.telescope = telescopes.KeckTelescopePar()
        self.camera = 'DEIMOS'

        # Don't instantiate these until they're needed
        self.grating = None
        self.optical_model = None
        self.detector_map = None
        self.amap = None
        self.bmap = None

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

        # Detector 1
        detector_dict1 = dict(
            binning         = binning,
            det             = 1,
            dataext         = 1,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.1185,
            darkcurr        = 4.19,
            saturation      = 65535., # ADU
            nonlinear       = 0.95,   # Changed by JFH from 0.86 to 0.95
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(1.226),
            ronoise         = np.atleast_1d(2.570),
            )
        # Detector 2
        detector_dict2 = detector_dict1.copy()
        detector_dict2.update(dict(
            det=2,
            dataext=2,
            darkcurr=3.46,
            gain=np.atleast_1d(1.188),
            ronoise=np.atleast_1d(2.491),
        ))
        # Detector 3
        detector_dict3 = detector_dict1.copy()
        detector_dict3.update(dict(
            det=3,
            dataext=3,
            darkcurr=4.03,
            gain=np.atleast_1d(1.248),
            ronoise=np.atleast_1d(2.618),
        ))
        # Detector 4
        detector_dict4 = detector_dict1.copy()
        detector_dict4.update(dict(
            det=4,
            dataext=4,
            darkcurr=3.80,
            gain=np.atleast_1d(1.220),
            ronoise=np.atleast_1d(2.557),
        ))
        # Detector 5
        detector_dict5 = detector_dict1.copy()
        detector_dict5.update(dict(
            det=5,
            dataext=5,
            darkcurr=4.71,
            gain=np.atleast_1d(1.184),
            ronoise=np.atleast_1d(2.482),
        ))
        # Detector 6
        detector_dict6 = detector_dict1.copy()
        detector_dict6.update(dict(
            det=6,
            dataext=6,
            darkcurr=4.28,
            gain=np.atleast_1d(1.177),
            ronoise=np.atleast_1d(2.469),
        ))
        # Detector 7
        detector_dict7 = detector_dict1.copy()
        detector_dict7.update(dict(
            det=7,
            dataext=7,
            darkcurr=3.33,
            gain=np.atleast_1d(1.201),
            ronoise=np.atleast_1d(2.518),
        ))
        # Detector 8
        detector_dict8 = detector_dict1.copy()
        detector_dict8.update(dict(
            det=8,
            dataext=8,
            darkcurr=3.69,
            gain=np.atleast_1d(1.230),
            ronoise=np.atleast_1d(2.580),
        ))
        detectors = [detector_dict1, detector_dict2, detector_dict3, detector_dict4,
                     detector_dict5, detector_dict6, detector_dict7, detector_dict8]
        # Return
        return detector_container.DetectorContainer(**detectors[det-1])


    def default_pypeit_par(self):
        """
        Set default parameters for Keck DEIMOS reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'keck_deimos'

        # Spectral flexure correction
        par['flexure']['spec_method'] = 'boxcar'
        # Set wave tilts order
        par['calibrations']['slitedges']['edge_thresh'] = 50.
        par['calibrations']['slitedges']['fit_order'] = 3
        par['calibrations']['slitedges']['minimum_slit_gap'] = 0.25
        par['calibrations']['slitedges']['minimum_slit_length_sci'] = 4.
#        par['calibrations']['slitedges']['sync_clip'] = False

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['ArI','NeI','KrI','XeI']
        #par['calibrations']['wavelengths']['nonlinear_counts'] \
        #        = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 3
        par['calibrations']['wavelengths']['match_toler'] = 2.5

        # Do not require bias frames
        turn_off = dict(use_biasimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Alter the method used to combine pixel flats
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['pixelflatframe']['process']['sig_lohi'] = [10.,10.]

        # Set the default exposure time ranges for the frame typing
#        par['calibrations']['biasframe']['exprng'] = [None, 2]
#        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
#        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
#        par['calibrations']['pixelflatframe']['exprng'] = [None, 30]
#        par['calibrations']['traceframe']['exprng'] = [None, 30]
#        par['scienceframe']['exprng'] = [30, None]
        
        # LACosmics parameters
        par['scienceframe']['process']['sigclip'] = 4.0
        par['scienceframe']['process']['objlim'] = 1.5

        # If telluric is triggered
        par['sensfunc']['IR']['telgridfile'] = resource_filename('pypeit', '/data/telluric/TelFit_MaunaKea_3100_26100_R20000.fits')

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

        headarr = self.get_headarr(scifile)

        # Turn PCA off for long slits
        # TODO: I'm a bit worried that this won't catch all
        # long-slits...
        if ('Long' in self.get_meta_value(headarr, 'decker')) or (
                'LVMslit' in self.get_meta_value(headarr, 'decker')):
            par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Turn on the use of mask design
        if 'Long' not in self.get_meta_value(headarr, 'decker'):
            par['calibrations']['slitedges']['use_maskdesign'] = True

        # Templates
        if self.get_meta_value(headarr, 'dispname') == '600ZD':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_deimos_600.fits'
            par['calibrations']['wavelengths']['lamps'] += ['CdI', 'ZnI', 'HgI']
        elif self.get_meta_value(headarr, 'dispname') == '830G':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_deimos_830G.fits'
        elif self.get_meta_value(headarr, 'dispname') == '1200G':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_deimos_1200G.fits'

        # FWHM
        binning = parse.parse_binning(self.get_meta_value(headarr, 'binning'))
        par['calibrations']['wavelengths']['fwhm'] = 6.0 / binning[1]

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
        meta['target'] = dict(ext=0, card='TARGNAME')
        meta['decker'] = dict(ext=0, card='SLMSKNAM')
        meta['binning'] = dict(card=None, compound=True)

        meta['mjd'] = dict(ext=0, card='MJD-OBS')
        meta['exptime'] = dict(ext=0, card='ELAPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')
        meta['dispname'] = dict(ext=0, card='GRATENAM')
        # Extras for config and frametyping
        meta['hatch'] = dict(ext=0, card='HATCHPOS')
        meta['dispangle'] = dict(card=None, compound=True, rtol=1e-5)
        # Image type
        meta['idname'] = dict(ext=0, card='OBSTYPE')
        # Lamps
        meta['lampstat01'] = dict(ext=0, card='LAMPS')

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
            binspatial, binspec = parse.parse_binning(headarr[0]['BINNING'])
            binning = parse.binning2string(binspec, binspatial)
            return binning
        elif meta_key == 'dispangle':
            if headarr[0]['GRATEPOS'] == 3:
                return headarr[0]['G3TLTWAV']
            elif headarr[0]['GRATEPOS'] == 4:
                return headarr[0]['G4TLTWAV']
            else:
                msgs.warn('This is probably a problem. Non-standard DEIMOS GRATEPOS={0}.'.format(headarr[0]['GRATEPOS']))
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
        return ['dispname', 'decker', 'binning', 'dispangle']


    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'Object') & (fitstbl['lampstat01'] == 'Off') \
                        & (fitstbl['hatch'] == 'open')
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'Bias') & (fitstbl['lampstat01'] == 'Off') \
                        & (fitstbl['hatch'] == 'closed')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            # Flats and trace frames are typed together
            is_flat = np.any(np.vstack(((fitstbl['idname'] == n) & (fitstbl['hatch'] == h)
                                    for n,h in zip(['IntFlat', 'DmFlat', 'SkyFlat'],
                                                   ['closed', 'open', 'open']))), axis=0)
            return good_exp & is_flat
        if ftype == 'pinhole':
            # Pinhole frames are never assigned for DEIMOS
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype == 'dark':
            return good_exp & (fitstbl['idname'] == 'Dark') & (fitstbl['lampstat01'] == 'Off') \
                        & (fitstbl['hatch'] == 'closed')
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'Line') & (fitstbl['hatch'] == 'closed')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    # TODO: We should aim to get rid of this...
    def idname(self, ftype):
        """
        Return the `idname` for the selected frame type for this instrument.

        Args:
            ftype (str):
                File type, which should be one of the keys in
                :class:`pypeit.core.framematch.FrameTypeBitMask`.

        Returns:
            str: The value of `idname` that should be available in the
            `PypeItMetaData` instance that identifies frames of this
            type.
        """
        # TODO: Fill in the rest of these.
        name = { 'arc': 'Line',
                 'tilt': None,
                 'bias': None,
                 'dark': None,
                 'pinhole': None,
                 'pixelflat': 'IntFlat',
                 'science': 'Object',
                 'standard': None,
                 'trace': 'IntFlat' }
        return name[ftype]

    def get_rawimage(self, raw_file, det):
        """
        Read a raw DEIMOS data frame (one or more detectors).

        Data are unpacked from the multi-extension HDU.  Function is
        based on :func:`pypeit.spectrographs.keck_lris.read_lris`, which
        was based on the IDL procedure ``readmhdufits.pro``.

        Parameters
        ----------
        raw_file : str
            Filename
        det : int or None
            if None, return all 8 detectors!

        Returns
        -------
        tuple
            See :func:`pypeit.spectrograph.spectrograph.get_rawimage`

        """
        # Check for file; allow for extra .gz, etc. suffix
        fil = glob.glob(raw_file + '*')
        if len(fil) != 1:
            msgs.error('Found {0} files matching {1}'.format(len(fil), raw_file + '*'))
        # Read
        try:
            msgs.info("Reading DEIMOS file: {:s}".format(fil[0]))
        except AttributeError:
            print("Reading DEIMOS file: {:s}".format(fil[0]))

        hdu = fits.open(fil[0])
        head0 = hdu[0].header

        # Get post, pre-pix values
        postpix = head0['POSTPIX']
        detlsize = head0['DETLSIZE']
        x0, x_npix, y0, y_npix = np.array(parse.load_sections(detlsize)).flatten()

        # Create final image
        if det is None:
            image = np.zeros((x_npix, y_npix + 4 * postpix))
            rawdatasec_img = np.zeros_like(image, dtype=int)
            oscansec_img = np.zeros_like(image, dtype=int)

        # get the x and y binning factors...
        binning = head0['BINNING']
        if binning != '1,1':
            msgs.error("This binning for DEIMOS might not work.  But it might..")

        # DEIMOS detectors
        nchip = 8

        if det is None:
            chips = range(nchip)
        else:
            chips = [det - 1]  # Indexing starts at 0 here
        # Loop
        for tt in chips:
            data, oscan = deimos_read_1chip(hdu, tt + 1)

            # One detector??
            if det is not None:
                image = np.zeros((data.shape[0], data.shape[1] + oscan.shape[1]))
                rawdatasec_img = np.zeros_like(image, dtype=int)
                oscansec_img = np.zeros_like(image, dtype=int)

            # Indexing
            x1, x2, y1, y2, o_x1, o_x2, o_y1, o_y2 = indexing(tt, postpix, det=det)

            # Fill
            image[y1:y2, x1:x2] = data
            rawdatasec_img[y1:y2, x1:x2] = 1 # Amp
            image[o_y1:o_y2, o_x1:o_x2] = oscan
            oscansec_img[o_y1:o_y2, o_x1:o_x2] = 1 # Amp

        # Return
        exptime = hdu[self.meta['exptime']['ext']].header[self.meta['exptime']['card']]
        return self.get_detector_par(hdu, det if det is not None else 1), \
               image, hdu, exptime, rawdatasec_img, oscansec_img

#    def load_raw_frame(self, raw_file, det=None):
#        """
#        Wrapper to the raw image reader for DEIMOS
#
#        Args:
#            raw_file:  str, filename
#            det: int, REQUIRED
#              Desired detector
#            **null_kwargs:
#              Captured and never used
#
#        Returns:
#            raw_img: ndarray
#              Raw image;  likely unsigned int
#            head0: Header
#        """
#        raw_img, hdu, _ = read_deimos(raw_file, det=det)
#
#        return raw_img, hdu

#    def get_image_section(self, inp=None, det=1, section='datasec'):
#        """
#        Return a string representation of a slice defining a section of
#        the detector image.
#
#        Overwrites base class function
#
#        Args:
#            inp (:obj:`str`, `astropy.io.fits.Header`_, optional):
#                String providing the file name to read, or the relevant
#                header object.  Default is None, meaning that the
#                detector attribute must provide the image section
#                itself, not the header keyword.
#            det (:obj:`int`, optional):
#                1-indexed detector number.
#            section (:obj:`str`, optional):
#                The section to return.  Should be either 'datasec' or
#                'oscansec', according to the
#                :class:`pypeitpar.DetectorPar` keywords.
#
#        Returns:
#            tuple: Returns three objects: (1) A list of string
#            representations for the image sections, one string per
#            amplifier.  The sections are *always* returned in PypeIt
#            order: spectral then spatial.  (2) Boolean indicating if the
#            slices are one indexed.  (3) Boolean indicating if the
#            slices should include the last pixel.  The latter two are
#            always returned as True following the FITS convention.
#        """
#        # Read the file
#        if inp is None:
#            msgs.error('Must provide Keck DEIMOS file or hdulist to get image section.')
#        # Read em
#        shape, datasec, oscansec, _ = deimos_image_sections(inp, det)
#        if section == 'datasec':
#            return datasec, False, False
#        elif section == 'oscansec':
#            return oscansec, False, False
#        else:
#            raise ValueError('Unrecognized keyword: {0}'.format(section))
#
#    def get_raw_image_shape(self, hdulist, det=None, **null_kwargs):
#        """
#        Overrides :class:`Spectrograph.get_image_shape` for LRIS images.
#
#        Must always provide a file.
#        """
#        # Do it
#        self._check_detector()
#        shape, datasec, oscansec, _ = deimos_image_sections(hdulist, det)
#        self.naxis = shape
#        return self.naxis

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Override parent bpm function with BPM specific to DEIMOS.

        .. todo::
            Allow for binning changes.

        Parameters
        ----------
        det : int, REQUIRED
        msbias : numpy.ndarray, required if the user wishes to generate a BPM based on a master bias
        **null_kwargs:
            Captured and never used

        Returns
        -------
        bpix : ndarray
          0 = ok; 1 = Mask

        """
        bpm_img = self.empty_bpm(filename, det, shape=shape)

        # Fill in bad pixels if a master bias frame is provided
        if msbias is not None:
            return self.bpm_frombias(msbias, det, bpm_img)

        if det == 1:
            bpm_img[:,1052:1054] = 1
        elif det == 2:
            bpm_img[:,0:4] = 1
            bpm_img[:,376:381] = 1
            bpm_img[:,489] = 1
            bpm_img[:,1333:1335] = 1
            bpm_img[:,2047] = 1
        elif det == 3:
            bpm_img[:,0:4] = 1
            bpm_img[:,221] = 1
            bpm_img[:,260] = 1
            bpm_img[:,366] = 1
            bpm_img[:,816:819] = 1
            bpm_img[:,851] = 1
            bpm_img[:,940] = 1
            bpm_img[:,1167] = 1
            bpm_img[:,1280] = 1
            bpm_img[:,1301:1303] = 1
            bpm_img[:,1744:1747] = 1
            bpm_img[:,-4:] = 1
        elif det == 4:
            bpm_img[:,0:4] = 1
            bpm_img[:,47] = 1
            bpm_img[:,744] = 1
            bpm_img[:,790:792] = 1
            bpm_img[:,997:999] = 1
        elif det == 5:
            bpm_img[:,25:27] = 1
            bpm_img[:,128:130] = 1
            bpm_img[:,1535:1539] = 1
        elif det == 7:
            bpm_img[:,426:428] = 1
            bpm_img[:,676] = 1
            bpm_img[:,1176:1178] = 1
        elif det == 8:
            bpm_img[:,440] = 1
            bpm_img[:,509:513] = 1
            bpm_img[:,806] = 1
            bpm_img[:,931:934] = 1

        return bpm_img

    def get_slitmask(self, filename):
        """
        Parse the slitmask data from a DEIMOS file into a
        :class:`pypeit.spectrographs.slitmask.SlitMask` object.

        Args:
            filename (:obj:`str`):
                Name of the file to read.
        """
        # Open the file
        hdu = fits.open(filename)

        # Build the object data
        #   - Find the index of the object IDs in the slit-object
        #     mapping that match the object catalog
        mapid = hdu['SlitObjMap'].data['ObjectID']
        catid = hdu['ObjectCat'].data['ObjectID']
        indx = index_of_x_eq_y(mapid, catid)
        #   - Pull out the slit ID, object ID, and object coordinates
        objects = np.array([hdu['SlitObjMap'].data['dSlitId'][indx].astype(float),
                            catid.astype(float),
                            hdu['ObjectCat'].data['RA_OBJ'],
                            hdu['ObjectCat'].data['DEC_OBJ']]).T
        #   - Only keep the objects that are in the slit-object mapping
        objects = objects[mapid[indx] == catid]

        # Match the slit IDs in DesiSlits to those in BluSlits
        indx = index_of_x_eq_y(hdu['DesiSlits'].data['dSlitId'], hdu['BluSlits'].data['dSlitId'],
                               strict=True)

        # Instantiate the slit mask object and return it
        self.slitmask = SlitMask(np.array([hdu['BluSlits'].data['slitX1'],
                                           hdu['BluSlits'].data['slitY1'],
                                           hdu['BluSlits'].data['slitX2'],
                                           hdu['BluSlits'].data['slitY2'],
                                           hdu['BluSlits'].data['slitX3'],
                                           hdu['BluSlits'].data['slitY3'],
                                           hdu['BluSlits'].data['slitX4'],
                                           hdu['BluSlits'].data['slitY4']]).T.reshape(-1,4,2),
                                 slitid=hdu['BluSlits'].data['dSlitId'],
                                 align=hdu['DesiSlits'].data['slitTyp'][indx] == 'A',
                                 science=hdu['DesiSlits'].data['slitTyp'][indx] == 'P',
                                 onsky=np.array([hdu['DesiSlits'].data['slitRA'][indx],
                                                 hdu['DesiSlits'].data['slitDec'][indx],
                                                 hdu['DesiSlits'].data['slitLen'][indx],
                                                 hdu['DesiSlits'].data['slitWid'][indx],
                                                 hdu['DesiSlits'].data['slitLPA'][indx]]).T,
                                 objects=objects)
        return self.slitmask

    def get_grating(self, filename):
        """
        Taken from xidl/DEEP2/spec2d/pro/deimos_omodel.pro and
        xidl/DEEP2/spec2d/pro/deimos_grating.pro
        """
        hdu = fits.open(filename)

        # Grating slider
        slider = hdu[0].header['GRATEPOS']
        # TODO: Add test for slider

        # Central wavelength, grating angle, and tilt position
        if slider == 3:
            central_wave = hdu[0].header['G3TLTWAV']
            # Not used
            #angle = (hdu[0].header['G3TLTRAW'] + 29094)/2500
            tilt = hdu[0].header['G3TLTVAL']
        elif slider in [2,4]:
            # Slider is 2 or 4
            central_wave = hdu[0].header['G4TLTWAV']
            # Not used
            #angle = (hdu[0].header['G4TLTRAW'] + 40934)/2500
            tilt = hdu[0].header['G4TLTVAL']
        else:
            raise ValueError('Slider has unknown value: {0}'.format(slider))

        # Ruling
        name = hdu[0].header['GRATENAM']
        if 'Mirror' in name:
            ruling = 0
        else:
            # Remove all non-numeric characters from the name and
            # convert to a floating point number
            ruling = float(re.sub('[^0-9]', '', name))
            # Adjust
            if abs(ruling-1200) < 0.5:
                ruling = 1200.06
            elif abs(ruling-831) <  2:
                ruling = 831.90

        # Get the orientation of the grating
        roll, yaw, tilt = KeckDEIMOSSpectrograph._grating_orientation(slider, ruling, tilt)

        self.grating = None if ruling == 0 else ReflectionGrating(ruling, tilt, roll, yaw,
                                                                  central_wave=central_wave)
        return self.grating

    def get_detector_map(self):
        if self.detector_map is None:
            self.detector_map = DEIMOSDetectorMap()
        return self.detector_map

    @staticmethod
    def _grating_orientation(slider, ruling, tilt):
        """
        Return the roll, yaw, and tilt of the grating.

        Numbers are hardwired.

        From xidl/DEEP2/spec2d/pro/omodel_params.pro
        """
        if slider == 2 and int(ruling) == 0:
            # Mirror in place of the grating
            return 0., 0., -19.423

        if slider == 2:
            raise ValueError('Ruling should be 0 if slider in position 2.')

        # Use the calibrated coefficients
        # These orientation coefficients are the newest ones and are meant for
        # observations obtained Post-2016 Servicing.
        # TODO: Figure out the impact of these coefficients on the slits identification.
        # We may not need to change them according to when the observations were taken
        _ruling = int(ruling) if int(ruling) in [600, 831, 900, 1200] else 'other'
        orientation_coeffs = {3: {    600: [ 0.145, -0.008, 5.6e-4, -0.146],
                                      831: [ 0.143,  0.000, 5.6e-4, -0.018],
                                      900: [ 0.141,  0.000, 5.6e-4, -0.118],
                                     1200: [ 0.145,  0.055, 5.6e-4, -0.141],
                                  'other': [ 0.145,  0.000, 5.6e-4, -0.141] },
                              4: {    600: [-0.065,  0.063, 6.9e-4, -0.108],
                                      831: [-0.034,  0.060, 6.9e-4, -0.038],
                                      900: [-0.064,  0.083, 6.9e-4, -0.060],
                                     1200: [-0.052,  0.122, 6.9e-4, -0.110],
                                  'other': [-0.050,  0.080, 6.9e-4, -0.110] } }

        # Orientation coefficients meant for observations taken Pre-2016 Servicing
        # orientation_coeffs = {3: {    600: [ 0.145, -0.008, 5.6e-4, -0.182],
        #                               831: [ 0.143,  0.000, 5.6e-4, -0.182],
        #                               900: [ 0.141,  0.000, 5.6e-4, -0.134],
        #                              1200: [ 0.145,  0.055, 5.6e-4, -0.181],
        #                           'other': [ 0.145,  0.000, 5.6e-4, -0.182] },
        #                       4: {    600: [-0.065,  0.063, 6.9e-4, -0.298],
        #                               831: [-0.034,  0.060, 6.9e-4, -0.196],
        #                               900: [-0.064,  0.083, 6.9e-4, -0.277],
        #                              1200: [-0.052,  0.122, 6.9e-4, -0.294],
        #                           'other': [-0.050,  0.080, 6.9e-4, -0.250] } }

        # Return calbirated roll, yaw, and tilt
        return orientation_coeffs[slider][_ruling][0], \
                orientation_coeffs[slider][_ruling][1], \
                tilt*(1-orientation_coeffs[slider][_ruling][2]) \
                    + orientation_coeffs[slider][_ruling][3]



    def get_amapbmap(self, filename):
        """
            Select the pre-grating (amap) and post-grating (bmap) maps according to the slider.

        Args:
            filename (:obj:`str`, optional):
                The filename to read the slider information from the header.

        Returns:
            Two attributes :attr:`amap` and :attr:`bmap`.

        """
        hdu = fits.open(filename)

        # Grating slider
        slider = hdu[0].header['GRATEPOS']

        mp_dir = resource_filename('pypeit', 'data/static_calibs/keck_deimos/')

        if slider in [3,4]:
            self.amap = fits.getdata(mp_dir+'amap.s{}.2003mar04.fits'.format(slider))
            self.bmap = fits.getdata(mp_dir+'bmap.s{}.2003mar04.fits'.format(slider))
        else:
            msgs.error('No amap/bmap available for slider {0}. Set `use_maskdesign = False`'.format(slider))
        #TODO: Figure out which amap and bmap to use for slider 2

        return self.amap, self.bmap



    def mask_to_pixel_coordinates(self, x=None, y=None, wave=None, order=1, filename=None,
                                  corners=False):
        r"""
        Convert the mask coordinates in mm to pixel coordinates on the
        DEIMOS detector.

        If not already instantiated, the :attr:`slitmask`,
        :attr:`grating`, :attr:`optical_model`, and :attr:`detector_map`
        attributes are instantiated.  If these are not instantiated, a
        file must be provided.  If no arguments are provided, the
        function expects these attributes to be set and will output the
        pixel coordinates for the centers of the slits in the
        :attr:`slitmask` at the central wavelength of the
        :attr:`grating`.

        Method generally expected to be executed in one of two modes:
            - Use the `filename` to read the slit mask and determine the
              detector positions at the central wavelength.
            - Specifically map the provided x, y, and wave values to the
              detector.

        If arrays are provided for both `x`, `y`, and `wave`, the
        returned objects have the shape :math:`N_\lambda\times S_x`,
        where :math:`S_x` is the shape of the x and y arrays.

        Args:
            x (array-like, optional):
                The x coordinates in the slit mask in mm.  Default is to
                use the center of the slits in the :attr:`slitmask`.
            y (array-like, optional):
                The y coordinates in the slit mask in mm.  Default is to
                use the center of the slits in the :attr:`slitmask`.
            wave (array-like, optional):
                The wavelengths in angstroms for the propagated
                coordinates.  If not provided, an array of wavelength
                covering the full DEIMOS wavelength range will be used.
            order (:obj:`int`, optional):
                The grating order.  Default is 1.
            filename (:obj:`str`, optional):
                The filename to use to (re)instantiate the
                :attr:`slitmask` and :attr:`grating`.  Default is to use
                previously instantiated attributes.
            corners (:obj:`bool`, optional):
                Instead of using the centers of the slits in the
                :attr:`slitmask`, return the detector pixel coordinates
                for the corners of all slits.

        Returns:
            numpy.ndarray: Returns 5 arrays: (1-2) the x and y
            coordinates in the image plane in mm, (3) the detector
            (1-indexed) where the slit should land at the provided
            wavelength(s), and (4-5) the pixel coordinates (1-indexed)
            in the relevant detector.

        Raises:
            ValueError:
                Raised if the user provides one but not both of the x
                and y coordinates, if no coordinates are provided or
                available within the :attr:`slitmask`, or if the
                :attr:`grating`, :attr:`amap` or :attr:`bmap` haven't been
                defined and not file is provided.
        """
        # Cannot provide just one of x or y
        if x is None and y is not None or x is not None and y is None:
            raise ValueError('Must provide both x and y or neither to use slit mask.')

        # Use the file to update the slitmask (if no x coordinates are
        # provided) and the grating
        if filename is not None:
            if x is None and y is None:
                # Reset the slit mask
                self.get_slitmask(filename)
            # Reset the grating
            self.get_grating(filename)
            # Load pre- and post-grating maps
            self.get_amapbmap(filename)

        if self.amap is None and self.bmap is None:
            raise ValueError('Must select amap and bmap; provide a file or use get_amapbmap()')

        # Check that any coordinates are available
        if x is None and y is None and self.slitmask is None:
            raise ValueError('No coordinates; Provide them directly or instantiate slit mask.')

        # Make sure the coordinates are numpy arrays
        _x = None if x is None else np.atleast_1d(x)
        _y = None if y is None else np.atleast_1d(y)
        if _x is None:
            # Use all the slit centers or corners
            _x = self.slitmask.corners[...,0].ravel() if corners else self.slitmask.center[:,0]
            _y = self.slitmask.corners[...,1].ravel() if corners else self.slitmask.center[:,1]

        # Check that the grating is defined
        if self.grating is None:
            raise ValueError('Must define a grating first; provide a file or use get_grating()')

        # Instantiate the optical model or reset it grating
        if self.optical_model is None:
            self.optical_model = DEIMOSOpticalModel(self.grating)
        else:
            self.optical_model.reset_grating(self.grating)

        # Instantiate the detector map, if necessary
        self.get_detector_map()

        # hard-coded for DEIMOS: wavelength array if wave is None
        if wave is None:
            npoints = 250
            wave = np.arange(npoints) * 24. + 4000.

        # Compute the detector image plane coordinates (in pixels)
        x_img, y_img = self.optical_model.mask_to_imaging_coordinates(_x, _y, self.amap, self.bmap,
                                                                      nslits=self.slitmask.nslits,
                                                                      wave=wave, order=order)
        # Reshape if computing the corner positions
        if corners:
            x_img = x_img.reshape(self.slitmask.corners.shape[:2])
            y_img = y_img.reshape(self.slitmask.corners.shape[:2])

        # Use the detector map to convert to the detector coordinates
        return (x_img, y_img) + self.detector_map.ccd_coordinates(x_img, y_img, in_mm=False)


class DEIMOSOpticalModel(OpticalModel):
    # TODO: Are focal_r_surface (!R_IMSURF) and focal_r_curvature
    # (!R_CURV) supposed to be the same?  If so, consolodate these into
    # a single number.
    def __init__(self, grating):
        super(DEIMOSOpticalModel, self).__init__(
                    20018.4,                # Pupil distance in mm (!PPLDIST, !D_1)
                    2133.6,                 # Radius of the image surface in mm (!R_IMSURF)
                    2124.71,                # Focal-plane radius of curvature in mm (!R_CURV)
                    2120.9,                 # Mask radius of curvature in mm (!M_RCURV)
                    np.radians(6.),         # Mask tilt angle in radians (!M_ANGLE)
                    128.803,                # Mask y zero point in mm (!ZPT_YM)
                    3.378,                  # Mask z zero-point in mm (!MASK_HT0)
                    2197.1,                 # Collimator distance in mm (sys.COL_DST)
                    4394.2,                 # Collimator radius of curvature in mm (!R_COLL)
                    -0.75,                  # Collimator curvature constant (!K_COLL)
                    np.radians(0.002),      # Collimator tilt error in radians (sys.COL_ERR)
                    0.0,                    # Collimator tilt phi angle in radians (sys.COL_PHI)
                    grating,                # DEIMOS grating object
                    np.radians(2.752),      # Camera angle in radians (sys.CAM_ANG)
                    np.pi/2,                # Camera tilt phi angle in radians (sys.CAM_PHI)
                    382.0,                  # Camera focal length in mm (sys.CAM_FOC)
                    DEIMOSCameraDistortion(),   # Object used to apply/remove camera distortions
                    np.radians(0.021),      # ICS rotation in radians (sys.MOS_ROT)
                    [-0.234, -3.822])       # Camera optical axis center in mm (sys.X_OPT,sys.Y_OPT)

        # Include tent mirror
        self.tent_theta = np.radians(71.5-0.5)  # Tent mirror theta angle (sys.TNT_ANG)
        self.tent_phi = np.radians(90.+0.081)   # Tent mirror phi angle (sys.TNT_PHI)

        #TENT MIRROR: this mirror is OK to leave in del-theta,phi
        self.tent_reflection \
                = OpticalModel.get_reflection_transform(self.tent_theta, self.tent_phi)

    def reset_grating(self, grating):
        self.grating = grating

    def mask_coo_to_grating_input_vectors(self, x, y):
        """
        Propagate rays from the mask plane to the grating.

        Taken from xidl/DEEP2/spec2d/pro/model/pre_grating.pro

        Need to override parent class to add tent mirror reflection.
        """
        r = super(DEIMOSOpticalModel, self).mask_coo_to_grating_input_vectors(x, y)
        # Reflect off the tent mirror and return
        return OpticalModel.reflect(r, self.tent_reflection)


class DEIMOSCameraDistortion:
    """Class to remove or apply DEIMOS camera distortion."""
    def __init__(self):
        self.c0 = 1.
        self.c2 = 0.0457563
        self.c4 = -0.3088123
        self.c6 = -14.917
    
        x = np.linspace(-0.6, 0.6, 1000)
        y = self.remove_distortion(x)
        self.interpolator = interpolate.interp1d(y, x)

    def remove_distortion(self, x):
        x2 = np.square(x)
        return x / (self.c0 + x2 * (self.c2 + x2 * (self.c4 + x2 * self.c6)))

    def apply_distortion(self, y):
        indx = (y > self.interpolator.x[0]) & (y < self.interpolator.x[-1])
        if not np.all(indx):
            warnings.warn('Some input angles outside of valid distortion interval!')
        x = np.zeros_like(y)
        x[indx] = self.interpolator(y[indx])
        return x


class DEIMOSDetectorMap(DetectorMap):
    """
    A map of the center coordinates and rotation of each CCD in DEIMOS.

    !! PIXEL COORDINATES ARE 1-INDEXED !!
    """
    def __init__(self):
        # Number of chips
        self.nccd = 8

        # Number of pixels for each chip in each dimension
        self.npix = np.array([2048, 4096])

        # The size of the CCD pixels in mm
        self.pixel_size = 0.015

        # Nominal gap between each CCD in each dimension in mm
        self.ccd_gap = np.array([1, 0.1])

        # Width of the CCD edge in each dimension in mm
        self.ccd_edge = np.array([0.154, 0.070])

        # Effective size of each chip in each dimension in pixels
        self.ccd_size = self.npix + (2*self.ccd_edge + self.ccd_gap)/self.pixel_size

        # Center coordinates
        origin = np.array([[-1.5,-0.5], [-0.5,-0.5], [ 0.5,-0.5], [ 1.5,-0.5],
                           [-1.5, 0.5], [-0.5, 0.5], [ 0.5, 0.5], [ 1.5, 0.5]])
        offset = np.array([[-20.05, 14.12], [-12.64, 7.25], [0.00, 0.00], [-1.34, -19.92],
                           [-19.02, 16.46], [ -9.65, 8.95], [1.88, 1.02], [ 4.81, -24.01]])
        self.ccd_center = origin * self.ccd_size[None,:] + offset
        
        # Construct the rotation matrix
        self.rotation = np.radians([-0.082, 0.030, 0.0, -0.1206, 0.136, -0.06, -0.019, -0.082])
        cosa = np.cos(self.rotation)
        sina = np.sin(self.rotation)
        self.rot_matrix = np.array([cosa, -sina, sina, cosa]).T.reshape(self.nccd,2,2)

        # ccd_geom.pro has offsets by sys.CN_XERR, but these are all 0.

#def deimos_image_sections(inp, det):
#    """
#    Parse the image for the raw image shape and data sections
#
#    Args:
#        inp (str or `astropy.io.fits.HDUList`_ object):
#        det (int):
#
#    Returns:
#        tuple:
#            shape, dsec, osec, ext_items
#            ext_items is a large tuple of bits and pieces for other methods
#                ext_items = hdu, chips, postpix, image
#    """
#    # Check for file; allow for extra .gz, etc. suffix
#    if isinstance(inp, str):
#        fil = glob.glob(inp + '*')
#        if len(fil) != 1:
#            msgs.error('Found {0} files matching {1}'.format(len(fil), inp + '*'))
#        # Read
#        try:
#            msgs.info("Reading DEIMOS file: {:s}".format(fil[0]))
#        except AttributeError:
#            print("Reading DEIMOS file: {:s}".format(fil[0]))
#        # Open
#        hdu = fits.open(fil[0])
#    else:
#        hdu = inp
#    head0 = hdu[0].header
#
#    # Get post, pre-pix values
#    precol = head0['PRECOL']
#    postpix = head0['POSTPIX']
#    preline = head0['PRELINE']
#    postline = head0['POSTLINE']
#    detlsize = head0['DETLSIZE']
#    x0, x_npix, y0, y_npix = np.array(parse.load_sections(detlsize)).flatten()
#
#
#    # Setup for datasec, oscansec
#    dsec = []
#    osec = []
#
#    # get the x and y binning factors...
#    binning = head0['BINNING']
#    if binning != '1,1':
#        msgs.error("This binning for DEIMOS might not work.  But it might..")
#
#    xbin, ybin = [int(ibin) for ibin in binning.split(',')]
#
#    # DEIMOS detectors
#    nchip = 8
#    if det is None:
#        chips = range(nchip)
#    else:
#        chips = [det-1] # Indexing starts at 0 here
#
#    for tt in chips:
#        x1, x2, y1, y2, o_x1, o_x2, o_y1, o_y2 = indexing(tt, postpix, det=det)
#        # Sections
#        idsec = '[{:d}:{:d},{:d}:{:d}]'.format(y1, y2, x1, x2)
#        iosec = '[{:d}:{:d},{:d}:{:d}]'.format(o_y1, o_y2, o_x1, o_x2)
#        dsec.append(idsec)
#        osec.append(iosec)
#
#    # Create final image (if the full image is requested)
#    if det is None:
#        image = np.zeros((x_npix,y_npix+4*postpix))
#        shape = image.shape
#    else:
#        image = None
#        head = hdu[chips[0]+1].header
#        shape = (head['NAXIS2'], head['NAXIS1']-precol)  # We don't load up the precol
#
#    # Pack up a few items for use elsewhere
#    ext_items = hdu, chips, postpix, image
#    # Return
#    return shape, dsec, osec, ext_items

def indexing(itt, postpix, det=None):
    """
    Some annoying book-keeping for instrument placement.

    Parameters
    ----------
    itt : int
    postpix : int
    det : int, optional

    Returns
    -------

    """
    # Deal with single chip
    if det is not None:
        tt = 0
    else:
        tt = itt
    ii = 2048
    jj = 4096
    # y indices
    if tt < 4:
        y1, y2 = 0, jj
    else:
        y1, y2 = jj, 2*jj
    o_y1, o_y2 = y1, y2

    # x
    x1, x2 = (tt%4)*ii, (tt%4 + 1)*ii
    if det is None:
        o_x1 = 4*ii + (tt%4)*postpix
    else:
        o_x1 = ii + (tt%4)*postpix
    o_x2 = o_x1 + postpix

    # Return
    return x1, x2, y1, y2, o_x1, o_x2, o_y1, o_y2

def deimos_read_1chip(hdu,chipno):
    """ Read one of the DEIMOS detectors

    Args:
        hdu (astropy.io.fits.HDUList):
        chipno (int):

    Returns:
        np.ndarray, np.ndarray:
            data, oscan
    """

    # Extract datasec from header
    datsec = hdu[chipno].header['DATASEC']
    detsec = hdu[chipno].header['DETSEC']
    postpix = hdu[0].header['POSTPIX']
    precol = hdu[0].header['PRECOL']

    x1_dat, x2_dat, y1_dat, y2_dat = np.array(parse.load_sections(datsec)).flatten()
    x1_det, x2_det, y1_det, y2_det = np.array(parse.load_sections(detsec)).flatten()

    # This rotates the image to be increasing wavelength to the top
    #data = np.rot90((hdu[chipno].data).T, k=2)
    #nx=data.shape[0]
    #ny=data.shape[1]


    # Science data
    fullimage = hdu[chipno].data
    data = fullimage[x1_dat:x2_dat,y1_dat:y2_dat]

    # Overscan
    oscan = fullimage[:,y2_dat:]

    # Flip as needed
    if x1_det > x2_det:
        data = np.flipud(data)
        oscan = np.flipud(oscan)
    if y1_det > y2_det:
        data = np.fliplr(data)
        oscan = np.fliplr(oscan)

    # Return
    return data, oscan


