'''
Implements DEIMOS-specific functions, including reading in slitmask design files.
'''
from __future__ import absolute_import, division, print_function

import glob
import re
import os
import numpy as np

from scipy import interpolate

from astropy.io import fits

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph

from pypeit.spectrographs.slitmask import SlitMask
from pypeit.spectrographs.opticalmodel import ReflectionGrating, OpticalModel, DetectorMap

from pypeit import debugger

class KeckDEIMOSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/DEIMOS specific code
    """
    def __init__(self):
        # Get it started
        super(KeckDEIMOSSpectrograph, self).__init__()
        self.spectrograph = 'keck_deimos'
        self.telescope = telescopes.KeckTelescopePar()
        self.camera = 'DEIMOS'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 1,
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.1185,
                            darkcurr        = 4.19,
                            saturation      = 65535.,
                            nonlinear       = 0.86,
                            numamplifiers   = 1,
                            gain            = 1.226,
                            ronoise         = 2.570,
                            datasec         = '',       # These are provided by read_deimos
                            oscansec        = '',
                            suffix          = '_01'
                            ),
                # Detector 2
                pypeitpar.DetectorPar(
                            dataext         = 2,
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.1185,
                            darkcurr        = 3.46,
                            saturation      = 65535.,
                            nonlinear       = 0.86,
                            numamplifiers   = 1,
                            gain            = 1.188,
                            ronoise         = 2.491,
                            datasec         = '',       # These are provided by read_deimos
                            oscansec        = '',
                            suffix          = '_02'
                            ),
                # Detector 3
                pypeitpar.DetectorPar(
                            dataext         = 3,
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.1185,
                            darkcurr        = 4.03,
                            saturation      = 65535.,
                            nonlinear       = 0.86,
                            numamplifiers   = 1,
                            gain            = 1.248,
                            ronoise         = 2.618,
                            datasec         = '',       # These are provided by read_deimos
                            oscansec        = '',
                            suffix          = '_03'
                            ),
                # Detector 4
                pypeitpar.DetectorPar(
                            dataext         = 4,
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.1185,
                            darkcurr        = 3.80,
                            saturation      = 65535.,
                            nonlinear       = 0.86,
                            numamplifiers   = 1,
                            gain            = 1.220,
                            ronoise         = 2.557,
                            datasec         = '',       # These are provided by read_deimos
                            oscansec        = '',
                            suffix          = '_04'
                            ),
                # Detector 5
                pypeitpar.DetectorPar(
                            dataext         = 5,
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.1185,
                            darkcurr        = 4.71,
                            saturation      = 65535.,
                            nonlinear       = 0.86,
                            numamplifiers   = 1,
                            gain            = 1.184,
                            ronoise         = 2.482,
                            datasec         = '',       # These are provided by read_deimos
                            oscansec        = '',
                            suffix          = '_05'
                            ),
                # Detector 6
                pypeitpar.DetectorPar(
                            dataext         = 6,
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.1185,
                            darkcurr        = 4.28,
                            saturation      = 65535.,
                            nonlinear       = 0.86,
                            numamplifiers   = 1,
                            gain            = 1.177,
                            ronoise         = 2.469,
                            datasec         = '',       # These are provided by read_deimos
                            oscansec        = '',
                            suffix          = '_06'
                            ),
                # Detector 7
                pypeitpar.DetectorPar(
                            dataext         = 7,
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.1185,
                            darkcurr        = 3.33,
                            saturation      = 65535.,
                            nonlinear       = 0.86,
                            numamplifiers   = 1,
                            gain            = 1.201,
                            ronoise         = 2.518,
                            datasec         = '',       # These are provided by read_deimos
                            oscansec        = '',
                            suffix          = '_07'),
                # Detector 8
                pypeitpar.DetectorPar(
                            dataext         = 8,
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1., 
                            platescale      = 0.1185,
                            darkcurr        = 3.69,
                            saturation      = 65535.,
                            nonlinear       = 0.86,
                            numamplifiers   = 1,
                            gain            = 1.230,
                            ronoise         = 2.580,
                            datasec         = '',       # These are provided by read_deimos
                            oscansec        = '',
                            suffix          = '_08'
                            )]
        self.numhead = 9
        # Uses default timeunit
        # Uses default primary_hdrext
        # self.sky_file ?

        # Don't instantiate these until they're needed
        self.grating = None
        self.optical_model = None
        self.detector_map = None

    def default_pypeit_par(self):
        """
        Set default parameters for Keck LRISb reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'keck_deimos'
        par['flexure']['method'] = 'boxcar'
        # Set wave tilts order
        par['calibrations']['slits']['sigdetect'] = 50.
        par['calibrations']['slits']['trace_npoly'] = 3

        # Overscan subtract the images
        par['calibrations']['biasframe']['useframe'] = 'overscan'

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['ArI','NeI','KrI','XeI']
        par['calibrations']['wavelengths']['nonlinear_counts'] \
                = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 3
        par['calibrations']['wavelengths']['match_toler'] = 2.5

        # Alter the method used to combine pixel flats
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['pixelflatframe']['process']['sig_lohi'] = [10.,10.]

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 2]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['pixelflatframe']['exprng'] = [None, 30]
        par['calibrations']['traceframe']['exprng'] = [None, 30]
        par['scienceframe']['exprng'] = [30, None]
        
        # LACosmics parameters
        par['scienceframe']['process']['sigclip'] = 4.0
        par['scienceframe']['process']['objlim'] = 1.5

        return par

    def config_specific_par(self, par, scifile):

        # Templates
        if self.get_meta_value(scifile, 'dispname') == '600ZD':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_deimos_600.fits'
            par['calibrations']['wavelengths']['lamps'] += ['CdI', 'ZnI', 'HgI']
        elif self.get_meta_value(scifile, 'dispname') == '830G':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_deimos_830G.fits'
        elif self.get_meta_value(scifile, 'dispname') == '1200G':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_deimos_1200G.fits'

        # FWHM
        binning = parse.parse_binning(self.get_meta_value(scifile, 'binning'))
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
                debugger.set_trace()
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
            return good_exp & (fitstbl['lampstat01'] == 'Off') & (fitstbl['hatch'] == 'open')
        if ftype == 'bias':
            return good_exp & (fitstbl['lampstat01'] == 'Off') & (fitstbl['hatch'] == 'closed')
        if ftype == 'pixelflat' or ftype == 'trace':
            # Flats and trace frames are typed together
            return good_exp & (fitstbl['lampstat01'] == 'Qz') & (fitstbl['hatch'] == 'closed')
        if ftype == 'pinhole' or ftype == 'dark':
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype == 'arc':
            return good_exp & (fitstbl['lampstat01'] == 'Kr Xe Ar Ne') & (fitstbl['hatch'] == 'closed')

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
                 'bias': None,
                 'dark': None,
                 'pinhole': None,
                 'pixelflat': 'IntFlat',
                 'science': 'Object',
                 'standard': None,
                 'trace': 'IntFlat' }
        return name[ftype]
  
    def load_raw_img_head(self, raw_file, det=None, **null_kwargs):
        """
        Wrapper to the raw image reader for DEIMOS

        Args:
            raw_file:  str, filename
            det: int, REQUIRED
              Desired detector
            **null_kwargs:
              Captured and never used

        Returns:
            raw_img: ndarray
              Raw image;  likely unsigned int
            head0: Header

        """
        raw_img, head0, _ = read_deimos(raw_file, det=det)

        return raw_img, head0

#    def get_match_criteria(self):
#        match_criteria = {}
#        for key in framematch.FrameTypeBitMask().keys():
#            match_criteria[key] = {}
#        # Standard
#        # Can be over-ruled by flux calibrate = False
#        match_criteria['standard']['match'] = {}
#        match_criteria['standard']['match']['decker'] = ''
#        match_criteria['standard']['match']['binning'] = ''
#        match_criteria['standard']['match']['filter1'] = ''
#        # Bias
#        match_criteria['bias']['match'] = {}
#        match_criteria['bias']['match']['binning'] = ''
#        # Pixelflat
#        match_criteria['pixelflat']['match'] = match_criteria['standard']['match'].copy()
#        # Traceflat
#        match_criteria['trace']['match'] = match_criteria['standard']['match'].copy()
#        # Arc
#        match_criteria['arc']['match'] = match_criteria['standard']['match'].copy()
#        # Return
#        return match_criteria

    def get_image_section(self, inp=None, det=1, section='datasec'):
        """
        Return a string representation of a slice defining a section of
        the detector image.

        Overwrites base class function to use :func:`read_deimos` to get
        the image sections.

        .. todo ::
            - It is really ineffiecient.  Can we parse
              :func:`read_deimos` into something that can give you the
              image section directly?

        This is done separately for the data section and the overscan
        section in case one is defined as a header keyword and the other
        is defined directly.
        
        Args:
            inp (:obj:`str`):
                String providing the file name to read.  Unlike the base
                class, a file name *must* be provided.
            det (:obj:`int`, optional):
                1-indexed detector number.
            section (:obj:`str`, optional):
                The section to return.  Should be either datasec or
                oscansec, according to the :class:`DetectorPar`
                keywords.

        Returns:
            list, bool: A list of string representations for the image
            sections, one string per amplifier, followed by three
            booleans: if the slices are one indexed, if the slices
            should include the last pixel, and if the slice should have
            their order transposed.
        """
        # Read the file
        if inp is None:
            msgs.error('Must provide Keck DEIMOS file to get image section.')
        elif not os.path.isfile(inp):
            msgs.error('File {0} does not exist!'.format(inp))
        temp, head0, secs = read_deimos(inp, det)
        if section == 'datasec':
            return secs[0], False, False, False
        elif section == 'oscansec':
            return secs[1], False, False, False
        else:
            raise ValueError('Unrecognized keyword: {0}'.format(section))

    def get_datasec_img(self, filename, det=1, force=True):
        """
        Create an image identifying the amplifier used to read each pixel.

        Args:
            filename (str):
                Name of the file from which to read the image size.
            det (:obj:`int`, optional):
                Detector number (1-indexed)
            force (:obj:`bool`, optional):
                Force the image to be remade

        Returns:
            `numpy.ndarray`: Integer array identifying the amplifier
            used to read each pixel.
        """
        if self.datasec_img is None or force:
            # Check the detector is defined
            self._check_detector()
            # Get the image shape
            raw_naxis = self.get_raw_image_shape(filename, det=det)

            # Binning is not required because read_deimos accounts for it
#            binning = self.get_meta_value(filename, 'binning')

            data_sections, one_indexed, include_end, transpose \
                    = self.get_image_section(filename, det, section='datasec')

            # Initialize the image (0 means no amplifier)
            self.datasec_img = np.zeros(raw_naxis, dtype=int)
            for i in range(self.detector[det-1]['numamplifiers']):
                # Convert the data section from a string to a slice
                datasec = parse.sec2slice(data_sections[i], one_indexed=one_indexed,
                                          include_end=include_end, require_dim=2,
                                          transpose=transpose) #, binning=binning)
                # Assign the amplifier
                self.datasec_img[datasec] = i+1
        return self.datasec_img

    # WARNING: Uses Spectrograph default get_image_shape.  If no file
    # provided it will fail.  Provide a function like in keck_lris.py
    # that forces a file to be provided?

    def bpm(self, shape=None, filename=None, det=None, **null_kwargs):
        """
        Override parent bpm function with BPM specific to DEIMOS.

        .. todo::
            Allow for binning changes.

        Parameters
        ----------
        det : int, REQUIRED
        **null_kwargs:
            Captured and never used

        Returns
        -------
        bpix : ndarray
          0 = ok; 1 = Mask

        """
        self.empty_bpm(shape=shape, filename=filename, det=det)
        if det == 1:
            self.bpm_img[:,1052:1054] = 1
        elif det == 2:
            self.bpm_img[:,0:4] = 1
            self.bpm_img[:,376:381] = 1
            self.bpm_img[:,489] = 1
            self.bpm_img[:,1333:1335] = 1
            self.bpm_img[:,2047] = 1
        elif det == 3:
            self.bpm_img[:,0:4] = 1
            self.bpm_img[:,221] = 1
            self.bpm_img[:,260] = 1
            self.bpm_img[:,366] = 1
            self.bpm_img[:,816:819] = 1
            self.bpm_img[:,851] = 1
            self.bpm_img[:,940] = 1
            self.bpm_img[:,1167] = 1
            self.bpm_img[:,1280] = 1
            self.bpm_img[:,1301:1303] = 1
            self.bpm_img[:,1744:1747] = 1
            self.bpm_img[:,-4:] = 1
        elif det == 4:
            self.bpm_img[:,0:4] = 1
            self.bpm_img[:,47] = 1
            self.bpm_img[:,744] = 1
            self.bpm_img[:,790:792] = 1
            self.bpm_img[:,997:999] = 1
        elif det == 5:
            self.bpm_img[:,25:27] = 1
            self.bpm_img[:,128:130] = 1
            self.bpm_img[:,1535:1539] = 1
        elif det == 7:
            self.bpm_img[:,426:428] = 1
            self.bpm_img[:,676] = 1
            self.bpm_img[:,1176:1178] = 1
        elif det == 8:
            self.bpm_img[:,440] = 1
            self.bpm_img[:,509:513] = 1
            self.bpm_img[:,806] = 1
            self.bpm_img[:,931:934] = 1

        return self.bpm_img

#    def setup_arcparam(self, arcparam, disperser=None, fitstbl=None, arc_idx=None,
#                       msarc_shape=None, **null_kwargs):
#        """
#
#        Args:
#            arcparam:
#            disperser:
#            fitstbl:
#            arc_idx:
#            msarc_shape:
#            binspectral:
#            **null_kwargs:
#
#        Returns:
#
#        """
#        arcparam['wv_cen'] = fitstbl['dispangle'][arc_idx]
#        # TODO -- Should set according to the lamps that were on
#        #arcparam['lamps'] = ['ArI','NeI','KrI','XeI']
#        # JFH Right now these are all hard wired to use det =1 numbers. Otherwise we will need a separate arcparam for each
#        # detector and there is no mechanism in place to create that yet
#
#        arcparam['nonlinear_counts'] = self.detector[0]['nonlinear']*self.detector[0]['saturation']
##        arcparam['min_nsig'] = 30.  # Minimum signififance
#        arcparam['sigdetect'] = 10.0      # Min significance for arc lines to be used
#        arcparam['wvmnx'] = [3000., 11000.]  # Guess at wavelength range
#        # These parameters influence how the fts are done by pypeit.core.wavecal.fitting.iterative_fitting
#        arcparam['match_toler'] = 3  # Matcing tolerance (pixels)
#        arcparam['func'] = 'legendre'  # Function for fitting
#        arcparam['n_first'] = 2  # Order of polynomial for first fit
#        arcparam['n_final'] = 4  # Order of polynomial for final fit
#        arcparam['nsig_rej'] = 2  # Number of sigma for rejection
#        arcparam['nsig_rej_final'] = 3.0  # Number of sigma for rejection (final fit)
#
#        arcparam['min_ampl'] = 1000.  # Lines tend to be very strong
#        arcparam['wvmnx'][0] = 4000.
#        arcparam['wvmnx'][1] = 11000.
#
#    #        if disperser == '830G': # Blaze 8640
#            arcparam['n_first']=2 # Too much curvature for 1st order
#            arcparam['disp']=0.47 # Ang per pixel (unbinned)
#            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0]
#            arcparam['wvmnx'][0] = 550.
#            arcparam['wvmnx'][1] = 11000.
#            arcparam['min_ampl'] = 3000.  # Lines tend to be very strong
#        elif disperser == '1200G': # Blaze 7760
#            arcparam['n_first']=2 # Too much curvature for 1st order
#            arcparam['disp']=0.32 # Ang per pixel (unbinned)
#            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0]
#            arcparam['wvmnx'][0] = 550.
#            arcparam['wvmnx'][1] = 11000.
#            arcparam['min_ampl'] = 2000.  # Lines tend to be very strong
#        else:
#            msgs.error('Not ready for this disperser {:s}!'.format(disperser))

    def get_slitmask(self, filename):
        hdu = fits.open(filename)
        corners = np.array([hdu['BluSlits'].data['slitX1'],
                            hdu['BluSlits'].data['slitY1'],
                            hdu['BluSlits'].data['slitX2'],
                            hdu['BluSlits'].data['slitY2'],
                            hdu['BluSlits'].data['slitX3'],
                            hdu['BluSlits'].data['slitY3'],
                            hdu['BluSlits'].data['slitX4'],
                            hdu['BluSlits'].data['slitY4']]).T.reshape(-1,4,2)
        self.slitmask = SlitMask(corners, slitid=hdu['BluSlits'].data['dSlitId'])
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
        _ruling = int(ruling) if int(ruling) in [600, 831, 900, 1200] else 'other'
        orientation_coeffs = {3: {    600: [ 0.145, -0.008, 5.6e-4, -0.182],
                                      831: [ 0.143,  0.000, 5.6e-4, -0.182],
                                      900: [ 0.141,  0.000, 5.6e-4, -0.134],
                                     1200: [ 0.145,  0.055, 5.6e-4, -0.181],
                                  'other': [ 0.145,  0.000, 5.6e-4, -0.182] },
                              4: {    600: [-0.065,  0.063, 6.9e-4, -0.298],
                                      831: [-0.034,  0.060, 6.9e-4, -0.196],
                                      900: [-0.064,  0.083, 6.9e-4, -0.277],
                                     1200: [-0.052,  0.122, 6.9e-4, -0.294],
                                  'other': [-0.050,  0.080, 6.9e-4, -0.250] } }

        # Return calbirated roll, yaw, and tilt
        return orientation_coeffs[slider][_ruling][0], \
                orientation_coeffs[slider][_ruling][1], \
                tilt*(1-orientation_coeffs[slider][_ruling][2]) \
                    + orientation_coeffs[slider][_ruling][3]

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
                coordinates.  Default is to use the central wavelength
                of the :attr:`grating`.
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
                :attr:`grating` hasn't been defined and not file is
                provided.
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

        # Compute the detector image plane coordinates (mm)
        x_img, y_img = self.optical_model.mask_to_imaging_coordinates(_x, _y, wave=wave,
                                                                      order=order)
        # Reshape if computing the corner positions
        if corners:
            x_img = x_img.reshape(self.slitmask.corners.shape[:2])
            y_img = y_img.reshape(self.slitmask.corners.shape[:2])

        # Use the detector map to convert to the detector coordinates
        return (x_img, y_img) + self.detector_map.ccd_coordinates(x_img, y_img)


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
   

def read_deimos(raw_file, det=None):
    """
    Read a raw DEIMOS data frame (one or more detectors)
    Packed in a multi-extension HDU
    Based on pypeit.arlris.read_lris...
       Based on readmhdufits.pro

    Parameters
    ----------
    raw_file : str
      Filename

    Returns
    -------
    array : ndarray
      Combined image
    header : FITS header
    sections : tuple
      List of datasec, oscansec sections
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
    precol = head0['PRECOL']
    postpix = head0['POSTPIX']
    preline = head0['PRELINE']
    postline = head0['POSTLINE']
    detlsize = head0['DETLSIZE']
    x0, x_npix, y0, y_npix = np.array(parse.load_sections(detlsize)).flatten()

    # Create final image
    if det is None:
        image = np.zeros((x_npix,y_npix+4*postpix))

    # Setup for datasec, oscansec
    dsec = []
    osec = []

    # get the x and y binning factors...
    binning = head0['BINNING']
    if binning != '1,1':
        msgs.error("This binning for DEIMOS might not work.  But it might..")

    xbin, ybin = [int(ibin) for ibin in binning.split(',')]

    # DEIMOS detectors
    nchip = 8


    if det is None:
        chips = range(nchip)
    else:
        chips = [det-1] # Indexing starts at 0 here
    # Loop
    for tt in chips:
        data, oscan = deimos_read_1chip(hdu, tt+1)


        #if n_elements(nobias) eq 0 then nobias = 0


        # One detector??
        if det is not None:
            image = np.zeros((data.shape[0],data.shape[1]+oscan.shape[1]))

        # Indexing
        x1, x2, y1, y2, o_x1, o_x2, o_y1, o_y2 = indexing(tt, postpix, det=det)

        # Fill
        image[y1:y2, x1:x2] = data
        image[o_y1:o_y2, o_x1:o_x2] = oscan

        # Sections
        idsec = '[{:d}:{:d},{:d}:{:d}]'.format(y1, y2, x1, x2)
        iosec = '[{:d}:{:d},{:d}:{:d}]'.format(o_y1, o_y2, o_x1, o_x2)
        dsec.append(idsec)
        osec.append(iosec)
    # Return
    return image, head0, (dsec,osec)


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

    Parameters
    ----------
    hdu : HDUList
    chipno : int

    Returns
    -------
    data : ndarray
    oscan : ndarray
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


