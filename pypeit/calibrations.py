""" Class for guiding calibration object generation in PypeIt
"""
from __future__ import absolute_import, division, print_function

import os
import numpy as np

from abc import ABCMeta

from astropy.table import Table

from pypeit import msgs
from pypeit.core import pixels
from pypeit import masterframe
from pypeit import arcimage
from pypeit import biasframe
from pypeit import bpmimage
from pypeit import flatfield
from pypeit import traceimage
from pypeit import traceslits
from pypeit import wavecalib
from pypeit import wavetilts
from pypeit import waveimage

from pypeit.metadata import PypeItMetaData

from pypeit.core import procimg
from pypeit.core import parse

from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph

from pypeit import debugger


class Calibrations(object):
    """
    This class is primarily designed to guide the generation of calibration images
    and objects in PypeIt

    Must be able to instantiate spectrograph.

    .. todo::
        Improve docstring...

    Parameters
    ----------
    redux_path : str, optional
      Path to location of PypeIt file (and reductions)
      Defaults to pwd

    Attributes
    ----------

    Inherited Attributes
    --------------------
    """
    __metaclass__ = ABCMeta

    # TODO: master_root is a bit of a kludge.  It could be defined
    # earlier and/or in par.
    def __init__(self, fitstbl, spectrograph=None, par=None, redux_path=None,
                 save_masters=True, write_qa=True, show = False):

        # Check the type of the provided fits table
        if not isinstance(fitstbl, PypeItMetaData):
            msgs.error('fitstbl must be an PypeItMetaData object')

        # Parameters unique to this Object
        self.fitstbl = fitstbl
        self.save_masters = save_masters
        self.write_qa = write_qa
        self.show = show

        # Spectrometer class
        _spectrograph = spectrograph
        if spectrograph is None:
            # Set spectrograph from FITS table instrument header
            # keyword.
            if par is not None and par['rdx']['spectrograph'] != fitstbl['instrume'][0]:
                msgs.error('Specified spectrograph does not match instrument in the fits table!')
            _spectrograph = fitstbl['instrume'][0]
        self.spectrograph = load_spectrograph(_spectrograph)

        # Instantiate the parameters
        # TODO: How far down through the other classes to we propagate
        # the spectrograph defaults as is done here...
        self.par = self.spectrograph.default_pypeit_par()['calibrations'] if par is None else par
        if not isinstance(self.par, pypeitpar.CalibrationsPar):
            raise TypeError('Input parameters must be a CalibrationsPar instance.')

        # Output dirs
        self.redux_path = os.getcwd() if redux_path is None else redux_path
        self.master_dir = masterframe.set_master_dir(self.redux_path, self.spectrograph, self.par)

        # Attributes
        self.calib_dict = {}
        self.det = None
        self.sci_ID = None
        self.setup = None

        # Steps
        self.steps = []

        # Internals
        self._reset_internals()

    def _reset_internals(self):
        """
        Reset all of the key internals to None

        Returns:

        """
        self.shape = None
        self.msarc = None
        self.msbias = None
        self.msbpm = None
        self.pixlocn = None
        self.tslits_dict = None
        self.maskslits = None
        self.wavecalib = None
        self.tilts_dict = None
        self.mspixflatnrm = None
        self.msillumflat = None
        self.mswave = None
        #self.datasec_img = None

    def reset(self, setup, det, sci_ID, par=None):
        """
        Specify the parameters of the Calibrations class and reset all
        the internals to None. The internal dict is left unmodified.

        Args:
            setup (str):
            det (int):
            sci_ID (int):
            par (CalibrationPar):

        Returns:

        """
        self.setup = setup
        self.det = det
        self.sci_ID = sci_ID
        if par is not None:
            self.par = par

        # Setup the calib_dict
        if self.setup not in self.calib_dict.keys():
            self.calib_dict[self.setup] = {}

        # Reset internals to None
        self._reset_internals()

    def get_arc(self):
        """
        Load or generate the bias frame/command

        Requirements:
          self.msbias
          setup, det, sci_ID, par

        Args:

        Returns:
            self.msarc: ndarray

        """
        # Check for existing data
        if not self._chk_objs(['msbias']):
            self.msarc = None
            return self.msarc

        # Check internals
        self._chk_set(['setup', 'det', 'sci_ID', 'par'])

        if 'arc' in self.calib_dict[self.setup].keys():
            # Previously calculated
            self.msarc = self.calib_dict[self.setup]['arc']
            return self.msarc

        self.arc_file_list = self.fitstbl.find_frame_files('arc', sci_ID=self.sci_ID)
        # Instantiate with everything needed to generate the image (in case we do)
        self.arcImage = arcimage.ArcImage(self.spectrograph, file_list = self.arc_file_list, det=self.det,msbias=self.msbias,
                                          par=self.par['arcframe'], setup=self.setup,
                                          master_dir=self.master_dir, mode=self.par['masters'])
        
        # Load the MasterFrame (if it exists and is desired)?
        self.msarc = self.arcImage.master()
        if self.msarc is None:  # Otherwise build it
            msgs.info("Preparing a master {0:s} frame".format(self.arcImage.frametype))
            self.msarc = self.arcImage.build_image()
            # Save to Masters
            if self.save_masters:
                self.arcImage.save_master(self.msarc, raw_files=self.arcImage.file_list,
                                          steps=self.arcImage.steps)

        # Save & return
        self.calib_dict[self.setup]['arc'] = self.msarc
        return self.msarc

    def get_bias(self):
        """
        Load or generate the bias frame/command

        Requirements:
           setup, det, sci_ID, par

        Returns:
            self.msbias: ndarray or str

        """
        # Check internals
        self._chk_set(['setup', 'det', 'sci_ID', 'par'])

        # Grab from internal dict?
        if 'bias' in self.calib_dict[self.setup].keys():
            self.msbias = self.calib_dict[self.setup]['bias']
            msgs.info("Reloading the bias from the internal dict")
            return self.msbias

        self.bias_file_list = self.fitstbl.find_frame_files('bias', sci_ID=self.sci_ID)
        # Instantiate
        self.biasFrame = biasframe.BiasFrame(self.spectrograph, file_list = self.bias_file_list, det=self.det,
                                             par=self.par['biasframe'], setup=self.setup,
                                             master_dir=self.master_dir, mode=self.par['masters'])

        # How are we treating biases: 1) No bias, 2) overscan, or 3) use bias subtraction. If use bias is there a master?
        self.msbias = self.biasFrame.determine_bias_mode()

        if self.msbias is None:  # Build it and save it
            self.msbias = self.biasFrame.build_image()
            if self.save_masters:
                self.biasFrame.save_master(self.msbias, raw_files=self.biasFrame.file_list,steps=self.biasFrame.steps)

        # Save & return
        self.calib_dict[self.setup]['bias'] = self.msbias
        return self.msbias

    def get_bpm(self):
        """
        Load or generate the bad pixel mask

        Requirements:
           Instrument dependent

        Returns:
            self.msbpm: ndarray

        """
        # Check internals
        self._chk_set(['sci_ID', 'par'])

        # Generate a bad pixel mask (should not repeat)
        if 'bpm' in self.calib_dict[self.setup].keys():
            self.msbpm = self.calib_dict[self.setup]['bpm']

        # Make sure shape is defined
        self._check_shape()

        # Always use the example file
        example_file = self.fitstbl.find_frame_files('science', sci_ID=self.sci_ID)[0]
        bpmImage = bpmimage.BPMImage(self.spectrograph, filename=example_file, det=self.det,
                                     msbias=self.msbias if self.par['badpix'] == 'bias' else None,
                                     trim=self.par['trim'])
        # Build, save, and return
        self.msbpm = bpmImage.build()
        self.calib_dict[self.setup]['bpm'] = self.msbpm
        return self.msbpm

    def get_flats(self, show=False):
        """
        Load or generate a normalized pixel flat
          and slit profile

        Requirements:
           tslits_dict
           tilts_dict
           det, sci_ID, par, setup

        Returns:
            self.mspixflatnrm: ndarray
            self.msillumflat: ndarray

        """

        if self.par['flatfield']['method'] is None:
            # User does not want to flat-field
            self.mspixflatnrm = None
            self.msillumflat = None
            msgs.warning('Parameter calibrations.flatfield.method is set to None. You are NOT '
                         'flatfielding your data!!!')
            return self.mspixflatnrm, self.msillumflat

        # Check for existing data necessary to build flats
        if not self._chk_objs(['tslits_dict', 'tilts_dict']):
            msgs.warning('Flats were requested, but there are quantities missing necessary to '
                         'create flats.  Proceeding without flat fielding....')
            # User cannot flat-field
            self.mspixflatnrm = None
            self.msillumflat = None
            return self.mspixflatnrm, self.msillumflat

        # Check internals
        self._chk_set(['setup', 'det', 'sci_ID', 'par'])

        # Return already generated data
        if np.all([k in self.calib_dict[self.setup].keys()
                   for k in ['normpixelflat', 'illumflat']]):
            self.mspixflatnrm = self.calib_dict[self.setup]['normpixelflat']
            self.msillumflat = self.calib_dict[self.setup]['illumflat']
            return self.mspixflatnrm, self.msillumflat

        # Instantiate
        pixflat_image_files = self.fitstbl.find_frame_files('pixelflat', sci_ID=self.sci_ID)
        self.flatField = flatfield.FlatField(self.spectrograph, file_list=pixflat_image_files,
                                             det=self.det, par=self.par['pixelflatframe'],
                                             setup=self.setup, master_dir=self.master_dir,
                                             mode=self.par['masters'],
                                             flatpar=self.par['flatfield'], msbias=self.msbias,
                                             tslits_dict=self.tslits_dict,
                                             tilts_dict=self.tilts_dict)

        # --- Pixel flats

        # 1)  Try to load a master frile from disk (MasterFrame)?
        self.mspixflatnrm = self.flatField.master()

        # 2) Did the user specify a flat? If so load it in  (e.g. LRISb with pixel flat)?
        if self.par['flatfield']['frame'] not in ['pixelflat']:
            # First try to find directly, then try to find it in the
            # masters directory, then fail
            if os.path.isfile(self.par['flatfield']['frame']):
                mspixelflat_name = self.par['flatfield']['frame']
            elif os.path.isfile(os.path.join(self.flatField.directory_path,
                                             self.par['flatfield']['frame'])):
                mspixelflat_name = os.path.join(self.flatField.directory_path,
                                                self.par['flatfield']['frame'])
            else:
                raise ValueError('Could not find user-defined flatfield master: {0}'.format(
                    self.par['flatfield']['frame']))
            msgs.info('Found user-defined file: {0}'.format(mspixelflat_name))
            self.mspixflatnrm = self.flatField.load_master(mspixelflat_name, exten=self.det)

        # 3) there is no master or no user supplied flat, generate the flat
        if self.mspixflatnrm is None and len(pixflat_image_files) != 0:
            # Run
            self.mspixflatnrm, self.msillumflat = self.flatField.run(show=self.show)

            # If we tweaked the slits, update the tilts_dict and
            # tslits_dict to reflect new slit edges
            if self.par['flatfield']['tweak_slits']:
                msgs.info('Using slit boundary tweaks from IllumFlat and updated tilts image')
                self.tslits_dict = self.flatField.tslits_dict
                self.tilts_dict = self.flatField.tilts_dict

            # Save to Masters
            if self.save_masters:
                self.flatField.save_master(self.mspixflatnrm, raw_files=pixflat_image_files,
                                           steps=self.flatField.steps)
                self.flatField.save_master(self.msillumflat, raw_files=pixflat_image_files,
                                           steps=self.flatField.steps,
                                           outfile=masterframe.master_name('illumflat', self.setup,self.master_dir))
                # If we tweaked the slits update the master files for tilts and slits
                if self.par['flatfield']['tweak_slits']:
                    msgs.info('Updating MasterTrace and MasterTilts using tweaked slit boundaries')
                    # Add tweaked boundaries to the MasterTrace file
                    self.traceSlits.lcen_tweak = self.flatField.tslits_dict['lcen']
                    self.traceSlits.rcen_tweak = self.flatField.tslits_dict['rcen']
                    self.traceSlits.save_master()
                    # Write the final_tilts using the new slit boundaries to the MasterTilts file
                    self.waveTilts.final_tilts = self.flatField.tilts_dict['tilts']
                    self.waveTilts.save_master()

        # 4) If we still don't have a pixel flat, then just use unity
        # everywhere and print out a warning
        if self.mspixflatnrm is None:
            self.mspixflatnrm = np.ones_like(self.tilts_dict['tilts'])
            msgs.warn('You are not pixel flat fielding your data!!!')

        # --- Illumination flats

        # 1) If we ran the flat field algorithm above, then the
        # illumination file was created. So check msillumflat is set
        if self.msillumflat is None:
            # 2) If no illumination file is set yet, try to read it in from a master
            self.msillumflat = self.flatField.load_master_illumflat()
            # 3) If there is no master file, then set illumflat to unit
            # and war user that they are not illumflatting their data
            if self.msillumflat is None:
                self.msillumflat = np.ones_like(self.tilts_dict['tilts'])
                msgs.warn('You are not illumination flat fielding your data!')

        # Save & return
        self.calib_dict[self.setup]['normpixelflat'] = self.mspixflatnrm
        self.calib_dict[self.setup]['illumflat'] = self.msillumflat

        return self.mspixflatnrm, self.msillumflat

    def get_slits(self, arms=True, redo=False):
        """
        Load or generate the slits.
        First, a trace flat image is generated

        .. todo::
            - arms is a parameter passed to traceSlits.  This may need
              to change if/when arms.py is replaced.

        Requirements:
           pixlocn
           det par setup sci_ID

        Returns:
            self.tslits_dict
            self.maskslits

        """
        # Check for existing data
        if not self._chk_objs(['pixlocn', 'msbpm']):
            self.tslits_dict = None
            self.maskslits = None
            return self.tslits_dict, self.maskslits

        # Check internals
        self._chk_set(['setup', 'det', 'sci_ID', 'par'])

        # Return already generated data
        if ('trace' in self.calib_dict[self.setup].keys()) and (not redo):
            self.tslits_dict = self.calib_dict[self.setup]['trace']
            self.maskslits = np.zeros(self.tslits_dict['lcen'].shape[1], dtype=bool)
            return self.tslits_dict, self.maskslits
                
        # Instantiate (without mstrace)
        self.traceSlits = traceslits.TraceSlits(None, self.pixlocn, par=self.par['slits'],
                                                det=self.det, setup=self.setup,
                                                master_dir=self.master_dir,
                                                redux_path=self.redux_path,
                                                mode=self.par['masters'], binbpx=self.msbpm)

        # Load via master, as desired
        if not self.traceSlits.master():
            # Build the trace image first
            self.trace_image_files = self.fitstbl.find_frame_files('trace', sci_ID=self.sci_ID)
            self.traceImage = traceimage.TraceImage(self.spectrograph,self.trace_image_files, det=self.det,
                                           par=self.par['traceframe'])
            # Load up and get ready
            self.traceSlits.mstrace = self.traceImage.process(bias_subtract=self.msbias,
                                                         trim=self.par['trim'], apply_gain=True)

            # Compute the plate scale in arcsec which is needed to trim short slits
            scidx = np.where(self.fitstbl.find_frames('science', sci_ID=self.sci_ID))[0][0]
            try:
                binspatial, binspectral = parse.parse_binning(self.fitstbl['binning'][scidx])
            except:
                binspatial, binspectral = 1,1
            ## Old code: binspatial, binspectral = parse.parse_binning(self.fitstbl['binning'][scidx])
            plate_scale = binspatial*self.spectrograph.detector[self.det-1]['platescale']

            # Now we go forth
            try:
                self.tslits_dict = self.traceSlits.run(arms=arms, plate_scale = plate_scale)
            except:
                self.traceSlits.save_master()
                msgs.error("Crashed out of finding the slits. Have saved the work done to disk but it needs fixing..")
            # No slits?
            if self.tslits_dict is None:
                self.maskslits = None
                return self.tslits_dict, self.maskslits
            # QA
            if self.write_qa:
                self.traceSlits._qa()
            # Save to disk
            if self.save_masters:
                # Master
                self.traceSlits.save_master()

        # Construct dictionary
        self.tslits_dict = self.traceSlits._fill_tslits_dict()

        # Save, initialize maskslits, and return
        self.calib_dict[self.setup]['trace'] = self.tslits_dict
        self.maskslits = np.zeros(self.tslits_dict['lcen'].shape[1], dtype=bool)
        return self.tslits_dict, self.maskslits

    def get_wave(self):
        """
        Load or generate a wavelength image

        Requirements:
           tilts_dict, tslits_dict, wv_calib, maskslits
           det, par, setup

        Returns:
            self.mswave: ndarray

        """
        # Check for existing data
        if not self._chk_objs(['tilts_dict','tslits_dict','wv_calib','maskslits']):
            self.mswave = None
            return self.mswave

        # Check internals
        self._chk_set(['setup', 'det', 'par'])

        # Return existing data
        if 'wave' in self.calib_dict[self.setup].keys():
            self.mswave = self.calib_dict[self.setup]['wave']
            return self.mswave

        # No wavelength calibration requested
        if self.par['wavelengths']['reference'] == 'pixel':
            self.mswave = self.tilts_dict['tilts'] * (self.tilts_dict['tilts'].shape[0]-1.0)
            self.calib_dict[self.setup]['wave'] = self.mswave
            return self.mswave

        # Instantiate
        # ToDO we are regenerating this mask a lot in this module. Could reduce that
        self.slitmask = self.spectrograph.slitmask(self.tslits_dict)
        self.waveImage = waveimage.WaveImage(self.slitmask,self.tilts_dict['tilts'], self.wv_calib,
                                             setup=self.setup, master_dir=self.master_dir,
                                             mode=self.par['masters'], maskslits=self.maskslits)
        # Attempt to load master
        self.mswave = self.waveImage.master()
        if self.mswave is None:
            self.mswave = self.waveImage._build_wave()
        # Save to hard-drive
        if self.save_masters:
            self.waveImage.save_master(self.mswave, steps=self.waveImage.steps)

        # Save & return
        self.calib_dict[self.setup]['wave'] = self.mswave
        return self.mswave

    def get_wv_calib(self):
        """
        Load or generate the 1D wavelength calibrations

        Requirements:
          msarc, msbpm, tslits_dict
          det, par, setup, sci_ID

        Returns:
            self.wv_calib: dict
            self.maskslits -- Updated
        """
        # Check for existing data
        if not self._chk_objs(['msarc', 'msbpm', 'tslits_dict']):
            msgs.error('dont have all the objects')
            self.wv_calib = None
            self.wv_maskslits = np.zeros_like(self.maskslits, dtype=bool)
            self.maskslits += self.wv_maskslits
            return self.wv_calib, self.maskslits

        # Check internals
        self._chk_set(['setup', 'det', 'sci_ID', 'par'])

        # Return existing data
        if 'wavecalib' in self.calib_dict[self.setup].keys():
            self.wv_calib = self.calib_dict[self.setup]['wavecalib']
            self.wv_maskslits = self.calib_dict[self.setup]['wvmask']
            self.maskslits += self.wv_maskslits
            return self.wv_calib, self.maskslits

        # No wavelength calibration requested
        if self.par['wavelengths']['reference'] == 'pixel':
            msgs.info("A wavelength calibration will not be performed")
            self.wv_calib = None
            self.wv_maskslits = np.zeros_like(self.maskslits, dtype=bool)
            self.maskslits += self.wv_maskslits
            return self.wv_calib, self.maskslits

        # Setup
        nonlinear = self.spectrograph.detector[self.det-1]['saturation'] \
                        * self.spectrograph.detector[self.det-1]['nonlinear']

        # Instantiate
        self.waveCalib = wavecalib.WaveCalib(self.msarc, spectrograph=self.spectrograph,det=self.det,
                                             par=self.par['wavelengths'], setup=self.setup, master_dir=self.master_dir,
                                             mode=self.par['masters'],redux_path=self.redux_path, bpm=self.msbpm)
        # Load from disk (MasterFrame)?
        self.wv_calib = self.waveCalib.master()
        # Build?
        if self.wv_calib is None:
            self.slitmask = self.spectrograph.slitmask(self.tslits_dict)
            self.wv_calib, _ = self.waveCalib.run(self.tslits_dict['lcen'],
                                                  self.tslits_dict['rcen'],
                                                  self.slitmask,
                                                  nonlinear=nonlinear, skip_QA=(not self.write_qa))
            # Save to Masters
            if self.save_masters:
                self.waveCalib.save_master(self.waveCalib.wv_calib)
        else:
            self.waveCalib.wv_calib = self.wv_calib

        # Create the mask
        self.wv_maskslits = self.waveCalib._make_maskslits(self.tslits_dict['lcen'].shape[1])

        # Save & return
        self.calib_dict[self.setup]['wavecalib'] = self.wv_calib
        self.calib_dict[self.setup]['wvmask'] = self.wv_maskslits
        self.maskslits += self.wv_maskslits
        return self.wv_calib, self.maskslits

    def get_pixlocn(self):
        """
        Generate the pixlocn image

        Requirements:
          spectrograph, shape

        Returns:
            self.pixlocn: ndarray
        """
        # Make sure shape is defined
        self._check_shape()
        # Check internals
        self._chk_set(['shape'])

        # Get the pixel locations
        xgap=self.spectrograph.detector[self.det-1]['xgap']
        ygap=self.spectrograph.detector[self.det-1]['ygap']
        ysize=self.spectrograph.detector[self.det-1]['ysize']
        self.pixlocn = pixels.gen_pixloc(self.shape, xgap=xgap, ygap=ygap, ysize=ysize)

        # Return
        return self.pixlocn

    def get_tilts(self):
        """
        Load or generate the tilts image

        Requirements:
           msarc, tslits_dict, pixlocn, wv_calib, maskslits
           det, par, setup, sci_ID, spectrograph

        Returns:
            self.tilts_dict: dictionary with tilts information (2D)
            self.maskslits: ndarray

        """
        # Check for existing data
        if not self._chk_objs(['msarc', 'msbpm', 'tslits_dict', 'wv_calib', 'maskslits']):
            msgs.error('dont have all the objects')
            self.tilts_dict = None
            self.wt_maskslits = np.zeros_like(self.maskslits, dtype=bool)
            self.maskslits += self.wt_maskslits
            return self.tilts_dict, self.maskslits

        # Check internals
        self._chk_set(['setup', 'det', 'sci_ID', 'par'])

        # Return existing data
        if 'tilts_dict' in self.calib_dict[self.setup].keys():
            self.tilts_dict = self.calib_dict[self.setup]['tilts_dict']
            self.wt_maskslits = self.calib_dict[self.setup]['wtmask']
            self.maskslits += self.wt_maskslits
            return self.tilts_dict, self.maskslits

        # Instantiate
        self.waveTilts = wavetilts.WaveTilts(self.msarc, spectrograph=self.spectrograph,
                                             par=self.par['tilts'], det=self.det,
                                             setup=self.setup, master_dir=self.master_dir,
                                             mode=self.par['masters'],
                                             tslits_dict=self.tslits_dict,
                                             redux_path=self.redux_path, bpm=self.msbpm)
        # Master
        self.tilts_dict = self.waveTilts.master()
        if self.tilts_dict is None:
            self.tilts_dict, self.wt_maskslits \
                    = self.waveTilts.run(maskslits=self.maskslits, wv_calib=self.wv_calib,
                                         doqa=self.write_qa)
            if self.save_masters:
                self.waveTilts.save_master()
        else:
            self.wt_maskslits = np.zeros_like(self.maskslits, dtype=bool)

        # Save & return
        self.calib_dict[self.setup]['tilts_dict'] = self.tilts_dict
        self.calib_dict[self.setup]['wtmask'] = self.wt_maskslits
        self.maskslits += self.wt_maskslits
        return self.tilts_dict, self.maskslits

    def run_the_steps(self):
        """
        Run full the full recipe of calibration steps

        Returns:

        """
        for step in self.steps:
            getattr(self, 'get_{:s}'.format(step))()

    # This is general to any attribute
    def _chk_set(self, items):
        for item in items:
            if getattr(self, item) is None:
                msgs.error("Use self.set to specify '{:s}' prior to generating XX".format(item))

    # This is specific to `self.ms*` attributes
    def _chk_objs(self, items):
        for obj in items:
            if getattr(self, obj) is None:
                msgs.warn("You need to generate {:s} prior to this calibration..".format(obj))
                if obj in ['tslits_dict','maskslits']:
                    msgs.warn("Use get_slits".format(obj))
                else:
                    # Strip ms
                    if obj[0:2] == 'ms':
                        iobj = obj[2:]
                    else:
                        iobj = obj
                    msgs.warn("Use get_{:s}".format(iobj))
                return False
        return True

    def _check_shape(self):
        """
        Check that the shape attribute is not None.  If it is use,
        define it using the shape of msarc

        .. warning::
            - This shape depends on if the images are trimmed or not!
        """
        # Check the shape is declared
        if self.shape is None and self.msarc is None:
            raise ValueError('Before calling BPM, must run get_arc to get image shape, or '
                             'provide shape directly.')
        if self.shape is None:
            self.shape = self.msarc.shape

    def show(self, obj):
        if isinstance(obj, np.ndarray):
            if len(obj.shape) == 2:
                debugger.show_image(obj)
        else:
            msgs.warn("Not ready for this type of object")

    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: setup={}, det={}, sci_ID={}'.format(self.__class__.__name__,
                                                          self.setup,
                                                          self.det,
                                                          self.sci_ID)
        txt += '>'
        return txt



class MultiSlitCalibrations(Calibrations):
    """
    Child of Calibrations class for performing multi-slit (and longslit)
    calibrations.
    """
    def __init__(self, fitstbl, spectrograph=None, par=None, redux_path=None, save_masters=True,
                 write_qa=True, show = False, steps=None):
        Calibrations.__init__(self, fitstbl, spectrograph=spectrograph, par=par,
                              redux_path=redux_path, save_masters=save_masters,
                              write_qa=write_qa, show = show)
        self.steps = MultiSlitCalibrations.default_steps() if steps is None else steps

    @staticmethod
    def default_steps():
        return ['bias', 'arc', 'bpm', 'pixlocn', 'slits', 'wv_calib', 'tilts',
                'flats', 'wave']


