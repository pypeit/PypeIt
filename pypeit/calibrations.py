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
from pypeit.core import trace_slits

from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph

from pypeit import debugger


class Calibrations(object):
    """
    This class is primarily designed to guide the generation of
    calibration images and objects in PypeIt.

    .. todo::
        Improve docstring...

    Args:
        fitstbl (:class:`pypeit.metadata.PypeItMetaData`):
            The class holding the metadata for all the frames in this
            PypeIt run.
        spectrograph (:obj:`str`, optional):
            The name of the spectrograph to reduce.  TODO: Not needed.
        par (:class:`pypeit.par.pypeitpar.PypeItPar`, optional):
            Parameter set defining optional parameters of PypeIt's
            low-level algorithms.  If None, defined using
            :func:`pypeit.spectrograph.spectrographs.Spectrograph.default_pypeit_par`.
        redux_path (:obj:`str`, optional):
            Top-level directory for PypeIt output.  If None, the current
            working directory is used.
        save_masters (:obj:`bool`, optional):
            Save Master files as they are generated for later use.
        write_qa (:obj:`bool`, optional):
            Create QA plots.
        show (:obj:`bool`, optional):
            Show plots of PypeIt's results as the code progesses.
            Requires interaction from the users.

    Attributes:
        fitstbl
        save_masters
        write_qa
        show
        spectrograph
        par
        redux_path
        master_dir
        calib_dict
        det
        frame (:obj:`int`):
            0-indexed row of the frame being calibrated in
            :attr:`fitstbl`.
        master_key


        
    """
    __metaclass__ = ABCMeta

    def __init__(self, fitstbl, spectrograph=None, par=None, redux_path=None, save_masters=True,
                 write_qa=True, show=False):

        # Check the type of the provided fits table
        if not isinstance(fitstbl, PypeItMetaData):
            msgs.error('fitstbl must be an PypeItMetaData object')

        # Parameters unique to this Object
        self.fitstbl = fitstbl
        self.save_masters = save_masters
        self.write_qa = write_qa
        self.show = show

        # Spectrometer class
        # TODO: the spectrograph is already defined in fitstbl
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
        self.frame = None

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

    def check_for_previous(self, ftype, master_key):
        """
        Check to see whether the master file(s) have been
        built during this run of PypeIt

        Args:
            ftype: str
            master_key: str

        Returns:
            previous: bool
               True = Built previously
        """
        previous = False
        if master_key in self.calib_dict.keys():
            if ftype in self.calib_dict[master_key].keys():
                previous = True
        else:  # Prep for saving
            self.calib_dict[master_key] = {}
        # Return
        return previous

    def set_config(self, frame, det, par=None):
        """
        Specify the parameters of the Calibrations class and reset all
        the internals to None. The internal dict is left unmodified.

        Args:
            master_key (str):
            det (int):
            par (CalibrationPar):

        Returns:

        """
        # TODO is the right behavior to just take the first one
        self.frame = frame
        self.calib_ID = int(self.fitstbl['calib'][frame])
        self.det = det
        if par is not None:
            self.par = par

        # Reset internals to None
        self._reset_internals()

    def get_arc(self):
        """
        Load or generate the bias frame/command

        Requirements:
          self.msbias
          master_key, det, par

        Args:

        Returns:
            self.msarc: ndarray

        """
        # Check for existing data
        ## JFH This check is wrong, if the user does not want to bias subtract, then it causes a crash
#        if not self._chk_objs(['msbias']):
#            self.msarc = None
#            return self.msarc        self.arc_file_list, arc_rows = self.fitstbl.find_frame_files('arc', calib_ID=self.calib_ID)

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        self.arc_file_list, arc_rows = self.fitstbl.find_frame_files('arc', calib_ID=self.calib_ID)
        self.arc_master_key = self.fitstbl.master_key(arc_rows[0], det=self.det)

        prev_build = self.check_for_previous('arc', self.arc_master_key)
        if prev_build:
            # Previously calculated
            self.msarc = self.calib_dict[self.arc_master_key]['arc']
            return self.msarc

        # Instantiate with everything needed to generate the image (in case we do)
        self.arcImage = arcimage.ArcImage(self.spectrograph, file_list = self.arc_file_list, det=self.det,msbias=self.msbias,
                                          par=self.par['arcframe'], master_key=self.arc_master_key,
                                          master_dir=self.master_dir, mode=self.par['masters'])

        # Load the MasterFrame (if it exists and is desired)?
        self.msarc = self.arcImage.master(force=prev_build)
        if self.msarc is None:  # Otherwise build it
            msgs.info("Preparing a master {0:s} frame".format(self.arcImage.frametype))
            self.msarc = self.arcImage.build_image()
            # Save to Masters
            if self.save_masters:
                self.arcImage.save_master(self.msarc, raw_files=self.arcImage.file_list,
                                          steps=self.arcImage.steps)

        # Save & return
        self.calib_dict[self.arc_master_key]['arc'] = self.msarc
        return self.msarc

    def get_bias(self):
        """
        Load or generate the bias frame/command

        Requirements:
           master_key, det, par

        Returns:
            self.msbias: ndarray or str

        """
        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        self.bias_file_list, bias_rows = self.fitstbl.find_frame_files('bias', calib_ID=self.calib_ID)
        if len(bias_rows) > 0:
            self.bias_master_key = self.fitstbl.master_key(bias_rows[0], det=self.det)
        else:  # Allow for user-supplied file (e.g. LRISb)
            self.bias_master_key = self.fitstbl.master_key(self.frame, det=self.det)

        # Grab from internal dict?
        prev_build = self.check_for_previous('bias', self.bias_master_key)
        if prev_build:
            self.msbias = self.calib_dict[self.bias_master_key]['bias']
            msgs.info("Reloading the bias from the internal dict")
            return self.msbias

        # Instantiate
        self.biasFrame = biasframe.BiasFrame(self.spectrograph, file_list = self.bias_file_list, det=self.det,
                                             par=self.par['biasframe'], master_key=self.bias_master_key,
                                             master_dir=self.master_dir, mode=self.par['masters'])

        # How are we treating biases: 1) No bias, 2) overscan, or 3) use bias subtraction. If use bias is there a master?
        self.msbias = self.biasFrame.determine_bias_mode(force=prev_build)
        # This could be made more elegant, like maybe msbias should be set to 'none' analgous to how overscan is treated???
        if (self.msbias is None) and (self.par['biasframe']['useframe'] != 'none'):  # Build it and save it
            self.msbias = self.biasFrame.build_image()
            if self.save_masters:
                self.biasFrame.save_master(self.msbias, raw_files=self.biasFrame.file_list,steps=self.biasFrame.steps)

        # Save & return
        self.calib_dict[self.bias_master_key]['bias'] = self.msbias
        return self.msbias

    def get_bpm(self):
        """
        Load or generate the bad pixel mask

        TODO -- Should consider doing this outside of calibrations as it is
        more specific to the science frame

        This needs to be for the *trimmed* image!

        Requirements:
           Instrument dependent

        Returns:
            self.msbpm: ndarray

        """
        # Check internals
        self._chk_set(['par'])


        # Generate a bad pixel mask (should not repeat)
        self.bpm_master_key = self.fitstbl.master_key(self.frame, det=self.det)
        prev_build = self.check_for_previous('bpm', self.bpm_master_key)
        if prev_build:
            self.msbpm = self.calib_dict[self.bpm_master_key]['bpm']
            return self.msbpm

        # Make sure shape is defined
        self._check_shape()

        # Always use the shape!
        #  But some instruments need the filename too, e.g. for binning
        sci_image_files = [self.fitstbl.frame_paths(self.frame)]
        dsec_img = self.spectrograph.get_datasec_img(sci_image_files[0], det=self.det)
        shape = procimg.trim_frame(dsec_img, dsec_img < 1).shape
        # Check it matches the processed arc;  if not we have issues..
        if not (self.shape == shape):
            msgs.error("You have an untrimmed arc!  We aren't prepared for this..")

        # Build it
        bpmImage = bpmimage.BPMImage(self.spectrograph, filename=sci_image_files[0],
                                     det=self.det, shape=self.shape,
                                     msbias=self.msbias if self.par['badpix'] == 'bias' else None)
        # Build, save, and return
        self.msbpm = bpmImage.build()
        self.calib_dict[self.bpm_master_key]['bpm'] = self.msbpm
        return self.msbpm

    def get_flats(self, show=False):
        """
        Load or generate a normalized pixel flat
          and slit profile

        Requirements:
           tslits_dict
           tilts_dict
           det, par

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
        self._chk_set(['det', 'calib_ID', 'par'])

        pixflat_image_files, pixflat_rows = self.fitstbl.find_frame_files('pixelflat', calib_ID=self.calib_ID)
        if len(pixflat_rows) > 0:
            self.pixflat_master_key = self.fitstbl.master_key(pixflat_rows[0], det=self.det)
        else:  # Allow for user-supplied file (e.g. LRISb)
            self.pixflat_master_key = self.fitstbl.master_key(self.frame, det=self.det)

        # Return already generated data
        prev_build1 = self.check_for_previous('normpixelflat', self.pixflat_master_key)
        prev_build2 = self.check_for_previous('illumflat', self.pixflat_master_key)
        if np.all([prev_build1, prev_build2]):
            self.mspixflatnrm = self.calib_dict[self.pixflat_master_key]['normpixelflat']
            self.msillumflat = self.calib_dict[self.pixflat_master_key]['illumflat']
            return self.mspixflatnrm, self.msillumflat

        # Instantiate
        self.flatField = flatfield.FlatField(self.spectrograph, file_list=pixflat_image_files,
                                             det=self.det, par=self.par['pixelflatframe'],
                                             master_key=self.pixflat_master_key, master_dir=self.master_dir,
                                             mode=self.par['masters'],
                                             flatpar=self.par['flatfield'], msbias=self.msbias,
                                             tslits_dict=self.tslits_dict,
                                             tilts_dict=self.tilts_dict)

        # --- Pixel flats

        # 1)  Try to load master files from disk (MasterFrame)?
        self.mspixflatnrm = self.flatField.master(force=prev_build1)
        if prev_build2:
            self.msillumflat = self.flatField.load_master_illumflat()

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
                                           outfile=masterframe.master_name('illumflat',
                                                                           self.pixflat_master_key,
                                                                           self.master_dir))
                # If we tweaked the slits update the master files for tilts and slits
                if self.par['flatfield']['tweak_slits']:
                    msgs.info('Updating MasterTrace and MasterTilts using tweaked slit boundaries')
                    # Add tweaked boundaries to the MasterTrace file
                    self.traceSlits.lcen_tweak = self.flatField.tslits_dict['lcen']
                    self.traceSlits.rcen_tweak = self.flatField.tslits_dict['rcen']
                    self.traceSlits.save_master()
                    # Write the final_tilts using the new slit boundaries to the MasterTilts file
                    self.waveTilts.final_tilts = self.flatField.tilts_dict['tilts']
                    self.waveTilts.tilts_dict = self.flatField.tilts_dict
                    self.waveTilts.save_master(self.flatField.tilts_dict, steps=self.waveTilts.steps)

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
        self.calib_dict[self.pixflat_master_key]['normpixelflat'] = self.mspixflatnrm
        self.calib_dict[self.pixflat_master_key]['illumflat'] = self.msillumflat

        return self.mspixflatnrm, self.msillumflat

    def get_slits(self, redo=False):
        """
        Load or generate the slits.
        First, a trace flat image is generated

        Requirements:
           pixlocn
           det par master_key

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
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        self.trace_image_files, trace_rows = self.fitstbl.find_frame_files('trace', calib_ID=self.calib_ID)
        self.trace_master_key = self.fitstbl.master_key(trace_rows[0], det=self.det)

        # Return already generated data
        prev_build = self.check_for_previous('trace', self.trace_master_key)
        if prev_build and (not redo):
            self.tslits_dict = self.calib_dict[self.trace_master_key]['trace']
            self.maskslits = np.zeros(self.tslits_dict['lcen'].shape[1], dtype=bool)
            return self.tslits_dict, self.maskslits
                
        # Instantiate (without mstrace)

        self.traceSlits = traceslits.TraceSlits(None, self.pixlocn, self.spectrograph,
                                                par=self.par['slits'],
                                                det=self.det, master_key=self.trace_master_key,
                                                master_dir=self.master_dir,
                                                redux_path=self.redux_path,
                                                mode=self.par['masters'], binbpx=self.msbpm)

        # Load via master, as desired
        if not self.traceSlits.master(force=prev_build):
            # Build the trace image first
            self.traceImage = traceimage.TraceImage(self.spectrograph,self.trace_image_files, det=self.det,
                                           par=self.par['traceframe'])
            # Load up and get ready
            self.traceSlits.mstrace = self.traceImage.process(bias_subtract=self.msbias,
                                                         trim=self.par['trim'], apply_gain=True)

            # Compute the plate scale in arcsec which is needed to trim short slits
            try:
                binspatial, binspectral = parse.parse_binning(self.fitstbl['binning'][self.frame])
            except:
                binspatial, binspectral = 1,1
            ## Old code: binspatial, binspectral = parse.parse_binning(self.fitstbl['binning'][scidx])
            plate_scale = binspatial*self.spectrograph.detector[self.det-1]['platescale']

            # User-defined slits??
            add_user_slits = trace_slits.parse_user_slits(self.par['slits']['add_slits'], self.det)
            rm_user_slits = trace_slits.parse_user_slits(self.par['slits']['rm_slits'], self.det, rm=True)

            # Now we go forth
            try:
                self.tslits_dict = self.traceSlits.run(plate_scale = plate_scale, show=self.show,
                                                       add_user_slits=add_user_slits, rm_user_slits=rm_user_slits)
            except:
                self.traceSlits.save_master()
                msgs.error("Crashed out of finding the slits. Have saved the work done to disk but it needs fixing..")
            # No slits?
            if self.tslits_dict is None:
                self.maskslits = None
                return self.tslits_dict, self.maskslits
            # Save to disk
            if self.save_masters:
                # Master
                self.traceSlits.save_master()
        else:
            msgs.info("TraceSlits master files loaded..")
            # Construct dictionary
            self.tslits_dict = self.traceSlits._fill_tslits_dict()

        # Save, initialize maskslits, and return
        self.calib_dict[self.trace_master_key]['trace'] = self.tslits_dict
        self.maskslits = np.zeros(self.tslits_dict['lcen'].shape[1], dtype=bool)
        return self.tslits_dict, self.maskslits

    def get_wave(self):
        """
        Load or generate a wavelength image

        Requirements:
           tilts_dict, tslits_dict, wv_calib, maskslits
           det, par, master_key

        Returns:
            self.mswave: ndarray

        """
        # Check for existing data
        if not self._chk_objs(['tilts_dict','tslits_dict','wv_calib','maskslits']):
            self.mswave = None
            return self.mswave

        # Check internals
        self._chk_set(['arc_master_key', 'det', 'par'])

        # Return existing data
        prev_build = self.check_for_previous('wave', self.arc_master_key)
        if prev_build:
            self.mswave = self.calib_dict[self.arc_master_key]['wave']
            return self.mswave

        # No wavelength calibration requested
        if self.par['wavelengths']['reference'] == 'pixel':
            self.mswave = self.tilts_dict['tilts'] * (self.tilts_dict['tilts'].shape[0]-1.0)
            self.calib_dict[self.arc_master_key]['wave'] = self.mswave
            return self.mswave

        # Instantiate
        # ToDO we are regenerating this mask a lot in this module. Could reduce that
        self.waveImage = waveimage.WaveImage(self.tslits_dict, self.tilts_dict['tilts'], self.wv_calib,self.spectrograph,
                                             master_key=self.arc_master_key, master_dir=self.master_dir,
                                             mode=self.par['masters'], maskslits=self.maskslits)
        # Attempt to load master
        self.mswave = self.waveImage.master(force=prev_build)
        if self.mswave is None:
            self.mswave = self.waveImage._build_wave()
        # Save to hard-drive
        if self.save_masters:
            self.waveImage.save_master(self.mswave, steps=self.waveImage.steps)

        # Save & return
        self.calib_dict[self.arc_master_key]['wave'] = self.mswave
        return self.mswave

    def get_wv_calib(self):
        """
        Load or generate the 1D wavelength calibrations

        Requirements:
          msarc, msbpm, tslits_dict
          det, par

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
        self._chk_set(['arc_master_key', 'det', 'calib_ID', 'par'])

        # Return existing data
        prev_build = self.check_for_previous('wavecalib', self.arc_master_key)
        if prev_build:
            self.wv_calib = self.calib_dict[self.arc_master_key]['wavecalib']
            self.wv_maskslits = self.calib_dict[self.arc_master_key]['wvmask']
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
                                             par=self.par['wavelengths'], master_key=self.arc_master_key, master_dir=self.master_dir,
                                             mode=self.par['masters'],redux_path=self.redux_path, bpm=self.msbpm)
        # Load from disk (MasterFrame)?
        self.wv_calib = self.waveCalib.master(force=prev_build)
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
        self.calib_dict[self.arc_master_key]['wavecalib'] = self.wv_calib
        self.calib_dict[self.arc_master_key]['wvmask'] = self.wv_maskslits
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
           det, par, arc_master_key, spectrograph

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
        self._chk_set(['arc_master_key', 'det', 'calib_ID', 'par'])

        # Return existing data
        prev_build = self.check_for_previous('tilts_dict', self.arc_master_key)
        if prev_build:
            self.tilts_dict = self.calib_dict[self.arc_master_key]['tilts_dict']
            self.wt_maskslits = self.calib_dict[self.arc_master_key]['wtmask']
            self.maskslits += self.wt_maskslits
            return self.tilts_dict, self.maskslits

        # Instantiate
        self.waveTilts = wavetilts.WaveTilts(self.msarc, self.tslits_dict, spectrograph=self.spectrograph,
                                             par=self.par['tilts'], wavepar = self.par['wavelengths'], det=self.det,
                                             master_key=self.arc_master_key, master_dir=self.master_dir,
                                             mode=self.par['masters'],
                                             redux_path=self.redux_path, bpm=self.msbpm)
        # Master
        self.tilts_dict = self.waveTilts.master(force=prev_build)
        if self.tilts_dict is None:
            # TODO still need to deal with syntax for LRIS ghosts. Maybe we don't need it
            self.tilts_dict, self.wt_maskslits \
                    = self.waveTilts.run(maskslits=self.maskslits,doqa=self.write_qa, show=self.show)
            if self.save_masters:
                self.waveTilts.save_master(self.tilts_dict, steps=self.waveTilts.steps)
        else:
            self.wt_maskslits = np.zeros_like(self.maskslits, dtype=bool)

        # Save & return
        self.calib_dict[self.arc_master_key]['tilts_dict'] = self.tilts_dict
        self.calib_dict[self.arc_master_key]['wtmask'] = self.wt_maskslits
        self.maskslits += self.wt_maskslits
        return self.tilts_dict, self.maskslits

    def run_the_steps(self):
        """
        Run full the full recipe of calibration steps

        Returns:

        """
        for step in self.steps:
            getattr(self, 'get_{:s}'.format(step))()
        msgs.info("Calibration complete!")

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
        txt = '<{:s}: frame={}, det={}, calib_ID={}'.format(self.__class__.__name__,
                                                          self.frame,
                                                          self.det,
                                                          self.calib_ID)
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


