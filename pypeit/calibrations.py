""" Class for guiding calibration object generation in PYPIT
"""
from __future__ import absolute_import, division, print_function

import os
import numpy as np

from abc import ABCMeta

from astropy.table import Table

from pypeit import msgs
from pypeit.core import pixels

from pypeit import arcimage
from pypeit import biasframe
from pypeit import bpmimage
from pypeit import flatfield
from pypeit import traceimage
from pypeit import traceslits
from pypeit import wavecalib
from pypeit import wavetilts
from pypeit import waveimage

from pypeit.core import fsort
from pypeit.core import masters
from pypeit.core import procimg
from pypeit.core import parse


from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
from pypeit.spectrographs.spectrograph import Spectrograph

from pypeit import debugger


class Calibrations(object):
    """
    This class is primarily designed to guide the generation of calibration images
    and objects in PYPIT

    Must be able to instantiate spectrograph.

    .. todo::
        Improve docstring...

    Parameters
    ----------

    Attributes
    ----------

    Inherited Attributes
    --------------------
    """
    __metaclass__ = ABCMeta

    # TODO: master_root is a bit of a kludge.  It could be defined
    # earlier and/or in par.
    def __init__(self, fitstbl, spectrograph=None, par=None, master_root=None, save_masters=True,
                 write_qa=True):

        # Check the type of the provided fits table
        if not isinstance(fitstbl, Table):
            msgs.error("fitstbl must be an astropy.Table")

        # Parameters unique to this Object
        self.fitstbl = fitstbl
        self.save_masters = save_masters
        self.write_qa = write_qa

        # Spectrometer class
        if spectrograph is None:
            # Set spectrograph from FITS table instrument header
            # keyword.
            if par is not None and par['rdx']['spectrograph'] != fitstbl['instrume'][0]:
                msgs.error('Specified spectrograph does not match instrument in the fits table!')
            self.spectrograph = load_spectrograph(spectrograph=fitstbl['instrume'][0])
        elif isinstance(spectrograph, str):
            self.spectrograph = load_spectrograph(spectrograph=spectrograph)
        elif isinstance(spectrograph, Spectrograph):
            self.spectrograph = spectrograph
        else:
            raise TypeError('Could not instantiate Spectrograph!')

        # Instantiate the parameters
        # TODO: How far down through the other classes to we propagate
        # the spectrograph defaults as is done here...
        self.par = self.spectrograph.default_pypeit_par()['calibrations'] if par is None else par
        if not isinstance(self.par, pypeitpar.CalibrationsPar):
            raise TypeError('Input parameters must be a CalibrationsPar instance.')

        # TODO: Kludge for now:  If master_root is provided, it
        # over-rides the default in CalibrationsPar.  Should just make
        # caldir and master_root the same...
        self.master_root = os.path.join(os.getcwd(), self.par['caldir']) \
                                if master_root is None else master_root
        self.master_dir = self.master_root+'_'+self.spectrograph.spectrograph

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
        self.mbpm = None
        self.pixlocn = None
        self.tslits_dict = None
        self.maskslits = None
        self.wavecalib = None
        self.mstilts = None
        self.mspixflatnrm = None
        self.slitprof = None
        self.mswave = None
        self.datasec_img = None

    def reset(self, setup, det, sci_ID, par):
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

        # Instantiate with everything needed to generate the image (in case we do)
        self.arcImage = arcimage.ArcImage(self.spectrograph, file_list=[], det=self.det,
                                          par=self.par['arcframe'], setup=self.setup,
                                          root_path=self.master_root, mode=self.par['masters'],
                                          fitstbl=self.fitstbl, sci_ID=self.sci_ID,
                                          msbias=self.msbias)
        
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

        # Instantiate
        self.biasFrame = biasframe.BiasFrame(self.spectrograph, det=self.det,
                                             par=self.par['biasframe'], setup=self.setup,
                                             root_path=self.master_root, mode=self.par['masters'],
                                             fitstbl=self.fitstbl, sci_ID=self.sci_ID)

        # Load the MasterFrame (if it exists and is desired) or the command (e.g. 'overscan')
        self.msbias = self.biasFrame.master()
        if self.msbias is None:  # Build it and save it
            self.msbias = self.biasFrame.build_image()
            if self.save_masters:
                self.biasFrame.save_master(self.msbias, raw_files=self.biasFrame.file_list,
                                           steps=self.biasFrame.steps)

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

        # TODO: This determination of scidx is done often.  It
        # should get moved to a function.
        scidx = np.where((self.fitstbl['sci_ID'] == self.sci_ID) \
                            & self.fitstbl['science'])[0][0]

        example_file = os.path.join(self.fitstbl['directory'][scidx],
                                    self.fitstbl['filename'][scidx])
        # Always use the example file
        bpmImage = bpmimage.BPMImage(self.spectrograph, filename=example_file, det=self.det,
                                     msbias=self.msbias if self.par['badpix'] == 'bias' else None,
                                     trim=self.par['trim'])
        # Build, save, and return
        self.msbpm = bpmImage.build()
        self.calib_dict[self.setup]['bpm'] = self.msbpm
        return self.msbpm

    def get_datasec_img(self):
        """
        Generate the datasec image

        Requirements:
           det, sci_ID, par

        Returns:
            self.datasec_img: ndarray

        """
        # Check internals
        self._chk_set(['det', 'sci_ID', 'par'])
        # Get an example science frame file name
        scidx = np.where((self.fitstbl['sci_ID'] == self.sci_ID) & self.fitstbl['science'])[0][0]
        scifile = os.path.join(self.fitstbl['directory'][scidx], self.fitstbl['filename'][scidx])
        # Generate the spectrograph-specific amplifier ID image
        self.datasec_img = self.spectrograph.get_datasec_img(scifile, self.det)
        if self.par['trim']:
            self.datasec_img = procimg.trim_frame(self.datasec_img, self.datasec_img < 1)
        return self.datasec_img

    def get_pixflatnrm(self):
        """
        Load or generate a normalized pixel flat
          and slit profile

        Requirements:
           tslits_dict
           mstilts
           datasec_img
           det, sci_ID, par, setup

        Returns:
            self.mspixflatnrm: ndarray
            self.slitprof: ndarray

        """
        if self.par['flatfield']['method'] is None:
            # User does not want to flat-field
            self.mspixflatnrm = None
            self.slitprof = None
            return self.mspixflatnrm, self.slitprof

        # Check for existing data
        if not self._chk_objs(['tslits_dict', 'mstilts', 'datasec_img']):
            # User cannot flat-field
            self.mspixflatnrm = None
            self.slitprof = None
            return self.mspixflatnrm, self.slitprof
       
        # Check internals
        self._chk_set(['setup', 'det', 'sci_ID', 'par'])

        # Return already generated data
        if np.all([ k in self.calib_dict[self.setup].keys() 
                            for k in ['normpixelflat','slitprof']]):
            self.mspixflatnrm = self.calib_dict[self.setup]['normpixelflat']
            self.slitprof = self.calib_dict[self.setup]['slitprof']
            return self.mspixflatnrm, self.slitprof

        # Instantiate
        pixflat_image_files = fsort.list_of_files(self.fitstbl, 'pixelflat', self.sci_ID)
        self.flatField = flatfield.FlatField(self.spectrograph, file_list=pixflat_image_files,
                                             det=self.det, par=self.par['pixelflatframe'],
                                             setup=self.setup, root_path=self.master_root,
                                             mode=self.par['masters'],
                                             flatpar=self.par['flatfield'], msbias=self.msbias,
                                             tslits_dict=self.tslits_dict, tilts=self.mstilts)

        # Load from disk (MasterFrame)?
        self.mspixflatnrm = self.flatField.master()
        self.slitprof = None

        # Load user supplied flat (e.g. LRISb with pixel flat)?
        if self.par['flatfield']['frame'] not in ['pixelflat', 'trace']:
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
            self.mspixflatnrm, head, _ = masters._load(mspixelflat_name, exten=self.det,
                                                         frametype=None, force=True)
            # TODO -- Handle slitprof properly, i.e.g from a slit flat for LRISb
            self.slitprof = np.ones_like(self.mspixflatnrm)

        if self.mspixflatnrm is None:
            # TODO -- Consider turning the following back on.  I'm regenerating the flat for now
            # Use mstrace if the indices are identical
            #if np.all(sort.ftype_indices(fitstbl,'trace',1) ==
            #                  sort.ftype_indices(fitstbl, 'pixelflat', 1))
            #            and (traceSlits.mstrace is not None):
            #    flatField.mspixelflat = traceSlits.mstrace.copy()
            # Run
            self.mspixflatnrm, self.slitprof = self.flatField.run(armed=False)
            # Save to Masters
            if self.save_masters:
                self.flatField.save_master(self.mspixflatnrm, raw_files=pixflat_image_files,
                                           steps=self.flatField.steps)
                self.flatField.save_master(self.slitprof, raw_files=pixflat_image_files,
                                           steps=self.flatField.steps,
                                           outfile=masters.master_name('slitprof', self.setup,
                                                        self.master_dir))
        elif self.slitprof is None:
            self.slitprof, _, _ = self.flatField.load_master_slitprofile()

        # Save & return
        self.calib_dict[self.setup]['normpixelflat'] = self.mspixflatnrm
        self.calib_dict[self.setup]['slitprof'] = self.slitprof
        return self.mspixflatnrm, self.slitprof

    def get_slits(self, arms=True):
        """
        Load or generate a normalized pixel flat
        First, a trace flat image is generated

        .. todo::
            - arms is a parameter passed to traceSlits.  This may need
              to change if/when arms.py is replaced.

        Requirements:
           pixlocn
           datasec_img
           det par setup sci_ID

        Returns:
            self.tslits_dict
            self.maskslits

        """
        # Check for existing data
        if not self._chk_objs(['pixlocn', 'datasec_img', 'msbpm']):
            self.tslits_dict = None
            self.maskslits = None
            return self.tslits_dict, self.maskslits

        # Check internals
        self._chk_set(['setup', 'det', 'sci_ID', 'par'])

        # Return already generated data
        if 'trace' in self.calib_dict[self.setup].keys():
            self.tslits_dict = self.calib_dict[self.setup]['trace']
            self.maskslits = np.zeros(self.tslits_dict['lcen'].shape[1], dtype=bool)
            return self.tslits_dict, self.maskslits
                
        # Instantiate (without mstrace)
        self.traceSlits = traceslits.TraceSlits(None, self.pixlocn, par=self.par['slits'],
                                                det=self.det, setup=self.setup,
                                                directory_path=self.master_dir,
                                                mode=self.par['masters'], binbpx=self.msbpm)

        # Load via master, as desired
        if not self.traceSlits.master():
            # Build the trace image first
            trace_image_files = fsort.list_of_files(self.fitstbl, 'trace', self.sci_ID)
            traceImage = traceimage.TraceImage(self.spectrograph,
                                           file_list=trace_image_files, det=self.det,
                                           par=self.par['traceframe'])
            # Load up and get ready
            self.traceSlits.mstrace = traceImage.process(bias_subtract=self.msbias,
                                                         trim=self.par['trim'], apply_gain=True)
            _ = self.traceSlits.make_binarr()

            # Compute the plate scale in arcsec which is needed to trim short slits
            scidx = np.where((self.fitstbl['sci_ID'] == self.sci_ID) & self.fitstbl['science'])[0][0]
            binspatial, binspectral = parse.parse_binning(self.fitstbl['binning'][scidx])
            plate_scale = binspatial*self.spectrograph.detector[self.det-1]['platescale']
            # Now we go forth
            self.tslits_dict = self.traceSlits.run(arms=arms, plate_scale = plate_scale)
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
           mstilts, tslits_dict, wv_calib, maskslits
           det, par, setup

        Returns:
            self.mswave: ndarray

        """
        # Check for existing data
        if not self._chk_objs(['mstilts','tslits_dict','wv_calib','maskslits']):
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
            self.mswave = self.mstilts * (self.mstilts.shape[0]-1.0)
            self.calib_dict[self.setup]['wave'] = self.mswave
            return self.mswave

        # Instantiate
        self.waveImage = waveimage.WaveImage(self.tslits_dict['slitpix'],
                                             self.mstilts, self.wv_calib,
                                             setup=self.setup, directory_path=self.master_dir,
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
          msarc, tslits_dict, pixlocn
          det, par, setup, sci_ID

        Returns:
            self.wv_calib: dict
            self.maskslits -- Updated
        """
        # Check for existing data
        if not self._chk_objs(['msarc', 'tslits_dict', 'pixlocn']):
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
        self.waveCalib = wavecalib.WaveCalib(self.msarc, spectrograph=self.spectrograph,
                                             par=self.par['wavelengths'], det=self.det,
                                             setup=self.setup, root_path=self.master_root,
                                             mode=self.par['masters'], fitstbl=self.fitstbl,
                                             sci_ID=self.sci_ID)
        # Load from disk (MasterFrame)?
        self.wv_calib = self.waveCalib.master()
        # Build?
        if self.wv_calib is None:
            self.wv_calib, _ = self.waveCalib.run(self.tslits_dict['lcen'],
                                                  self.tslits_dict['rcen'], self.pixlocn,
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
            self.mstilts: ndarray (2D)
            self.maskslits: ndarray

        """
        # Check for existing data
        if not self._chk_objs(['msarc', 'tslits_dict', 'pixlocn', 'wv_calib', 'maskslits']):
            msgs.error('dont have all the objects')
            self.mstilts = None
            self.wt_maskslits = np.zeros_like(self.maskslits, dtype=bool)
            self.maskslits += self.wt_maskslits
            return self.mstilts, self.maskslits

        # Check internals
        self._chk_set(['setup', 'det', 'sci_ID', 'par'])

        # Return existing data
        if 'tilts' in self.calib_dict[self.setup].keys():
            self.mstilts = self.calib_dict[self.setup]['tilts']
            self.wt_maskslits = self.calib_dict[self.setup]['wtmask']
            self.maskslits += self.wt_maskslits
            return self.mstilts, self.maskslits

        # Instantiate
        self.waveTilts = wavetilts.WaveTilts(self.msarc, spectrograph=self.spectrograph,
                                             par=self.par['tilts'], det=self.det,
                                             setup=self.setup, root_path=self.master_root,
                                             mode=self.par['masters'], pixlocn=self.pixlocn,
                                             tslits_dict=self.tslits_dict)
        # Master
        self.mstilts = self.waveTilts.master()
        if self.mstilts is None:
            self.mstilts, self.wt_maskslits \
                    = self.waveTilts.run(maskslits=self.maskslits, wv_calib=self.wv_calib,
                                         doqa=self.write_qa)
            if self.save_masters:
                self.waveTilts.save_master()
        else:
            self.wt_maskslits = np.zeros_like(self.maskslits, dtype=bool)

        # Save & return
        self.calib_dict[self.setup]['tilts'] = self.mstilts
        self.calib_dict[self.setup]['wtmask'] = self.wt_maskslits
        self.maskslits += self.wt_maskslits
        return self.mstilts, self.maskslits

    def run_the_steps(self):
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
    def __init__(self, fitstbl, spectrograph=None, par=None, master_root=None, save_masters=True,
                 write_qa=True, steps=None):
        Calibrations.__init__(self, fitstbl, spectrograph=spectrograph, par=par,
                              master_root=master_root, save_masters=save_masters,
                              write_qa=write_qa)
        self.steps = MultiSlitCalibrations.default_steps() if steps is None else steps

    @staticmethod
    def default_steps():
        return ['datasec_img', 'bias', 'arc', 'bpm', 'pixlocn', 'slits', 'wv_calib', 'tilts',
                'pixflatnrm', 'wave']

