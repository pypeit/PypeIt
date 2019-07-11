"""
Class for guiding calibration object generation in PypeIt

.. _numpy.ndarray: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html

"""
import os
import numpy as np

from abc import ABCMeta

from astropy.io import fits
from astropy.table import Table
from IPython import embed

from pypeit import msgs
from pypeit import masterframe
from pypeit import arcimage
from pypeit import tiltimage
from pypeit import biasframe
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
from pypeit.spectrographs.spectrograph import Spectrograph

class Calibrations(object):
    """
    This class is primarily designed to guide the generation of
    calibration images and objects in PypeIt.

    To avoid rebuilding MasterFrames that were generated during this execution
    of PypeIt, the class performs book-keeping of these master frames and
    holds that info in self.calib_dict

    Args:
        fitstbl (:class:`pypeit.metadata.PypeItMetaData`):
            The class holding the metadata for all the frames in this
            PypeIt run.
        par (:class:`pypeit.par.pypeitpar.PypeItPar`):
            Parameter set defining optional parameters of PypeIt's
            low-level algorithms.  Needs to specifically be a
            CalibrationsPar child.
        spectrograph (:obj:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph object
        caldir (:obj:`str`, optional):
            Path to write the output calibrations.  If None, calibration
            data are not saved.
        qadir (:obj:`str, optional):
            Path for quality assessment output.  If not provided, no QA
            plots are saved.
        save_masters (:obj:`bool`, optional):
            Save the calibration frames to disk.
        reuse_masters (:obj:`bool`, optional):
            Load calibration files from disk if they exist
        show (:obj:`bool`, optional):
            Show plots of PypeIt's results as the code progesses.
            Requires interaction from the users.

    Attributes:
        TODO: Fix these
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
        calib_ID (:obj:`int`):
            calib group ID of the current frame

    """
    __metaclass__ = ABCMeta

    # TODO: I added back save_masters as a parameter because if you
    # provide a caldir, you may just want to be reusing the masters.  I
    # think the code won't save masters if they're reused, but allowing
    # save_masters as an argument allows us to make this explicit.
    def __init__(self, fitstbl, par, spectrograph, caldir=None, qadir=None, save_masters=True,
                 reuse_masters=False, show=False):

        # Check the types
        if not isinstance(fitstbl, PypeItMetaData):
            msgs.error('fitstbl must be an PypeItMetaData object')
        if not isinstance(par, pypeitpar.CalibrationsPar):
            msgs.error('Input parameters must be a CalibrationsPar instance.')
        if not isinstance(spectrograph, Spectrograph):
            msgs.error('Must provide Spectrograph instance to Calibrations.')

        # Required inputs
        self.fitstbl = fitstbl
        self.par = par
        self.spectrograph = spectrograph

        # Control flow
        self.reuse_masters = reuse_masters
        self.master_dir = caldir
        self.save_masters = save_masters
        self.qa_path = qadir
        self.write_qa = qadir is not None
        self.show = show

        # Check that the masters can be reused and/or saved
        if self.master_dir is None:
            if self.save_masters:
                # TODO: Default to current directory instead?
                msgs.warn('To save masters, must provide the directory (caldir).  '
                          'Masters will not be saved!')
                self.save_masters = False
            if self.reuse_masters:
                # TODO: Default to current directory instead?
                msgs.warn('To reuse masters, must provide the directory (caldir).  '
                          'Masters will not be reused!')
                self.reuse_masters = False

        # Check the directories exist
        # TODO: This should be done when the masters are saved
        if self.save_masters and not os.path.isdir(self.master_dir):
            os.makedirs(self.master_dir)
        # TODO: This should be done when the qa plots are saved
        if self.write_qa and not os.path.isdir(os.path.join(self.qa_path, 'PNGs')):
            os.makedirs(os.path.join(self.qa_path, 'PNGs'))

        # Attributes
        self.calib_dict = {}
        self.det = None
        self.frame = None
        self.binning = None

        # Steps
        self.steps = []

        # Internals
        self._reset_internals()

    def _reset_internals(self):
        """
        Reset all of the key internals to None or an empty object

        """
        # TODO: I think all of these should change to the relevant class
        # objects, particularly if we're going to allow a class (e.g.,
        # FlatField) to alter the internals of another classe (e.g.,
        # TraceSlits)
        self.shape = None
        self.msarc = None
        self.msbias = None
        self.msbpm = None
        self.mstrace = None
        self.tslits_dict = None
        self.wavecalib = None
        self.tilts_dict = None
        self.mspixelflat = None
        self.msillumflat = None
        self.mswave = None
        self.calib_ID = None
        self.master_key_dict = {}

    def _update_cache(self, master_key, master_type, data):
        """
        Update or add new cached data held in memory.

        Fundamentally just executes::

            self.calib_dict[self.master_key_dict[master_key]][master_type] = data

        Args:
            master_key (:obj:`str`):
                Keyword used to select the master key from
                :attr:`master_key_dict`.
            master_type (:obj:`str`, :obj:`tuple`):
                One or more keywords setting the type of master frame
                being saved to :attr:`calib_dict`.  E.g.
                `master_type=bpm` for the data saved to
                `self.calib_dict['A_01_1']['bpm'].
            data (object, :obj:`tuple`):
                One or more data objects to save to :attr:`calib_dict`.
        """
        # Handle a single entry
        _master_type = master_type if isinstance(master_type, tuple) else (master_type,)
        _data = data if isinstance(data, tuple) else (data,)
        # Save the data
        # TODO: Allow for copy option?  Do we now whether or not this is
        # actually being done correctly as it is?
        for key, d in zip(_master_type, _data):
            self.calib_dict[self.master_key_dict[master_key]][key] = d

    def _cached(self, master_type, master_key):
        """
        Check if the calibration frame data has been cached in memory.

        The image data is saved in memory as, e.g.::

           self.calib_dict[master_key][master_type] = {}

        If the master has not yet been generated, an empty dict is
        prepared::

           self.calib_dict[master_key] = {}

        Args:
            master_type (str): Type of calibration frame
            master_key (str): Master key naming

        Returns:
             bool: True = Built previously
        """
        if master_key not in self.calib_dict.keys():
            # Masters are not available for any frames in this
            # configuration + calibration group + detector
            self.calib_dict[master_key] = {}
            return False
        if master_key in self.calib_dict.keys():
            if master_type in self.calib_dict[master_key].keys():
                # Found the previous master in memory
                msgs.info('Using {0} for {1} found in cache.'.format(master_type, master_key))
                return True
        # Master key exists but no master in memory for this specific
        # type
        self.calib_dict[master_key][master_type] = {}
        return False

    def set_config(self, frame, det, par=None):
        """
        Specify the parameters of the Calibrations class and reset all
        the internals to None. The internal dict is left unmodified.

        Args:
            frame (int): Frame index in the fitstbl
            det (int): Detector number
            par (:class:`pypeit.par.pypeitpar.CalibrationPar`):

        """
        # Reset internals to None
        # NOTE: This sets empties calib_ID and master_key_dict so must
        # be done here first before these things are initialized below.
        self._reset_internals()

        # Initialize for this setup
        self.frame = frame
        self.calib_ID = int(self.fitstbl['calib'][frame])
        self.det = det
        if par is not None:
            self.par = par
        # Deal with binning
        self.binning = self.fitstbl['binning'][self.frame]

        # Initialize the master key dict for this science/standard frame
        self.master_key_dict['frame'] = self.fitstbl.master_key(frame, det=det)

    def get_arc(self):
        """
        Load or generate the Arc image

        Requirements:
          master_key, det, par

        Args:

        Returns:
            ndarray: :attr:`msarc` image

        """
        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        arc_rows = self.fitstbl.find_frames('arc', calib_ID=self.calib_ID, index=True)
        self.arc_files = self.fitstbl.frame_paths(arc_rows)
        self.master_key_dict['arc'] \
                = self.fitstbl.master_key(arc_rows[0] if len(arc_rows) > 0 else self.frame,
                                          det=self.det)

        if self._cached('arc', self.master_key_dict['arc']):
            # Previously calculated
            self.msarc = self.calib_dict[self.master_key_dict['arc']]['arc']
            return self.msarc

        # Instantiate with everything needed to generate the image (in case we do)
        self.arcImage = arcimage.ArcImage(self.spectrograph, files=self.arc_files,
                                          det=self.det, msbias=self.msbias,
                                          par=self.par['arcframe'],
                                          master_key=self.master_key_dict['arc'],
                                          master_dir=self.master_dir,
                                          reuse_masters=self.reuse_masters)

        # Load the MasterFrame (if it exists and is desired)?
        self.msarc = self.arcImage.load()
        if self.msarc is None:  # Otherwise build it
            msgs.info("Preparing a master {0:s} frame".format(self.arcImage.frametype))
            self.msarc = self.arcImage.build_image(bias=self.msbias, bpm=self.msbpm)
            # Save to Masters
            if self.save_masters:
                self.arcImage.save()

        # Save & return
        self._update_cache('arc', 'arc', self.msarc)
        return self.msarc


    def get_tiltimg(self):
        """
        Load or generate the Tilt image

        Requirements:
          master_key, det, par

        Args:

        Returns:
            ndarray: :attr:`mstilt` image

        """
        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        tilt_rows = self.fitstbl.find_frames('tilt', calib_ID=self.calib_ID, index=True)
        if len(tilt_rows) == 0:
            msgs.error('Must identify tilt frames to construct tilt image.')
        self.tilt_files = self.fitstbl.frame_paths(tilt_rows)
        self.master_key_dict['tilt'] \
                = self.fitstbl.master_key(tilt_rows[0] if len(tilt_rows) > 0 else self.frame,
                                          det=self.det)

        if self._cached('tiltimg', self.master_key_dict['tilt']):
            # Previously calculated
            self.mstilt = self.calib_dict[self.master_key_dict['tilt']]['tiltimg']
            return self.mstilt

        # Instantiate with everything needed to generate the image (in case we do)
        self.tiltImage = tiltimage.TiltImage(self.spectrograph, files=self.tilt_files,
                                          det=self.det, msbias=self.msbias,
                                          par=self.par['tiltframe'],
                                          master_key=self.master_key_dict['tilt'],
                                          master_dir=self.master_dir,
                                          reuse_masters=self.reuse_masters)

        # Load the MasterFrame (if it exists and is desired)?
        self.mstilt = self.tiltImage.load()
        if self.mstilt is None:  # Otherwise build it
            msgs.info("Preparing a master {0:s} frame".format(self.tiltImage.frametype))
            self.mstilt = self.tiltImage.build_image(bias=self.msbias, bpm=self.msbpm)
            # JFH Add a cr_masking option here. The image processing routines are not ready for it yet.

            # Save to Masters
            if self.save_masters:
                self.tiltImage.save()

        # Save & return
        self._update_cache('tilt', 'tiltimg', self.mstilt)
        # TODO in the future add in a tilt_inmask
        #self._update_cache('tilt', 'tilt_inmask', self.mstilt_inmask)
        return self.mstilt

    def get_bias(self):
        """
        Load or generate the bias frame/command

        Requirements:
           master_key, det, par

        Returns:
            ndarray or str: :attr:`bias`

        """

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        bias_rows = self.fitstbl.find_frames('bias', calib_ID=self.calib_ID, index=True)
        self.bias_files = self.fitstbl.frame_paths(bias_rows)

        self.master_key_dict['bias'] \
                = self.fitstbl.master_key(bias_rows[0] if len(bias_rows) > 0 else self.frame,
                                          det=self.det)

        # Grab from internal dict (or hard-drive)?
        if self._cached('bias', self.master_key_dict['bias']):
            self.msbias = self.calib_dict[self.master_key_dict['bias']]['bias']
            msgs.info("Reloading the bias from the internal dict")
            return self.msbias

        # Instantiate
        self.biasFrame = biasframe.BiasFrame(self.spectrograph, files=self.bias_files,
                                             det=self.det, par=self.par['biasframe'],
                                             master_key=self.master_key_dict['bias'],
                                             master_dir=self.master_dir,
                                             reuse_masters=self.reuse_masters)

        # Try to load the master bias
        self.msbias = self.biasFrame.load()
        if self.msbias is None:
            # Build it and save it
            self.msbias = self.biasFrame.build_image()
            if self.save_masters:
                self.biasFrame.save()

        # Save & return
        self._update_cache('bias', 'bias', self.msbias)
        return self.msbias

    def get_bpm(self):
        """
        Load or generate the bad pixel mask

        TODO -- Should consider doing this outside of calibrations as it is
        more specific to the science frame

        This needs to be for the *trimmed* and correctly oriented image!

        Requirements:
           Instrument dependent

        Returns:
            ndarray: :attr:`msbpm` image of bad pixel mask

        """
        # Check internals
        self._chk_set(['par', 'det'])

        # Generate a bad pixel mask (should not repeat)
        self.master_key_dict['bpm'] = self.fitstbl.master_key(self.frame, det=self.det)

        if self._cached('bpm', self.master_key_dict['bpm']):
            self.msbpm = self.calib_dict[self.master_key_dict['bpm']]['bpm']
            return self.msbpm

        # Build the data-section image
        sci_image_file = self.fitstbl.frame_paths(self.frame)
        rdsec_img = self.spectrograph.get_rawdatasec_img(sci_image_file, det=self.det)

        # Instantiate the shape here, based on the shape of the science
        # image. This is the shape of most calibrations, although we are
        # allowing for arcs of different shape because of X-shooter etc.
        trim = procimg.trim_frame(rdsec_img, rdsec_img < 1)
        orient = self.spectrograph.orient_image(trim, self.det)
        self.shape = orient.shape

        # Build it
        self.msbpm = self.spectrograph.bpm(shape=self.shape, filename=sci_image_file, det=self.det)

        # Record it
        self._update_cache('bpm', 'bpm', self.msbpm)
        # Return
        return self.msbpm

    def get_flats(self):
        """
        Load or generate a normalized pixel flat and slit illumination
        flat.

        Requires :attr:`tslits_dict`, :attr:`tilts_dict`, :attr:`det`,
        :attr:`par`.

        Returns:
            `numpy.ndarray`_: Two arrays are returned, the normalized
            pixel flat image (:attr:`mspixelflat`) and the slit
            illumination flat (:attr:`msillumflat`).  If the user
            requested the field flattening be skipped
            (`FlatFieldPar['method'] == 'skip'`) or if the slit and tilt
            traces are not provided, the function returns two None
            objects instead.
        """
        # Check for existing data
        if not self._chk_objs(['msarc', 'msbpm', 'tslits_dict', 'wv_calib']):
            msgs.error('dont have all the objects')

        if self.par['flatfield']['method'] is 'skip':
            # User does not want to flat-field
            self.mspixelflat = None
            self.msillumflat = None
            msgs.warning('Parameter calibrations.flatfield.method is set to skip. You are NOT '
                         'flatfielding your data!!!')
            # TODO: Why does this not return unity arrays, like what's
            # done below?
            return self.mspixelflat, self.msillumflat

        # Slit and tilt traces are required to flat-field the data
        if not self._chk_objs(['tslits_dict', 'tilts_dict']):
            msgs.warning('Flats were requested, but there are quantities missing necessary to '
                         'create flats.  Proceeding without flat fielding....')
            # User cannot flat-field
            self.mspixelflat = None
            self.msillumflat = None
            # TODO: Why does this not return unity arrays, like what's
            # done below?
            return self.mspixelflat, self.msillumflat

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        pixflat_rows = self.fitstbl.find_frames('pixelflat', calib_ID=self.calib_ID, index=True)
        # TODO: Why aren't these set to self
        #   KBW: They're kept in self.flatField.files
        pixflat_image_files = self.fitstbl.frame_paths(pixflat_rows)
        # Allow for user-supplied file (e.g. LRISb)
        self.master_key_dict['flat'] \
                = self.fitstbl.master_key(pixflat_rows[0] if len(pixflat_rows) > 0 else self.frame,
                                          det=self.det)

        # Return already generated data
        if self._cached('pixelflat', self.master_key_dict['flat']) \
                and self._cached('illumflat', self.master_key_dict['flat']):
            self.mspixelflat = self.calib_dict[self.master_key_dict['flat']]['pixelflat']
            self.msillumflat = self.calib_dict[self.master_key_dict['flat']]['illumflat']
            return self.mspixelflat, self.msillumflat

        # Instantiate
        # TODO: This should automatically attempt to load and instatiate
        # from a file if it exists.
        self.flatField = flatfield.FlatField(self.spectrograph, self.par['pixelflatframe'],
                                             files=pixflat_image_files, det=self.det,
                                             master_key=self.master_key_dict['flat'],
                                             master_dir=self.master_dir,
                                             reuse_masters=self.reuse_masters,
                                             flatpar=self.par['flatfield'],
                                             msbias=self.msbias,
                                             # TODO: msbpm is not passed?
                                             tslits_dict=self.tslits_dict,
                                             tilts_dict=self.tilts_dict)

        # --- Pixel flats

        # 1)  Try to load master files from disk (MasterFrame)?
        _, self.mspixelflat, self.msillumflat = self.flatField.load()

        # 2) Did the user specify a flat? If so load it in  (e.g. LRISb with pixel flat)?
        # TODO: We need to document this format for the user!
        if self.par['flatfield']['frame'] != 'pixelflat':
            # - Name is explicitly correct?
            if os.path.isfile(self.par['flatfield']['frame']):
                flat_file = self.par['flatfield']['frame']
            # - Is it in the master directory?
            elif os.path.isfile(os.path.join(self.flatField.master_dir,
                                             self.par['flatfield']['frame'])):
                flat_file = os.path.join(self.flatField.master_dir, self.par['flatfield']['frame'])
            else:
                msgs.error('Could not find user-defined flatfield file: {0}'.format(
                           self.par['flatfield']['frame']))
            msgs.info('Using user-defined file: {0}'.format(flat_file))
            with fits.open(flat_file) as hdu:
                self.mspixelflat = hdu[self.det].data
            self.msillumflat = None

        # 3) there is no master or no user supplied flat, generate the flat
        if self.mspixelflat is None and len(pixflat_image_files) != 0:
            # Run
            self.mspixelflat, self.msillumflat = self.flatField.run(show=self.show,
                                                                    maskslits=self.tslits_dict['maskslits'])

            # If we tweaked the slits, update the tilts_dict and
            # tslits_dict to reflect new slit edges
            if self.par['flatfield']['tweak_slits']:
                msgs.info('Using slit boundary tweaks from IllumFlat and updated tilts image')
                self.tslits_dict = self.flatField.tslits_dict
                self.tilts_dict = self.flatField.tilts_dict

            # Save to Masters
            if self.save_masters:
                self.flatField.save()

                # If we tweaked the slits update the master files for tilts and slits
                # TODO: These should be saved separately
                if self.par['flatfield']['tweak_slits']:
                    msgs.info('Updating MasterTrace and MasterTilts using tweaked slit boundaries')
                    # Add tweaked boundaries to the MasterTrace file
                    self.traceSlits.tslits_dict = self.flatField.tslits_dict
                    try:
                        self.traceSlits.save(traceImage=self.traceImage)
                    except:
                        self.traceSlits.save(traceImage=self.mstrace)
                    # Write the final_tilts using the new slit boundaries to the MasterTilts file
                    self.waveTilts.final_tilts = self.flatField.tilts_dict['tilts']
                    self.waveTilts.tilts_dict = self.flatField.tilts_dict
                    self.waveTilts.save()

        # 4) If either of the two flats are still None, use unity
        # everywhere and print out a warning
        # TODO: These will barf if self.tilts_dict['tilts'] isn't
        # defined.
        if self.mspixelflat is None:
            self.mspixelflat = np.ones_like(self.tilts_dict['tilts'])
            msgs.warn('You are not pixel flat fielding your data!!!')
        if self.msillumflat is None:
            self.msillumflat = np.ones_like(self.tilts_dict['tilts'])
            msgs.warn('You are not illumination flat fielding your data!')

        # Save & return
        self._update_cache('flat', ('pixelflat','illumflat'), (self.mspixelflat,self.msillumflat))
        return self.mspixelflat, self.msillumflat

    # TODO: if write_qa need to provide qa_path!
    # TODO: why do we allow redo here?
    def get_slits(self, redo=False, write_qa=True):
        """
        Load or generate the definition of the slit boundaries.

        Internals that must be available are :attr:`fitstbl`,
        :attr:`calib_ID`, :attr:`det`.

        Args:
            redo (bool): Redo
            write_qa (bool, optional):
              Generate the QA?  Turn off for testing..

        Returns:
            Returns the trace-slits dictionary (also kept internally as
            :attr:`tslits_dict`) and the slit mask array (numpy.ndarray;
            also kept internally as :attr:`maskslits`)

        """
        # Check for existing data
        if not self._chk_objs(['msbpm']):
            self.tslits_dict = None
            return self.tslits_dict

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        trace_rows = self.fitstbl.find_frames('trace', calib_ID=self.calib_ID, index=True)
        self.trace_image_files = self.fitstbl.frame_paths(trace_rows)
        self.master_key_dict['trace'] \
                = self.fitstbl.master_key(trace_rows[0] if len(trace_rows) > 0 else self.frame,
                                          det=self.det)

        # Return already generated data
        if self._cached('trace', self.master_key_dict['trace']) and not redo:
            self.tslits_dict = self.calib_dict[self.master_key_dict['trace']]['trace']
            return self.tslits_dict

        # Instantiate
        self.traceSlits = traceslits.TraceSlits(self.spectrograph, self.par['slits'], det=self.det,
                                                master_key=self.master_key_dict['trace'],
                                                master_dir=self.master_dir, qa_path=self.qa_path,
                                                reuse_masters=self.reuse_masters, msbpm=self.msbpm)

        # Load the MasterFrame (if it exists and is desired)?
        self.tslits_dict, _ = self.traceSlits.load()
        if self.tslits_dict is None:
            # Build the trace image
            self.traceImage = traceimage.TraceImage(self.spectrograph,
                                                    files=self.trace_image_files, det=self.det,
                                                    par=self.par['traceframe'],
                                                    bias=self.msbias)
            self.traceImage.build_image(bias=self.msbias, bpm=self.msbpm)

            # Compute the plate scale in arcsec which is needed to trim short slits
            binspectral, binspatial = parse.parse_binning(self.binning)
            plate_scale = binspatial*self.spectrograph.detector[self.det-1]['platescale']

            # JFH Why is this stuff on user defined slits here and not
            # in the class?  User-defined slits??
            # TODO: this should be done inside TraceSlits so that the
            # call to run() or whatever has the same format as what the
            # user sees in TraceSlitsPar
            add_user_slits = None if self.par['slits']['add_slits'] is None \
                                else trace_slits.parse_user_slits(self.par['slits']['add_slits'],
                                                                  self.det)
            rm_user_slits = None if self.par['slits']['rm_slits'] is None \
                                else trace_slits.parse_user_slits(self.par['slits']['rm_slits'],
                                                                  self.det, rm=True)
            # Now we go forth
            try:
                self.tslits_dict = self.traceSlits.run(self.traceImage.image,
                                                       self.binning,
                                                       add_user_slits=add_user_slits,
                                                       rm_user_slits=rm_user_slits,
                                                       plate_scale=plate_scale,
                                                       show=self.show,
                                                       write_qa=write_qa)
            except:
                self.traceSlits.save(traceImage=self.traceImage)
                msgs.error('Crashed out of finding the slits. Have saved the work done to disk '
                           'but it needs fixing.')

            # No slits?
            if self.tslits_dict is None:
                return self.tslits_dict

            # Save to disk
            if self.save_masters:
                self.traceSlits.save(traceImage=self.traceImage)
        # Save, initialize maskslits, and return
        # TODO: We're not caching self.mstrace.  And actually there is
        # no mstrace in Calibrations anymore; only in TraceSlits?
        self._update_cache('trace', 'trace', self.tslits_dict)
        return self.tslits_dict

    def get_wave(self):
        """
        Load or generate a wavelength image

        Requirements:
           tilts_dict, tslits_dict, wv_calib
           det, par, master_key

        Returns:
            ndarray: :attr:`mswave` wavelength image

        """
        # Check for existing data
        if not self._chk_objs(['tilts_dict', 'tslits_dict', 'wv_calib']):
            self.mswave = None
            return self.mswave

        # Check internals
        self._chk_set(['det', 'par'])

        # Return existing data
        if self._cached('wave', self.master_key_dict['arc']):
            self.mswave = self.calib_dict[self.master_key_dict['arc']]['wave']
            return self.mswave

        # No wavelength calibration requested
        if self.par['wavelengths']['reference'] == 'pixel':
            msgs.warn('No wavelength calibration performed!')
            self.mswave = self.tilts_dict['tilts'] * (self.tilts_dict['tilts'].shape[0]-1.0)
            self.calib_dict[self.master_key_dict['arc']]['wave'] = self.mswave
            return self.mswave

        # Instantiate
        # TODO we are regenerating this mask a lot in this module. Could reduce that
        self.waveImage = waveimage.WaveImage(self.tslits_dict, self.tilts_dict['tilts'],
                                             self.wv_calib, self.spectrograph, self.det,
                                             self.tslits_dict['maskslits'],
                                             master_key=self.master_key_dict['arc'],
                                             master_dir=self.master_dir,
                                             reuse_masters=self.reuse_masters)

        # Attempt to load master
        self.mswave = self.waveImage.load()
        if self.mswave is None:
            self.mswave = self.waveImage.build_wave()
            # Save to hard-drive
            if self.save_masters:
                self.waveImage.save()

        # Save & return
        self._update_cache('arc', 'wave', self.mswave)
        return self.mswave

    def get_wv_calib(self):
        """
        Load or generate the 1D wavelength calibrations

        Requirements:
          msarc, msbpm, tslits_dict
          det, par

        Returns:
            dict, ndarray: :attr:`wv_calib` calibration dict and the updated slit mask array
        """
        # Check for existing data
        if not self._chk_objs(['msarc', 'msbpm', 'tslits_dict']):
            msgs.error('dont have all the objects')

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])
        if 'arc' not in self.master_key_dict.keys():
            msgs.error('Arc master key not set.  First run get_arc.')

        # Return existing data
        if self._cached('wavecalib', self.master_key_dict['arc']) \
                and self._cached('wvmask', self.master_key_dict['arc']):
            self.wv_calib = self.calib_dict[self.master_key_dict['arc']]['wavecalib']
            self.wv_maskslits = self.calib_dict[self.master_key_dict['arc']]['wvmask']
            self.tslits_dict['maskslits'] += self.wv_maskslits
            return self.wv_calib

        # No wavelength calibration requested
        if self.par['wavelengths']['reference'] == 'pixel':
            msgs.info("A wavelength calibration will not be performed")
            self.wv_calib = None
            self.wv_maskslits = np.zeros_like(self.maskslits, dtype=bool)
            self.tslits_dict['maskslits'] += self.wv_maskslits
            return self.wv_calib

        # Grab arc binning (may be different from science!)
        arc_rows = self.fitstbl.find_frames('arc', calib_ID=self.calib_ID, index=True)
        self.arc_files = self.fitstbl.frame_paths(arc_rows)
        binspec, binspat = parse.parse_binning(self.spectrograph.get_meta_value(self.arc_files[0],
                                                                                'binning'))
        # Instantiate
        self.waveCalib = wavecalib.WaveCalib(self.msarc, self.tslits_dict, self.spectrograph,
                                             self.par['wavelengths'], binspectral=binspec,
                                             det=self.det, master_key=self.master_key_dict['arc'],
                                             master_dir=self.master_dir,
                                             reuse_masters=self.reuse_masters,
                                             qa_path=self.qa_path, msbpm=self.msbpm)
        # Load from disk (MasterFrame)?
        self.wv_calib = self.waveCalib.load()
        if self.wv_calib is None:
            self.wv_calib, _ = self.waveCalib.run(skip_QA=(not self.write_qa))
            # Save to Masters
            if self.save_masters:
                self.waveCalib.save()

        # Create the mask (needs to be done here in case wv_calib was loaded from Masters)
        # TODO: This should either be done here or save as part of the
        # master frame file.  As it is, if not loaded from the master
        # frame file, mask_maskslits is run twice, once in run above and
        # once here...
        self.wv_maskslits = self.waveCalib.make_maskslits(self.tslits_dict['slit_left'].shape[1])
        self.tslits_dict['maskslits'] += self.wv_maskslits

        # Save & return
        self._update_cache('arc', ('wavecalib','wvmask'), (self.wv_calib,self.wv_maskslits))
        # Return
        return self.wv_calib

    def get_tilts(self):
        """
        Load or generate the tilts image

        Requirements:
           mstilt, tslits_dict, wv_calib
           det, par, spectrograph

        Returns:
            dict, ndarray: :attr:`tilts_dict` dictionary with tilts information (2D)
            and the updated slit mask array

        """
        # Check for existing data
        #TODO add mstilt_inmask to this list when it gets implemented.
        if not self._chk_objs(['mstilt', 'msbpm', 'tslits_dict', 'wv_calib']):
            msgs.error('dont have all the objects')
            self.tilts_dict = None
            self.wt_maskslits = np.zeros_like(self.maskslits, dtype=bool)
            self.tslits_dict['maskslits'] += self.wt_maskslits
            return self.tilts_dict

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])
        if 'tilt' not in self.master_key_dict.keys():
            msgs.error('Tilt master key not set.  First run get_tiltimage.')

        # Return existing data
        if self._cached('tilts_dict', self.master_key_dict['tilt']) \
                and self._cached('wtmask', self.master_key_dict['tilt']):
            self.tilts_dict = self.calib_dict[self.master_key_dict['tilt']]['tilts_dict']
            self.wt_maskslits = self.calib_dict[self.master_key_dict['tilt']]['wtmask']
            self.tslits_dict['maskslits'] += self.wt_maskslits
            return self.tilts_dict

        # Instantiate
        self.waveTilts = wavetilts.WaveTilts(self.mstilt, self.tslits_dict, self.spectrograph,
                                             self.par['tilts'], self.par['wavelengths'],
                                             det=self.det, master_key=self.master_key_dict['tilt'],
                                             master_dir=self.master_dir,
                                             reuse_masters=self.reuse_masters,
                                             qa_path=self.qa_path, msbpm=self.msbpm)
        # Master
        self.tilts_dict = self.waveTilts.load()
        if self.tilts_dict is None:
            # TODO still need to deal with syntax for LRIS ghosts. Maybe we don't need it
            self.tilts_dict, self.wt_maskslits \
                    = self.waveTilts.run(maskslits=self.tslits_dict['maskslits'], doqa=self.write_qa,
                                         show=self.show)
            if self.save_masters:
                self.waveTilts.save()
        else:
            self.wt_maskslits = np.zeros_like(self.tslits_dict['maskslits'], dtype=bool)

        # Save & return
        self._update_cache('tilt', ('tilts_dict','wtmask'), (self.tilts_dict,self.wt_maskslits))
        self.tslits_dict['maskslits'] += self.wt_maskslits
        return self.tilts_dict

    def run_the_steps(self):
        """
        Run full the full recipe of calibration steps
        """
        for step in self.steps:
            getattr(self, 'get_{:s}'.format(step))()
        msgs.info("Calibration complete!")

    def _chk_set(self, items):
        """
        Check whether a needed attribute has previously been set

        Args:
            items (list): Attributes to check

        """
        for item in items:
            if getattr(self, item) is None:
                msgs.error("Use self.set to specify '{:s}' prior to generating XX".format(item))

    # This is specific to `self.ms*` attributes
    def _chk_objs(self, items):
        """
        Check that the input items exist internally as attributes

        Args:
            items (list):

        Returns:
            bool: True if all exist

        """
        for obj in items:
            if getattr(self, obj) is None:
                msgs.warn("You need to generate {:s} prior to this calibration..".format(obj))
                if obj in ['tslits_dict']:
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
    calibrations.  See :class:`pypeit.calibrations.Calibrations` for
    arguments.
    """
    def __init__(self, fitstbl, par, spectrograph, caldir=None, qadir=None, reuse_masters=False,
                 show=False, steps=None):
        super(MultiSlitCalibrations, self).__init__(fitstbl, par, spectrograph, caldir=caldir,
                                                    qadir=qadir, reuse_masters=reuse_masters,
                                                    show=show)
        self.steps = MultiSlitCalibrations.default_steps() if steps is None else steps

    @staticmethod
    def default_steps():
        """
        This defines the steps for calibrations and their order

        Returns:
            list: Calibration steps, in order of execution

        """
        # Order matters!
        return ['bpm', 'bias', 'arc', 'tiltimg', 'slits', 'wv_calib', 'tilts', 'flats', 'wave']

    # Notes on order of steps if flexure compensation is implemented
    #  ['bpm', 'bias', 'arc', 'tiltimg', 'slits', 'wv_calib', 'tilts', 'flats', 'wave']