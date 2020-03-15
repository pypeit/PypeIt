"""
Class for guiding calibration object generation in PypeIt

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os

from abc import ABCMeta

from IPython import embed

import numpy as np

from astropy.io import fits

from pypeit import msgs
from pypeit import arcimage
from pypeit import tiltimage
from pypeit import biasframe
from pypeit import flatfield
from pypeit import edgetrace
from pypeit import masterframe
from pypeit import slittrace
from pypeit import traceimage
from pypeit import wavecalib
from pypeit import wavetilts
from pypeit import waveimage
from pypeit.metadata import PypeItMetaData
from pypeit.core import parse
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
        fitstbl (:class:`pypeit.metadata.PypeItMetaData`, None):
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
        qadir (:obj:`str`, optional):
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
        # TODO -- Remove this None option once we have data models for all the Calibrations
        #  outputs and use them to feed Reduce instead of the Calibrations object
        if not isinstance(fitstbl, PypeItMetaData) and fitstbl is not None:
            msgs.error('fitstbl must be an PypeItMetaData object')
        if not isinstance(par, pypeitpar.CalibrationsPar):
            msgs.error('Input parameters must be a CalibrationsPar instance.')
        if not isinstance(spectrograph, Spectrograph):
            msgs.error('Must provide Spectrograph instance to Calibrations.')

        # Required inputs
        self.fitstbl = fitstbl
        self.par = par
        self.spectrograph = spectrograph

        # Masters
        self.reuse_masters = reuse_masters
        self.master_dir = caldir
        self.save_masters = save_masters

        # QA
        self.qa_path = qadir
        self.write_qa = qadir is not None
        self.show = show

        # Check that the masters can be reused and/or saved
        #if caldir is None:
        #    if self.save_masters:
        #        # TODO: Default to current directory instead?
        #        msgs.warn('To save masters, must provide the directory (caldir).  '
        #                  'Masters will not be saved!')
        #        self.save_masters = False
        #    if self.reuse_masters:
        #        # TODO: Default to current directory instead?
        #        msgs.warn('To reuse masters, must provide the directory (caldir).  '
        #                  'Masters will not be reused!')
        #        self.reuse_masters = False

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

    def _prep_calibrations(self, ctype):
        # Grab rows and files
        rows = self.fitstbl.find_frames(ctype, calib_ID=self.calib_ID, index=True)
        image_files = self.fitstbl.frame_paths(rows)
        # Update the internal dict
        #   Kludge for flats
        if ctype == 'pixelflat':
            _ctype = 'flat'
        else:
            _ctype = ctype
        self.master_key_dict[_ctype] \
            = self.fitstbl.master_key(rows[0] if len(rows) > 0 else self.frame, det=self.det)
        # Return
        return image_files

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
        self.slits = None
        self.wavecalib = None
        self.wavetilts = None
        #self.mspixelflat = None
        #self.msillumflat = None
        self.flatimages = None
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
                being saved to :attr:`calib_dict`. E.g.
                ``master_type=bpm`` for the data saved to
                ``self.calib_dict['A_01_1']['bpm']``.
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
        # Master key exists but no master in memory for this specific type
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
        # Initialize the master dict for input, output

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
        arc_files = self._prep_calibrations('arc')
        masterframe_name = masterframe.construct_file_name(
            arcimage.ArcImage, self.master_key_dict['arc'], master_dir=self.master_dir)

        # Previously calculated?  If so, reuse
        if self._cached('arc', self.master_key_dict['arc']):
            self.msarc = self.calib_dict[self.master_key_dict['arc']]['arc']
            return self.msarc

        # Reuse master frame?
        if os.path.isfile(masterframe_name) and self.reuse_masters:
            self.msarc = arcimage.ArcImage.from_master_file(masterframe_name)
        else:  # Build it
            msgs.info("Preparing a master {0:s} frame".format(arcimage.ArcImage.frametype))
            buildArcImage = arcimage.BuildArcImage(self.spectrograph, self.det,
                                                        self.par['arcframe']['process'],
                                                        arc_files,
                                                        bias=self.msbias)
            self.msarc = buildArcImage.build_image(bias=self.msbias, bpm=self.msbpm)

            # Save
            if self.save_masters:
                self.msarc.to_master_file(self.master_dir, self.master_key_dict['arc'],  # Naming
                                          self.spectrograph.spectrograph,  # Header
                                          steps=buildArcImage.process_steps,
                                          raw_files=arc_files)
        # Cache
        self._update_cache('arc', 'arc', self.msarc)
        # Return
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
        tilt_files = self._prep_calibrations('tilt')
        masterframe_name = masterframe.construct_file_name(
            tiltimage.TiltImage, self.master_key_dict['tilt'], master_dir=self.master_dir)

        # Previously calculated?  If so, reuse
        if self._cached('tiltimg', self.master_key_dict['tilt']):
            self.mstilt = self.calib_dict[self.master_key_dict['tilt']]['tiltimg']
            return self.mstilt

        # Reuse master frame?
        if os.path.isfile(masterframe_name) and self.reuse_masters:
            self.mstilt = tiltimage.TiltImage.from_master_file(masterframe_name)
        else: # Build
            msgs.info("Preparing a master {0:s} frame".format(tiltimage.TiltImage.frametype))
            buildtiltImage = tiltimage.BuildTiltImage(self.spectrograph, self.det,
                                                       self.par['tiltframe']['process'],
                                                       tilt_files, bias=self.msbias)
            self.mstilt = buildtiltImage.build_image(bias=self.msbias, bpm=self.msbpm)

            # Save to Masters
            if self.save_masters:
                self.mstilt.to_master_file(self.master_dir, self.master_key_dict['tilt'],  # Naming
                                      self.spectrograph.spectrograph,  # Header
                                      steps=buildtiltImage.process_steps,
                                      raw_files=tilt_files)

        # Cache
        self._update_cache('tilt', 'tiltimg', self.mstilt)
        # TODO in the future add in a tilt_inmask
        #self._update_cache('tilt', 'tilt_inmask', self.mstilt_inmask)

        # Return
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
        bias_files = self._prep_calibrations('bias')

        # Grab from internal dict (or hard-drive)?
        if self._cached('bias', self.master_key_dict['bias']):
            self.msbias = self.calib_dict[self.master_key_dict['bias']]['bias']
            msgs.info("Reloading the bias from the internal dict")
            return self.msbias

        # Instantiate
        self.biasFrame = biasframe.BiasFrame(self.spectrograph, self.par['biasframe'],
                                             self.det, files=bias_files)

        # Construct the name, in case we need it
        masterframe_name = masterframe.construct_file_name(biasframe.BiasImage,
                                                           self.master_key_dict['bias'],
                                                           master_dir=self.master_dir)
        # Try to load the master bias
        self.msbias = self.biasFrame.load(masterframe_name, reuse_masters=self.reuse_masters)
        if self.msbias is None:
            # Build it and save it
            self.msbias = self.biasFrame.build_image()
            if self.save_masters:
                self.msbias.to_master_file(self.master_dir, self.master_key_dict['bias'],  # Naming
                                          self.spectrograph.spectrograph,  # Header
                                          steps=self.biasFrame.process_steps,
                                          raw_files=bias_files)

        # Save & return
        self._update_cache('bias', 'bias', self.msbias)

        return self.msbias

    def get_bpm(self):
        """
        Load or generate the bad pixel mask

        TODO -- Should consider doing this outside of calibrations as it is
        more specific to the science frame - unless we want to generate a BPM
        from the bias frame.

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

        # Check if a bias frame exists, and if a BPM should be generated
        msbias = None
        if self.par['bpm_usebias'] and self._cached('bias', self.master_key_dict['bias']):
            msbias = self.msbias
        # Build it
        self.msbpm = self.spectrograph.bpm(sci_image_file, self.det, msbias=msbias)
        self.shape = self.msbpm.shape

        # Record it
        self._update_cache('bpm', 'bpm', self.msbpm)
        # Return
        return self.msbpm

    def get_flats(self):
        """
        Load or generate a normalized pixel flat and slit illumination
        flat.

        Requires :attr:`slits`, :attr:`wavetilts`, :attr:`det`,
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
        if not self._chk_objs(['msarc', 'msbpm', 'slits', 'wv_calib']):
            msgs.error('Must have the arc, bpm, slits, and wv_calib defined to proceed!')

        if self.par['flatfield']['method'] is 'skip':
            # User does not want to flat-field
            msgs.warn('Parameter calibrations.flatfield.method is set to skip. You are NOT '
                      'flatfielding your data!!!')
            # TODO: Why does this not return unity arrays, like what's
            # done below?
            return flatfield.FlatImages(None, None, None, None)

        # Slit and tilt traces are required to flat-field the data
        if not self._chk_objs(['slits', 'wavetilts']):
            # TODO: Why doesn't this fault?
            msgs.warn('Flats were requested, but there are quantities missing necessary to '
                      'create flats.  Proceeding without flat fielding....')
            return flatfield.FlatImages(None, None, None, None)

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        pixflat_image_files = self._prep_calibrations('pixelflat')

        # Return cached images
        if self._cached('flatimages', self.master_key_dict['flat']):
            self.flatimages = self.calib_dict[self.master_key_dict['flat']]['flatimages']
            return self.flatimages

        masterframe_name = masterframe.construct_file_name(flatfield.FlatImages,
                                                           self.master_key_dict['flat'], master_dir=self.master_dir)
        # The following if-elif-else does:
        #   1.  First try to load a user-supplied image
        #   2.  Try to load a MasterFrame (if reuse_masters is True)
        #   3.  Build from scratch
        # TODO: We need to document this format for the user!
        if self.par['flatfield']['frame'] != 'pixelflat':
            # - Name is explicitly correct?
            if os.path.isfile(self.par['flatfield']['frame']):
                flat_file = self.par['flatfield']['frame']
            # - Or is it in the master directory?
            elif os.path.isfile(os.path.join(self.flatField.master_dir,
                                             self.par['flatfield']['frame'])):
                flat_file = os.path.join(self.flatField.master_dir, self.par['flatfield']['frame'])
            else:
                msgs.error('Could not find user-defined flatfield file: {0}'.format(
                    self.par['flatfield']['frame']))
            # Load
            msgs.info('Using user-defined file: {0}'.format(flat_file))
            with fits.open(flat_file) as hdu:
                pixelflat = hdu[self.det].data
            # Build
            self.flatimages = flatfield.FlatImages(None, pixelflat, None, None)
        elif os.path.isfile(masterframe_name) and self.reuse_masters:
            # Load MasterFrame
            self.flatimages = flatfield.FlatImages.from_file(masterframe_name)
        else:
            # Process/combine the input pixelflat frames
            buildflatImage = flatfield.BuildFlatImage(self.spectrograph, self.det,
                                                      self.par['pixelflatframe']['process'],
                                                      pixflat_image_files, bias=self.msbias)
            pixflatimage = buildflatImage.build_image(bias=self.msbias, bpm=self.msbpm)

            # Normalize and illumination
            flatField = flatfield.FlatField(pixflatimage, self.spectrograph, self.par['flatfield'],
                det=self.det, slits=self.slits, wavetilts=self.wavetilts)
            # Run
            self.flatimages = flatField.run(show=self.show)

            # Save to Masters
            if self.save_masters:
                self.flatimages.to_master_file(self.master_dir, self.master_key_dict['flat'],  # Naming
                                           self.spectrograph.spectrograph,  # Header
                                           steps=flatField.steps)

                # If slits were tweaked by the slit illumination
                # profile, re-write them so that the tweaked slits are
                # included.
                if self.par['flatfield']['tweak_slits']:
                    # Update the SlitTraceSet master
                    self.slits.to_master_file(self.master_dir, self.master_key_dict['trace'],  # Naming
                                              self.spectrograph.spectrograph)
                    #self.slits.to_master()
                    # TODO: The waveTilts datamodel needs to be improved
                    # Objects should point to the same data
                    # TODO: Remove this line once we're sure the coding
                    # is correct so that they're not tripped.
                    # assert self.waveTilts.tilts_dict is self.flatField.tilts_dict
                    # Update the WaveTilts master
                    self.wavetilts['tilts'] = flatField.wavetilts['tilts'].copy()
                    self.wavetilts.to_master_file(self.master_dir, self.master_key_dict['tilt'],
                                                  self.spectrograph.spectrograph)

        # 4) If either of the two flats are still None, use unity
        # everywhere and print out a warning
        # TODO: These will barf if self.wavetilts['tilts'] isn't
        # defined.
        if self.flatimages.pixelflat is None:
            msgs.warn('You are not pixel flat fielding your data!!!')
        if self.flatimages.illumflat is None or not self.par['flatfield']['illumflatten']:
            msgs.warn('You are not illumination flat fielding your data!')

        # Cache & return
        self._update_cache('flat', 'flatimages', self.flatimages)
        return self.flatimages

    # TODO: why do we allow redo here?
    def get_slits(self, redo=False):
        """
        Load or generate the definition of the slit boundaries.

        Internals that must be available are :attr:`fitstbl`,
        :attr:`calib_ID`, :attr:`det`.

        Args:
            redo (bool): Redo

        Returns:
            Returns the :class:`SlitTraceSet` object (also kept
            internally as :attr:`slits`) and the slit mask array
            (numpy.ndarray; also kept internally as
            :attr:`maskslits`)

        """
        # Check for existing data
        if not self._chk_objs(['msbpm']):
            self.slits = None
            return self.slits

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        trace_image_files = self._prep_calibrations('trace')

        # Previously calculated?  If so, reuse
        if self._cached('trace', self.master_key_dict['trace']) and not redo:
            self.slits = self.calib_dict[self.master_key_dict['trace']]['trace']
            return self.slits

        # Reuse master frame?
        slit_masterframe_name = masterframe.construct_file_name(slittrace.SlitTraceSet,
                                                           self.master_key_dict['trace'],
                                                           master_dir=self.master_dir)
        if os.path.isfile(slit_masterframe_name) and self.reuse_masters:
            self.slits = slittrace.SlitTraceSet.from_file(slit_masterframe_name)
            return self.slits

        # Slits don't exist or we're not resusing them
        edge_masterframe_name = masterframe.construct_file_name(edgetrace.EdgeTraceSet,
                                                           self.master_key_dict['trace'],
                                                           master_dir=self.master_dir)
        # Reuse master frame?
        if os.path.isfile(edge_masterframe_name) and self.reuse_masters:
            self.edges = edgetrace.EdgeTraceSet.from_file(edge_masterframe_name)
        else:
            # Build me
            self.edges = edgetrace.EdgeTraceSet(self.spectrograph, self.par['slitedges'],
                                                files=trace_image_files)
            # Build the trace image
            buildtraceImage = traceimage.BuildTraceImage(self.spectrograph, self.det,
                                                         self.par['traceframe']['process'],
                                                         trace_image_files, bias=self.msbias)
            self.traceImage = buildtraceImage.build_image(bias=self.msbias, bpm=self.msbpm)

            try:
                self.edges.auto_trace(self.traceImage, bpm=self.msbpm, det=self.det,
                                      save=False) # self.save_masters) #, debug=True, show_stages=True)
            except:
                self.edges.save(edge_masterframe_name, master_dir=self.master_dir,
                                master_key=self.master_key_dict['trace'])
                msgs.error('Crashed out of finding the slits. Have saved the work done to '
                           'disk but it needs fixing.')
                return None
            else:  # This is not elegant..
                self.edges.save(edge_masterframe_name, master_dir=self.master_dir,
                                master_key=self.master_key_dict['trace'])

            # Show the result if requested
            if self.show:
                self.edges.show(thin=10, in_ginga=True)

        # Get the slits from the result of the edge tracing, delete
        # the edges object, and save the slits, if requested
        self.slits = self.edges.get_slits()
        self.edges = None
        if self.save_masters:
            self.slits.to_master_file(self.master_dir, self.master_key_dict['trace'],  # Naming
                                      self.spectrograph.spectrograph,  # Header
                                      raw_files=trace_image_files)

        # Save, initialize maskslits, and return
        self._update_cache('trace', 'trace', self.slits)
        return self.slits

    def get_wave(self):
        """
        Load or generate a wavelength image

        Requirements:
           wavetilts, slits, wv_calib, det, par, master_key

        Returns:
            ndarray: :attr:`mswave` wavelength image

        """
        # Check for existing data
        if not self._chk_objs(['wavetilts', 'slits', 'wv_calib']):
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
            self.mswave = waveimage.WaveImage(self.wavetilts['tilts'] * (self.wavetilts['tilts'].shape[0]-1.0))
            self.calib_dict[self.master_key_dict['arc']]['wave'] = self.mswave
            return self.mswave

        # Load?
        masterframe_name = masterframe.construct_file_name(
            waveimage.WaveImage, self.master_key_dict['arc'], master_dir=self.master_dir)
        if os.path.isfile(masterframe_name) and self.reuse_masters:
            self.mswave = waveimage.WaveImage.from_file(masterframe_name)
        else:  # Build
            # Instantiate
            # TODO we are regenerating this mask a lot in this module. Could reduce that
            buildwaveImage = waveimage.BuildWaveImage(self.slits, self.wavetilts['tilts'], self.wv_calib,
                                             self.spectrograph, self.det)
            self.mswave = buildwaveImage.build_wave()
            # Save to hard-drive
            if self.save_masters:
                self.mswave.to_master_file(self.master_dir, self.master_key_dict['tilts'],  # Naming
                                          self.spectrograph.spectrograph,  # Header
                                          steps=buildwaveImage.steps)

        # Cache & return
        self._update_cache('arc', 'wave', self.mswave)
        return self.mswave

    def get_wv_calib(self):
        """
        Load or generate the 1D wavelength calibrations

        Requirements:
          msarc, msbpm, slits, det, par

        Returns:
            dict, ndarray: :attr:`wv_calib` calibration dict and the updated slit mask array
        """
        # Check for existing data
        if not self._chk_objs(['msarc', 'msbpm', 'slits']):
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
            self.slits.mask |= self.wv_maskslits
            return self.wv_calib

        # No wavelength calibration requested
        if self.par['wavelengths']['reference'] == 'pixel':
            msgs.info("A wavelength calibration will not be performed")
            self.wv_calib = None
            self.wv_maskslits = np.zeros_like(self.maskslits, dtype=bool)
            self.slits.mask |= self.wv_maskslits
            return self.wv_calib

        # Grab arc binning (may be different from science!)
        #arc_rows = self.fitstbl.find_frames('arc', calib_ID=self.calib_ID, index=True)
        #self.arc_files = self.fitstbl.frame_paths(arc_rows)
        #binspec, binspat = parse.parse_binning(self.spectrograph.get_meta_value(self.arc_files[0],
        #                                                                        'binning'))
        # TODO : Do this internally when we have a wv_calib DataContainer
        binspec, binspat = parse.parse_binning(self.msarc.detector.binning)

        # Instantiate
        self.waveCalib = wavecalib.WaveCalib(self.msarc, self.slits, self.spectrograph,
                                             self.par['wavelengths'], binspectral=binspec,
                                             det=self.det,
                                             master_key=self.master_key_dict['arc'],  # For QA naming
                                             qa_path=self.qa_path, msbpm=self.msbpm)
                                             #master_dir=self.master_dir,
                                             #reuse_masters=self.reuse_masters,
        # Load from disk (MasterFrame)?
        masterframe_name = masterframe.construct_file_name(wavecalib.WaveCalib, self.master_key_dict['arc'],
                                                           master_dir=self.master_dir)
        if os.path.isfile(masterframe_name) and self.reuse_masters:
            self.wv_calib = self.waveCalib.load(masterframe_name)
        else:
            self.wv_calib, _ = self.waveCalib.run(skip_QA=(not self.write_qa))
            # Save to Masters
            if self.save_masters:
                self.waveCalib.save(outfile=masterframe_name)

        # Create the mask (needs to be done here in case wv_calib was loaded from Masters)
        # TODO: This should either be done here or save as part of the
        # master frame file.  As it is, if not loaded from the master
        # frame file, mask_maskslits is run twice, once in run above and
        # once here...
        self.wv_maskslits = self.waveCalib.make_maskslits(self.slits.nslits)
        self.slits.mask |= self.wv_maskslits

        # Save & return
        self._update_cache('arc', ('wavecalib','wvmask'), (self.wv_calib,self.wv_maskslits))
        # Return
        return self.wv_calib

    def get_tilts(self):
        """
        Load or generate the tilts image

        Requirements:
           mstilt, slits, wv_calib
           det, par, spectrograph

        Returns:
            dict, ndarray: :attr:`wavetilts` dictionary with tilts information (2D)
            and the updated slit mask array

        """
        # Check for existing data
        #TODO add mstilt_inmask to this list when it gets implemented.
        if not self._chk_objs(['mstilt', 'msbpm', 'slits', 'wv_calib']):
            msgs.error('dont have all the objects')
            self.wavetilts = None
            self.wt_maskslits = np.zeros_like(self.maskslits, dtype=bool)
            self.slits.mask |= self.wt_maskslits
            return self.wavetilts

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])
        if 'tilt' not in self.master_key_dict.keys():
            msgs.error('Tilt master key not set.  First run get_tiltimage.')

        # Return existing data
        if self._cached('wavetilts', self.master_key_dict['tilt']) \
                and self._cached('wtmask', self.master_key_dict['tilt']):
            self.wavetilts = self.calib_dict[self.master_key_dict['tilt']]['wavetilts']
            self.wt_maskslits = self.calib_dict[self.master_key_dict['tilt']]['wtmask']
            self.slits.mask |= self.wt_maskslits
            return self.wavetilts

        # Load up?
        masterframe_name = masterframe.construct_file_name(wavetilts.WaveTilts, self.master_key_dict['tilt'],
                                                           master_dir=self.master_dir)
        if os.path.isfile(masterframe_name) and self.reuse_masters:
            self.wavetilts = wavetilts.WaveTilts.from_file(masterframe_name)
            self.wt_maskslits = np.zeros(self.slits.nslits, dtype=bool)
        else: # Build
            buildwaveTilts = wavetilts.BuildWaveTilts(
                self.mstilt, self.slits, self.spectrograph, self.par['tilts'],
                self.par['wavelengths'], det=self.det, qa_path=self.qa_path,
                msbpm=self.msbpm, master_key=self.master_key_dict['tilt'])

            # TODO still need to deal with syntax for LRIS ghosts. Maybe we don't need it
            self.wavetilts, self.wt_maskslits \
                    = buildwaveTilts.run(maskslits=self.slits.mask, doqa=self.write_qa, show=self.show)
            # Save?
            if self.save_masters:
                self.wavetilts.to_master_file(self.master_dir, self.master_key_dict['tilt'],
                    self.spectrograph.spectrograph, steps=buildwaveTilts.steps)

        # Save & return
        self._update_cache('tilt', ('wavetilts','wtmask'), (self.wavetilts, self.wt_maskslits))
        self.slits.mask |= self.wt_maskslits
        return self.wavetilts

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
                # Strip ms
                iobj = obj[2:] if obj[0:2] == 'ms' else obj
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

    NOTE: Echelle uses this same class.  It had been possible there would be
    a different order of the default_steps

    ..todo:: Rename this child or eliminate altogether
    """
    def __init__(self, fitstbl, par, spectrograph, caldir=None, qadir=None, reuse_masters=False,
                 show=False, steps=None, save_masters=True):
        super(MultiSlitCalibrations, self).__init__(fitstbl, par, spectrograph, caldir=caldir,
                                                    qadir=qadir, reuse_masters=reuse_masters,
                                                    show=show, save_masters=save_masters)
        self.steps = MultiSlitCalibrations.default_steps() if steps is None else steps

    @staticmethod
    def default_steps():
        """
        This defines the steps for calibrations and their order

        Returns:
            list: Calibration steps, in order of execution

        """
        # Order matters!
        return ['bias', 'bpm', 'arc', 'tiltimg', 'slits', 'wv_calib', 'tilts', 'flats', 'wave']

    # TODO For flexure compensation add a method adjust_flexure to calibrations which will get called from extract_one
    # Notes on order of steps if flexure compensation is implemented
    #  ['bpm', 'bias', 'arc', 'tiltimg', 'slits', 'wv_calib', 'tilts', 'flats', 'wave']

