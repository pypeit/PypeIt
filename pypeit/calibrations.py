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
from pypeit import alignframe
from pypeit import flatfield
from pypeit import edgetrace
from pypeit import masterframe
from pypeit import slittrace
from pypeit import wavecalib
from pypeit import wavetilts
from pypeit.images import buildimage
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
        par (:class:`pypeit.par.pypeitpar.CalibrationsPar`):
            Parameter set defining optional parameters of PypeIt's algorithms
            for Calibrations
        spectrograph (:obj:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph object
        caldir (:obj:`str`, optional):
            Path to write the output calibrations.  If None, calibration
            data are not saved.
        qadir (:obj:`str`, optional):
            Path for quality assessment output.  If not provided, no QA
            plots are saved.
        reuse_masters (:obj:`bool`, optional):
            Load calibration files from disk if they exist
        show (:obj:`bool`, optional):
            Show plots of PypeIt's results as the code progesses.
            Requires interaction from the users.

    .. todo: Fix these

    Attributes:
        wavetilts (:class:`pypeit.wavetilts.WaveTilts`):
        mstilt (:class:`pypeit.images.buildimage.TiltImage`):
        write_qa
        show
        spectrograph
        par (:class:`pypeit.par.pypeitpar.CalibrationsPar`):
        full_par (:class:`pypeit.par.pypeitpar.PypeItPar`):
        redux_path
        master_dir
        calib_dict
        det
        frame (:obj:`int`):
            0-indexed row of the frame being calibrated in
            :attr:`fitstbl`.
        calib_ID (:obj:`int`):
            calib group ID of the current frame
        slitspat_num (:obj:`str` or :obj:`list, optional):
            Identifies a slit or slits to restrict the analysis on
            Used in :func:`get_slits` and propagated beyond

    """
    __metaclass__ = ABCMeta

    def __init__(self, fitstbl, par, spectrograph, caldir, qadir=None,
                 reuse_masters=False, show=False, slitspat_num=None):

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

        # Restrict on slits?
        self.slitspat_num = slitspat_num

        # QA
        self.qa_path = qadir
        self.write_qa = qadir is not None
        self.show = show

        # Check the directories exist
        # TODO: This should be done when the masters are saved
        if not os.path.isdir(self.master_dir):
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
        """
        Parse self.fitstbl for rows matching the calibration type
        and initialize the self.master_key_dict

        Args:
            ctype (str):
                Calibration type, e.g. 'flat', 'arc', 'bias'

        Returns:
            tuple:  2 objects
               - list:  List of image files matching input type
               - str:  master_key

        """
        # Grab rows and files
        rows = self.fitstbl.find_frames(ctype, calib_ID=self.calib_ID, index=True)
        image_files = self.fitstbl.frame_paths(rows)
        # Return
        return image_files, self.fitstbl.master_key(rows[0] if len(rows) > 0 else self.frame, det=self.det)

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
        self.msalign = None
        self.alignment = None
        self.msbias = None
        self.msbpm = None
        self.slits = None
        self.wavecalib = None
        self.wavetilts = None
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
        #  actually being done correctly as it is?
        #  We should *not* copy.
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
            master_type (str): Type of calibration frame, e.g. 'bias', 'arc', ..
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
                # Found the previous master in memory (may be None)
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
        arc_files, self.master_key_dict['arc'] = self._prep_calibrations('arc')
        masterframe_name = masterframe.construct_file_name(
            buildimage.ArcImage, self.master_key_dict['arc'], master_dir=self.master_dir)

        # Previously calculated?  If so, reuse
        if self._cached('arc', self.master_key_dict['arc']):
            self.msarc = self.calib_dict[self.master_key_dict['arc']]['arc']
            return self.msarc

        # Reuse master frame?
        if os.path.isfile(masterframe_name) and self.reuse_masters:
            self.msarc = buildimage.ArcImage.from_file(masterframe_name)
        else:  # Build it
            msgs.info("Preparing a master {0:s} frame".format(buildimage.ArcImage.master_type))
            self.msarc = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                        self.par['arcframe'], arc_files,
                                                        bias=self.msbias, bpm=self.msbpm)
            # Save
            self.msarc.to_master_file(masterframe_name)
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
        tilt_files, self.master_key_dict['tilt'] = self._prep_calibrations('tilt')
        masterframe_name = masterframe.construct_file_name(
            buildimage.TiltImage, self.master_key_dict['tilt'], master_dir=self.master_dir)

        # Previously calculated?  If so, reuse
        if self._cached('tiltimg', self.master_key_dict['tilt']):
            self.mstilt = self.calib_dict[self.master_key_dict['tilt']]['tiltimg']
            return self.mstilt

        # Reuse master frame?
        if os.path.isfile(masterframe_name) and self.reuse_masters:
            self.mstilt = buildimage.TiltImage.from_file(masterframe_name)
        else: # Build
            msgs.info("Preparing a master {0:s} frame".format(buildimage.TiltImage.master_type))
            self.mstilt = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                self.par['tiltframe'],
                                                tilt_files, bias=self.msbias, bpm=self.msbpm,
                                                         slits=self.slits)  # For flexure

            # Save to Masters
            self.mstilt.to_master_file(masterframe_name)

        # Cache
        self._update_cache('tilt', 'tiltimg', self.mstilt)
        # TODO in the future add in a tilt_inmask
        #self._update_cache('tilt', 'tilt_inmask', self.mstilt_inmask)

        # Return
        return self.mstilt

    def get_align(self):
        """
        Load or generate the alignment frame

        Requirements:
           master_key, det, par

        Returns:
            ndarray or str: :attr:`align`

        """
        # Check for existing data
        if not self._chk_objs(['msbpm', 'tslits_dict']):
            msgs.error("Don't have all the objects")

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        align_files = self._prep_calibrations('align')
        #align_rows = self.fitstbl.find_frames('align', calib_ID=self.calib_ID, index=True)
        #self.align_files = self.fitstbl.frame_paths(align_rows)
        #self.master_key_dict['align'] \
        #        = self.fitstbl.master_key(align_rows[0] if len(align_rows) > 0 else self.frame,
        #                                  det=self.det)
        masterframe_name = masterframe.construct_file_name(
            buildimage.AlignImage, self.master_key_dict['align'], master_dir=self.master_dir)

        # Previously cahded?
        if self._cached('align', self.master_key_dict['align']):
            # Previously calculated
            self.msalign = self.calib_dict[self.master_key_dict['align']]['align']
        elif os.path.isfile(masterframe_name) and self.reuse_masters:
            self.msalign = buildimage.AlignImage.from_file(masterframe_name)
        else:
            # Instantiate with everything needed to generate the image (in case we do)
            #self.alignFrame = alignframe.AlignFrame(self.spectrograph, files=self.align_files,
            #                                  det=self.det, msbias=self.msbias,
            #                                  par=self.par['alignframe'],
            #                                  master_key=self.master_key_dict['align'],
            #                                  master_dir=self.master_dir,
            #                                  reuse_masters=self.reuse_masters)
            self.align = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                         self.par['alignframe'],
                                                         align_files, bias=self.msbias, bpm=self.msbpm)

            # Load the MasterFrame (if it exists and is desired)?
            #self.msalign = self.alignFrame.load()
            #if self.msalign is None:  # Otherwise build it
            #    msgs.info("Preparing a master {0:s} frame".format(self.alignFrame.master_type))
            #    self.msalign = self.alignFrame.build_image(bias=self.msbias, bpm=self.msbpm)
            #    # Need to set head0 here, since a master align frame loaded from file will have head0 set.
            #    self.msalign.head0 = self.alignFrame.build_master_header(steps=self.alignFrame.process_steps,
            #                                                         raw_files=self.alignFrame.file_list)
            #   # Save to Masters
            # Save to Masters
            self.msalign.to_master_file(self.master_dir, self.master_key_dict['align'],  # Naming
                                       self.spectrograph.spectrograph,  # Header
                                       steps=self.msalign.process_steps,
                                       raw_files=align_files)

            # Store the alignment frame
            self._update_cache('align', 'align', self.msalign)

        # Check if the alignment dictionary exists
        if self._cached('align_dict', self.master_key_dict['align']) \
                and self._cached('wtmask', self.master_key_dict['align']):
            self.align_dict = self.calib_dict[self.master_key_dict['align']]['align_dict']
        else:
            # Extract some header info needed by the algorithm
            binning = self.spectrograph.get_meta_value(self.align_files[0], 'binning')

            # Instantiate
            self.alignment = alignframe.Alignment(self.msalign, self.tslits_dict, self.spectrograph,
                                                  self.par['alignment'],
                                                  det=self.det, binning=binning,
                                                  master_key=self.master_key_dict['align'],
                                                  master_dir=self.master_dir,
                                                  reuse_masters=self.reuse_masters,
                                                  qa_path=self.qa_path, msbpm=self.msbpm)

            # Master
            self.align_dict = self.alignment.load()
            if self.align_dict is None:
                self.align_dict = self.alignment.run(self.show)
                self.alignment.save()

            # Save & return
            self._update_cache('align', 'align_dict', self.align_dict)

        return self.msalign, self.align_dict

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
        bias_files, self.master_key_dict['bias'] = self._prep_calibrations('bias')
        # Construct the name, in case we need it
        masterframe_name = masterframe.construct_file_name(buildimage.BiasImage,
                                                           self.master_key_dict['bias'],
                                                           master_dir=self.master_dir)

        # This needs to come after prep or the code crashes when saving as master_key_dict['bias'] is not set
        if self.par['biasframe']['useframe'].lower() == 'none':
            self.msbias = None
            return self.msbias


        # Grab from internal dict (or hard-drive)?
        if self._cached('bias', self.master_key_dict['bias']):
            self.msbias = self.calib_dict[self.master_key_dict['bias']]['bias']
            msgs.info("Reloading the bias from the internal dict")
            return self.msbias

        # Try to load?
        if os.path.isfile(masterframe_name) and self.reuse_masters:
            self.msbias = buildimage.BiasImage.from_file(masterframe_name)
        else:
            # Without files, we are stuck
            if len(bias_files) == 0:
                self.msbias = None
            else:
                # Build it
                self.msbias = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                    self.par['biasframe'], bias_files)
                # Save it?
                self.msbias.to_master_file(masterframe_name)
                        #self.master_dir, self.master_key_dict['bias'],  # Naming
                        #                  self.spectrograph.spectrograph,  # Header
                        #                  steps=self.msbias.process_steps,
                        #                  raw_files=bias_files)
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
            :class:`pypeit.flatfield.FlatImages`:
        """
        # Check for existing data
        if not self._chk_objs(['msarc', 'msbpm', 'slits', 'wv_calib']):
            msgs.error('Must have the arc, bpm, slits, and wv_calib defined to proceed!')

        if self.par['flatfield']['method'] == 'skip':
            # User does not want to flat-field
            msgs.warn('Parameter calibrations.flatfield.method is set to skip. You are NOT '
                      'flatfielding your data!!!')
            # TODO: Why does this not return unity arrays, like what's
            # done below?
            return flatfield.FlatImages(None)

        # Slit and tilt traces are required to flat-field the data
        if not self._chk_objs(['slits', 'wavetilts']):
            # TODO: Why doesn't this fault?
            msgs.warn('Flats were requested, but there are quantities missing necessary to '
                      'create flats.  Proceeding without flat fielding....')
            return flatfield.FlatImages(None)

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        trace_image_files, self.master_key_dict['flat'] = self._prep_calibrations('trace')

        # Return cached images
        if self._cached('flatimages', self.master_key_dict['flat']):
            self.flatimages = self.calib_dict[self.master_key_dict['flat']]['flatimages']
            self.flatimages.is_synced(self.slits)
            self.slits.mask_flats(self.flatimages)
            return self.flatimages

        masterframe_filename = masterframe.construct_file_name(flatfield.FlatImages,
                                                           self.master_key_dict['flat'], master_dir=self.master_dir)
        # The following if-elif-else does:
        #   1.  Try to load a MasterFrame (if reuse_masters is True)
        #   2.  Build from scratch
        #   3.  Replace the built pixel-flat with user supplied (e.g. LRISb)

        if os.path.isfile(masterframe_filename) and self.reuse_masters:
            # Load MasterFrame
            self.flatimages = flatfield.FlatImages.from_file(masterframe_filename)
            self.flatimages.is_synced(self.slits)
            self.slits.mask_flats(self.flatimages)
        elif len(trace_image_files) > 0:
            # Process/combine the input pixelflat frames
            # TODO -- Include an illum frametype eventually
            # TODO This is incorrect nad is only a hack for LRIS-Blue where we want to only run this code
            # for the illumination flat construction. There needs to be a pixelflat and an illumflat, and
            # in cases where the files are different, we need to run the flat field code twice.
            stacked_traceflat = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                    self.par['traceframe'],
                                                    trace_image_files,
                                                    bias=self.msbias, bpm=self.msbpm)
            # Normalize and illumination
            flatField = flatfield.FlatField(stacked_traceflat, self.spectrograph, self.par['flatfield'],
                self.slits, self.wavetilts)
            # Run
            self.flatimages = flatField.run(show=self.show)

            # Save to Masters
            self.flatimages.to_master_file(masterframe_filename)
            # Save slits too, in case they were tweaked
            self.slits.to_master_file()
        else:
            self.flatimages = flatfield.FlatImages(None, None, None, None)
            msgs.warn("No pixelflats provided")

        # TODO: We need to document this format for the user!
        # User supplied pixelf flat??
        if self.par['flatfield']['frame'] != 'pixelflat':
            # - Name is explicitly correct?
            if os.path.isfile(self.par['flatfield']['frame']):
                flat_file = self.par['flatfield']['frame']
            else:
                msgs.error('Could not find user-defined flatfield file: {0}'.format(
                    self.par['flatfield']['frame']))
            # Load
            msgs.info('Using user-defined file: {0}'.format(flat_file))
            with fits.open(flat_file) as hdu:
                user_pixelflat = hdu[self.det].data
            self.flatimages.pixelflat = user_pixelflat

        # 4) If flat is still None, print out a warning
        if self.flatimages.pixelflat is None:
            msgs.warn('You are not pixel flat fielding your data!!!')

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
        trace_image_files, self.master_key_dict['trace'] = self._prep_calibrations('trace')

        # Previously calculated?  If so, reuse
        if self._cached('trace', self.master_key_dict['trace']) and not redo:
            self.slits = self.calib_dict[self.master_key_dict['trace']]['trace']
            # Reset the bitmask
            self.slits.mask = self.slits.mask_init.copy()
            return self.slits

        # Reuse master frame?
        slit_masterframe_name = masterframe.construct_file_name(slittrace.SlitTraceSet,
                                                           self.master_key_dict['trace'],
                                                           master_dir=self.master_dir)
        if os.path.isfile(slit_masterframe_name) and self.reuse_masters:
            self.slits = slittrace.SlitTraceSet.from_file(slit_masterframe_name)
            # Reset the bitmask
            self.slits.mask = self.slits.mask_init.copy()
        else:
            # Slits don't exist or we're not resusing them
            edge_masterframe_name = masterframe.construct_file_name(edgetrace.EdgeTraceSet,
                                                               self.master_key_dict['trace'],
                                                               master_dir=self.master_dir)
            # Reuse master frame?
            if os.path.isfile(edge_masterframe_name) and self.reuse_masters:
                self.edges = edgetrace.EdgeTraceSet.from_file(edge_masterframe_name)
            else:
                # Build the trace image
                self.traceImage = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                        self.par['traceframe'], trace_image_files,
                                                        bias=self.msbias, bpm=self.msbpm)
                # Build me
                self.edges = edgetrace.EdgeTraceSet(self.traceImage, self.spectrograph, self.par['slitedges'],
                                                    files=trace_image_files)

                try:
                    self.edges.auto_trace(bpm=self.msbpm, det=self.det, save=False)
                except:
                    self.edges.save(edge_masterframe_name, master_dir=self.master_dir,
                                    master_key=self.master_key_dict['trace'])
                    msgs.error('Crashed out of finding the slits. Have saved the work done to '
                               'disk but it needs fixing.')
                    return None
                else:
                    self.edges.save(edge_masterframe_name, master_dir=self.master_dir,
                                    master_key=self.master_key_dict['trace'])

                # Show the result if requested
                if self.show:
                    self.edges.show(thin=10, in_ginga=True)

            # Get the slits from the result of the edge tracing, delete
            # the edges object, and save the slits, if requested
            self.slits = self.edges.get_slits()
            self.edges = None
            self.slits.to_master_file(slit_masterframe_name)

        # User mask?
        if self.slitspat_num is not None:
            self.slits.user_mask(self.det, self.slitspat_num)

        # Save, initialize maskslits, and return
        self._update_cache('trace', 'trace', self.slits)
        return self.slits

#    def get_wave(self):
#        """
#        Load or generate a wavelength image
#
#        Requirements:
#           wavetilts, slits, wv_calib, det, par, master_key
#
#        Returns:
#            `numpy.ndarray`_: :attr:`mswave` wavelength image
#
#        """
#        msgs.error("NO LONGER USED.  GENERATE ON-THE-SPOT with code in pypeit.wavecalib")
#        # Check for existing data
#        if not self._chk_objs(['wavetilts', 'slits', 'wv_calib']):
#            self.mswave = None
#            return self.mswave
#
#        # Check internals
#        self._chk_set(['det', 'par'])
#
#        # Return existing data
#        if self._cached('wave', self.master_key_dict['arc']):
#            self.mswave = self.calib_dict[self.master_key_dict['arc']]['wave']
#            return self.mswave
#
#        # No wavelength calibration requested
#        if self.par['wavelengths']['reference'] == 'pixel':
#            msgs.warn('No wavelength calibration performed!')
#            self.mswave = waveimage.WaveImage(self.wavetilts['tilts'] * (self.wavetilts['tilts'].shape[0]-1.0))
#            self.calib_dict[self.master_key_dict['arc']]['wave'] = self.mswave
#            return self.mswave
#
#        # Load?
#        masterframe_name = masterframe.construct_file_name(
#            waveimage.WaveImage, self.master_key_dict['arc'], master_dir=self.master_dir)
#        if os.path.isfile(masterframe_name) and self.reuse_masters:
#            self.mswave = waveimage.WaveImage.from_file(masterframe_name)
#        else:  # Build
#            # Instantiate
#            # TODO we are regenerating this mask a lot in this module. Could reduce that
#            buildwaveImage = waveimage.BuildWaveImage(self.slits, self.wavetilts['tilts'], self.wv_calib,
#                                             self.spectrograph, self.det)
#            self.mswave = buildwaveImage.build_wave()
#            # Save to hard-drive
#            self.mswave.to_master_file(masterframe_name)
#                    #self.master_dir, self.master_key_dict['arc'],  # Naming
#                    #                      self.spectrograph.spectrograph,  # Header
#                    #                      steps=buildwaveImage.steps)
#
#        # Cache & return
#        self._update_cache('arc', 'wave', self.mswave)
#        return self.mswave

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
            msgs.error('Not enough information to load/generate the wavelength calibration')

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])
        if 'arc' not in self.master_key_dict.keys():
            msgs.error('Arc master key not set.  First run get_arc.')

        # Return existing data
        if self._cached('wavecalib', self.master_key_dict['arc']):
            self.wv_calib = self.calib_dict[self.master_key_dict['arc']]['wavecalib']
            self.slits.mask_wvcalib(self.wv_calib)
            return self.wv_calib

        # No wavelength calibration requested
        if self.par['wavelengths']['reference'] == 'pixel':
            msgs.info("A wavelength calibration will not be performed")
            self.wv_calib = None
            return self.wv_calib

        # Grab arc binning (may be different from science!)
        # TODO : Do this internally when we have a wv_calib DataContainer
        binspec, binspat = parse.parse_binning(self.msarc.detector.binning)

        # Instantiate
        self.waveCalib = wavecalib.WaveCalib(self.msarc, self.slits, self.spectrograph,
                                             self.par['wavelengths'], binspectral=binspec,
                                             det=self.det,
                                             master_key=self.master_key_dict['arc'],  # For QA naming
                                             qa_path=self.qa_path, msbpm=self.msbpm)
        # Load from disk (MasterFrame)?
        masterframe_name = masterframe.construct_file_name(wavecalib.WaveCalib, self.master_key_dict['arc'],
                                                           master_dir=self.master_dir)
        if os.path.isfile(masterframe_name) and self.reuse_masters:
            # Load from disk
            self.wv_calib = self.waveCalib.load(masterframe_name)
            self.slits.mask_wvcalib(self.wv_calib)
        else:
            self.wv_calib = self.waveCalib.run(skip_QA=(not self.write_qa))
            # Save to Masters
            self.waveCalib.save(outfile=masterframe_name)

        # Save & return
        self._update_cache('arc', 'wavecalib', self.wv_calib)
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

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])
        if 'tilt' not in self.master_key_dict.keys():
            msgs.error('Tilt master key not set.  First run get_tiltimage.')

        # Return existing data
        if self._cached('wavetilts', self.master_key_dict['tilt']):
            self.wavetilts = self.calib_dict[self.master_key_dict['tilt']]['wavetilts']
            self.wavetilts.is_synced(self.slits)
            self.slits.mask_wavetilts(self.wavetilts)
            return self.wavetilts

        # Load up?
        masterframe_name = masterframe.construct_file_name(wavetilts.WaveTilts, self.master_key_dict['tilt'],
                                                           master_dir=self.master_dir)
        if os.path.isfile(masterframe_name) and self.reuse_masters:
            self.wavetilts = wavetilts.WaveTilts.from_file(masterframe_name)
            self.wavetilts.is_synced(self.slits)
            self.slits.mask_wavetilts(self.wavetilts)
        else: # Build
            # Flexure
            _spat_flexure = self.mstilt.spat_flexure \
                if self.par['tiltframe']['process']['spat_flexure_correct'] else None
            # Instantiate
            buildwaveTilts = wavetilts.BuildWaveTilts(
                self.mstilt, self.slits, self.spectrograph, self.par['tilts'],
                self.par['wavelengths'], det=self.det, qa_path=self.qa_path,
                master_key=self.master_key_dict['tilt'], spat_flexure=_spat_flexure)

            # TODO still need to deal with syntax for LRIS ghosts. Maybe we don't need it
            self.wavetilts = buildwaveTilts.run(doqa=self.write_qa, show=self.show)
            # Save?
            self.wavetilts.to_master_file(masterframe_name)

        # Save & return
        self._update_cache('tilt', 'wavetilts', self.wavetilts)
        return self.wavetilts

    def run_the_steps(self):
        """
        Run full the full recipe of calibration steps

        """
        for step in self.steps:
            getattr(self, 'get_{:s}'.format(step))()
        msgs.info("Calibration complete!")
        msgs.info("#######################################################################")

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
    def __init__(self, fitstbl, par, spectrograph, caldir, **kwargs):
        super(MultiSlitCalibrations, self).__init__(fitstbl, par, spectrograph, caldir, **kwargs)
        self.steps = MultiSlitCalibrations.default_steps()

    @staticmethod
    def default_steps():
        """
        This defines the steps for calibrations and their order

        Returns:
            list: Calibration steps, in order of execution

        """
        # Order matters!
        return ['bias', 'bpm', 'slits', 'arc', 'tiltimg', 'wv_calib', 'tilts', 'flats']


class IFUCalibrations(Calibrations):
    """
    Child of Calibrations class for performing IFU calibrations.
    See :class:`pypeit.calibrations.Calibrations` for arguments.

    """

    def __init__(self, fitstbl, par, spectrograph, caldir, **kwargs):
        super(IFUCalibrations, self).__init__(fitstbl, par, spectrograph, caldir, **kwargs)
        self.steps = IFUCalibrations.default_steps()

    @staticmethod
    def default_steps():
        """
        This defines the steps for calibrations and their order

        Returns:
            list: Calibration steps, in order of execution

        """
        # Order matters!
        return ['bias', 'bpm', 'arc', 'tiltimg', 'slits', 'wv_calib', 'tilts', 'align', 'flats', 'wave']

