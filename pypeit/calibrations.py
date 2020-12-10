"""
Class for guiding calibration object generation in PypeIt

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os

from abc import ABCMeta
from collections import Counter

from IPython import embed

import numpy as np

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
from pypeit import io
from pypeit import utils


class Calibrations(object):
    """
    This class is primarily designed to guide the generation of
    calibration images and objects in PypeIt.

    To avoid rebuilding MasterFrames that were generated during this execution
    of PypeIt, the class performs book-keeping of these master frames

    Args:
        fitstbl (:class:`pypeit.metadata.PypeItMetaData`, None):
            The class holding the metadata for all the frames in this PypeIt run.
            If None, we are using this class as a glorified dict to hold the objects.
        par (:class:`pypeit.par.pypeitpar.CalibrationsPar`):
            Parameter set defining optional parameters of PypeIt's algorithms
            for Calibrations
        spectrograph (:obj:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph object
        caldir (:obj:`str`, None):
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
        fitstbl (:class:`pypeit.metadata.PypeItMetaData`):
            Table with metadata for all fits files to reduce.
        wavetilts (:class:`pypeit.wavetilts.WaveTilts`):
        mstilt (:class:`pypeit.images.buildimage.TiltImage`):
        flatimages (:class:`pypeit.flatfield.FlatImages`):
        msbias (:class:`pypeit.images.buildimage.BiasImage`):
        msdark (:class:`pypeit.images.buildimage.DarkImage`):
        msbpm (`numpy.ndarray`_):
        msarc (:class:`pypeit.images.buildimage.ArcImage`):
            Master arc-lamp image.
        alignments (:class:`pypeit.alignframe.Alignments`):
        wv_calib (:class:`pypeit.wavecalib.WaveCalib`):
        slits (:class:`pypeit.slittrace.SlitTraceSet`):

        write_qa
        show
        spectrograph
        par (:class:`pypeit.par.pypeitpar.CalibrationsPar`):
        full_par (:class:`pypeit.par.pypeitpar.PypeItPar`):
        redux_path
        master_dir
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

    @classmethod
    def get_instance(cls, fitstbl, par, spectrograph, caldir, qadir=None,
                     reuse_masters=False, show=False, slitspat_num=None):
        """
        """
        pypeline = spectrograph.pypeline
        if spectrograph.pypeline == 'Echelle':
            pypeline = 'MultiSlit'
        return next(c for c in cls.__subclasses__()
                    if c.__name__ == (pypeline + 'Calibrations'))(
            fitstbl, par, spectrograph, caldir, qadir=qadir,
                     reuse_masters=reuse_masters, show=show, slitspat_num=slitspat_num)

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
        if caldir is not None and not os.path.isdir(self.master_dir):
            os.makedirs(self.master_dir)
        # TODO: This should be done when the qa plots are saved
        if self.write_qa and not os.path.isdir(os.path.join(self.qa_path, 'PNGs')):
            os.makedirs(os.path.join(self.qa_path, 'PNGs'))

        # Attributes
        self.det = None
        self.frame = None
        self.binning = None

        self.shape = None

        self.msarc = None
        self.mstilt = None
        self.alignments = None
        self.msbias = None
        self.msdark = None
        self.msbpm = None
        self.wv_calib = None
        self.slits = None

        self.wavecalib = None
        self.wavetilts = None
        self.flatimages = None
        self.calib_ID = None
        self.master_key_dict = {}

        # Steps
        self.steps = []

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

    def set_config(self, frame, det, par=None):
        """
        Specify the parameters of the Calibrations class and reset all
        the internals to None. The internal dict is left unmodified.

        Args:
            frame (int): Frame index in the fitstbl
            det (int): Detector number
            par (:class:`pypeit.par.pypeitpar.CalibrationsPar`):

        """
        # Reset internals to None
        # NOTE: This sets empties calib_ID and master_key_dict so must
        # be done here first before these things are initialized below.

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
            `numpy.ndarray`_: :attr:`msarc` image

        """
        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        arc_files, self.master_key_dict['arc'] = self._prep_calibrations('arc')
        masterframe_name = masterframe.construct_file_name(
            buildimage.ArcImage, self.master_key_dict['arc'], master_dir=self.master_dir)

        # Reuse master frame?
        if os.path.isfile(masterframe_name) and self.reuse_masters:
            self.msarc = buildimage.ArcImage.from_file(masterframe_name)
        elif len(arc_files) == 0:
            msgs.warn("No frametype=arc files to build arc")
            return
        else:  # Build it
            msgs.info("Preparing a master {0:s} frame".format(buildimage.ArcImage.master_type))
            self.msarc = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                        self.par['arcframe'], arc_files,
                                                        bias=self.msbias, bpm=self.msbpm)
            # Save
            self.msarc.to_master_file(masterframe_name)

        # Return
        return self.msarc

    def get_tiltimg(self):
        """
        Load or generate the Tilt image

        Requirements:
          master_key, det, par

        Args:

        Returns:
            `numpy.ndarray`_: :attr:`mstilt` image

        """
        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        tilt_files, self.master_key_dict['tilt'] = self._prep_calibrations('tilt')
        masterframe_name = masterframe.construct_file_name(
            buildimage.TiltImage, self.master_key_dict['tilt'], master_dir=self.master_dir)

        # Reuse master frame?
        if os.path.isfile(masterframe_name) and self.reuse_masters:
            self.mstilt = buildimage.TiltImage.from_file(masterframe_name)
        elif len(tilt_files) == 0:
            msgs.warn("No frametype=tilt files to build tiltimg")
            return
        else: # Build
            msgs.info("Preparing a master {0:s} frame".format(buildimage.TiltImage.master_type))
            self.mstilt = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                self.par['tiltframe'],
                                                tilt_files, bias=self.msbias, bpm=self.msbpm,
                                                         slits=self.slits)  # For flexure

            # Save to Masters
            self.mstilt.to_master_file(masterframe_name)

        # TODO in the future add in a tilt_inmask
        #self._update_cache('tilt', 'tilt_inmask', self.mstilt_inmask)

        # Return
        return self.mstilt

    def get_align(self):
        """
        Load or generate the alignment frame

        Requires: :attr:`slits`, :attr:`det`, :attr:`par`

        Returns:
            :class:`pypeit.alignframe.Alignments`:

        """
        # Check for existing data
        if not self._chk_objs(['msbpm', 'slits']):
            msgs.error('Must have the bpm and slits to make the alignments!')

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        align_files, self.master_key_dict['align'] = self._prep_calibrations('align')

        masterframe_filename = masterframe.construct_file_name(alignframe.Alignments,
                                                               self.master_key_dict['align'],
                                                               master_dir=self.master_dir)

        # Reuse master frame?
        if os.path.isfile(masterframe_filename) and self.reuse_masters:
            self.alignments = alignframe.Alignments.from_file(masterframe_filename)
            self.alignments.is_synced(self.slits)
            return self.alignments

        msalign = buildimage.buildimage_fromlist(self.spectrograph, self.det, self.par['alignframe'], align_files,
                                                 bias=self.msbias, bpm=self.msbpm)

        # Extract some header info needed by the algorithm
        binning = self.spectrograph.get_meta_value(align_files[0], 'binning')

        # Instantiate
        alignment = alignframe.TraceAlignment(msalign, self.slits, self.spectrograph, self.par['alignment'],
                                              det=self.det, binning=binning, qa_path=self.qa_path, msbpm=self.msbpm)
        # Run
        self.alignments = alignment.run(show=self.show)
        # Save to Masters
        self.alignments.to_master_file(masterframe_filename)

        return self.alignments

    def get_bias(self):
        """
        Load or generate the bias frame/command

        Requirements:
           master_key, det, par

        Returns:
            :class:`pypeit.images.buildimage.BiasImage`:

        """

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        bias_files, self.master_key_dict['bias'] = self._prep_calibrations('bias')
        # Construct the name, in case we need it
        masterframe_name = masterframe.construct_file_name(buildimage.BiasImage,
                                                           self.master_key_dict['bias'],
                                                           master_dir=self.master_dir)

        if self.par['biasframe']['useframe'] is not None:
            msgs.error("Not ready to load from disk")

        # Try to load?
        if os.path.isfile(masterframe_name) and self.reuse_masters:
            self.msbias = buildimage.BiasImage.from_file(masterframe_name)
        elif len(bias_files) == 0:
            self.msbias = None
        else:
            # Build it
            self.msbias = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                         self.par['biasframe'], bias_files)
            # Save it?
            self.msbias.to_master_file(masterframe_name)

        # Return
        return self.msbias

    def get_dark(self):
        """
        Load or generate the dark image

        Requirements:
           master_key, det, par

        Returns:
            :class:`pypeit.images.buildimage.DarkImage`:

        """

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        dark_files, self.master_key_dict['dark'] = self._prep_calibrations('dark')
        # Construct the name, in case we need it
        masterframe_name = masterframe.construct_file_name(buildimage.DarkImage,
                                                           self.master_key_dict['dark'],
                                                           master_dir=self.master_dir)

        # Try to load?
        if os.path.isfile(masterframe_name) and self.reuse_masters:
            self.msdark = buildimage.DarkImage.from_file(masterframe_name)
        elif len(dark_files) == 0:
            self.msdark = None
        else:
            # TODO: Should this include the bias?
            # Build it
            self.msdark = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                    self.par['darkframe'], dark_files)
            # Save it?
            self.msdark.to_master_file(masterframe_name)

        # Return
        return self.msdark


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
            `numpy.ndarray`_: :attr:`msbpm` image of bad pixel mask

        """
        # Check internals
        self._chk_set(['par', 'det'])

        # Generate a bad pixel mask (should not repeat)
        self.master_key_dict['bpm'] = self.fitstbl.master_key(self.frame, det=self.det)

        # Build the data-section image
        sci_image_file = self.fitstbl.frame_paths(self.frame)

        # Check if a bias frame exists, and if a BPM should be generated
        msbias = None
        if self.par['bpm_usebias']:
            msbias = self.msbias
        # Build it
        self.msbpm = self.spectrograph.bpm(sci_image_file, self.det, msbias=msbias)
        self.shape = self.msbpm.shape

        # Return
        return self.msbpm

    def get_flats(self):
        """
        Load or generate a normalized pixel flat and slit illumination
        flat.

        Requires :attr:`slits`, :attr:`wavetilts`, :attr:`det`,
        :attr:`par`.

        Constructs :attr:`flatimages`.

        """
        # Check for existing data
        if not self._chk_objs(['msarc', 'msbpm', 'slits', 'wv_calib']):
            msgs.warn('Must have the arc, bpm, slits, and wv_calib defined to make flats!  Skipping and may crash down the line')
            self.flatimages = flatfield.FlatImages()
            return

        # Slit and tilt traces are required to flat-field the data
        if not self._chk_objs(['slits', 'wavetilts']):
            # TODO: Why doesn't this fault?
            msgs.warn('Flats were requested, but there are quantities missing necessary to '
                      'create flats.  Proceeding without flat fielding....')
            self.flatimages = flatfield.FlatImages()
            return

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        illum_image_files, self.master_key_dict['flat'] = self._prep_calibrations('illumflat')
        pixflat_image_files, self.master_key_dict['flat'] = self._prep_calibrations('pixelflat')

        masterframe_filename = masterframe.construct_file_name(flatfield.FlatImages,
                                                           self.master_key_dict['flat'], master_dir=self.master_dir)
        # The following if-elif-else does:
        #   1.  Try to load a MasterFrame (if reuse_masters is True).  If successful, pass it back
        #   2.  Build from scratch
        #   3.  Load any user-supplied images to over-ride any built

        # Load MasterFrame?
        if os.path.isfile(masterframe_filename) and self.reuse_masters:
            flatimages = flatfield.FlatImages.from_file(masterframe_filename)
            flatimages.is_synced(self.slits)
            # Load user defined files
            if self.par['flatfield']['pixelflat_file'] is not None:
                # Load
                msgs.info('Using user-defined file: {0}'.format('pixelflat_file'))
                with io.fits_open(self.par['flatfield']['pixelflat_file']) as hdu:
                    nrm_image = flatfield.FlatImages(pixelflat_norm=hdu[self.det].data)
                    flatimages = flatfield.merge(flatimages, nrm_image)
            self.flatimages = flatimages
            # update slits
            self.slits.mask_flats(self.flatimages)
            return self.flatimages

        # Generate the image
        pixelflatImages, illumflatImages = None, None
        # Check if the image files are the same
        pix_is_illum = Counter(illum_image_files) == Counter(pixflat_image_files)
        if len(pixflat_image_files) > 0:
            pixel_flat = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                        self.par['pixelflatframe'],
                                                        pixflat_image_files, dark=self.msdark,
                                                        bias=self.msbias, bpm=self.msbpm)
            # Initialise the pixel flat
            pixelFlatField = flatfield.FlatField(pixel_flat, self.spectrograph,
                                                 self.par['flatfield'], self.slits, self.wavetilts, self.wv_calib)
            # Generate
            pixelflatImages = pixelFlatField.run(show=self.show)

        # Only build illum_flat if the input files are different from the pixel flat
        if (not pix_is_illum) and len(illum_image_files) > 0:
            illum_flat = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                        self.par['illumflatframe'],
                                                        illum_image_files, dark=self.msdark,
                                                        bias=self.msbias, bpm=self.msbpm)
            # Initialise the pixel flat
            illumFlatField = flatfield.FlatField(illum_flat, self.spectrograph,
                                                 self.par['flatfield'], self.slits, self.wavetilts,
                                                 self.wv_calib, spat_illum_only=True)
            # Generate
            illumflatImages = illumFlatField.run(show=self.show)

        # Merge the illum flat with the pixel flat
        if pixelflatImages is not None:
            # Combine the pixelflat and illumflat parameters into flatimages.
            # This will merge the attributes of pixelflatImages that are not None
            # with the attributes of illflatImages that are not None. Default is
            # to take pixelflatImages.
            flatimages = flatfield.merge(pixelflatImages, illumflatImages)
        else:
            # No pixel flat, but there might be an illumflat. This will mean that
            # the attributes prefixed with 'pixelflat_' will all be None.
            flatimages = illumflatImages

        # Save flat images
        if flatimages is not None:
            flatimages.to_master_file(masterframe_filename)
            # Save slits too, in case they were tweaked
            self.slits.to_master_file()

        # 3) Load user-supplied images
        #  NOTE:  This is the *final* images, not just a stack
        #  And it will over-ride what is generated below (if generated)
        if self.par['flatfield']['pixelflat_file'] is not None:
            # Load
            msgs.info('Using user-defined file: {0}'.format('pixelflat_file'))
            with io.fits_open(self.par['flatfield']['pixelflat_file']) as hdu:
                flatimages = flatfield.merge(flatimages, flatfield.FlatImages(pixelflat_norm=hdu[self.det].data))

        self.flatimages = flatimages
        # Return
        return self.flatimages

    def get_slits(self):
        """
        Load or generate the definition of the slit boundaries.

        Internals that must be available are :attr:`fitstbl`,
        :attr:`calib_ID`, :attr:`det`.

        Returns:
            :class:`pypeit.slittrace.SlitTraceSet`: Traces of the
            slit edges; also kept internally as :attr:`slits`.

        """
        # Check for existing data
        if not self._chk_objs(['msbpm']):
            return

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        trace_image_files, self.master_key_dict['trace'] = self._prep_calibrations('trace')

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
            elif len(trace_image_files) == 0:
                msgs.warn("No frametype=trace files to build slits")
                return None
            else:
                # Build the trace image
                self.traceImage = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                        self.par['traceframe'], trace_image_files,
                                                        bias=self.msbias, bpm=self.msbpm,
                                                        dark=self.msdark)
                self.edges = edgetrace.EdgeTraceSet(self.traceImage, self.spectrograph,
                                                    self.par['slitedges'], bpm=self.msbpm,
                                                    auto=True)
                self.edges.to_master_file(edge_masterframe_name)

                # Show the result if requested
                if self.show:
                    self.edges.show(in_ginga=True)

            # Get the slits from the result of the edge tracing, delete
            # the edges object, and save the slits, if requested
            self.slits = self.edges.get_slits()
            self.edges = None
            self.slits.to_master_file(slit_masterframe_name)

        # User mask?
        if self.slitspat_num is not None:
            self.slits.user_mask(self.det, self.slitspat_num)

        return self.slits

    def get_wv_calib(self):
        """
        Load or generate the 1D wavelength calibrations

        Requirements:
          msarc, msbpm, slits, det, par

        Returns:
            dict: :attr:`wv_calib` calibration dict and the updated slit mask array
        """
        # Check for existing data
        if not self._chk_objs(['msarc', 'msbpm', 'slits']):
            msgs.warn('Not enough information to load/generate the wavelength calibration. Skipping and may crash down the line')
            return None

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])
        if 'arc' not in self.master_key_dict.keys():
            msgs.error('Arc master key not set.  First run get_arc.')

        # No wavelength calibration requested
        if self.par['wavelengths']['reference'] == 'pixel':
            msgs.info("A wavelength calibration will not be performed")
            self.wv_calib = None
            return self.wv_calib

        # Grab arc binning (may be different from science!)
        # TODO : Do this internally when we have a wv_calib DataContainer
        binspec, binspat = parse.parse_binning(self.msarc.detector.binning)

        masterframe_name = masterframe.construct_file_name(wavecalib.WaveCalib,
                                                           self.master_key_dict['arc'],
                                                           master_dir=self.master_dir)
        if os.path.isfile(masterframe_name) and self.reuse_masters:
            self.wv_calib = wavecalib.WaveCalib.from_file(masterframe_name)
            self.wv_calib.chk_synced(self.slits)
            self.slits.mask_wvcalib(self.wv_calib)
        else:
            # Instantiate
            self.waveCalib = wavecalib.BuildWaveCalib(self.msarc, self.slits, self.spectrograph,
                                             self.par['wavelengths'], binspectral=binspec,
                                             det=self.det,
                                             master_key=self.master_key_dict['arc'],  # For QA naming
                                             qa_path=self.qa_path, msbpm=self.msbpm)
            self.wv_calib = self.waveCalib.run(skip_QA=(not self.write_qa))
            # Save to Masters
            self.wv_calib.to_master_file(masterframe_name)

        # Return
        return self.wv_calib

    def get_tilts(self):
        """
        Load or generate the tilts image

        Requirements:
           mstilt, slits, wv_calib
           det, par, spectrograph

        Returns:
            :class:`pypeit.wavetilts.WaveTilts`:

        """
        # Check for existing data
        #TODO add mstilt_inmask to this list when it gets implemented.
        if not self._chk_objs(['mstilt', 'msbpm', 'slits', 'wv_calib']):
            msgs.warn('dont have all the objects for tilts.  Skipping and may crash down the line..')
            return None

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])
        if 'tilt' not in self.master_key_dict.keys():
            msgs.error('Tilt master key not set.  First run get_tiltimage.')

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
        return ['bias', 'dark', 'bpm', 'slits', 'arc', 'tiltimg', 'wv_calib', 'tilts', 'flats']


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
        return ['bias', 'dark', 'bpm', 'arc', 'tiltimg', 'slits', 'wv_calib', 'tilts', 'align', 'flats']


def check_for_calibs(par, fitstbl, raise_error=True, cut_cfg=None):
    """
    Perform a somewhat quick and dirty check to see if the user
    has provided all of the calibration frametype's to reduce
    the science frames

    Args:
        par (:class:`pypeit.par.pypeitpar.PypeItPar`):
        fitstbl (:class:`pypeit.metadata.PypeItMetaData`, None):
            The class holding the metadata for all the frames in this
            PypeIt run.
        raise_error (:obj:`bool`, optional):
            If True, crash out
        cut_cfg (`numpy.ndarray`_, optional):
            Also cut on this restricted configuration (mainly for chk_calibs)

    Returns:
        bool: True if we passed all checks
    """
    if cut_cfg is None:
        cut_cfg = np.ones(len(fitstbl), dtype=bool)
    pass_calib = True
    # Find the science frames
    is_science = fitstbl.find_frames('science')
    # Frame indices
    frame_indx = np.arange(len(fitstbl))

    for i in range(fitstbl.n_calib_groups):
        in_grp = fitstbl.find_calib_group(i)
        grp_science = frame_indx[is_science & in_grp & cut_cfg]
        u_combid = np.unique(fitstbl['comb_id'][grp_science])
        for j, comb_id in enumerate(u_combid):
            frames = np.where(fitstbl['comb_id'] == comb_id)[0]
            calib_ID = int(fitstbl['calib'][frames[0]])

            # Arc, tilt, science
            for ftype in ['arc', 'tilt', 'science', 'trace']:
                rows = fitstbl.find_frames(ftype, calib_ID=calib_ID, index=True)
                if len(rows) == 0:
                    # Fail
                    msg = "No frames of type={} provided. Add them to your PypeIt file if this is a standard run!".format(ftype)
                    pass_calib = False
                    if raise_error:
                        msgs.error(msg)
                    else:
                        msgs.warn(msg)

            # Explore science frame
            for key, ftype in zip(['use_biasimage', 'use_darkimage', 'use_pixelflat', 'use_illumflat'],
                                  ['bias', 'dark', 'pixelflat', 'illumflat']):
                if par['scienceframe']['process'][key]:
                    rows = fitstbl.find_frames(ftype, calib_ID=calib_ID, index=True)
                    if len(rows) == 0:
                        # Allow for pixelflat inserted
                        if ftype == 'pixelflat' and par['calibrations']['flatfield']['pixelflat_file'] is not None:
                            continue
                        # Otherwise fail
                        msg = "No frames of type={} provide for the *{}* processing step. Add them to your PypeIt file!".format(ftype, key)
                        pass_calib = False
                        if raise_error:
                            msgs.error(msg)
                        else:
                            msgs.warn(msg)

    if pass_calib:
        msgs.info("Congrats!!  You passed the calibrations inspection!!")
    return pass_calib
