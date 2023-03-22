"""
Class for guiding calibration object generation in PypeIt.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from pathlib import Path
from datetime import datetime
from copy import deepcopy
from abc import ABCMeta
from collections import Counter

from IPython import embed

import numpy as np
import yaml

from pypeit import __version__
from pypeit import msgs
from pypeit import alignframe
from pypeit import flatfield
from pypeit import edgetrace
from pypeit import slittrace
from pypeit import wavecalib
from pypeit import wavetilts
from pypeit.calibframe import CalibFrame
from pypeit.images import buildimage
from pypeit.metadata import PypeItMetaData
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs.spectrograph import Spectrograph
from pypeit import io
from pypeit import utils


class Calibrations:
    """
    Class designed to guide the generation of calibration images and objects in
    PypeIt.

    It also keeps track of existing processed calibration frames to avoid
    rebuilding those that were generated during this execution of PypeIt.

    Args:
        fitstbl (:class:`pypeit.metadata.PypeItMetaData`, None):
            The class holding the metadata for all the frames in this PypeIt run.
            If None, we are using this class as a glorified dict to hold the objects.
        par (:class:`pypeit.par.pypeitpar.CalibrationsPar`):
            Parameter set defining optional parameters of PypeIt's algorithms
            for Calibrations
        spectrograph (:obj:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph object
        caldir (:obj:`str`, `Path`_):
            Path to write the output calibrations.
        qadir (:obj:`str`, optional):
            Path for quality assessment output.  If not provided, no QA
            plots are saved.
        reuse (:obj:`bool`, optional):
            Instead of reprocessing them, load existing calibration files from
            disk if they exist.
        show (:obj:`bool`, optional):
            Show plots of PypeIt's results as the code progresses.  Requires
            interaction from the users.
        user_slits (:obj:`dict`, optional):
            A limited set of slits selected by the user for analysis.  See
            :func:`~pypeit.slittrace.SlitTraceSet.user_mask`.

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
            Arc-lamp calibration frame.
        alignments (:class:`pypeit.alignframe.Alignments`):
        wv_calib (:class:`pypeit.wavecalib.WaveCalib`):
        slits (:class:`pypeit.slittrace.SlitTraceSet`):

        write_qa
        show
        spectrograph
        par (:class:`pypeit.par.pypeitpar.CalibrationsPar`):
        full_par (:class:`pypeit.par.pypeitpar.PypeItPar`):
        redux_path
        calib_dir
        det
        frame (:obj:`int`):
            0-indexed row of the frame being calibrated in
            :attr:`fitstbl`.
        calib_ID (:obj:`int`):
            calib group ID of the current frame
        user_slits (:obj:`list, optional):
            Identifies a slit or slits to restrict the analysis on
            Used in :func:`get_slits` and propagated beyond

    """
    __metaclass__ = ABCMeta

    @staticmethod
    def get_instance(fitstbl, par, spectrograph, caldir, **kwargs):
        """
        Get the instance of the appropriate subclass of :class:`Calibrations` to
        use for reducing data from the provided ``spectrograph``.  For argument
        descriptions, see :class:`Calibrations`.
        """
        calibclass = MultiSlitCalibrations if spectrograph.pypeline in ['MultiSlit', 'Echelle'] \
                        else IFUCalibrations
        return calibclass(fitstbl, par, spectrograph, caldir, **kwargs)

    def __init__(self, fitstbl, par, spectrograph, caldir, qadir=None,
                 reuse=False, show=False, user_slits=None):

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

        # Calibrations
        self.reuse = reuse
        self.calib_dir = Path(caldir).resolve()
        if not self.calib_dir.exists():
            self.calib_dir.mkdir(parents=True)

        # QA
        self.qa_path = None if qadir is None else Path(qadir).resolve()
        if self.qa_path is not None:
            # TODO: This should only be defined in one place!  Where?...
            qa_png_path = self.qa_path / 'PNGs'
        self.write_qa = self.qa_path is not None
        if self.write_qa and not qa_png_path.exists():
            qa_png_path.mkdir(parents=True)

        # Debugging
        self.show = show

        # Restrict on slits?
        self.user_slits = user_slits

        # Attributes
        self.det = None
        self.frame = None
        self.binning = None

#        self.shape = None

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

        # Steps
        self.steps = []
        self.success = False
        self.failed_step = None

    def _prep_calibrations(self, frametype):
        """
        Find all frames matching a given calibration type and return the
        necessary calibration identifiers.

        Args:
            frametype (:obj:`str`):
                Calibration frame type.  Must be a valid frame type; see
                :func:`~pypeit.core.framematch.valid_frametype`.

        Returns:
            :obj:`tuple`:  Returns a :obj:`list` of image files matching the
            input type, the setup/configuration string, and the calibration group string.
        """
        # NOTE: This will raise an exception if the frametype is not valid!
        valid = framematch.valid_frametype(frametype, raise_error=True)

        # Grab rows and files
        rows = self.fitstbl.find_frames(frametype, calib_ID=self.calib_ID, index=True)
        if len(rows) == 0:
            return [], None, None

        # Grab the setup and calibration id(s)
        if self.par[f'{frametype}frame']['process']['calib_setup_and_bit'] is None:
            setup = self.fitstbl['setup'][rows[0]]
            # Here calib_id is a string of comma-separated values
            calib_id = self.fitstbl['calib'][rows[0]]
        else:
            setup, calib_id \
                    = self.par[f'{frametype}frame']['process']['calib_setup_and_bit'].split('_')
            # Here calib_id is a list of strings
            calib_id = calib_id.split('-')

        return self.fitstbl.frame_paths(rows), setup, CalibFrame.ingest_calib_id(calib_id)

    def set_config(self, frame, det, par=None):
        """
        Specify the parameters of the Calibrations class and reset all
        the internals to None.

        Args:
            frame (:obj:`int`):
                Frame index in the fitstbl
            det (:obj:`int`):
                Detector number
            par (:class:`~pypeit.par.pypeitpar.CalibrationsPar`):
                Parameters used by the calibration procedures.

        """
        # Reset internals to None
        # NOTE: This sets calib_ID so must
        # be done here first before these things are initialized below.

        # Initialize for this setup
        self.frame = frame
        # Find the calibration groups associated with this frame.  Note
        # find_frame_calib_groups *always* returns a list.  Science frames only
        # have one calibration group, but calibration frames can have many.  So
        # for both science and calibration frames, we just set the calibration
        # group to the first in the returned list.
        # TODO: May need to revisit this...
        self.calib_ID = self.fitstbl.find_frame_calib_groups(self.frame)[0]
#        self.calib_ID = int(self.fitstbl['calib'][frame])
        self.det = det
        if par is not None:
            self.par = par
        # Deal with binning
        self.binning = self.fitstbl['binning'][self.frame]

    def get_arc(self):
        """
        Load or generate the arc calibration frame.

        Returns:
            :class:`pypeit.images.buildimage.ArcImage`: The processed
            calibration image.
        """
        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        raw_files, setup, calib_id = self._prep_calibrations('arc')

        if len(raw_files) == 0:
            # There are no arc files, so we're done!
            self.msarc = None
            msgs.warn('No frametype=arc files available!')
            return self.msarc

        # Construct the detector/mosaic identifier
        detname = self.spectrograph.get_det_name(self.det)

        # Construct the expected calibration frame file name
        arc_file = Path(buildimage.ArcImage.construct_file_name(
                            CalibFrame.construct_calib_key(setup, calib_id, detname),
                            calib_dir=self.calib_dir)).resolve()

        # If it exists and we want to reuse it, do so:
        if arc_file.exists() and self.reuse:
            self.msarc = buildimage.ArcImage.from_file(arc_file)
            return self.msarc

        # Otherwise, create the processed file.
        msgs.info(f'Preparing a {buildimage.ArcImage.calib_type} calibration frame.')
        self.msarc = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                    self.par['arcframe'], raw_files,
                                                    bias=self.msbias, bpm=self.msbpm,
                                                    dark=self.msdark, calib_dir=self.calib_dir,
                                                    setup=setup, calib_id=calib_id)
        # Save the result
        self.msarc.to_file()
        # Return it
        return self.msarc

    def get_tiltimg(self):
        """
        Load or generate the tilt calibration frame.

        Returns:
            :class:`pypeit.images.buildimage.TiltImage`: The processed
            calibration image.
        """
        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        raw_files, setup, calib_id = self._prep_calibrations('tilt')

        if len(raw_files) == 0:
            # There are no arc files, so we're done!
            self.mstilt = None
            msgs.warn('No frametype=tilt files available!')
            return self.mstilt

        # Construct the detector/mosaic identifier
        detname = self.spectrograph.get_det_name(self.det)

        # Construct the expected calibration frame file name
        tilt_file = Path(buildimage.TiltImage.construct_file_name(
                            CalibFrame.construct_calib_key(setup, calib_id, detname),
                            calib_dir=self.calib_dir)).resolve()

        # If it exists and we want to reuse it, do so:
        if tilt_file.exists() and self.reuse:
            self.mstilt = buildimage.TiltImage.from_file(tilt_file)
            return self.mstilt

        # Otherwise, create the processed file.
        msgs.info(f'Preparing a {buildimage.TiltImage.calib_type} calibration frame.')
        self.mstilt = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                     self.par['tiltframe'], raw_files,
                                                     bias=self.msbias, bpm=self.msbpm,
                                                     dark=self.msdark, slits=self.slits,
                                                     calib_dir=self.calib_dir, setup=setup,
                                                     calib_id=calib_id)
        # Save the result
        self.mstilt.to_file()
        # Return it
        return self.mstilt

    def get_align(self):
        """
        Load or generate the alignment calibration frame.

        Requires: :attr:`slits`, :attr:`det`, :attr:`par`

        Returns:
            :class:`pypeit.alignframe.Alignments`: The processed alignment
            image.
        """
        # Check for existing data
        if not self._chk_objs(['msbpm', 'slits']):
            msgs.error('Must have the bpm and slits to make the alignments!')

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        raw_files, setup, calib_id = self._prep_calibrations('align')

        if len(raw_files) == 0:
            # There are no arc files, so we're done!
            self.alignments = None
            return self.alignments

        # Construct the detector/mosaic identifier
        detname = self.spectrograph.get_det_name(self.det)

        # Construct the expected calibration frame file name
        align_file = Path(alignframe.Alignments.construct_file_name(
                            CalibFrame.construct_calib_key(setup, calib_id, detname),
                            calib_dir=self.calib_dir)).resolve()

        # If it exists and we want to reuse it, do so:
        if align_file.exists() and self.reuse:
            self.alignments = alignframe.Alignments.from_file(align_file)
            self.alignments.is_synced(self.slits)
            return self.alignments

        # Otherwise, create the processed file.
        msgs.info(f'Preparing a {alignframe.Alignments.calib_type} calibration frame.')
        msalign = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                 self.par['alignframe'], raw_files,
                                                 bias=self.msbias, bpm=self.msbpm,
                                                 dark=self.msdark, calib_dir=self.calib_dir,
                                                 setup=setup, calib_id=calib_id)

        # Extract some header info needed by the algorithm
        # TODO: We have the fitstbl, can we get it from there, instead of having
        # to reread the file?
        binning = self.spectrograph.get_meta_value(raw_files[0], 'binning')

        # Instantiate
        # TODO: From JFH: Do we need the bpm here?  Check that this was in the previous code.
        alignment = alignframe.TraceAlignment(msalign, self.slits, self.spectrograph,
                                              self.par['alignment'], det=self.det, binning=binning,
                                              qa_path=self.qa_path, msbpm=self.msbpm)
        self.alignments = alignment.run(show=self.show)
        self.alignments.to_file()
        return self.alignments

    def get_bias(self):
        """
        Load or generate the bias calibration frame.

        Returns:
            :class:`pypeit.images.buildimage.BiasImage`: The processed
            calibration image.
        """
        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        raw_files, setup, calib_id = self._prep_calibrations('bias')

        if len(raw_files) == 0:
            # There are no bias files, so we're done!
            self.msbias = None
            return self.msbias

        # Construct the detector/mosaic identifier
        detname = self.spectrograph.get_det_name(self.det)

        # Construct the expected calibration frame file name
        bias_file = Path(buildimage.BiasImage.construct_file_name(
                            CalibFrame.construct_calib_key(setup, calib_id, detname),
                            calib_dir=self.calib_dir)).resolve()

        # If it exists and we want to reuse it, do so:
        if bias_file.exists() and self.reuse:
            self.msbias = buildimage.BiasImage.from_file(bias_file)
            return self.msbias

        # Otherwise, create the processed file.
        self.msbias = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                     self.par['biasframe'], raw_files,
                                                     calib_dir=self.calib_dir, setup=setup,
                                                     calib_id=calib_id)
        # Save the result
        self.msbias.to_file()
        # Return it
        return self.msbias

    def get_dark(self):
        """
        Load or generate the dark calibration frame.

        Returns:
            :class:`pypeit.images.buildimage.DarkImage`: The processed
            calibration image.
        """
        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        raw_files, setup, calib_id = self._prep_calibrations('dark')

        if len(raw_files) == 0:
            # There are no dark files, so we're done!
            self.msdark = None
            return self.msdark

        # Construct the detector/mosaic identifier
        detname = self.spectrograph.get_det_name(self.det)

        # Construct the expected calibration frame file name
        dark_file = Path(buildimage.DarkImage.construct_file_name(
                            CalibFrame.construct_calib_key(setup, calib_id, detname),
                            calib_dir=self.calib_dir)).resolve()

        # If it exists and we want to reuse it, do so:
        if dark_file.exists() and self.reuse:
            self.msdark = buildimage.DarkImage.from_file(dark_file)
            return self.msdark

        # TODO: If a bias has been constructed and it will be subtracted from
        # the science images, it should also be subtracted from this image.  If
        # it isn't, subtracting the dark will effectively lead to subtracting
        # the bias twice.

        # TODO: The order is such that the bpm doesn't exist yet.  But calling
        # buildimage_fromlist will create the bpm if it isn't passed.  So
        # calling get_dark then get_bpm unnecessarily creates the bpm twice.  Is
        # there any reason why creation of the bpm should come after the dark,
        # or can we change the order?

        # Otherwise, create the processed file.
        self.msdark = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                     self.par['darkframe'], raw_files,
                                                     bias=self.msbias, calib_dir=self.calib_dir,
                                                     setup=setup, calib_id=calib_id)
        # Save the result
        self.msdark.to_file()
        # Return it
        return self.msdark

    def get_bpm(self):
        """
        Load or generate the bad pixel mask.

        This is primarily a wrapper for
        :func:`~pypeit.spectrographs.spectrograph.Spectrograph.bpm`.

        Returns:
            `numpy.ndarray`_: The bad pixel mask, which should match the shape
            and orientation of a *trimmed* and PypeIt-oriented science image!
        """
        # Check internals
        self._chk_set(['par', 'det'])
        # Build it
        self.msbpm = self.spectrograph.bpm(self.fitstbl.frame_paths(self.frame), self.det,
                                           msbias=self.msbias if self.par['bpm_usebias'] else None)
        # Return
        return self.msbpm

    def get_flats(self):
        """
        Load or generate the flat-field calibration images.

        Returns:
            :class:`~pypeit.flatfield.FlatImages`: The processed calibration
            image.
        """
        # Check for existing data
        if not self._chk_objs(['msarc', 'msbpm', 'slits', 'wv_calib']):
            msgs.warn('Must have the arc, bpm, slits, and wv_calib defined to make flats!  '
                      'Skipping and may crash down the line')
            # TODO: Why is this an empty object and not None?
            self.flatimages = flatfield.FlatImages()
            return self.flatimages

        # Slit and tilt traces are required to flat-field the data
        if not self._chk_objs(['slits', 'wavetilts']):
            # TODO: Why doesn't this fault?
            msgs.warn('Flats were requested, but there are quantities missing necessary to '
                      'create flats.  Proceeding without flat fielding....')
            # TODO: Why is this an empty object and not None?
            self.flatimages = flatfield.FlatImages()
            return self.flatimages

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        raw_illum_files, setup, calib_id = self._prep_calibrations('illumflat')
        # Using _prep_calibrations is a bit overkill here...
        raw_pixel_files = self._prep_calibrations('pixelflat')[0]
        raw_lampoff_files = self._prep_calibrations('lampoffflats')[0]

        # Construct the detector/mosaic identifier
        detname = self.spectrograph.get_det_name(self.det)

        # Construct the expected calibration frame file name
        flat_file = Path(flatfield.FlatImages.construct_file_name(
                            CalibFrame.construct_calib_key(setup, calib_id, detname),
                            calib_dir=self.calib_dir)).resolve()

        # If it exists and we want to reuse it, do so:
        if flat_file.exists() and self.reuse:
            self.flatimages = flatfield.FlatImages.from_file(flat_file)
            self.flatimages.is_synced(self.slits)
            # Load user defined files
            if self.par['flatfield']['pixelflat_file'] is not None:
                # Load
                msgs.info(f'Using user-defined file: {self.par["flatfield"]["pixelflat_file"]}')
                with io.fits_open(self.par['flatfield']['pixelflat_file']) as hdu:
                    nrm_image = flatfield.FlatImages(pixelflat_norm=hdu[self.det].data)
                    self.flatimages = flatfield.merge(self.flatimages, nrm_image)
            # update slits
            self.slits.mask_flats(self.flatimages)
            return self.flatimages

        # Generate the image
        pixelflatImages, illumflatImages = None, None
        lampoff_flat = None
        # Check if the image files are the same
        pix_is_illum = Counter(illum_image_files) == Counter(pixflat_image_files)
        if len(raw_pixel_files) > 0:
            msgs.info('Creating pixel-flat calibration frame using files: ')
            for f in pixflat_image_files:
                msgs.prindent(f'{Path(f).name}')
            pixel_flat = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                        self.par['pixelflatframe'],
                                                        raw_pixel_files, dark=self.msdark,
                                                        bias=self.msbias, bpm=self.msbpm,
                                                        calib_dir=self.calib_dir, setup=setup,
                                                        calib_id=calib_id)
            if len(raw_lampoff_files) > 0:
                msgs.info('Subtracting lamp off flats using files: ')
                for f in raw_lampoff_files:
                    msgs.prindent(f'{Path(f).name}')
                lampoff_flat = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                              self.par['lampoffflatsframe'],
                                                              raw_lampoff_files,
                                                              dark=self.msdark, bias=self.msbias,
                                                              bpm=self.msbpm,
                                                              calib_dir=self.calib_dir, setup=setup,
                                                              calib_id=calib_id)
                pixel_flat = pixel_flat.sub(lampoff_flat)

            # Initialise the pixel flat
            pixelFlatField = flatfield.FlatField(pixel_flat, self.spectrograph,
                                                 self.par['flatfield'], self.slits, self.wavetilts,
                                                 self.wv_calib, qa_path=self.qa_path,
                                                 calib_key=pixel_flat.calib_key)
            # Generate
            pixelflatImages = pixelFlatField.run(doqa=self.write_qa, show=self.show)
            # Set flatimages in case we want to apply the pixel-to-pixel sensitivity corrections to the illumflat
            self.flatimages = pixelflatImages

        # Only build illum_flat if the input files are different from the pixel flat
        if not pix_is_illum and len(raw_illum_files) > 0:
            msgs.info('Creating illumination flat calibration frame using files: ')
            for f in illum_image_files:
                msgs.prindent(f'{Path(f).name}')
            illum_flat = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                        self.par['illumflatframe'], raw_illum_files,
                                                        dark=self.msdark, bias=self.msbias,
                                                        flatimages=self.flatimages, bpm=self.msbpm,
                                                        calib_dir=self.calib_dir, setup=setup,
                                                        calib_id=calib_id)
            if len(raw_lampoff_files) > 0:
                msgs.info('Subtracting lamp off flats using files: ')
                for f in raw_lampoff_files:
                    msgs.prindent(f'{Path(f).name}')
                if lampoff_flat is None:
                    lampoff_flat = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                                  self.par['lampoffflatsframe'],
                                                                  raw_lampoff_files,
                                                                  dark=self.msdark,
                                                                  bias=self.msbias, bpm=self.msbpm,
                                                                  calib_dir=self.calib_dir,
                                                                  setup=setup, calib_id=calib_id)
                illum_flat = illum_flat.sub(lampoff_flat)

            # Initialise the pixel flat
            illumFlatField = flatfield.FlatField(illum_flat, self.spectrograph,
                                                 self.par['flatfield'], self.slits, self.wavetilts,
                                                 self.wv_calib, spat_illum_only=True,
                                                 qa_path=self.qa_path,
                                                 calib_key=illum_flat.calib_key)
            # Generate
            illumflatImages = illumFlatField.run(doqa=self.write_qa, show=self.show)

        # Merge the illum flat with the pixel flat
        if pixelflatImages is not None:
            # Combine the pixelflat and illumflat parameters into flatimages.
            # This will merge the attributes of pixelflatImages that are not None
            # with the attributes of illflatImages that are not None. Default is
            # to take pixelflatImages.
            self.flatimages = flatfield.merge(pixelflatImages, illumflatImages)
        else:
            # No pixel flat, but there might be an illumflat. This will mean that
            # the attributes prefixed with 'pixelflat_' will all be None.
            self.flatimages = illumflatImages

        if self.flatimages is not None:
            # Save flat images
            self.flatimages.to_file()
            # Save slits too, in case they were tweaked
            self.slits.to_file()

        # 3) Load user-supplied images
        # NOTE: These are the *final* images, not just a stack, and it will
        # over-ride what is generated below (if generated).

        # TODO: Why is this done after writing the image above?  If we did,  we
        # wouldn't need to re-read the user-provided file when ingesting the
        # existing flat file.  Is this to allow the user to change the pixel
        # flat file?  Should we allow the user to change the file name?
        if self.par['flatfield']['pixelflat_file'] is not None:
            # Load
            msgs.info(f'Using user-defined file: {self.par["flatfield"]["pixelflat_file"]}')
            with io.fits_open(self.par['flatfield']['pixelflat_file']) as hdu:
                self.flatimages = flatfield.merge(self.flatimages,
                                        flatfield.FlatImages(pixelflat_norm=hdu[self.det].data))

        return self.flatimages

    def get_slits(self):
        """
        Load or generate the definition of the slit boundaries.

        Returns:
            :class:`pypeit.slittrace.SlitTraceSet`: Traces of the
            slit edges; also kept internally as :attr:`slits`.
        """
        # Check for existing data
        if not self._chk_objs(['msbpm']):
            return None

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        raw_trace_files, setup, calib_id = self._prep_calibrations('trace')
        raw_lampoff_files = self._prep_calibrations('lampoffflats')[0]

        # Construct the detector/mosaic identifier
        detname = self.spectrograph.get_det_name(self.det)

        # Construct the expected calibration frame file name
        slits_file = Path(slittrace.SlitTraceSet.construct_file_name(
                            CalibFrame.construct_calib_key(setup, calib_id, detname),
                            calib_dir=self.calib_dir)).resolve()

        # If it exists and we want to reuse it, do so:
        if slits_file.exists() and self.reuse:
            self.slits = slittrace.SlitTraceSet.from_file(slits_file)
            self.slits.mask = self.slits.mask_init.copy()
            if self.user_slits is not None:
                self.slits.user_mask(detname, self.user_slits)            
            return self.slits

        if len(raw_trace_files) == 0:
            self.slits = None
            msgs.warn('No frametype=trace files to build slits')
            return self.slits

        # Slits don't exist or we're not resusing them
        edges_file = Path(edgetrace.EdgeTraceSet.construct_file_name(
                            CalibFrame.construct_calib_key(setup, calib_id, detname),
                            calib_dir=self.calib_dir)).resolve()
        # Reuse an exising calibration frame?
        if edges_file.exists() and self.reuse:
            self.slits = edgetrace.EdgeTraceSet.from_file(edges_file).get_slits()
            self.slits.to_file()
            if self.user_slits is not None:
                self.slits.user_mask(detname, self.user_slits)            
            return self.slits

        # Build the trace image
        msgs.info('Creating edge tracing calibration frame using files: ')
        for f in raw_trace_files:
            msgs.prindent(f'{Path(f).name}')
        self.traceImage = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                         self.par['traceframe'], raw_trace_files,
                                                         bias=self.msbias, bpm=self.msbpm,
                                                         dark=self.msdark,
                                                         calib_dir=self.calib_dir,
                                                         setup=setup, calib_id=calib_id)
        if len(raw_lampoff_files) > 0:
            msgs.info('Subtracting lamp off flats using files: ')
            for f in raw_lampoff_files:
                msgs.prindent(f'{Path(f).name}')
            lampoff_flat = buildimage.buildimage_fromlist(self.spectrograph, self.det,
                                                          self.par['lampoffflatsframe'],
                                                          raw_lampoff_files, dark=self.msdark,
                                                          bias=self.msbias, bpm=self.msbpm,
                                                          calib_dir=self.calib_dir, setup=setup,
                                                          calib_id=calib_id)
            self.traceImage = self.traceImage.sub(lampoff_flat)

        self.edges = edgetrace.EdgeTraceSet(self.traceImage, self.spectrograph,
                                            self.par['slitedges'], auto=True)
        if not self.edges.success:
            self.success = False
            return None
        self.edges.to_file()

        # Show the result if requested
        if self.show:
            self.edges.show(in_ginga=True)

        # Get the slits from the result of the edge tracing, delete
        # the edges object, and save the slits, if requested
        self.slits = self.edges.get_slits()
        self.edges = None
        self.slits.to_file()
        if self.user_slits is not None:
            self.slits.user_mask(detname, self.user_slits)            

        # User mask?
#        if self.user_slits is not None:
#            # Parse the DET/MSC name
#            if isinstance(self.det, tuple):
#                detname = self.spectrograph.list_detectors(
#                mosaic=True)[self.spectrograph.allowed_mosaics.index(self.det)]
#            elif isinstance(self.det, int):
#                detname = self.spectrograph.list_detectors()[self.det-1]
#            else:
#                msgs.error("Bad type for self.det")
#
#            self.slits.user_mask(detname, self.user_slits)

        return self.slits

    def get_wv_calib(self):
        """
        Load or generate the 1D wavelength calibrations

        Returns:
            :class:`~pypeit.wavecalib.WaveCalib`: Object containing wavelength
            calibrations and the updated slit mask array.
        """
        # No wavelength calibration requested
        if self.par['wavelengths']['reference'] == 'pixel':
            msgs.info('Wavelength "reference" parameter set to "pixel"; no wavelength '
                      'calibration will be performed.')
            self.wv_calib = None
            return self.wv_calib

        # Check for existing data
        if not self._chk_objs(['msarc', 'msbpm', 'slits']):
            msgs.warn('Not enough information to load/generate the wavelength calibration. '
                      'Skipping and may crash down the line')
            return None

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Prep
        _, setup, calib_id = self._prep_calibrations('arc')

        # Construct the detector/mosaic identifier
        detname = self.spectrograph.get_det_name(self.det)

        # Grab arc binning (may be different from science!)
        # TODO : Do this internally when we have a wv_calib DataContainer
        binspec, binspat = parse.parse_binning(self.msarc.detector.binning)

        # Construct the expected calibration frame file name
        wvcalib_file = Path(wavecalib.WaveCalib.construct_file_name(
                            CalibFrame.construct_calib_key(setup, calib_id, detname),
                            calib_dir=self.calib_dir)).resolve()

        # If it exists and we want to reuse it, do so:
        if wvcalib_file.exists() and self.reuse:
            self.wv_calib = wavecalib.WaveCalib.from_file(wvcalib_file)
            self.wv_calib.chk_synced(self.slits)
            self.slits.mask_wvcalib(self.wv_calib)
            return self.wv_calib

        # Determine lamp list to use for wavecalib
        # Find all the arc frames in this calibration group
        is_arc = self.fitstbl.find_frames('arc', calib_ID=self.calib_ID)
        lamps = self.spectrograph.get_lamps(self.fitstbl[is_arc]) \
                    if self.par['wavelengths']['lamps'] == ['use_header'] \
                    else self.par['wavelengths']['lamps']
        meta_dict = dict(self.fitstbl[is_arc][0]) \
                    if self.spectrograph.pypeline == 'Echelle' \
                        and not self.spectrograph.ech_fixed_format else None
        # Instantiate
        # TODO: Pull out and pass only the necessary parts of meta_dict to
        # this, or include the relevant parts as parameters.  See comments
        # in PRs #1454 and #1476 on this.
        msgs.info(f'Preparing a {wavecalib.WaveCalib.calib_type} calibration frame.')
        self.waveCalib = wavecalib.BuildWaveCalib(self.msarc, self.slits, self.spectrograph,
                                                  self.par['wavelengths'], lamps,
                                                  meta_dict=meta_dict, binspectral=binspec,
                                                  det=self.det, qa_path=self.qa_path)
        self.wv_calib = self.waveCalib.run(skip_QA=(not self.write_qa))
        # If orders were found, save slits to disk
        if self.spectrograph.pypeline == 'Echelle' and not self.spectrograph.ech_fixed_format:
            self.slits.to_file()
        # Save calibration frame
        self.wv_calib.to_file()

        # Return
        return self.wv_calib

    def get_tilts(self):
        """
        Load or generate the wavelength tilts calibration frame

        Returns:
            :class:`~pypeit.wavetilts.WaveTilts`: Object containing the
            wavelength tilt calibration.
        """
        # Check for existing data
        # TODO: add mstilt_inmask to this list when it gets implemented.
        if not self._chk_objs(['mstilt', 'msbpm', 'slits', 'wv_calib']):
            msgs.warn('Do not have all the necessary objects for tilts.  Skipping and may crash '
                      'down the line.')
            return None

        # Check internals
        self._chk_set(['det', 'calib_ID', 'par'])

        # Load up?
        _, setup, calib_id = self._prep_calibrations('tilt')
        detname = self.spectrograph.get_det_name(self.det)
        tilts_file = Path(wavetilts.WaveTilts.construct_file_name(
                            CalibFrame.construct_calib_key(setup, calib_id, detname),
                            calib_dir=self.calib_dir)).resolve()

        if tilts_file.exists() and self.reuse:
            self.wavetilts = wavetilts.WaveTilts.from_file(tilts_file)
            self.wavetilts.is_synced(self.slits)
            self.slits.mask_wavetilts(self.wavetilts)
            return self.wavetilts

        # Get flexure
        _spat_flexure = self.mstilt.spat_flexure \
            if self.par['tiltframe']['process']['spat_flexure_correct'] else None

        # Build
        buildwaveTilts = wavetilts.BuildWaveTilts(
            self.mstilt, self.slits, self.spectrograph, self.par['tilts'],
            self.par['wavelengths'], det=self.det, qa_path=self.qa_path,
            spat_flexure=_spat_flexure)

        # TODO still need to deal with syntax for LRIS ghosts. Maybe we don't need it
        self.wavetilts = buildwaveTilts.run(doqa=self.write_qa, show=self.show)
        self.wavetilts.to_file()
        return self.wavetilts

    def run_the_steps(self):
        """
        Run full the full recipe of calibration steps.
        """
        self.success = True
        for step in self.steps:
            getattr(self, f'get_{step}')()
            if not self.success:
                self.failed_step = f'get_{step}'
                return
        msgs.info("Calibration complete and/or fully loaded!")
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

    @staticmethod
    def association_summary(ofile, fitstbl, spectrograph, caldir, subset=None, det=None,
                            overwrite=False):
        """
        Write a file listing the associations between the processed calibration
        frames and their source raw files for every setup and every calibration
        group.

        Args:
            ofile (:obj:`str`, `Path`_):
                Full path to the output file.
            fitstbl (:class:`pypeit.metadata.PypeItMetaData`):
                The class holding the metadata for all the frames to process.
            spectrograph (:obj:`pypeit.spectrographs.spectrograph.Spectrograph`):
                Spectrograph object
            caldir (:obj:`str`, `Path`_):
                Path for the processed calibration frames.
            subset (`numpy.ndarray`_, optional):
                A boolean array selecting a subset of rows from ``fitstbl`` for
                output.
            det (:obj:`int`, :obj:`tuple`, optional):
                The specific detector (or mosaic) to use when constructing the
                output processed calibration group file names.  If None, a
                placeholder is used.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file of the same name.
        """
        if fitstbl.calib_groups is None:
            msgs.error('Calibration groups have not been defined!')

        _ofile = Path(ofile).resolve()
        if _ofile.exists() and not overwrite:
            msgs.error(f'{_ofile} exists!  To overwrite, set overwrite=True.')

        _caldir = str(Path(caldir).resolve())

        # This defines the classes used by each frametype that results in an
        # output calibration frame:
        frame_calibrations = {'align': [alignframe.Alignments],
                              'arc': [buildimage.ArcImage, wavecalib.WaveCalib],
                              'bias': [buildimage.BiasImage],
                              'dark': [buildimage.DarkImage],
                              'pixelflat': [flatfield.FlatImages],
                              'illumflat': [flatfield.FlatImages],
                              'lampoffflats': [flatfield.FlatImages], 
                              'trace': [edgetrace.EdgeTraceSet, slittrace.SlitTraceSet],
                              'tilt': [buildimage.TiltImage, wavetilts.WaveTilts]
                             }

        # Get the name of the detector/mosaic
        detname = spectrograph.get_det_name(1 if det is None else det)

        # Find the unique configuations in the metaddata
        setups = fitstbl.unique_configurations(copy=True, rm_none=True)

        # Subset to output
        if subset is None:
            subset = np.ones(len(fitstbl), dtype=bool)

        cfg = {}
        # Iterate through each setup
        for setup in setups.keys():
            cfg[setup] = {}
            cfg[setup]['--'] = deepcopy(setups[setup])
            in_setup = fitstbl.find_configuration(setup) & subset
            if not any(in_setup):
                continue
            # Iterate through each calibration group
            for calib_ID in fitstbl.calib_groups:
                # Find all the frames in this calibration group
                in_grp = fitstbl.find_calib_group(calib_ID) & in_setup
                if not any(in_grp):
                    continue

                cfg[setup][calib_ID] = {}

                # First list the science and standard (and sky?) frames in this
                # calibration group
                for frametype in ['science', 'standard']:
                    indx = fitstbl.find_frames(frametype) & in_grp
                    if not any(indx):
                        continue
                    cfg[setup][calib_ID][frametype] = fitstbl.frame_paths(indx)

                # Iterate through each frame type
                for frametype, calib_classes in frame_calibrations.items():
                    indx = fitstbl.find_frames(frametype) & in_grp
                    if not any(indx):
                        continue
                    if not all(fitstbl['calib'][indx] == fitstbl['calib'][indx][0]):
                        msgs.error(f'CODING ERROR: All {frametype} frames in group {calib_ID} '
                                   'are not all associated with the same subset of calibration '
                                   'groups; calib for the first file is '
                                   f'{fitstbl["calib"][indx][0]}.')
                    calib_key = CalibFrame.construct_calib_key(setup, fitstbl['calib'][indx][0],
                                                               detname)
                    cfg[setup][calib_ID][frametype] = {}
                    cfg[setup][calib_ID][frametype]['raw'] = fitstbl.frame_paths(indx)
                    cfg[setup][calib_ID][frametype]['proc'] \
                            = [calib_class.construct_file_name(calib_key, calib_dir=_caldir)
                                for calib_class in frame_calibrations[frametype]]

        # Write it
        with open(_ofile, 'w') as ff:
            ff.write('# Auto-generated calibration association file using PypeIt version: '
                     f' {__version__}\n')
            ff.write(f'# UTC {datetime.utcnow().isoformat(timespec="milliseconds")}\n')
            if det is None:
                ff.write(f'# NOTE: {detname} is a placeholder for the reduced detectors/mosaics\n')
            ff.write(yaml.dump(utils.yamlify(cfg)))
        msgs.info(f'Calibration association file written to: {_ofile}')


class MultiSlitCalibrations(Calibrations):
    """
    Calibration class for performing multi-slit calibrations (and also long-slit
    and echelle).  See :class:`~pypeit.calibrations.Calibrations` for arguments.

    .. note::

        Calibrations are not sufficiently different yet to warrant a different
        class for echelle reductions.  This may change if a different order is
        eventually required for the set of processing steps (see
        :func:`default_steps`).
    """
    def __init__(self, fitstbl, par, spectrograph, caldir, **kwargs):
        super().__init__(fitstbl, par, spectrograph, caldir, **kwargs)
        self.steps = MultiSlitCalibrations.default_steps()

    @staticmethod
    def default_steps():
        """
        This defines the calibration steps and their order.

        Returns:
            :obj:`list`: Calibration steps, in order of execution.
        """
        # Order matters!  And the name must match a viable "get_{step}" method
        # in Calibrations.
        # TODO: Does the bpm need to be done after the dark?
        return ['bias', 'dark', 'bpm', 'slits', 'arc', 'tiltimg', 'wv_calib', 'tilts', 'flats']


class IFUCalibrations(Calibrations):
    """
    Child of Calibrations class for performing IFU calibrations.
    See :class:`pypeit.calibrations.Calibrations` for arguments.
    """

    def __init__(self, fitstbl, par, spectrograph, caldir, **kwargs):
        super().__init__(fitstbl, par, spectrograph, caldir, **kwargs)
        self.steps = IFUCalibrations.default_steps()

    @staticmethod
    def default_steps():
        """
        This defines the steps for calibrations and their order

        Returns:
            list: Calibration steps, in order of execution
        """
        # Order matters!
        return ['bias', 'dark', 'bpm', 'arc', 'tiltimg', 'slits', 'wv_calib', 'tilts', 'align',
                'flats']


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

    for calib_ID in fitstbl.calib_groups:
        in_grp = fitstbl.find_calib_group(calib_ID)
        if not np.any(is_science & in_grp & cut_cfg):
            continue
        grp_science = frame_indx[is_science & in_grp & cut_cfg]
        u_combid = np.unique(fitstbl['comb_id'][grp_science])
        for j, comb_id in enumerate(u_combid):
            frames = np.where(fitstbl['comb_id'] == comb_id)[0]
            # Arc, tilt, science
            for ftype in ['arc', 'tilt', 'science', 'trace']:
                rows = fitstbl.find_frames(ftype, calib_ID=calib_ID, index=True)
                if len(rows) == 0:
                    # Fail
                    msg = f'No frames of type={ftype} provided. Add them to your PypeIt file ' \
                          'if this is a standard run!'
                    pass_calib = False
                    if raise_error:
                        msgs.error(msg)
                    else:
                        msgs.warn(msg)

            # Explore science frame
            for key, ftype in zip(['use_biasimage', 'use_darkimage', 'use_pixelflat',
                                   'use_illumflat'], ['bias', 'dark', 'pixelflat', 'illumflat']):
                if par['scienceframe']['process'][key]:
                    rows = fitstbl.find_frames(ftype, calib_ID=calib_ID, index=True)
                    if len(rows) == 0:
                        # Allow for pixelflat inserted
                        if ftype == 'pixelflat' \
                                and par['calibrations']['flatfield']['pixelflat_file'] is not None:
                            continue
                        # Otherwise fail
                        msg = f'No frames of type={ftype} provide for the *{key}* processing ' \
                              'step. Add them to your PypeIt file!'
                        pass_calib = False
                        if raise_error:
                            msgs.error(msg)
                        else:
                            msgs.warn(msg)

    if pass_calib:
        msgs.info("Congrats!!  You passed the calibrations inspection!!")
    return pass_calib



