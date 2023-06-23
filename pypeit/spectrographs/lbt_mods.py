"""
Module for LBT/MODS specific methods.

.. include:: ../include/links.rst
"""
import glob
import numpy as np
from astropy.io import fits

from pypeit import msgs
from pypeit import telescopes
from pypeit import io
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container

# TODO: FW: test MODS1B and MODS2B

class LBTMODSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Shane/Kast specific code
    """
    ndet = 1
    telescope = telescopes.LBTTelescopePar()
    url = 'https://scienceops.lbto.org/mods/'

#    def __init__(self):
#        super().__init__()
#        self.timeunit = 'isot'

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Scienceimage default parameters
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['pixelflatframe']['exprng'] = [0, None]
        par['calibrations']['traceframe']['exprng'] = [0, None]
        par['calibrations']['arcframe']['exprng'] = [None, None]
        par['calibrations']['standardframe']['exprng'] = [1, 200]
        par['scienceframe']['exprng'] = [200, None]

        # Do not sigmaclip the arc frames for better Arc and better wavecalib
        par['calibrations']['arcframe']['process']['clip'] = False
        # Do not sigmaclip the tilt frames
        par['calibrations']['tiltframe']['process']['clip'] = False

        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='OBJRA')
        self.meta['dec'] = dict(ext=0, card='OBJDEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(ext=0, card='MASKNAME')
        self.meta['binning'] = dict(card=None, compound=True)
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        self.meta['dispname'] = dict(ext=0, card='GRATNAME')
        self.meta['dichroic'] = dict(ext=0, card='FILTNAME')
        self.meta['idname'] = dict(ext=0, card='IMAGETYP')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')

    def compound_meta(self, headarr, meta_key):
        """
        Methods to generate metadata requiring interpretation of the header
        data, instead of simply reading the value of a header card.

        Args:
            headarr (:obj:`list`):
                List of `astropy.io.fits.Header`_ objects.
            meta_key (:obj:`str`):
                Metadata keyword to construct.

        Returns:
            object: Metadata value read from the header(s).
        """
        if meta_key == 'binning':
            binspatial, binspec = parse.parse_binning(np.array([headarr[0]['CCDXBIN'], headarr[0]['CCDYBIN']]))
            binning = parse.binning2string(binspatial, binspec)
            return binning
        msgs.error("Not ready for this compound meta")

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        This list is used by :class:`~pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            :obj:`list`: List of keywords of data pulled from file headers
            and used to constuct the :class:`~pypeit.metadata.PypeItMetaData`
            object.
        """
        # decker is not included because standards are usually taken with a 5" slit and arc using 0.8" slit
        return ['dispname', 'binning' ]

    def raw_header_cards(self):
        """
        Return additional raw header cards to be propagated in
        downstream output files for configuration identification.

        The list of raw data FITS keywords should be those used to populate
        the :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.configuration_keys`
        or are used in :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.config_specific_par`
        for a particular spectrograph, if different from the name of the
        PypeIt metadata keyword.

        This list is used by :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.subheader_for_spec`
        to include additional FITS keywords in downstream output files.

        Returns:
            :obj:`list`: List of keywords from the raw data files that should
            be propagated in output files.
        """
        return ['GRATNAME', 'CCDXBIN', 'CCDYBIN']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.

        Args:
            ftype (:obj:`str`):
                Type of frame to check. Must be a valid frame type; see
                frame-type :ref:`frame_type_defs`.
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for one or more frames to check.
            exprng (:obj:`list`, optional):
                Range in the allowed exposure time for a frame of type
                ``ftype``. See
                :func:`pypeit.core.framematch.check_frame_exptime`.

        Returns:
            `numpy.ndarray`_: Boolean array with the flags selecting the
            exposures in ``fitstbl`` that are ``ftype`` type frames.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['science', 'standard']:
            return good_exp & (fitstbl['idname'] == 'OBJECT') & (fitstbl['ra'] != 'none') \
                   & (fitstbl['dispname'] != 'Flat')
        if ftype == 'bias':
            return good_exp  & (fitstbl['idname'] == 'BIAS')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            # Flats and trace frames are typed together
            return good_exp  & (fitstbl['idname'] == 'FLAT') & (fitstbl['decker'] != 'Imaging')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'COMP') & (fitstbl['dispname'] != 'Flat')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces
        that are key for image processing.

        Parameters
        ----------
        raw_file : :obj:`str`
            File to read
        det : :obj:`int`
            1-indexed detector to read

        Returns
        -------
        detector_par : :class:`pypeit.images.detector_container.DetectorContainer`
            Detector metadata parameters.
        raw_img : `numpy.ndarray`_
            Raw image for this detector.
        hdu : `astropy.io.fits.HDUList`_
            Opened fits file
        exptime : :obj:`float`
            Exposure time read from the file header
        rawdatasec_img : `numpy.ndarray`_
            Data (Science) section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.
        oscansec_img : `numpy.ndarray`_
            Overscan section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.
        """
        # Check for file; allow for extra .gz, etc. suffix
        fil = glob.glob(raw_file + '*')
        if len(fil) != 1:
            msgs.error("Found {:d} files matching {:s}".format(len(fil)))

        # Read
        msgs.info("Reading LBT/MODS file: {:s}".format(fil[0]))
        hdu = io.fits_open(fil[0])
        head = hdu[0].header

        # TODO These parameters should probably be stored in the detector par

        # Number of amplifiers (could pull from DetectorPar but this avoids needing the spectrograph, e.g. view_fits)
        detector_par = self.get_detector_par(det if det is not None else 1, hdu=hdu)
        numamp = detector_par['numamplifiers']

        # get the x and y binning factors...
        xbin, ybin = head['CCDXBIN'], head['CCDYBIN']

        datasize = head['DETSIZE'] # Unbinned size of detector full array
        _, nx_full, _, ny_full = np.array(parse.load_sections(datasize, fmt_iraf=False)).flatten()

        # Determine the size of the output array...
        nx, ny = int(nx_full / xbin), int(ny_full / ybin)
        nbias1 = 48
        nbias2 = 8240

        # allocate output array...
        array = hdu[0].data.T * 1.0 ## Convert to float in order to get it processed with procimg.py
        rawdatasec_img = np.zeros_like(array, dtype=int)
        oscansec_img = np.zeros_like(array, dtype=int)

        ## allocate datasec and oscansec to the image
        # apm 1
        rawdatasec_img[int(nbias1/xbin):int(nx/2), :int(ny/2)] = 1
        oscansec_img[1:int(nbias1/xbin), :int(ny/2)] = 1 # exclude the first pixel since it always has problem

        # apm 2
        rawdatasec_img[int(nx/2):int(nbias2/xbin), :int(ny/2)] = 2
        oscansec_img[int(nbias2/xbin):nx-1, :int(ny/2)] = 2 # exclude the last pixel since it always has problem

        # apm 3
        rawdatasec_img[int(nbias1/xbin):int(nx/2), int(ny/2):] = 3
        oscansec_img[1:int(nbias1/xbin), int(ny/2):] = 3 # exclude the first pixel since it always has problem

        # apm 4
        rawdatasec_img[int(nx/2):int(nbias2/xbin), int(ny/2):] = 4
        oscansec_img[int(nbias2/xbin):nx-1, int(ny/2):] = 4 # exclude the last pixel since it always has problem

        # Need the exposure time
        exptime = hdu[self.meta['exptime']['ext']].header[self.meta['exptime']['card']]
        # Return, transposing array back to orient the overscan properly
        return detector_par,np.flipud(array), hdu, exptime, np.flipud(rawdatasec_img), np.flipud(oscansec_img)


class LBTMODS1RSpectrograph(LBTMODSSpectrograph):
    """
    Child to handle LBT/MODS1R specific code
    """
    name = 'lbt_mods1r'
    camera = 'MODS1R'
    header_name = 'MODS1R'
    supported = True
    comment = 'MODS-I red spectrometer'

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        # Binning
        binning = '1,1' if hdu is None \
                    else f"{hdu[0].header['CCDXBIN']},{hdu[0].header['CCDYBIN']}"

        # Detector 1
        detector_dict = dict(
            binning= binning,
            det=1,
            dataext         = 0,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.123,
            darkcurr        = 0.4,
            saturation      = 65535.,
            nonlinear       = 0.99,
            mincounts       = -1e10,
            numamplifiers   = 4,
            gain            = np.atleast_1d([2.38,2.50,2.46,2.81]),
            ronoise         = np.atleast_1d([3.78,4.04,4.74,4.14]),
# TODO: The raw image reader sets these up by hand
#            datasec         = np.atleast_1d('[:,:]'),
#            oscansec        = np.atleast_1d('[:,:]')
            )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        par['flexure']['spec_method'] = 'boxcar'

        # 1D wavelength solution
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['rms_threshold'] = 0.4
        par['calibrations']['wavelengths']['fwhm'] = 10.
        #par['calibrations']['wavelengths']['lamps'] = ['XeI','ArII','ArI','NeI','KrI']]
        par['calibrations']['wavelengths']['lamps'] = ['ArI','NeI','KrI','XeI']
        #par['calibrations']['wavelengths']['lamps'] = ['OH_MODS']
        par['calibrations']['wavelengths']['n_first'] = 3
        par['calibrations']['wavelengths']['match_toler'] = 2.5

        # slit
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['edge_thresh'] = 100.

        # Set wave tilts order
        par['calibrations']['tilts']['spat_order'] = 5
        par['calibrations']['tilts']['spec_order'] = 5
        par['calibrations']['tilts']['maxdev_tracefit'] = 0.02
        par['calibrations']['tilts']['maxdev2d'] = 0.02

        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            scifile (:obj:`str`):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        par = super().config_specific_par(scifile, inp_par=inp_par)

        if self.get_meta_value(scifile, 'dispname') == 'G670L':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'lbt_mods1r_red.fits'

        return par

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask.

        Even though they are both optional, either the precise shape for
        the image (``shape``) or an example file that can be read to get
        the shape (``filename`` using :func:`get_image_shape`) *must* be
        provided.

        Args:
            filename (:obj:`str` or None):
                An example file to use to get the image shape.
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.
            shape (tuple, optional):
                Processed image shape
                Required if filename is None
                Ignored if filename is not None
            msbias (`numpy.ndarray`_, optional):
                Processed bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)

        msgs.info("Using hard-coded BPM for  MODS1R")

        # TODO: Fix this
        # Get the binning
        hdu = io.fits_open(filename)
        header = hdu[0].header
        xbin, ybin = header['CCDXBIN'], header['CCDYBIN']
        hdu.close()

        # Apply the mask
        bpm_img[6278//xbin:6289//xbin, 1544//ybin:1634//ybin] = 1
        bpm_img[4202//xbin:4204//xbin, 1474//ybin:1544//ybin] = 1
        bpm_img[3551//xbin:3558//xbin, 2391//ybin:2903//ybin] = 1
        bpm_img[3553//xbin:3558//xbin, 1454//ybin:1544//ybin] = 1

        bpm_img[5650//xbin, 1280//ybin:1544//ybin] = 1
        bpm_img[4780//xbin, 1406//ybin:1536//ybin] = 1
        bpm_img[3554//xbin, 1544//ybin:2392//ybin] = 1
        bpm_img[163//xbin, 1544//ybin:1963//ybin] = 1

        return bpm_img


class LBTMODS1BSpectrograph(LBTMODSSpectrograph):
    """
    Child to handle LBT/MODS1R specific code
    """

    name = 'lbt_mods1b'
    camera = 'MODS1B'
    header_name = 'MODS1B'
    supported = True
    comment = 'MODS-I blue spectrometer'

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        binning = '1,1' if hdu is None \
                    else f"{hdu[0].header['CCDXBIN']},{hdu[0].header['CCDYBIN']}"

        # Detector 1
        detector_dict = dict(
            binning= binning,
            det=1,
            dataext         = 0,
            specaxis        = 0,
            specflip        = True,
            spatflip        = False,
            platescale      = 0.120,
            darkcurr        = 0.5,
            saturation      = 65535.,
            nonlinear       = 0.99,
            mincounts       = -1e10,
            numamplifiers   = 4,
            gain            = np.atleast_1d([2.55,1.91,2.09,2.02]),
            ronoise         = np.atleast_1d([3.41,2.93,2.92,2.76]),
# TODO: The raw image reader sets these up by hand
#            datasec         = np.atleast_1d('[:,:]'),
#            oscansec        = np.atleast_1d('[:,:]')
            )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        par['flexure']['spec_method'] = 'boxcar'

        # 1D wavelength solution
        par['calibrations']['wavelengths']['sigdetect'] = 10.
        par['calibrations']['wavelengths']['rms_threshold'] = 0.4
        par['calibrations']['wavelengths']['lamps'] = ['XeI','KrI','ArI','HgI']

        # slit
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['edge_thresh'] = 100.

        # Set wave tilts order
        par['calibrations']['tilts']['spat_order'] = 5
        par['calibrations']['tilts']['spec_order'] = 5
        par['calibrations']['tilts']['maxdev_tracefit'] = 0.02
        par['calibrations']['tilts']['maxdev2d'] = 0.02

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
        # Start with instrument wide
        par = super(LBTMODS1BSpectrograph, self).config_specific_par(scifile, inp_par=inp_par)

        if self.get_meta_value(scifile, 'dispname') == 'G400L':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'lbt_mods1b_blue.fits'

        return par

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask.

        Even though they are both optional, either the precise shape for
        the image (``shape``) or an example file that can be read to get
        the shape (``filename`` using :func:`get_image_shape`) *must* be
        provided.

        Args:
            filename (:obj:`str` or None):
                An example file to use to get the image shape.
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.
            shape (tuple, optional):
                Processed image shape
                Required if filename is None
                Ignored if filename is not None
            msbias (`numpy.ndarray`_, optional):
                Processed bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)
        msgs.info("Using hard-coded BPM for  MODS1B")

        # Get the binning
        hdu = io.fits_open(filename)
        header = hdu[0].header
        xbin, ybin = header['CCDXBIN'], header['CCDYBIN']
        hdu.close()

        # Apply the mask
        bpm_img[6390//xbin:6392//xbin, 1262//ybin:1545//ybin] = 1

        bpm_img[3064//xbin, 1437//ybin:1937//ybin] = 1
        bpm_img[6490//xbin, 1161//ybin:1545//ybin] = 1
        bpm_img[7306//xbin, 783//ybin:1531//ybin] = 1

        return bpm_img


class LBTMODS2RSpectrograph(LBTMODSSpectrograph):
    """
    Child to handle LBT/MODS1R specific code
    """
    name = 'lbt_mods2r'
    camera = 'MODS2R'
    header_name = 'MODS2R'
    supported = True
    comment = 'MODS-II red spectrometer'

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        # Binning
        binning = '1,1' if hdu is None \
                    else f"{hdu[0].header['CCDXBIN']},{hdu[0].header['CCDYBIN']}"

        # Detector 1
        detector_dict = dict(
            binning= binning,
            det=1,
            dataext         = 0,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.123,
            darkcurr        = 0.4,
            saturation      = 65535.,
            nonlinear       = 0.99,
            mincounts       = -1e10,
            numamplifiers   = 4,
            gain            = np.atleast_1d([1.70,1.67,1.66,1.66]),
            ronoise         = np.atleast_1d([2.95,2.65,2.78,2.87]),
# TODO: The raw image reader sets these up by hand
#            datasec         = np.atleast_1d('[:,:]'),
#            oscansec        = np.atleast_1d('[:,:]')
            )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        par['flexure']['spec_method'] = 'boxcar'

        # 1D wavelength solution
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['rms_threshold'] = 1.0
        par['calibrations']['wavelengths']['fwhm'] = 10.
        #par['calibrations']['wavelengths']['lamps'] = ['XeI','ArII','ArI','NeI','KrI']]
        par['calibrations']['wavelengths']['lamps'] = ['ArI','NeI','KrI','XeI']
        #par['calibrations']['wavelengths']['lamps'] = ['OH_MODS']
        par['calibrations']['wavelengths']['n_first'] = 3
        par['calibrations']['wavelengths']['match_toler'] = 2.5

        # slit
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['edge_thresh'] = 300.

        # Set wave tilts order
        par['calibrations']['tilts']['spat_order'] = 5
        par['calibrations']['tilts']['spec_order'] = 5
        par['calibrations']['tilts']['maxdev_tracefit'] = 0.02
        par['calibrations']['tilts']['maxdev2d'] = 0.02

        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            scifile (:obj:`str`):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        par = super().config_specific_par(scifile, inp_par=inp_par)
        if self.get_meta_value(scifile, 'dispname') == 'G670L':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'lbt_mods2r_red.fits'
        return par

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask.

        Even though they are both optional, either the precise shape for
        the image (``shape``) or an example file that can be read to get
        the shape (``filename`` using :func:`get_image_shape`) *must* be
        provided.

        Args:
            filename (:obj:`str` or None):
                An example file to use to get the image shape.
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.
            shape (tuple, optional):
                Processed image shape
                Required if filename is None
                Ignored if filename is not None
            msbias (`numpy.ndarray`_, optional):
                Processed bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)
        msgs.info("Using hard-coded BPM for  MODS2R")

        # Get the binning
        hdu = io.fits_open(filename)
        header = hdu[0].header
        xbin, ybin = header['CCDXBIN'], header['CCDYBIN']
        hdu.close()

        # Apply the mask
        bpm_img[6148//xbin:6150//xbin, 1333//ybin:1544//ybin] = 1
        bpm_img[6207//xbin:6209//xbin, 1396//ybin:1544//ybin] = 1

        bpm_img[6101//xbin, 1342//ybin:1544//ybin] = 1
        bpm_img[6159//xbin, 1399//ybin:1544//ybin] = 1
        bpm_img[6189//xbin, 1316//ybin:1544//ybin] = 1
        bpm_img[7552//xbin, 1544//ybin:2771//ybin] = 1
        bpm_img[7504//xbin, 1544//ybin:2774//ybin] = 1
        bpm_img[4203//xbin, 0:1544//ybin] = 1
        bpm_img[4155//xbin, 0:1544//ybin] = 1

        return bpm_img



class LBTMODS2BSpectrograph(LBTMODSSpectrograph):
    """
    Child to handle LBT/MODS1R specific code
    """
    name = 'lbt_mods2b'
    camera = 'MODS2B'
    header_name = 'MODS2B'
    supported = True
    comment = 'MODS-II blue spectrometer'

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        # Binning
        binning = '1,1' if hdu is None \
                    else f"{hdu[0].header['CCDXBIN']},{hdu[0].header['CCDYBIN']}"

        # Detector 1
        detector_dict = dict(
            binning= binning,
            det=1,
            dataext         = 0,
            specaxis        = 0,
            specflip        = True,
            spatflip        = False,
            platescale      = 0.120,
            darkcurr        = 0.5,
            saturation      = 65535.,
            nonlinear       = 0.99,
            mincounts       = -1e10,
            numamplifiers   = 4,
            gain            = np.atleast_1d([1.99,2.06,1.96,2.01]),
            ronoise         = np.atleast_1d([3.66,3.62,3.72,3.64]),
# TODO: The raw image reader sets these up by hand
#            datasec         = np.atleast_1d('[:,:]'),
#            oscansec        = np.atleast_1d('[:,:]')
            )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        par['flexure']['spec_method'] = 'boxcar'

        # 1D wavelength solution
        par['calibrations']['wavelengths']['sigdetect'] = 10.
        par['calibrations']['wavelengths']['rms_threshold'] = 0.4
        par['calibrations']['wavelengths']['lamps'] = ['XeI','KrI','ArI','HgI']

        # slit
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['edge_thresh'] = 100.

        # Set wave tilts order
        par['calibrations']['tilts']['spat_order'] = 5
        par['calibrations']['tilts']['spec_order'] = 5
        par['calibrations']['tilts']['maxdev_tracefit'] = 0.02
        par['calibrations']['tilts']['maxdev2d'] = 0.02

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
        # Start with instrument wide
        par = super(LBTMODS2BSpectrograph, self).config_specific_par(scifile, inp_par=inp_par)

        if self.get_meta_value(scifile, 'dispname') == 'G400L':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'lbt_mods1b_blue.fits'

        return par


    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask.

        Even though they are both optional, either the precise shape for
        the image (``shape``) or an example file that can be read to get
        the shape (``filename`` using :func:`get_image_shape`) *must* be
        provided.

        Args:
            filename (:obj:`str` or None):
                An example file to use to get the image shape.
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.
            shape (tuple, optional):
                Processed image shape
                Required if filename is None
                Ignored if filename is not None
            msbias (`numpy.ndarray`_, optional):
                Processed bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)
        msgs.info("Using hard-coded BPM for  MODS2B")

        # Get the binning
        hdu = io.fits_open(filename)
        header = hdu[0].header
        xbin, ybin = header['CCDXBIN'], header['CCDYBIN']
        hdu.close()

        # Apply the mask
        bpm_img[5176//xbin:5179//xbin, 549//ybin:1544//ybin] = 1
        bpm_img[5176//xbin:5179//xbin, 1544//ybin:1628//ybin] = 1
        bpm_img[4408//xbin:4410//xbin, 1544//ybin:2661//ybin] = 1
        bpm_img[4408//xbin:4411//xbin, 2660//ybin:2663//ybin] = 1
        bpm_img[2495//xbin:2499//xbin, 1326//ybin:1544//ybin] = 1
        bpm_img[2391//xbin:2394//xbin, 1048//ybin:1051//ybin] = 1
        bpm_img[1974//xbin:1980//xbin, 806//ybin:1544//ybin] = 1
        bpm_img[1975//xbin:1980//xbin, 1544//ybin:1607//ybin] = 1
        bpm_img[1972//xbin:1974//xbin, 1587//ybin:1589//ybin] = 1
        bpm_img[274//xbin:278//xbin, 1341//ybin:1544//ybin] = 1
        bpm_img[275//xbin:278//xbin, 1251//ybin:1341//ybin] = 1
        bpm_img[276//xbin:278//xbin, 1242//ybin:1251//ybin] = 1
        bpm_img[274//xbin:277//xbin, 1544//ybin:3066//ybin] = 1

        bpm_img[2392//xbin, 1051//ybin:1544//ybin] = 1
        bpm_img[275//xbin, 1220//ybin:1242//ybin] = 1

        return bpm_img

