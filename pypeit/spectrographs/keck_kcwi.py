""" Implements KCWI-specific functions.
"""

import glob
import re
import os
import numpy as np
import warnings

from scipy import interpolate

from astropy.io import fits

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph

from pypeit.utils import index_of_x_eq_y

from pypeit.spectrographs.slitmask import SlitMask


class KeckKCWISpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/KCWI specific code
    """

    def __init__(self):
        # Get it started
        super(KeckKCWISpectrograph, self).__init__()
        self.spectrograph = 'keck_kcwi_base'
        self.telescope = telescopes.KeckTelescopePar()

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Keck KCWI reductions.
        """
        par = pypeitpar.PypeItPar()
        # Set wave tilts order
        #par['calibrations']['slitedges']['edge_thresh'] = 15.
        #par['calibrations']['slitedges']['fit_order'] = 3
        #par['calibrations']['slitedges']['sync_center'] = 'gap'
        #par['calibrations']['slitedges']['minimum_slit_length'] = 6
        # 1D wavelengths
        #par['calibrations']['wavelengths']['rms_threshold'] = 0.20  # Might be grism dependent
        # Always sky subtract, starting with default parameters
        #par['scienceimage'] = pypeitpar.ScienceImagePar()

        # Always flux calibrate, starting with default parameters
        #par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Always correct for flexure, starting with default parameters
        #par['flexure']['method'] = 'boxcar'

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [1, None]
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['pixelflatframe']['exprng'] = [None, 30]  # This may be too low for KCWI
        par['calibrations']['traceframe']['exprng'] = [None, 30]
        par['scienceframe']['exprng'] = [29, None]
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
        par = self.default_pypeit_par() if inp_par is None else inp_par

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
        meta['decker'] = dict(ext=0, card='IFUNAM')
        meta['binning'] = dict(card=None, compound=True)

        meta['mjd'] = dict(ext=0, card='MJD')
        meta['exptime'] = dict(ext=0, card='ELAPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        meta['hatch'] = dict(ext=0, card='HATNUM')
        meta['dispname'] = dict(ext=0, card='BGRATNAM')
        meta['dispangle'] = dict(ext=0, card='BGRANGLE', rtol=1e-2)

        # Lamps
        lamp_names = ['LMP0', 'LMP1', 'LMP2', 'LMP3']  # FeAr, ThAr, Aux, Continuum
        for kk, lamp_name in enumerate(lamp_names):
            meta['lampstat{:02d}'.format(kk + 1)] = dict(ext=0, card=lamp_name+'STAT')
        for kk, lamp_name in enumerate(lamp_names):
            if lamp_name == 'LMP3':
                # There is no shutter on LMP3
                continue
            meta['lampshst{:02d}'.format(kk + 1)] = dict(ext=0, card=lamp_name+'SHST')
        # Ingest
        self.meta = meta

    def compound_meta(self, headarr, meta_key):
        if meta_key == 'binning':
            binspatial, binspec = parse.parse_binning(headarr[0]['BINNING'])
            binning = parse.binning2string(binspec, binspatial)
            return binning
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
            list: List of keywords of data pulled from file headers and
            used to constuct the :class:`pypeit.metadata.PypeItMetaData`
            object.
        """
        return ['dispname', 'decker', 'binning']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype == 'science':
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == 1)  #hatch=1,0=open,closed
        if ftype == 'bias':
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == 0)
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & self.lamps(fitstbl, 'dome') & (fitstbl['hatch'] == 1)
        if ftype in ['dark']:
            # dark
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == 0)
        if ftype in ['pinhole']:
            # Don't type pinhole frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & self.lamps(fitstbl, 'arcs') & (fitstbl['hatch'] == 0)

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def lamps(self, fitstbl, status):
        """
        Check the lamp status.

        Args:
            fitstbl (:obj:`astropy.table.Table`):
                The table with the fits header meta data.
            status (:obj:`str`):
                The status to check.  Can be `off`, `arcs`, or `dome`.

        Returns:
            numpy.ndarray: A boolean array selecting fits files that
            meet the selected lamp status.

        Raises:
            ValueError:
                Raised if the status is not one of the valid options.
        """
        if status == 'off':
            # Check if all are off
            lampstat = np.array([(fitstbl[k] == 0) | (fitstbl[k] == 'None')
                                    for k in fitstbl.keys() if 'lampstat' in k])
            lampshst = np.array([(fitstbl[k] == 0) | (fitstbl[k] == 'None')
                                    for k in fitstbl.keys() if 'lampshst' in k])
            return np.all(lampstat|lampshst, axis=0)  # i.e. either the shutter is closed or the lamp is off
        if status == 'arcs':
            # Check if any arc lamps are on (FeAr | ThAr)
            arc_lamp_stat = ['lampstat{0:02d}'.format(i) for i in range(1, 3)]
            arc_lamp_shst = ['lampshst{0:02d}'.format(i) for i in range(1, 3)]
            lamp_stat = np.array([fitstbl[k] == 1 for k in fitstbl.keys()
                                    if k in arc_lamp_stat])
            lamp_shst = np.array([fitstbl[k] == 1 for k in fitstbl.keys()
                                    if k in arc_lamp_shst])
            return np.any(lamp_stat&lamp_shst, axis=0)  # i.e. lamp on and shutter open
        if status == 'dome':
            # Check if any dome lamps are on (Continuum) - Ignore lampstat03 (Aux) - not sure what this is used for
            dome_lamp_stat = ['lampstat{0:02d}'.format(i) for i in range(4, 5)]
            lamp_stat = np.array([fitstbl[k] == 1 for k in fitstbl.keys()
                                    if k in dome_lamp_stat])
            return np.any(lamp_stat, axis=0)  # i.e. lamp on
        raise ValueError('No implementation for status = {0}'.format(status))

    def get_rawimage(self, raw_file, det):
        """
        Read a raw KCWI data frame

        Parameters
        ----------
        raw_file : str
          Filename
        det (int or None):
          Detector number
        Returns
        -------
        array : ndarray
          Combined image
        hdu : HDUList
        sections : list
          List of datasec, oscansec, ampsec sections
          datasec, oscansec needs to be for an *unbinned* image as per standard convention
        """
        # Check for file; allow for extra .gz, etc. suffix
        fil = glob.glob(raw_file + '*')
        if len(fil) != 1:
            msgs.error("Found {:d} files matching {:s}".format(len(fil), raw_file))

        # Read
        msgs.info("Reading KCWI file: {:s}".format(fil[0]))
        hdu = fits.open(fil[0])
        head0 = hdu[0].header
        raw_img = hdu[self.detector[det-1]['dataext']].data.astype(float)

        # Some properties of the image
        numamps = head0['NVIDINP']
        # Exposure time (used by ProcessRawImage)
        headarr = self.get_headarr(hdu)
        exptime = self.get_meta_value(headarr, 'exptime')

        # get the x and y binning factors...
        binning = head0['BINNING']
        xbin, ybin = [int(ibin) for ibin in binning.split(',')]
        binning_raw = binning

        # Always assume normal FITS header formatting
        one_indexed = True
        include_last = True
        for section in ['DSEC', 'BSEC']:

            # Initialize the image (0 means no amplifier)
            pix_img = np.zeros(raw_img.shape, dtype=int)
            for i in range(numamps):
                # Get the data section
                sec = head0[section+"{0:1d}".format(i+1)]

                # Convert the data section from a string to a slice
                datasec = parse.sec2slice(sec, one_indexed=one_indexed,
                                          include_end=include_last, require_dim=2,
                                          binning=binning_raw)
                # Assign the amplifier
                pix_img[datasec] = i+1
            # Finish
            if section == 'DSEC':
                rawdatasec_img = pix_img.copy()
            elif section == 'BSEC':
                oscansec_img = pix_img.copy()

        # Return
        return raw_img, hdu, exptime, rawdatasec_img, oscansec_img


class KeckKCWIBSpectrograph(KeckKCWISpectrograph):
    """
    Child to handle Keck/KCWI specific code
    """
    def __init__(self):
        # Get it started
        super(KeckKCWISpectrograph, self).__init__()
        self.spectrograph = 'keck_kcwi_blue'
        self.telescope = telescopes.KeckTelescopePar()
        self.camera = 'KCWIb'
        self.detector = [pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = ,
                            darkcurr        = ,
                            saturation      = 65535.,
                            nonlinear       = 0.95,   # For lack of a better number!
                            numamplifiers   = 0,
                            gain            = 0,
                            ronoise         = 0,
                            datasec         = '',       # These are provided by get_rawimage
                            oscansec        = '',
                            suffix          = '_01'
                            )]
        self.numhead = 1
        # Uses default timeunit
        # Uses default primary_hdrext
        # self.sky_file ?

        # Don't instantiate these until they're needed
        self.grating = None
        self.optical_model = None
        self.detector_map = None

    def default_pypeit_par(self):
        """
        Set default parameters for Keck DEIMOS reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'keck_kcwi_blue'
        par['flexure']['method'] = 'boxcar'
        # Set wave tilts order
        par['calibrations']['slitedges']['edge_thresh'] = 50.
        par['calibrations']['slitedges']['fit_order'] = 3
        par['calibrations']['slitedges']['minimum_slit_gap'] = 0.25
        par['calibrations']['slitedges']['minimum_slit_length'] = 4.
        par['calibrations']['slitedges']['sync_clip'] = False

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
        par = self.default_pypeit_par() if inp_par is None else inp_par

        headarr = self.get_headarr(scifile)

        # Templates
        if self.get_meta_value(headarr, 'dispname') == '600ZD':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_deimos_600.fits'
            par['calibrations']['wavelengths']['lamps'] += ['CdI', 'ZnI', 'HgI']
        elif self.get_meta_value(headarr, 'dispname') == '830G':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_deimos_830G.fits'
        elif self.get_meta_value(headarr, 'dispname') == '1200G':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_deimos_1200G.fits'

        # FWHM
        binning = parse.parse_binning(self.get_meta_value(headarr, 'binning'))
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
        # Image type
        meta['idname'] = dict(ext=0, card='OBSTYPE')
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
                msgs.warn('This is probably a problem. Non-standard DEIMOS GRATEPOS={0}.'.format(headarr[0]['GRATEPOS']))
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
            #return good_exp & (fitstbl['lampstat01'] == 'Off') & (fitstbl['hatch'] == 'open')
            return good_exp & (fitstbl['lampstat01'] == 'Off') & (fitstbl['hatch'] == 'open')
        if ftype == 'bias':
            return good_exp & (fitstbl['lampstat01'] == 'Off') & (fitstbl['hatch'] == 'closed')
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & (fitstbl['idname'] == 'IntFlat') & (fitstbl['hatch'] == 'closed')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'Line') & (fitstbl['hatch'] == 'closed')

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
                 'tilt': None,
                 'bias': None,
                 'dark': None,
                 'pinhole': None,
                 'pixelflat': 'IntFlat',
                 'science': 'Object',
                 'standard': None,
                 'trace': 'IntFlat' }
        return name[ftype]

    def get_rawimage(self, raw_file, det):
        """
        Read a raw KCWI data frame (one or more detectors).

        Data are unpacked from the multi-extension HDU.  Function is
        based :func:`pypeit.spectrographs.keck_lris.read_lris`, which
        was based on the IDL procedure ``readmhdufits.pro``.

        Parameters
        ----------
        raw_file : str
            Filename

        Returns
        -------
        array : ndarray
            Combined image
        hdu: HDUList
        sections : tuple
            List of datasec, oscansec sections

        """
        # Check for file; allow for extra .gz, etc. suffix
        fil = glob.glob(raw_file + '*')
        if len(fil) != 1:
            msgs.error('Found {0} files matching {1}'.format(len(fil), raw_file + '*'))
        # Read
        try:
            msgs.info("Reading KCWI file: {:s}".format(fil[0]))
        except AttributeError:
            print("Reading KCWI file: {:s}".format(fil[0]))

        hdu = fits.open(fil[0])
        head0 = hdu[0].header

        # Get post, pre-pix values
        postpix = head0['POSTPIX']
        detlsize = head0['DETLSIZE']
        x0, x_npix, y0, y_npix = np.array(parse.load_sections(detlsize)).flatten()

        # Create final image
        if det is None:
            image = np.zeros((x_npix, y_npix + 4 * postpix))
            rawdatasec_img = np.zeros_like(image, dtype=int)
            oscansec_img = np.zeros_like(image, dtype=int)

        # get the x and y binning factors...
        binning = head0['BINNING']

        # DEIMOS detectors
        nchip = 8

        if det is None:
            chips = range(nchip)
        else:
            chips = [det - 1]  # Indexing starts at 0 here
        # Loop
        for tt in chips:
            data, oscan = deimos_read_1chip(hdu, tt + 1)

            # One detector??
            if det is not None:
                image = np.zeros((data.shape[0], data.shape[1] + oscan.shape[1]))
                rawdatasec_img = np.zeros_like(image, dtype=int)
                oscansec_img = np.zeros_like(image, dtype=int)

            # Indexing
            x1, x2, y1, y2, o_x1, o_x2, o_y1, o_y2 = indexing(tt, postpix, det=det)

            # Fill
            image[y1:y2, x1:x2] = data
            rawdatasec_img[y1:y2, x1:x2] = 1 # Amp
            image[o_y1:o_y2, o_x1:o_x2] = oscan
            oscansec_img[o_y1:o_y2, o_x1:o_x2] = 1 # Amp

        # Return
        exptime = hdu[self.meta['exptime']['ext']].header[self.meta['exptime']['card']]
        return image, hdu, exptime, rawdatasec_img, oscansec_img
        #return image, hdu, (dsec, osec)

    '''
    def load_raw_frame(self, raw_file, det=None):
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
        raw_img, hdu, _ = read_deimos(raw_file, det=det)

        return raw_img, hdu
    '''

    '''
    def get_image_section(self, inp=None, det=1, section='datasec'):
        """
        Return a string representation of a slice defining a section of
        the detector image.

        Overwrites base class function

        Args:
            inp (:obj:`str`, `astropy.io.fits.Header`_, optional):
                String providing the file name to read, or the relevant
                header object.  Default is None, meaning that the
                detector attribute must provide the image section
                itself, not the header keyword.
            det (:obj:`int`, optional):
                1-indexed detector number.
            section (:obj:`str`, optional):
                The section to return.  Should be either 'datasec' or
                'oscansec', according to the
                :class:`pypeitpar.DetectorPar` keywords.

        Returns:
            tuple: Returns three objects: (1) A list of string
            representations for the image sections, one string per
            amplifier.  The sections are *always* returned in PypeIt
            order: spectral then spatial.  (2) Boolean indicating if the
            slices are one indexed.  (3) Boolean indicating if the
            slices should include the last pixel.  The latter two are
            always returned as True following the FITS convention.
        """
        # Read the file
        if inp is None:
            msgs.error('Must provide Keck DEIMOS file or hdulist to get image section.')
        # Read em
        shape, datasec, oscansec, _ = deimos_image_sections(inp, det)
        if section == 'datasec':
            return datasec, False, False
        elif section == 'oscansec':
            return oscansec, False, False
        else:
            raise ValueError('Unrecognized keyword: {0}'.format(section))

    def get_raw_image_shape(self, hdulist, det=None, **null_kwargs):
        """
        Overrides :class:`Spectrograph.get_image_shape` for LRIS images.

        Must always provide a file.
        """
        # Do it
        self._check_detector()
        shape, datasec, oscansec, _ = deimos_image_sections(hdulist, det)
        self.naxis = shape
        return self.naxis
    '''

    def bpm(self, filename, det, shape=None):
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
        bpm_img = self.empty_bpm(filename, det, shape=shape)
        if det == 1:
            bpm_img[:,1052:1054] = 1
        elif det == 2:
            bpm_img[:,0:4] = 1
            bpm_img[:,376:381] = 1
            bpm_img[:,489] = 1
            bpm_img[:,1333:1335] = 1
            bpm_img[:,2047] = 1
        elif det == 3:
            bpm_img[:,0:4] = 1
            bpm_img[:,221] = 1
            bpm_img[:,260] = 1
            bpm_img[:,366] = 1
            bpm_img[:,816:819] = 1
            bpm_img[:,851] = 1
            bpm_img[:,940] = 1
            bpm_img[:,1167] = 1
            bpm_img[:,1280] = 1
            bpm_img[:,1301:1303] = 1
            bpm_img[:,1744:1747] = 1
            bpm_img[:,-4:] = 1
        elif det == 4:
            bpm_img[:,0:4] = 1
            bpm_img[:,47] = 1
            bpm_img[:,744] = 1
            bpm_img[:,790:792] = 1
            bpm_img[:,997:999] = 1
        elif det == 5:
            bpm_img[:,25:27] = 1
            bpm_img[:,128:130] = 1
            bpm_img[:,1535:1539] = 1
        elif det == 7:
            bpm_img[:,426:428] = 1
            bpm_img[:,676] = 1
            bpm_img[:,1176:1178] = 1
        elif det == 8:
            bpm_img[:,440] = 1
            bpm_img[:,509:513] = 1
            bpm_img[:,806] = 1
            bpm_img[:,931:934] = 1

        return bpm_img

    def get_slitmask(self, filename):
        """
        Parse the slitmask data from a DEIMOS file into a
        :class:`pypeit.spectrographs.slitmask.SlitMask` object.

        Args:
            filename (:obj:`str`):
                Name of the file to read.
        """
        # Open the file
        hdu = fits.open(filename)

        # Build the object data
        #   - Find the index of the object IDs in the slit-object
        #     mapping that match the object catalog
        mapid = hdu['SlitObjMap'].data['ObjectID']
        catid = hdu['ObjectCat'].data['ObjectID']
        indx = index_of_x_eq_y(mapid, catid)
        #   - Pull out the slit ID, object ID, and object coordinates
        objects = np.array([hdu['SlitObjMap'].data['dSlitId'][indx].astype(float),
                            catid.astype(float), hdu['ObjectCat'].data['RA_OBJ'],
                            hdu['ObjectCat'].data['DEC_OBJ']]).T
        #   - Only keep the objects that are in the slit-object mapping
        objects = objects[mapid[indx] == catid]

        # Match the slit IDs in DesiSlits to those in BluSlits
        indx = index_of_x_eq_y(hdu['DesiSlits'].data['dSlitId'], hdu['BluSlits'].data['dSlitId'],
                               strict=True)

        # Instantiate the slit mask object and return it
        self.slitmask = SlitMask(np.array([hdu['BluSlits'].data['slitX1'],
                                           hdu['BluSlits'].data['slitY1'],
                                           hdu['BluSlits'].data['slitX2'],
                                           hdu['BluSlits'].data['slitY2'],
                                           hdu['BluSlits'].data['slitX3'],
                                           hdu['BluSlits'].data['slitY3'],
                                           hdu['BluSlits'].data['slitX4'],
                                           hdu['BluSlits'].data['slitY4']]).T.reshape(-1,4,2),
                                 slitid=hdu['BluSlits'].data['dSlitId'],
                                 align=hdu['DesiSlits'].data['slitTyp'][indx] == 'A',
                                 science=hdu['DesiSlits'].data['slitTyp'][indx] == 'P',
                                 onsky=np.array([hdu['DesiSlits'].data['slitRA'][indx],
                                                 hdu['DesiSlits'].data['slitDec'][indx],
                                                 hdu['DesiSlits'].data['slitLen'][indx],
                                                 hdu['DesiSlits'].data['slitWid'][indx],
                                                 hdu['DesiSlits'].data['slitLPA'][indx]]).T,
                                 objects=objects)
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
