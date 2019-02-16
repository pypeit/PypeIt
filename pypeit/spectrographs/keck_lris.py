""" Module for LRIS specific codes
"""
from __future__ import absolute_import, division, print_function

import glob
import os
import numpy as np
from astropy.io import fits

from pkg_resources import resource_filename

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph

from pypeit import debugger

class KeckLRISSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/LRIS specific code
    """
    def __init__(self):
        # Get it started
        super(KeckLRISSpectrograph, self).__init__()
        self.spectrograph = 'keck_lris_base'
        self.telescope = telescopes.KeckTelescopePar()

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Keck LRISr reductions.
        """
        par = pypeitpar.PypeItPar()
        # Set wave tilts order
        par['calibrations']['slits']['sigdetect'] = 30.
        # 1D wavelengths
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20  # Might be grism dependent
        # Always sky subtract, starting with default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()

        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Always correct for flexure, starting with default parameters
        par['flexure']['method'] = 'boxcar'

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['pixelflatframe']['exprng'] = [None, 30]    # This may be too low for LRISb
        par['calibrations']['traceframe']['exprng'] = [None, 30]
        par['scienceframe']['exprng'] = [29, None]
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
        meta['decker'] = dict(ext=0, card='SLITNAME')
        meta['binning'] = dict(card=None, compound=True)

        meta['mjd'] = dict(ext=0, card='MJD-OBS')
        meta['exptime'] = dict(ext=0, card='ELAPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        meta['dichroic'] = dict(ext=0, card='DICHNAME')
        meta['hatch'] = dict(ext=0, card='TRAPDOOR')
        # Red only, but grabbing here
        meta['dispangle'] = dict(ext=0, card='GRANGLE', rtol=1e-2)

        # Lamps
        lamp_names = ['MERCURY', 'NEON', 'ARGON', 'CADMIUM', 'ZINC', 'KRYPTON', 'XENON',
                      'FEARGON', 'DEUTERI', 'FLAMP1', 'FLAMP2', 'HALOGEN']
        for kk,lamp_name in enumerate(lamp_names):
            meta['lampstat{:02d}'.format(kk+1)] = dict(ext=0, card=lamp_name)
        # Ingest
        self.meta = meta

    def compound_meta(self, headarr, meta_key):
        if meta_key == 'binning':
#            return '1,1'
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
        return ['dispname', 'dichroic', 'decker', 'binning']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype == 'science':
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == 'open')
        if ftype == 'bias':
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == 'closed')
        if ftype == 'pixelflat' or ftype == 'trace':
            # Flats and trace frames are typed together
            return good_exp & self.lamps(fitstbl, 'dome') & (fitstbl['hatch'] == 'open')
        if ftype == 'pinhole' or ftype == 'dark':
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype == 'arc':
            return good_exp & self.lamps(fitstbl, 'arcs') & (fitstbl['hatch'] == 'closed')

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
            return np.all(np.array([ (fitstbl[k] == 'off') | (fitstbl[k] == 'None')
                                        for k in fitstbl.keys() if 'lampstat' in k]), axis=0)
        if status == 'arcs':
            # Check if any arc lamps are on
            arc_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(1,9) ]
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys()
                                            if k in arc_lamp_stat]), axis=0)
        if status == 'dome':
            # Check if any dome lamps are on
            # Warning 9, 10 are FEARGON and DEUTERI
            dome_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(9,13) ]
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys()
                                            if k in dome_lamp_stat]), axis=0)
        raise ValueError('No implementation for status = {0}'.format(status))
        
    def load_raw_img_head(self, raw_file, det=None, **null_kwargs):
        """
        Wrapper to the raw image reader for LRIS

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
        raw_img, head0, _ = read_lris(raw_file, det=det)

        return raw_img, head0

    def get_image_section(self, inp=None, det=1, section='datasec'):
        """
        Return a string representation of a slice defining a section of
        the detector image.

        Overwrites base class function to use :func:`read_lris` to get
        the image sections.

        .. todo ::
            - It is really ineffiecient.  Can we parse
              :func:`read_deimos` into something that can give you the
              image section directly?

        This is done separately for the data section and the overscan
        section in case one is defined as a header keyword and the other
        is defined directly.
        
        Args:
            inp (:obj:`str`):
                String providing the file name to read.  Unlike the base
                class, a file name *must* be provided.
            det (:obj:`int`, optional):
                1-indexed detector number.
            section (:obj:`str`, optional):
                The section to return.  Should be either datasec or
                oscansec, according to the :class:`DetectorPar`
                keywords.

        Returns:
            list, bool: A list of string representations for the image
            sections, one string per amplifier, followed by three
            booleans: if the slices are one indexed, if the slices
            should include the last pixel, and if the slice should have
            their order transposed.
        """
        # Read the file
        if inp is None:
            msgs.error('Must provide Keck LRIS file to get image section.')
        elif not os.path.isfile(inp):
            msgs.error('File {0} does not exist!'.format(inp))
        temp, head0, secs = read_lris(inp, det)
        if section == 'datasec':
            return secs[0], False, False, False
        elif section == 'oscansec':
            return secs[1], False, False, False
        else:
            raise ValueError('Unrecognized keyword: {0}'.format(section))

    '''
    def get_datasec_img(self, filename, det=1, force=True):
        """
        Create an image identifying the amplifier used to read each pixel.

        Args:
            filename (str):
                Name of the file from which to read the image size.
            det (:obj:`int`, optional):
                Detector number (1-indexed)
            force (:obj:`bool`, optional):
                Force the image to be remade

        Returns:
            `numpy.ndarray`: Integer array identifying the amplifier
            used to read each pixel.
        """
        if self.datasec_img is None or force:
            # Check the detector is defined
            self._check_detector()
            # Get the image shape
            raw_naxis = self.get_raw_image_shape(filename, det=det)

            # Binning is not required because read_lris accounts for it
#            binning = self.get_meta_value(filename, 'binning')

            data_sections, one_indexed, include_end, transpose \
                    = self.get_image_section(filename, det, section='datasec')

            # Initialize the image (0 means no amplifier)
            self.datasec_img = np.zeros(raw_naxis, dtype=int)
            for i in range(self.detector[det-1]['numamplifiers']):
                # Convert the data section from a string to a slice
                datasec = parse.sec2slice(data_sections[i], one_indexed=one_indexed,
                                          include_end=include_end, require_dim=2,
                                          transpose=transpose) #, binning=binning)
                # Assign the amplifier
                self.datasec_img[datasec] = i+1
        return self.datasec_img
    '''

    def get_image_shape(self, filename=None, det=None, **null_kwargs):
        """
        Overrides :class:`Spectrograph.get_image_shape` for LRIS images.

        Must always provide a file.
        """
        # Cannot be determined without file
        if filename is None:
            raise ValueError('Must provide a file to determine the shape of an LRIS image.')

        # Use a file
        self._check_detector()
        self.naxis = (self.load_raw_frame(filename, det=det)[0]).shape
        return self.naxis

    def get_match_criteria(self):
        match_criteria = {}
        for key in framematch.FrameTypeBitMask().keys():
            match_criteria[key] = {}
        #
        match_criteria['standard']['match'] = {}
        match_criteria['standard']['match']['dispname'] = ''
        match_criteria['standard']['match']['dichroic'] = ''
        match_criteria['standard']['match']['binning'] = ''
        match_criteria['standard']['match']['decker'] = ''
        # Bias
        match_criteria['bias']['match'] = {}
        match_criteria['bias']['match']['binning'] = ''
        # Pixelflat
        match_criteria['pixelflat']['match'] = match_criteria['standard']['match'].copy()
        # Traceflat
        match_criteria['trace']['match'] = match_criteria['standard']['match'].copy()
        # Arc
        match_criteria['arc']['match'] = match_criteria['standard']['match'].copy()

        # Return
        return match_criteria


class KeckLRISBSpectrograph(KeckLRISSpectrograph):
    """
    Child to handle Keck/LRISb specific code
    """
    def __init__(self):
        # Get it started
        super(KeckLRISBSpectrograph, self).__init__()
        self.spectrograph = 'keck_lris_blue'
        self.camera = 'LRISb'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 1,
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.135,
                            darkcurr        = 0.0,
                            saturation      = 65535.,
                            nonlinear       = 0.86,
                            numamplifiers   = 2,
                            gain            = [1.55, 1.56],
                            ronoise         = [3.9, 4.2],
                            datasec         = ['',''],      # These are provided by read_lris
                            oscansec        = ['',''],
                            suffix          = '_01blue'
                            ),
                #Detector 2
                pypeitpar.DetectorPar(
                            dataext         = 2,
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.135,
                            darkcurr        = 0.,
                            saturation      = 65535.,
                            nonlinear       = 0.86,
                            numamplifiers   = 2,
                            gain            = [1.63, 1.70],
                            ronoise         = [3.6, 3.6],
                            datasec         = ['',''],      # These are provided by read_lris
                            oscansec        = ['',''],
                            suffix          = '_02blue'
                            )]
        self.numhead = 5
        # Uses default timeunit
        # Uses default primary_hdrext
        self.sky_file = 'sky_LRISb_600.fits'

    def default_pypeit_par(self):
        """
        Set default parameters for Keck LRISr reductions.
        """
        par = KeckLRISSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'keck_lris_blue'
        # 1D wavelength solution -- Additional parameters are grism dependent
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20  # Might be grism dependent..
        par['calibrations']['wavelengths']['sigdetect'] = 10.0

        par['calibrations']['wavelengths']['lamps'] = ['NeI', 'ArI', 'CdI', 'KrI', 'XeI', 'ZnI', 'HgI']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 3
        par['calibrations']['wavelengths']['match_toler'] = 2.5
        par['calibrations']['wavelengths']['method'] = 'full_template'


        return par

    def config_specific_par(self, par, scifile):
        """
        Set par values according to the specific frame

        Args:
            par:  ParSet
            scifile: str
              Name of the science file to use

        Returns:
            par

        """
        # Wavelength calibrations
        if self.get_meta_value(scifile, 'dispname') == '300/5000':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_blue_300_d680.fits'
            par['flexure']['spectrum'] = os.path.join(resource_filename('pypeit', 'data/sky_spec/'),
                                                      'sky_LRISb_400.fits')
        elif self.get_meta_value(scifile, 'dispname') == '400/3400':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_blue_400_d560.fits'
            par['flexure']['spectrum'] = os.path.join(resource_filename('pypeit', 'data/sky_spec/'),
                                                  'sky_LRISb_400.fits')
        elif self.get_meta_value(scifile, 'dispname') == '600/4000':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_blue_600_d560.fits'
            par['flexure']['spectrum'] = os.path.join(resource_filename('pypeit', 'data/sky_spec/'),
                                                      'sky_LRISb_600.fits')
        elif self.get_meta_value(scifile, 'dispname') == '1200/3400':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_blue_1200_d460.fits'
            par['flexure']['spectrum'] = os.path.join(resource_filename('pypeit', 'data/sky_spec/'),
                                                      'sky_LRISb_600.fits')

        # FWHM
        binning = parse.parse_binning(self.get_meta_value(scifile, 'binning'))
        par['calibrations']['wavelengths']['fwhm'] = 8.0 / binning[0]

        # Slit tracing
        # Reduce the slit parameters because the flux does not span the full detector
        #   It is primarily on the upper half of the detector (usually)
        if self.get_meta_value(scifile, 'dispname') == '300/5000':
            par['calibrations']['slits']['mask_frac_thresh'] = 0.45
            par['calibrations']['slits']['smash_range'] = [0.5, 1.]

        # Return
        return par

    '''
    def check_headers(self, headers):
        """
        Check headers match expectations for an LRISb exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = { '0.INSTRUME': 'LRISBLUE',
                               '1.NAXIS': 2,
                               '2.NAXIS': 2,
                               '3.NAXIS': 2,
                               '4.NAXIS': 2,
                             '1.CCDGEOM': 'e2v (Marconi) CCD44-82',
                             '1.CCDNAME': '00151-14-1' }
        super(KeckLRISBSpectrograph, self).check_headers(headers, expected_values=expected_values)
    '''

    '''
    def header_keys(self):
        hdr_keys = super(KeckLRISBSpectrograph, self).header_keys()
        hdr_keys[0]['filter1'] = 'BLUFILT'
        return hdr_keys
    '''

    def init_meta(self):
        """
        Meta data specific to Keck LRIS red

        Returns:

        """
        super(KeckLRISBSpectrograph, self).init_meta()
        # Add the name of the dispersing element
        self.meta['dispname'] = dict(ext=0, card='GRISNAME')

    def bpm(self, shape=None, filename=None, det=None, **null_kwargs):
        """ Generate a BPM

        Parameters
        ----------
        shape : tuple, REQUIRED
        filename : str,
        det : int, REQUIRED
        **null_kwargs:
           Captured and never used

        Returns
        -------
        badpix : ndarray

        """
        # Get the empty bpm: force is always True
        self.empty_bpm(shape=shape, filename=filename, det=det)

        # Only defined for det=1
        if det == 1:
            msgs.info("Using hard-coded BPM for det=1 on LRISb")
            self.bpm_img[:, 0:2] = 1

        return self.bpm_img


class KeckLRISRSpectrograph(KeckLRISSpectrograph):
    """
    Child to handle Keck/LRISr specific code
    """
    def __init__(self):
        # Get it started
        super(KeckLRISRSpectrograph, self).__init__()
        self.spectrograph = 'keck_lris_red'
        self.camera = 'LRISr'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         =1,
                            specaxis        =0,
                            specflip        = False,
                            xgap            =0.,
                            ygap            =0.,
                            ysize           =1.,
                            platescale      =0.135,
                            darkcurr        =0.0,
                            saturation      =65535.,
                            nonlinear       =0.76,
                            numamplifiers   =2,
                            gain            =[1.255, 1.18],
                            ronoise         =[4.64, 4.76],
                            datasec         = ['',''],      # These are provided by read_lris
                            oscansec        = ['',''],
                            suffix          ='_01red'
                            ),
                #Detector 2
                pypeitpar.DetectorPar(
                            dataext         =2,
                            specaxis        =0,
                            specflip        = False,
                            xgap            =0.,
                            ygap            =0.,
                            ysize           =1.,
                            platescale      =0.135,
                            darkcurr        =0.,
                            saturation      =65535., 
                            nonlinear       =0.76,
                            numamplifiers   =2,
                            gain            =[1.191, 1.162],
                            ronoise         =[4.54, 4.62],
                            datasec         = ['',''],      # These are provided by read_lris
                            oscansec        = ['',''],
                            suffix          ='_02red'
                            )]
        self.numhead = 5
        # Uses default timeunit

    def default_pypeit_par(self):
        """
        Set default parameters for Keck LRISr reductions.
        """
        par = KeckLRISSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'keck_lris_red'
        #
        par['calibrations']['slits']['sigdetect'] = 50.

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['NeI', 'ArI', 'CdI', 'KrI', 'XeI', 'ZnI', 'HgI']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['sigdetect'] = 10.0
        # Tilts
        # These are the defaults
        par['calibrations']['tilts']['tracethresh'] = 25
        par['calibrations']['tilts']['spat_order'] = 4
        par['calibrations']['tilts']['spec_order'] = 7
        par['calibrations']['tilts']['maxdev2d'] = 1.0
        par['calibrations']['tilts']['maxdev_tracefit'] = 1.0
        par['calibrations']['tilts']['sigrej2d'] = 5.0

        # Scienceimage
        par['scienceimage']['bspline_spacing'] = 0.8

        # Defaults for anything other than 1,1 binning
        #  Rest config_specific_par below if binning is (1,1)
        par['scienceframe']['process']['sigclip'] = 5.
        par['scienceframe']['process']['objlim'] = 5.

        # reidentification stuff
        #par['calibrations']['wavelengths']['method'] = 'reidentify'
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_red_400_8500_d560.json'
        return par

    def config_specific_par(self, par, scifile):
        """
        Set par values according to the specific frame

        Args:
            par:  ParSet
            scifile: str
              Name of the science file to use

        Returns:
            par

        """
        # Lacosmic CR settings
        #   Grab the defaults for LRISr
        binning = self.get_meta_value(scifile, 'binning')
        # Unbinned LRISr needs very aggressive LACosmics parameters for 1x1 binning
        if binning == '1,1':
            sigclip = 3.0
            objlim = 0.5
            par['scienceframe']['process']['sigclip'] = sigclip
            par['scienceframe']['process']['objlim'] = objlim

        # Wavelength calibrations
        if self.get_meta_value(scifile, 'dispname') == '400/8500':  # This is basically a reidentify
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_red_400.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['sigdetect'] = 20.0
            par['calibrations']['wavelengths']['nsnippet'] = 1
        elif self.get_meta_value(scifile, 'dispname') == '1200/9000':  # This is basically a reidentify
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_red_1200_9000.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'

        # FWHM
        binning = parse.parse_binning(self.get_meta_value(scifile, 'binning'))
        par['calibrations']['wavelengths']['fwhm'] = 8.0 / binning[0]


        # Return
        return par

    '''
    def check_headers(self, headers):
        """
        Check headers match expectations for an LRISr exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = { '0.INSTRUME': 'LRIS',
                               '1.NAXIS': 2,
                               '2.NAXIS': 2,
                               '3.NAXIS': 2,
                               '4.NAXIS': 2,
                             '1.CCDGEOM': 'LBNL Thick High-Resistivity',
                             '1.CCDNAME': '19-3',
                             '3.CCDNAME': '19-2' }
        super(KeckLRISRSpectrograph, self).check_headers(headers, expected_values=expected_values)
    '''

    '''
    def header_keys(self):
        hdr_keys = super(KeckLRISRSpectrograph, self).header_keys()
        hdr_keys[0]['filter1'] = 'REDFILT'
        return hdr_keys
    '''

    def init_meta(self):
        """
        Meta data specific to Keck LRIS red

        Returns:

        """
        super(KeckLRISRSpectrograph, self).init_meta()
        # Add the name of the dispersing element
        self.meta['dispname'] = dict(ext=0, card='GRANAME')

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
        cfg_keys = super(KeckLRISRSpectrograph, self).configuration_keys()
        # Add grating tilt
        return cfg_keys+['dispangle']

    def bpm(self, shape=None, filename=None, det=None, **null_kwargs):
        """ Generate a BPM

        Parameters
        ----------
        shape : tuple, REQUIRED
        filename : str, REQUIRED for binning
        det : int, REQUIRED
        **null_kwargs:
           Captured and never used

        Returns
        -------
        badpix : ndarray

        """
        # Get the empty bpm: force is always True
        self.empty_bpm(shape=shape, filename=filename, det=det)
        
        # Only defined for det=2
        if det == 2:
            msgs.info("Using hard-coded BPM for det=2 on LRISr")

            # Get the binning
            hdu = fits.open(filename)
            binning = hdu[0].header['BINNING']
            hdu.close()

            # Apply the mask
            xbin = int(binning.split(',')[0])
            badc = 16//xbin
            self.bpm_img[:, 0:badc] = 1

        return self.bpm_img


class KeckLRISRLSpectrograph(KeckLRISRSpectrograph):
    """
    Child to handle Keck/LRISr in Long-slit readout mode (Vid1, Vid4)
    """
    def __init__(self):
        # Get it started
        super(KeckLRISRSpectrograph, self).__init__()
        self.spectrograph = 'keck_lris_red_longonly'
        self.camera = 'LRISr'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         =1,
                            specaxis        =0,
                            specflip        = False,
                            xgap            =0.,
                            ygap            =0.,
                            ysize           =1.,
                            platescale      =0.135,
                            darkcurr        =0.0,
                            saturation      =65535.*1.255,  # Gain applied
                            nonlinear       =0.86,          # Modified by JXP to go higher
                            numamplifiers   =1,
                            gain            =[1.255],
                            ronoise         =[4.64],
                            datasec         = ['',''],      # These are provided by read_lris
                            oscansec        = ['',''],
                            suffix          ='_01red'
                            ),
                #Detector 2
                pypeitpar.DetectorPar(
                            dataext         =2,
                            specaxis        =0,
                            specflip        = False,
                            xgap            =0.,
                            ygap            =0.,
                            ysize           =1.,
                            platescale      =0.135,
                            darkcurr        =0.,
                            saturation      =65535.*1.162,  # Gain applied
                            nonlinear       =0.86,
                            numamplifiers   =1,
                            gain            =[1.162],
                            ronoise         =[4.62],
                            datasec         = ['',''],      # These are provided by read_lris
                            oscansec        = ['',''],
                            suffix          ='_02red'
                            )]
        self.numhead = 3
        # Uses default timeunit

    def load_raw_frame(self, raw_file, det=None):
        """
        Wrapper to the raw image reader for LRIS

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
        hdu = fits.open(raw_file)
        header = hdu[det].header

        # Grab data (this includes flips as needed)
        data, predata, postdata, x1, y1 = lris_read_amp(hdu, det)
        # Pack
        raw_img = np.zeros((data.shape[0]+predata.shape[0]+postdata.shape[0], data.shape[1]))
        raw_img[:predata.shape[0],:] = predata
        raw_img[predata.shape[0]:predata.shape[0]+data.shape[0],:] = data
        raw_img[-postdata.shape[0]:,:] = postdata

        # Return
        return raw_img.T, header

    def get_image_section(self, inp=None, det=1, section='datasec'):
        #
        hdu = fits.open(inp)
        head0 = hdu[0].header
        binning = head0['BINNING']
        xbin, ybin = [int(ibin) for ibin in binning.split(',')]

        # Get post, pre-pix values
        precol = head0['PRECOL']
        postpix = head0['POSTPIX']
        preline = head0['PRELINE']
        postline = head0['POSTLINE']

        if section == 'datasec':
            datsec = hdu[det].header['DATASEC']  # THIS IS BINNED
            x1, x2, y1, y2 = np.array(parse.load_sections(datsec, fmt_iraf=False)).flatten()
            dy = (y2-y1)+1
            section = '[{:d}:{:d},{:d}:{:d}]'.format(preline*ybin, preline*ybin+(dy)*ybin, x1*xbin, x2*xbin)  # Eliminate lines
        elif section == 'oscansec':
            nx = hdu[det].data.shape[1]
            section = '[:,{:d}:{:d}]'.format(nx*2-postpix, nx*2)
        #
        return [section], False, False, False


def read_lris(raw_file, det=None, TRIM=False):
    """
    Read a raw LRIS data frame (one or more detectors)
    Packed in a multi-extension HDU
    Based on readmhdufits.pro

    Parameters
    ----------
    raw_file : str
      Filename
    det : int, optional
      Detector number; Default = both
    TRIM : bool, optional
      Trim the image?
      This doesn't work....

    Returns
    -------
    array : ndarray
      Combined image 
    header : FITS header
    sections : list
      List of datasec, oscansec, ampsec sections
      datasec, oscansec needs to be for an *unbinned* image as per standard convention
    """

    # Check for file; allow for extra .gz, etc. suffix
    fil = glob.glob(raw_file+'*') 
    if len(fil) != 1:
        msgs.error("Found {:d} files matching {:s}".format(len(fil)))

    # Read
    msgs.info("Reading LRIS file: {:s}".format(fil[0]))
    hdu = fits.open(fil[0])
    head0 = hdu[0].header

    # Get post, pre-pix values
    precol = head0['PRECOL']
    postpix = head0['POSTPIX']
    preline = head0['PRELINE']
    postline = head0['POSTLINE']

    # Setup for datasec, oscansec
    dsec = []
    osec = []

    # get the x and y binning factors...
    binning = head0['BINNING']
    xbin, ybin = [int(ibin) for ibin in binning.split(',')]

    # First read over the header info to determine the size of the output array...
    n_ext = len(hdu)-1  # Number of extensions (usually 4)
    xcol = []
    xmax = 0
    ymax = 0
    xmin = 10000
    ymin = 10000
    for i in np.arange(1, n_ext+1):
        theader = hdu[i].header
        detsec = theader['DETSEC']
        if detsec != '0':
            # parse the DETSEC keyword to determine the size of the array.
            x1, x2, y1, y2 = np.array(parse.load_sections(detsec, fmt_iraf=False)).flatten()

            # find the range of detector space occupied by the data
            # [xmin:xmax,ymin:ymax]
            xt = max(x2, x1)
            xmax = max(xt, xmax)
            yt =  max(y2, y1)
            ymax = max(yt, ymax)

            # find the min size of the array
            xt = min(x1, x2)
            xmin = min(xmin, xt)
            yt = min(y1, y2)
            ymin = min(ymin, yt)
            # Save
            xcol.append(xt)

    # determine the output array size...
    nx = xmax - xmin + 1
    ny = ymax - ymin + 1

    # change size for binning...
    nx = nx // xbin
    ny = ny // ybin

    # Update PRECOL and POSTPIX
    precol = precol // xbin
    postpix = postpix // xbin

    # Deal with detectors
    if det in [1,2]:
        nx = nx // 2
        n_ext = n_ext // 2
        det_idx = np.arange(n_ext, dtype=np.int) + (det-1)*n_ext
        ndet = 1
    elif det is None:
        ndet = 2
        det_idx = np.arange(n_ext).astype(int)
    else:
        raise ValueError('Bad value for det')

    # change size for pre/postscan...
    if not TRIM:
        nx += n_ext*(precol+postpix)
        ny += preline + postline

    # allocate output array...
    array = np.zeros( (nx, ny) )
    order = np.argsort(np.array(xcol))

    # insert extensions into master image...
    for kk, i in enumerate(order[det_idx]):

        # grab complete extension...
        data, predata, postdata, x1, y1 = lris_read_amp(hdu, i+1)
                            #, linebias=linebias, nobias=nobias, $
                            #x1=x1, x2=x2, y1=y1, y2=y2, gaindata=gaindata)
        # insert components into output array...
        if not TRIM:
            # insert predata...
            buf = predata.shape
            nxpre = buf[0]
            xs = kk*precol
            xe = xs + nxpre
            '''
            if keyword_set(VERBOSITY) then begin
                section = '['+stringify(xs)+':'+stringify(xe)+',*]'
                message, 'inserting extension '+stringify(i)+ $
                         ' predata  in '+section, /info
            endif 
            '''
            array[xs:xe, :] = predata

            # insert data...
            buf = data.shape
            nxdata = buf[0]
            nydata = buf[1]
            xs = n_ext*precol + kk*nxdata #(x1-xmin)/xbin
            xe = xs + nxdata
            #section = '[{:d}:{:d},{:d}:{:d}]'.format(preline,nydata-postline, xs, xe)  # Eliminate lines
            section = '[{:d}:{:d},{:d}:{:d}]'.format(preline*ybin, (nydata-postline)*ybin, xs*xbin, xe*xbin)  # Eliminate lines
            dsec.append(section)
            #print('data',xs,xe)
            array[xs:xe, :] = data   # Include postlines

            #; insert postdata...
            buf = postdata.shape
            nxpost = buf[0]
            xs = nx - n_ext*postpix + kk*postpix
            xe = xs + nxpost 
            section = '[:,{:d}:{:d}]'.format(xs*xbin, xe*xbin)
            osec.append(section)
            '''
            if keyword_set(VERBOSITY) then begin
                section = '['+stringify(xs)+':'+stringify(xe)+',*]'
                message, 'inserting extension '+stringify(i)+ $
                         ' postdata in '+section, /info
            endif 
            '''
            array[xs:xe, :] = postdata
        else:
            buf = data.shape
            nxdata = buf[0]
            nydata = buf[1]

            xs = (x1-xmin)//xbin
            xe = xs + nxdata 
            ys = (y1-ymin)//ybin
            ye = ys + nydata - postline

            yin1 = preline
            yin2 = nydata - postline 

            '''
            if keyword_set(VERBOSITY) then begin
                section = '['+stringify(xs)+':'+stringify(xe)+ $
                          ','+stringify(ys)+':'+stringify(ye)+']'
                message, 'inserting extension '+stringify(i)+ $
                         ' data     in '+section, /info
            endif 
            '''
            array[xs:xe, ys:ye] = data[:, yin1:yin2]

    # make sure BZERO is a valid integer for IRAF
    obzero = head0['BZERO']
    head0['O_BZERO'] = obzero
    head0['BZERO'] = 32768-obzero

    # Return, transposing array back to goofy Python indexing
    return array.T, head0, (dsec, osec)


def lris_read_amp(inp, ext):
    """
    Read one amplifier of an LRIS multi-extension FITS image

    Parameters
    ----------
    inp: tuple 
      (str,int) filename, extension
      (hdu,int) FITS hdu, extension

    Returns
    -------
    data
    predata
    postdata
    x1
    y1

    ;------------------------------------------------------------------------
    function lris_read_amp, filename, ext, $
      linebias=linebias, nobias=nobias, $
      predata=predata, postdata=postdata, header=header, $
      x1=x1, x2=x2, y1=y1, y2=y2, GAINDATA=gaindata
    ;------------------------------------------------------------------------
    ; Read one amp from LRIS mHDU image
    ;------------------------------------------------------------------------
    """
    # Parse input
    if isinstance(inp, str):
        hdu = fits.open(inp)
    else:
        hdu = inp

    # Get the pre and post pix values
    # for LRIS red POSTLINE = 20, POSTPIX = 80, PRELINE = 0, PRECOL = 12
    head0 = hdu[0].header
    precol = head0['precol']
    postpix = head0['postpix']

    # Deal with binning
    binning = head0['BINNING']
    xbin, ybin = [int(ibin) for ibin in binning.split(',')]
    precol = precol//xbin
    postpix = postpix//xbin

    # get entire extension...
    temp = hdu[ext].data.transpose() # Silly Python nrow,ncol formatting
    tsize = temp.shape
    nxt = tsize[0]

    # parse the DETSEC keyword to determine the size of the array.
    header = hdu[ext].header
    detsec = header['DETSEC']
    x1, x2, y1, y2 = np.array(parse.load_sections(detsec, fmt_iraf=False)).flatten()

    # parse the DATASEC keyword to determine the size of the science region (unbinned)
    datasec = header['DATASEC']
    xdata1, xdata2, ydata1, ydata2 = np.array(parse.load_sections(datasec, fmt_iraf=False)).flatten()

    # grab the components...
    predata = temp[0:precol, :]
    # datasec appears to have the x value for the keywords that are zero
    # based. This is only true in the image header extensions
    # not true in the main header.  They also appear inconsistent between
    # LRISr and LRISb!
    #data     = temp[xdata1-1:xdata2-1,*]
    #data = temp[xdata1:xdata2+1, :]
    if (xdata1-1) != precol:
        msgs.error("Something wrong in LRIS datasec or precol")
    xshape = 1024 // xbin
    if (xshape+precol+postpix) != temp.shape[0]:
        msgs.warn("Unexpected size for LRIS detector.  We expect you did some windowing...")
        xshape = temp.shape[0] - precol - postpix
    data = temp[precol:precol+xshape,:]
    postdata = temp[nxt-postpix:nxt, :]

    # flip in X as needed...
    if x1 > x2:
        xt = x2
        x2 = x1
        x1 = xt
        data = np.flipud(data) #reverse(temporary(data),1)

    # flip in Y as needed...
    if y1 > y2:
        yt = y2
        y2 = y1
        y1 = yt
        data = np.fliplr(data)
        predata = np.fliplr(predata)
        postdata = np.fliplr(postdata)

    '''
    #; correct gain if requested...
    if keyword_set(GAINDATA) then begin
        gain = gainvalue( gaindata, header)
        data = FLOAT(temporary(data)) * gain
        predata = FLOAT(temporary(predata)) * gain
        postdata = FLOAT(temporary(postdata)) * gain
    endif
    '''

    '''
    ;; optional bias subtraction...
    if ~ keyword_set(NOBIAS) then begin
        if keyword_set( LINEBIAS) then begin
            ;; compute a bias for each line...
            bias = median( postdata, dim=1)

            ;; subtract for data...
            buf = size(data)
            nx = buf[1]
            ny = buf[2]
            data2 = fltarr(nx,ny)
            for i=0,nx-1 do begin
                data2[i,*] = float(data[i,*]) - bias
            endfor 
            data = data2
        endif else begin
            ;; compute a scalar bias....
            bias = median( postdata)
            data -= bias
        endelse
    endif
    '''

    return data, predata, postdata, x1, y1



def convert_lowredux_pixelflat(infil, outfil):
    """ Convert LowRedux pixelflat to PYPIT format
    Returns
    -------

    """
    # Read
    hdu = fits.open(infil)
    data = hdu[0].data

    #
    prihdu = fits.PrimaryHDU()
    hdus = [prihdu]
    prihdu.header['FRAMETYP'] = 'pixelflat'

    # Detector 1
    img1 = data[:,:data.shape[1]//2]
    hdu = fits.ImageHDU(img1)
    hdu.name = 'DET1'
    prihdu.header['EXT0001'] = 'DET1-pixelflat'
    hdus.append(hdu)

    # Detector 2
    img2 = data[:,data.shape[1]//2:]
    hdu = fits.ImageHDU(img2)
    hdu.name = 'DET2'
    prihdu.header['EXT0002'] = 'DET2-pixelflat'
    hdus.append(hdu)

    # Finish
    hdulist = fits.HDUList(hdus)
    hdulist.writeto(outfil, clobber=True)
    print('Wrote {:s}'.format(outfil))

