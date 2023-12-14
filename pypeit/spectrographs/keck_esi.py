"""
Module for Keck/ESI specific methods.

.. include:: ../include/links.rst
"""
import os

from IPython import embed

import numpy as np

from astropy.time import Time
import datetime

from pypeit import msgs
from pypeit import telescopes
from pypeit import io
from pypeit.core import framematch
from pypeit.core import parse
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container


class KeckESISpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/ESI specific code
    """
    ndet = 1
    name = 'keck_esi'
    camera = 'ESI'
    header_name = 'ESI'
    #url = 'https://www.lco.cl/?epkb_post_type_1=mage'
    ech_fixed_format = True
    telescope = telescopes.KeckTelescopePar()
    pypeline = 'Echelle'
    supported = True
    #comment = 'See :doc:`mage`'

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
        # TODO: Could this be detector dependent??
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')

        # Detector 1
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            # plate scale in arcsec/pixel
            platescale      = 0.1542,
            # electrons/pixel/hour. 
            darkcurr        = 2.10, # e/pixel/hour... Note : Could be updated
            saturation      = 65535.,
            # CCD is linear to better than 0.5 per cent up to digital saturation (65,536 DN including bias) in the Fast readout mode.
            nonlinear       = 0.99,
            mincounts       = -1e10,
            numamplifiers   = 2,
            gain            = np.atleast_1d([1.3, 1.3]), 
            ronoise         = np.atleast_1d([2.5, 2.5]),
            )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        # Bias
        #par['calibrations']['biasframe']['useframe'] = 'overscan'
        # Wavelengths
        # 1D wavelength solution
        # This is for 1x1
        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.103
        par['calibrations']['wavelengths']['fwhm'] = 2.9
        #
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['lamps'] = ['CuI', 'ArI', 'NeI', 'HgI', 'XeI', 'ArII']

        par['calibrations']['wavelengths']['method'] = 'reidentify'
        par['calibrations']['wavelengths']['cc_thresh'] = 0.50
        par['calibrations']['wavelengths']['cc_local_thresh'] = 0.50

        # Reidentification parameters
        par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_esi_ECH.fits'
        #par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] = 'nothing'

        # Set slits and tilts parameters
        par['calibrations']['tilts']['tracethresh'] = 10. #[10]*self.norders
        par['calibrations']['slitedges']['fit_order'] = 5
        par['calibrations']['slitedges']['max_shift_adj'] = 3.
        par['calibrations']['slitedges']['left_right_pca'] = True

        par['calibrations']['slitedges']['edge_thresh'] = 5.0
        par['calibrations']['slitedges']['det_min_spec_length'] = 0.2
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.4
        par['calibrations']['slitedges']['pca_sigrej'] = 1.5
        par['calibrations']['slitedges']['pca_order'] = 3
        par['calibrations']['slitedges']['add_missed_orders'] = True
        # Find object parameters
        par['reduce']['findobj']['find_trim_edge'] = [4,4]    # Slit is too short to trim 5,5 especially with 2x binning
        par['reduce']['findobj']['maxnumber_sci'] = 2  # Slit is narrow so allow one object per order
        par['reduce']['findobj']['maxnumber_std'] = 1  # Slit is narrow so allow one object per order
        par['reduce']['extraction']['model_full_slit'] = True  # local sky subtraction operates on entire slit

        # Scattered light
        par['calibrations']['pixelflatframe']['process']['subtract_scattlight'] = True
        par['calibrations']['illumflatframe']['process']['subtract_scattlight'] = True
        par['scienceframe']['process']['subtract_scattlight'] = True

        # Always flux calibrate, starting with default parameters
        # Do not correct for flexure
        par['flexure']['spec_method'] = 'skip'
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 60]
        par['calibrations']['arcframe']['exprng'] = [300, None] # Allow for CuAr which can be quite long
        par['calibrations']['darkframe']['exprng'] = [1, None]
        par['scienceframe']['exprng'] = [60, None]
        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='TARGNAME')
        
        self.meta['decker'] = dict(ext=0, card='SLMSKNAM')
        self.meta['binning'] = dict(card=None, compound=True)
        self.meta['mjd'] = dict(card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='ELAPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(card=None, compound=True)
        self.meta['idname'] = dict(ext=0, card='OBSTYPE')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')

        # Lamps -- Have varied in time..
        # From Jim Lyke in the PypeIt slack:
        # HgNe = LAMPAR1 = on
        # CuAr = LAMPCU1 = on
        # Xe   = LAMPNE1 = on

        self.meta['lampstat01'] = dict(ext=0, card='LAMPAR1')
        self.meta['lampstat02'] = dict(ext=0, card='LAMPCU1')
        self.meta['lampstat03'] = dict(ext=0, card='LAMPNE1')
        self.meta['lampstat04'] = dict(ext=0, card='LAMPNE2')
        self.meta['lampstat05'] = dict(ext=0, card='LAMPQTZ1')
        self.meta['lampstat06'] = dict(ext=0, card='FLIMAGIN')
        self.meta['lampstat07'] = dict(ext=0, card='FLSPECTR')

        # Hatch
        self.meta['hatch'] = dict(ext=0, card='HATCHPOS')

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
            binspatial, binspec = parse.parse_binning(headarr[0]['BINNING'])
            return parse.binning2string(binspec, binspatial)
        elif meta_key == 'mjd':
            mjd = headarr[0].get('MJD-OBS', None)
            if mjd is not None:
                mjd_time = Time(mjd,format="mjd")
                # The MJD header value is often invalid in a way that gives it a year > 9000, So sanity check it
                try:
                    mjd_year = mjd_time.to_value(format="decimalyear") 
                    if mjd_year >= datetime.MINYEAR and mjd_year < 9000:
                        return mjd_time.mjd
                except Exception as e:
                    # A problem parsing the MJD, we'll try DATE-OBS and UT
                    msgs.warn("Problem parsing MJD-OBS, trying DATE-OBS and UT instead.")
                    pass             
            return Time('{}T{}'.format(headarr[0]['DATE-OBS'], headarr[0]['UT'])).mjd
        elif meta_key == 'dispname':
            if headarr[0]['PRISMNAM'] == 'in':
                dname = 'Echellette'
            else: # TODO -- Figure out prism and imaging modes
                dname = 'UNKNWN'
            return dname
        else:
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
        return []

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
        if ftype in ['pinhole', 'dark']:
            # No pinhole or pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['bias']:
            return fitstbl['idname'] == 'Bias'
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            ans = np.zeros(len(fitstbl), dtype=bool)
            for kk, idnm in enumerate(fitstbl['idname']):
                if idnm in ['DmFlat', 'IntFlat', 'SkyFlat']:
                    ans[kk] = True
            return ans
        if ftype in ['arc', 'tilt']:
            return fitstbl['idname'] == 'Line'
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'Object') 
        if ftype == 'standard':
            return good_exp & (fitstbl['idname'] == 'Object') 
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

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
                Master bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)


        # Get the binning 
        msgs.info("Custom bad pixel mask for ESI")
        hdu = io.fits_open(filename)
        binspatial, binspec = parse.parse_binning(hdu[0].header['BINNING'])
        hdu.close()

        # Binning Independent Masking 
        bpm_img[:,0:(2//binspatial)] = 1
        # Mask out the 'hotspot' on the upper left coner 
        bpm_img[(3842//binspec):(3944//binspec), (19//binspatial):(161//binspatial)] = 1
        # Mask out the 'bad columns' on the upper left coner
        bpm_img[(2642//binspec):, (418//binspatial):(442//binspatial)] = 1

        # Return
        return bpm_img

    def scattered_light_archive(self, binning, dispname):
        """Archival model parameters for the scattered light. These are based on best fits to currently available data.

        Parameters
        ----------
        binning : :obj:`str`, optional
            Comma-separated binning along the spectral and spatial directions; e.g., ``2,1``
        dispname : :obj:`str`, optional
            Name of the disperser

        Returns
        -------
        x0 : `numpy.ndarray`_
            A 1D array containing the best-fitting model parameters
        bounds : :obj:`tuple`
            A tuple of two elements, containing two `numpy.ndarray`_ of the same length as x0. These
            two arrays contain the lower (first element of the tuple) and upper (second element of the tuple)
            bounds to consider on the scattered light model parameters.
        """
        # Grab the binning for convenience
        specbin, spatbin = parse.parse_binning(binning)

        # Get some starting parameters (these were determined by fitting spectra,
        # and should be close to the final fitted values to reduce computational time)
        # Note :: These values need to be originally based on data that uses 1x1 binning,
        # and are now scaled here according to the binning of the current data to be analysed.
        # These parameters give a cost of 8.0517e+08 with the science frame used as scattlight (1x1 binning, pad=5)
        x0 = np.array([272.33958742493064 / specbin, 115.501464689107 / spatbin,  # Gaussian kernel widths
                       272.3418000034377 / specbin, 168.0591427733949 / spatbin,  # Lorentzian kernel widths
                       -141.2552517318941 / specbin, 79.25936221285629 / spatbin,  # pixel offsets
                       1.0877734248786808, 1.0562808322123667,  # Zoom factor (spec, spat)
                       5.876311151022701,  # constant flux offset
                       0.0444248025888341,  # kernel angle
                       0.6090358292193677,  # Relative kernel scale (>1 means the kernel is more Gaussian, >0 but <1 makes the profile more lorentzian)
                       0.135392229831296, -0.16167521454188258, # Polynomial terms (coefficients of "spat" and "spat*spec")
                       0.06148093592863097, 0.10305719952486242])  # Polynomial terms (coefficients of spec**index)

        # Now set the bounds of the fitted parameters
        bounds = ([# Lower bounds:
                      1, 1,  # Gaussian kernel widths
                      1, 1,  # Lorentzian kernel widths
                      -200 / specbin, -200 / spatbin,  # pixel offsets
                      0, 0,  # Zoom factor (spec, spat)
                      -1000, -(10 / 180) * np.pi, 0.0,  # constant flux offset, kernel angle, relative kernel scale
                      -10, -10, -10, -10],  # Polynomial terms
                  # Upper bounds
                     [600 / specbin, 600 / spatbin,  # Gaussian kernel widths
                      600 / specbin, 600 / spatbin,  # Lorentzian kernel widths
                      200 / specbin, 200 / spatbin,  # pixel offsets
                      2, 2,  # Zoom factor (spec, spat)
                      1000.0, +(10 / 180) * np.pi, 1000.0,  # constant flux offset, kernel angle, relative kernel scale
                      10, 10, 10, 10])  # Polynomial terms

        # Return the best-fitting archival parameters and the bounds
        return x0, bounds

    @property
    def norders(self):
        """
        Number of orders for this spectrograph. Should only defined for
        echelle spectrographs, and it is undefined for the base class.
        """
        return 10   # 15-6

    @property
    def order_spat_pos(self):
        """
        Return the expected spatial position of each echelle order.
        """
        return np.array([0.115, 0.245, 0.362, 0.465, 0.558, 0.642, 0.719, 0.791, 0.861, 0.933])

    @property
    def order_spat_width(self):
        """
        Return the expected spatial width of each slit trace for each order,
        relative to the spatial size of the detector.
        """
        return np.array([0.0879, 0.0818, 0.0779, 0.0747, 0.0720, 0.0696, 0.0676, 0.0658, 0.0640,
                         0.0617])

    @property
    def orders(self):
        """
        Return the order number for each echelle order.
        """
        return np.arange(15, 5, -1, dtype=int)

    @property
    def spec_min_max(self):
        """
        Return the minimum and maximum spectral pixel expected for the
        spectral range of each order.
        """
        spec_max = np.full(self.norders, np.inf)
        spec_min = np.full(self.norders, -np.inf)
        return np.vstack((spec_min, spec_max))

    def order_platescale(self, order_vec, binning=None):
        """
        Return the platescale for each echelle order.

        This routine is only defined for echelle spectrographs, and it is
        undefined in the base class.

        Args:
            order_vec (`numpy.ndarray`_):
                The vector providing the order numbers.
            binning (:obj:`str`, optional):
                The string defining the spectral and spatial binning.

        Returns:
            `numpy.ndarray`_: An array with the platescale for each order
            provided by ``order``.
        """
        norders = len(order_vec)
        binspatial, binspec = parse.parse_binning(binning)
        # Plate scales
        unbinned_pscale = [0.120, #15 
                           0.127, 0.134, 0.137, 0.144, 0.149, 0.153, 
                           0.158, 0.163, 
                           0.168, # 6 
                           ]

        return np.array(unbinned_pscale)*binspatial



    def get_rawimage(self, raw_file, det, spectrim=None):
        """ Read the image
        """
        # Check for file; allow for extra .gz, etc. suffix
        if not os.path.isfile(raw_file):
            msgs.error(f'{raw_file} not found!')
        hdu = io.fits_open(raw_file)
        head0 = hdu[0].header

        # Number of AMPS
        namp = head0['NUMAMPS']

        # Get post, pre-pix values
        prepix = head0['PREPIX']
        postpix = head0['POSTPIX']
        preline = head0['PRELINE']
        postline = head0['POSTLINE']

        # Grab the data
        full_image = hdu[0].data.astype(float)
        rawdatasec_img = np.zeros_like(full_image, dtype=int)
        oscansec_img = np.zeros_like(full_image, dtype=int)

        # 
        nspat = int(head0['WINDOW'].split(',')[3]) // namp
        for amp in range(namp):
            col0 = prepix*2 + nspat*amp
            # Data
            rawdatasec_img[:, col0:col0+nspat] = amp+1
            # Overscan
            o0 = prepix*2 + nspat*namp + postpix*amp
            oscansec_img[:, o0:o0+postpix] = amp+1

        #embed(header='435 of keck_esi.py')

        return self.get_detector_par(1, hdu=hdu), \
                full_image, hdu, head0['ELAPTIME'], rawdatasec_img, oscansec_img