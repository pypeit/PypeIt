"""
Module for Keck/NIRES specific methods.

.. include:: ../include/links.rst
"""
import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container

from IPython import embed


class KeckNIRESSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/NIRES specific code
    """
    ndet = 1
    name = 'keck_nires'
    telescope = telescopes.KeckTelescopePar()
    camera = 'NIRES'
    url = 'https://www2.keck.hawaii.edu/inst/nires/'
    header_name = 'NIRES'
    pypeline = 'Echelle'
    ech_fixed_format = True
    supported = True
    comment = 'see :doc:`keck_nires`'


    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Args:
            det (:obj:`int`):
                1-indexed detector number.  This is not used because NIRES only
                has one detector!
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        # Detector 1
        detector_dict = dict(
            binning='1,1',
            det=1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = True,
            spatflip=False,
            platescale      = 0.15,
            darkcurr        = 0.01,
            saturation      = 1e6, # I'm not sure we actually saturate with the DITs???
            nonlinear       = 0.76,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(3.8),
            ronoise         = np.atleast_1d(5.0),
            datasec         = np.atleast_1d('[:,:]'),
            oscansec        = None, #np.atleast_1d('[980:1024,:]')  # Is this a hack??
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

        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20 #0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect']=5.0
        par['calibrations']['wavelengths']['fwhm']= 5.0
        par['calibrations']['wavelengths']['n_final']= [3,4,4,4,4]
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'reidentify'
        # Reidentification parameters
        par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_nires.fits'
#        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 6
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.4
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['fwhm_gaussian'] = 4.0

        # Tilt parameters
        par['calibrations']['tilts']['tracethresh'] =  10.0
        #par['calibrations']['tilts']['spat_order'] =  3
        #par['calibrations']['tilts']['spec_order'] =  3

        # Processing steps
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0

        # Flexure
        par['flexure']['spec_method'] = 'skip'

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] ='nothing'
        par['reduce']['extraction']['boxcar_radius'] = 0.75  # arcsec

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 60]
        par['calibrations']['arcframe']['exprng'] = [61, None]
        par['calibrations']['tiltframe']['exprng'] = [61, None]
        par['scienceframe']['exprng'] = [61, None]

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 8
        par['sensfunc']['IR']['maxiter'] = 2
        par['sensfunc']['IR']['telgridfile'] = 'TelFit_MaunaKea_3100_26100_R20000.fits'

        # COADD2D
        # set offsets for coadd2d
        par['coadd2d']['offsets'] = 'header'

        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='TARGNAME')
        self.meta['decker'] = dict(ext=0, card=None, default='0.55 slit')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')

        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='ITIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='INSTR')
        self.meta['idname'] = dict(ext=0, card='OBSTYPE')
        self.meta['frameno'] = dict(ext=0, card='FRAMENUM')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')

        # Dithering
        self.meta['dithpat'] = dict(ext=0, card='DPATNAME')
        self.meta['dithpos'] = dict(card=None, compound=True)
        self.meta['dithoff'] = dict(ext=0, card='YOFFSET')

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
        if meta_key == 'dithpos':
            # the dither positions in NIRES are expressed in numbers 1,2,3,4.
            # We want to convert those into A,B,C. We need to know the dither pattern to do so.
            dpat = headarr[0].get('DPATNAME')
            dpos = headarr[0].get('DPATIPOS')
            if dpos is not None and ((dpos in [1, 2, 3, 4] and dpat in ["ABBA", "ABBAprime", "ABAB"]) or
                                     (dpos in [1, 2, 3] and dpat == 'ABC') or
                                     (dpos in [1, 2] and dpat == 'ABpat')):
                return dpat[dpos-1]
            else:
                return None
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
        return ['dispname']

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
        return ['INSTR']

    def get_comb_group(self, fitstbl):
        """
        Automatically assign combination groups and background images by parsing
        known dither patterns.

        This method is used in
        :func:`~pypeit.metadata.PypeItMetaData.set_combination_groups`, and
        directly modifies the ``comb_id`` and ``bkg_id`` columns in the provided
        table.

        Specifically here for NIRES, since it's likely to have one set of flat/dark frames for
        different targets, this method sets calib = "all" for the flat and dark frames and
        assigns different calib values to the science/standard frames of different targets.

        Moreover, this method parses from the header the dither pattern of the
        science/standard frames in a given calibration group and assigns to each
        of them a default ``comb_id`` and ``bkg_id``. The dither patterns used
        here are: "ABAB", "ABBA", "ABpat", and "ABC".  Note that the frames in
        the same dither positions (A positions or B positions) of each "ABAB" or
        "ABBA" sequence are 2D coadded  (without optimal weighting) before the
        background subtraction, while for the other dither patterns (e.g.,
        "ABpat"), the frames in the same dither positions are not coadded.  The
        ``comb_id`` and ``bkg_id`` will *not* assigned if:

            - the dither offset is zero for every frame in the dither sequence

            - the dither pattern recorded in the header is not recognized or set
              to NONE or MANUAL.

        Args:
            fitstbl(`astropy.table.Table`_):
                The table with the metadata for all the frames.

        Returns:
            `astropy.table.Table`_: modified fitstbl.
        """
        #TODO incorporate parse_dither_pattern() here.

        if 'calib' in fitstbl.keys():
            # find index of fitstbl that contains dark, pixelflat, trace, or lampoffflats frames
            flat_idx = np.array(['pixelflat' in _tab for _tab in fitstbl['frametype']]) | \
                       np.array(['trace' in _tab for _tab in fitstbl['frametype']]) | \
                       np.array(['lampoffflats' in _tab for _tab in fitstbl['frametype']]) | \
                       np.array(['dark' in _tab for _tab in fitstbl['frametype']])
            # set calib for those frames to "all" since it's likely that the same flats are used for different targets
            fitstbl['calib'][flat_idx] = 'all'
            # initialize target calib
            targ_calib = 0

        # find index of fitstbl that contains science and standard frames
        # where science
        sci_idx = np.array(['science' in _tab for _tab in fitstbl['frametype']])
        # where standard
        std_idx = np.array(['standard' in _tab for _tab in fitstbl['frametype']])

        sci_std_idx = [sci_idx, std_idx]
        # loop over the science and standard frames
        for idx in sci_std_idx:
            setups = np.unique(fitstbl[idx]['setup'])
            # loop over the setups (generally there is only one setup, but we check anyway)
            for setup in setups:
                in_cfg = idx & np.array([setup in _set for _set in fitstbl['setup']])
                if len(fitstbl[in_cfg]) == 1:
                    continue
                # generally there is only one setup, but different targets. We want to separate those.
                # how may targets are in this setup?
                targets = np.unique(fitstbl[in_cfg]['target'])
                # loop through targets
                for targ in targets:
                    # where this targ
                    targ_idx = in_cfg & (fitstbl['target'] == targ)
                    if 'calib' in fitstbl.keys():
                        # set different calib for different targs
                        if 'science' in fitstbl['frametype'][targ_idx][0] or \
                           ('standard' in fitstbl['frametype'][targ_idx][0] and 'arc' in fitstbl['frametype'][targ_idx][0]):
                            fitstbl['calib'][targ_idx] = str(targ_calib)
                        elif 'standard' in fitstbl['frametype'][targ_idx]:
                            # find the science frames
                            sci_in_cfg = sci_idx & np.array([setup in _set for _set in fitstbl['setup']])
                            if len(fitstbl[sci_in_cfg]) > 0:
                                # find the closest (in time) science frame to the standard target
                                close_idx = np.argmin(np.absolute(fitstbl[sci_in_cfg]['mjd'] - fitstbl[targ_idx]['mjd'][0]))
                                fitstbl['calib'][targ_idx] = fitstbl['calib'][sci_in_cfg][close_idx]
                        targ_calib += 1

                    # how many dither patterns are used for the selected science/standard target?
                    uniq_dithpats = np.unique(fitstbl[targ_idx]['dithpat'])
                    # loop through the dither patterns
                    for dpat in uniq_dithpats:
                        if dpat in ['NONE', 'none', 'MANUAL']:
                            continue
                        # where this dpat
                        dpat_idx = targ_idx & (fitstbl['dithpat'] == dpat)
                        # get doff
                        doff = fitstbl[dpat_idx]['dithoff'].data.astype(int)

                        # compute comb_id
                        if len(fitstbl[dpat_idx]) > 1 and np.any(doff != 0):
                            # get default combid and bkgid
                            combid = np.copy(fitstbl['comb_id'][dpat_idx].data)
                            bkgid = np.copy(fitstbl['bkg_id'][dpat_idx].data)
                            dpos = fitstbl[dpat_idx]['dithpos'].data

                            if dpat == "ABAB":
                                # find the starting index of the ABAB sequence
                                dpos_idx = np.where((dpos == "A") & (np.roll(dpos, -1) == "B") &
                                                    (np.roll(dpos, -2) == "A") & (np.roll(dpos, -3) == "B"))[0]
                                for i in dpos_idx:
                                    # make sure that that dither offsets are correct
                                    if i < len(dpos) - 3 and doff[i] == doff[i+2] and doff[i+1] == doff[i+3]:
                                        bkgid[i] = combid[i+1]
                                        bkgid[i+1] = combid[i]
                                        combid[i+2] = combid[i]
                                        bkgid[i+2] = bkgid[i]
                                        combid[i+3] = combid[i+1]
                                        bkgid[i+3] = bkgid[i+1]

                            elif dpat in ["ABBA", "ABBAprime"]:
                                # find the starting index of the ABBA sequence
                                dpos_idx = np.where((dpos == "A") & (np.roll(dpos, -1) == "B") &
                                                    (np.roll(dpos, -2) == "B") & (np.roll(dpos, -3) == "A"))[0]
                                for i in dpos_idx:
                                    # make sure that that dither offsets are correct
                                    if i < len(dpos) - 3 and doff[i] == doff[i+3] and doff[i+1] == doff[i+2]:
                                        bkgid[i] = combid[i+1]
                                        bkgid[i+1] = combid[i]
                                        combid[i+2] = combid[i+1]
                                        bkgid[i+2] = bkgid[i+1]
                                        combid[i+3] = combid[i]
                                        bkgid[i+3] = bkgid[i]

                            elif dpat == "ABpat":
                                # find the starting index of the ABpat sequence
                                dpos_idx = np.where((dpos == "A") & (np.roll(dpos, -1) == "B"))[0]
                                for i in dpos_idx:
                                    # exclude when np.roll counts the 1st element of dpos to be in a
                                    # sequence with the last element
                                    if i < len(dpos)-1:
                                        bkgid[i] = combid[i+1]
                                        bkgid[i+1] = combid[i]

                            elif dpat == "ABC":
                                # find the starting index of the ABC sequence
                                dpos_idx = np.where((dpos == "A") & (np.roll(dpos, -1) == "B") &
                                                    (np.roll(dpos, -2) == "C"))[0]
                                for i in dpos_idx:
                                    # exclude when np.roll counts the 1st element of dpos to be in a
                                    # sequence with the last element
                                    if i < len(dpos) - 2:
                                        bkgid[i] = combid[i+1]
                                        bkgid[i+1] = combid[i+2]
                                        bkgid[i+2] = combid[i+1]

                            # assign bkgid for files that are part of an incomplete sequence
                            for i in range(len(fitstbl[dpat_idx])):
                                # initialize pos_idx
                                pos_idx = None
                                # if A frame doesn't have bkgid assigned
                                if bkgid[i] == -1 and (fitstbl[dpat_idx]['dithpos'][i] == "A"):
                                    # find B frames to subtract from this A
                                    pos_idx = fitstbl[dpat_idx]['dithpos'] == "B"
                                # if B frame doesn't have bkgid assigned
                                elif bkgid[i] == -1 and (fitstbl[dpat_idx]['dithpos'][i] == "B"):
                                    # find A frames to subtract from this B
                                    pos_idx = fitstbl[dpat_idx]['dithpos'] == "A"
                                # if C frame doesn't have bkgid assigned
                                elif bkgid[i] == -1 and (fitstbl[dpat_idx]['dithpos'][i] == "C"):
                                    # find A or B frames to subtract from this C
                                    pos_idx = (fitstbl[dpat_idx]['dithpos'] == "A") | (fitstbl[dpat_idx]['dithpos'] == "B")
                                # assign bkgid
                                if pos_idx is not None and np.any(pos_idx):
                                    # pick closest (in mjd)
                                    close_idx = np.argmin(np.absolute(fitstbl[dpat_idx][pos_idx]['mjd'] - fitstbl[dpat_idx]['mjd'][i]))
                                    bkgid[i] = combid[pos_idx][close_idx]

                            fitstbl['bkg_id'][dpat_idx] = bkgid
                            fitstbl['comb_id'][dpat_idx] = combid

        return fitstbl

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        pypeit_keys = super().pypeit_file_keys()
        pypeit_keys += ['dithpat', 'dithpos', 'dithoff','frameno']
        return pypeit_keys

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
        if ftype in ['pinhole', 'bias', 'dark']:
            # No pinhole or bias frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype == 'standard':
            return good_exp & ((fitstbl['idname'] == 'object') | (fitstbl['idname'] == 'Object') |
                               (fitstbl['idname'] == 'standard') | (fitstbl['idname'] == 'telluric')) \
                   & (fitstbl['target'] != 'DOME PHLAT')
        if ftype == 'lampoffflats':
            return good_exp & ((fitstbl['idname'] == 'dark') | (fitstbl['idname'] == 'Dark'))
        if ftype in ['pixelflat', 'trace']:
            return (fitstbl['idname'] == 'domeflat') | (fitstbl['idname'] == 'DomeFlat')
        if ftype in 'science':
            return good_exp & ((fitstbl['idname'] == 'object') | (fitstbl['idname'] == 'Object')) \
                   & (fitstbl['target'] != 'DOME PHLAT')
        if ftype in ['arc', 'tilt']:
            return good_exp & ((fitstbl['idname'] == 'object') | (fitstbl['idname'] == 'Object')) \
                   & (fitstbl['target'] != 'DOME PHLAT')
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
                Processed bias frame used to identify bad pixels.

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        msgs.info("Custom bad pixel mask for NIRES")
        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)

        if det == 1:
            bpm_img[:, :20] = 1.
            bpm_img[:, 1000:] = 1.

        return bpm_img

    @property
    def norders(self):
        """
        Number of orders for this spectograph. Should only defined for
        echelle spectrographs, and it is undefined for the base class.
        """
        return 5

    @property
    def order_spat_pos(self):
        """
        Return the expected spatial position of each echelle order.
        """
        return np.array([0.22773035, 0.40613574, 0.56009658, 0.70260714, 0.86335914])

    @property
    def orders(self):
        """
        Return the order number for each echelle order.
        """
        return np.arange(7, 2, -1, dtype=int)

    @property
    def spec_min_max(self):
        """
        Return the minimum and maximum spectral pixel expected for the
        spectral range of each order.
        """
        spec_max = np.asarray([np.inf]*self.norders)
        spec_min = np.asarray([1024, -np.inf, -np.inf, -np.inf, -np.inf])
        return np.vstack((spec_min, spec_max))

    def order_platescale(self, order_vec, binning=None):
        """
        Return the platescale for each echelle order.

        Note that NIRES has no binning.

        Args:
            order_vec (`numpy.ndarray`_):
                The vector providing the order numbers.
            binning (:obj:`str`, optional):
                The string defining the spectral and spatial binning. **This
                is always ignored.**

        Returns:
            `numpy.ndarray`_: An array with the platescale for each order
            provided by ``order``.
        """
        return np.full(order_vec.size, 0.15)

