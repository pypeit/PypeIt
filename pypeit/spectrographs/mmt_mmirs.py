"""
Module for MMT MMIRS

.. include:: ../include/links.rst
"""
from pathlib import Path

import numpy as np
from scipy.signal import savgol_filter

from astropy.table import Table
from astropy.time import Time
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from pypeit import log
from pypeit import PypeItError
from pypeit import telescopes
from pypeit import utils
from pypeit import io
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.core import fitramp
from pypeit.images import detector_container
from pypeit.spectrographs import spectrograph
from pypeit.par import parset


class MMTMMIRSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle MMT/MMIRS specific code
    """
    ndet = 1
    name = 'mmt_mmirs'
    telescope = telescopes.MMTTelescopePar()
    camera = 'MMIRS'
    url = 'https://lweb.cfa.harvard.edu/mmti/mmirs.html'
    header_name = 'mmirs'
    supported = True

    # Up-the-ramp fitting configuration
    ramp_sig_guess = 15.0
    """Initial guess for the single-read noise in electrons."""
    ramp_sig_range = (3.0, 50.0)
    """Allowed range for the calibrated single-read noise (electrons)."""
    ramp_min_dark_groups = 10
    """Minimum number of reads for a dark to be usable for noise calibration."""
    _ramp_dark_files = None
    _ramp_sigma = None

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=1, card='RA')
        self.meta['dec'] = dict(ext=1, card='DEC')
        self.meta['target'] = dict(ext=1, card='OBJECT')
        self.meta['decker'] = dict(ext=1, card='APERTURE')
        self.meta['dichroic'] = dict(ext=1, card='FILTER')
        self.meta['binning'] = dict(ext=1, card=None, default='1,1')

        self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['exptime'] = dict(ext=1, card='EXPTIME')
        self.meta['airmass'] = dict(ext=1, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=1, card='DISPERSE')
        self.meta['idname'] = dict(ext=1, card='IMAGETYP')
        self.meta['instrument'] = dict(ext=1, card='INSTRUME')

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
        # TODO: This should be how we always deal with timeunit = 'isot'. Are
        # we doing that for all the relevant spectrographs?
        if meta_key == 'mjd':
            time = headarr[1]['DATE-OBS']
            ttime = Time(time, format='isot')
            return ttime.mjd
        raise PypeItError("Not ready for this compound meta")

    def cache_metadata(self, fitstbl):
        """
        Record the dark frames in the metadata table for later use in
        calibrating the up-the-ramp single-read noise.

        Cheap and idempotent; the darks are only opened (lazily) by
        :func:`get_ramp_sigma`.

        Args:
            fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
                The class holding the metadata for all the frames.
        """
        tbl = fitstbl.table
        if 'directory' not in tbl.colnames or 'filename' not in tbl.colnames:
            return
        if 'frametype' in tbl.colnames:
            indx = np.array([ft is not None and 'dark' in str(ft)
                             for ft in tbl['frametype']])
        elif 'idname' in tbl.colnames:
            indx = np.array([str(idn).strip().lower() == 'dark'
                             for idn in tbl['idname']])
        else:
            return
        self._ramp_dark_files = [Path(str(d)) / str(f)
                                 for d, f in zip(tbl['directory'][indx],
                                                 tbl['filename'][indx])
                                 if not str(f).startswith('#')]

    def _best_ramp_dark(self):
        """
        Select the recorded dark frame with the most reads.

        Returns:
            `Path`_: Path to the best dark, or None if no recorded dark
            exists with at least :attr:`ramp_min_dark_groups` reads.
        """
        best, best_n = None, 0
        for f in self._ramp_dark_files or []:
            try:
                with io.fits_open(f) as dhdu:
                    n = sum(1 for h in dhdu if h.header.get('NAXIS') == 2
                            and h.header.get('NAXIS1', 0) > 0)
            except (OSError, PypeItError):
                log.warning(f'Could not open recorded dark frame {f}; skipping it '
                            'for read-noise calibration.')
                continue
            if n >= self.ramp_min_dark_groups and n > best_n:
                best, best_n = f, n
        return best

    def get_ramp_sigma(self, diffs, covar):
        """
        Determine the single-read noise for up-the-ramp fitting.

        Preferentially calibrates the noise from a dark frame recorded by
        :func:`cache_metadata` (the result is cached for the rest of the
        run).  If no suitable dark is available, self-calibrates from the
        provided ramp differences (not cached).  The result is clamped to
        :attr:`ramp_sig_range` in both cases.

        Args:
            diffs (`numpy.ndarray`_):
                Scaled resultant differences of the frame being processed,
                shape ``(ndiffs, ny, nx)``, in electrons.
            covar (:class:`~pypeit.core.fitramp.Covar`):
                Covariance object matching ``diffs``.

        Returns:
            :obj:`float`: Single-read noise in electrons.
        """
        if self._ramp_sigma is not None:
            return self._ramp_sigma
        dark = self._best_ramp_dark()
        if dark is not None:
            log.info(f'Calibrating MMIRS single-read noise from dark: {dark.name}')
            with io.fits_open(dark) as dhdu:
                dreads, dhead = mmirs_load_ramp(dhdu)
            ngroups = dreads.shape[0]
            dcovar = fitramp.Covar([dhead['GRPTIME'] * (i + 1)
                                    for i in range(ngroups)])
            ddiffs = mmirs_ramp_diffs(dreads * dhead['GAIN'], dcovar)
            sig = mmirs_calibrate_sigma(ddiffs, dcovar,
                                        sig_guess=self.ramp_sig_guess)
            self._ramp_sigma = float(np.clip(sig, *self.ramp_sig_range))
            log.info(f'Calibrated single-read noise: {self._ramp_sigma:.2f} e-')
            return self._ramp_sigma
        log.info('No suitable dark listed; self-calibrating MMIRS single-read '
                 'noise from the frame itself')
        sig = mmirs_calibrate_sigma(diffs, covar, sig_guess=self.ramp_sig_guess)
        sig = float(np.clip(sig, *self.ramp_sig_range))
        log.info(f'Self-calibrated single-read noise: {sig:.2f} e-')
        return sig

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
        return ['DISPERSE']

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
        # Detector 1
        detector_dict = dict(
            binning='1,1',
            det             = 1,
            dataext         = 1,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.2012,
            darkcurr        = 36.0,  # e-/pixel/hour  (=0.01 e-/pixel/s)
            saturation      = 700000., #155400.,
            nonlinear       = 1.0,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(0.95),
            ronoise         = np.atleast_1d(3.14),
            datasec         = np.atleast_1d('[:,:]'),
            oscansec        = None, #np.atleast_1d('[:,:]')
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

        # Image processing steps
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)
        #par['calibrations']['traceframe']['process']['use_darkimage'] = True
        #par['calibrations']['pixelflatframe']['process']['use_darkimage'] = True
        #par['calibrations']['illumflatframe']['process']['use_darkimage'] = True
        #par['scienceframe']['process']['use_darkimage'] = True
        par['scienceframe']['process']['use_illumflat'] = True

        # Wavelengths
        # 1D wavelength solution with arc lines
        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.125
        par['calibrations']['wavelengths']['sigdetect']=5
        par['calibrations']['wavelengths']['fwhm'] = 4.
        par['calibrations']['wavelengths']['n_first']=2
        par['calibrations']['wavelengths']['n_final']=4
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        par['calibrations']['wavelengths']['match_toler']=5.0

        # Set slits and tilts parameters
        par['calibrations']['tilts']['tracethresh'] = 5
        par['calibrations']['tilts']['spat_order'] = 7
        par['calibrations']['tilts']['spec_order'] = 5
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['edge_thresh'] = 100.
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.4
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['bound_detector'] = True

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 60]
        par['calibrations']['tiltframe']['exprng'] = [60, None]
        par['calibrations']['arcframe']['exprng'] = [60, None]
        par['calibrations']['darkframe']['exprng'] = [30, None]
        par['scienceframe']['exprng'] = [30, None]

        # dark
        # TODO: This is now the default.
        par['calibrations']['darkframe']['process']['apply_gain'] = True

        # cosmic ray rejection
        par['scienceframe']['process']['sigclip'] = 5.0
        par['scienceframe']['process']['objlim'] = 2.0
        par['scienceframe']['process']['grow'] = 0.5

        # Science reduction
        par['reduce']['findobj']['snr_thresh'] = 5.0
        par['reduce']['skysub']['sky_sigrej'] = 5.0
        par['reduce']['findobj']['find_trim_edge'] = [5,5]
        # Object tracing: bound_detector=True gives synthetic straight slit
        # edges, so the initial object trace is a straight line while the real
        # trace is inclined; a linear fit with loose rejection lets the
        # centroids reach and follow the real trace.
        par['reduce']['findobj']['trace_npoly'] = 1
        par['reduce']['findobj']['trace_maxdev'] = 50.
        par['reduce']['findobj']['trace_maxshift'] = 20.
        # Do not correct for flexure
        par['flexure']['spec_method'] = 'skip'

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 8
        # ToDo: replace the telluric grid file for MMT site.
        par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_26000_R10000.fits'

        return par

    def config_specific_par(
            self,
            inp:str|list|Path|fits.Header|Table,
            inp_par:parset.ParSet|None=None
        ) -> parset.ParSet:
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            inp (:obj:`str`, :obj:`list`, `Path`_, `astropy.io.fits.Header`_, `astropy.table.Table`_):
                Input filename, an `astropy.io.fits.Header`_ object, or a list
                of `astropy.io.fits.Header`_ objects.  Or a row from the
                metadata table.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        # Start with instrument-wide parameters
        par = super().config_specific_par(inp, inp_par=inp_par)

       # Adjust parameters based on grating & dichroic used
        grating = self.get_meta_value(inp, 'dispname')
        dichroic = self.get_meta_value(inp, 'dichroic')

        if (grating=='HK') and (dichroic=='zJ'):
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'mmt_mmirs_HK_zJ.fits'
        elif (grating=='K3000') and (dichroic=='Kspec'):
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'mmt_mmirs_K3000_Kspec.fits'
        elif (grating=='J') and (dichroic=='zJ'):
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'mmt_mmirs_J_zJ.fits'

        return par

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
        if ftype in ['pinhole', 'bias']:
            # No pinhole or bias frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (fitstbl['idname'] == 'flat')
        if ftype == 'standard':
            return good_exp & (fitstbl['idname'] == 'object')
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'object')
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'object')
        if ftype == 'dark':
            return good_exp & (fitstbl['idname'] == 'dark')
        log.debug('Cannot determine if frames are of type {0}.'.format(ftype))
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
                Processed bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)

        log.info("Using hard-coded BPM for det=1 on MMIRS")

        # Get the binning
        hdu = io.fits_open(filename)
        binning = hdu[1].header['CCDSUM']
        hdu.close()

        # Apply the mask
        xbin, ybin = int(binning.split(' ')[0]), int(binning.split(' ')[1])
        bpm_img[:, 187 // ybin] = 1

        return bpm_img

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

        Frames with three or more non-destructive reads are combined using
        up-the-ramp fitting (see :func:`_ramp_fit_image`); frames with two
        (or one) reads use correlated double sampling, as before.
        """
        fil = utils.find_single_file(f'{raw_file}*', required=True)

        # Read
        log.info(f'Reading MMIRS file: {fil}')
        hdu = io.fits_open(fil)
        head1 = hdu[1].header

        detector_par = self.get_detector_par(det if det is not None else 1, hdu=hdu)

        # get the x and y binning factors...
        binning = head1['CCDSUM']
        xbin, ybin = [int(ibin) for ibin in binning.split(' ')]

        # Need the exposure time
        exptime = hdu[self.meta['exptime']['ext']].header[self.meta['exptime']['card']]

        # Number of non-destructive reads stored in the file
        n_reads = sum(1 for h in hdu if h.header.get('NAXIS') == 2
                      and h.header.get('NAXIS1', 0) > 0)

        if n_reads >= 3:
            # Up-the-ramp fitting with jump detection
            array, eff_ronoise = self._ramp_fit_image(hdu, detector_par, exptime)
            detector_par['ronoise'] = np.atleast_1d(eff_ronoise)
        else:
            # Correlated double sampling (first minus last read)
            datasec = head1['DATASEC']
            x1, x2, y1, y2 = np.array(parse.load_sections(datasec,
                                                          fmt_iraf=False)).flatten()
            if len(hdu) > 2:
                data = mmirs_read_amp(hdu[1].data.astype('float64')) \
                        - mmirs_read_amp(hdu[2].data.astype('float64'))
            else:
                data = mmirs_read_amp(hdu[1].data.astype('float64'))
            array = data[x1-1:x2, y1-1:y2]

        ## ToDo: This is a hack. Need to solve this issue. I cut at 998 due to the HK zero order contaminating
        ## the blue part of the zJ+HK spectrum. For other setup, you do not need to cut the detector.
        if (head1['FILTER']=='zJ') and (head1['DISPERSE']=='HK'):
            array = array[:int(998/ybin),:]
        rawdatasec_img = np.ones_like(array,dtype='int')
        # NOTE: If there is no overscan, must be set to 0s
        oscansec_img = np.zeros_like(array,dtype='int')

        # Return, transposing array back to orient the overscan properly
        return detector_par, np.flipud(array), hdu, exptime, np.flipud(rawdatasec_img),\
               np.flipud(np.flipud(oscansec_img))

    def _ramp_fit_image(self, hdu, detector_par, exptime):
        """
        Perform up-the-ramp fitting of a multi-read MMIRS frame.

        Args:
            hdu (`astropy.io.fits.HDUList`_):
                Opened raw file with at least 3 non-destructive reads.
            detector_par (:class:`~pypeit.images.detector_container.DetectorContainer`):
                Detector parameters; provides the gain.
            exptime (:obj:`float`):
                Exposure time in seconds, used to scale the fitted count
                rate to total counts.

        Returns:
            :obj:`tuple`: The ramp-fitted image in ADU (`numpy.ndarray`_,
            shape of the trimmed data section) and the effective read noise
            in electrons (:obj:`float`).
        """
        reads, head1 = mmirs_load_ramp(hdu)
        ngroups = reads.shape[0]
        gain = detector_par['gain'][0]
        reads *= gain      # ADU -> electrons
        covar = fitramp.Covar([head1['GRPTIME'] * (i + 1)
                               for i in range(ngroups)])
        diffs = mmirs_ramp_diffs(reads, covar)
        del reads
        sig = self.get_ramp_sigma(diffs, covar)
        log.info(f'Up-the-ramp fitting {ngroups} reads '
                 f'(single-read noise {sig:.2f} e-)')
        countrate = mmirs_fit_ramp(diffs, covar, sig)
        eff_ronoise = mmirs_effective_ronoise(sig, ngroups)
        log.info(f'Effective read noise: {eff_ronoise:.2f} e-')
        return countrate * exptime / gain, eff_ronoise

def mmirs_read_amp(img, namps=32):
    """
    MMIRS has 32 reading out channels. Need to deal with this issue a little
    bit. The pypeit overscan subtraction is not used; reference-pixel correction
    is done here instead.

    Imported from MMIRS IDL pipeline refpix.pro
    """

    # number of channels for reading out
    if namps is None:
        namps = 32

    data_shape = np.shape(img)
    ampsize = int(data_shape[0] / namps)

    refpix1 = np.array([1, 2, 3])
    refpix2 = np.arange(4) + data_shape[0] - 4
    refpix_all = np.hstack([[0, 1, 2, 3], np.arange(4) + data_shape[0] - 4])
    refvec = np.sum(img[:, refpix_all], axis=1) / np.size(refpix_all)
    svec = savgol_filter(refvec, 11, polyorder=5)

    refvec_2d = np.reshape(np.repeat(svec, data_shape[0], axis=0), data_shape)
    img_out = img - refvec_2d

    for amp in range(namps):
        img_out_ref = img_out[np.hstack([refpix1, refpix2]), :]
        ref1, _, _ = sigma_clipped_stats(
            img_out_ref[:, amp * ampsize + 2 * np.arange(int(ampsize / 2))], sigma=3
        )
        ref2, _, _ = sigma_clipped_stats(
            img_out_ref[:, amp * ampsize + 2 * np.arange(int(ampsize / 2)) + 1], sigma=3
        )
        ref12 = (ref1 + ref2) / 2.
        img_out[:, amp * ampsize:(amp + 1) * ampsize] -= ref12

    return img_out


def mmirs_load_ramp(hdu, namps=32):
    """
    Load the non-destructive reads of an MMIRS ramp in time order.

    Image extensions are sorted by ``EXTVER`` (in raw MMIRS files, ext 1 is
    the *final* read and holds the complete metadata).  Each read is
    reference-pixel corrected with :func:`mmirs_read_amp` on the full frame
    and then trimmed to ``DATASEC``.

    Parameters
    ----------
    hdu : `astropy.io.fits.HDUList`_
        Opened raw MMIRS file.
    namps : :obj:`int`, optional
        Number of readout amplifiers passed to :func:`mmirs_read_amp`.

    Returns
    -------
    reads : `numpy.ndarray`_
        Float64 array with shape ``(ngroups, ny, nx)`` holding the
        reference-pixel-corrected reads in ADU, trimmed to ``DATASEC``.
    head1 : `astropy.io.fits.Header`_
        Header of extension 1 (the metadata-complete final read).
    """
    head1 = hdu[1].header
    img_hdus = sorted([h for h in hdu if h.header.get('NAXIS') == 2
                       and h.header.get('NAXIS1', 0) > 0],
                      key=lambda h: h.header['EXTVER'])
    x1, x2, y1, y2 = np.array(parse.load_sections(head1['DATASEC'],
                                                  fmt_iraf=False)).flatten()
    ngroups = len(img_hdus)
    reads = np.empty((ngroups, x2 - x1 + 1, y2 - y1 + 1), dtype=np.float64)
    for i, h in enumerate(img_hdus):
        frame = mmirs_read_amp(h.data.astype(np.float64), namps=namps)
        reads[i] = frame[x1-1:x2, y1-1:y2]
    return reads, head1


def mmirs_ramp_diffs(reads, covar):
    """
    Compute scaled resultant differences from ramp reads.

    Parameters
    ----------
    reads : `numpy.ndarray`_
        Reads in time order, shape ``(ngroups, ny, nx)``, in electrons.
    covar : :class:`~pypeit.core.fitramp.Covar`
        Covariance object providing the time intervals ``delta_t``.

    Returns
    -------
    diffs : `numpy.ndarray`_
        ``(reads[i+1] - reads[i]) / covar.delta_t[i]``, shape
        ``(ngroups-1, ny, nx)``, in e-/s.
    """
    diffs = np.diff(reads, axis=0)
    diffs /= np.asarray(covar.delta_t)[:, None, None]
    return diffs


def mmirs_calibrate_sigma(diffs, covar, sig_guess=15.0, nrows=200):
    """
    Calibrate the single-read noise from ramp differences.

    Fits a subsample of rows without jump detection (two-pass, with the
    count-rate guess clamped to non-negative values to debias the second
    pass) and rescales ``sig_guess`` so that the median chi-squared matches
    the expected degrees of freedom (``ngroups - 2``).

    Parameters
    ----------
    diffs : `numpy.ndarray`_
        Scaled resultant differences, shape ``(ndiffs, ny, nx)``, electrons.
    covar : :class:`~pypeit.core.fitramp.Covar`
        Covariance object matching ``diffs``.
    sig_guess : :obj:`float`, optional
        Initial guess for the single-read noise in electrons.
    nrows : :obj:`int`, optional
        Number of evenly spaced rows (from the central 80% of the detector)
        to include in the calibration.

    Returns
    -------
    :obj:`float`
        Calibrated single-read noise in electrons (unclamped).
    """
    ndiffs, ny, nx = diffs.shape
    margin = int(ny * 0.10)
    row_candidates = np.arange(margin, ny - margin)
    nrows = min(nrows, len(row_candidates))
    indices = np.linspace(0, len(row_candidates) - 1, nrows, dtype=int)
    sig_row = np.full(nx, sig_guess, dtype=np.float64)
    all_chisq = []
    for row in row_candidates[indices]:
        result = fitramp.fit_ramps(diffs[:, row, :], covar, sig_row)
        guess = result.countrate * (result.countrate > 0)
        result = fitramp.fit_ramps(diffs[:, row, :], covar, sig_row,
                                   countrateguess=guess)
        all_chisq.append(result.chisq)
    median_chisq = float(np.median(np.concatenate(all_chisq)))
    expected_chisq = float(ndiffs - 1)
    return sig_guess * np.sqrt(median_chisq / expected_chisq)


def mmirs_fit_ramp(diffs, covar, sig):
    """
    Fit all pixels of a ramp, row by row, with jump detection.

    Parameters
    ----------
    diffs : `numpy.ndarray`_
        Scaled resultant differences, shape ``(ndiffs, ny, nx)``, electrons.
    covar : :class:`~pypeit.core.fitramp.Covar`
        Covariance object matching ``diffs``.
    sig : :obj:`float`
        Single-read noise in electrons.

    Returns
    -------
    countrate : `numpy.ndarray`_
        Fitted count rates in e-/s, shape ``(ny, nx)``.
    """
    ndiffs, ny, nx = diffs.shape
    countrate = np.empty((ny, nx), dtype=np.float64)
    sig_row = np.full(nx, sig, dtype=np.float64)
    for row in range(ny):
        diffs2use, guess = fitramp.mask_jumps(diffs[:, row, :], covar, sig_row)
        result = fitramp.fit_ramps(diffs[:, row, :], covar, sig_row,
                                   diffs2use=diffs2use,
                                   countrateguess=guess * (guess > 0))
        countrate[row] = result.countrate
    return countrate


def mmirs_effective_ronoise(sig, ngroups):
    """
    Effective read noise of an up-the-ramp-fitted image.

    For ``N`` uniformly spaced reads with single-read noise ``sig``, the
    read-noise contribution to the total-count uncertainty of the fitted
    slope is ``sig * sqrt(12 (N-1) / (N (N+1)))`` (Brandt 2024a).

    Parameters
    ----------
    sig : :obj:`float`
        Single-read noise in electrons.
    ngroups : :obj:`int`
        Number of reads in the ramp.

    Returns
    -------
    :obj:`float`
        Effective read noise in electrons.
    """
    return sig * np.sqrt(12. * (ngroups - 1) / (ngroups * (ngroups + 1)))
