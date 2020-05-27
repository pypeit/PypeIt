"""
Module for MMT MMIRS
"""
import glob
import numpy as np
from scipy.signal import savgol_filter
from astropy.time import Time
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.images import detector_container
from pypeit.spectrographs import spectrograph
from pkg_resources import resource_filename


class MMTMMIRSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle MMT/MMIRS specific code


    """
    ndet = 1

    def __init__(self):
        # Get it started
        super(MMTMMIRSSpectrograph, self).__init__()
        self.spectrograph = 'mmt_mmirs'
        self.telescope = telescopes.MMTTelescopePar()

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for the reductions.
        """
        par = pypeitpar.PypeItPar()
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
        meta['ra'] = dict(ext=1, card='RA')
        meta['dec'] = dict(ext=1, card='DEC')
        meta['target'] = dict(ext=1, card='OBJECT')
        meta['decker'] = dict(ext=1, card='APERTURE')
        meta['dichroic'] = dict(ext=1, card='FILTER')
        meta['binning'] = dict(ext=1, card=None, default='1,1')

        meta['mjd'] = dict(ext=0, card=None, compound=True)
        meta['exptime'] = dict(ext=1, card='EXPTIME')
        meta['airmass'] = dict(ext=1, card='AIRMASS')
        # Extras for config and frametyping
        meta['dispname'] = dict(ext=1, card='DISPERSE')
        meta['idname'] = dict(ext=1, card='IMAGETYP')

        # Ingest
        self.meta = meta

    def compound_meta(self, headarr, meta_key):
        if meta_key == 'mjd':
            time = headarr[1]['DATE-OBS']
            ttime = Time(time, format='isot')
            return ttime.mjd
        else:
            msgs.error("Not ready for this compound meta")

    def get_detector_par(self, hdu, det):
        """
        Return a DectectorContainer for the current image

        Args:
            hdu (`astropy.io.fits.HDUList`):
                HDUList of the image of interest.
                Ought to be the raw file, or else..
            det (int):

        Returns:
            :class:`pypeit.images.detector_container.DetectorContainer`:

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
            darkcurr        = 0.01,
            saturation      = 700000., #155400.,
            nonlinear       = 1.0,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(0.95),
            ronoise         = np.atleast_1d(3.14),
            datasec         = np.atleast_1d('[:,:]'),
            oscansec        = np.atleast_1d('[:,:]')
            )
        return detector_container.DetectorContainer(**detector_dict)

    def default_pypeit_par(self):
        """
        Set default parameters for the reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'mmt_mmirs'


        # Image processing steps
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False, use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)
        #par['calibrations']['traceframe']['process']['use_darkimage'] = True
        #par['calibrations']['pixelflatframe']['process']['use_darkimage'] = True
        #par['calibrations']['illumflatframe']['process']['use_darkimage'] = True
        #par['scienceframe']['process']['use_darkimage'] = True
        par['scienceframe']['process']['use_illumflat'] = True

        # Wavelengths
        # 1D wavelength solution with arc lines
        par['calibrations']['wavelengths']['rms_threshold'] = 0.5
        par['calibrations']['wavelengths']['sigdetect']=5
        par['calibrations']['wavelengths']['fwhm'] = 5
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

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 60]
        par['calibrations']['tiltframe']['exprng'] = [60, None]
        par['calibrations']['arcframe']['exprng'] = [60, None]
        par['calibrations']['darkframe']['exprng'] = [30, None]
        par['scienceframe']['exprng'] = [30, None]

        # dark
        par['calibrations']['darkframe']['process']['apply_gain'] = True

        # cosmic ray rejection
        par['scienceframe']['process']['sigclip'] = 5.0
        par['scienceframe']['process']['objlim'] = 2.0
        par['scienceframe']['process']['grow'] = 0.5

        # Science reduction
        par['reduce']['findobj']['sig_thresh'] = 5.0
        par['reduce']['skysub']['sky_sigrej'] = 5.0
        par['reduce']['findobj']['find_trim_edge'] = [5,5]
        # Do not correct for flexure
        par['flexure']['spec_method'] = 'skip'

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 8
        # ToDo: replace the telluric grid file for MMT site.
        par['sensfunc']['IR']['telgridfile'] = resource_filename('pypeit', '/data/telluric/TelFit_MaunaKea_3100_26100_R20000.fits')

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

        if (self.get_meta_value(scifile, 'dispname')=='HK') and (self.get_meta_value(scifile, 'dichroic')=='zJ'):
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'mmt_mmirs_HK_zJ.fits'
        elif (self.get_meta_value(scifile, 'dispname')=='K3000') and (self.get_meta_value(scifile, 'dichroic')=='Kspec'):
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'mmt_mmirs_K3000_Kspec.fits'
        elif (self.get_meta_value(scifile, 'dispname')=='J') and (self.get_meta_value(scifile, 'dichroic')=='zJ'):
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'mmt_mmirs_J_zJ.fits'

        return par

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
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
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def bpm(self, filename, det, shape=None, msbias=None):
        """ Generate a BPM

        Parameters
        ----------
        det : int, REQUIRED
        **null_kwargs:
           Captured and never used

        Returns
        -------
        badpix : ndarray

        """
        # Get the empty bpm: force is always True
        bpm_img = self.empty_bpm(filename, det, shape=shape)

        # Fill in bad pixels if a master bias frame is provided
        if msbias is not None:
            return self.bpm_frombias(msbias, det, bpm_img)

        msgs.info("Using hard-coded BPM for det=1 on MMIRS")

        # Get the binning
        hdu = fits.open(filename)
        binning = hdu[1].header['CCDSUM']
        hdu.close()

        # Apply the mask
        xbin, ybin = int(binning.split(' ')[0]), int(binning.split(' ')[1])
        bpm_img[:, 187 // ybin] = 1

        return bpm_img

    def get_rawimage(self, raw_file, det):
        """
        Load up the raw image and generate a few other bits and pieces
        that are key for image processing

        Args:
            raw_file (str):
            det (int):

        Returns:
            tuple:
                raw_img (np.ndarray) -- Raw image for this detector
                hdu (astropy.io.fits.HDUList)
                exptime (float)
                rawdatasec_img (np.ndarray)
                oscansec_img (np.ndarray)

        """
        # Check for file; allow for extra .gz, etc. suffix
        fil = glob.glob(raw_file + '*')
        if len(fil) != 1:
            msgs.error("Found {:d} files matching {:s}".format(len(fil)))

        # Read
        msgs.info("Reading MMIRS file: {:s}".format(fil[0]))
        hdu = fits.open(fil[0])
        head1 = fits.getheader(fil[0],1)

        detector_par = self.get_detector_par(hdu, det if det is None else 1)

        # get the x and y binning factors...
        binning = head1['CCDSUM']
        xbin, ybin = [int(ibin) for ibin in binning.split(' ')]

        # First read over the header info to determine the size of the output array...
        datasec = head1['DATASEC']
        x1, x2, y1, y2 = np.array(parse.load_sections(datasec, fmt_iraf=False)).flatten()

        # ToDo: I am currently using the standard double correlated frame, that is a difference between
        # the first and final read-outs. In the future need to explore up-the-ramp fitting.
        if len(hdu)>2:
            data = mmirs_read_amp(hdu[1].data.astype('float64')) - mmirs_read_amp(hdu[2].data.astype('float64'))
        else:
            data = mmirs_read_amp(hdu[1].data.astype('float64'))
        array = data[x1-1:x2,y1-1:y2]

        ## ToDo: This is a hack. Need to solve this issue. I cut at 998 due to the HK zero order contaminating
        ## the blue part of the zJ+HK spectrum. For other setup, you do not need to cut the detector.
        if (head1['FILTER']=='zJ') and (head1['DISPERSE']=='HK'):
            array = array[:int(998/ybin),:]
        rawdatasec_img = np.ones_like(array,dtype='int')
        oscansec_img = np.ones_like(array,dtype='int')

        # Need the exposure time
        exptime = hdu[self.meta['exptime']['ext']].header[self.meta['exptime']['card']]
        # Return, transposing array back to orient the overscan properly
        return detector_par, np.flipud(array), hdu, exptime, np.flipud(rawdatasec_img),\
               np.flipud(np.flipud(oscansec_img))

def mmirs_read_amp(img, namps=32):
    '''
    MMIRS has 32 reading out channels. Need to deal with this issue a little bit. I am not using the pypeit overscan
    subtraction since we need to do the up-the-ramp fitting in the future.
    Imported from MMIRS IDL pipeline refpix.pro
    '''

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
        ref1, med1, std1 = sigma_clipped_stats(img_out_ref[:, amp * ampsize + 2 * np.arange(int(ampsize / 2))],
                                               sigma=3)
        ref2, med2, std2 = sigma_clipped_stats(img_out_ref[:, amp * ampsize + 2 * np.arange(int(ampsize / 2)) + 1],
                                               sigma=3)
        ref12 = (ref1 + ref2) / 2.
        img_out[:, amp * ampsize:(amp + 1) * ampsize] -= ref12

    return img_out