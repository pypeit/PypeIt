""" Module for Keck/NIRSPEC specific codes
"""
import numpy as np

from pypeit import msgs
from pypeit.images import detector_container
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pkg_resources import resource_filename

from pypeit import debugger

class KeckNIRSPECSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/NIRSPEC specific code
    """
    ndet = 1

    def __init__(self):
        # Get it started
        super(KeckNIRSPECSpectrograph, self).__init__()
        self.telescope = telescopes.KeckTelescopePar()
        self.camera = 'NIRSPEC'

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
        detector_dict = dict(
            det=1,
            binning         ='1,1',  # No binning allowed
            dataext         = 0,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.193,
            darkcurr        = 0.8,
            saturation      = 100000.,
            nonlinear       = 1.00,  # docs say linear to 90,000 but our flats are usually higher
            numamplifiers   = 1,
            mincounts       = -1e10,
            gain            = np.atleast_1d(5.8),
            ronoise         = np.atleast_1d(23.),
            datasec         = np.atleast_1d('[:,:]'),
            oscansec        = np.atleast_1d('[:,:]')
            )
        return detector_container.DetectorContainer(**detector_dict)

    def default_pypeit_par(self):
        """
        Set default parameters for Keck/NIRSPEC
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'keck_nirspec_low'
        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20 #0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect']=5.0
        par['calibrations']['wavelengths']['fwhm']= 5.0
        par['calibrations']['wavelengths']['n_final']= 4
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'
        # Reidentification parameters
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_nires.fits'
        par['calibrations']['slitedges']['edge_thresh'] = 200.
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Flats
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.80

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0

        # Flexure
        par['flexure']['spec_method'] = 'skip'

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] ='nothing'


        # Overscan but not bias
        #  This seems like a kludge of sorts
        par['calibrations']['biasframe']['useframe'] = 'none'
        # No overscan
        par['scienceframe']['process']['overscan'] ='none'
        for key in par['calibrations'].keys():
            if 'frame' in key:
                par['calibrations'][key]['process']['overscan'] = 'none'


        # The settings below enable NIRSPEC dark subtraction from the traceframe and pixelflatframe, but enforce
        # that this bias won't be subtracted from other images. It is a hack for now, because eventually we want to
        # perform this operation with the dark frame class, and we want to attach individual sets of darks to specific
        # images.
        par['calibrations']['biasframe']['useframe'] = 'bias'
        par['calibrations']['traceframe']['process']['bias'] = 'force'
        par['calibrations']['pixelflatframe']['process']['bias'] = 'force'
        par['calibrations']['arcframe']['process']['bias'] = 'skip'
        par['calibrations']['tiltframe']['process']['bias'] = 'skip'
        par['calibrations']['standardframe']['process']['bias'] = 'skip'
        par['scienceframe']['process']['bias'] = 'skip'


        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 20]
        par['calibrations']['arcframe']['exprng'] = [20, None]
        par['calibrations']['darkframe']['exprng'] = [20, None]
        par['scienceframe']['exprng'] = [20, None]

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 8
        par['sensfunc']['IR']['telgridfile'] = resource_filename('pypeit', '/data/telluric/TelFit_MaunaKea_3100_26100_R20000.fits')
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
        meta['binning'] = dict(ext=0, card=None, default='1,1')

        meta['mjd'] = dict(ext=0, card='MJD-OBS')
        meta['exptime'] = dict(ext=0, card='ELAPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        meta['dispname'] = dict(ext=0, card='DISPERS')
        meta['hatch'] = dict(ext=0, card='CALMPOS')
        meta['idname'] = dict(ext=0, card='IMAGETYP')
        # Lamps
        lamp_names = ['NEON', 'ARGON', 'KRYPTON', 'XENON', 'ETALON', 'FLAT']
        for kk,lamp_name in enumerate(lamp_names):
            meta['lampstat{:02d}'.format(kk+1)] = dict(ext=0, card=lamp_name)
        # Ingest
        self.meta = meta

    def configuration_keys(self):
        return ['decker', 'dispname']

    def pypeit_file_keys(self):
        pypeit_keys = super(KeckNIRSPECSpectrograph, self).pypeit_file_keys()
        pypeit_keys += ['calib', 'comb_id', 'bkg_id']
        return pypeit_keys

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['science', 'standard']:
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == 0) \
                        & (fitstbl['idname'] == 'object')
        if ftype in ['bias', 'dark']:
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == 0) \
                        & (fitstbl['idname'] == 'dark')
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & self.lamps(fitstbl, 'dome') & (fitstbl['hatch'] == 1) \
                        & (fitstbl['idname'] == 'flatlamp')
        if ftype == 'pinhole':
            # Don't type pinhole frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            # TODO: This is a kludge.  Allow science frames to also be
            # classified as arcs
            is_arc = self.lamps(fitstbl, 'arcs') & (fitstbl['hatch'] == 1) \
                            & (fitstbl['idname'] == 'arclamp')
            is_obj = self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == 0) \
                        & (fitstbl['idname'] == 'object')
            return good_exp & (is_arc | is_obj)
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
            return np.all(np.array([fitstbl[k] == 0 for k in fitstbl.keys() if 'lampstat' in k]),
                          axis=0)
        if status == 'arcs':
            # Check if any arc lamps are on
            arc_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(1,6) ]
            return np.any(np.array([ fitstbl[k] == 1 for k in fitstbl.keys()
                                            if k in arc_lamp_stat]), axis=0)
        if status == 'dome':
            return fitstbl['lampstat06'] == 1

        raise ValueError('No implementation for status = {0}'.format(status))

    # TODO: This function is unstable to shape...
    def bpm(self, filename, det, shape=None, msbias=None):
        """ Generate a BPM
        Parameters
        ----------
        shape : tuple, REQUIRED
        msbias : numpy.ndarray, required if the user wishes to generate a BPM based on a master bias

        Returns
        -------
        badpix : ndarray
        """
        bpm_img = self.empty_bpm(filename, det)

        # Fill in bad pixels if a master bias frame is provided
        if msbias is not None:
            return self.bpm_frombias(msbias, det, bpm_img)

        # Edges of the detector are junk
        msgs.info("Custom bad pixel mask for NIRSPEC")
        bpm_img[:, :20] = 1.
        bpm_img[:, 1000:] = 1.

        return bpm_img

class KeckNIRSPECLowSpectrograph(KeckNIRSPECSpectrograph):
    """
    Child to handle NIRSPEC low-dispersion specific code
    """

    def __init__(self):
        # Get it started
        super(KeckNIRSPECLowSpectrograph, self).__init__()
        self.spectrograph = 'keck_nirspec_low'


    @property
    def telluric_grid_file(self):
        """Return the grid of HITRAN atmosphere models for telluric correctinos"""
        return resource_filename('pypeit', '/data/telluric/TelFit_MaunaKea_3100_26100_R20000.fits')
