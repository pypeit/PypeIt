"""
Module for VLT/SINFONI specific methods.

.. include:: ../include/links.rst
"""

from IPython import embed

import numpy as np
from astropy.io import fits
from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container

class VLTSINFONISpectrograph(spectrograph.Spectrograph):
    """
    Child to handle VLT/SINFONI specific code
    """
    ndet = 1
    name = 'vlt_sinfoni'
    telescope = telescopes.VLTTelescopePar()
    camera = 'SINFONI'
    url = 'https://www.eso.org/sci/facilities/paranal/decommissioned/sinfoni.html'
    header_name = 'SINFONI'
    supported = True
    comment = 'Gratings tested: K'

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
            binning         = '1,1',
            det             = 1,
            dataext         = 0,
            specaxis        = 0,
            specflip        = True,
            spatflip        = False,
            platescale      = 0.0125,
            darkcurr        = 0.15,
            saturation      = 1e9, # ADU, this is hacked for now
            nonlinear       = 1.00,  # docs say linear to 90,000 but our flats are usually higher
            numamplifiers   = 1,
            mincounts       = -1e10,
            gain            = np.atleast_1d(2.42),
            ronoise         = np.atleast_1d(7.0),
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

        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.30 #0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect']=5.0
        par['calibrations']['wavelengths']['fwhm']= 5.0
        par['calibrations']['wavelengths']['n_final']= 4
        par['calibrations']['wavelengths']['lamps'] = ['OH_FIRE_Echelle']
        #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        #par['calibrations']['wavelengths']['method'] = 'holy-grail'
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_sinfoni_K.fits'
        par['calibrations']['wavelengths']['nsnippet'] = 1

        # Reidentification parameters
        par['calibrations']['slitedges']['edge_thresh'] = 50.
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['rm_slits'] = '1:1024:983' # Remove the center slit that is not illuminated

        # Tilts parameters
        par['calibrations']['tilts']['tracethresh'] = 5.0

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 20]
        par['calibrations']['arcframe']['exprng'] = [20, None]
        par['calibrations']['darkframe']['exprng'] = [20, None]
        par['scienceframe']['exprng'] = [20, None]


        # TODO: We need to implement dark subtraction for the arcframe and
        # tiltframe. Currently the pypeit file won't let me do this.
        par['calibrations']['arcframe']['process']['sigclip'] = 20.0
        #par['calibrations']['arcframe']['process']['combine'] = 'median'
        par['calibrations']['arcframe']['process']['mask_cr'] = True


        par['calibrations']['tiltframe']['process']['sigclip'] = 20.0
        #par['calibrations']['tiltframe']['process']['combine'] = 'median'
        par['calibrations']['tiltframe']['process']['mask_cr'] = True

        par['calibrations']['skyframe']['process']['sigclip'] = 20.0
        #par['calibrations']['skyframe']['process']['combine'] = 'median'
        par['calibrations']['skyframe']['process']['mask_cr'] = True


        # Flats
        turn_off = dict(use_biasimage=False, use_overscan=False, use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # # Require dark images to be subtracted from the flat images used for tracing, pixelflats, and illumflats
        # par['calibrations']['pixelflatframe']['process']['use_darkimage'] = True
        # par['calibrations']['illumflatframe']['process']['use_darkimage'] = True
        # par['calibrations']['traceframe']['process']['use_darkimage'] = True
        # TODO: `mask_cr` now defaults to True for darks.  Should this be turned off?

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.9
        par['reduce']['extraction']['sn_gauss'] = 5.0
        par['reduce']['extraction']['model_full_slit'] = True  # local sky subtraction operates on entire slit

        # Object finding
        par['reduce']['findobj']['find_fwhm'] = 10
        par['reduce']['findobj']['skip_second_find'] = True


        # Sky subtraction
        par['reduce']['skysub']['global_sky_std']  = False # Do not perform global sky subtraction for standard stars

        # Flexure
        par['flexure']['spec_method'] = 'skip'


        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] ='nothing'

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 20]
        par['calibrations']['arcframe']['exprng'] = [20, None]
        par['calibrations']['darkframe']['exprng'] = [20, None]
        par['scienceframe']['exprng'] = [20, None]

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 7
        par['sensfunc']['IR']['telgridfile'] = 'TelFit_Paranal_NIR_9800_25000_R25000.fits'

        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA', required_ftypes=['science', 'standard'])  # Need to convert to : separated
        self.meta['dec'] = dict(ext=0, card='DEC', required_ftypes=['science', 'standard'])
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='HIERARCH ESO TEL AIRM START', required_ftypes=['science', 'standard'])
        # Extras for config and frametyping
        self.meta['decker'] = dict(ext=0, card='HIERARCH ESO INS OPTI1 NAME')
        self.meta['filter1'] = dict(ext=0, card='HIERARCH ESO INS FILT1 NAME')
        self.meta['dispname'] = dict(ext=0, card='HIERARCH ESO INS GRAT1 NAME')
        self.meta['idname'] = dict(ext=0, card='HIERARCH ESO OCS DET IMGNAME')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')
        # self.meta['idname'] = dict(ext=0, card='HIERARCH ESO DPR CATG')
        # Dithering
        self.meta['dithoff'] = dict(ext=0, card='HIERARCH ESO SEQ CUMOFFSETY',
                                   required_ftypes=['science', 'standard'])

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
        if meta_key == 'decker':
            try:  # Science
                decker = headarr[0]['HIERARCH ESO INS SLIT NAME']
            except KeyError:  # Standard!
                try:
                    decker = headarr[0]['HIERARCH ESO SEQ SPEC TARG']
                except KeyError:
                    return None
            return decker
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
        return ['decker', 'dispname', 'filter1']

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
        return ['HIERARCH ESO INS OPTI1 NAME', 'HIERARCH ESO INS GRAT1 NAME',
                'HIERARCH ESO INS FILT1 NAME']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        pypeit_keys = super().pypeit_file_keys()
        # TODO: Why are these added here? See
        # pypeit.metadata.PypeItMetaData.set_pypeit_cols
        pypeit_keys += ['calib', 'comb_id', 'bkg_id']
        return pypeit_keys + ['dithoff']



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
        # TODO: Allow for 'sky' frame type, for now include sky in
        # 'science' category
        if ftype == 'science':
            return good_exp & ((fitstbl['idname'] == 'SINFONI_IFS_OBS')
                                | (fitstbl['target'] == 'STD,TELLURIC')
                                | (fitstbl['target'] == 'SKY,STD'))
        if ftype == 'standard':
            return good_exp & ((fitstbl['target'] == 'STD') | (fitstbl['target'] == 'SKY,STD'))
        #if ftype == 'bias':
        #    return good_exp & (fitstbl['target'] == 'BIAS')
        if ftype == 'dark':
            return good_exp & (fitstbl['target'] == 'DARK')
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & (fitstbl['target'] == 'FLAT,LAMP')
        #if ftype == 'pinhole':
        #    # Don't type pinhole
        #    return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & ((fitstbl['target'] == 'WAVE,LAMP') | (fitstbl['idname'] == 'SINFONI_IFS_OBS') |
                               (fitstbl['idname'] == 'SINFONI_IFS_SKY'))
        # Putting this in now in anticipation of the sky class
        if ftype in ['sky']:
            return good_exp & (fitstbl['idname'] == 'SINFONI_IFS_SKY')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


    def lamps(self, fitstbl, status):
        """
        Check the lamp status.

        Args:
            fitstbl (`astropy.table.Table`_):
                The table with the fits header meta data.
            status (:obj:`str`):
                The status to check. Can be ``'off'``, ``'arcs'``, or
                ``'dome'``.

        Returns:
            `numpy.ndarray`_: A boolean array selecting fits files that meet
            the selected lamp status.

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
            return fitstbl['lampstat01'] == '1'

        raise ValueError('No implementation for status = {0}'.format(status))



