""" Module for Keck/NIRSPEC specific codes
"""
import numpy as np

from pypeit import msgs
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph

from pypeit import debugger

class KeckNIRSPECSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/NIRSPEC specific code
    """
    def __init__(self):
        # Get it started
        super(KeckNIRSPECSpectrograph, self).__init__()
        self.telescope = telescopes.KeckTelescopePar()
        self.camera = 'NIRSPEC'
        self.detector = [
                # Detector 1
            pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.193,
                            darkcurr        = 0.8,
                            saturation      = 100000.,
                            nonlinear       = 1.00,  # docs say linear to 90,000 but our flats are usually higher
                            numamplifiers   = 1,
                            gain            = 5.8,
                            ronoise         = 23,
                            datasec         = '[:,:]',
                            oscansec        = '[:,:]'
                            )]
        self.numhead = 1
        # Uses default timeunit
        # Uses default primary_hdrext
        # self.sky_file = ?


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
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'
        # Reidentification parameters
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_nires.fits'
        par['calibrations']['slitedges']['edge_thresh'] = 200.
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Flats
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.80
        par['calibrations']['flatfield']['illumflatten'] = True

        # Extraction
        par['scienceimage']['skysub']['bspline_spacing'] = 0.8
        par['scienceimage']['extraction']['sn_gauss'] = 4.0

        # Flexure
        par['flexure']['method'] = 'skip'

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
        return par

    # JFH Replaced with updated values based on experienced with MOSFIRE above.
    #
    # @staticmethod
    # def default_pypeit_par():
    #     """
    #     Set default parameters for NIRSPEC reductions
    #     """
    #     par = pypeitpar.PypeItPar()
    #     # Frame numbers
    #     par['calibrations']['standardframe']['number'] = 1
    #     par['calibrations']['biasframe']['number'] = 0
    #     par['calibrations']['pixelflatframe']['number'] = 5
    #     par['calibrations']['traceframe']['number'] = 5
    #     par['calibrations']['arcframe']['number'] = 1
    #     # Do not flux calibrate
    #     # NIRSPEC uses sky lines to wavelength calibrate; no need for flexure correction
    #     par['flexure'] = pypeitpar.FlexurePar()
    #     par['flexure']['method'] = 'skip'
    #     # Set the default exposure time ranges for the frame typing
    #     par['calibrations']['arcframe']['exprng'] = [1, None]
    #     par['calibrations']['biasframe']['exprng'] = [None, 2]
    #     par['calibrations']['darkframe']['exprng'] = [None, 5]
    #     par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
    #     par['calibrations']['pixelflatframe']['exprng'] = [0, None]
    #     par['calibrations']['traceframe']['exprng'] = [0, None]
    #     par['calibrations']['standardframe']['exprng'] = [None,5]
    #     par['scienceframe']['exprng'] = [1, None]
    #     # Lower the default threshold for tilts
    #     par['calibrations']['tilts']['tracethresh'] = 10.
    #     # Slits
    #     par['calibrations']['slitedges']['edge_thresh'] = 200.
    #     # 1D wavelength solution
    #     par['calibrations']['wavelengths']['lamps']  = ['OH_R24000']
    #     par['calibrations']['wavelengths']['rms_threshold'] = 0.20  # Good for NIRSPEC-1
    #     par['calibrations']['wavelengths']['sigdetect'] = 5.      # Good for NIRSPEC-1
    #
    #     return par
    #

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
        expected_values = { '0.INSTRUME': 'NIRSPEC',
                               '0.NAXIS': 2,
                              '0.NAXIS1': 1024,
                              '0.NAXIS2': 1024 }
        super(KeckNIRSPECSpectrograph, self).check_headers(headers,
                                                           expected_values=expected_values)

    def header_keys(self):
        """
        Return a dictionary with the header keywords to read from the
        fits file.
        Returns:
            dict: A nested dictionary with the header keywords to read.
            The first level gives the extension to read and the second
            level gives the common name for header values that is passed
            on to the PypeItMetaData object.
        """
        hdr_keys = {}
        hdr_keys[0] = {}
        hdr_keys[0]['idname'] = 'IMAGETYP'
        #hdr_keys[0]['date'] = 'DATE-OBS'
        hdr_keys[0]['utc'] = 'UTC'
        hdr_keys[0]['target'] = 'OBJECT'
        hdr_keys[0]['time'] = 'MJD-OBS'
        hdr_keys[0]['ra'] = 'RA'
        hdr_keys[0]['dec'] = 'DEC'
        hdr_keys[0]['airmass'] = 'AIRMASS'
        hdr_keys[0]['decker'] = 'SLITNAME'
        hdr_keys[0]['echellepos'] = 'ECHLPOS'
        hdr_keys[0]['crosspos'] = 'DISPPOS'
        hdr_keys[0]['naxis0'] = 'NAXIS2'
        hdr_keys[0]['naxis1'] = 'NAXIS1'
        hdr_keys[0]['filter1'] = 'FILNAME'
        hdr_keys[0]['dispname'] = 'DISPERS'
        hdr_keys[0]['hatch'] = 'CALMPOS'
        hdr_keys[0]['slitwid'] = 'SLITWIDT'
        hdr_keys[0]['slitlen'] = 'SLITLEN'

        # 'ELAPTIME' is added by KOA, but otherwise would need to do 'ITIME' * 'COADDS'
        hdr_keys[0]['exptime'] = 'ELAPTIME'

        # Lamp names and statuses
        lamp_names = ['NEON', 'ARGON', 'KRYPTON', 'XENON', 'ETALON', 'FLAT']
        for kk,lamp_name in enumerate(lamp_names):
            hdr_keys[0]['lampstat{:02d}'.format(kk+1)] = lamp_name

        return hdr_keys
    '''
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

    '''
    def metadata_keys(self):
        return super(KeckNIRSPECSpectrograph, self).metadata_keys() \
                    + ['echellepos', 'crosspos', 'idname']
    '''
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
    def bpm(self, filename, det, shape=None):
        """ Generate a BPM
        Parameters
        ----------
        shape : tuple, REQUIRED
        Returns
        -------
        badpix : ndarray
        """
        bpm_img = self.empty_bpm(filename, det)
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


    '''
    def check_header(self, headers):
        """Validate elements of the header."""
        chk_dict = {}
        # chk_dict is 1-indexed!
        chk_dict[1] = {}
        # THIS CHECK IS A MUST! It performs a standard check to make sure the data are 2D.
        chk_dict[1]['NAXIS'] = 2
        return chk_dict
    '''


