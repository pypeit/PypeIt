""" Module for Keck/NIRSPEC specific codes
"""
from __future__ import absolute_import, division, print_function

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
                            dispaxis        = 0,
                            dispflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.193,
                            darkcurr        = 0.8,
                            saturation      = 65535.,
                            nonlinear       = 0.76,
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

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for NIRSPEC reductions
        """
        par = pypeitpar.PypeItPar()
        # Frame numbers
        par['calibrations']['standardframe']['number'] = 1
        par['calibrations']['biasframe']['number'] = 0
        par['calibrations']['pixelflatframe']['number'] = 5
        par['calibrations']['traceframe']['number'] = 5
        par['calibrations']['arcframe']['number'] = 0
        # Set wave tilts order
        par['calibrations']['tilts']['order'] = 2
        # Scienceimage default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()
        # Do not flux calibrate
        # NIRSPEC uses sky lines to wavelength calibrate; no need for flexure correction
        par['flexure'] = pypeitpar.FlexurePar()
        par['flexure']['method'] = 'skip'
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 2]
        par['calibrations']['darkframe']['exprng'] = [None, 5]
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['pixelflatframe']['exprng'] = [0, None]
        par['calibrations']['traceframe']['exprng'] = [0, None]
        par['calibrations']['standardframe']['exprng'] = [None,5]
        par['scienceframe']['exprng'] = [1, None]
        # Lower the default threshold for tilts
        par['calibrations']['tilts']['tracethresh'] = 10.
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20  # Good for NIRSPEC-1
        par['calibrations']['wavelengths']['lowest_nsig'] = 5.      # Good for NIRSPEC-1
        par['calibrations']['wavelengths']['min_nsig'] = 5.      # Good for NIRSPEC-1

        return par

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

    def metadata_keys(self):
        return super(KeckNIRSPECSpectrograph, self).metadata_keys() \
                    + ['echellepos', 'crosspos', 'idname']

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
        if ftype == 'arc':
            return good_exp & self.lamps(fitstbl, 'arcs') & (fitstbl['hatch'] == 1) \
                        & (fitstbl['idname'] == 'arclamp')

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

    def get_match_criteria(self):
        """Set the general matching criteria for Keck NIRSPEC."""
        match_criteria = {}
        for key in framematch.FrameTypeBitMask().keys():
            match_criteria[key] = {}

        match_criteria['standard']['match'] = {}
        match_criteria['standard']['match']['naxis0'] = '=0'
        match_criteria['standard']['match']['naxis1'] = '=0'

        match_criteria['bias']['match'] = {}
        match_criteria['bias']['match']['naxis0'] = '=0'
        match_criteria['bias']['match']['naxis1'] = '=0'

        match_criteria['pixelflat']['match'] = {}
        match_criteria['pixelflat']['match']['naxis0'] = '=0'
        match_criteria['pixelflat']['match']['naxis1'] = '=0'
        match_criteria['pixelflat']['match']['decker'] = ''
        match_criteria['pixelflat']['match']['dispname'] = ''

        match_criteria['trace']['match'] = {}
        match_criteria['trace']['match']['naxis0'] = '=0'
        match_criteria['trace']['match']['naxis1'] = '=0'
        match_criteria['trace']['match']['decker'] = ''
        match_criteria['trace']['match']['dispname'] = ''

        match_criteria['arc']['match'] = {}
        match_criteria['arc']['match']['naxis0'] = '=0'
        match_criteria['arc']['match']['naxis1'] = '=0'
        match_criteria['arc']['match']['dispname'] = ''

        return match_criteria

    # TODO: This function is unstable to shape...
    def bpm(self, shape=None, **null_kwargs):
        """ Generate a BPM

        Parameters
        ----------
        shape : tuple, REQUIRED

        Returns
        -------
        badpix : ndarray

        """
        if shape is None:
            raise ValueError('Must provide shape for Keck NIRSPEC bpm.')
        # Edges of the detector are junk
        msgs.info("Custom bad pixel mask for NIRSPEC")
        self.bpm_img = np.zeros(shape, dtype=np.int8)
        self.bpm_img[:, :20] = 1.
        self.bpm_img[:, 1000:] = 1.

        return self.bpm_img

class KeckNIRSPECLowSpectrograph(KeckNIRSPECSpectrograph):
    """
    Child to handle NIRSPEC low-dispersion specific code
    """

    def __init__(self):
        # Get it started
        super(KeckNIRSPECLowSpectrograph, self).__init__()
        self.spectrograph = 'keck_nirspec_low'


    def check_header(self, headers):
        """Validate elements of the header."""
        chk_dict = {}
        # chk_dict is 1-indexed!
        chk_dict[1] = {}
        # THIS CHECK IS A MUST! It performs a standard check to make sure the data are 2D.
        chk_dict[1]['NAXIS'] = 2
        return chk_dict

    '''
    def header_keys(self):
        """
        Header keys specific to keck_nirspec

        Returns:

        """
        head_keys = self.nirspec_header_keys()
        # Add the name of the filter used
        head_keys[0]['filter'] = 'FILNAME'
        return head_keys
    '''

    def setup_arcparam(self, arcparam, disperser=None, fitstbl=None, arc_idx=None, msarc_shape=None,
                       binspectral=None, **null_kwargs):
        """

        Args:
            arcparam:
            disperser:
            fitstbl:
            arc_idx:
            msarc_shape:
            **null_kwargs:

        Returns:

        """
        arcparam['lamps'] = ['OH_R24000']
        arcparam['nonlinear_counts'] = self.detector[0]['nonlinear']*self.detector[0]['saturation']
        arcparam['min_ampl'] = 1000.       # Minimum amplitude

        if fitstbl['filter1'][arc_idx] == 'NIRSPEC-1':
            arcparam['n_first'] = 2  # Too much curvature for 1st order
            arcparam['disp'] = 2.1093  # Ang per pixel for Low-Res, NIRSPEC-1 filter
            arcparam['b1'] = 1. / arcparam['disp'] / msarc_shape[0]
            arcparam['wvmnx'][0] = 9400.  # Min wavelength
            arcparam['wvmnx'][1] = 11300.  # Max wavelength
            arcparam['wv_cen'] = 10000.  # Central wavelength

