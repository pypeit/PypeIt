""" Module for Keck/NIRSPEC specific codes
"""
import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph

from IPython import embed

class KeckMOSFIRESpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/NIRSPEC specific code
    """
    def __init__(self):
        # Get it started
        super(KeckMOSFIRESpectrograph, self).__init__()
        self.telescope = telescopes.KeckTelescopePar()
        self.camera = 'MOSFIRE'
        self.detector = [
                # Detector 1
            pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 1,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.193,
                            darkcurr        = 0.8,
                            saturation      = 33000., # DN
                            nonlinear       = 1.00,  # docs say linear to 90,000 but our flats are usually higher
                            numamplifiers   = 1,
                            gain            = 5.65685,  # Taken from a random Header
                            ronoise         = 23,
                            datasec         = '[:,:]',
                            oscansec        = '[:,:]'
                            )]
        self.numhead = 1

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for NIRSPEC reductions
        """
        par = pypeitpar.PypeItPar()
        # No bias subtraction
        par['calibrations']['biasframe']['useframe'] = 'none'
        # Scienceimage default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()
        # Do not flux calibrate
        # NIRSPEC uses sky lines to wavelength calibrate; no need for flexure correction
        par['flexure'] = pypeitpar.FlexurePar()
        par['flexure']['method'] = 'skip'
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['arcframe']['exprng'] = [1, None]
        par['calibrations']['biasframe']['exprng'] = [None, 2]
        par['calibrations']['darkframe']['exprng'] = [None, 5]
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['pixelflatframe']['exprng'] = [0, None]
        par['calibrations']['traceframe']['exprng'] = [0, None]
        par['calibrations']['standardframe']['exprng'] = [None,5]
        par['scienceframe']['exprng'] = [1, None]
        # Lower the default threshold for tilts
        par['calibrations']['tilts']['tracethresh'] = 10.
        # Slits
        par['calibrations']['slits']['sigdetect'] = 200.
        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps']  = ['OH_R24000']
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20
        par['calibrations']['wavelengths']['sigdetect'] = 5.

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
        meta['decker'] = dict(ext=0, card='MASKNAME')
        meta['binning'] = dict(ext=0, card=None, default='1,1')

        meta['mjd'] = dict(ext=0, card='MJD-OBS')
        meta['exptime'] = dict(ext=0, card='TRUITIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        meta['dispname'] = dict(ext=0, card='GRATMODE')
        #meta['hatch'] = dict(ext=0, card='CALMPOS')
        meta['idname'] = dict(ext=0, card='KOAIMTYP')
        # Filter
        meta['filter1'] = dict(ext=0, card='FILTER')
        # Lamps
        lamp_names = ['FLATSPEC']
        for kk,lamp_name in enumerate(lamp_names):
            meta['lampstat{:02d}'.format(kk+1)] = dict(ext=0, card=lamp_name)
        # Ingest
        self.meta = meta

    def configuration_keys(self):
        return ['decker', 'dispname', 'filter1']

    def pypeit_file_keys(self):
        pypeit_keys = super(KeckMOSFIRESpectrograph, self).pypeit_file_keys()
        pypeit_keys += ['calib', 'comb_id', 'bkg_id']
        return pypeit_keys

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['science', 'standard']:
            return good_exp & (fitstbl['idname'] == 'object')
        if ftype in ['bias', 'dark']:
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['idname'] == 'dark')
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & self.lamps(fitstbl, 'dome') & (fitstbl['idname'] == 'flatlamp')
        if ftype == 'pinhole':
            # Don't type pinhole frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            # TODO: This is a kludge.  Allow science frames to also be
            # classified as arcs
            is_arc = self.lamps(fitstbl, 'arcs') & (fitstbl['idname'] == 'arclamp')
            is_obj = self.lamps(fitstbl, 'off') & (fitstbl['idname'] == 'object')
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
            return fitstbl['lampstat01'] == '1'

        raise ValueError('No implementation for status = {0}'.format(status))

    def bpm(self, shape=None, **null_kwargs):
        """ Generate a BPM
        Parameters
        ----------
        shape : tuple, REQUIRED
        Returns
        -------
        badpix : ndarray
        """
        # Edges of the detector are junk
        msgs.info("Custom bad pixel mask for MOSFIRE")
        self.bpm_img = np.zeros(shape, dtype=np.int8)
        #self.bpm_img[:, :20] = 1.
        #self.bpm_img[:, 1000:] = 1.

        return self.bpm_img
