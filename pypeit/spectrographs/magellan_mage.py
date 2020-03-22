""" Module for Magellan/MAGE specific codes
"""
import numpy as np

from astropy.time import Time
from astropy.io import fits

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.core import parse
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.core import pixels
from pypeit import debugger

from IPython import embed

class MagellanMAGESpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Magellan/MAGE specific code
    """
    def __init__(self):
        # Get it started
        super(MagellanMAGESpectrograph, self).__init__()
        self.spectrograph = 'magellan_mage'
        self.camera = 'MagE'
        self.telescope = telescopes.MagellanTelescopePar()
        self.numhead = 1
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 1,
                            specflip        = True,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            # plate scale in arcsec/pixel
                            platescale      = 0.3,
                            # electrons/pixel/hour. From: http://www.lco.cl/telescopes-information/magellan/instruments/mage/the-mage-spectrograph-user-manual
                            darkcurr        = 1.00,
                            saturation      = 65535.,
                            # CCD is linear to better than 0.5 per cent up to digital saturation (65,536 DN including bias) in the Fast readout mode.
                            nonlinear       = 0.99,
                            numamplifiers   = 1,
                            gain            = 1.02, # depends on the readout
                            ronoise         = 2.9, # depends on the readout
                            datasec         = '[1:1024, 1:2048]',
                            oscansec        = '[1:1024, 2049:2176]',
                            )]
        # Taken from the MASE paper: https://arxiv.org/pdf/0910.1834.pdf
        #self.norders = 15
        # 20-6
        # Uses default timeunit
        # Uses default primary_hdrext
        # self.sky_file = ?

    @property
    def pypeline(self):
        return 'Echelle'

    def default_pypeit_par(self):
        """
        Set default parameters for magellan MagE reduction.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'magellan_mage'
        # Bias
        #par['calibrations']['biasframe']['useframe'] = 'overscan'
        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['lamps'] = ['ThAr_MagE']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']

        par['calibrations']['wavelengths']['method'] = 'reidentify'
        par['calibrations']['wavelengths']['cc_thresh'] = 0.50
        par['calibrations']['wavelengths']['cc_local_thresh'] = 0.50

        # Reidentification parameters
        par['calibrations']['wavelengths']['reid_arxiv'] = 'magellan_mage.fits'
        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        # Always correct for flexure, starting with default parameters
        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] = 'nothing'

        # Set slits and tilts parameters
        par['calibrations']['tilts']['tracethresh'] = [10]*self.norders
        par['calibrations']['slitedges']['fit_order'] = 5
        par['calibrations']['slitedges']['max_shift_adj'] = 3.
        par['calibrations']['slitedges']['edge_thresh'] = 10.  # Tough to get the bluest orders
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.3  # Allow for a short detected blue order
        # Find object parameters
        par['reduce']['findobj']['find_trim_edge'] = [4,4]    # Slit is too short to trim 5,5 especially with 2x binning
        # Always flux calibrate, starting with default parameters
        # Do not correct for flexure
        par['flexure']['method'] = 'skip'
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 20]
        par['calibrations']['arcframe']['exprng'] = [20, None]
        par['calibrations']['darkframe']['exprng'] = [20, None]
        par['scienceframe']['exprng'] = [20, None]
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
        meta['target'] = dict(ext=0, card='OBJECT')
        #TODO: Check decker is correct
        meta['decker'] = dict(ext=0, card='SLITNAME')
        meta['binning'] = dict(card=None, compound=True)
#        self.meta['binning'] = dict(ext=0, card='BINNING')
        meta['mjd'] = dict(ext=0, card=None, compound=True)
        meta['exptime'] = dict(ext=0, card='EXPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        meta['dispname'] = dict(ext=0, card='INSTRUME')
        meta['idname'] = dict(ext=0, card='EXPTYPE')

        # Ingest
        self.meta = meta

    def compound_meta(self, headarr, meta_key):
        """

        Args:
            headarr: list
            meta_key: str

        Returns:
            value

        """
        if meta_key == 'binning':
            binspatial, binspec = parse.parse_binning(headarr[0]['BINNING'])
            return parse.binning2string(binspec, binspatial)
        elif meta_key == 'mjd':
            time = '{:s}T{:s}'.format(headarr[0]['UT-DATE'], headarr[0]['UT-TIME'])
            ttime = Time(time, format='isot')
            return ttime.mjd
        else:
            msgs.error("Not ready for this compound meta")

    def configuration_keys(self):
        return []

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        if ftype in ['pinhole', 'dark']:
            # No pinhole or bias or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        elif ftype in ['bias']:
            return fitstbl['idname'] == 'Bias'
        elif ftype in ['pixelflat', 'trace']:
            return fitstbl['idname'] == 'Flat'
        elif ftype in ['arc']:
            return fitstbl['idname'] == 'ThAr-Lamp'
        else:
            return (fitstbl['idname'] == 'Object') \
                        & framematch.check_frame_exptime(fitstbl['exptime'], exprng)

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Override parent bpm function with BPM specific to X-Shooter VIS.

        .. todo::
            Allow for binning changes.

        Parameters
        ----------
        det : int, REQUIRED
        msbias : numpy.ndarray, required if the user wishes to generate a BPM based on a master bias
        **null_kwargs:
            Captured and never used

        Returns
        -------
        bpix : ndarray
          0 = ok; 1 = Mask

        """
        msgs.info("Custom bad pixel mask for MAGE")
        bpm_img = self.empty_bpm(filename, det, shape=shape)

        # Fill in bad pixels if a master bias frame is provided
        if msbias is not None:
            return self.bpm_frombias(msbias, det, bpm_img)

        # Get the binning
        hdu = fits.open(filename)
        binspatial, binspec = parse.parse_binning(hdu[0].header['BINNING'])
        hdu.close()
        # Do it
        bpm_img[:, :10//binspatial] = 1.  # Setting BPM on the edge of the detector often leads to false edges
        bpm_img[:, 1020//binspatial:] = 1.
        # Return
        return bpm_img

# TODO: Not sure if this was ever used.
#    @staticmethod
#    def slitmask(tslits_dict, pad=None, binning=None):
#        """
#         Generic routine ton construct a slitmask image from a tslits_dict. Children of this class can
#         overload this function to implement instrument specific slitmask behavior, for example setting
#         where the orders on an echelle spectrograph end
#
#         Parameters
#         -----------
#         tslits_dict: dict
#            Trace slits dictionary with slit boundary information
#
#         Optional Parameters
#         pad: int or float
#            Padding of the slit boundaries
#         binning: tuple
#            Spectrograph binning in spectral and spatial directions
#
#         Returns
#         -------
#         slitmask: ndarray int
#            Image with -1 where there are no slits/orders, and an integer where there are slits/order with the integer
#            indicating the slit number going from 0 to nslit-1 from left to right.
#
#         """
#
#        # These lines are always the same
#        pad = tslits_dict['pad'] if pad is None else pad
#        slitmask = pixels.slit_pixels(tslits_dict['lcen'], tslits_dict['rcen'], tslits_dict['nspat'], pad=pad)
#
#        return slitmask

    @property
    def norders(self):
        return 15   # 20-6

    @property
    def order_spat_pos(self):
        ord_spat_pos =  np.array([0.3157, 0.3986, 0.47465896, 0.5446689, 0.60911287, 0.66850584, 0.72341316,
               0.77448156, 0.82253604, 0.86875753, 0.91512689, 0.96524312])
        return ord_spat_pos

    @property
    def orders(self):
        return  np.arange(17, 5, -1, dtype=int)


    @property
    def spec_min_max(self):
        spec_max = np.full(self.norders, np.inf)
        spec_min = np.full(self.norders, -np.inf)
        return np.vstack((spec_min, spec_max))

    def order_platescale(self, order_vec, binning=None):
        """
        Returns the plate scale in arcseconds for each order

        Args:
            order_vec (np.ndarray): Order numbers
            binning (optional):

        Returns:
            np.ndarray: Platescale

        """
        norders = len(order_vec)
        binspatial, binspec = parse.parse_binning(binning)
        return np.full(norders, 0.30*binspatial)


