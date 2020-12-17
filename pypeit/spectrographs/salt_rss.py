"""
Spectrograph parameters module for the Robert Stobie Spectrograph  on SALT
"""
import numpy as np

from astropy.io import fits
from astropy.coordinates import Angle
from astropy import units as u
from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container

from typing import List
from pkg_resources import resource_filename

def flip_fits_slice(s: str) -> str:
    return '[' + ','.join(s.strip('[]').split(',')[::-1]) + ']'

class SALTRSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Robert Stobie Spectrograph on SALT specific code
    """
    ndet = 6

    def __init__(self):
        # Get it started
        super(SALTRSSpectrograph, self).__init__()
        self.spectrograph = 'salt_rss'
        self.telescope = telescopes.SALTelescopePar()

    def configuration_keys(self):
        """
        Return the metadata keys that defines a unique instrument
        configuration.

        This list is used by :class:`pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            list: List of keywords of data pulled from meta
        """
        #return ['gainset', 'readout_speed', 'decker', 'dispangle', 'dispname', 'binning', 'filter1']
        # return ['decker', 'dispangle', 'dispname', 'binning', 'filter1']
        #return ['decker', 'dispname', 'binning', 'filter1']
        return ['decker', 'dispname', 'binning']

    def init_meta(self):
        """
        Generate the meta data dict
        Note that the children can add to this

        Returns:
            self.meta: dict (generated in place)

        """
        meta = {}
        # Required (core)
        meta['ra'] = dict(ext=0, card='RA', required_ftypes=['science', 'standard'])
        meta['dec'] = dict(ext=0, card='DEC', required_ftypes=['science', 'standard'])
        meta['target'] = dict(ext=0, card='OBJECT')

        meta['dispname'] = dict(ext=0, card='GRATING')
        #meta['readout_speed'] = dict(ext=0, card='ROSPEED') #TODO: include this
        #meta['gainset'] = dict(ext=0, card='GAINSET') #TODO: include this
        meta['binning'] = dict(card=None, compound=True)

        meta['mjd'] = dict(card=None, compound=True)
        meta['exptime'] = dict(ext=0, card='EXPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS', required_ftypes=['science', 'standard'])


        # Extras for config and frametyping
        meta['dispangle'] = dict(ext=0, card='GRTILT', rtol=1e-2)
        meta['filter1'] = dict(ext=0, card='CALFILT')
        meta['decker'] = dict(ext=0, card='MASKID')
        #meta['masktype'] = dict(card=None, compound=True) #TODO: include this

        meta['idname'] = dict(ext=0, card='CCDTYPE')

        # Lamps
        meta['lampstat01'] = dict(ext=0, card='LAMPID')

        # Ingest
        self.meta = meta

    def compound_meta(self, headarr: List[fits.Header], meta_key: str):
        """
        Methods to generate meta in a more complex manner than simply
        reading from the header.

        Args:
            headarr: List[fits.Header]
              List of headers
            meta_key: str

        Returns:
            value:

        """

        if meta_key == 'binning':
            binspatial, binspec = headarr[0]['CCDSUM'].split(' ')
            return parse.binning2string(binspec, binspatial)
        elif meta_key == 'mjd':
            return Time(headarr[0]['JD'], format='jd').mjd

        # elif meta_key == 'masktype':
        #     if headarr[0]['MASKTYPE'] != 'LONGSLIT' :
        #         msgs.error('only longslit data is currently supported on SALT RSS!')
        #     return headarr[0]['MASKTYPE']

        else:
            return None

    def pypeit_file_keys(self):
        pypeit_keys = super(SALTRSSpectrograph, self).pypeit_file_keys()
        pypeit_keys += ['masktype']
        return pypeit_keys

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['science', 'standard']:
            return good_exp & (fitstbl['lampstat01'] == 'NONE') & (fitstbl['idname'] == 'OBJECT')
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'ZERO')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (fitstbl['idname'] == 'FLAT')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['lampstat01'] != 'NONE') & (fitstbl['idname'] == 'ARC')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))

        return np.zeros(len(fitstbl), dtype=bool)


class SALTRSSVisiblepectrograph(SALTRSSpectrograph):
    """
    Child to handle SALT/RSS visible beam specific code
    """
    # TODO: LOOK UP ON RSS
    def __init__(self):
        # Get it started
        super(SALTRSSVisiblepectrograph, self).__init__()
        self.spectrograph = 'salt_rss_visible'
        self.camera = 'RSSv'

    def get_detector_par(self, hdu: fits.HDUList, det: int):
        """
        Return a DectectorContainer for the current image

        Args:
            hdu (`astropy.io.fits.HDUList`):
                HDUList of the image of interest.
                Ought to be the raw file, or else..
            det (int): 1-indexed detector number

        Returns:
            :class:`pypeit.images.detector_container.DetectorContainer`:

        """
        # from: http://pysalt.salt.ac.za/proposal_calls/current/ProposalCall.html#h.47382lsidgt4
        # The dark current is typically less than 1 electron per pixel per hour.
        # #Full well depth is on the order of 170k electrons.

        # these should be the same for all detectors/not in the header
        this_detector = dict(
            binning         = self.get_meta_value(self.get_headarr(hdu), 'binning'),
            specaxis        = 0, #check
            specflip        = True,  #check
            spatflip        = False, # check
            platescale      = 0.138, # unbinned - assuming that's what's meant?
            darkcurr        = 0.5, # "less than 1" is what docs say  ¯\_(ツ)_/¯
            nonlinear       = .96, #guess, needs verification!
            mincounts       = -1e10, #guess, needs verification!
            numamplifiers   = 1,
            )


        this_detector['det'] = det
        this_detector['dataext'] = det

        detheader = hdu[det].header
        this_detector['gain'] = np.atleast_1d(detheader['GAIN'])
        this_detector['ronoise'] = np.atleast_1d(detheader['RDNOISE'])
        this_detector['saturation'] = detheader['SATURATE']
        this_detector['datasec'] = np.atleast_1d(flip_fits_slice(detheader['DATASEC']))
        this_detector['oscansec'] = np.atleast_1d(flip_fits_slice(detheader['BIASSEC']))

        return detector_container.DetectorContainer(**this_detector)


    def default_pypeit_par(self):
        """
        Set default parameters for salt rss reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'salt_rss_visible'

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.5 # check

        par['scienceframe']['process']['use_overscan'] = True
        # Make a bad pixel mask
        par['calibrations']['bpm_usebias'] = True

        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['pixelflatframe']['process']['sig_lohi'] = [10.,10.]

        # Change the wavelength calibration method
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['lamps'] = ['XeI'] # check

        # Do not flux calibrate
        par['fluxcalib'] = None

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, 120] #check
        par['calibrations']['standardframe']['exprng'] = [None, 120] #check
        par['scienceframe']['exprng'] = [90, None] #check

        par['sensfunc']['UVIS']['nresln'] = 5 #check

        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

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

        # check - this is adaptred from p200, and is definitely only right for longslits...

        angle = Angle(self.get_meta_value(scifile, 'dispangle'), unit=u.deg).rad
        # an example bar code is PL0150N001, and I *think* the 0150 is the slit width? # check
        slitwidth = (float(self.get_meta_value(scifile, 'decker').split('N')[0].replace('PL', ''))/100) * u.arcsec
        lines_mm = float(self.get_meta_value(scifile, 'dispname').replace('PG', '').strip()) / u.mm

        theta_m = 38.5 * 2*np.pi / 360. - angle
        order = 2. if lines_mm == 158. * u.mm else 1.
        platescale = 0.389 * u.arcsec / u.pix
        pix_size = 15 * u.um

        disp = np.cos(theta_m)/(lines_mm * 9 * u.imperial.inch) * 1e7 * u.AA / u.mm
        cen_wv = np.abs(1/lines_mm * (np.sin(theta_m) - np.sin(angle)) / order)
        dlam = slitwidth / platescale * pix_size / u.pix * disp

        resolving_power = cen_wv / dlam

        par['sensfunc']['UVIS']['resolution'] = resolving_power.decompose().value

        return par


class SALTRSSNIRpectrograph(SALTRSSpectrograph):
    """
    Child to handle SALT/RSS NIR beam specific code
    """
    def __init__(self):
        super(SALTRSSNIRpectrograph, self).__init__()
        self.spectrograph = 'salt_rss_nir'
        self.camera = 'RSSir'

        raise NotImplementedError('NIR beam not yet implemented')
