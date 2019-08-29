""" Module for the SpecObjs and SpecObj classes
"""
import copy
import re

import numpy as np

from scipy import interpolate

from astropy import units
from astropy.table import Table
from astropy.units import Quantity
from astropy.utils import isiterable

from linetools.spectra import xspectrum1d

from pypeit import msgs
from pypeit.core import parse

naming_model = {}
for skey in ['SPAT', 'SLIT', 'DET', 'SCI','OBJ', 'ORDER']:
    naming_model[skey.lower()] = skey

# Data model -- Put here to be able to reach it without instantiating the class
#  These are outward facing items, i.e. items that the user will receive and use.
data_model = {
    'idx': dict(otype=str, desc='Name of the object in PypeIt convention'),
    'trace_spat': dict(otype=np.ndarray, atype=float, desc='Object trace along the spec (spatial pixel)'),
    'fwhm': dict(otype=float, desc='Spatial FWHM of the object (pixels)'),
    'fwhmfit': dict(otype=np.ndarray, desc='Spatial FWHM across the detector (pixels)'),
    'OPT_WAVE': dict(otype=np.ndarray, atype=float, desc='Optimal Wavelengths (Angstroms)'),
    'OPT_COUNTS': dict(otype=np.ndarray, atype=float, desc='Optimal flux (counts)'),
    'OPT_COUNTS_IVAR': dict(otype=np.ndarray, atype=float,
                            desc='Inverse variance of optimally extracted flux using modelivar image (counts^2)'),
    'OPT_COUNTS_SIG': dict(otype=np.ndarray, atype=float,
                             desc='Optimally extracted noise from IVAR (counts)'),
    'OPT_COUNTS_NIVAR': dict(otype=np.ndarray, atype=float,
                             desc='Optimally extracted noise variance, sky+read noise only (counts^2)'),
    'OPT_MASK': dict(otype=np.ndarray, atype=bool, desc='Mask for optimally extracted flux'),
    'OPT_COUNTS_SKY': dict(otype=np.ndarray, atype=float, desc='Optimally extracted sky (counts)'),
    'OPT_COUNTS_RN': dict(otype=np.ndarray, atype=float, desc='Optimally extracted RN squared (counts)'),
    'OPT_FRAC_USE': dict(otype=np.ndarray, atype=float,
                         desc='Fraction of pixels in the object profile subimage used for this extraction'),
    'OPT_CHI2': dict(otype=np.ndarray, atype=float,
                     desc='Reduced chi2 of the model fit for this spectral pixel'),
    #
    'BOX_WAVE': dict(otype=np.ndarray, atype=float, desc='Optimal Wavelengths (Angstroms)'),
    'BOX_COUNTS': dict(otype=np.ndarray, atype=float, desc='Optimal flux (counts)'),
    'BOX_COUNTS_IVAR': dict(otype=np.ndarray, atype=float,
                            desc='Inverse variance of optimally extracted flux using modelivar image (counts^2)'),
    'BOX_COUNTS_SIG': dict(otype=np.ndarray, atype=float,
                           desc='Optimally extracted noise from IVAR (counts)'),
    'BOX_COUNTS_NIVAR': dict(otype=np.ndarray, atype=float,
                             desc='Optimally extracted noise variance, sky+read noise only (counts^2)'),
    'BOX_MASK': dict(otype=np.ndarray, atype=bool, desc='Mask for optimally extracted flux'),
    'BOX_COUNTS_SKY': dict(otype=np.ndarray, atype=float, desc='Optimally extracted sky (counts)'),
    'BOX_COUNTS_RN': dict(otype=np.ndarray, atype=float, desc='Optimally extracted RN squared (counts)'),
    'BOX_FRAC_USE': dict(otype=np.ndarray, atype=float,
                         desc='Fraction of pixels in the object profile subimage used for this extraction'),
    'BOX_CHI2': dict(otype=np.ndarray, atype=float,
                     desc='Reduced chi2 of the model fit for this spectral pixel'),
    'BOX_RADIUS': dict(otype=float, desc='Size of boxcar radius (pixels)'),
}


class SpecObj(dict):
    """Class to handle object spectra from a single exposure
    One generates one of these Objects for each spectrum in the exposure. They are instantiated by the object
    finding routine, and then all spectral extraction information for the object are assigned as attributes

    Args:
        slit_spat_pos (tuple): tuple of floats (spat_left,spat_right)
            The spatial pixel location of the left and right slit trace arrays evaluated at slit_spec_pos (see below). These
            will be in the range (0,nspat)
        slit_spec_pos (float):
            The midpoint of the slit location in the spectral direction. This will typically be nspec/2, but must be in the
            range (0,nspec)
        det (int): Detector number
        objtype (str, optional)
           Type of object ('unknown', 'standard', 'science')
        slitid (int, optional):
           Identifier for the slit (max=9999)

    Attributes:
        slitcen (float): Center of slit in fraction of total (trimmed) detector size at ypos
        objid (int): Identifier for the object (max=999)
        flex_shift (float): Flexure correction in pixels

    Extraction dict's
        'WAVE' : wave_opt  # Optimally extracted wavelengths
        'COUNTS' : flux_opt  # Optimally extracted flux
        'COUNTS_IVAR' : mivar_opt  # Inverse variance of optimally extracted flux using modelivar image
        'COUNTS_NIVAR' : nivar_opt  # Optimally extracted noise variance (sky + read noise) only
        'MASK' : mask_opt  # Mask for optimally extracted flux
        'COUNTS_SKY' : sky_opt  # Optimally extracted sky
        'COUNTS_RN' : rn_opt  # Square root of optimally extracted read noise squared
        'FRAC_USE' : frac_use  # Fraction of pixels in the object profile subimage used for this extraction
        'CHI2' : chi2  # Reduced chi2 of the model fit for this spectral pixel
    """

    def __init__(self, slit_spat_pos, slit_spec_pos, det=1, idx=None,
                 indict=None, slitid=999, orderindx=999, objtype='unknown',
                 pypeline='unknown', spat_pixpos=None):

        # For copying the object
        if indict is not None:
            if '_SpecObj_initialised' in indict:
                indict.pop('_SpeObj_initialised')
            self.__dict__ = indict
        else:
            # set any attributes here - before initialisation
            # these remain as normal attributes
            # We may wish to eliminate *all* of these
            self.pypeline = pypeline
            self.slit_spat_pos = slit_spat_pos
            self.slit_spec_pos = slit_spec_pos
            self.det = det
            self.spat_pixpos = spat_pixpos   # Position on the image in pixels at the midpoint of the slit in spectral direction
            self.slitid = slitid
            self.objid = 999

            # Object finding
            self.spat_fracpos = None
            self.smash_peakflux = None
            self.smash_nsig = None
            self.maskwidth = None

            # Object profile
            self.prof_nsigma = None
            self.sign = 1.0
            self.min_spat = None
            self.max_spat = None

            # Trace
            self.trace_spec = None  # Only for debuggin, internal plotting

            # Echelle specific
            self.ech_orderindx = orderindx

        # after initialisation, setting attributes is the same as setting an item
        self.__initialised = True

        # Set index
        if idx is None:
            self.set_idx()
        else:
            self.idx = idx

    def __getattr__(self, item):
        """Maps values to attributes.
        Only called if there *isn't* an attribute with this name
        """
        try:
            return self.__getitem__(item)
        except KeyError:
            raise AttributeError(item)

    def __setattr__(self, item, value):

        if not '_SpecObj__initialised' in self.__dict__:  # this test allows attributes to be set in the __init__ method
            return dict.__setattr__(self, item, value)
        elif item in self.__dict__:       # any normal attributes are handled normally
            dict.__setattr__(self, item, value)
        else:
            if item not in data_model.keys():
                raise IOError("Cannot set this attribute.  It is not in the data model")
            if not isinstance(value, data_model[item]['otype']):
                print("Wrong data type for attribute: {}".format(item))
                print("Allowed type(s) are: {}".format(data_model[item]['otype']))
                raise IOError("Try again")
            # Special checking for arrays?
            self.__setitem__(item, value)

    @staticmethod
    def sobjs_key():
        """
        This function returns the dictionary that defines the mapping between specobjs attributes and the fits header
        cards

        Returns:
            dict:

        """
        sobjs_key_dict = dict(det='DET',
                              objid='OBJID',
                              slitid='SLITID',
                              ech_objid='ECHOBJID',
                              ech_orderindx='ECHOINDX',
                              ech_order='ECHORDER',
                              pypeline='PYPELINE')

        return sobjs_key_dict

    def set_idx(self):
        """
        Generate a unique index for this spectrum based on the
        slit/order, its position and for multi-slit the detector.

        Sets self.idx internally

        Returns:
            str: :attr:`self.idx`

        """
        # Detector string
        sdet = parse.get_dnum(self.det, prefix=False)

        if 'Echelle' in self.pypeline:
            # ObjID
            self.idx = naming_model['obj']
            if self.ech_objid is None:
                self.idx += '----'
            else:
                self.idx += '{:04d}'.format(self.ech_objid)
            self.idx += '-'+naming_model['order']
            # Order
            if self.ech_orderindx is None:
                self.idx += '----'
            else:
                self.idx += '{:04d}'.format(self.ech_order)
        else:
            # Spat
            self.idx = naming_model['spat']
            if self.spat_pixpos is None:
                self.idx += '----'
            else:
                self.idx += '{:04d}'.format(int(np.rint(self.spat_pixpos)))
            # Slit
            self.idx += '-'+naming_model['slit']
            if self.slitid is None:
                self.idx += '----'
            else:
                self.idx += '{:04d}'.format(self.slitid)

        self.idx += '-{:s}{:s}'.format(naming_model['det'], sdet)
        # Return
        return self.idx

    def copy(self):
        """
        Generate a copy of this object

        Returns:
            SpecObj

        """
        sobj_copy = SpecObj(self.slit_spat_pos, self.slit_spec_pos,
                            indict=copy.deepcopy(self.__dict__))  # Instantiate
        for key, item in self.items():
            sobj_copy[key] = item
        return sobj_copy

    def flexure_interp(self, sky_wave, fdict):
        """
        Apply interpolation with the flexure dict

        Args:
            sky_wave (np.ndarray): Wavelengths of the extracted sky
            fdict (dict): Holds the various flexure items

        Returns:
            np.ndarray:  New sky spectrum (mainly for QA)

        """
        # Simple interpolation to apply
        npix = len(sky_wave)
        x = np.linspace(0., 1., npix)
        # Apply
        for attr in ['boxcar', 'optimal']:
            if not hasattr(self, attr):
                continue
            if 'WAVE' in getattr(self, attr).keys():
                msgs.info("Applying flexure correction to {0:s} extraction for object:".format(attr) +
                          msgs.newline() + "{0:s}".format(str(self)))
                f = interpolate.interp1d(x, sky_wave, bounds_error=False, fill_value="extrapolate")
                getattr(self, attr)['WAVE'] = f(x + fdict['shift'] / (npix - 1)) * units.AA
        # Shift sky spec too
        cut_sky = fdict['sky_spec']
        x = np.linspace(0., 1., cut_sky.npix)
        f = interpolate.interp1d(x, cut_sky.wavelength.value, bounds_error=False, fill_value="extrapolate")
        twave = f(x + fdict['shift'] / (cut_sky.npix - 1)) * units.AA
        new_sky = xspectrum1d.XSpectrum1D.from_tuple((twave, cut_sky.flux))
        # Save
        self.flex_shift = fdict['shift']
        # Return
        return new_sky

    def to_xspec1d(self, extraction='optimal'):
        """
        Convert the SpecObj to an XSpectrum1D object

        Args:
            extraction (str): Extraction method to convert

        Returns:
            linetools.spectra.xspectrum1d.XSpectrum1D:  Spectrum object

        """
        extract = getattr(self, extraction)
        if len(extract) == 0:
            msgs.warn("This object has not been extracted with extract={}".format(extraction))
        if 'FLAM' in extract:
            flux = extract['FLAM']
            sig = extract['FLAM_SIG']
        else:
            flux = extract['COUNTS']
            sig = np.zeros_like(flux)
            gdc = extract['COUNTS_IVAR'] > 0.
            sig[gdc] = 1./np.sqrt(extract['COUNTS_IVAR'][gdc])
        # Create
        xspec = xspectrum1d.XSpectrum1D.from_tuple((extract['WAVE'], flux, sig))
        return xspec

    def show(self, extraction='optimal'):
        """
        Show the spectrum by converting it to a XSpectrum1D object

        Args:
            extraction (str): Extraction option 'optimal' or 'boxcar'

        Returns:

        """
        extract = getattr(self, extraction)
        # Generate an XSpec
        xspec = self.to_xspec1d(extraction=extraction)
        if xspec is None:
            return
        xspec.plot(xspec=True)

    #def __getitem__(self, key):
    #    # Access the DB groups
    #    return getattr(self, key)

    #def __repr__(self):
    #    # Create a single summary table for one object, so that the representation is always the same
    #    sobjs = SpecObjs(specobjs=[self])
    #    return sobjs.summary.__repr__()



def objnm_to_dict(objnm):
    """ Convert an object name or list of them into a dict

    Parameters
    ----------
    objnm : str or list of str

    Returns
    -------
    odict : dict
      Object value or list of object values
    """
    if isinstance(objnm, list):
        tdict = {}
        for kk,iobj in enumerate(objnm):
            idict = objnm_to_dict(iobj)
            if kk == 0:
                for key in idict.keys():
                    tdict[key] = []
            # Fill
            for key in idict.keys():
                tdict[key].append(idict[key])
        # Generate the Table
        return tdict
    # Generate the dict
    prs = objnm.split('-')
    odict = {}
    for iprs in prs:
        # Find first character that is an integer
        idig = re.search("\d", iprs).start()
        odict[iprs[:idig]] = int(iprs[idig:])
    # Return
    return odict


def mtch_obj_to_objects(iobj, objects, stol=50, otol=10, **kwargs):
    """
    Parameters
    ----------
    iobj : str
      Object identifier in format O###-S####-D##
    objects : list
      List of object identifiers
    stol : int
      Tolerance in slit matching
    otol : int
      Tolerance in object matching

    Returns
    -------
    matches : list
      indices of matches in objects
      None if none
    indcies : list

    """
    # Parse input object
    odict = objnm_to_dict(iobj)
    # Generate a Table of the objects
    tbl = Table(objnm_to_dict(objects))

    # Logic on object, slit and detector [ignoring sciidx for now]
    gdrow = (np.abs(tbl[naming_model['spat']]-odict[naming_model['spat']]) < otol) & (
            np.abs(tbl[naming_model['slit']]-odict[naming_model['slit']]) < stol) & (
            tbl[naming_model['det']] == odict[naming_model['det']])
    if np.sum(gdrow) == 0:
        return None
    else:
        return np.array(objects)[gdrow].tolist(), np.where(gdrow)[0].tolist()



def dummy_specobj(shape, det=1, extraction=True):
    """ Generate dummy specobj classes
    Parameters
    ----------
    shape : tuple
      naxis1, naxis0
    Returns
    sobj_list: list
      Pair of SpecObj objects
    -------

    """
    config = 'AA'
    scidx = 5 # Could be wrong
    xslit = (0.3,0.7) # Center of the detector
    ypos = 0.5
    xobjs = [0.4, 0.6]
    sobj_list = []
    for jj,xobj in enumerate(xobjs):
        specobj = SpecObj(shape, 1240, xslit, spat_pixpos=900, det=det, config=config)
        specobj.slitid = jj+1
        #specobj = SpecObj(shape, config, scidx, det, xslit, ypos, xobj)
        # Dummy extraction?
        if extraction:
            npix = 2001
            specobj.boxcar['WAVE'] = np.linspace(4000., 6000., npix)*units.AA
            specobj.boxcar['COUNTS'] = 50.*(specobj.boxcar['WAVE'].value/5000.)**-1.
            specobj.boxcar['COUNTS_IVAR']  = 1./specobj.boxcar['COUNTS'].copy()
        # Append
        sobj_list.append(specobj)
    # Return
    return sobj_list




def unravel_specobjs(specobjs):
    """
    Likely to be Deprecated

    Method to unwrap nested specobjs objects into a single list

    Args:
        specobjs (list of lists or list of SpecObj):

    Returns:
        list: list of SpecObj

    """
    # Wrapped is all None and lists
    ans = [isinstance(ispec, (list, type(None))) for ispec in specobjs]
    if np.all(ans):
        all_specobj = []
        for det in range(len(specobjs)):           # detector loop
            if specobjs[det] is None:
                continue
            for sl in range(len(specobjs[det])):   # slit loop
                for spobj in specobjs[det][sl]:    # object loop
                    all_specobj.append(spobj)
    else:
        all_specobj = specobjs
    # Return
    return all_specobj
