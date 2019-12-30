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
for key in ['SPAT', 'SLIT', 'DET', 'SCI','OBJ', 'ORDER']:
    naming_model[key.lower()] = key

extraction_data_model = {
    'BOX_WAVE': dict(otype=np.ndarray, dtype=float, desc='Boxcar Wavelengths (Angstroms)'),
                      }


class SpecObj(object):
    """Class to handle object spectra from a single exposure
    One generates one of these Objects for each spectrum in the exposure. They are instantiated by the object
    finding routine, and then all spectral extraction information for the object are assigned as attributes

    Args:
        shape (tuple): nspec, nspat
           dimensions of the spectral image that the object is identified on
        slit_spat_pos (tuple): tuple of floats (spat_left,spat_right)
            The spatial pixel location of the left and right slit trace arrays evaluated at slit_spec_pos (see below). These
            will be in the range (0,nspat)
        slit_spec_pos (float):
            The midpoint of the slit location in the spectral direction. This will typically be nspec/2, but must be in the
            range (0,nspec)
        det (int): Detector number
        config (str, optional): Instrument configuration
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
    # Attributes

    def __init__(self, shape, slit_spat_pos, slit_spec_pos, det=1, setup=None, idx=None,
                 slitid=999, orderindx=999, objtype='unknown', pypeline='unknown', spat_pixpos=None, config=None):


        #Assign from init parameters
        self.shape = shape
        self.slit_spat_pos = slit_spat_pos
        self.slit_spec_pos = slit_spec_pos
        self.setup = setup
        self.slitid = slitid
        self.det = det
        self.objtype = objtype
        self.config = config
        self.pypeline = pypeline

        # ToDo add all attributes here and to the documentaiton

        # Object finding attributes
        self.sign = 1.0
        self.objid = 999
        self.spat_fracpos = None
        self.smash_peakflux = None
        self.fwhm = None
        self.trace_spat = None
        self.spat_pixpos = spat_pixpos # Position on the image in pixels at the midpoint of the slit in spectral direction
        self.maskwidth = None
        self.min_spat = None
        self.max_spat = None
        self.prof_nsigma = None
        self.fwhmfit = None
        self.smash_nsig = None

        # Wavelength items
        self.flex_shift = 0.

        # Some things for echelle functionality
        self.ech_order = 0 # Needs a default value
        self.ech_orderindx = orderindx
        self.ech_objid = 999
        self.ech_snr = None
        self.ech_fracpos = None
        self.ech_frac_was_fit = None
        self.ech_usepca = False

        # Attributes for HAND apertures, which are object added to the extraction by hand
        self.hand_extract_spec = None
        self.hand_extract_spat = None
        self.hand_extract_det = None
        self.hand_extract_fwhm = None
        self.hand_extract_flag = False


        # Dictionaries holding boxcar and optimal extraction parameters
        self.boxcar = {}   # Boxcar extraction 'wave', 'counts', 'var', 'sky', 'mask', 'flam', 'flam_var'
        self.optimal = {}  # Optimal extraction 'wave', 'counts', 'var', 'sky', 'mask', 'flam', 'flam_var'



        # Generate IDs
        #self.slitid = int(np.round(self.slitcen*1e4))
        #self.objid = int(np.round(xobj*1e3))

        # Set index
        if idx is None:
            self.set_idx()
        else:
            self.idx = idx
        #

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
        Generate a unique index for this spectrum

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

    def check_trace(self, trace, toler=1.):
        """
        Check that the input trace matches the defined specobjexp

        Args:
            trace (ndarray): Trace of the object
            toler (float): Tolerance for matching, in pixels

        Returns:
            bool:  True = match within tolerance

        """
        # Trace
        yidx = int(np.round(self.ypos*trace.size))
        obj_trc = trace[yidx]
        # self
        nslit = self.shape[1]*(self.xslit[1]-self.xslit[0])
        xobj_pix = self.shape[1]*self.xslit[0] + nslit*self.xobj
        # Check
        if np.abs(obj_trc-xobj_pix) < toler:
            return True
        else:
            return False

    def copy(self):
        """
        Generate a copy of this object

        Returns:
            SpecObj

        """
        sobj_copy = SpecObj(self.shape, self.slit_spat_pos, self.slit_spec_pos) # Instantiate
#        sobj_copy.__dict__ = self.__dict__.copy() # Copy over all attributes
#        sobj_copy.boxcar = self.boxcar.copy() # Copy boxcar and optimal dicts
#        sobj_copy.optimal = self.optimal.copy()
        sobj_copy.__dict__ = copy.deepcopy(self.__dict__)
        sobj_copy.boxcar = copy.deepcopy(self.boxcar) # Copy boxcar and optimal dicts
        sobj_copy.optimal = copy.deepcopy(self.optimal)
        # These attributes are numpy arrays that don't seem to copy from the lines above??
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

    def __getitem__(self, key):
        # Access the DB groups
        return getattr(self, key)

    def __repr__(self):
        # Create a single summary table for one object, so that the representation is always the same
        sobjs = SpecObjs(specobjs=[self])
        return sobjs.summary.__repr__()


class SpecObjs(object):
    """
    Object to hold a set of SpecObj objects

    Args:
        specobjs (ndarray or list, optional):  One or more SpecObj objects

    Internals:
        summary (astropy.table.Table):


    __getitem__ is overloaded to allow one to pull an attribute or a
                portion of the SpecObjs list
        Args:
            item (str or int (or slice)
        Returns:
            item (object, SpecObj or SpecObjs):  Depends on input item..

    __setitem__ is over-loaded using our custom set() method
        Args:
            name (str):  Item to set
            value (anything) : Value of the item
        Returns:
    __getattr__ is overloaded to generate an array of attribute 'k' from the specobjs
        First attempts to grab data from the Summary table, then the list
    """

    def __init__(self, specobjs=None):
        if specobjs is None:
            self.specobjs = np.array([])
        else:
            if isinstance(specobjs, (list, np.ndarray)):
                specobjs = np.array(specobjs)
            self.specobjs = specobjs

        # Internal summary Table
        self.build_summary()

    @property
    def nobj(self):
        """
        Return the number of SpecObj objects

        Returns:
            int

        """
        return self.specobjs.size

    def get_std(self):
        """
        Return the standard star from this Specobjs. For MultiSlit this
        will be a single specobj in SpecObjs container, for Echelle it
        will be the standard for all the orders.

        Args:

        Returns:
            SpecObj or SpecObjs

        """
        # Is this MultiSlit or Echelle
        pypeline = (self.pypeline)[0]
        if 'MultiSlit' in pypeline:
            nspec = self[0].optimal['COUNTS'].size
            SNR = np.zeros(self.nobj)
            # Have to do a loop to extract the counts for all objects
            for iobj in range(self.nobj):
                SNR[iobj] = np.median(self[iobj].optimal['COUNTS']*np.sqrt(self[iobj].optimal['COUNTS_IVAR']))
            istd = SNR.argmax()
            return SpecObjs(specobjs=[self[istd]])
        elif 'Echelle' in pypeline:
            uni_objid = np.unique(self.ech_objid)
            uni_order = np.unique(self.ech_orderindx)
            nobj = len(uni_objid)
            norders = len(uni_order)
            SNR = np.zeros((norders, nobj))
            for iobj in range(nobj):
                for iord in range(norders):
                    ind = (self.ech_objid == uni_objid[iobj]) & (self.ech_orderindx == uni_order[iord])
                    spec = self[ind]
                    SNR[iord, iobj] = np.median(spec[0].optimal['COUNTS']*np.sqrt(spec[0].optimal['COUNTS_IVAR']))
            SNR_all = np.sqrt(np.sum(SNR**2,axis=0))
            objid_std = uni_objid[SNR_all.argmax()]
            indx = self.ech_objid == objid_std
            return SpecObjs(specobjs=self[indx])
        else:
            msgs.error('Unknown pypeline')


    def append_neg(self, sobjs_neg):
        """
        Append negative objects and change the sign of their objids for IR reductions

        Args:
            sobjs_neg (SpecObjs):

        """

        # Assign the sign and the objids
        for spec in sobjs_neg:
            spec.sign = -1.0
            try:
                spec.objid = -spec.objid
            except TypeError:
                pass
            try:
                spec.ech_objid = -spec.ech_objid
            except TypeError:
                pass

        self.add_sobj(sobjs_neg)

        # Sort objects according to their spatial location. Necessary for the extraction to properly work
        if self.nobj > 0:
            spat_pixpos = self.spat_pixpos
            self.specobjs = self.specobjs[spat_pixpos.argsort()]

    def purge_neg(self):
        """
        Purge negative objects from specobjs for IR reductions

        """
        # Assign the sign and the objids
        if self.nobj > 0:
            index = (self.objid < 0) | (self.ech_objid < 0)
            self.remove_sobj(index)


    def add_sobj(self, sobj):
        """
        Add one or more SpecObj
        The summary table is rebuilt

        Args:
            sobj (SpecObj or list or ndarray):  On or more SpecObj objects

        Returns:


        """
        if isinstance(sobj, SpecObj):
            self.specobjs = np.append(self.specobjs, [sobj])
        elif isinstance(sobj, (np.ndarray,list)):
            self.specobjs = np.append(self.specobjs, sobj)
        elif isinstance(sobj, SpecObjs):
            self.specobjs = np.append(self.specobjs, sobj)

        # Rebuild summary table
        self.build_summary()

    def build_summary(self):
        """
        Build the internal Summary Table

        Returns:
            Builds self.summary Table internally

        """
        # Dummy?
        if len(self.specobjs) == 0:
            self.summary = Table()
            return
        #
        atts = self.specobjs[0].__dict__.keys()
        uber_dict = {}
        for key in atts:
            uber_dict[key] = []
            for sobj in self.specobjs:
                uber_dict[key] += [getattr(sobj, key)]
        # Build it
        self.summary = Table(uber_dict)

    def remove_sobj(self, index):
        """
        Remove an object

        Args:
            index: int

        Returns:

        """
        msk = np.ones(self.specobjs.size, dtype=bool)
        msk[index] = False
        # Do it
        self.specobjs = self.specobjs[msk]
        # Update
        self.build_summary()


    def copy(self):
        """
        Generate a copy of self

        Returns:
            SpecObjs

        """
        sobj_copy = SpecObjs()
        for sobj in self.specobjs:
            sobj_copy.add_sobj(sobj.copy())
        sobj_copy.build_summary()
        return sobj_copy

    def set_idx(self):
        """
        Set the idx in all the SpecObj
        Update the summary Table

        Returns:

        """
        for sobj in self.specobjs:
            sobj.set_idx()
        self.build_summary()


    def __getitem__(self, item):
        if isinstance(item, str):
            return self.__getattr__(item)
        elif isinstance(item, (int, np.integer)):
            return self.specobjs[item] # TODO Is this using pointers or creating new data????
        elif (isinstance(item, slice) or  # Stolen from astropy.table
            isinstance(item, np.ndarray) or
            isinstance(item, list) or
            isinstance(item, tuple) and all(isinstance(x, np.ndarray) for x in item)):
            # here for the many ways to give a slice; a tuple of ndarray
            # is produced by np.where, as in t[np.where(t['a'] > 2)]
            # For all, a new table is constructed with slice of all columns
            return SpecObjs(specobjs=self.specobjs[item])


    # TODO this code fails for assignments of this nature sobjs[:].attribute = np.array(5)
    def __setitem__(self, name, value):
        self.set(slice(0,self.nobj), name, value)

    def set(self, islice, attr, value):
        """
        Set the attribute for a slice of the specobjs

        Args:
            islice (int, ndarray of bool, slice):  Indicates SpecObj to affect
            attr (str):
            value (anything) : Value of the item

        Returns:

        """
        sub_sobjs = self.specobjs[islice]
        if isiterable(value):
            if sub_sobjs.size == len(value):  # Assume you want each paired up
                for kk,sobj in enumerate(sub_sobjs):
                    setattr(sobj, attr, value[kk])
                    return
        # Assuming scalar assignment
        if isinstance(sub_sobjs, SpecObj):
            setattr(sub_sobjs, attr, value)
        else:
            for sobj in sub_sobjs:
                setattr(sobj, attr, value)
        return


    def __getattr__(self, k):
        # Overloaded
        self.build_summary()
        # Special case(s)
        if k in self.summary.keys():  # _data
            lst = self.summary[k]
        else:
            lst = None
        # specobjs last!
        if lst is None:
            if len(self.specobjs) == 0:
                raise ValueError("Attribute not available!")
            try:
                lst = [getattr(specobj, k) for specobj in self.specobjs]
            except ValueError:
                raise ValueError("Attribute does not exist")
        # Recast as an array
        return lst_to_array(lst)

    # Printing
    def __repr__(self):
        return self.summary.__repr__()

    def __len__(self):
        return len(self.specobjs)

    def keys(self):
        self.build_summary()
        return self.summary.keys()



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


def lst_to_array(lst, mask=None):
    """
    Simple method to convert a list to an array

    Allows for a list of Quantity objects

    Args:
        lst : list
          Should be number or Quantities
        mask (ndarray of bool, optional):  Limit to a subset of the list.  True=good

    Returns:
        ndarray or Quantity array:  Converted list

    """
    if mask is None:
        mask = np.array([True]*len(lst))
    if isinstance(lst[0], Quantity):
        return Quantity(lst)[mask]
    else:
        return np.array(lst)[mask]


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
