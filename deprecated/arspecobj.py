# Module for handling extracted spectra
#  Includes ArSpecObj class
from __future__ import absolute_import, division, print_function

import copy
from collections import OrderedDict

import numpy as np

from astropy import units
from astropy.table import Table

from pypit import msgs
from pypit import arparse
from pypit.core import artraceslits
from pypit import ardebug as debugger

class SpecObjExp(object):
    """Class to handle object spectra from a single exposure
    One generates one of these Objects for each spectrum in the exposure

    Parameters:
    ----------
    shape: tuple
       row,col of the frame
    config: str
       Instrument configuration
    scidx: int
       Exposure index (max=9999)
    det: int
       Detector index (max=99)
    xslit: tuple
       float (0-1), float (0-1)
       left, right of slit in fraction of total (trimmed) detector size defined at ypos
    ypos: float
       ypos of slit in fraction of total (trimmed) detector size 
    xobj: float
       float (0-1)
       Position of object in fraction of total slit at same ypos that defines the slit
    objtype: str, optional
       Type of object ('unknown', 'standard', 'science')

    Attributes:
    ----------
    slitcen: float
       Center of slit in fraction of total (trimmed) detector size at ypos
    slitid: int
       Identifier for the slit (max=9999) 
    objid: int
       Identifier for the object (max=999)
    """
    # Attributes
    # Init
    def __init__(self, shape, config, scidx, det, xslit, ypos, xobj, objtype='unknown'):
        self.shape = shape
        self.config = config
        self.scidx = copy.deepcopy(scidx)
        self.det = det
        self.xslit = xslit
        self.ypos = ypos
        self.slitcen = np.mean([xslit[0], xslit[1]])
        self.xobj = xobj
        self.objtype = objtype

        # Generate IDs
        self.slitid = int(np.round(self.slitcen*1e4))
        self.objid = int(np.round(xobj*1e3))

        # Set index
        self.set_idx()

        # Items that are generally filled
        self.boxcar = {}   # Boxcar extraction 'wave', 'counts', 'var', 'sky', 'mask', 'flam', 'flam_var'
        self.optimal = {}  # Optimal extraction 'wave', 'counts', 'var', 'sky', 'mask', 'flam', 'flam_var'
        #

    def set_idx(self):
        # Generate a unique index for this exposure
        #self.idx = '{:02d}'.format(self.setup)
        self.idx = 'O{:03d}'.format(self.objid)
        self.idx += '-S{:04d}'.format(self.slitid)
        sdet = arparse.get_dnum(self.det, prefix=False)
        self.idx += '-D{:s}'.format(sdet)
        self.idx += '-I{:04d}'.format(self.scidx)

    def check_trace(self, trace, toler=1.):
        """Check that the input trace matches the defined specobjexp

        Parameters:
        ----------
        trace: ndarray
          Trace of the object
        toler: float, optional
          Tolerance for matching, in pixels
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
        slf = SpecObjExp(self.shape, self.config, self.scidx, self.det, self.xslit, self.ypos, self.xobj,
                       objtype=self.objtype)
        slf.boxcar = self.boxcar.copy()
        slf.optimal = self.optimal.copy()
        return slf

    # Printing
    def __repr__(self):
        # Generate sets string
        sdet = arparse.get_dnum(self.det, prefix=False)
        return ('<SpecObjExp: {:s} == Setup {:s} Object at {:g} in Slit at {:g} with det={:s}, scidx={:d} and objtype={:s}>'.format(
                self.idx, self.config, self.xobj, self.slitcen, sdet, self.scidx, self.objtype))


def init_exp(lordloc, rordloc, shape, maskslits, det, scidx, fitstbl, tracelist, binning=None,
             ypos=0.5, objtype='unknown'): #**kwargs):
    """ Generate a list of SpecObjExp objects for a given exposure

    Parameters
    ----------
    self
       Instrument "setup" (min=10,max=99)
    scidx : int
       Index of file
    det : int
       Detector index 
    tracelist : list of dict
       Contains trace info
    ypos : float, optional [0.5]
       Row on trimmed detector (fractional) to define slit (and object)

    Returns
    -------
    specobjs : list
      List of SpecObjExp objects
    """

    # Init
    specobjs = []
    if fitstbl is None:
        fitsrow = None
    else:
        fitsrow = fitstbl[scidx]
    config = instconfig(fitsrow=fitsrow, binning=binning)
    slits = range(len(tracelist))
    gdslits = np.where(~maskslits)[0]

    # Loop on slits
    for sl in slits:
        specobjs.append([])
        # Analyze the slit?
        if sl not in gdslits:
            specobjs[sl].append(None)
            continue
        # Object traces
        if tracelist[sl]['nobj'] != 0:
            # Loop on objects
            #for qq in range(trc_img[sl]['nobj']):
            for qq in range(tracelist[sl]['traces'].shape[1]):
                slitid, slitcen, xslit = artraceslits.get_slitid(shape, lordloc, rordloc,
                                                                 sl, ypos=ypos)
                # xobj
                _, xobj = get_objid(lordloc, rordloc, sl, qq, tracelist, ypos=ypos)
                # Generate
                spec_obj_shape = (tracelist[0]['object'].shape[:2]) \
                                        if tracelist[sl]['object'] is None \
                                        else (tracelist[sl]['object'].shape[:2])
                specobj = SpecObjExp(spec_obj_shape, config, scidx, det, xslit, ypos, xobj,
                                     objtype=objtype) #**kwargs)
                # Add traces
                specobj.trace = tracelist[sl]['traces'][:, qq]
                # Append
                specobjs[sl].append(copy.deepcopy(specobj))
        else:
            msgs.warn("No objects for slit {0:d}".format(sl+1))
            specobjs[sl].append(None)
    # Return
    return specobjs


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
        odict[iprs[0]] = int(iprs[1:])
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
    gdrow = (np.abs(tbl['O']-odict['O']) < otol) & (np.abs(tbl['S']-odict['S']) < stol) & (tbl['D'] == odict['D'])
    if np.sum(gdrow) == 0:
        return None
    else:
        return np.array(objects)[gdrow].tolist(), np.where(gdrow)[0].tolist()



def get_objid(lordloc, rordloc, islit, iobj, trc_img, ypos=0.5):
    """ Convert slit position to a slitid
    Parameters
    ----------
    det : int
    islit : int
    iobj : int
    trc_img : list of dict
    ypos : float, optional

    Returns
    -------
    objid : int
    xobj : float
    """
    yidx = int(np.round(ypos*lordloc.shape[0]))
    pixl_slit = lordloc[yidx, islit]
    pixr_slit = rordloc[yidx, islit]
    #
    xobj = (trc_img[islit]['traces'][yidx,iobj]-pixl_slit) / (pixr_slit-pixl_slit)
    objid= int(np.round(xobj*1e3))
    # Return
    return objid, xobj


def instconfig(fitsrow=None, binning=None):
    """ Returns a unique config string

    Parameters
    ----------
    fitsrow : Row
    binnings : str, optional

    Returns
    -------
    config : str
    """

    config_dict = OrderedDict()
    config_dict['S'] = 'slitwid'
    config_dict['D'] = 'dichroic'
    config_dict['G'] = 'dispname'
    config_dict['T'] = 'dispangle'
    #
    config = ''
    for key in config_dict.keys():
        try:
            comp = str(fitsrow[config_dict[key]])
        except (KeyError, TypeError):
            comp = '0'
        #
        val = ''
        for s in comp:
            if s.isdigit():
                val += s
        config = config + key+'{:s}-'.format(val)
    # Binning
    if binning is None:
        msgs.warn("Assuming 1x1 binning for your detector")
        binning = '1x1'
    val = ''
    for s in binning:
        if s.isdigit():
            val = val + s
    config += 'B{:s}'.format(val)
    # Return
    return config


def dummy_specobj(fitstbl, det=1, extraction=True):
    """ Generate dummy specobj classes
    Parameters
    ----------
    fitstbl : Table
      Expecting the fitsdict from dummy_fitsdict
    Returns
    -------

    """
    shape = fitstbl['naxis1'][0], fitstbl['naxis0'][0]
    config = 'AA'
    scidx = 5 # Could be wrong
    xslit = (0.3,0.7) # Center of the detector
    ypos = 0.5
    xobjs = [0.4, 0.6]
    specobjs = []
    for xobj in xobjs:
        specobj = SpecObjExp(shape, config, scidx, det, xslit, ypos, xobj)
        # Dummy extraction?
        if extraction:
            npix = 2001
            specobj.boxcar['wave'] = np.linspace(4000., 6000., npix)*units.AA
            specobj.boxcar['counts'] = 50.*(specobj.boxcar['wave'].value/5000.)**-1.
            specobj.boxcar['var']  = specobj.boxcar['counts'].copy()
        # Append
        specobjs.append(specobj)
    # Return
    return specobjs




