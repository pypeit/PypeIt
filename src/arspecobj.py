# Module for handling extracted spectra
#  Includes ArSpecObj class
import numpy as np
import copy

import armsgs

# Logging
msgs = armsgs.get_logger()


class SpecObjExp(object):
    '''Class to handle object spectra from a single exposure 
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
       x0, x1 of slit in fraction of total (trimmed) detector size defined at ypos
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
    '''
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
        self.slitid= int(np.round(self.slitcen*1e4))
        self.objid= int(np.round(xobj*1e3))

        # Generate a unique index for this expsoure
        #self.idx = '{:02d}'.format(self.setup)
        self.idx = 'O{:03d}'.format(self.objid)
        self.idx += '-S{:04d}'.format(self.slitid)
        self.idx += '-D{:02d}'.format(self.det)
        self.idx += '-I{:04d}'.format(self.scidx)

        # Items that are generally filled
        self.boxcar = {}   # Boxcar extraction 'wave', 'counts', 'var', 'sky', 'mask', 'flam', 'flam_var'
        self.optimal = {}  # Optimal extraction 'wave', 'counts', 'var', 'sky', 'mask', 'flam', 'flam_var'
        #
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

    # Printing
    def __repr__(self):
        # Generate sets string
        return ('[SpecObjExp: {:s} == Setup {:s} Object at {:g} in Slit at {:g} with det={:d}, scidx={:d} and objtype={:s}]'.format(
                self.idx, self.config, self.xobj, self.slitcen, self.det, self.scidx, self.objtype))

def init_exp(slf, scidx, det, fitsdict, trc_img=None, ypos=0.5, **kwargs):
    """Generate a list of SpecObjExp objects for a given exposure

    Parameters
    ----------
    self
       Instrument "setup" (min=10,max=99)
    scidx : int
       Index of file
    det : int
       Detector index 
    ypos : float, optional [0.5]
       Row on trimmed detector (fractional) to define slit (and object)
    """
    from armlsd import instconfig

    # Init
    specobjs = []
    config = instconfig(slf, det, scidx, fitsdict)
    yidx = int(np.round(ypos*slf._lordloc[det-1].shape[0]))
    pixl_slits = slf._lordloc[det-1][yidx, :]
    pixr_slits = slf._rordloc[det-1][yidx, :]
    #
    if trc_img is not None: # Object traces
        for qq in xrange(trc_img['traces'].shape[1]): # Loop on objects
            # Find the slit
            gds = np.where( (trc_img['traces'][yidx,qq]>pixl_slits) & 
                (trc_img['traces'][yidx,qq]<pixr_slits))[0]
            if len(gds) != 1:
                msgs.error('arspecobj.init_exp: Problem finding the slit')
            else:
                islit = gds[0]
                pixl_slit = pixl_slits[islit]
                pixr_slit = pixr_slits[islit]
                xl_slit = pixl_slit/trc_img['object'].shape[1]
                xr_slit = pixr_slit/trc_img['object'].shape[1]
            # xobj
            xobj = (trc_img['traces'][yidx,qq]-pixl_slit) / (pixr_slit-pixl_slit)
            # Generate 
            specobj = SpecObjExp((trc_img['object'].shape[:2]), config, scidx, det, (xl_slit,xr_slit),ypos, xobj, **kwargs)
            # Add traces
            specobj.trace = trc_img['traces'][:,qq]
            # Append
            specobjs.append(specobj)
    else:
        msgs.error("No way given to generate the list of SpecObjExp")
    # Return
    return specobjs



