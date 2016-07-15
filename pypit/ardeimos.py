'''
Implements DEIMOS-specific functions, including reading in slitmask design files and lamp ID.
'''
from __future__ import absolute_import, division, print_function

import numpy as np
from astropy.io import fits
from pypit import armsgs
# import armsgs
try:
    from xastropy.xutils import xdebug as debugger
except:
    import pdb as debugger

# Logging
msgs = armsgs.get_logger()

class DEIMOS_slits(object):
    '''
    DEIMOS slit mask design.  Coordinates are stored in pixel units on the 
    detector, such that the x-axis is along the spatial axis and the y-axis is 
    along the dispersion axis.

    DEIMOS multi-extension fits files should have at least the following named 
    binary tables:
        MaskDesign - a single line containing the mask name and other details
        ObjectCat - a record of accepted objects from the DSIMULATOR input file
        DesiSlits - design parameters for each slit
        SlitObjMap - a record of the mapping between slit IDs and object IDs
        BluSlits - contains the slit corners in arcseconds away from the center
    '''

    def __init__(self, hdulist, plate_scale=0.1185):
        '''
        Parameters
        ----------
        hdulist : either a string, in which case it corresponds to 
                  the fits filename, or an HDUList instance from
                  astropy's fits module
        plate_scale : the DEIMOS plate scale, arcseconds per pixel
        '''

        if isinstance(hdulist, str):
            hdulist = fits.open(hdulist)
        else:
            assert isinstance(hdulist, fits.hdu.hdulist.HDUList)
        
        # mask attributes
        self.plate_scale = plate_scale
        
        error_str = ' bin table not found, verify the FITS file is in the correct format!'
        try:
            mask_design = hdulist['MaskDesign'].data
        except KeyError:
            msgs.error('MaskDesign' + error_str)
        try:
            obj_cat = hdulist['ObjectCat'].data
        except KeyError:
            msgs.error('ObjectCat' + error_str)
        try:
            desi_slits = hdulist['DesiSlits'].data
        except KeyError:
            msgs.error('DesiSlits' + error_str)
        try:
            slit_obj_map = hdulist['SlitObjMap'].data
        except KeyError:
            msgs.error('SlitObjMap' + error_str)
        try:
            blu_slits = hdulist['BluSlits'].data
        except KeyError:
            msgs.error('BluSlits' + error_str)
            
        self.mask_name = mask_design[0]['DesName']
        self.mask_pa = mask_design[0]['PA_PNT']        # in degrees ccw of North
        self.nslits = len(desi_slits)

        # slit attributes, each is an array of size nslits
        self.slit_ra = desi_slits['slitRA']            # in degrees
        self.slit_dec = desi_slits['slitDec']          # in degrees
        self.slit_pa = desi_slits['slitLPA']           # in degrees ccw of North
        self.slit_length = desi_slits['slitLen']       # in arcsec
        self.slit_width = desi_slits['slitWid']        # in arcsec
        self.slit_isalign = (desi_slits['slitTyp'] == 'A')  # boolean array, true if alignment slit
        # assign slit names based on the DSIM input names
        assert len(desi_slits) == len(slit_obj_map)
        idx = np.array([(i in slit_obj_map['ObjectId']) for i in obj_cat['ObjectId']])
        self.slit_name = obj_cat[idx]['OBJECT']

        # [b]ottom/[t]op, [l]eft/[r]ight corners, given as [nslits x 2] arrays
        # with (x, y) coordinates in the latter axis.
        # Coordinates are in arcseconds relative the mask center, with a
        # left-handed orientation, such that increasing y moves up, increasing
        # x moves left.
        self.slit_brc = np.array(zip(blu_slits['SlitX1'], blu_slits['SlitY1']))
        self.slit_blc = np.array(zip(blu_slits['SlitX2'], blu_slits['SlitY2']))
        self.slit_tlc = np.array(zip(blu_slits['SlitX3'], blu_slits['SlitY3']))
        self.slit_trc = np.array(zip(blu_slits['SlitX4'], blu_slits['SlitY4']))

    def __repr__(self):
        return '<DEIMOS_slits: ' + self.mask_name + ' (' + str(self.nslits) + ' slits)>'

    def __len__(self):
        return self.nslits

    def get_edges(self, det):
        '''
        Calculates the pixel values of slit edges for the associated detector.

        The DEIMOS chip layout is shown in the engineering drawings found here:
        http://www.ucolick.org/~sla/fits/mosaic/d0307j.pdf

        Parameters
        ----------
        det : int, detector index, measured from 1, matching the 'detXX' syntax 
              in the settings.deimos file

        Returns
        -------
        edges : [n by 2] array, where n is the number of slits on the detector
                Each row in the array gives the left and right pixel values of
                the slit.
        '''
        pass

    
