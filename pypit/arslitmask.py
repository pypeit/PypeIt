'''
Module for dealing with multi-object spectroscopy slitmasks
'''
from __future__ import (print_function, absolute_import,
                        division, unicode_literals)

import numpy as np
from astropy.io import fits as pyfits

from pypit import armsgs

try:
    from xastropy.xutils import xdebug as debugger
except:
    import pdb as debugger

# Logging
msgs = armsgs.get_logger()

class Slit(object):
    '''
    Slit with mask coordinates
    '''

    def __init__(self, left_edge, right_edge, name):
        '''
        Parameters
        ----------
        left_edge : float, 
        right_edge : 
        name : str, name of slit
        '''
        self.left_edge = left_edge
        self.right_edge = right_edge
        self.name = name

    def __repr__(self):
        return '<' + self.__name__ + ': ' + self.name + '>'

        
class Slitmask(object):
    '''
    Generic slitmask class, should be sub-classed for specific instruments.
    '''

    def __init__(self):
        self.mask_name = None
        self.mask_pa = None
        self.slits = []

    def __len__(self):
        return len(self.slits)
    
    def __repr__(self):
        s = ('<' + self.__name__ + ': ' + self.mask_name +  ' (' +
             str(len(self)) + ' slits)>')

class DEIMOS_slitmask(Slitmask):
    pass

class LRIS_slitmask(Slitmask):
    pass

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

    def __init__(self, hdulist, plate_scale=0.1185, pix_per_chip=2048):
        '''
        Parameters
        ----------
        hdulist : either a string, in which case it corresponds to 
                  the fits filename, or an HDUList instance from
                  astropy's fits module
        plate_scale : the DEIMOS plate scale, arcseconds per pixel
        pix_per_chip : number of pixels along the spatial direction of CCDs
        '''

        if isinstance(hdulist, str):
            hdulist = fits.open(hdulist)
        else:
            assert isinstance(hdulist, fits.hdu.hdulist.HDUList)
        
        # mask attributes
        self.plate_scale = plate_scale
        self.pix_per_chip = pix_per_chip
        
        error_str = (' bin table not found, verify the FITS file '
                     'is in the correct format!')
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
        self.slit_isalign = (desi_slits['slitTyp'] == 'A')  # boolean array
        # assign slit names based on the DSIM input names
        assert len(desi_slits) == len(slit_obj_map)
        idx = np.array([(i in slit_obj_map['ObjectId'])
                        for i in obj_cat['ObjectId']])
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

        # check that slits have paired left and paired right corneres
        assert np.all(np.isclose(self.slit_brc[:, 0],
                                 self.slit_trc[:, 0], atol=0.01))
        assert np.all(np.isclose(self.slit_blc[:, 0],
                                 self.slit_tlc[:, 0], atol=0.01))
        
    def __repr__(self):
        return ('<DEIMOS_slits: ' + self.mask_name +
                ' (' + str(self.nslits) + ' slits)>')

    def __len__(self):
        return self.nslits

    def mosaic_location(self, det):
        '''
        Locates the chip within the DEIMOS detector mosaic.

        Parameters
        ----------
        det : int, detector index, measured from 1, matching the 'detXX' syntax 
              in the settings.deimos file

        Returns
        -------
        row : int, bottow (0) is blue side, top (1) is red side
        col : int, left to right
        '''
        # check that the requested detector(s) make(s) sense
        try:
            good = det in range(1, 9)
        except ValueError:
            # det should be iterable
            good = all([d in range(1, 9) for d in det])
            det = np.array(det)
        assert good
        return divmod(det - 1, 4)

    def mask_to_det_coords(self, x, det):
        '''
        Transform mask spatial arcsec coordinates to pixel coordinates of the 
        detector.

        Parameters
        ----------
        x : float, arcsec in mask coord, can be array
        det : int, detector index, measured from 1, matching the 'detXX' syntax 
              in the settings.deimos file

        Returns
        -------
        pix_x : float, detector pixel value of mask coord, will be array if 
                input is an array
        '''
        row, col = self.mosaic_location(det)
        
        # bottom row detectors are oriented with origin in blc right-handed
        pix_x = x / self.plate_scale - (col - 2) * self.pix_per_chip
        if row == 1:
            # top row detectors are oriented with origin in trc, right-handed
            pix_x = self.pix_per_chip - pix_x
        return pix_x

    def det_to_mask_coords(self, pix_x, det):
        '''

        Transforms detector pixel coordinate to mask coordinate in arcsec.

        Parameters
        ----------
        pix_x : float, detector pixel value of mask coord, can be array
        det : int, detector index, measured from 1, matching the 'detXX' syntax 
              in the settings.deimos file

        Returns
        -------
        x : float, arcsec in mask coord, will be array if pix_x is one
        '''

        row, col = self.mosaic_location(det)
        
        if row == 1:
            pix_x = self.pix_per_chip - pix_x
        x = (pix_x + (col - 2) * self.pix_per_chip) * self.plate_scale
        return x
            
    def slits_in_det(self, det):
        '''
        Gets a boolean array corresponding to slits falling on the specified
        detector.

        Parameters
        ----------
        det : int, detector index, measured from 1, matching the 'detXX' syntax 
              in the settings.deimos file
        
        Returns
        -------
        slits : boolean array of size nslits
        '''
        row, col = self.mosaic_location(det)

        if row == 0:
            left = 0
            right = self.pix_per_chip
        else:
            left = self.pix_per_chip
            right = 0
        chip_left_edge = self.det_to_mask_coords(left, det)
        chip_right_edge = self.det_to_mask_coords(right, det)

        left_edges = self.slit_blc[:, 0]
        right_edges = self.slit_brc[:, 0]
        slits = (chip_left_edge < left_edges) & (right_edges < chip_right_edge)
        return slits
        

    def get_edges(self, det, sort=True):
        '''
        Calculates the pixel values of slit edges for the associated detector.

        The DEIMOS chip layout is shown in the engineering drawings found here:
        http://www.ucolick.org/~sla/fits/mosaic/d0307j.pdf

        Parameters
        ----------
        det : int, detector index, measured from 1, matching the 'detXX' syntax 
              in the settings.deimos file
        sort : boolean, if true, sort in ascending order

        Returns
        -------
        edges : [n by 2] array, where n is the number of slits on the detector
                Each row in the array gives the left and right pixel values of
                the slit.
        '''
        slits = self.slits_in_det(det)

        left_edges = self.mask_to_det_coords(self.slit_blc[:, 0][slits], det)
        right_edges = self.mask_to_det_coords(self.slit_brc[:, 0][slits], det)
        names = self.slit_name[slits]
        
        if sort:
            idx = np.argsort(left_edges)
            left_edges = left_edges[idx]
            right_edges = right_edges[idx]
            names = names[idx]
            
        return names, np.column_stack([left_edges, right_edges])
        

    
