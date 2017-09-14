'''
Implements DEIMOS-specific functions, including reading in slitmask design files.
'''
from __future__ import absolute_import, division, print_function

import glob
import numpy as np
#from astropy.io import fits
import astropy.io.fits as pyfits


from pypit import armsgs
from pypit.arparse import load_sections
from pypit import ardebug as debugger
from IPython import embed


# Logging
msgs = armsgs.get_logger()

def read_deimos(raw_file, det, trim=False):
    """
    Read a raw DEIMOS data frame (one or more detectors)
    Packed in a multi-extension HDU
    Based on pypit.arlris.read_lris...
       Based on readmhdufits.pro

    Parameters
    ----------
    raw_file : str
      Filename
    det : int
      detector index starting at 1
    trim : bool, optional
      Trim the image?

    Returns
    -------
    array : ndarray
      Combined image 
    header : FITS header
    sections : list
      List of datasec, oscansec, ampsec sections
    """

    # Check for file; allow for extra .gz, etc. suffix
    fil = glob.glob(raw_file+'*') 
    if len(fil) != 1:
        msgs.error("Found {:d} files matching {:s}".format(len(fil)))

    # Read
    msgs.info("Reading DEIMOS file: {:s}".format(fil[0]))
    hdu = pyfits.open(fil[0])
    head0 = hdu[0].header

    # Get post, pre-pix values
    precol = head0['PRECOL']
    postpix = head0['POSTPIX']
    preline = head0['PRELINE']
    postline = head0['POSTLINE']

    # Setup for datasec, oscansec
    dsec = []
    osec = []

    # get the x and y binning factors...
    binning = head0['BINNING']
    xbin, ybin = [int(ibin) for ibin in binning.split(',')]

    # First read over the header info to determine the size of the output array...
    #n_ext = len(hdu)-1  # Number of extensions (usually 4)
    n_ext = 8
    xcol = []
    xmax = -np.inf
    ymax = -np.inf
    xmin = np.inf
    ymin = np.inf
    for i in np.arange(1, n_ext+1):
        theader = hdu[i].header
        detsec = theader['DETSEC']
        if detsec != '0':
            # parse the DETSEC keyword to determine the size of the array.
            x1, x2, y1, y2 = np.array(load_sections(detsec)).flatten()

            # find the range of detector space occupied by the data
            # [xmin:xmax,ymin:ymax]
            xt = max(x2, x1)
            xmax = max(xt, xmax)
            yt =  max(y2, y1)
            ymax = max(yt, ymax)

            # find the min size of the array
            xt = min(x1, x2)
            xmin = min(xmin, xt)
            yt = min(y1, y2)
            ymin = min(ymin, yt)
            # Save
            xcol.append(xt)

    # determine the output array size...
    nx = xmax - xmin + 1
    ny = ymax - ymin + 1

    # change size for binning...
    nx = nx // xbin
    ny = ny // ybin

    # Update PRECOL and POSTPIX
    precol = precol // xbin
    postpix = postpix // xbin

    # Deal with detectors
    assert det in np.arange(1, 9)
    det_idx = det - 1
    ndet = 1
    
    # if det in [1,2]:
    #     nx = nx // 2
    #     n_ext = n_ext // 2
    #     det_idx = np.arange(n_ext, dtype=np.int) + (det-1)*n_ext
    #     ndet = 1
    # elif det is None:
    #     ndet = 2
    #     det_idx = np.arange(n_ext).astype(int)
    # else:
    #     raise ValueError('Bad value for det')

    # change size for pre/postscan...
    if not trim:
        nx += n_ext*(precol+postpix)
        ny += preline + postline

    # allocate output array...
    array = np.zeros( (nx, ny) )
    order = np.argsort(np.array(xcol))


    # insert extensions into master image...
    for kk, i in enumerate(order[det_idx]):

        # grab complete extension...
        data, predata, postdata, x1, y1 = deimos_read_amp(hdu, i+1)
                            #, linebias=linebias, nobias=nobias, $
                            #x1=x1, x2=x2, y1=y1, y2=y2, gaindata=gaindata)
        # insert components into output array...
        if not trim:
            # insert predata...
            buf = predata.shape
            nxpre = buf[0]
            xs = kk*precol
            xe = xs + nxpre
            '''
            if keyword_set(VERBOSITY) then begin
                section = '['+stringify(xs)+':'+stringify(xe)+',*]'
                message, 'inserting extension '+stringify(i)+ $
                         ' predata  in '+section, /info
            endif 
            '''
            array[xs:xe, :] = predata

            # insert data...
            buf = data.shape
            nxdata = buf[0]
            nydata = buf[1]
            xs = n_ext*precol + kk*nxdata #(x1-xmin)/xbin
            xe = xs + nxdata
            # Data section
            section = '[{:d}:{:d},{:d}:{:d}]'.format(preline,nydata-postline, xs, xe)  # Eliminate lines
            dsec.append(section)
            #print('data',xs,xe)
            array[xs:xe, :] = data   # Include postlines

            #; insert postdata...
            buf = postdata.shape
            nxpost = buf[0]
            xs = nx - n_ext*postpix + kk*postpix
            xe = xs + nxpost 
            section = '[:,{:d}:{:d}]'.format(xs, xe)
            osec.append(section)
            '''
            if keyword_set(VERBOSITY) then begin
                section = '['+stringify(xs)+':'+stringify(xe)+',*]'
                message, 'inserting extension '+stringify(i)+ $
                         ' postdata in '+section, /info
            endif 
            '''
            array[xs:xe, :] = postdata
        else:
            buf = data.shape
            nxdata = buf[0]
            nydata = buf[1]

            xs = (x1-xmin)//xbin
            xe = xs + nxdata 
            ys = (y1-ymin)//ybin
            ye = ys + nydata - postline

            yin1 = preline
            yin2 = nydata - postline 

            '''
            if keyword_set(VERBOSITY) then begin
                section = '['+stringify(xs)+':'+stringify(xe)+ $
                          ','+stringify(ys)+':'+stringify(ye)+']'
                message, 'inserting extension '+stringify(i)+ $
                         ' data     in '+section, /info
            endif 
            '''
            array[xs:xe, ys:ye] = data[:, yin1:yin2]

    # make sure BZERO is a valid integer for IRAF
    obzero = head0['BZERO']
    head0['O_BZERO'] = obzero
    head0['BZERO'] = 32768-obzero

    # Return, transposing array back to goofy Python indexing
    return array.T, head0, (dsec, osec)


def deimos_read_amp(inp, ext):
    """
    Read one amplifier of a DEIMOS multi-extension FITS image
    Copied from arlris.lris_read_amp

    TODO: review this function, is it necessary?

    Parameters
    ----------
    inp: tuple 
      (str,int) filename, extension
      (hdu,int) FITS hdu, extension

    Returns
    -------
    data
    predata
    postdata
    x1
    y1
    """
    # Parse input
    if isinstance(inp, basestring):
        hdu = pyfits.open(inp)
    else:
        hdu = inp

    # Get the pre and post pix values
    # for LRIS red POSTLINE = 20, POSTPIX = 80, PRELINE = 0, PRECOL = 12
    head0 = hdu[0].header
    precol = head0['precol']
    postpix = head0['postpix']

    # Deal with binning
    binning = head0['BINNING']
    xbin, ybin = [int(ibin) for ibin in binning.split(',')]
    precol = precol//xbin
    postpix = postpix//xbin

    # get entire extension...
    temp = hdu[ext].data.transpose() # Silly Python nrow,ncol formatting
    tsize = temp.shape
    nxt = tsize[0]
    # parse the DETSEC keyword to determine the size of the array.
    header = hdu[ext].header
    detsec = header['DETSEC']
    x1, x2, y1, y2 = np.array(load_sections(detsec)).flatten()
    # parse the DATASEC keyword to determine the size of the science region (unbinned)
    datasec = header['DATASEC']
    xdata1, xdata2, ydata1, ydata2 = np.array(load_sections(datasec)).flatten()
    # grab the components...
    predata = temp[0:precol, :]
    # datasec appears to have the x value for the keywords that are zero
    # based. This is only true in the image header extensions
    # not true in the main header.  They also appear inconsistent between
    # LRISr and LRISb!
    #data     = temp[xdata1-1:xdata2-1,*]
    #data = temp[xdata1:xdata2+1, :]
    if (xdata1-1) != precol:
        msgs.error("Something wrong in LRIS datasec or precol")
    xshape = 1024 // xbin
    if (xshape+precol+postpix) != temp.shape[0]:
        msgs.error("Wrong size for in LRIS detector somewhere.  Funny binning?")
    data = temp[precol:precol+xshape,:]
    postdata = temp[nxt-postpix:nxt, :]
    # flip in X as needed...
    if x1 > x2:
        xt = x2
        x2 = x1
        x1 = xt
        data = np.flipud(data) #reverse(temporary(data),1)
    # flip in Y as needed...
    if y1 > y2:
        yt = y2
        y2 = y1
        y1 = yt
        data = np.fliplr(data)
        predata = np.fliplr(predata)
        postdata = np.fliplr(postdata)
    return data, predata, postdata, x1, y1


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
            hdulist = pyfits.open(hdulist)
        else:
            assert isinstance(hdulist, pyfits.hdu.hdulist.HDUList)
        
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
        

    
