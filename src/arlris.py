# Module for LRIS spectific codes
import numpy as np
import glob
import armsgs as msgs
import astropy.io.fits as pyfits

import arcyextract
import arcyutils
import arcyproc
import arload
import artrace
import arutils
import arplot

from xastropy.xutils import xdebug as xdb

def read_lris(raw_file, TRIM=False):
    ''' Read a full raw LRIS data frame (both detectors)
    Packed in a multi-extension HDU
    Based on readmhdufits.pro

    Parameters:
    ----------
    raw_file: str
      Filename
    TRIM: bool, optional
      Trim the image?

    Returns:
    --------
    array: ndarray
      Combined image 
    header: FITS header
    '''
    # Check for file; allow for extra .gz, etc. suffix
    fil = glob.glob(raw_file+'*') 
    if len(fil) != 1:
        msgs.error("Found {:d} files matching {:s}".format(len(fil)))

    # Read
    msgs.info("Reading LRIS file: {:s}".format(fil[0]))
    hdu = pyfits.open(fil[0])
    head0 = hdu[0].header

    # Get post, pre-pix values
    precol = head0['PRECOL']
    postpix = head0['POSTPIX']
    preline = head0['PRELINE']
    postline = head0['POSTLINE']

    # get the x and y binning factors...
    binning = head0['BINNING']
    xbin,ybin = [int(ibin) for ibin in binning.split(',')]

    #; First read over the header info to determine the size of the output
    #; array...
    n_ext = len(hdu)
    xcol  = []
    xmax = 0
    ymax = 0
    xmin = 10000
    ymin = 10000
    for i in range(1,n_ext):
        theader = hdu[i].header
        detsec = theader['DETSEC']
        if detsec != '0':
            #;parse the DETSEC keyword to determine the size of the array.
            x1,x2,y1,y2 = np.array(arload.load_sections(detsec)).flatten()

            #; find the range of detector space occupied by the data 
            #; [xmin:xmax,ymin:ymax]
            xt = max(x2,x1)
            xmax = max(xt,xmax) 
            yt =  max(y2,y1)
            ymax = max(yt,ymax) 

            #; find the min size of the array
            xt = min(x1,x2)
            xmin = min(xmin,xt)
            yt = min(y1,y2)
            ymin = min(ymin,yt)
            # Save
            xcol.append(xt)

    #; determine the output array size...
    nx = xmax - xmin + 1
    ny = ymax - ymin + 1

    #; change size for binning...
    nx = nx / xbin
    ny = ny / ybin

    #; Update PRECOL and POSTPIX
    precol = precol / xbin
    postpix = postpix / xbin

    #; change size for pre/postscan...
    if not TRIM:
        nx += n_ext*(precol+postpix)
        ny += preline + postline

    #; allocate output array...
    array = np.zeros( (nx, ny) )
    order = np.argsort(np.array(xcol))

    #; insert extensions into master image...
    for i in range(1,n_ext):

        #; grab complete extension...
        data, predata, postdata, x1, y1 = lris_read_amp(hdu, i)
                            #, linebias=linebias, nobias=nobias, $
                            #x1=x1, x2=x2, y1=y1, y2=y2, gaindata=gaindata)
        #; insert components into output array...
        if not TRIM:
            #; insert predata...
            buf = predata.shape
            nxpre = buf[0]
            xs = order[i-1]*precol
            xe = xs + nxpre 
            '''
            if keyword_set(VERBOSE) then begin
                section = '['+stringify(xs)+':'+stringify(xe)+',*]'
                message, 'inserting extension '+stringify(i)+ $
                         ' predata  in '+section, /info
            endif 
            '''
            array[xs:xe,:] = predata

            #; insert data...
            buf = data.shape
            nxdata = buf[0]
            xs = n_ext*precol + (x1-xmin)/xbin
            xe = xs + nxdata 
            '''
            if keyword_set(VERBOSE) then begin
                section = '['+stringify(xs)+':'+stringify(xe)+',*]'

                message, 'inserting extension '+stringify(i)+ $
                         ' data     in '+section, /info
            endif 
            '''
            array[xs:xe,:] = data

            #; insert postdata...
            buf = postdata.shape
            nxpost = buf[0]
            xs = nx - n_ext*postpix + order[i-1]*postpix
            xe = xs + nxpost 
            '''
            if keyword_set(VERBOSE) then begin
                section = '['+stringify(xs)+':'+stringify(xe)+',*]'
                message, 'inserting extension '+stringify(i)+ $
                         ' postdata in '+section, /info
            endif 
            '''
            array[xs:xe,:] = postdata
        else:
            buf = data.shape
            nxdata = buf[0]
            nydata = buf[1]

            xs = (x1-xmin)/xbin
            xe = xs + nxdata 
            ys = (y1-ymin)/ybin
            ye = ys + nydata - postline

            yin1 = preline
            yin2 = nydata - postline 

            '''
            if keyword_set(VERBOSE) then begin
                section = '['+stringify(xs)+':'+stringify(xe)+ $
                          ','+stringify(ys)+':'+stringify(ye)+']'
                message, 'inserting extension '+stringify(i)+ $
                         ' data     in '+section, /info
            endif 
            '''
           
            array[xs:xe,ys:ye] = data[:,yin1:yin2]

    #; make sure BZERO is a valid integer for IRAF
    obzero = head0['BZERO']
    head0['O_BZERO'] = obzero
    head0['BZERO'] = 32768-obzero

    # Return, transposing array back to goofy Python indexing
    return array.transpose(), head0

def lris_read_amp(inp, ext):
    ''' Read one amplifier of an LRIS multi-extension FITS image

    Parameters:
    ----------
    inp: tuple 
      (str,int) filename, extension
      (hdu,int) FITS hdu, extension

    Returns:
    ----------
    data
    predata
    postdata
    x1
    y1

    ;------------------------------------------------------------------------
    function lris_read_amp, filename, ext, $
      linebias=linebias, nobias=nobias, $
      predata=predata, postdata=postdata, header=header, $
      x1=x1, x2=x2, y1=y1, y2=y2, GAINDATA=gaindata
    ;------------------------------------------------------------------------
    ; Read one amp from LRIS mHDU image
    ;------------------------------------------------------------------------
    '''
    # Parse input
    if isinstance(inp,basestring):
        hdu = fits.open(inp)
    else:
        hdu = inp

    #; Get the pre and post pix values
    #; for LRIS red POSTLINE = 20, POSTPIX = 80, PRELINE = 0, PRECOL = 12
    head0 = hdu[0].header
    precol   = head0['precol']
    postpix  = head0['postpix']

    #; Deal with binning
    binning = head0['BINNING']
    xbin,ybin = [int(ibin) for ibin in binning.split(',')]
    precol = precol/xbin
    postpix = postpix/xbin

    #; get entire extension...
    temp = hdu[ext].data.transpose() # Silly Python nrow,ncol formatting
    tsize = temp.shape
    nxt = tsize[0]

    #; parse the DETSEC keyword to determine the size of the array.
    header = hdu[ext].header
    detsec = header['DETSEC']
    x1,x2,y1,y2 = np.array(arload.load_sections(detsec)).flatten()

    #; parse the DATASEC keyword to determine the size of the science region
    datasec = header['DATASEC']
    xdata1,xdata2,ydata1,ydata2 = np.array(arload.load_sections(datasec)).flatten()

    #; grab the components...
    predata  = temp[0:precol-1,:]
    # datasec appears to have the x value for the keywords that are zero
    # based. This is only true in the image header extensions
    # not true in the main header.
    #data     = temp[xdata1-1:xdata2-1,*]
    data     = temp[xdata1:xdata2,:]
    postdata = temp[nxt-postpix:nxt-1,:]

    #; flip in X as needed...
    if x1 > x2:
        xt=x2
        x2=x1
        x1=xt
        data = np.flipud(data) #reverse(temporary(data),1)

    #; flip in Y as needed...
    if y1 > y2:
        yt=y2
        y2=y1
        y1=yt
        data  = np.fliplr(data)
        predata  = np.fliplr(predata)
        postdata = np.fliplr(postdata)

    '''
    #; correct gain if requested...
    if keyword_set(GAINDATA) then begin
        gain = gainvalue( gaindata, header)
        data = FLOAT(temporary(data)) * gain
        predata = FLOAT(temporary(predata)) * gain
        postdata = FLOAT(temporary(postdata)) * gain
    endif
    '''

    '''
    ;; optional bias subtraction...
    if ~ keyword_set(NOBIAS) then begin
        if keyword_set( LINEBIAS) then begin
            ;; compute a bias for each line...
            bias = median( postdata, dim=1)

            ;; subtract for data...
            buf = size(data)
            nx = buf[1]
            ny = buf[2]
            data2 = fltarr(nx,ny)
            for i=0,nx-1 do begin
                data2[i,*] = float(data[i,*]) - bias
            endfor 
            data = data2
        endif else begin
            ;; compute a scalar bias....
            bias = median( postdata)
            data -= bias
        endelse
    endif
    '''

    return data, predata, postdata, x1, y1

