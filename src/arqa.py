# Module for QA in PYPIT
import os
import astropy.io.fits as pyfits
import armsgs as msgs
import numpy as np
from arplot import zscale

import matplotlib
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages

try:
    from xastropy.xutils import xdebug as xdb
except:
    pass

def trace_qa(slf, frame, ltrace, rtrace, root='trace', outfil=None, normalize=True):
    ''' Generate a QA plot for the traces
    Parameters:
    ------------
    frame: ndarray
      image
    ltrace: ndarray
      Left edge traces
    rtrace: ndarray
      Right edge traces
    root: str, optional
      Root name for generating output file, e.g. msflat_01blue_000.fits
    outfil: str, optional
      Output file
    '''
    # Outfil
    if outfil is None:
        if 'fits' in root: # Expecting name of msflat FITS file
            outfil = root.replace('.fits', '_trc.pdf')
            outfil = outfil.replace('MasterFrames', 'Plots')
        else:
            outfil = root+'.pdf'
    ntrc = ltrace.shape[1]
    ycen = np.arange(frame.shape[0])
    # Normalize flux in the traces
    if normalize:
        nrm_frame = np.zeros_like(frame)
        for ii in xrange(ntrc):
            xtrc = (ltrace[:,ii] + rtrace[:,ii])/2.
            ixtrc = np.round(xtrc).astype(int)
            # Simple 'extraction'
            dumi = np.zeros( (frame.shape[0],3) )
            for jj in xrange(3):
                dumi[:,jj] = frame[ycen,ixtrc-1+jj]
            trc = np.median(dumi, axis=1)
            # Find portion of the image and normalize
            for yy in ycen:
                xi = max(0,int(ltrace[yy,ii])-3)
                xe = min(frame.shape[1],int(rtrace[yy,ii])+3)
                # Fill + normalize
                nrm_frame[yy,xi:xe] = frame[yy,xi:xe] / trc[yy]
        sclmin, sclmax = 0.4, 1.1
    else:
        nrm_frame = frame.copy()
        sclmin, sclmax = zscale(nrm_frame)

    # Plot
    pp = PdfPages(outfil)
    plt.clf()
    fig = plt.figure(dpi=1200)
    #fig.set_size_inches(10.0,6.5)

    plt.rcParams['font.family']= 'times new roman'
    ticks_font = matplotlib.font_manager.FontProperties(family='times new roman', 
       style='normal', size=16, weight='normal', stretch='normal')
    ax = plt.gca()
    for label in ax.get_yticklabels() :
        label.set_fontproperties(ticks_font)
    for label in ax.get_xticklabels() :
        label.set_fontproperties(ticks_font)
    cmm = cm.Greys_r
    mplt = plt.imshow(nrm_frame,origin='lower', cmap=cmm, extent=(0., frame.shape[1]-1, 0., frame.shape[0]-1))
    mplt.set_clim(vmin=sclmin, vmax=sclmax)

    # Axes
    plt.xlim(0., frame.shape[1]-1)
    plt.ylim(0., frame.shape[0]-1)

    # Traces
    for ii in xrange(ntrc):
        # Left
        plt.plot(ltrace[:,ii], ycen, 'r--',alpha=0.7)
        # Right
        plt.plot(rtrace[:,ii], ycen, 'g--',alpha=0.7)
        # Label
        iy = int(frame.shape[0]/2.)
        plt.text(ltrace[iy,ii], ycen[iy], '{:d}'.format(ii+1), color='red', ha='center')
        plt.text(rtrace[iy,ii], ycen[iy], '{:d}'.format(ii+1), color='green', ha='center')

    pp.savefig(bbox_inches='tight')
    pp.close()
    plt.close()
