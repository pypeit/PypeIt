# Module for QA in PYPIT
import os
import astropy.io.fits as pyfits
import armsgs as msgs
import numpy as np

import matplotlib
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages

from xastropy.xutils import xdebug as xdb

def trace_qa(slf, flat, lordloc, rordloc, root='trace', outfil=None):
    ''' Generate a QA plot for the traces
    Parameters:
    ------------
    flat: ndarray
      Flat field image
    lordloc: ndarray
      Left edge traces
    rordloc: ndarray
      Left edge traces
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
    ntrc = lordloc.shape[1]
    ycen = np.arange(flat.shape[0])
    # Normalize flux in the traces
    nrm_flat = np.zeros_like(flat)
    for ii in range(ntrc):
        xtrc = (lordloc[:,ii] + rordloc[:,ii])/2.
        ixtrc = np.round(xtrc).astype(int)
        # Simple 'extraction'
        dumi = np.zeros( (flat.shape[0],3) )
        for jj in range(3):
            dumi[:,jj] = flat[ycen,ixtrc-1+jj]
        trc = np.median(dumi, axis=1)
        # Find portion of the image and normalize
        for yy in ycen:
            xi = max(0,int(lordloc[yy,ii])-3)
            xe = min(flat.shape[1],int(rordloc[yy,ii])+3)
            # Fill + normalize
            nrm_flat[yy,xi:xe] = flat[yy,xi:xe] / trc[yy]

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
    mplt = plt.imshow(nrm_flat,origin='lower', cmap=cmm, extent=(0., flat.shape[1], 0., flat.shape[0]))
    mplt.set_clim(vmin=0.4, vmax=1.1)

    # Axes
    plt.xlim(0., flat.shape[1])
    plt.ylim(0., flat.shape[0])

    # Traces
    for ii in range(ntrc):
        # Left
        plt.plot(lordloc[:,ii], ycen, 'r--',alpha=0.7)
        # Right
        plt.plot(rordloc[:,ii], ycen, 'g--',alpha=0.7)
        # Label
        iy = int(flat.shape[0]/2.)
        plt.text(lordloc[iy,ii], ycen[iy], '{:d}'.format(ii+1), color='red', ha='center')
        plt.text(rordloc[iy,ii], ycen[iy], '{:d}'.format(ii+1), color='green', ha='center')

    pp.savefig(bbox_inches='tight')
    pp.close()
    plt.close()
