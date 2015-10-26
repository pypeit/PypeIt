# Module for QA in PYPIT
import os
import astropy.io.fits as pyfits
import armsgs as msgs
import arutils
import numpy as np

import matplotlib
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages

try:
    from xastropy.xutils import xdebug as xdb
except:
    pass


def arc_fit_qa(fit, arc_spec, outroot=None, outfil=None):
    '''QA for Arc spectrum
    Parameters:
    -----------
    outfil: str, optional
      Name of output file
    '''
    import matplotlib.gridspec as gridspec
    if outfil is None:
        if outroot is None:
            outfil = 'Plots/arc_qa.pdf'
        else:
            outfil = outroot.replace('.fits','_fit.pdf')
            outfil = outfil.replace('MasterFrames', 'Plots')

    # Begin
    pp = PdfPages(outfil)
    plt.figure(figsize=(8, 4.0))
    plt.clf()
    gs = gridspec.GridSpec(2,2)

    # Simple spectrum plot
    ax_spec = plt.subplot(gs[:,0])
    ax_spec.plot(np.arange(len(arc_spec)), arc_spec)
    ymin, ymax = 0., np.max(arc_spec)
    ysep = ymax*0.03
    for kk,x in enumerate(fit['xfit']*fit['xnorm']):
        yline = np.max(arc_spec[int(x)-2:int(x)+2])
        # Tick mark
        ax_spec.plot([x,x], [yline+ysep*0.25, yline+ysep], 'g-')
        # label
        ax_spec.text(x, yline+ysep*1.3, 
            '{:g}'.format(fit['yfit'][kk]), ha='center', va='bottom',
            size='xx-small', rotation=90., color='green')
    ax_spec.set_xlim(0., len(arc_spec))
    ax_spec.set_ylim(ymin, ymax*1.2)
    ax_spec.set_xlabel('Pixel')
    ax_spec.set_ylabel('Flux')

    # Arc Fit
    ax_fit = plt.subplot(gs[0,1])
    # Points
    ax_fit.scatter(fit['xfit']*fit['xnorm'], fit['yfit'], marker='x')
    if len(fit['xrej']) > 0:
        ax_fit.scatter(fit['xrej']*fit['xnorm'], fit['yrej'], marker='o',
            edgecolor='gray', facecolor='none')
    # Solution
    xval = np.arange(len(arc_spec))
    wave = arutils.func_val(fit['fitc'], xval/fit['xnorm'], 'legendre', 
        min=fit['fmin'], max=fit['fmax'])
    ax_fit.plot(xval, wave, 'r-')
    xmin, xmax = 0., len(arc_spec)
    ax_fit.set_xlim(xmin, xmax)
    ymin,ymax = np.min(wave)*.95,  np.max(wave)*1.05
    ax_fit.set_ylim(np.min(wave)*.95,  np.max(wave)*1.05)
    ax_fit.set_ylabel('Wavelength')
    ax_fit.get_xaxis().set_ticks([]) # Suppress labeling
    # Stats
    wave_fit = arutils.func_val(fit['fitc'], fit['xfit'], 'legendre', 
        min=fit['fmin'], max=fit['fmax'])
    dwv_pix = np.median(np.abs(wave-np.roll(wave,1)))
    rms = np.sqrt(np.sum((fit['yfit']-wave_fit)**2)/len(fit['xfit'])) # Ang
    ax_fit.text(0.1*len(arc_spec), 0.90*ymin+(ymax-ymin),
        r'$\Delta\lambda$={:.3f}$\AA$ (per pix)'.format(dwv_pix), size='small')
    ax_fit.text(0.1*len(arc_spec), 0.80*ymin+(ymax-ymin),
        'RMS={:.3f} (pixels)'.format(rms/dwv_pix), size='small')
    # Arc Residuals
    ax_res = plt.subplot(gs[1,1])
    res = fit['yfit']-wave_fit
    ax_res.scatter(fit['xfit']*fit['xnorm'], res/dwv_pix, marker='x')
    ax_res.plot([xmin,xmax], [0.,0], 'k--')
    ax_res.set_xlim(xmin, xmax)
    ax_res.set_xlabel('Pixel')
    ax_res.set_ylabel('Residuals (Pix)')

    # Finish
    plt.tight_layout(pad=0.2,h_pad=0.0,w_pad=0.0)
    pp.savefig(bbox_inches='tight')
    pp.close()
    plt.close()

def slit_trace_qa(slf, flat, lordloc, rordloc, root='trace', outfil=None):
    ''' Generate a QA plot for the slit traces
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
