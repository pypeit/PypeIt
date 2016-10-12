""" Module for QA in PYPIT
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

from astropy import units as u

import os
import numpy as np
from pypit.arplot import zscale
from pypit import armsgs
from pypit import arutils

import matplotlib
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages

msgs = armsgs.get_logger()

plt.rcParams['font.family']= 'times new roman'
ticks_font = matplotlib.font_manager.FontProperties(family='times new roman',
                                                    style='normal', size=16, weight='normal', stretch='normal')

try:
    from xastropy.xutils import xdebug as debugger
except:
    import pdb as debugger

def arc_fit_qa(slf, fit, outfil=None, ids_only=False, title=None):
    """
    QA for Arc spectrum

    Parameters
    ----------
    fit : Wavelength fit
    arc_spec : ndarray
      Arc spectrum
    outfil : str, optional
      Name of output file
    """
    arc_spec = fit['spec']
    if outfil is not None:
        pp = PdfPages(outfil)

    # Begin
    if not ids_only:
        plt.figure(figsize=(8, 4.0))
        plt.clf()
        gs = gridspec.GridSpec(2, 2)
        idfont = 'xx-small'
    else:
        plt.figure(figsize=(11, 8.5))
        plt.clf()
        gs = gridspec.GridSpec(1, 1)
        idfont = 'small'

    # Simple spectrum plot
    ax_spec = plt.subplot(gs[:,0])
    ax_spec.plot(np.arange(len(arc_spec)), arc_spec)
    ymin, ymax = 0., np.max(arc_spec)
    ysep = ymax*0.03
    for kk, x in enumerate(fit['xfit']*fit['xnorm']):
        yline = np.max(arc_spec[int(x)-2:int(x)+2])
        # Tick mark
        ax_spec.plot([x,x], [yline+ysep*0.25, yline+ysep], 'g-')
        # label
        ax_spec.text(x, yline+ysep*1.3, 
            '{:s} {:g}'.format(fit['ions'][kk], fit['yfit'][kk]), ha='center', va='bottom',
            size=idfont, rotation=90., color='green')
    ax_spec.set_xlim(0., len(arc_spec))
    ax_spec.set_ylim(ymin, ymax*1.2)
    ax_spec.set_xlabel('Pixel')
    ax_spec.set_ylabel('Flux')
    if title is not None:
        ax_spec.text(0.04, 0.93, title, transform=ax_spec.transAxes,
                     size='x-large', ha='left')#, bbox={'facecolor':'white'})
    if ids_only:
        plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
        if outfil is not None:
            pp.savefig(bbox_inches='tight')
            pp.close()
        plt.close()
        return

    # Arc Fit
    ax_fit = plt.subplot(gs[0, 1])
    # Points
    ax_fit.scatter(fit['xfit']*fit['xnorm'], fit['yfit'], marker='x')
    if len(fit['xrej']) > 0:
        ax_fit.scatter(fit['xrej']*fit['xnorm'], fit['yrej'], marker='o',
            edgecolor='gray', facecolor='none')
    # Solution
    xval = np.arange(len(arc_spec))
    wave = arutils.func_val(fit['fitc'], xval/fit['xnorm'], 'legendre', 
        minv=fit['fmin'], maxv=fit['fmax'])
    ax_fit.plot(xval, wave, 'r-')
    xmin, xmax = 0., len(arc_spec)
    ax_fit.set_xlim(xmin, xmax)
    ymin,ymax = np.min(wave)*.95,  np.max(wave)*1.05
    ax_fit.set_ylim(np.min(wave)*.95,  np.max(wave)*1.05)
    ax_fit.set_ylabel('Wavelength')
    ax_fit.get_xaxis().set_ticks([]) # Suppress labeling
    # Stats
    wave_fit = arutils.func_val(fit['fitc'], fit['xfit'], 'legendre', 
        minv=fit['fmin'], maxv=fit['fmax'])
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
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
    slf._qa.savefig(bbox_inches='tight')
    plt.close()
    return


def flexure(slf, det, flex_dict, slit_cen=False):
    """ QA on flexure measurement

    Parameters
    ----------
    slf
    det
    flex_dict
    slit_cen : bool, optional
      QA on slit center instead of objects

    Returns
    -------

    """
    # Setup
    if slit_cen:
        nobj = 1
        ncol = 1
    else:
        nobj = len(slf._specobjs[det-1])
        if nobj == 0:
            return
        ncol = min(3,nobj)
    #
    nrow = nobj // ncol + ((nobj%ncol) > 0)


    plt.figure(figsize=(8, 5.0))
    plt.clf()
    gs = gridspec.GridSpec(nrow, ncol)

    # Correlation QA
    for o in range(nobj):
        ax = plt.subplot(gs[o//ncol, o%ncol])
        # Fit
        fit = flex_dict['polyfit'][o]
        xval = np.linspace(-10., 10, 100) + flex_dict['corr_cen'][o] #+ flex_dict['shift'][o]
        #model = (fit[2]*(xval**2.))+(fit[1]*xval)+fit[0]
        model = arutils.func_val(fit, xval, 'polynomial')
        mxmod = np.max(model)
        ylim = [np.min(model/mxmod), 1.3]
        ax.plot(xval-flex_dict['corr_cen'][o], model/mxmod, 'k-')
        # Measurements
        ax.scatter(flex_dict['subpix'][o]-flex_dict['corr_cen'][o],
                   flex_dict['corr'][o]/mxmod, marker='o')
        # Final shift
        ax.plot([flex_dict['shift'][o]]*2, ylim, 'g:')
        # Label
        if slit_cen:
            ax.text(0.5, 0.25, 'Slit Center', transform=ax.transAxes, size='large', ha='center')
        else:
            ax.text(0.5, 0.25, '{:s}'.format(slf._specobjs[det-1][o].idx), transform=ax.transAxes, size='large', ha='center')
        ax.text(0.5, 0.15, 'flex_shift = {:g}'.format(flex_dict['shift'][o]),
                transform=ax.transAxes, size='large', ha='center')#, bbox={'facecolor':'white'})
        # Axes
        ax.set_ylim(ylim)
        ax.set_xlabel('Lag')

    # Finish
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
    slf._qa.savefig(bbox_inches='tight')
    plt.close()

    # Sky line QA (just one object)
    if slit_cen:
        o=0
    else:
        o=0
        specobj = slf._specobjs[det-1][o]
    sky_spec = flex_dict['sky_spec'][o]
    arx_spec = flex_dict['arx_spec'][o]

    # Sky lines
    sky_lines = np.array([3370.0, 3914.0, 4046.56, 4358.34, 5577.338, 6300.304,
              7340.885, 7993.332, 8430.174, 8919.610, 9439.660,
              10013.99, 10372.88])*u.AA
    dwv = 20.*u.AA
    gdsky = np.where((sky_lines>sky_spec.wvmin) & (sky_lines < sky_spec.wvmax))[0]
    if len(gdsky) == 0:
        msgs.warn("No sky lines for Flexure QA")
        return
    if len(gdsky) > 6:
        idx = np.array([0,1,len(gdsky)//2,len(gdsky)//2+1,-2,-1])
        gdsky = gdsky[idx]

    # Figure
    plt.figure(figsize=(8, 5.0))
    plt.clf()
    nrow, ncol = 2, 3
    gs = gridspec.GridSpec(nrow, ncol)
    if slit_cen:
        plt.suptitle('Sky Comparison for Slit Center',y=1.05)
    else:
        plt.suptitle('Sky Comparison for {:s}'.format(specobj.idx),y=1.05)

    for ii,igdsky in enumerate(gdsky):
        skyline = sky_lines[igdsky]
        ax = plt.subplot(gs[ii//ncol, ii%ncol])
        # Norm
        pix = np.where(np.abs(sky_spec.wavelength-skyline) < dwv)[0]
        f1 = np.sum(sky_spec.flux[pix])
        f2 = np.sum(arx_spec.flux[pix])
        norm = f1/f2
        # Plot
        ax.plot(sky_spec.wavelength[pix], sky_spec.flux[pix], 'k-', label='Obj',
                drawstyle='steps-mid')
        pix2 = np.where(np.abs(arx_spec.wavelength-skyline) < dwv)[0]
        ax.plot(arx_spec.wavelength[pix2], arx_spec.flux[pix2]*norm, 'r-', label='Arx',
                drawstyle='steps-mid')
        # Axes
        ax.xaxis.set_major_locator(plt.MultipleLocator(dwv.value))
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Counts')

    # Legend
    legend = plt.legend(loc='upper left', scatterpoints=1, borderpad=0.3,
                        handletextpad=0.3, fontsize='small', numpoints=1)

    # Finish
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
    slf._qa.savefig(bbox_inches='tight')
    plt.close()

    return


def obj_trace_qa(slf, frame, ltrace, rtrace, root='trace', outfil=None, normalize=True):
    """ Generate a QA plot for the object trace

    Parameters
    ----------
    frame : ndarray
      image
    ltrace : ndarray
      Left edge traces
    rtrace : ndarray
      Right edge traces
    root : str, optional
      Root name for generating output file, e.g. msflat_01blue_000.fits
    outfil : str, optional
      Output file
    normalize : bool, optional
      Normalize the flat?  If not, use zscale for output
    """
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
        for ii in range(ntrc):
            xtrc = (ltrace[:,ii] + rtrace[:,ii])/2.
            ixtrc = np.round(xtrc).astype(int)
            # Simple 'extraction'
            dumi = np.zeros( (frame.shape[0],3) )
            for jj in range(3):
                dumi[:,jj] = frame[ycen,ixtrc-1+jj]
            trc = np.median(dumi, axis=1)
            # Find portion of the image and normalize
            for yy in ycen:
                xi = max(0, int(ltrace[yy,ii])-3)
                xe = min(frame.shape[1],int(rtrace[yy,ii])+3)
                # Fill + normalize
                nrm_frame[yy, xi:xe] = frame[yy,xi:xe] / trc[yy]
        sclmin, sclmax = 0.4, 1.1
    else:
        nrm_frame = frame.copy()
        sclmin, sclmax = zscale(nrm_frame)

    # Plot
    plt.clf()
    fig = plt.figure(dpi=1200)

    plt.rcParams['font.family']= 'times new roman'
    ticks_font = matplotlib.font_manager.FontProperties(family='times new roman', 
       style='normal', size=16, weight='normal', stretch='normal')
    ax = plt.gca()
    for label in ax.get_yticklabels() :
        label.set_fontproperties(ticks_font)
    for label in ax.get_xticklabels() :
        label.set_fontproperties(ticks_font)
    cmm = cm.Greys_r
    mplt = plt.imshow(nrm_frame,origin='lower', cmap=cmm, extent=(0., frame.shape[1], 0., frame.shape[0]))
    mplt.set_clim(vmin=sclmin, vmax=sclmax)

    # Axes
    plt.xlim(0., frame.shape[1])
    plt.ylim(0., frame.shape[0])
    plt.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off')

    # Traces
    for ii in range(ntrc):
        # Left
        plt.plot(ltrace[:,ii]+0.5, ycen, 'r--',alpha=0.7)
        # Right
        plt.plot(rtrace[:,ii]+0.5, ycen, 'g--',alpha=0.7)
        # Label
        iy = int(frame.shape[0]/2.)
        plt.text(ltrace[iy,ii], ycen[iy], '{:d}'.format(ii+1), color='red', ha='center')
        plt.text(rtrace[iy,ii], ycen[iy], '{:d}'.format(ii+1), color='green', ha='center')

    slf._qa.savefig(bbox_inches='tight')
    plt.close()

def obj_profile_qa(slf, specobjs, scitrace):
    """ Generate a QA plot for the object spatial profile
    Parameters
    ----------
    """
    # Setup
    nobj = scitrace['traces'].shape[1]
    ncol = min(3,nobj)
    nrow = nobj // ncol + ((nobj%ncol) > 0)
    # Plot
    plt.figure(figsize=(8, 5.0))
    plt.clf()
    gs = gridspec.GridSpec(nrow, ncol)

    # Plot
    for o in range(nobj):
        fdict = scitrace['opt_profile'][o]
        if 'param' not in fdict.keys():  # Not optimally extracted
            continue
        ax = plt.subplot(gs[o//ncol,o%ncol])

        # Data
        gdp = fdict['mask'] == 0
        ax.scatter(fdict['slit_val'][gdp], fdict['flux_val'][gdp], marker='.',
                   s=0.5, edgecolor='none')

        # Fit
        mn = np.min(fdict['slit_val'][gdp])
        mx = np.max(fdict['slit_val'][gdp])
        xval = np.linspace(mn,mx,1000)
        fit = arutils.func_val(fdict['param'], xval, fdict['func'])
        ax.plot(xval, fit, 'r')
        # Axes
        ax.set_xlim(mn,mx)
        # Label
        ax.text(0.02, 0.90, 'Obj={:s}'.format(specobjs[o].idx),
                transform=ax.transAxes, size='large', ha='left')

    slf._qa.savefig(bbox_inches='tight')
    plt.close()


def slit_trace_qa(slf, frame, ltrace, rtrace, extslit, desc="", root='trace', outfil=None, normalize=True):
    """
    Generate a QA plot for the traces

    Parameters
    ----------
    slf : class
      An instance of the Science Exposure Class
    frame : ndarray
      trace image
    ltrace : ndarray
      Left slit edge traces
    rtrace : ndarray
      Right slit edge traces
    extslit : ndarray
      Mask of extrapolated slits (True = extrapolated)
    desc : str, optional
      A description to be used as a title for the page
    root : str, optional
      Root name for generating output file, e.g. msflat_01blue_000.fits
    outfil : str, optional
      Output file
    normalize: bool, optional
      Normalize the flat?  If not, use zscale for output
    """
    # Outfil
    # if outfil is None:
    #     if '.fits' in root: # Expecting name of msflat FITS file
    #         outfil = root.replace('.fits', '_trc.pdf')
    #         outfil = outfil.replace('MasterFrames', 'Plots')
    #     else:
    #         outfil = root+'.pdf'
    ntrc = ltrace.shape[1]
    ycen = np.arange(frame.shape[0])
    # Normalize flux in the traces
    if normalize:
        nrm_frame = np.zeros_like(frame)
        for ii in range(ntrc):
            xtrc = (ltrace[:, ii] + rtrace[:, ii])/2.
            ixtrc = np.round(xtrc).astype(int)
            # Simple 'extraction'
            dumi = np.zeros((frame.shape[0], 3))
            for jj in range(3):
                dumi[:, jj] = frame[ycen, ixtrc-1+jj]
            trc = np.median(dumi, axis=1)
            # Find portion of the image and normalize
            for yy in ycen:
                xi = max(0, int(ltrace[yy, ii])-3)
                xe = min(frame.shape[1], int(rtrace[yy, ii])+3)
                # Fill + normalize
                nrm_frame[yy, xi:xe] = frame[yy, xi:xe] / trc[yy]
        sclmin, sclmax = 0.4, 1.1
    else:
        nrm_frame = frame.copy()
        sclmin, sclmax = zscale(nrm_frame)

    # Plot
    plt.clf()

    ax = plt.gca()
    set_fonts(ax)
    for label in ax.get_yticklabels():
        label.set_fontproperties(ticks_font)
    for label in ax.get_xticklabels():
        label.set_fontproperties(ticks_font)
    cmm = cm.Greys_r
    mplt = plt.imshow(nrm_frame, origin='lower', cmap=cmm, interpolation=None,
                      extent=(0., frame.shape[1], 0., frame.shape[0]))
    mplt.set_clim(vmin=sclmin, vmax=sclmax)

    # Axes
    plt.xlim(0., frame.shape[1])
    plt.ylim(0., frame.shape[0])
    plt.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off',
                    labelbottom='off', labelleft='off')

    # Traces
    iy = int(frame.shape[0]/2.)
    for ii in range(ntrc):
        if extslit[ii] is True:
            ptyp = ':'
        else:
            ptyp = '--'
        # Left
        plt.plot(ltrace[:, ii]+0.5, ycen, 'r'+ptyp, alpha=0.7)
        # Right
        plt.plot(rtrace[:, ii]+0.5, ycen, 'c'+ptyp, alpha=0.7)
        # Label
        #plt.text(ltrace[iy, ii], ycen[iy], '{0:d}'.format(ii+1), color='red', ha='left')
        plt.text(0.5*(ltrace[iy, ii]+rtrace[iy, ii]), ycen[iy], '{0:d}'.format(ii+1), color='green', ha='center')
    if desc != "":
        plt.suptitle(desc)

    slf._qa.savefig(dpi=1200, orientation='portrait', bbox_inches='tight')
    #pp.savefig()
    #pp.close()
    plt.close()


def set_fonts(ax):
    """ Set axes fonts
    Parameters
    ----------
    plt
    Returns
    -------
    """
    for label in ax.get_yticklabels():
        label.set_fontproperties(ticks_font)
    for label in ax.get_xticklabels():
        label.set_fontproperties(ticks_font)
