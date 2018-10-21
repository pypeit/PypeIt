""" QA for arclines.holy
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np

from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.backends.backend_pdf import PdfPages

from pypeit import utils


def arc_fit_qa(fit, outfile, ids_only=False, title=None):
    """
    QA for Arc spectrum

    Parameters
    ----------
    fit : Wavelength fit
    arc_spec : ndarray
      Arc spectrum
    outfile : str, optional
      Name of output file
    """

    plt.rcdefaults()
    plt.rcParams['font.family']= 'times new roman'

    arc_spec = fit['spec']

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
        plt.savefig(outfile, dpi=800)
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
    wave = utils.func_val(fit['fitc'], xval/fit['xnorm'], 'legendre',
                            minv=fit['fmin'], maxv=fit['fmax'])
    ax_fit.plot(xval, wave, 'r-')
    xmin, xmax = 0., len(arc_spec)
    ax_fit.set_xlim(xmin, xmax)
    ymin,ymax = np.min(wave)*.95,  np.max(wave)*1.05
    ax_fit.set_ylim(np.min(wave)*.95,  np.max(wave)*1.05)
    ax_fit.set_ylabel('Wavelength')
    ax_fit.get_xaxis().set_ticks([]) # Suppress labeling
    # Stats
    wave_fit = utils.func_val(fit['fitc'], fit['xfit'], 'legendre',
                                minv=fit['fmin'], maxv=fit['fmax'])
    rms = np.sqrt(np.sum((fit['yfit']-wave_fit)**2)/len(fit['xfit'])) # Ang
    dwv_pix = np.median(np.abs(wave-np.roll(wave,1)))
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
    plt.savefig(outfile, dpi=800)
    plt.close()

    plt.rcdefaults()

    return


def match_qa(arc_spec, tcent, line_list, IDs, scores, outfile = None, title=None, path=None):
    """
    Parameters
    ----------
    arc_spec
    tcent
    line_list
    IDs
    scores
    outfile
    title
    path

    Returns
    -------

    """


    # Plot
    plt.figure(figsize=(11, 8.5))
    plt.clf()
    gs = gridspec.GridSpec(1, 1)
    idfont = 'small'

    # Simple spectrum plot
    ax_spec = plt.subplot(gs[0])
    ax_spec.plot(np.arange(len(arc_spec)), arc_spec, 'k')
    ymin, ymax = 0., np.max(arc_spec)
    ysep = ymax*0.03
    mn_yline = 1e9

    # Standard IDs
    clrs = dict(Perfect='green', Good='orange', Ok='red')
    clrs['Very Good'] = 'blue'
    for kk, score in enumerate(scores):
        x = tcent[kk]
        # Color
        try:
            clr = clrs[score]
        except KeyError:
            clr = 'gray'
        yline = np.max(arc_spec[int(x)-2:int(x)+2])
        mn_yline = min(mn_yline, yline)
        # Tick mark
        ax_spec.plot([x,x], [yline+ysep*0.25, yline+ysep], '-', color=clr)
        if score in ['Good', 'Ok', 'Perfect', 'Very Good']:
            # Label
            imin = np.argmin(np.abs(line_list['wave']-IDs[kk]))
            row = line_list[imin]
            lbl = '{:s} {:.4f}'.format(row['ion'], row['wave'])
            # label
            ax_spec.text(x, yline+ysep*1.3, '{:s}'.format(lbl), ha='center', va='bottom',
                size=idfont, rotation=90., color=clr)
    # Axes
    ax_spec.set_xlim(0., len(arc_spec))
    ax_spec.set_ylim(ymin, ymax*1.3)
    ax_spec.set_xlabel('Pixel')
    ax_spec.minorticks_on()
    ax_spec.set_ylabel('Counts')
    plt.legend()
    if title is not None:
        ax_spec.text(0.04, 0.93, title, transform=ax_spec.transAxes,
                     size='x-large', ha='left')#, bbox={'facecolor':'white'})
    # Finish
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
    if outfile is None:
        plt.show()
    else:
        pp = PdfPages(outfile)
        pp.savefig(bbox_inches='tight')
        pp.close()

    plt.close()
    return
