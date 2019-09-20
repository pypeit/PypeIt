""" QA for arclines.holy
"""
import numpy as np

from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch

from pypeit import utils


def arc_fit_qa(fit, outfile=None, ids_only=False, title=None):
    """
    QA for Arc spectrum

    Parameters
    ----------
    setup: str
      For outfile
    fit : dict
    arc_spec : ndarray
      Arc spectrum
    outfile : str, optional
      Name of output file
      or 'show' to show on screen
    """

    plt.rcdefaults()
    plt.rcParams['font.family']= 'times new roman'

    arc_spec = fit['spec']

    # Begin
    plt.close('all')
    if ids_only:
        nrows, ncols = 1,1
        figsize =(11,8.5)
        idfont = 'small'
    else:
        nrows, ncols = 2,2
        if outfile is None:
            figsize = (16,8)
            idfont = 'small'
        else:
            figsize = (8,4)
            idfont = 'xx-small'
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(nrows,ncols)#, figure = fig)


    # Simple spectrum plot
    ax_spec = plt.subplot(gs[:,0])
    ax_spec.plot(np.arange(len(arc_spec)), arc_spec)
    ymin, ymax = np.min(arc_spec), np.max(arc_spec)
    ysep = ymax*0.03
    mask = fit['mask']
    pixel_fit = fit['pixel_fit'][mask]
    wave_fit = fit['wave_fit'][mask]
    ions = fit['ions'][mask]
    xnorm = fit['xnorm']
    for kk, x in enumerate(pixel_fit):
        ind_left = np.fmax(int(x)-2, 0)
        ind_righ = np.fmin(int(x)+2,arc_spec.size-1)
        yline = np.max(arc_spec[ind_left:ind_righ])
        # Tick mark
        ax_spec.plot([x,x], [yline+ysep*0.25, yline+ysep], 'g-')
        # label
        ax_spec.text(x, yline+ysep*1.3,'{:s} {:g}'.format(ions[kk], wave_fit[kk]), ha='center', va='bottom',size=idfont,
                     rotation=90., color='green')
    ax_spec.set_xlim(0., len(arc_spec))
    ax_spec.set_ylim(1.05*ymin, ymax*1.2)
    ax_spec.set_xlabel('Pixel')
    ax_spec.set_ylabel('Flux')
    if title is not None:
        ax_spec.text(0.04, 0.93, title, transform=ax_spec.transAxes,
                     size='x-large', ha='left')#, bbox={'facecolor':'white'})
    if ids_only:
        plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
        if outfile is None:
            plt.show()
        else:
            plt.savefig(outfile, dpi=800)
        plt.close()
        return

    # Arc Fit
    ax_fit = plt.subplot(gs[0, 1])
    # Points
    ax_fit.scatter(pixel_fit,wave_fit, marker='x')
    if len(fit['xrej']) > 0:
        ax_fit.scatter(fit['xrej'], fit['yrej'], marker='o',
            edgecolor='gray', facecolor='none')
    # Solution
    xval = np.arange(len(arc_spec))
    wave_soln = fit['wave_soln'] #utils.func_val(fit['fitc'], xval, 'legendre',minx=fit['fmin'], maxx=fit['fmax'])
    ax_fit.plot(xval, wave_soln, 'r-')
    xmin, xmax = 0., len(arc_spec)
    ax_fit.set_xlim(xmin, xmax)
    ymin,ymax = np.min(wave_soln)*.95,  np.max(wave_soln)*1.05
    ax_fit.set_ylim((ymin, ymax))
    ax_fit.set_ylabel('Wavelength')
    ax_fit.get_xaxis().set_ticks([]) # Suppress labeling
    # Stats
    wave_soln_fit = utils.func_val(fit['fitc'], pixel_fit/xnorm, 'legendre',minx=fit['fmin'], maxx=fit['fmax'])
    rms = np.sqrt(np.sum((wave_fit-wave_soln_fit)**2)/len(pixel_fit)) # Ang
    dwv_pix = np.median(np.abs(wave_soln-np.roll(wave_soln,1)))
    ax_fit.text(0.1*len(arc_spec), 0.90*ymin+(ymax-ymin),r'$\Delta\lambda$={:.3f}$\AA$ (per pix)'.format(dwv_pix), size='small')
    ax_fit.text(0.1*len(arc_spec), 0.80*ymin+(ymax-ymin),'RMS={:.3f} (pixels)'.format(rms/dwv_pix), size='small')
    # Arc Residuals
    ax_res = plt.subplot(gs[1,1])
    res = wave_fit-wave_soln_fit
    ax_res.scatter(pixel_fit, res/dwv_pix, marker='x')
    ax_res.plot([xmin,xmax], [0.,0], 'k--')
    ax_res.set_xlim(xmin, xmax)
    ax_res.set_xlabel('Pixel')
    ax_res.set_ylabel('Residuals (Pix)')

    # Finish
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
    if outfile is None:
        plt.show()
    else:
        plt.savefig(outfile, dpi=400)
    plt.close('all')

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
    ymin, ymax = np.min(arc_spec), np.max(arc_spec)
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
    # Overplot the line classification legened
    clrs['Not reidentified'] ='gray'
    legend_elements = []
    for key, clr in clrs.items():
        legend_elements.append(Patch(facecolor=clr, edgecolor=clr,label=key))

    # Axes
    ax_spec.set_xlim(0., len(arc_spec))
    ax_spec.set_ylim(ymin*1.05, ymax*1.3)
    ax_spec.set_xlabel('Pixel')
    ax_spec.minorticks_on()
    ax_spec.set_ylabel('Counts')
    plt.legend(handles=legend_elements)
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
