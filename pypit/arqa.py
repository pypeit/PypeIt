""" Module for QA in PYPIT
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import inspect

import numpy as np
import glob

from astropy import units as u

from pypit.arplot import zscale
from pypit import armsgs
from pypit import arutils

import matplotlib
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages

try:
    basestring
except NameError:  # For Python 3
    basestring = str

msgs = armsgs.get_logger()

# Force the default matplotlib plotting parameters
plt.rcdefaults()
plt.rcParams['font.family']= 'times new roman'
ticks_font = matplotlib.font_manager.FontProperties(family='times new roman',
                                                    style='normal', size=16, weight='normal', stretch='normal')

from pypit import ardebug as debugger


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
    # Outfil
    module = inspect.stack()[0][3]
    outfile = set_qa_filename(slf.setup, module)
    #
    arc_spec = fit['spec']
    if outfil is not None:  # Will deprecate
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
        else:
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
    if slf is not None:
        if False:
            slf._qa.savefig(bbox_inches='tight')
        else:
            plt.savefig(outfile, dpi=800)
            plt.close()
    else:  # FOR ARCLINES MATCH script
        pp.savefig(bbox_inches='tight')
        pp.close()
    return


def coaddspec_qa(ispectra, rspec, spec1d, qafile=None):
    """  QA plot for 1D coadd of spectra

    Parameters
    ----------
    ispectra : XSpectrum1D
      Multi-spectra object
    rspec : XSpectrum1D
      Rebinned spectra with updated variance
    spec1d : XSpectrum1D
      Final coadd

    """
    from pypit.arcoadd import get_std_dev as gsd
    from scipy.stats import norm
    from astropy.stats import sigma_clip

    if qafile is not None:
        pp = PdfPages(qafile)

    plt.clf()
    plt.figure()
    gs = gridspec.GridSpec(1,2)

    # Deviate
    std_dev, dev_sig = gsd(rspec, spec1d)
    #dev_sig = (rspec.data['flux'] - spec1d.flux) / (rspec.data['sig']**2 + spec1d.sig**2)
    #std_dev = np.std(sigma_clip(dev_sig, sigma=5, iters=2))
    flat_dev_sig = dev_sig.flatten()

    xmin = -10
    xmax = 10
    n_bins = 100

    # Deviation
    ax = plt.subplot(gs[0])
    hist, edges = np.histogram(flat_dev_sig, range=(xmin, xmax), bins=n_bins)
    area = len(flat_dev_sig)*((xmax-xmin)/float(n_bins))
    xppf = np.linspace(norm.ppf(0.0001), norm.ppf(0.9999), 100)
    ax.plot(xppf, area*norm.pdf(xppf), color='black', linewidth=2.0)
    ax.bar(edges[:-1], hist, width=((xmax-xmin)/float(n_bins)), alpha=0.5)

    # Coadd on individual
    ax = plt.subplot(gs[1])
    for idx in range(ispectra.nspec):
        ispectra.select = idx
        ax.plot(ispectra.wavelength, ispectra.flux, alpha=0.5)#, label='individual exposure')

    ax.plot(spec1d.wavelength, spec1d.flux, color='black', label='coadded spectrum')
    debug=False
    if debug:
        ax.set_ylim(0., 180.)
        ax.set_xlim(3840, 3860.)
    plt.legend()
    plt.title('Coadded + Original Spectra')

    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.2)
    if qafile is not None:
        pp.savefig(bbox_inches='tight')
        pp.close()
        msgs.info("Wrote coadd QA: {:s}".format(qafile))
    plt.close()
    return


def flexure(slf, det, flex_list, slit_cen=False):
    """ QA on flexure measurement

    Parameters
    ----------
    slf
    det
    flex_list : list
      list of dict containing flexure results
    slit_cen : bool, optional
      QA on slit center instead of objects

    Returns
    -------

    """
    # Outfile
    module = inspect.stack()[0][3]
    #
    for sl in range(len(slf._specobjs[det-1])):
        # Setup
        if slit_cen:
            nobj = 1
            ncol = 1
        else:
            nobj = len(slf._specobjs[det-1][sl])
            if nobj == 0:
                continue
            ncol = min(3, nobj)
        #
        nrow = nobj // ncol + ((nobj % ncol) > 0)

        # Get the flexure dictionary
        flex_dict = flex_list[sl]

        # Outfile
        outfile = set_qa_filename(slf._basename, module+'_corr', det=det,
            slit=slf._specobjs[det-1][sl][0].slitid)

        plt.figure(figsize=(8, 5.0))
        plt.clf()
        gs = gridspec.GridSpec(nrow, ncol)

        # Correlation QA
        for o in range(nobj):
            ax = plt.subplot(gs[o//ncol, o % ncol])
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
                ax.text(0.5, 0.25, '{:s}'.format(slf._specobjs[det-1][sl][o].idx), transform=ax.transAxes, size='large', ha='center')
            ax.text(0.5, 0.15, 'flex_shift = {:g}'.format(flex_dict['shift'][o]),
                    transform=ax.transAxes, size='large', ha='center')#, bbox={'facecolor':'white'})
            # Axes
            ax.set_ylim(ylim)
            ax.set_xlabel('Lag')

        # Finish
        plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
        if False:
            slf._qa.savefig(bbox_inches='tight')
        plt.savefig(outfile, dpi=600)
        plt.close()

        # Sky line QA (just one object)
        if slit_cen:
            o = 0
        else:
            o = 0
            specobj = slf._specobjs[det-1][sl][o]
        sky_spec = flex_dict['sky_spec'][o]
        arx_spec = flex_dict['arx_spec'][o]

        # Sky lines
        sky_lines = np.array([3370.0, 3914.0, 4046.56, 4358.34, 5577.338, 6300.304,
                  7340.885, 7993.332, 8430.174, 8919.610, 9439.660,
                  10013.99, 10372.88])*u.AA
        dwv = 20.*u.AA
        gdsky = np.where((sky_lines > sky_spec.wvmin) & (sky_lines < sky_spec.wvmax))[0]
        if len(gdsky) == 0:
            msgs.warn("No sky lines for Flexure QA")
            return
        if len(gdsky) > 6:
            idx = np.array([0, 1, len(gdsky)//2, len(gdsky)//2+1, -2, -1])
            gdsky = gdsky[idx]

        # Outfile
        outfile = set_qa_filename(slf._basename, module+'_sky', det=det,
                                  slit=slf._specobjs[det-1][sl][0].slitid)
        # Figure
        plt.figure(figsize=(8, 5.0))
        plt.clf()
        nrow, ncol = 2, 3
        gs = gridspec.GridSpec(nrow, ncol)
        if slit_cen:
            plt.suptitle('Sky Comparison for Slit Center', y=1.05)
        else:
            plt.suptitle('Sky Comparison for {:s}'.format(specobj.idx), y=1.05)

        for ii, igdsky in enumerate(gdsky):
            skyline = sky_lines[igdsky]
            ax = plt.subplot(gs[ii//ncol, ii % ncol])
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
        plt.legend(loc='upper left', scatterpoints=1, borderpad=0.3,
                   handletextpad=0.3, fontsize='small', numpoints=1)

        # Finish
        plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
        if False:
            slf._qa.savefig(bbox_inches='tight')
        plt.savefig(outfile, dpi=800)
        plt.close()
        #plt.close()

    return



def obj_trace_qa(slf, frame, ltrace, rtrace, objids, det,
                 root='trace', normalize=True, desc=""):
    """ Generate a QA plot for the object trace

    Parameters
    ----------
    frame : ndarray
      image
    ltrace : ndarray
      Left edge traces
    rtrace : ndarray
      Right edge traces
    objids : list
    det : int
    desc : str, optional
      Title
    root : str, optional
      Root name for generating output file, e.g. msflat_01blue_000.fits
    normalize : bool, optional
      Normalize the flat?  If not, use zscale for output
    """
    module = inspect.stack()[0][3]
    outfile = set_qa_filename(slf._basename, module, det=det)
    #
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

    plt.rcParams['font.family'] = 'times new roman'
    ticks_font = matplotlib.font_manager.FontProperties(family='times new roman', 
       style='normal', size=16, weight='normal', stretch='normal')
    ax = plt.gca()
    for label in ax.get_yticklabels() :
        label.set_fontproperties(ticks_font)
    for label in ax.get_xticklabels() :
        label.set_fontproperties(ticks_font)
    cmm = cm.Greys_r
    mplt = plt.imshow(nrm_frame, origin='lower', cmap=cmm, extent=(0., frame.shape[1], 0., frame.shape[0]))
    mplt.set_clim(vmin=sclmin, vmax=sclmax)

    # Axes
    plt.xlim(0., frame.shape[1])
    plt.ylim(0., frame.shape[0])
    plt.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off')

    # Traces
    iy_mid = int(frame.shape[0] / 2.)
    iy = np.linspace(iy_mid,frame.shape[0]*0.9,ntrc).astype(int)
    for ii in range(ntrc):
        # Left
        plt.plot(ltrace[:, ii]+0.5, ycen, 'r--', alpha=0.7)
        # Right
        plt.plot(rtrace[:, ii]+0.5, ycen, 'c--', alpha=0.7)
        if objids is not None:
            # Label
            # plt.text(ltrace[iy,ii], ycen[iy], '{:d}'.format(ii+1), color='red', ha='center')
            lbl = 'O{:03d}'.format(objids[ii])
            plt.text((ltrace[iy[ii], ii]+rtrace[iy[ii], ii])/2., ycen[iy[ii]],
                lbl, color='green', ha='center')
    # Title
    tstamp = gen_timestamp()
    if desc == "":
        plt.suptitle(tstamp)
    else:
        plt.suptitle(desc+'\n'+tstamp)

    if False:
        slf._qa.savefig(dpi=1200, orientation='portrait', bbox_inches='tight')
    plt.savefig(outfile, dpi=800)
    plt.close()


def obj_profile_qa(slf, specobjs, scitrace, det):
    """ Generate a QA plot for the object spatial profile
    Parameters
    ----------
    """
    module = inspect.stack()[0][3]
    for sl in range(len(specobjs)):
        # Setup
        nobj = scitrace[sl]['traces'].shape[1]
        ncol = min(3, nobj)
        nrow = nobj // ncol + ((nobj % ncol) > 0)
        # Outfile
        outfile = set_qa_filename(slf._basename, module, det=det, slit=specobjs[sl][0].slitid)
        # Plot
        plt.figure(figsize=(8, 5.0))
        plt.clf()
        gs = gridspec.GridSpec(nrow, ncol)

        # Plot
        for o in range(nobj):
            fdict = scitrace[sl]['opt_profile'][o]
            if 'param' not in fdict.keys():  # Not optimally extracted
                continue
            ax = plt.subplot(gs[o//ncol, o % ncol])

            # Data
            gdp = fdict['mask'] == 0
            ax.scatter(fdict['slit_val'][gdp], fdict['flux_val'][gdp], marker='.',
                       s=0.5, edgecolor='none')

            # Fit
            mn = np.min(fdict['slit_val'][gdp])
            mx = np.max(fdict['slit_val'][gdp])
            xval = np.linspace(mn, mx, 1000)
            fit = arutils.func_val(fdict['param'], xval, fdict['func'])
            ax.plot(xval, fit, 'r')
            # Axes
            ax.set_xlim(mn,mx)
            # Label
            ax.text(0.02, 0.90, 'Obj={:s}'.format(specobjs[sl][o].idx),
                    transform=ax.transAxes, ha='left', size='small')

        if False:
            slf._qa.savefig(bbox_inches='tight')
        plt.savefig(outfile, dpi=500)
        plt.close()


def plot_orderfits(slf, model, ydata, xdata=None, xmodl=None, textplt="Slit",
                   maxp=4, desc="", maskval=-999999.9):
    """ Generate a QA plot for the blaze function fit to each slit
    Or the arc line tilts

    Parameters
    ----------
    slf : class
      Science Exposure class
    model : ndarray
      (m x n) 2D array containing the model blaze function (m) of a flat frame for each slit (n)
    ydata : ndarray
      (m x n) 2D array containing the extracted 1D spectrum (m) of a flat frame for each slit (n)
    xdata : ndarray, optional
      x values of the data points
    xmodl : ndarry, optional
      x values of the model points
    textplt : str, optional
      A string printed above each panel
    maxp : int, (optional)
      Maximum number of panels per page
    desc : str, (optional)
      A description added to the top of each page
    maskval : float, (optional)
      Value used in arrays to indicate a masked value
    """
    # Outfil
    module = inspect.stack()[0][3]
    if 'Arc' in desc:
        module += '_Arc'
    elif 'Blaze' in desc:
        module += '_Blaze'
    else:
        msgs.error("Should not get here")
    outroot = set_qa_filename(slf.setup, module)
    #
    npix, nord = ydata.shape
    pages, npp = get_dimen(nord, maxp=maxp)
    if xdata is None: xdata = np.arange(npix).reshape((npix, 1)).repeat(nord, axis=1)
    if xmodl is None: xmodl = np.arange(model.shape[0])
    # Loop through all pages and plot the results
    ndone = 0
    axesIdx = True
    for i in range(len(pages)):
        f, axes = plt.subplots(pages[i][1], pages[i][0])
        ipx, ipy = 0, 0
        for j in range(npp[i]):
            if pages[i][0] == 1 and pages[i][1] == 1: axesIdx = False
            elif pages[i][1] == 1: ind = (ipx,)
            elif pages[i][0] == 1: ind = (ipy,)
            else: ind = (ipy, ipx)
            if axesIdx:
                axes[ind].plot(xdata[:,ndone+j], ydata[:,ndone+j], 'bx', drawstyle='steps')
                axes[ind].plot(xmodl, model[:,ndone+j], 'r-')
            else:
                axes.plot(xdata[:,ndone+j], ydata[:,ndone+j], 'bx', drawstyle='steps')
                axes.plot(xmodl, model[:,ndone+j], 'r-')
            ytmp = ydata[:,ndone+j]
            gdy = ytmp != maskval
            ytmp = ytmp[gdy]
            if ytmp.size != 0:
                amn = min(np.min(ytmp), np.min(model[gdy,ndone+j]))
            else:
                amn = np.min(model[:,ndone+j])
            if ytmp.size != 0:
                amx = max(np.max(ytmp), np.max(model[gdy,ndone+j]))
            else: amx = np.max(model[:,ndone+j])
            # Restrict to good pixels
            xtmp = xdata[:,ndone+j]
            gdx = xtmp != maskval
            xtmp = xtmp[gdx]
            if xtmp.size == 0:
                xmn = np.min(xmodl)
                xmx = np.max(xmodl)
            else:
                xmn = np.min(xtmp)
                xmx = np.max(xtmp)
                #xmn = min(np.min(xtmp), np.min(xmodl))
                #xmx = max(np.max(xtmp), np.max(xmodl))
            if axesIdx:
                axes[ind].axis([xmn, xmx, amn-1, amx+1])
                axes[ind].set_title("{0:s} {1:d}".format(textplt, ndone+j+1))
            else:
                axes.axis([xmn, xmx, amn, amx])
                axes.set_title("{0:s} {1:d}".format(textplt, ndone+j+1))
            ipx += 1
            if ipx == pages[i][0]:
                ipx = 0
                ipy += 1
        # Delete the unnecessary axes
        if axesIdx:
            for j in range(npp[i], axes.size):
                if pages[i][1] == 1: ind = (ipx,)
                elif pages[i][0] == 1: ind = (ipy,)
                else: ind = (ipy, ipx)
                f.delaxes(axes[ind])
                if ipx == pages[i][0]:
                    ipx = 0
                    ipy += 1
        ndone += npp[i]
        # Save the figure
        if axesIdx: axsz = axes.size
        else: axsz = 1.0
        if pages[i][1] == 1 or pages[i][0] == 1: ypngsiz = 11.0/axsz
        else: ypngsiz = 11.0*axes.shape[0]/axes.shape[1]
        f.set_size_inches(11.0, ypngsiz)
        if desc != "":
            pgtxt = ""
            if len(pages) != 1:
                pgtxt = ", page {0:d}/{1:d}".format(i+1, len(pages))
            f.suptitle(desc + pgtxt, y=1.02, size=16)
        f.tight_layout()
        if False:
            slf._qa.savefig(dpi=200, orientation='landscape', bbox_inches='tight')
        else:
            outfile = outroot+'{:03d}.png'.format(i)
            plt.savefig(outfile, dpi=200)
            plt.close()
        f.clf()
    del f
    return


def slit_profile(slf, mstrace, model, lordloc, rordloc, msordloc, textplt="Slit", maxp=16, desc=""):
    """ Generate a QA plot for the slit profile of each slit

    Parameters
    ----------
    slf : class
      Science Exposure class
    mstrace : ndarray
      trace frame
    model : ndarray
      model of slit profiles, same shape as frame.
    lordloc : ndarray
      left edge locations of all slits
    rordloc : ndarray
      right edge locations of all slits
    msordloc : ndarray
      An array the same size as frame that determines which pixels contain a given order.
    textplt : str, optional
      A string printed above each panel
    maxp : int, (optional)
      Maximum number of panels per page
    desc : str, (optional)
      A description added to the top of each page
    """

    npix, nord = lordloc.shape
    nbins = 40
    pages, npp = get_dimen(nord, maxp=maxp)
    # Loop through all pages and plot the results
    ndone = 0
    axesIdx = True
    for i in range(len(pages)):
        f, axes = plt.subplots(pages[i][1], pages[i][0])
        ipx, ipy = 0, 0
        for j in range(npp[i]):
            if pages[i][0] == 1 and pages[i][1] == 1: axesIdx = False
            elif pages[i][1] == 1: ind = (ipx,)
            elif pages[i][0] == 1: ind = (ipy,)
            else: ind = (ipy, ipx)
            # Get data to be plotted
            word = np.where(msordloc == ndone+j+1)
            if word[0].size == 0:
                msgs.warn("There are no pixels in slit {0:d}".format(ndone + j + 1))
                # Delete the axis
                if pages[i][1] == 1: ind = (ipx,)
                elif pages[i][0] == 1: ind = (ipy,)
                else: ind = (ipy, ipx)
                f.delaxes(axes[ind])
                ipx += 1
                if ipx == pages[i][0]:
                    ipx = 0
                    ipy += 1
                continue
            spatval = (word[1] + 0.5 - lordloc[:, ndone+j][word[0]]) / (rordloc[:, ndone+j][word[0]] - lordloc[:, ndone+j][word[0]])
            fluxval = mstrace[word]
            mxval = np.max(fluxval)
            modvals = np.zeros(nbins)
            if axesIdx:
                cnts, xedges, yedges, null = axes[ind].hist2d(spatval, fluxval, bins=nbins, cmap=plt.cm.Greys)
                groups = np.digitize(spatval, xedges)
                modelw = model[word]
                for mm in range(1, xedges.size):
                    modvals[mm-1] = modelw[groups == mm].mean()
                axes[ind].plot(0.5*(xedges[1:]+xedges[:-1]), modvals, 'g-', linewidth=2.0)
                axes[ind].plot([0.0, 0.0], [0.0, mxval], 'r-')
                axes[ind].plot([1.0, 1.0], [0.0, mxval], 'r-')
            else:
                cnts, xedges, yedges, null = axes.hist2d(spatval, fluxval, bins=nbins, cmap=plt.cm.Greys)
                groups = np.digitize(spatval, xedges)
                modelw = model[word]
                for mm in range(1, xedges.size):
                    modvals[mm-1] = modelw[groups == mm].mean()
                axes.plot(0.5*(xedges[1:]+xedges[:-1]), modvals, 'g-', linewidth=2.0)
                axes.plot([0.0, 0.0], [0.0, mxval], 'r-')
                axes.plot([1.0, 1.0], [0.0, mxval], 'r-')
            if axesIdx:
                axes[ind].axis([xedges[0], xedges[-1], 0.0, 1.1*mxval])
                axes[ind].set_title("{0:s} {1:d}".format(textplt, ndone+j+1))
                axes[ind].tick_params(labelsize=10)
            else:
                axes.axis([xedges[0], xedges[-1], 0.0, 1.1*mxval])
                axes.set_title("{0:s} {1:d}".format(textplt, ndone+j+1))
                axes.tick_params(labelsize=10)
            ipx += 1
            if ipx == pages[i][0]:
                ipx = 0
                ipy += 1
        # Delete the unnecessary axes
        if axesIdx:
            for j in range(npp[i], axes.size):
                if pages[i][1] == 1: ind = (ipx,)
                elif pages[i][0] == 1: ind = (ipy,)
                else: ind = (ipy, ipx)
                f.delaxes(axes[ind])
                ipx += 1
                if ipx == pages[i][0]:
                    ipx = 0
                    ipy += 1
        ndone += npp[i]
        # Save the figure
        if axesIdx: axsz = axes.size
        else: axsz = 1.0
        if pages[i][1] == 1 or pages[i][0] == 1: ypngsiz = 11.0/axsz
        else: ypngsiz = 11.0*axes.shape[0]/axes.shape[1]
        f.set_size_inches(11.0, ypngsiz)
        if desc != "":
            pgtxt = ""
            if len(pages) != 1:
                pgtxt = ", page {0:d}/{1:d}".format(i+1, len(pages))
            f.suptitle(desc + pgtxt, y=1.02, size=16)
        f.tight_layout()
        slf._qa.savefig(dpi=200, orientation='landscape', bbox_inches='tight')
        #plt.close()
        f.clf()
    del f
    return


def slit_trace_qa(slf, frame, ltrace, rtrace, extslit, desc="",
                  root='trace', outfil=None, normalize=True, use_slitid=None):
    """ Generate a QA plot for the slit traces

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
    module = inspect.stack()[0][3]
    outfile = set_qa_filename(slf.setup, module)
    # if outfil is None:
    #     if '.fits' in root: # Expecting name of msflat FITS file
    #         outfil = root.replace('.fits', '_trc.pdf')
    #         outfil = outfil.replace('MasterFrames', 'Plots')
    #     else:
    #         outfil = root+'.pdf'
    from pypit.arspecobj import get_slitid
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
        nrm_frame[frame > 0.0] = np.sqrt(nrm_frame[frame > 0.0])
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
        plt.plot(ltrace[:, ii]+0.5, ycen, 'r'+ptyp, linewidth=0.3, alpha=0.7)
        # Right
        plt.plot(rtrace[:, ii]+0.5, ycen, 'c'+ptyp, linewidth=0.3, alpha=0.7)
        # Label
        if use_slitid:
            slitid, _, _ = get_slitid(slf, use_slitid, ii, ypos=0.5)
            lbl = 'S{:04d}'.format(slitid)
        else:
            lbl = '{0:d}'.format(ii+1)
        plt.text(0.5*(ltrace[iy, ii]+rtrace[iy, ii]), ycen[iy], lbl, color='green', ha='center', size='small')
    # Title
    tstamp = gen_timestamp()
    if desc == "":
        plt.suptitle(tstamp)
    else:
        plt.suptitle(desc+'\n'+tstamp)

    # Write
    if False:
        slf._qa.savefig(dpi=1200, orientation='portrait', bbox_inches='tight')
    plt.savefig(outfile, dpi=800)
    plt.close()


def pca_plot(slf, inpar, ofit, prefix, maxp=25, pcadesc="", addOne=True):
    """ Saves quality control plots for a PCA analysis
    Parameters
    ----------
    slf
    inpar
    ofit
    prefix : str
      prefix for the filenames
    maxp
    pcadesc
    addOne

    Returns
    -------

    """
    # Setup
    module = inspect.stack()[0][3]
    outroot = set_qa_filename(slf.setup, module, prefix=prefix)
    npc = inpar['npc']+1
    pages, npp = get_dimen(npc, maxp=maxp)
    #
    x0 = inpar['x0']
    ordernum = inpar['x0in']
    x0fit = inpar['x0fit']
    usetrc = inpar['usetrc']
    hidden = inpar['hidden']
    high_fit = inpar['high_fit']
    nc = np.max(ordernum[usetrc])
    # Loop through all pages and plot the results
    ndone = 0
    for i in range(len(pages)):
        plt.clf()
        f, axes = plt.subplots(pages[i][1], pages[i][0])
        ipx, ipy = 0, 0
        if i == 0:
            if pages[i][1] == 1: ind = (0,)
            elif pages[i][0] == 1: ind = (0,)
            else: ind = (0,0)
            axes[ind].plot(ordernum[usetrc], x0[usetrc], 'bx')
            axes[ind].plot(ordernum, x0fit, 'k-')
            amn, amx = np.min(x0fit), np.max(x0fit)
            diff = x0[usetrc]-x0fit[usetrc]
            tdiffv = np.median(diff)
            mdiffv = 1.4826*np.median(np.abs(tdiffv-diff))
            amn -= 2.0*mdiffv
            amx += 2.0*mdiffv
            mval = amn-0.15*(amx-amn)
            dmin, dmax = tdiffv-2.0*mdiffv, tdiffv+2.0*mdiffv
            diff = mval + diff*0.20*(amx-amn)/(dmax-dmin)
            wign = np.where(np.abs(diff-np.median(diff))<4.0*1.4826*np.median(np.abs(diff-np.median(diff))))[0]
            dmin, dmax = np.min(diff[wign]), np.max(diff[wign])
            axes[ind].plot(ordernum[usetrc], diff, 'rx')
            if addOne:
                axes[ind].plot([0, nc+1], [mval,mval], 'k-')
                axes[ind].axis([0, nc+1, dmin-0.5*(dmax-dmin), amx + 0.05*(amx-amn)])
            else:
                axes[ind].plot([0, nc], [mval, mval], 'k-')
                axes[ind].axis([0, nc, dmin-0.5*(dmax-dmin), amx + 0.05*(amx-amn)])
            axes[ind].set_title("Mean Value")
            ipx += 1
            if ipx == pages[i][0]:
                ipx = 0
                ipy += 1
            npp[0] -= 1
        for j in range(npp[i]):
            if pages[i][1] == 1: ind = (ipx,)
            elif pages[i][0] == 1: ind = (ipy,)
            else: ind = (ipy, ipx)
            axes[ind].plot(ordernum[usetrc], hidden[j+ndone,:], 'bx')
            axes[ind].plot(ordernum, high_fit[:,j+ndone], 'k-')
            vmin, vmax = np.min(hidden[j+ndone,:]), np.max(hidden[j+ndone,:])
            if ofit[1+j+ndone] != -1:
                cmn, cmx = np.min(high_fit[:,j+ndone]), np.max(high_fit[:,j+ndone])
                diff = hidden[j+ndone,:]-high_fit[:,j+ndone][usetrc]
                tdiffv = np.median(diff)
                mdiffv = 1.4826*np.median(np.abs(tdiffv-diff))
                cmn -= 2.0*mdiffv
                cmx += 2.0*mdiffv
                mval = cmn-0.15*(cmx-cmn)
                dmin, dmax = tdiffv-2.0*mdiffv, tdiffv+2.0*mdiffv
                #dmin, dmax = np.min(diff), np.max(diff)
                diff = mval + diff*0.20*(cmx-cmn)/(dmax-dmin)
                wign = np.where(np.abs(diff-np.median(diff))<4.0*1.4826*np.median(np.abs(diff-np.median(diff))))[0]
                dmin, dmax = np.min(diff[wign]), np.max(diff[wign])
                #vmin, vmax = np.min(hidden[j+ndone,:][wign]), np.max(hidden[j+ndone,:][wign])
                axes[ind].plot(ordernum[usetrc], diff, 'rx')
                axes[ind].plot([0, 1+nc], [mval, mval], 'k-')
#				ymin = np.min([(3.0*dmin-dmax)/2.0,vmin-0.1*(vmax-dmin),dmin-0.1*(vmax-dmin)])
#				ymax = np.max([np.max(high_fit[:,j+ndone]),vmax+0.1*(vmax-dmin),dmax+0.1*(vmax-dmin)])
                ymin = dmin-0.5*(dmax-dmin)
                ymax = cmx + 0.05*(cmx-cmn)
                if addOne: axes[ind].axis([0, nc+1, ymin, ymax])
                else: axes[ind].axis([0, nc, ymin, ymax])
            else:
                if addOne: axes[ind].axis([0, nc+1, vmin-0.1*(vmax-vmin), vmax+0.1*(vmax-vmin)])
                else: axes[ind].axis([0, nc, vmin-0.1*(vmax-vmin), vmax+0.1*(vmax-vmin)])
            axes[ind].set_title("PC {0:d}".format(j+ndone))
            axes[ind].tick_params(labelsize=8)
            ipx += 1
            if ipx == pages[i][0]:
                ipx = 0
                ipy += 1
        if i == 0: npp[0] = npp[0] + 1
        # Delete the unnecessary axes
        for j in range(npp[i], axes.size):
            if pages[i][1] == 1: ind = (ipx,)
            elif pages[i][0] == 1: ind = (ipy,)
            else: ind = (ipy, ipx)
            f.delaxes(axes[ind])
            ipx += 1
            if ipx == pages[i][0]:
                ipx = 0
                ipy += 1
        ndone += npp[i]
        # Save the figure
        if pages[i][1] == 1 or pages[i][0] == 1: ypngsiz = 11.0/axes.size
        else: ypngsiz = 11.0*axes.shape[0]/axes.shape[1]
        f.set_size_inches(11.0, ypngsiz)
        if pcadesc != "":
            pgtxt = ""
            if len(pages) != 1:
                pgtxt = ", page {0:d}/{1:d}".format(i+1, len(pages))
            f.suptitle(pcadesc + pgtxt, y=1.02, size=16)
        f.tight_layout()
        if False:
            slf._qa.savefig(dpi=200, orientation='landscape', bbox_inches='tight')
        outfile = outroot+'{:02d}.png'.format(i)
        f.savefig(outfile, dpi=200)
        plt.close()
        f.clf()
    del f
    return


def pca_arctilt(slf, tiltang, centval, tilts, maxp=25, maskval=-999999.9):
    """ Generate a QA plot for the blaze function fit to each slit

    Parameters
    ----------
    slf : class
      Science Exposure class
    tiltang : ndarray
      (m x n) 2D array containing the measured tilts (m) for each slit (n)
    centval : ndarray
      (m x n) 2D array containing the pixel location (in the spectral direction) of the measured tilts (m) for each slit (n)
    tilts : ndarray
      (m x n) 2D array containing the model tilts (m) for each slit (n)
    maxp : int, (optional)
      Maximum number of panels per page
    maskval : float, (optional)
      Value used in arrays to indicate a masked value
    """
    module = inspect.stack()[0][3]
    outroot = set_qa_filename(slf.setup, module)
    #
    npc = tiltang.shape[1]
    pages, npp = get_dimen(npc, maxp=maxp)
    x0 = np.arange(tilts.shape[0])
    # First calculate the min and max values for the plotting axes, to make sure they are all the same
    w = np.where(tiltang != maskval)
    medv = np.median(tiltang[w])
    madv = 1.4826 * np.median(np.abs(medv - tiltang[w]))
    ymin, ymax = medv - 3.0 * madv, medv + 3.0 * madv
    ymin = min(ymin, np.min(tilts))
    ymax = max(ymax, np.max(tilts))
    # Check that ymin and ymax are set, if not, return without plotting
    if ymin is None or ymax is None:
        msgs.warn("Arc tilt fits were not plotted")
        return
    # Generate the plots
    ndone = 0
    for i in range(len(pages)):
        f, axes = plt.subplots(pages[i][1], pages[i][0])
        ipx, ipy = 0, 0
        for j in range(npp[i]):
            if pages[i][1] == 1:
                ind = (ipx,)
            elif pages[i][0] == 1:
                ind = (ipy,)
            else:
                ind = (ipy, ipx)
            w = np.where(tiltang[:, ndone] != maskval)
            if np.size(w[0]) != 0:
                axes[ind].plot(centval[:, ndone][w], tiltang[:, ndone][w], 'bx')
                axes[ind].plot(x0, tilts[:, ndone], 'r-')
            axes[ind].axis([0, tilts.shape[0]-1, ymin, ymax])
            axes[ind].set_title("Slit {0:d}".format(1+ndone))
            axes[ind].tick_params(labelsize=8)
            ipx += 1
            if ipx == pages[i][0]:
                ipx = 0
                ipy += 1
            ndone += 1
        # Delete the unnecessary axes
        for j in range(npp[i], axes.size):
            if pages[i][1] == 1:
                ind = (ipx,)
            elif pages[i][0] == 1:
                ind = (ipy,)
            else:
                ind = (ipy, ipx)
            f.delaxes(axes[ind])
            ipx += 1
            if ipx == pages[i][0]:
                ipx = 0
                ipy += 1
        # Save the figure
        if pages[i][1] == 1 or pages[i][0] == 1:
            ypngsiz = 11.0/axes.size
        else:
            ypngsiz = 11.0*axes.shape[0]/axes.shape[1]
        f.set_size_inches(11.0, ypngsiz)
        f.tight_layout()
        if False:
            slf._qa.savefig(dpi=200, orientation='landscape', bbox_inches='tight')
        outfile = outroot+'{:02d}.png'.format(i)
        f.savefig(outfile, dpi=200)
        plt.close()
        f.clf()
        del f
    return



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


def set_qa_filename(root, module, det=None, slit=None, prefix=None):
    """
    Parameters
    ----------
    root : str
      Root name
    module : str
      Describes the QA routine
    det : int, optional
    slit : int, optional
    prefix : str, optional
      start the name of the QA file

    Returns
    -------
    outfile : str
      Filename
    """
    if module == 'slit_trace_qa':
        outfile = 'QA/PNGs/Slit_Trace_{:s}.png'.format(root)
    elif module == 'arc_fit_qa':
        outfile = 'QA/PNGs/Arc_1dfit_{:s}.png'.format(root)
    elif module == 'plot_orderfits_Arc':  # This is root for multiple PNGs
        outfile = 'QA/PNGs/Arc_tilts_{:s}_'.format(root)
    elif module == 'pca_plot':  # This is root for multiple PNGs
        outfile = 'QA/PNGs/{:s}_pca_{:s}_'.format(prefix, root)
    elif module == 'pca_arctilt':  # This is root for multiple PNGs
        outfile = 'QA/PNGs/Arc_pca_{:s}_'.format(root)
    elif module == 'plot_orderfits_Blaze':  # This is root for multiple PNGs
        outfile = 'QA/PNGs/Blaze_{:s}_'.format(root)
    elif module == 'obj_trace_qa':
        outfile = 'QA/PNGs/{:s}_D{:02d}_obj_trace.png'.format(root, det)
    elif module == 'obj_profile_qa':
        outfile = 'QA/PNGs/{:s}_D{:02d}_S{:04d}_obj_prof.png'.format(root, det, slit)
    elif module == 'flexure_corr':
        outfile = 'QA/PNGs/{:s}_D{:02d}_S{:04d}_flex_corr.png'.format(root, det, slit)
    elif module == 'flexure_sky':
        outfile = 'QA/PNGs/{:s}_D{:02d}_S{:04d}_flex_sky.png'.format(root, det, slit)
    else:
        msgs.error("NOT READY FOR THIS QA: {:s}".format(module))
    # Return
    return outfile


def get_dimen(x, maxp=25):
    """ Assign the plotting dimensions to be the "most square"

    Parameters
    ----------
    x : int
      An integer that equals the number of panels to be plot
    maxp : int (optional)
      The maximum number of panels to plot on a single page

    Returns
    -------
    pages : list
      The number of panels in the x and y direction on each page
    npp : list
      The number of panels on each page
    """
    pages, npp = [], []
    xr = x
    while xr > 0:
        if xr > maxp:
            xt = maxp
        else:
            xt = xr
        ypg = int(np.sqrt(np.float(xt)))
        if int(xt) % ypg == 0:
            xpg = int(xt)/ypg
        else:
            xpg = 1 + int(xt)/ypg
        pages.append([int(xpg), int(ypg)])
        npp.append(int(xt))
        xr -= xt
    return pages, npp


def gen_timestamp():
    """ Generate a simple time stamp including the current user
    Returns
    -------
    timestamp : str
      user_datetime
    """
    import datetime
    tstamp = datetime.datetime.today().strftime('%Y-%b-%d-T%Hh%Mm%Ss')
    import getpass
    user = getpass.getuser()
    # Return
    return '{:s}_{:s}'.format(user, tstamp)


def html_header(title):
    """ 
    Parameters
    ----------
    title : str, optional

    Returns
    -------

    """
    head = '<?xml version="1.0" encoding="UTF-8"?>\n'
    head += '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n'


    head += '<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">\n'
    head += '\n'
    head += '<head>\n'
    head += '\n'
    head += '<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />\n'
    head += '<title>{:s}</title>\n'.format(title)
    head += '<meta name="keywords" content="" />\n'
    head += '<meta name="description" content="" />\n'
    head += '<script type="text/javascript" src="jquery/jquery-1.4.2.min.js"></script>\n'
    head += '<script type="text/javascript" src="jquery/jquery.slidertron-0.1.js"></script>\n'
    head += '<link href="style.css" rel="stylesheet" type="text/css" media="screen" />\n'
    head += '\n'
    head += '</head>\n'

    # Begin the Body
    head += '<body>\n'
    head += '<h1>{:s}</h1>\n'.format(title)
    head += '<hr>\n'

    return head

def html_end(f, body, links=None):
    """ Fill in the HTML file and end it
    Parameters
    ----------
    f : file
    body : str
    links : str, optional
    """
    # Write links
    if links is not None:
        f.write(links)
        f.write('</ul>\n')
        f.write('<hr>\n')
    # Write body
    f.write(body)
    # Finish
    end = '</body>\n'
    end += '</html>\n'
    f.write(end)

    return end


def html_init(f, title):
    head = html_header(title)
    f.write(head)
    # Init links
    links = '<h2>Quick Links</h2>\n'
    links += '<ul>\n'
    return links


def html_mf_pngs(setup, cbset, det):
    """ Geneate HTML for MasterFrame PNGs
    Parameters
    ----------
    setup : str
    cbset : str
    det : int

    Returns
    -------
    links : str
    body : str

    """
    links = ''
    body = ''
    # QA root
    # Search for PNGs
    idval = '{:s}_{:02d}_{:s}'.format(setup, det, cbset)

    # Organize the outputs
    html_dict = {}
    html_dict['strace'] = dict(fname='slit_trace_qa', ext='',
        href='strace', label='Slit Trace')
    html_dict['blaze'] = dict(fname='plot_orderfits_Blaze', ext='*.png',
                               href='blaze', label='Blaze')
    html_dict['arc_fit'] = dict(fname='arc_fit_qa', ext='',
                              href='arc_fit', label='Arc 1D Fit')
    html_dict['arc_tilt'] = dict(fname='plot_orderfits_Arc', ext='*.png',
                              href='arc_tilts', label='Arc Tilts')
    html_dict['arc_pca'] = dict(fname='pca_arctilt', ext='*.png',
                                 href='arc_pca', label='Arc Tilt PCA')

    # Generate HTML
    for key in ['strace', 'blaze', 'arc_fit', 'arc_pca', 'arc_tilt']:
        png_root = set_qa_filename(idval, html_dict[key]['fname'])
        pngs = glob.glob(png_root+html_dict[key]['ext'])
        if len(pngs) > 0:
            href="{:s}_{:s}".format(html_dict[key]['href'], idval)
            # Link
            links += '<li><a class="reference internal" href="#{:s}">{:s} {:s}</a></li>\n'.format(href, html_dict[key]['label'], idval)
            # Body
            body += '<hr>\n'
            body += '<div class="section" id="{:s}">\n'.format(href)
            body += '<h2> {:s} {:s} </h2>\n'.format(html_dict[key]['label'], idval)
            for png in pngs:
                # Remove QA
                ifnd = png.find('QA/')
                if ifnd < 0:
                    msgs.error("QA is expected to be in the path!")
                body += '<img class ="research" src="{:s}" width="100%" height="auto"/>\n'.format(png[ifnd+3:])
            body += '</div>\n'

    # Return
    return links, body


def html_exp_pngs(exp_name, det):
    """ 
    Parameters
    ----------
    exp_name : str
    det : int

    Returns
    -------
    links : str
    body : str

    """
    import os
    links = ''
    body = ''
    # QA root

    # Organize the outputs
    html_dict = {}
    html_dict['trace'] = dict(fname='obj_trace_qa', ext='', slit=False,
                               href='otrace', label='Object Traces')
    html_dict['prof'] = dict(fname='obj_profile_qa', ext='', slit=True,
                              href='oprofile', label='Object Profiles')
    html_dict['flex_corr'] = dict(fname='flexure_corr', ext='', slit=True,
                             href='flex_corr', label='Flexure Cross Correlation')
    html_dict['flex_sky'] = dict(fname='flexure_sky', ext='', slit=True,
                                  href='flex_sky', label='Flexure Sky')

    # Generate HTML
    for key in ['trace', 'prof', 'flex_corr', 'flex_sky']:
        png_root = set_qa_filename(exp_name, html_dict[key]['fname'], det=det, slit=9999)
        if html_dict[key]['slit']:  # Kludge to handle multiple slits
            png_root = png_root.replace('S9999', 'S*')
        pngs = glob.glob(png_root+html_dict[key]['ext'])
        if len(pngs) > 0:
            href="{:s}_{:02d}".format(html_dict[key]['href'], det)
            # Link
            links += '<li><a class="reference internal" href="#{:s}">{:s} {:02d}</a></li>\n'.format(href, html_dict[key]['label'], det)
            # Body
            body += '<hr>\n'
            body += '<div class="section" id="{:s}">\n'.format(href)
            body += '<h2> {:s} {:02d} </h2>\n'.format(html_dict[key]['label'], det)
            for png in pngs:
                # Remove QA
                ifnd = png.find('QA/')
                if ifnd < 0:
                    msgs.error("QA is expected to be in the path!")
                body += '<img class ="research" src="{:s}" width="100%" height="auto"/>\n'.format(png[ifnd+3:])
            body += '</div>\n'

    # Return
    return links, body


