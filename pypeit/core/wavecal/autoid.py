""" Module for finding patterns in arc line spectra

.. include:: ../include/links.rst
"""
import copy
import itertools

import astropy.stats
import astropy.table
import numpy as np
import scipy.interpolate
import scipy.ndimage
import scipy.spatial

from linetools import utils as ltu

from IPython import embed


from pypeit.par import pypeitpar
from pypeit.core.wavecal import kdtree_generator
from pypeit.core.wavecal import waveio
from pypeit.core.wavecal import patterns
from pypeit.core.wavecal import wv_fitting
from pypeit.core.wavecal import wvutils
from pypeit.core import arc

from pypeit.core import pca
from pypeit import utils

from pypeit import msgs

from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch


def arc_fit_qa(waveFit, outfile=None, ids_only=False, title=None,
               log=True):
    """
    QA for Arc spectrum

    Args:
        waveFit (:class:`pypeit.core.wavecal.wv_fitting.WaveFit`):
        outfile (:obj:`str`, optional): Name of output file or 'show' to show on screen
        ids_only (bool, optional):
        title (:obj:`str`, optional):
        log (:obj:`bool`, optional):
            If True, use log scaling for the spectrum

    Returns:

    """
    plt.rcdefaults()
    plt.rcParams['font.family']= 'serif'

    arc_spec = waveFit['spec']

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

    # log is True by default, but if a large part of spectrum is < 0, the log plot will look very bad
    neg_values = np.where(arc_spec < 0)[0]
    if neg_values.size > 0.3 * len(arc_spec):
        log = False


    # Simple spectrum plot
    ax_spec = plt.subplot(gs[:,0])
    ax_spec.minorticks_on()
    ax_spec.plot(np.arange(len(arc_spec)), arc_spec)
    ymin, ymax = np.min(arc_spec), np.max(arc_spec)
    if log:
        ymax *= 4
        ymin = max(1., ymin)
    ysep = ymax*0.03
    yscl = (1.2, 1.5, 1.7)

    # Label all found lines
    for kk, x in enumerate(waveFit.tcent):
        ind_left = np.fmax(int(x)-2, 0)
        ind_righ = np.fmin(int(x)+2,arc_spec.size-1)
        yline = np.max(arc_spec[ind_left:ind_righ])
        # Tick mark
        if log:
            ax_spec.plot([x,x], [yline*yscl[0], yline*yscl[1]], '-', color='gray')
        else:
            ax_spec.plot([x,x], [yline+ysep*0.25, yline+ysep], '-', color='gray')

    # Label the ID'd lines
    for kk, x in enumerate(waveFit.pixel_fit):
        ind_left = np.fmax(int(x)-2, 0)
        ind_righ = np.fmin(int(x)+2,arc_spec.size-1)
        yline = np.max(arc_spec[ind_left:ind_righ])
        # Tick mark
        if log:
            ax_spec.plot([x,x], [yline*yscl[0], yline*yscl[1]], 'g-')
        else:
            ax_spec.plot([x,x], [yline+ysep*0.25, yline+ysep], 'g-')
        # label
        if log:
            ypos = yline*yscl[2]
        else:
            ypos = yline+ysep*1.3
        ax_spec.text(x, ypos, '{:s} {:g}'.format(waveFit.ions[kk],
                                                          waveFit.wave_fit[kk]),
                     ha='center', va='bottom',size=idfont,
                     rotation=90., color='green')

    # Axes
    ax_spec.set_xlim(0., len(arc_spec))
    if not log:
        ax_spec.set_ylim(1.05*ymin, ymax*1.2)
    else:
        ax_spec.set_ylim(ymin, ymax)
    ax_spec.set_xlabel('Pixel')
    ax_spec.set_ylabel('Flux')
    if log:
        ax_spec.set_yscale('log')

    # Title
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
    ax_fit.scatter(waveFit.pixel_fit,waveFit.wave_fit, marker='x')
    # Rejections?
    gpm = waveFit.pypeitfit.bool_gpm
    bpm = np.logical_not(gpm)
    if np.any(bpm):
        xrej = waveFit.pixel_fit[bpm]
        yrej = waveFit.wave_fit[bpm]
        ax_fit.scatter(xrej, yrej, marker='o', edgecolor='gray', facecolor='none')
    # Solution
    xval = np.arange(len(arc_spec))
    ax_fit.plot(xval, waveFit.wave_soln, 'r-')
    xmin, xmax = 0., len(arc_spec)
    ax_fit.set_xlim(xmin, xmax)
    ymin,ymax = np.min(waveFit.wave_soln)*.95,  np.max(waveFit.wave_soln)*1.05
    ax_fit.set_ylim((ymin, ymax))
    ax_fit.set_ylabel('Wavelength')
    ax_fit.get_xaxis().set_ticks([]) # Suppress labeling
    ax_fit.minorticks_on()
    ax_fit.tick_params(axis="y", which='both', right=True)

    # Stats
    wave_soln_fit = waveFit.pypeitfit.eval(waveFit.pixel_fit/waveFit.xnorm)#, 'legendre',minx=fit['fmin'], maxx=fit['fmax'])
    ax_fit.text(0.1*len(arc_spec), 0.90*ymin+(ymax-ymin),r'$\Delta\lambda$={:.3f}$\AA$ (per pix)'.format(waveFit.cen_disp), size='small')
    ax_fit.text(0.1*len(arc_spec), 0.80*ymin+(ymax-ymin),'RMS={:.3f} (pixels)'.format(waveFit.rms), size='small')
    # Arc Residuals
    ax_res = plt.subplot(gs[1,1])
    res = waveFit.wave_fit-wave_soln_fit
    ax_res.scatter(waveFit.pixel_fit[gpm], res[gpm]/waveFit.cen_disp, marker='x')
    ax_res.plot([xmin,xmax], [0.,0], 'k--')
    ax_res.set_xlim(xmin, xmax)
    ax_res.set_xlabel('Pixel')
    ax_res.set_ylabel('Residuals (Pix)')
    ax_res.minorticks_on()
    ax_res.tick_params(axis="y", which='both', right=True)

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





def reidentify(spec, spec_arxiv_in, wave_soln_arxiv_in, line_list, nreid_min, det_arxiv=None, detections=None, cc_thresh=0.8,cc_local_thresh = 0.8,
               match_toler=2.0, nlocal_cc=11, nonlinear_counts=1e10,sigdetect=5.0,fwhm=4.0,
               debug_xcorr=False, debug_reid=False, debug_peaks = False):
    """ Determine  a wavelength solution for a set of spectra based on archival wavelength solutions

    Parameters
    ----------
    spec:  float ndarray shape (nspec)
       Arc spectrum for which wavelength identifications are desired.

    spec_arxiv:  float ndarray shape (nspec, narxiv) or (nspec)
       Collection of archival arc spectra for which wavelength solution and line identifications are known

    wave_soln_arxiv:  float ndarray shape (nspec, narxiv) or (nspec)
       Wavelength solutions for the archival arc spectra spec_arxiv

    line_list: astropy table
       The arc line list used for thew wavelength solution in pypeit format.

    nreid_min: int
       Minimum number of times that a given candidate reidentified line must be properly matched with a line in the arxiv
       to be considered a good reidentification. If there is a lot of duplication in the arxiv of the spectra in question
       (i.e. multislit) set this to a number like 2-4. For echelle this depends on the number of solutions in the arxiv.
       For fixed format echelle (ESI, X-SHOOTER, NIRES) set this 1. For an echelle with a tiltable grating, it will depend
       on the number of solutions in the arxiv.

    Optional Parameters
    -------------------

    det_arxiv (optional):  dict, the dict has narxiv keys which are '0','1', ... up to str(narxiv-1). det_arxiv['0'] points to an
                an ndarray of size determined by the number of lines that were detected.

       Arc line pixel locations in the spec_arxiv spectra that were used in combination with line identifications from the
       line list to determine the wavelength solution wave_soln_arxiv.

    detections: float ndarray, default = None
       An array containing the pixel centroids of the lines in the arc as computed by the pypeit.core.arc.detect_lines
       code. If this is set to None, the line detection will be run inside the code.

    cc_thresh: float, default = 0.8
       Threshold for the *global* cross-correlation coefficient between an input spectrum and member of the archive required to
       attempt reidentification. Spectra from the archive with a lower cross-correlation are not used for reidentification

    cc_local_thresh: float, default = 0.8
       Threshold for the *local* cross-correlation coefficient, evaluated at each reidentified line,  between an input
       spectrum and the shifted and stretched archive spectrum above which a line must be to be considered a good line for
       reidentification. The local cross-correlation is evaluated at each candidate reidentified line
       (using a window of nlocal_cc), and is then used to score the the reidentified lines to arrive at the final set of
       good reidentifications

    match_toler: float, default = 2.0
       Matching tolerance in pixels for a line reidentification. A good line match must match within this tolerance to the
       the shifted and stretched archive spectrum, and the archive wavelength solution at this match must be within
       match_toler dispersion elements from the line in line list.

    n_local_cc: int, defualt = 11
       Size of pixel window used for local cross-correlation computation for each arc line. If not an odd number one will
       be added to it to make it odd.


    debug_xcorr: bool, default = False
       Show plots useful for debugging the cross-correlation used for shift/stretch computation

    debug_reid: bool, default = False
       Show plots useful for debugging the line reidentification

    Returns
    -------
    (detections, spec_cont_sub, patt_dict)

    detections: ndarray,
       Pixel locations of arc lines detected.
    spec_cont_sub: ndarray
       Array of continuum subtracted arc spectra
    patt_dict: dict
       Arc lines pattern dictionary with some information about the IDs as well as the cross-correlation values

    Revision History
    ----------------
    November 2018 by J.F. Hennawi. Built from an initial version of cross_match code written by Ryan Cooke.
    """

    # Determine the seed for scipy.optimize.differential_evolution optimizer. Just take the sum of all the elements
    # and round that to an integer
    
    seed = np.fmin(int(np.abs(np.sum(spec[np.isfinite(spec)]))),2**32-1)
    random_state = np.random.RandomState(seed = seed)

    nlocal_cc_odd = nlocal_cc + 1 if nlocal_cc % 2 == 0 else nlocal_cc
    window = 1.0/nlocal_cc_odd* np.ones(nlocal_cc_odd)

    # Generate the wavelengths from the line list and sort
    wvdata = np.array(line_list['wave'].data)  # Removes mask if any
    wvdata.sort()
    # Determine whether wavelengths correlate or anti-correlation with pixels for patt_dicts. This is not used
    # but just comptued for compatibility

    # Do some input checking
    if spec.ndim == 1:
        nspec = spec.size
    else:
        msgs.error('spec must be a one dimensional numpy array ')

    if spec_arxiv_in.ndim != wave_soln_arxiv_in.ndim:
        msgs.error('spec arxiv and wave_soln_arxiv must have the same dimensions')

    if spec_arxiv_in.ndim == 1:
        spec_arxiv1 = spec_arxiv_in.reshape(spec_arxiv_in.size,1)
        wave_soln_arxiv1 = wave_soln_arxiv_in.reshape(wave_soln_arxiv_in.size,1)
    elif spec_arxiv_in.ndim == 2:
        spec_arxiv1 = spec_arxiv_in.copy()
        wave_soln_arxiv1 = wave_soln_arxiv_in.copy()
    else:
        msgs.error('Unrecognized shape for spec_arxiv. It must be either a one dimensional or two dimensional numpy array')

    # TODO: JFH I would like to take these calls out. This reidentify code should only ever be run by comparing
    # data with the same binning. That would then allow me to drop the requirement that this code operate
    # on arrays of same number of pixels. I'm a big confused on how that interacts with stretch though so postponing
    # these changes for now.
    spec_arxiv = arc.resize_spec(spec_arxiv1, nspec)
    wave_soln_arxiv = arc.resize_spec(wave_soln_arxiv1, nspec)

    nspec_arxiv, narxiv = spec_arxiv.shape

    this_soln = wave_soln_arxiv[:,0]
    sign = 1 if (this_soln[this_soln.size // 2] > this_soln[this_soln.size // 2 - 1]) else -1

    xrng = np.arange(nspec)
    if nspec_arxiv != nspec:
        msgs.error('Spectrum sizes do not match. Something is very wrong!')

    # Search for lines no matter what to continuum subtract the input arc
    tcent, ecent, cut_tcent, icut, spec_cont_sub = wvutils.arc_lines_from_spec(
        spec, sigdetect=sigdetect,nonlinear_counts=nonlinear_counts, fwhm = fwhm, debug = debug_peaks)
    # If the detections were not passed in measure them
    if detections is None:
        detections = tcent[icut]

    spec_arxiv_cont_sub = np.zeros_like(spec_arxiv)

    # Search for lines no matter what to continuum subtract the arxiv arc, also determine the central wavelength and
    # dispersion of wavelength arxiv
    det_arxiv1 = {}
    for iarxiv in range(narxiv):
        tcent_arxiv, ecent_arxiv, cut_tcent_arxiv, icut_arxiv, spec_cont_sub_now = wvutils.arc_lines_from_spec(
            spec_arxiv[:,iarxiv], sigdetect=sigdetect,nonlinear_counts=nonlinear_counts, fwhm = fwhm, debug = debug_peaks)
        spec_arxiv_cont_sub[:,iarxiv] = spec_cont_sub_now
        det_arxiv1[str(iarxiv)] = tcent_arxiv[icut_arxiv]

    if det_arxiv is None:
        det_arxiv = det_arxiv1

    wvc_arxiv = np.zeros(narxiv, dtype=float)
    disp_arxiv = np.zeros(narxiv, dtype=float)
    # Determine the central wavelength and dispersion of wavelength arxiv
    for iarxiv in range(narxiv):
        wvc_arxiv[iarxiv] = wave_soln_arxiv[nspec//2, iarxiv]
        igood = wave_soln_arxiv[:,iarxiv] > 1.0
        disp_arxiv[iarxiv] = np.median(wave_soln_arxiv[igood,iarxiv] - np.roll(wave_soln_arxiv[igood,iarxiv], 1))

    marker_tuple = ('o','v','<','>','8','s','p','P','*','X','D','d','x')
    color_tuple = ('black','green','red','cyan','magenta','blue','darkorange','yellow','dodgerblue','purple','lightgreen','cornflowerblue')
    marker = itertools.cycle(marker_tuple)
    colors = itertools.cycle(color_tuple)

    # Cross-correlate with each arxiv spectrum to identify lines
    line_indx = np.array([], dtype=int)
    det_indx = np.array([], dtype=int)
    line_cc = np.array([], dtype=float)
    line_iarxiv = np.array([], dtype=int)
    wcen = np.zeros(narxiv)
    disp = np.zeros(narxiv)
    shift_vec = np.zeros(narxiv)
    stretch_vec = np.zeros(narxiv)
    ccorr_vec = np.zeros(narxiv)
    for iarxiv in range(narxiv):
        msgs.info('Cross-correlating with arxiv slit # {:d}'.format(iarxiv))
        this_det_arxiv = det_arxiv[str(iarxiv)]
        # Match the peaks between the two spectra. This code attempts to compute the stretch if cc > cc_thresh
        success, shift_vec[iarxiv], stretch_vec[iarxiv], ccorr_vec[iarxiv], _, _ = \
            wvutils.xcorr_shift_stretch(spec_cont_sub, spec_arxiv[:, iarxiv],
                                        cc_thresh=cc_thresh, fwhm=fwhm, seed=random_state,
                                        debug=debug_xcorr)
        # If cc < cc_thresh or if this optimization failed, don't reidentify from this arxiv spectrum
        if success != 1:
            continue
        # Estimate wcen and disp for this slit based on its shift/stretch relative to the archive slit
        disp[iarxiv] = disp_arxiv[iarxiv] / stretch_vec[iarxiv]
        wcen[iarxiv] = wvc_arxiv[iarxiv] - shift_vec[iarxiv]*disp[iarxiv]
        # For each peak in the arxiv spectrum, identify the corresponding peaks in the input spectrum. Do this by
        # transforming these arxiv slit line pixel locations into the (shifted and stretched) input spectrum frame
        det_arxiv_ss = this_det_arxiv*stretch_vec[iarxiv] + shift_vec[iarxiv]
        spec_arxiv_ss = wvutils.shift_and_stretch(spec_arxiv[:, iarxiv], shift_vec[iarxiv], stretch_vec[iarxiv])

        if debug_xcorr:
            plt.figure(figsize=(14, 6))
            tampl_slit = np.interp(detections, xrng, spec_cont_sub)
            plt.plot(xrng, spec_cont_sub, color='red', drawstyle='steps-mid', label='input arc',linewidth=1.0, zorder=10)
            plt.plot(detections, tampl_slit, 'r.', markersize=10.0, label='input arc lines', zorder=10)
            tampl_arxiv = np.interp(this_det_arxiv, xrng, spec_arxiv[:, iarxiv])
            plt.plot(xrng, spec_arxiv[:, iarxiv], color='black', drawstyle='steps-mid', linestyle=':',
                     label='arxiv arc', linewidth=0.5)
            plt.plot(this_det_arxiv, tampl_arxiv, 'k+', markersize=8.0, label='arxiv arc lines')
            # tampl_ss = np.interp(gsdet_ss, xrng, gdarc_ss)
            for iline in range(det_arxiv_ss.size):
                plt.plot([this_det_arxiv[iline], det_arxiv_ss[iline]], [tampl_arxiv[iline], tampl_arxiv[iline]],
                         color='cornflowerblue', linewidth=1.0)
            plt.plot(xrng, spec_arxiv_ss, color='black', drawstyle='steps-mid', label='arxiv arc shift/stretch',linewidth=1.0)
            plt.plot(det_arxiv_ss, tampl_arxiv, 'k.', markersize=10.0, label='predicted arxiv arc lines')
            plt.title(
                'Cross-correlation of input slit and arxiv slit # {:d}'.format(iarxiv + 1) +
                ': ccor = {:5.3f}'.format(ccorr_vec[iarxiv]) +
                ', shift = {:6.1f}'.format(shift_vec[iarxiv]) +
                ', stretch = {:5.4f}'.format(stretch_vec[iarxiv]) +
                ', wv_cen = {:7.1f}'.format(wcen[iarxiv]) +
                ', disp = {:5.3f}'.format(disp[iarxiv]))
            plt.ylim(1.2*spec_cont_sub.min(), 1.5 *spec_cont_sub.max())
            plt.legend()
            plt.show()


        # Calculate wavelengths for all of the this_det_arxiv detections. This step could in principle be done more accurately
        # with the polynomial solution itself, but the differences are 1e-12 of a pixel, and this interpolate of the tabulated
        # solution makes the code more general.
        wvval_arxiv = (scipy.interpolate.interp1d(xrng, wave_soln_arxiv[:, iarxiv], kind='cubic'))(this_det_arxiv)

        # Compute a "local" zero lag correlation of the slit spectrum and the shifted and stretch arxiv spectrum over a
        # a nlocal_cc_odd long segment of spectrum. We will then uses spectral similarity as a further criteria to
        # decide which lines are good matches
        prod_smooth = scipy.ndimage.convolve1d(spec_cont_sub*spec_arxiv_ss, window)
        spec2_smooth = scipy.ndimage.convolve1d(spec_cont_sub**2, window)
        arxiv2_smooth = scipy.ndimage.convolve1d(spec_arxiv_ss**2, window)
        denom = np.sqrt(spec2_smooth*arxiv2_smooth)
        corr_local = np.zeros_like(denom)
        corr_local[denom > 0] = prod_smooth[denom > 0]/denom[denom > 0]
        corr_local[denom == 0.0] = -1.0

        # Loop over the current slit line pixel detections and find the nearest arxiv spectrum line
        # JFH added this if statement to prevent crashes for cases where no arc lines where found. This is because
        # full_template keeps passing in tiny snippets of mostly junk padded spectra that cause all kind of crashes.
        # A better approach would be to fix full_template so as to not enter reidentify unless the "arxiv_arcs"
        # are not almost entirely zero padded snippets.
        if det_arxiv_ss.size > 0:
            for iline in range(detections.size):
                # match to pixel in shifted/stretch arxiv spectrum
                pdiff = np.abs(detections[iline] - det_arxiv_ss)
                bstpx = np.argmin(pdiff)
                # If a match is found within 2 pixels, consider this a successful match
                if pdiff[bstpx] < match_toler:
                    # Using the arxiv arc wavelength solution, search for the nearest line in the line list
                    bstwv = np.abs(wvdata - wvval_arxiv[bstpx])
                    # This is a good wavelength match if it is within match_toler disperion elements
                    if bstwv[np.argmin(bstwv)] < match_toler*disp_arxiv[iarxiv]:
                        line_indx = np.append(line_indx, np.argmin(bstwv))  # index in the line list array wvdata of this match
                        det_indx = np.append(det_indx, iline)             # index of this line in the detected line array detections
                        line_cc = np.append(line_cc,np.interp(detections[iline],xrng,corr_local)) # local cross-correlation at this match
                        line_iarxiv = np.append(line_iarxiv,iarxiv)

    narxiv_used = np.sum(wcen != 0.0)

    # Initialise the patterns dictionary, sigdetect not used anywhere
    if (narxiv_used == 0) or (len(np.unique(line_indx)) < 3):
        patt_dict_slit = patterns.empty_patt_dict(detections.size)
        patt_dict_slit['sigdetect'] = sigdetect
        return detections, spec_cont_sub, patt_dict_slit


    # Finalize the best guess of each line
    patt_dict_slit = patterns.solve_xcorr(detections, wvdata, det_indx, line_indx, line_cc,nreid_min=nreid_min,cc_local_thresh=cc_local_thresh)
    patt_dict_slit['sign'] = sign # This is not used anywhere
    patt_dict_slit['bwv'] = np.median(wcen[wcen != 0.0])
    patt_dict_slit['bdisp'] = np.median(disp[disp != 0.0])
    patt_dict_slit['sigdetect'] = sigdetect


    if debug_reid:
        plt.figure(figsize=(14, 6))
        # Plot a summary of the local x-correlation values for each line on each slit
        for iarxiv in range(narxiv):
            # Only plot those that we actually tried to reidentify (i.e. above cc_thresh)
            if wcen[iarxiv] != 0.0:
                this_iarxiv = line_iarxiv == iarxiv
                plt.plot(wvdata[line_indx[this_iarxiv]], line_cc[this_iarxiv], marker=next(marker), color=next(colors),
                         linestyle='', markersize=5.0, label='arxiv slit={:d}'.format(iarxiv))

        plt.hlines(cc_local_thresh, wvdata[line_indx].min(), wvdata[line_indx].max(), color='red', linestyle='--',
                   label='Local xcorr threshhold')
        plt.title('Local x-correlation for reidentified lines from narxiv_used={:d}'.format(narxiv_used) +
                  ' arxiv slits. Requirement: nreid_min={:d}'.format(nreid_min) + ' matches > threshold')
        plt.xlabel('wavelength from line list')
        plt.ylabel('Local x-correlation coefficient')
        # plt.ylim((0.0, 1.2))
        plt.legend()
        plt.show()
        # QA Plot ofthe reidentifications
        match_qa(spec_cont_sub, detections, line_list, patt_dict_slit['IDs'], patt_dict_slit['scores'])

    # Use only the perfect IDs
    iperfect = np.array(patt_dict_slit['scores']) != 'Perfect'
    patt_dict_slit['mask'][iperfect] = False
    patt_dict_slit['nmatch'] = np.sum(patt_dict_slit['mask'])
    if patt_dict_slit['nmatch'] < 3:
        patt_dict_slit['acceptable'] = False

    return detections, spec_cont_sub, patt_dict_slit


def measure_fwhm(spec, sigdetect=10., fwhm=5.):
    """
    Measure the arc lines FWHM, i.e, approximate spectral resolution

    Args:
        spec (`numpy.ndarray`_):
            Arc spectrum from a single slit.
        sigdetect (:obj:`float`, optional):
            Sigma threshold above fluctuations for arc-line detection.
            Used by :func:`pypeit.core.arc.detect_lines`.
        fwhm (:obj:`float`, optional):
            Number of pixels per fwhm resolution element.
            Used by :func:`pypeit.core.arc.detect_lines`.

    Returns:
        :obj:`float`: Measured arc lines FWHM in binned pixels of the input arc image
    """

    # Determine the lines FWHM, i.e, approximate spectral resolution
    #  This may only be recorded and not used by the algorithms
    _, _, _, wdth, _, best, _, nsig = arc.detect_lines(spec, sigdetect=sigdetect, fwhm=fwhm)
    # 1sigma Gaussian widths of the line detections
    wdth = wdth[best]
    # significance of each line detected
    nsig = nsig[best]
    # Nsigma (significance) threshold. We use only lines that have the highest significance
    # We start with nsig_thrshd of 500 and iteratively reduce it if there are not more than 6 lines
    nsig_thrshd = 500.
    measured_fwhm = None
    while nsig_thrshd >= sigdetect:
        if wdth[nsig >= nsig_thrshd].size > 6:
            # compute average `wdth`
            mean, med, _ = astropy.stats.sigma_clipped_stats(
                wdth[nsig >= nsig_thrshd], sigma_lower=2.0, sigma_upper=2.0
            )
            # FWHM in pixels
            measured_fwhm = med * (2 * np.sqrt(2 * np.log(2)))
            break
        nsig_thrshd -= sigdetect/2.

    return measured_fwhm


def set_fwhm(par, measured_fwhm=None):
    """
    Set the value of the arc lines FWHM by choosing between the provided parset
    and the measured_fwhm

    Args:
        par (:class:`~pypeit.par.pypeitpar.WaveSolutionPar`):
            Key parameters that drive the behavior of the
            wavelength-solution algorithms.
        measured_fwhm (:obj:`float`):
            Measured arc lines FWHM in binned pixels of the input arc image

    Returns:
       :obj:`float`: Chosen arc lines FWHM in binned pixels of the input arc image
    """

    # Set FWHM for the methods that follow
    if par['fwhm_fromlines'] is False:
        fwhm = par['fwhm']
        msgs.info(f"User-provided arc lines FWHM: {fwhm:.1f} pixels")
    elif measured_fwhm is None:
        fwhm = par['fwhm']
        msgs.warn(f"Assumed arc lines FWHM: {fwhm:.1f} pixels")
    else:
        fwhm = measured_fwhm
        msgs.info(f"Measured arc lines FWHM: {fwhm:.1f} pixels")

    return fwhm


def full_template(spec, lamps, par, ok_mask, det, binspectral, nsnippet=2,
                  measured_fwhms=None, debug_xcorr=False, debug_reid=False,
                  x_percentile=50., template_dict=None, debug=False, 
                  nonlinear_counts=1e10):
    """
    Method of wavelength calibration using a single, comprehensive template spectrum

    The steps are:
      1. Load the template and rebin, as necessary
      2. Cross-correlate input spectrum and template to find the shift between the two
      3. Loop on snippets of the input spectrum to ID lines using reidentify()
      4. Fit with fitting.iterative_fitting()

    Args:
        spec: ndarray (nspec, nslit)
          Spectra to be calibrated
        lamps : :obj:`list`
            List of arc lamps to be used for wavelength calibration.
            E.g., ['ArI','NeI','KrI','XeI']
        par: WavelengthSolutionPar ParSet
          Calibration parameters
        ok_mask: ndarray, bool
          Mask of indices of good slits
        det: int
          Detector index
        binspectral: int
          Binning of the input arc in the spectral dimension
        nsnippet: int, optional
          Number of snippets to chop the input spectrum into when ID'ing lines
          This deals with differences due to non-linearity between the template
          and input spectrum.
        measured_fwhms: ndarray, optional
            Array of FWHM (in pixels) measured from the arc lines. Shape (nslit,)
        x_percentile: float, optional
          Passed to reidentify to reduce the dynamic range of arc line amplitudes
        template_dict (dict, optional): Dict containing tempmlate items, largely for development

    Returns:
        wvcalib: dict
          Dict of wavelength calibration solutions

    """
    #debug = True
    #debug_xcorr = True
    #debug_reid = True
    # Load line lists
    if 'ThAr' in lamps:
        line_lists_all = waveio.load_line_lists(lamps)
        line_lists = line_lists_all[np.where(line_lists_all['ion'] != 'UNKNWN')]
    else:
        line_lists = waveio.load_line_lists(lamps)

    # Load template
    if template_dict is None:
        temp_wv, temp_spec, temp_bin = waveio.load_template(
            par['reid_arxiv'], det, wvrng=par['wvrng_arxiv'])
    else:
        temp_wv = template_dict['wave']
        temp_spec = template_dict['spec']
        temp_bin = template_dict['bin']

    # Deal with binning (not yet tested)
    if binspectral != temp_bin:
        msgs.info("Resizing the template due to different binning.")
        new_npix = int(temp_wv.size * temp_bin / binspectral)
        temp_wv = arc.resize_spec(temp_wv, new_npix)
        temp_spec = arc.resize_spec(temp_spec, new_npix)

    # Dimensions
    if spec.ndim == 2:
        nspec, nslits = spec.shape
    elif spec.ndim == 1:
        nspec = spec.size
        nslits = 1
        spec = np.reshape(spec, (nspec,1))

    # Loop on slits
    wvcalib = {}
    for slit in range(nslits):
        # Sigdetect
        sigdetect = wvutils.parse_param(par, 'sigdetect', slit)
        # Check
        if slit not in ok_mask:
            wvcalib[str(slit)] = None
            continue
        msgs.info("Processing slit {}".format(slit))
        msgs.info("Using sigdetect = {}".format(sigdetect))
        # Grab the observed arc spectrum
        obs_spec_i = spec[:,slit]
        # get FWHM for this slit
        fwhm = set_fwhm(par, measured_fwhm=measured_fwhms[slit])

        # Find the shift
        ncomb = temp_spec.size
        # Remove the continuum before adding the padding to obs_spec_i
        _, _, _, _, obs_spec_cont_sub = wvutils.arc_lines_from_spec(obs_spec_i)
        _, _, _, _, templ_spec_cont_sub = wvutils.arc_lines_from_spec(temp_spec)
        # Pad
        pad_spec = np.zeros_like(temp_spec)
        nspec = len(obs_spec_i)
        npad = ncomb - nspec
        if npad > 0:    # Pad the input spectrum
            pad_spec[npad // 2:npad // 2 + len(obs_spec_i)] = obs_spec_cont_sub
            tspec = templ_spec_cont_sub
        elif npad < 0:  # Pad the template!
            pad_spec = obs_spec_cont_sub
            npad *= -1
            tspec = np.zeros(nspec)
            tspec[npad // 2:npad // 2 + ncomb] = templ_spec_cont_sub
        else:  # No padding necessary
            pad_spec = obs_spec_cont_sub
            tspec = templ_spec_cont_sub
        # Cross-correlate
        shift_cc, corr_cc = wvutils.xcorr_shift(tspec, pad_spec, debug=debug, fwhm=fwhm, percent_ceil=x_percentile)
        #shift_cc, corr_cc = wvutils.xcorr_shift(temp_spec, pspec, debug=debug, percent_ceil=x_percentile)
        msgs.info("Shift = {}; cc = {}".format(shift_cc, corr_cc))
        if debug:
            xvals = np.arange(tspec.size)
            plt.clf()
            ax = plt.gca()
            #
            ax.plot(xvals, tspec, label='template')  # Template
            ax.plot(xvals, np.roll(pspec, int(shift_cc)), 'k', label='input')  # Input
            ax.legend()
            plt.show()
            #embed(header='909 autoid')
        i0 = npad // 2 + int(shift_cc)

        # Generate the template snippet
        if i0 < 0: # Pad?
            mspec = np.concatenate([np.zeros(-1*i0), temp_spec[0:i0+nspec]])
            mwv = np.concatenate([np.zeros(-1*i0), temp_wv[0:i0+nspec]])
        elif (i0+nspec) > temp_spec.size: # Pad?
            mspec = np.concatenate([temp_spec[i0:], np.zeros(nspec-temp_spec.size+i0)])
            mwv = np.concatenate([temp_wv[i0:], np.zeros(nspec-temp_spec.size+i0)])
        else: # Don't pad
            mspec = temp_spec[i0:i0 + nspec]
            mwv = temp_wv[i0:i0 + nspec]

        # Loop on snippets
        nsub = obs_spec_i.size // nsnippet
        sv_det, sv_IDs = [], []
        for kk in range(nsnippet):
            # Construct
            j0 = nsub * kk
            j1 = min(nsub*(kk+1), obs_spec_i.size)
            tsnippet = obs_spec_i[j0:j1]
            msnippet = mspec[j0:j1]
            mwvsnippet = mwv[j0:j1]
            # TODO: JFH This continue statement deals with the case when the msnippet derives from *entirely* zero-padded
            #  pixels, and allows the code to continue with crashing. This code is constantly causing reidentify to crash
            #  by passing in these junk snippets that are almost entirely zero-padded for large shifts. We should
            #  be checking for this intelligently rather than constantly calling reidentify with basically junk arxiv
            #  spectral snippets.
            if not np.any(msnippet):
                continue
            # TODO -- JXP
            #  should we use par['cc_thresh'] instead of hard-coding cc_thresh??
            # Run reidentify
            detections, spec_cont_sub, patt_dict = reidentify(tsnippet, msnippet, mwvsnippet,
                                                              line_lists, 1, debug_xcorr=debug_xcorr,
                                                              sigdetect=sigdetect,
                                                              nonlinear_counts=nonlinear_counts,
                                                              debug_reid=debug_reid,  # verbose=True,
                                                              match_toler=par['match_toler'],
                                                              cc_thresh=0.1, fwhm=fwhm)
            # Deal with IDs
            sv_det.append(j0 + detections)
            try:
                sv_IDs.append(patt_dict['IDs'])
            except KeyError:
                msgs.warn("Failed to perform wavelength calibration in reidentify..")
                sv_IDs.append(np.zeros_like(detections))
            else:
                # Save now in case the next one barfs
                bdisp = patt_dict['bdisp']

        # Collate and proceed
        dets = np.concatenate(sv_det)
        IDs = np.concatenate(sv_IDs)
        gd_det = np.where(IDs > 0.)[0]
        if len(gd_det) < 4:
            msgs.warn("Not enough useful IDs")
            wvcalib[str(slit)] = None
            continue
        # Fit
        try:
            final_fit = wv_fitting.iterative_fitting(obs_spec_i, dets, gd_det,
                                              IDs[gd_det], line_lists, bdisp,
                                              verbose=False, n_first=par['n_first'],
                                              match_toler=par['match_toler'],
                                              func=par['func'],
                                              n_final=par['n_final'],
                                              sigrej_first=par['sigrej_first'],
                                              sigrej_final=par['sigrej_final'])
        except TypeError:
            #embed(header='974 of autoid')
            wvcalib[str(slit)] = None
        else:
            wvcalib[str(slit)] = copy.deepcopy(final_fit)
    # Finish
    return wvcalib

def echelle_wvcalib(spec, orders, spec_arxiv, wave_arxiv, lamps, par, ok_mask=None, use_unknowns=True, debug_all=False,
                    debug_peaks=False, debug_xcorr=False, debug_reid=False, debug_fits=False, nonlinear_counts=1e10):
    r"""
    Algorithm to wavelength calibrate echelle data based on a predicted or archived wavelength solution

    Parameters
    ----------
    spec :  `numpy.ndarray_`, shape=(nspec, norders)
        Array of arc spectra for each order for which wavelength solutions are
        desired.
    orders : `numpy.ndarray_`
        Order numbers for the provided spectra. Used to match against
        the relevant archived spectrum for echelle spectrographs.
        Shape must be :math:`(N_{\rm orders},)`
    spec_arxiv :  `numpy.ndarray_`, shape=(nspec, narxiv) or (nspec)
        Collection of archival arc spectra for which wavelength solution and line identifications are known
    wave_arxiv:  float ndarray shape (nspec, narxiv) or (nspec)
        Wavelength solutions for the archival arc spectra spec_arxiv
    lamps : :obj:`list`
        List of arc lamps to be used for wavelength calibration.
        E.g., ['ArI','NeI','KrI','XeI']
    par : :class:`~pypeit.par.pypeitpar.WaveSolutionPar`
        Key parameters that drive the behavior of the
        wavelength-solution algorithms.
    ok_mask : `numpy.ndarray`, optional
        Integer array with the list of valid spectra ``spec`` to use.
        If None, all spectra are used.
    use_unknowns : bool, default = True, optional
        If True, arc lines that are known to be present in the
        spectra, but have not been attributed to an element+ion, will
        be included in the fit.
    debug_all: :obj:`bool`, optional
        Convenience parameter that turns on all debugging. Setting
        ``debug_all`` to True is equivalent to setting
        ``debug_peaks``, ``debug_xcorr``, ``debug_reid``, and
        ``debug_fits`` to True.
    debug_peaks : :obj:`bool`, optional
        Debug the line identification in the arcline spectra. See
        ``debug`` parameter in
        func:`pypeit.core.wavecal.wvutils.arc_lines_from_spec`.
    debug_xcorr: bool, default = False, optional
        Show plots useful for debugging the cross-correlation used
        for shift/stretch computation.
    debug_reid: bool, default = False, optional
        Show plots useful for debugging the line reidentification
    debug_fits : :obj:`bool`, optional
        Show the arc-line fit debugging plot. See :func:`arc_fit_qa`.
    nonlinear_counts: float, default = 1e10
        For arc line detection: Arc lines above this saturation
        threshold are not used in wavelength solution fits because
        they cannot be accurately centroided

    Returns
    -------
    all_patt_dict: dict
       Arc lines pattern dictionary with some information about the IDs as well as the cross-correlation values
    wv_calib: dict
       Dictionary containing the wavelength solution for each order

    """

    # TODO: Perform detailed checking of the input

    # Check input
    if not isinstance(par, pypeitpar.WavelengthSolutionPar):
        msgs.error('Input parameters must be provided by a WavelengthSolutionPar object.')


    if spec.ndim != 2:
        msgs.error('Input spec must be a 2D numpy array!')

    nspec, norders = spec.shape

    if orders.size != norders:
        msgs.error('Number of provided orders does not match the number of provided spectra.')

    # Mask info
    ok_mask = np.arange(norders) if ok_mask is None else ok_mask
    if np.amax(ok_mask) >= norders:
        msgs.error('Spectrum selected by ok_mask is beyond the limits of the provided '
                   'spec array.')

    # Load the line lists
    if 'ThAr' in lamps:
        line_lists_all = waveio.load_line_lists(lamps)
        line_lists = line_lists_all[np.where(line_lists_all['ion'] != 'UNKNWN')]
        unknwns = line_lists_all[np.where(line_lists_all['ion'] == 'UNKNWN')]
    else:
        line_lists = waveio.load_line_lists(lamps)
        unknwns = waveio.load_unknown_list(lamps)

    tot_line_list = astropy.table.vstack([line_lists, unknwns]) if use_unknowns else line_lists

    # Array to hold continuum subtracted arcs
    spec_cont_sub = np.zeros_like(spec)

    # These are the final outputs
    all_patt_dict = {}
    detections = {}
    wv_calib = {}
    bad_orders = np.array([], dtype=int)
    # Reidentify each slit, and perform a fit
    for iord in range(norders):
        # ToDO should we still be populating wave_calib with an empty dict here?
        if iord not in ok_mask:
            wv_calib[str(iord)] = None
            continue
        msgs.info('Reidentifying and fitting Order = {0:d}, which is {1:d}/{2:d}'.format(orders[iord], iord, norders - 1))
        sigdetect = wvutils.parse_param(par, 'sigdetect', iord)
        cc_thresh = wvutils.parse_param(par, 'cc_thresh', iord)
        rms_threshold = wvutils.parse_param(par, 'rms_threshold', iord)
        msgs.info("Using sigdetect =  {}".format(sigdetect))
        msgs.info("Using rms_threshold =  {}".format(rms_threshold))
        detections[str(iord)], spec_cont_sub[:, iord], all_patt_dict[str(iord)] = reidentify(
            spec[:, iord], spec_arxiv[:, iord], wave_arxiv[:, iord], tot_line_list, par['nreid_min'],
            cc_thresh=cc_thresh, match_toler=par['match_toler'],
            cc_local_thresh=par['cc_local_thresh'], nlocal_cc=par['nlocal_cc'],
            nonlinear_counts=nonlinear_counts, sigdetect=sigdetect, fwhm=par['fwhm'],
            debug_peaks=(debug_peaks or debug_all),
            debug_xcorr=(debug_xcorr or debug_all),
            debug_reid=(debug_reid or debug_all))
        # Check if an acceptable reidentification solution was found
        if not all_patt_dict[str(iord)]['acceptable']:
            wv_calib[str(iord)] = None
            bad_orders = np.append(bad_orders, iord)
            continue

        # Perform the fit
        n_final = wvutils.parse_param(par, 'n_final', iord)
        final_fit = wv_fitting.fit_slit(spec_cont_sub[:, iord], all_patt_dict[str(iord)], detections[str(iord)],
                                        tot_line_list, match_toler=par['match_toler'], func=par['func'],
                                        n_first=par['n_first'],
                                        sigrej_first=par['sigrej_first'], n_final=n_final, sigrej_final=par['sigrej_final'])

        # Did the fit succeed?
        if final_fit is None:
            # This pattern wasn't good enough
            wv_calib[str(iord)] = None
            bad_orders = np.append(bad_orders, iord)
            continue
        # Is the RMS below the threshold?
        if final_fit['rms'] > rms_threshold:
            msgs.warn('---------------------------------------------------' + msgs.newline() +
                      'Reidentify report for slit {0:d}/{1:d}:'.format(iord, norders - 1) + msgs.newline() +
                      '  Poor RMS ({0:.3f})! Need to add additional spectra to arxiv to improve fits'.format(
                          final_fit['rms']) + msgs.newline() +
                      '---------------------------------------------------')
            bad_orders = np.append(bad_orders, iord)
            # Note this result in new_bad_orders, but store the solution since this might be the best possible

        # Add the patt_dict and wv_calib to the output dicts
        wv_calib[str(iord)] = copy.deepcopy(final_fit)
        if (debug_fits or debug_all):
            arc_fit_qa(wv_calib[str(iord)], title='Silt: {}'.format(str(iord)))

    # Print the final report of all lines
    report_final(norders, all_patt_dict, detections, wv_calib, ok_mask, bad_orders)

    return all_patt_dict, wv_calib


def report_final(nslits, all_patt_dict, detections, wv_calib, ok_mask, bad_slits):
    """
    Print out the final report for wavelength calibration

    Args:
        nslits (int):
            Number of slits or ders
        all_patt_dict (dict):
            Dictionary containing reidentification information.
        detections (dict):
            Dictionary containing the lines that were detected.
        wv_calib (dict):
            Dictionary holding the wavelength solutions for each slit/orders
        ok_mask (ndarray, bool):
            Mask of indices of good slits
        bad_slits (ndarray, bool):
            List of slits that are bad
    """
    for slit in range(nslits):
        # Prepare a message for bad wavelength solutions
        badmsg = '---------------------------------------------------' + msgs.newline() + \
                 'Final report for slit {0:d}/{1:d}:'.format(slit, nslits) + msgs.newline() + \
                 '  Wavelength calibration not performed!'
        if slit not in ok_mask or slit in bad_slits or all_patt_dict[str(slit)] is None:
            msgs.warn(badmsg)
            continue
        st = str(slit)
        if wv_calib[st] is None:  # or len(wv_calib[st]) == 0:
            msgs.warn('Bad solution for slit: {}'.format(st))
            continue
        if all_patt_dict[st]['sign'] == +1:
            signtxt = 'correlate'
        else:
            signtxt = 'anitcorrelate'
        # Report
        cen_wave = wv_calib[st]['cen_wave']
        cen_disp = wv_calib[st]['cen_disp']
        msgs.info(msgs.newline() +
                  '---------------------------------------------------' + msgs.newline() +
                  'Final report for slit {0:d}/{1:d}:'.format(slit, nslits - 1) + msgs.newline() +
                  '  Pixels {:s} with wavelength'.format(signtxt) + msgs.newline() +
                  '  Number of lines detected      = {:d}'.format(detections[st].size) + msgs.newline() +
                  '  Number of lines that were fit = {:d}'.format(
                      len(wv_calib[st]['pixel_fit'])) + msgs.newline() +
                  '  Central wavelength            = {:g}A'.format(cen_wave) + msgs.newline() +
                  '  Central dispersion            = {:g}A/pix'.format(cen_disp) + msgs.newline() +
                  '  Central wave/disp             = {:g}'.format(cen_wave / cen_disp) + msgs.newline() +
                  '  Final RMS of fit              = {:g}'.format(wv_calib[st]['rms']))


class ArchiveReid:
    r"""
    Algorithm to wavelength calibrate spectroscopic data based on an
    archive of wavelength solutions.

    Parameters
    ----------
    spec :  float ndarray shape of (nspec, nslits) or (nspec)
        Array of arc spectra for which wavelength solutions are
        desired.
    lamps : :obj:`list`
        List of arc lamps to be used for wavelength calibration.
        E.g., ['ArI','NeI','KrI','XeI']
    par : :class:`~pypeit.par.pypeitpar.WaveSolutionPar`
        Key parameters that drive the behavior of the
        wavelength-solution algorithms.
    ech_fixed_format: bool
        Set to True if this is a fixed format echelle spectrograp. The code will then
        align the archive_arc and the extracted arc for each order for the reidentification.
    ok_mask : `numpy.ndarray`, optional
        Integer array with the list of valid spectra ``spec`` to use.
        If None, all spectra are used.
    measured_fwhms: ndarray, optional
        Array of FWHM (in pixels) measured from the arc lines. Shape (nslit,)
    use_unknowns : bool, default = True, optional
        If True, arc lines that are known to be present in the
        spectra, but have not been attributed to an element+ion, will
        be included in the fit.
    debug_all: :obj:`bool`, optional
        Convenience parameter that turns on all debugging. Setting
        ``debug_all`` to True is equivalent to setting
        ``debug_peaks``, ``debug_xcorr``, ``debug_reid``, and
        ``debug_fits`` to True.
    debug_peaks : :obj:`bool`, optional
        Debug the line identification in the arcline spectra. See
        ``debug`` parameter in
        func:`pypeit.core.wavecal.wvutils.arc_lines_from_spec`.
    debug_xcorr: bool, default = False, optional
        Show plots useful for debugging the cross-correlation used
        for shift/stretch computation.
    debug_reid: bool, default = False, optional
        Show plots useful for debugging the line reidentification
    debug_fits : :obj:`bool`, optional
        Show the arc-line fit debugging plot. See :func:`arc_fit_qa`.
    orders : `numpy.ndarray`, optional
        Order numbers for the provided spectra. Used to match against
        the relevant archived spectrum for echelle spectrographs.
        Shape must be :math:`(N_{\rm spec},)` and these *must* be
        provided if ech_fixed_format is True.
    nonlinear_counts: float, default = 1e10
        For arc line detection: Arc lines above this saturation
        threshold are not used in wavelength solution fits because
        they cannot be accurately centroided

    Attributes
    ----------
    debug_peaks : :obj:`bool`
        Debug the peak finding.

    .. todo::
        - Fill in the rest of the attributes.

    """
    # TODO: Because we're passing orders directly, we no longer need spectrograph...
    def __init__(self, spec, lamps, par, ech_fixed_format=False, ok_mask=None,  measured_fwhms=None, use_unknowns=True, debug_all=False,
                 debug_peaks=False, debug_xcorr=False, debug_reid=False, debug_fits=False,
                 orders=None, nonlinear_counts=1e10):

        # TODO: Perform detailed checking of the input

        # Check input
        if not isinstance(par, pypeitpar.WavelengthSolutionPar):
            msgs.error('Input parameters must be provided by a WavelengthSolutionPar object.')
        # TODO: Do we need ech_fix_format if we have
        # spectrograph.pypeline, assuming we keep passing spectrograph?
        if ech_fixed_format and orders is None:
            msgs.error('If the specrograph is a fixed-format echelle (ech_fix_format is True), '
                       'the orders must be provided.')

        # TODO: What does and does not need to be an attribute?

        # Debugging
        self.debug_peaks = debug_peaks or debug_all
        self.debug_xcorr = debug_xcorr or debug_all
        self.debug_reid = debug_reid or debug_all
        self.debug_fits = debug_fits or debug_all

        self.spec = spec
        if spec.ndim == 2:
            self.nspec, self.nslits = spec.shape
        elif spec.ndim == 1:
            self.nspec = spec.size
            self.nslits = 1
        else:
            msgs.error('Input spec must be a 1D or 2D numpy array!')

        if orders is not None and orders.size != self.nslits:
            msgs.error('Number of provided orders does not match the number of provided spectra.')

        self.par = par
        self.lamps = lamps
        self.use_unknowns = use_unknowns

        # Mask info
        self.ok_mask = np.arange(self.nslits) if ok_mask is None else ok_mask
        if np.amax(ok_mask) >= self.nslits:
            msgs.error('Spectrum selected by ok_mask is beyond the limits of the provided '
                       'spec array.')
        # List of bad slits
        self.bad_slits = []

        # Pull paramaters out of the parset
        # TODO: Why are we doing this?
        # Parameters for arc line detction
        self.nonlinear_counts = nonlinear_counts # self.par['nonlinear_counts']
        # Paramaters that govern reidentification
        self.reid_arxiv = self.par['reid_arxiv']
        self.nreid_min = self.par['nreid_min']
        self.nlocal_cc = self.par['nlocal_cc']
        self.cc_thresh = self.par['cc_thresh']
        self.cc_local_thresh = self.par['cc_local_thresh']

        # Paramters that govern wavelength solution fitting
        self.match_toler = self.par['match_toler']
        self.func = self.par['func']
        self.n_first= self.par['n_first']
        self.sigrej_first= self.par['sigrej_first']
        self.sigrej_final= self.par['sigrej_final']

        # Load the line lists
        if 'ThAr' in self.lamps:
            line_lists_all = waveio.load_line_lists(self.lamps)
            self.line_lists = line_lists_all[np.where(line_lists_all['ion'] != 'UNKNWN')]
            self.unknwns = line_lists_all[np.where(line_lists_all['ion'] == 'UNKNWN')]
        else:
            self.line_lists = waveio.load_line_lists(self.lamps)
            self.unknwns = waveio.load_unknown_list(self.lamps)

        self.tot_line_list = astropy.table.vstack([self.line_lists, self.unknwns]) if self.use_unknowns \
                                else self.line_lists

        # Read in the wv_calib_arxiv and pull out some relevant quantities
        # ToDO deal with different binnings!
        self.wv_calib_arxiv, self.par_arxiv = waveio.load_reid_arxiv(self.reid_arxiv)

        # Determine the number of spectra in the arxiv, check that it
        # matches nslits if this is fixed format.
        narxiv = len(self.wv_calib_arxiv)
        for key in self.wv_calib_arxiv.keys():
            try:
                test = int(key)
            except ValueError:
                narxiv -=1

        #if self.ech_fix_format and (self.nslits != narxiv):
        #    msgs.error('You have set ech_fix_format = True, but nslits={:d} != narxiv={:d}'.format(self.nslits,narxiv) + '.' +
        #               msgs.newline() + 'The number of orders identified does not match the number of solutions in the arxiv')
        #

        # Array to hold continuum subtracted arcs
        self.spec_cont_sub = np.zeros_like(self.spec)

        nspec_arxiv = self.wv_calib_arxiv['0']['spec'].size
        self.spec_arxiv = np.zeros((nspec_arxiv, narxiv))
        self.wave_soln_arxiv = np.zeros((nspec_arxiv, narxiv))
        self.det_arxiv = {}
        for iarxiv in range(narxiv):
            self.spec_arxiv[:, iarxiv] = self.wv_calib_arxiv[str(iarxiv)]['spec']
            self.wave_soln_arxiv[:, iarxiv] = self.wv_calib_arxiv[str(iarxiv)]['wave_soln']
        # arxiv orders (echelle only)
        if ech_fixed_format:
            arxiv_orders = []
            for iarxiv in range(narxiv):
                arxiv_orders.append(self.wv_calib_arxiv[str(iarxiv)]['order'])
#            orders, _ = self.spectrograph.slit2order(slit_spat_pos)

        ind_arxiv = np.arange(narxiv, dtype=int)
        # These are the final outputs
        self.all_patt_dict = {}
        self.detections = {}
        self.wv_calib = {}
        self.bad_slits = np.array([], dtype=int)
        # Reidentify each slit, and perform a fit
        for slit in range(self.nslits):
            # ToDO should we still be populating wave_calib with an empty dict here?
            if slit not in self.ok_mask:
                self.wv_calib[str(slit)] = None
                continue
            msgs.info('Reidentifying and fitting slit # {0:d}/{1:d}'.format(slit,self.nslits-1))
            # If this is a fixed format echelle, arxiv has exactly the same orders as the data and so
            # we only pass in the relevant arxiv spectrum to make this much faster
            ind_sp = arxiv_orders.index(orders[slit]) if ech_fixed_format else ind_arxiv
            if ech_fixed_format:
                msgs.info(f'Order: {orders[slit]}')
            sigdetect = wvutils.parse_param(self.par, 'sigdetect', slit)
            cc_thresh = wvutils.parse_param(self.par, 'cc_thresh', slit)
            rms_threshold = wvutils.parse_param(self.par, 'rms_threshold', slit)
            msgs.info("Using sigdetect =  {}".format(sigdetect))
            msgs.info("Using rms_threshold =  {}".format(rms_threshold))
            # get FWHM for this slit
            fwhm = set_fwhm(self.par, measured_fwhm=measured_fwhms[slit])
            self.detections[str(slit)], self.spec_cont_sub[:,slit], self.all_patt_dict[str(slit)] = \
                reidentify(self.spec[:,slit], self.spec_arxiv[:,ind_sp], self.wave_soln_arxiv[:,ind_sp],
                           self.tot_line_list, self.nreid_min, cc_thresh=cc_thresh, match_toler=self.match_toler,
                           cc_local_thresh=self.cc_local_thresh, nlocal_cc=self.nlocal_cc, nonlinear_counts=self.nonlinear_counts,
                           sigdetect=sigdetect, fwhm=fwhm, debug_peaks=self.debug_peaks, debug_xcorr=self.debug_xcorr,
                           debug_reid=self.debug_reid)
            # Check if an acceptable reidentification solution was found
            if not self.all_patt_dict[str(slit)]['acceptable']:
                self.wv_calib[str(slit)] = None
                self.bad_slits = np.append(self.bad_slits, slit)
                continue

            # Perform the fit
            n_final = wvutils.parse_param(self.par, 'n_final', slit)
            final_fit = wv_fitting.fit_slit(self.spec_cont_sub[:, slit], self.all_patt_dict[str(slit)],
                                         self.detections[str(slit)],
                                         self.tot_line_list, match_toler=self.match_toler,func=self.func, n_first=self.n_first,
                                         sigrej_first=self.sigrej_first, n_final=n_final,sigrej_final=self.sigrej_final)

            # Did the fit succeed?
            if final_fit is None:
                # This pattern wasn't good enough
                self.wv_calib[str(slit)] = None
                self.bad_slits = np.append(self.bad_slits, slit)
                continue
            # Is the RMS below the threshold?
            if final_fit['rms'] > rms_threshold:
                msgs.warn('---------------------------------------------------' + msgs.newline() +
                          'Reidentify report for slit {0:d}/{1:d}:'.format(slit, self.nslits-1) + msgs.newline() +
                          '  Poor RMS ({0:.3f})! Need to add additional spectra to arxiv to improve fits'.format(
                              final_fit['rms']) + msgs.newline() +
                          '---------------------------------------------------')
                self.bad_slits = np.append(self.bad_slits, slit)
                # Note this result in new_bad_slits, but store the solution since this might be the best possible

            # Add the patt_dict and wv_calib to the output dicts
            self.wv_calib[str(slit)] = copy.deepcopy(final_fit)
            if self.debug_fits:
                arc_fit_qa(self.wv_calib[str(slit)], title='Slit: {}'.format(str(slit)))

        # Print the final report of all lines
        report_final(self.nslits, self.all_patt_dict, self.detections, self.wv_calib, self.ok_mask, self.bad_slits)


    def get_results(self):
        return copy.deepcopy(self.all_patt_dict), copy.deepcopy(self.wv_calib)





class HolyGrail:
    """ General algorithm to wavelength calibrate spectroscopic data

    Parameters
    ----------
    spec : ndarray
        2D array of arcline spectra (nspec,nslit)
    lamps : :obj:`list`
        List of arc lamps to be used for wavelength calibration.
        E.g., ['ArI','NeI','KrI','XeI']
    par : ParSet or dict, default = default parset, optional
        This is the parset par['calibrations']['wavelengths']. A
        dictionary with the corresponding parameter names also works.
    ok_mask : ndarray, optional
        Array of good slits
    islinelist : bool, optional
        Is lines a linelist (True), or a list of ions (False)
    outroot : str, optional
        Name of output file
    debug : bool, optional
        Used to debug the algorithm
    verbose : bool, optional
        If True, the final fit will print out more detail as the RMS is
        refined, and lines are rejected. This is mostly helpful for
        developing the algorithm.
    binw : ndarray, optional
        Set the wavelength grid when identifying the best solution
    bind : ndarray, optional
        Set the dispersion grid when identifying the best solution
    nstore : int, optional
        The number of "best" initial solutions to consider
    use_unknowns : bool, optional
        If True, arc lines that are known to be present in the spectra,
        but have not been attributed to an element+ion, will be included
        in the fit.
    spectrograph : str, optional
        Spectrograph name

    Returns
    -------
    all_patt_dict : list of dicts
        A list of dictionaries, which contain the results from the
        preliminary pattern matching algorithm providing the first guess
        at the ID lines
    all_final_fit : list of dicts
        A list of dictionaries, which contain the full fitting results
        and final best guess of the line IDs

    """

    def __init__(self, spec, lamps, par=None, ok_mask=None, islinelist=False,
                 outroot=None, debug=False, verbose=False,
                 binw=None, bind=None, nstore=1, use_unknowns=True, 
                 nonlinear_counts=None, spectrograph=None):

        # Set some default parameters
        self._spec = spec
        self._par = pypeitpar.WavelengthSolutionPar() if par is None else par
        self._lamps = lamps
        self._npix, self._nslit = spec.shape
        self._nstore = nstore
        self._binw = binw
        self._bind = bind

        # Mask info
        if ok_mask is None:
            self._ok_mask = np.arange(self._nslit)
        else:
            self._ok_mask = ok_mask
        self._bad_slits = []  # List of bad slits

        # Set the input parameters
        self._nonlinear_counts = nonlinear_counts
        # JFH I'm not convinced that the codea actually does anything except use the lowest nsig, but am not sure
        self._match_toler = self._par['match_toler']
        self._func = self._par['func']
        self._n_first= self._par['n_first']
        self._sigrej_first= self._par['sigrej_first']
        self._sigrej_final= self._par['sigrej_final']

        self._use_unknowns = use_unknowns
        self._islinelist = islinelist

        self._outroot = outroot

        self._debug = debug
        self._verbose = verbose

        # Load the linelist to be used for pattern matching
        if self._islinelist:
            self._line_lists = self._lamps
            self._unknwns = self._lamps[:0].copy()
        else:
            if 'ThAr' in self._lamps:
                line_lists_all = waveio.load_line_lists(self._lamps)
                self._line_lists = line_lists_all[np.where(line_lists_all['ion'] != 'UNKNWN')]
                self._unknwns = line_lists_all[np.where(line_lists_all['ion'] == 'UNKNWN')]
            else:
                restrict = spectrograph if self._par['use_instr_flag'] else None
                self._line_lists = waveio.load_line_lists(
                    self._lamps, restrict_on_instr=restrict)
                self._unknwns = waveio.load_unknown_list(self._lamps)

        if self._use_unknowns:
            self._tot_list = astropy.table.vstack([self._line_lists, self._unknwns])
        else:
            self._tot_list = self._line_lists

        # Generate the final linelist and sort
        self._wvdata = np.array(self._tot_list['wave'].data)  # Removes mask if any
        self._wvdata.sort()

        # Find the wavelength solution!
        # KD Tree algorithm only works for ThAr - check first that this is what is being used
        self._thar = False
        if 'ThAr' in self._lamps and len(self._lamps) == 1:
            self._thar = True
            # Set up the grids to be used for pattern matching
            self.set_grids(ngridw=5000, ngridd=1000)
            msgs.info("Using KD Tree pattern matching algorithm to wavelength calibrate")
            self.run_kdtree()
        else:
            # Set up the grids to be used for pattern matching
            self.set_grids()
            msgs.info("Using brute force pattern matching algorithm to wavelength calibrate")
            self.run_brute()

    def get_results(self):
        return copy.deepcopy(self._all_patt_dict), copy.deepcopy(self._all_final_fit)

    def set_grids(self, ngridw = 300, ngridd=3000): #ngridw = 200, ngridd=2000):
        # Set the wavelength grid
        if self._binw is None:
            # Ideally, you want binw to roughly sample the A/pix of the spectrograph
            self._ngridw = ngridw
            self._binw = np.linspace(np.min(self._wvdata), np.max(self._wvdata), self._ngridw)
        else:
            self._ngridw = self._binw.size
        # Set the dispersion grid
        if self._bind is None:
            self._ngridd = ngridd
            #self._bind = np.linspace(-3.0, 1.0, self._ngridd)
            # JFH I have no idea why this goes down as low as -3.0. 3000/10^(-3.0) would be R ~ 3e6. No spectrograph
            # has a dispersion that high. I'm changing this to be -1.5 which would be R ~ 100,000 at 3000A. In this regime
            # one would anyway use the ThAr routine. I'm rasing
            # the upper limit to be 2.0 to handle low-resolution data (i.e in the near-IR 2.5e4/100 = R ~ 250

            self._bind = np.linspace(-1.5, 2.0, self._ngridd)
        else:
            self._ngridd = self._bind.size
        return

    def run_brute_loop(self, slit, tcent_ecent, wavedata=None):
        # Set the parameter space that gets searched
        rng_poly = [3, 4]            # Range of algorithms to check (only trigons+tetragons are supported)
        rng_list = range(3, 6)       # Number of lines to search over for the linelist
        rng_detn = range(3, 6)       # Number of lines to search over for the detected lines
        rng_pixt = [1.0]             # Pixel tolerance
        idthresh = 0.5               # Criteria for early return (at least this fraction of lines must have
                                     # an ID on either side of the spectrum)

        rms_threshold = wvutils.parse_param(self._par, 'rms_threshold', slit)
        msgs.info("Using rms_threshold =  {}".format(rms_threshold))
        best_patt_dict, best_final_fit = None, None
        # Loop through parameter space
        for poly in rng_poly:
            for detsrch in rng_detn:
                for lstsrch in rng_list:
                    for pix_tol in rng_pixt:
                        # JFH Note that results_brute and solve_slit are running on the same set of detections. I think this is the way
                        # it should be.
                        psols, msols = self.results_brute(tcent_ecent,poly=poly, pix_tol=pix_tol,
                                                          detsrch=detsrch, lstsrch=lstsrch,wavedata=wavedata)
                        patt_dict, final_fit = self.solve_slit(slit, psols, msols,tcent_ecent)
                        if final_fit is None:
                            # This is not a good solution
                            continue
                        # Test if this solution is better than the currently favoured solution
                        if best_patt_dict is None:
                            # First time a fit is found
                            best_patt_dict, best_final_fit = copy.deepcopy(patt_dict), copy.deepcopy(final_fit)
                            continue
                        elif final_fit['rms'] < rms_threshold:
                            # Has a better fit been identified (i.e. more lines identified)?
                            if len(final_fit['pixel_fit']) > len(best_final_fit['pixel_fit']):
                                best_patt_dict, best_final_fit = copy.deepcopy(patt_dict), copy.deepcopy(final_fit)
                            # Decide if an early return is acceptable
                            nlft = np.sum(best_final_fit['tcent'] < best_final_fit['spec'].size/2.0)
                            nrgt = best_final_fit['tcent'].size-nlft
                            if np.sum(best_final_fit['pixel_fit'] < 0.5)/nlft > idthresh and\
                                np.sum(best_final_fit['pixel_fit'] >= 0.5) / nrgt > idthresh:
                                # At least half of the lines on either side of the spectrum have been identified
                                return best_patt_dict, best_final_fit

        return best_patt_dict, best_final_fit

    def run_brute(self, min_nlines=10):
        """Run through the parameter space and determine the best solution
        """

        # ToDo This code appears to use the weak lines for everything throughout
        self._all_patt_dict = {}
        self._all_final_fit = {}
        good_fit = np.zeros(self._nslit, dtype=bool)
        self._det_weak = {}
        self._det_stro = {}
        for slit in range(self._nslit):
            msgs.info("Working on slit: {}".format(slit))
            if slit not in self._ok_mask:
                self._all_final_fit[str(slit)] = None
                continue
            # TODO Pass in all the possible params for detect_lines to arc_lines_from_spec, and update the parset
            # Detect lines, and decide which tcent to use
            sigdetect = wvutils.parse_param(self._par, 'sigdetect', slit)
            msgs.info("Using sigdetect =  {}".format(sigdetect))
            self._all_tcent, self._all_ecent, self._cut_tcent, self._icut, _  =\
                wvutils.arc_lines_from_spec(self._spec[:, slit].copy(), sigdetect=sigdetect, nonlinear_counts = self._nonlinear_counts)
            self._all_tcent_weak, self._all_ecent_weak, self._cut_tcent_weak, self._icut_weak, _  =\
                wvutils.arc_lines_from_spec(self._spec[:, slit].copy(), sigdetect=sigdetect, nonlinear_counts = self._nonlinear_counts)

            # Were there enough lines?  This mainly deals with junk slits
            if self._all_tcent.size < min_nlines:
                msgs.warn("Not enough lines to identify in slit {0:d}!".format(slit))
                self._det_weak[str(slit)] = [None,None]
                self._det_stro[str(slit)] = [None,None]
                # Remove from ok mask
                oklist = self._ok_mask.tolist()
                oklist.pop(slit)
                self._ok_mask = np.array(oklist)
                self._all_final_fit[str(slit)] = None
                continue
            # Setup up the line detection dicts
            self._det_weak[str(slit)] = [self._all_tcent_weak[self._icut_weak].copy(),self._all_ecent_weak[self._icut_weak].copy()]
            self._det_stro[str(slit)] = [self._all_tcent[self._icut].copy(),self._all_ecent[self._icut].copy()]

            # Run brute force algorithm on the weak lines
            best_patt_dict, best_final_fit = self.run_brute_loop(slit, self._det_weak[str(slit)])

            # Print preliminary report
            good_fit[slit] = self.report_prelim(slit, best_patt_dict, best_final_fit)

        # Now that all slits have been inspected, cross match to generate a
        # list of all lines in every slit, and refit all spectra
        if self._nslit > 1:
            msgs.info('Checking wavelength solution by cross-correlating with all slits')

            msgs.info('Cross-correlation iteration #1')
            obad_slits = self.cross_match(good_fit, self._det_weak)
            cntr = 2
            while obad_slits.size > 0:
                msgs.info('Cross-correlation iteration #{:d}'.format(cntr))
                good_fit = np.ones(self._nslit, dtype=bool)
                good_fit[obad_slits] = False
                bad_slits = self.cross_match(good_fit,self._det_weak)
                if np.array_equal(bad_slits, obad_slits):
                    break
                obad_slits = bad_slits.copy()
                cntr += 1
                if cntr > 10:
                    msgs.warn("Breaking while loop before convergence. Check the wavelength solution!")
                    break

        # With these updates to the fits of each slit, determine the final fit.
        self.finalize_fit(self._det_weak)

        # Print the final report of all lines
        self.report_final()
        return

    def run_kdtree(self, polygon=4, detsrch=7, lstsrch=10, pixtol=5):
        """
        KD Tree algorithm to wavelength calibrate spectroscopic data.
        Currently, this is only designed for ThAr lamp spectra. See the
        'run_brute' function if you want to calibrate longslit spectra.

        Parameters
        ----------
        polygon : int
          Number of sides to the polygon used in pattern matching.  For example:

            - polygon=3  -->  trigon (two anchor lines and one floating line)
            - polygon=4  -->  tetragon (two anchor lines and two floating lines)
            - polygon=5  -->  pentagon (two anchor lines and three floating lines)

        detsrch : int
            Number of consecutive detected lines used to generate a
            pattern. For example, if detsrch is 4, then for a trigon,
            the following patterns will be generated (assuming line #1
            is the left anchor):

                - 1 2 3:  (in this case line #3 is the right anchor)
                - 1 2 4:  (in this case line #4 is the right anchor)
                - 1 3 4:  (in this case line #4 is the right anchor)

        lstsrch : int
            Number of consecutive lines in the linelist used to generate
            a pattern.  See example above for detsrch
        pixtol : float
            Tolerance used to find good patterns. An acceptable match if
            the closest distance to a pattern is < pixtol/npix, where
            npix is the number of pixels in the spectral direction.
            Ideally, this should depend on the pattern...

        """

        # Load the linelist KD Tree
        lsttree, lindex = waveio.load_tree(polygon=polygon, numsearch=lstsrch)

        # Set the search error to be 5 pixels
        err = pixtol / self._npix

        self._all_patt_dict = {}
        self._all_final_fit = {}
        good_fit = np.zeros(self._nslit, dtype=bool)
        self._det_weak = {}
        self._det_stro = {}
        for slit in range(self._nslit):
            if slit not in self._ok_mask:
                self._all_final_fit[str(slit)] = {}
                continue
            # Detect lines, and decide which tcent to use
            sigdetect = wvutils.parse_param(self._par, 'sigdetect', slit)
            self._all_tcent, self._all_ecent, self._cut_tcent, self._icut, _ =\
                wvutils.arc_lines_from_spec(self._spec[:, slit], sigdetect=sigdetect, nonlinear_counts = self._nonlinear_counts)
            self._all_tcent_weak, self._all_ecent_weak, self._cut_tcent_weak, self._icut_weak, _ =\
                wvutils.arc_lines_from_spec(self._spec[:, slit], sigdetect=sigdetect, nonlinear_counts = self._nonlinear_counts)
            if self._all_tcent.size == 0:
                msgs.warn("No lines to identify in slit {0:d}!".format(slit+ 1))
                continue

            # Save the detections
            self._det_weak[str(slit)] = [self._all_tcent_weak[self._icut_weak].copy(),self._all_ecent_weak[self._icut_weak].copy()]
            self._det_stro[str(slit)] = [self._all_tcent[self._icut].copy(),self._all_ecent[self._icut].copy()]

            use_tcentp, use_ecentp = self.get_use_tcent(1, self._det_weak[str(slit)])
            use_tcentm, use_ecentm = self.get_use_tcent(-1, self._det_weak[str(slit)])
            if use_tcentp.size < detsrch:
                if self._verbose:
                    msgs.info("Not enough lines to test this solution, will attempt another.")
                return None, None

            # Create a detlines KD Tree
            maxlinear = 0.5*self._npix
            if polygon == 3:
                msgs.info("Generating patterns for a trigon")
                patternp, indexp = kdtree_generator.trigon(use_tcentp, detsrch, maxlinear)
                patternm, indexm = kdtree_generator.trigon(use_tcentm, detsrch, maxlinear)
            elif polygon == 4:
                msgs.info("Generating patterns for a tetragon")
                patternp, indexp = kdtree_generator.tetragon(use_tcentp, detsrch, maxlinear)
                patternm, indexm = kdtree_generator.tetragon(use_tcentm, detsrch, maxlinear)
            elif polygon == 5:
                msgs.info("Generating patterns for a pentagon")
                patternp, indexp = kdtree_generator.pentagon(use_tcentp, detsrch, maxlinear)
                patternm, indexm = kdtree_generator.pentagon(use_tcentm, detsrch, maxlinear)
            elif polygon == 6:
                msgs.info("Generating patterns for a hexagon")
                patternp, indexp = kdtree_generator.hexagon(use_tcentp, detsrch, maxlinear)
                patternm, indexm = kdtree_generator.hexagon(use_tcentm, detsrch, maxlinear)
            else:
                msgs.warn("Patterns can only be generated with 3 <= polygon <= 6")
                return None

            dettreep = scipy.spatial.cKDTree(patternp, leafsize=30)
            dettreem = scipy.spatial.cKDTree(patternm, leafsize=30)

            # Query the detections tree
            msgs.info("Querying KD tree patterns (slit {0:d}/{1:d})".format(slit+1, self._nslit))
            resultp = dettreep.query_ball_tree(lsttree, r=err)
            resultm = dettreem.query_ball_tree(lsttree, r=err)

            msgs.info("Identifying wavelengths for each pattern")
            # First flatten the KD Tree query results so numba can handle the input array
            flatresp = [item for sublist in resultp for item in sublist]
            flatresm = [item for sublist in resultm for item in sublist]
            flatidxp = [ii for ii, sublist in enumerate(resultp) for item in sublist]
            flatidxm = [ii for ii, sublist in enumerate(resultm) for item in sublist]
            # Obtain the correlate and anti-correlate solutions
            psols = results_kdtree_nb(use_tcentp, self._wvdata, flatresp, flatidxp, indexp,
                                      lindex, indexp.shape[1], self._npix)
            msols = results_kdtree_nb(use_tcentm, self._wvdata, flatresm, flatidxm, indexm,
                                      lindex, indexm.shape[1], self._npix)

            msgs.info("Identifying the best solution")
            patt_dict, final_fit = self.solve_slit(slit, psols, msols, self._det_weak[str(slit)], nselw=1, nseld=2)

            # Print preliminary report
            good_fit[slit] = self.report_prelim(slit, patt_dict, final_fit)

        # Now that all slits have been inspected, attempt to find a better
        # solution for orders that were not fit well, by estimating the
        # wavelength coverage of that slit
        #self.cross_match_order(good_fit)

        # With the updates to the fits of each slit, determine the final fit, and save the QA
        self.finalize_fit()

        # Print the final report of all lines
        self.report_final()


    # TODO This routine should be replaced with a new version based on my reidentify code
    def cross_match(self, good_fit, detections):
        """Cross-correlate the spectra across all slits to ID all of the lines.

        Parameters
        ----------
        good_fit : ndarray (bool)
            Indicates which slits are deemed to be a good fit (although,
            sometimes a bad fit can be labelled as a good fit). To
            remedy this, the true good fits are determined in this
            routine.

        """
        # Steps:
        # Check that all of the "good" slits are indeed good
        # For all of the bad slits, cross-correlate against all of the good slits to label each line
        # For all newly labeled lines, create a patt_dict of these labeled lines
        # Perform a final fit on these lines

        #self._debug = True
        # First, sort spectra according to increasing central wavelength
        ngd = good_fit.sum()
        idx_gd = np.zeros(ngd, dtype=int)
        wvc_gd = np.zeros(ngd, dtype=float)
        dsp_gd = np.zeros(ngd, dtype=float)
        wvc_gd_jfh = np.zeros(ngd, dtype=float)
        dsp_gd_jfh = np.zeros(ngd, dtype=float)
        xrng = np.arange(self._npix)
        xnpixmin1 = float(self._npix-1)
        cntr = 0
        for slit in range(self._nslit):
            # Masked?
            if slit not in self._ok_mask:
                continue
            if good_fit[slit]:
                idx_gd[cntr] = slit
                # TODO JFH We can get rid of this and thus not need patt_dict
                wvc_gd[cntr] = self._all_patt_dict[str(slit)]["bwv"]
                dsp_gd[cntr] = self._all_patt_dict[str(slit)]["bdisp"]
                # JFH stuff
                wave_soln = self._all_final_fit[str(slit)].pypeitfit.eval(xrng/xnpixmin1)
                wvc_gd_jfh[cntr] = wave_soln[self._npix//2]
                dsp_gd_jfh[cntr]= np.median(wave_soln - np.roll(wave_soln,1))
                cntr += 1
        srt = np.argsort(wvc_gd_jfh)
        sort_idx = idx_gd[srt]
        sort_wvc = wvc_gd[srt]
        sort_dsp = dsp_gd[srt]
        sort_wvc_jfh = wvc_gd_jfh[srt]
        sort_dsp_jfh = dsp_gd_jfh[srt]

        # Cross correlate all good spectra with each other, in order of wavelength
        ncrco = ngd*(ngd-1)//2
        ccorr_val = np.zeros(ncrco)
        shift_val = np.zeros(ncrco)
        dwvc_val = np.zeros(ncrco)
        slit_ids = np.zeros((ncrco, 2), dtype=int)
        cntr = 0
        # JFH Consider adding something in here that takes advantage of the
        for gd in range(0, sort_idx.size-1):
            for gc in range(gd+1, sort_idx.size):
                #corr = scipy.signal.correlate(self._spec[:, sort_idx[gd]], self._spec[:, sort_idx[gc]], mode='same')
                #amax = np.argmax(corr)
                # dwvc_val[cntr] = (sort_wvc[gc]-sort_wvc[gd]) / (0.5*(sort_dsp[gc]+sort_dsp[gd])) - (amax - self._spec.shape[0] // 2)
                # JFH replaced with more robust xcorr
                shift_val[cntr], ccorr_val[cntr]= wvutils.xcorr_shift(self._spec[:, sort_idx[gd]],self._spec[:, sort_idx[gc]],
                                                                      percent_ceil=50.0)
                #dwvc_val[cntr] = (sort_wvc[gc]-sort_wvc[gd]) / (0.5*(sort_dsp[gc]+sort_dsp[gd])) - shift
                # JFH TESTING
                dwvc_val[cntr] = (sort_wvc_jfh[gc]-sort_wvc_jfh[gd]) / (0.5*(sort_dsp_jfh[gc]+sort_dsp_jfh[gd])) - shift_val[cntr]
                slit_ids[cntr, 0] = gd
                slit_ids[cntr, 1] = gc
                cntr += 1


        # TODO Replace this code below with code based on either sigma_clipped_stats or djs_reject
        # Identify the good slits as those for which the cross-correlation is consistent with the mad of all the slits.
        # Bad slits are then the outliers.
        sigrej = 3.0
        mad = 1.4826 * np.median(np.abs(dwvc_val))
        gdmsk = np.abs(dwvc_val) < sigrej * mad
        for ii in range(100):  # Limit to 100 iterations - this will likely never be reached...
            ogdmsk = gdmsk.copy()
            mad = 1.4826 * np.median(np.abs(dwvc_val[gdmsk]))
            gdmsk = np.abs(dwvc_val) < sigrej*mad
            if np.array_equal(gdmsk, ogdmsk):
                break

        if self._debug:
            # TODO Add something here indicating slit indices? Like plot the cc pairs next to the points?
            xplt = np.arange(dwvc_val.size)
            plt.plot(xplt[~gdmsk], dwvc_val[~gdmsk], 'rx',label ='bad slit')
            plt.plot(xplt[gdmsk], dwvc_val[gdmsk], 'bo',label = 'good slit')
            plt.hlines(0,xplt.min(),xplt.max(), color='black',linestyle='--')
            plt.xticks(xplt)
            plt.legend()
            plt.show()

        # Catalogue the good and bad slits.

        # ToDO Basically a slit needs to have a bad cross-correlation with every other slit in order
        # to be classified as a bad slit here. Is this the behavior we want?? Maybe we should be more
        # conservative and call a bad any slit which results in an outlier here?
        good_slits = np.sort(sort_idx[np.unique(slit_ids[gdmsk, :].flatten())])
        bad_slits = np.setdiff1d(np.arange(self._nslit), good_slits, assume_unique=True)
        nbad = bad_slits.size
        if nbad > 0:
            msgs.info('Working on {:d}'.format(nbad) + ' bad slits: {:}'.format(bad_slits + 1))

        # Get the sign (i.e. if pixels correlate/anticorrelate with wavelength)
        # and dispersion (A/pix). Assume these are the same for all slits

        # JFH Changed this to take the median which is more robust. Could even reject outliers
        disp_good = np.zeros(good_slits.size,dtype=float)
        sign_good = np.zeros(good_slits.size,dtype=int)
        wvc_good  = np.zeros(good_slits.size,dtype=float)
        for islit in range(good_slits.size):
            # JFH ToDO Could just use the good wavelength solutions and then we would not need this sign and hence all_patt_ict
            sign_good[islit] =  self._all_patt_dict[str(good_slits[islit])]['sign']
            # JFH stuff
            wave_soln = self._all_final_fit[str(good_slits[islit])].pypeitfit.eval(xrng / xnpixmin1)
            wvc_good[islit] = wave_soln[self._npix // 2]
            disp_good[islit] = np.median(wave_soln - np.roll(wave_soln, 1))


        disp_med = np.median(disp_good)
        sign = np.median(sign_good)
        #disp = self._all_patt_dict[str(good_slits[0])]['bdisp']
        #sign = self._all_patt_dict[str(good_slits[0])]['sign']

        # For all of the bad slits, estimate some line wavelengths
        new_bad_slits = np.array([], dtype=int)
        for bs in bad_slits:
            if bs not in self._ok_mask:
                continue
            if detections[str(bs)][0] is None:  # No detections at all; slit is hopeless
                msgs.warn('Slit {:d}'.format(bs) + ' has no arc line detections.  Likely this slit is junk!')
                self._bad_slits.append(bs)
                continue
            bsdet, _ = self.get_use_tcent(sign, detections[str(bs)])
            lindex = np.array([], dtype=int)
            dindex = np.array([], dtype=int)
            wcen = np.zeros(good_slits.size)
            disp = np.zeros(good_slits.size)
            shift_vec = np.zeros(good_slits.size)
            stretch_vec = np.zeros(good_slits.size)
            ccorr_vec = np.zeros(good_slits.size)
            for cntr, gs in enumerate(good_slits):
                msgs.info('Cross-correlating bad slit # {:d}'.format(bs + 1) + ' with good slit # {:d}'.format(gs + 1))
                # Match the peaks between the two spectra.
                # spec_gs_adj is the stretched spectrum
                success, shift_vec[cntr], stretch_vec[cntr], ccorr_vec[cntr], _, _ =  \
                    wvutils.xcorr_shift_stretch(self._spec[:, bs],self._spec[:, gs], debug = self._debug)
                if not success:
                    continue
                # ToDo Put in a cut on the cross-correlation value here in this logic so that we only consider slits that are sufficiently similar

                # Estimate wcen and disp for this bad slit based on its shift/stretch relative to the good slit
                disp[cntr] = disp_good[cntr]/stretch_vec[cntr]
                wcen[cntr] = wvc_good[cntr] - shift_vec[cntr]*disp[cntr]

                # For each peak in the gs spectrum, identify the corresponding peaks in the bs spectrum. Do this by
                # transforming these good slit line pixel locations into the (shifted and stretched) bs frame
                gsdet, _ = self.get_use_tcent(sign, detections[str(gs)])
                gsdet_ss = gsdet*stretch_vec[cntr] + shift_vec[cntr]
                if self._debug:
                    plt.figure(figsize=(14, 6))
                    xrng = np.arange(self._npix)
                    tampl_bs = np.interp(bsdet, xrng, self._spec[:, bs])
                    plt.plot(xrng, self._spec[:, bs], color='red', drawstyle='steps-mid', label='bad slit arc', linewidth=1.0, zorder= 10)
                    plt.plot(bsdet, tampl_bs, 'r.', markersize=10.0, label='bad slit lines', zorder= 10)
                    tampl_gs = np.interp(gsdet, xrng, self._spec[:, gs])
                    plt.plot(xrng, self._spec[:, gs], color='black', drawstyle='steps-mid', linestyle=':',
                             label='good slit arc', linewidth=0.5)
                    plt.plot(gsdet, tampl_gs, 'k+', markersize=8.0, label='good slit lines')
                    gdarc_ss = wvutils.shift_and_stretch(self._spec[:, gs], shift_vec[cntr], stretch_vec[cntr])
                    #tampl_ss = np.interp(gsdet_ss, xrng, gdarc_ss)
                    for iline in range(gsdet_ss.size):
                        plt.plot([gsdet[iline],gsdet_ss[iline]],[tampl_gs[iline], tampl_gs[iline]], color='cornflowerblue', linewidth = 1.0)
                    plt.plot(xrng, gdarc_ss, color='black', drawstyle='steps-mid', label='good slit arc shift/stretch', linewidth=1.0)
                    plt.plot(gsdet_ss, tampl_gs, 'k.', markersize=10.0, label='predicted good slit lines')
                    plt.title('Cross-correlation of bad slit # {:d}'.format(bs+1) + ' and good slit # {:d}'.format(gs+1) +
                              ': ccor = {:5.3f}'.format(ccorr_vec[cntr]) +
                              ', shift = {:6.1f}'.format(shift_vec[cntr]) +
                              ', stretch = {:5.4f}'.format(stretch_vec[cntr]) +
                              ', wv_cen = {:7.1f}'.format(wcen[cntr]) +
                              ', disp = {:5.3f}'.format(disp[cntr]))
                    plt.ylim(-5.0, 1.5*self._spec[:, bs].max())
                    plt.legend()
                    plt.show()

                # Calculate wavelengths for all of the gsdet detections
                wvval = self._all_final_fit[str(gs)].pypeitfit.eval(xrng / xnpixmin1)
                # Loop over the bad slit line pixel detections and find the nearest good slit line
                for dd in range(bsdet.size):
                    pdiff = np.abs(bsdet[dd]-gsdet_ss)
                    bstpx = np.argmin(pdiff)
                    # If a match is found within 2 pixels, consider this a successful match
                    if pdiff[bstpx] < 2.0:
                        # Using the good slit wavelength solution, search for the nearest line in the line list
                        bstwv = np.abs(self._wvdata - wvval[bstpx])
                        # This is probably not a good match
                        if bstwv[np.argmin(bstwv)] > 2.0*disp_med:
                            continue
                        lindex = np.append(lindex, np.argmin(bstwv)) # index in the line list self._wvdata
                        dindex = np.append(dindex, dd)               # index in the array of pixel detections bsdet
            # Finalize the best guess of each line
            # Initialise the patterns dictionary
            patt_dict = dict(acceptable=False, nmatch=0, ibest=-1, bwv=0.,
                             sigdetect= wvutils.parse_param(self._par, 'sigdetect', bs),
                             mask=np.zeros(bsdet.size, dtype=bool), scores = None)
            patt_dict['sign'] = sign
            patt_dict['bwv'] = np.median(wcen[wcen != 0.0])
            patt_dict['bdisp'] = np.median(disp[disp != 0.0])
            patterns.solve_triangles(bsdet, self._wvdata, dindex, lindex, patt_dict = patt_dict)

            if self._debug:
                tmp_list = astropy.table.vstack([self._line_lists, self._unknwns])
                match_qa(self._spec[:, bs], bsdet, tmp_list,patt_dict['IDs'], patt_dict['scores'])

            # Use only the perfect IDs
            iperfect = np.array(patt_dict['scores']) != 'Perfect'
            patt_dict['mask'][iperfect] = False
            patt_dict['nmatch'] = np.sum(patt_dict['mask'])
            if patt_dict['nmatch'] < 3:
                patt_dict['acceptable'] = False

            # Check if a solution was found
            if not patt_dict['acceptable']:
                new_bad_slits = np.append(new_bad_slits, bs)
                continue
            final_fit = wv_fitting.fit_slit(self._spec[:, bs], patt_dict, bsdet, self._line_lists)
            #final_fit = self.fit_slit(bs, patt_dict, bsdet)
            if final_fit is None:
                # This pattern wasn't good enough
                new_bad_slits = np.append(new_bad_slits, bs)
                continue

            if final_fit['rms'] >  wvutils.parse_param(self._par, 'rms_threshold', bs):
                msgs.warn('---------------------------------------------------' + msgs.newline() +
                          'Cross-match report for slit {0:d}/{1:d}:'.format(bs + 1, self._nslit) + msgs.newline() +
                          '  Poor RMS ({0:.3f})! Will try cross matching iteratively'.format(final_fit['rms']) + msgs.newline() +
                          '---------------------------------------------------')
                # Store this result in new_bad_slits, so the iteration can be performed,
                # but make sure to store the result, as this might be the best possible.
                new_bad_slits = np.append(new_bad_slits, bs)
            self._all_patt_dict[str(bs)] = copy.deepcopy(patt_dict)
            self._all_final_fit[str(bs)] = copy.deepcopy(final_fit)
            if self._debug:
                xplt = np.linspace(0.0, 1.0, self._npix)
                yplt = utils.func_val(final_fit['fitc'], xplt, 'legendre', minx=0.0, maxx=1.0)
                plt.plot(final_fit['pixel_fit'], final_fit['wave_fit'], 'bx')
                plt.plot(xrng, yplt, 'r-')
                plt.show()
        return new_bad_slits

    def cross_match_order(self, good_fit):
        """Using the solutions of all orders, identify the good solutions, and refit the bad ones!

        TODO: This function needs work... The first few lines of code successfully pick up the good orders,
        but we need a new routine that (based on an estimated central wavelength and dispersion) can successfully
        ID all of the lines.
        """
        # DEPRECATED (NOT USED)

        # First determine the central wavelength and dispersion of every slit, using the known good solutions
        xplt = np.arange(self._nslit)
        yplt, dplt = np.zeros(self._nslit), np.zeros(self._nslit)
        imsk = np.ones(self._nslit, dtype=int)
        for slit in range(self._nslit):
            if good_fit[slit]:
                yplt[slit] = self._all_patt_dict[str(slit)]['bwv']
                dplt[slit] = self._all_patt_dict[str(slit)]['bdisp']
                imsk[slit] = 0

        mask, fit = utils.robust_polyfit(xplt, yplt, 2, function='polynomial', sigma=2,
                                         initialmask=imsk, forceimask=True)
        good_fit[mask == 1] = False
        wavemodel = utils.func_val(fit, xplt, 'polynomial')
        disp = np.median(dplt[good_fit])

        # TODO: maybe rethink the model at this point? Using the derived
        # central wavelength and dispersion identify liens in all orders?

        if self._debug:
            plt.subplot(211)
            plt.plot(xplt, wavemodel, 'r-')
            ww = np.where(mask==0)
            plt.plot(xplt[ww], yplt[ww], 'bx')
            ww = np.where(mask==1)
            plt.plot(xplt[ww], yplt[ww], 'rx')
            plt.subplot(212)
            plt.plot(xplt, dplt, 'bx')
            plt.show()
            #embed()

        fact_nl = 1.2  # Non linear factor
        new_good_fit = np.zeros(self._nslit, dtype=bool)
        for slit in range(self._nslit):
            wmin = wavemodel[slit] - fact_nl*disp*self._npix/2
            wmax = wavemodel[slit] + fact_nl*disp*self._npix/2
            ww = np.where((self._wvdata > wmin) & (self._wvdata < wmax))
            wavedata = self._wvdata[ww]
            msgs.info('Brute force ID for slit {0:d}/{1:d}'.format(slit+1, self._nslit))
            best_patt_dict, best_final_fit =\
                self.run_brute_loop(slit, arrerr=self._det_weak[str(slit)], wavedata=wavedata)

            self._all_patt_dict[str(slit)] = copy.deepcopy(best_patt_dict)
            self._all_final_fit[str(slit)] = copy.deepcopy(best_final_fit)
            new_good_fit[slit] = self.report_prelim(slit, best_patt_dict, best_final_fit)
        return new_good_fit


        # Set some fitting parameters
        if self._n_final is None:
            order = 4
        else:
            order = self._n_final

        ofit = [5, 3, 1, 0]
        lnpc = len(ofit) - 1

        # Prepare the fitting coefficients
        xv = np.arange(self._npix)/(self._npix-1)
        ords = np.arange(self._nslit)
        xcen = xv[:, np.newaxis].repeat(self._nslit, axis=1)
        extrapord = ~good_fit
        maskord = np.where(extrapord)[0]

        coeffs = None
        waves = np.zeros(xcen.shape, dtype=float)
        for slit in range(self._nslit):
            if good_fit[slit]:
                func = self._all_final_fit[str(slit)]['function']
                fmin = self._all_final_fit[str(slit)]['fmin']
                fmax = self._all_final_fit[str(slit)]['fmax']
                fitc = self._all_final_fit[str(slit)]['fitc']
                if coeffs is None:
                    coeffs = np.zeros((fitc.size, self._nslit))
                coeffs[:, slit] = fitc.copy()
                waves[:, slit] = utils.func_val(fitc, xv, func, minx=fmin, maxx=fmax)

        msgs.info("Performing a PCA on the order wavelength solutions")
        #embed()
        pca_wave, outpar = pca.basis(xcen, waves, coeffs, lnpc, ofit, x0in=ords, mask=maskord, skipx0=False, function=func)

        # Report the QA
        # TODO: fix setup passing
        setup = "BLAH"
        pca.pca_plot(setup, outpar, ofit, "wave_cross_match", pcadesc="Wavelength calibration PCA")


        # Extrapolate the remaining orders requested
        #extrap_wave, outpar = pca.extrapolate(outpar, ords)

        # Determine if pixels correlate and anti-correlate with wavelength
        signs = np.zeros(self._nslit, dtype=int)
        for slit in range(self._nslit):
            wvval = pca_wave[:, slit]
            if wvval[wvval.size//2] > wvval[wvval.size//2-1]:
                signs[slit] = 1
            else:
                signs[slit] = -1
        sign = 1
        if np.sum(signs) < 0:
            sign = -1

        new_bad_slits = np.array([], dtype=int)
        # Using the first guesses at the wavelength solution, identify lines
        for slit in range(self._nslit):
            # Get the detections
            dets, _ = self.get_use_tcent(sign, self._det_weak[str(slit)])
            lindex = np.array([], dtype=int)
            dindex = np.array([], dtype=int)
            # Calculate wavelengths for the gsdet detections
            wvval = pca_wave[:, slit]
            wvcen = wvval[wvval.size//2]
            disp = abs(wvval[wvval.size//2] - wvval[wvval.size//2-1])
            for dd in range(dets.size):
                pdiff = np.abs(dets[dd] - xv)
                bstpx = np.argmin(pdiff)
                bstwv = np.abs(self._wvdata - wvval[bstpx])
                if bstwv[np.argmin(bstwv)] > 10.0 * disp:
                    # This is probably not a good match
                    continue
                lindex = np.append(lindex, np.argmin(bstwv))
                dindex = np.append(dindex, dd)

            # Finalize the best guess of each line
            # Initialise the patterns dictionary
            patt_dict = dict(acceptable=False, nmatch=0, ibest=-1, bwv=0.,
                             sigdetect=wvutils.parse_param(self._par, 'sigdetect', slit),
                             mask=np.zeros(dets.size, dtype=bool))
            patt_dict['sign'] = sign
            patt_dict['bwv'] = wvcen
            patt_dict['bdisp'] = disp

            patterns.solve_triangles(dets, self._wvdata, dindex, lindex, patt_dict)
            # Check if a solution was found
            if not patt_dict['acceptable']:
                new_bad_slits = np.append(new_bad_slits, slit)
                msgs.warn('---------------------------------------------------' + msgs.newline() +
                          'Cross-match report for slit {0:d}/{1:d}:'.format(slit, self._nslit-1) + msgs.newline() +
                          '  Lines could not be identified! Will try cross matching iteratively' + msgs.newline() +
                          '---------------------------------------------------')
                continue
            final_fit = self.fit_slit(slit, patt_dict, dets)
            if final_fit is None:
                # This pattern wasn't good enough
                new_bad_slits = np.append(new_bad_slits, slit)
                msgs.warn('---------------------------------------------------' + msgs.newline() +
                          'Cross-match report for slit {0:d}/{1:d}:'.format(slit, self._nslit-1) + msgs.newline() +
                          '  Fit was not good enough! Will try cross matching iteratively' + msgs.newline() +
                          '---------------------------------------------------')
                continue
            if final_fit['rms'] > wvutils.parse_param(self._par, 'rms_threshold', slit):
                msgs.warn('---------------------------------------------------' + msgs.newline() +
                          'Cross-match report for slit {0:d}/{1:d}:'.format(slit, self._nslit-1) + msgs.newline() +
                          '  Poor RMS ({0:.3f})! Will try cross matching iteratively'.format(final_fit['rms']) + msgs.newline() +
                          '---------------------------------------------------')
                # Store this result in new_bad_slits, so the iteration can be performed,
                # but make sure to store the result, as this might be the best possible.
                new_bad_slits = np.append(new_bad_slits, slit)
            self._all_patt_dict[str(slit)] = copy.deepcopy(patt_dict)
            self._all_final_fit[str(slit)] = copy.deepcopy(final_fit)
            if self._debug:
                xplt = np.linspace(0.0, 1.0, self._npix)
                yplt = utils.func_val(final_fit['fitc'], xplt, 'legendre', minx=0.0, maxx=1.0)
                plt.plot(final_fit['pixel_fit'], final_fit['wave_fit'], 'bx')
                plt.plot(xplt, yplt, 'r-')
                plt.show()
                #embed()

        # debugging
        if self._debug:
            # First determine the central wavelength and dispersion of every slit, using the known good solutions
            xplt = np.arange(self._nslit)
            yplt, dplt = np.zeros(self._nslit), np.zeros(self._nslit)
            imsk = np.ones(self._nslit, dtype=int)
            for slit in range(self._nslit):
                if good_fit[slit]:
                    yplt[slit] = self._all_patt_dict[str(slit)]['bwv']
                    dplt[slit] = self._all_patt_dict[str(slit)]['bdisp']
                    imsk[slit] = 0

            mask, fit = utils.robust_polyfit(xplt, yplt, 2, function='polynomial', sigma=2,
                                             initialmask=imsk, forceimask=True)

            ymodel = utils.func_val(fit, xplt, 'polynomial')
            plt.subplot(211)
            plt.plot(xplt, ymodel, 'r-')
            ww = np.where(mask==0)
            plt.plot(xplt[ww], yplt[ww], 'bx')
            ww = np.where(mask==1)
            plt.plot(xplt[ww], yplt[ww], 'rx')
            plt.subplot(212)
            plt.plot(xplt, dplt, 'bx')
            plt.show()
            #embed()

        return new_bad_slits

    def get_use_tcent_old(self, corr, cut=True, arr_err=None, weak=False):
        """
        Grab the lines to use

        Args:
            corr:  int
                Set if pixels correlate with wavelength (corr==1) or
                anticorrelate (corr=-1)
            arr_err:
                A list [tcent, ecent] indicating which detection list
                should be used. Note that if arr_err is set then the
                weak keyword is ignored.
            weak: bool, optional
                If True, return the weak lines
            cut: bool, optional
                Cut on the lines according to significance

        Returns:
            tuple: arr, err

        """
        # Decide which array to use
        if arr_err is None:
            if weak:
                if cut:
                    arr = self._all_tcent_weak.copy()[self._icut_weak]
                    err = self._all_ecent_weak.copy()[self._icut_weak]
                else:
                    msgs.error('CODING ERROR: Cut must be True')
            else:
                if cut:
                    arr = self._all_tcent.copy()[self._icut]
                    err = self._all_ecent.copy()[self._icut]
                else:
                    msgs.error('CODING ERROR: Cut must be True')
        else:
            arr, err = arr_err[0], arr_err[1]
        # Return the appropriate tcent
        if corr == 1:
            return arr, err
        else:
            return (self._npix - 1.0) - arr[::-1], err[::-1]

    def get_use_tcent(self, corr, tcent_ecent):
        """
        Grab the lines to use

        Args:
            corr:  int
              Set if pixels correlate with wavelength (corr==1) or
              anticorrelate (corr=-1)
            tcent_ecent:
              A list [tcent, ecent] indicating which detection list
              should be used. Note that if arr_err is set then the weak
              keyword is ignored.

        Returns:
            tuple: arr, err
        """
        # Return the appropriate tcent
        tcent, ecent = tcent_ecent[0], tcent_ecent[1]
        if corr == 1:
            return tcent, ecent
        else:
            return (self._npix - 1.0) - tcent[::-1], ecent[::-1]


    def results_brute(self, tcent_ecent, poly=3, pix_tol=0.5, detsrch=5, lstsrch=5, wavedata=None):
        """
        Need some docs here. I think this routine generates the
        patterns, either triangles are quadrangles.

        Parameters
        ----------
        tcent_ecent: list of ndarrays, [tcent, ecent]
        poly, optional:
            algorithms to use for pattern matching. Only triangles (3)
            and quadrangles (4) are supported
        pix_tol, optional:
            tolerance that is used to determine if a pattern match is
            successful (in units of pixels)
        detsrch, optional:
            Number of lines to search over for the detected lines
        lstsrch, optional:
            Number of lines to search over for the detected lines
        wavedata, optional:
        arrerr, optional:

        """
        # Import the pattern matching algorithms
        if poly == 3:
            from pypeit.core.wavecal.patterns import triangles as generate_patterns
        elif poly == 4:
            from pypeit.core.wavecal.patterns import quadrangles as generate_patterns
        else:
            msgs.warn("Pattern matching is only available for trigons and tetragons.")
            return None, None

        if wavedata is None:
            wavedata = self._wvdata

        # Test if there are enough lines to generate a solution
        use_tcent, _ = self.get_use_tcent(1, tcent_ecent)
        if use_tcent.size < lstsrch or use_tcent.size < detsrch:
            if self._verbose:
                msgs.info("Not enough lines to test this solution, will attempt another.")
            return None, None

        if self._verbose:
            msgs.info("Begin pattern matching")

        # First run pattern recognition assuming pixels correlate with wavelength
        dindexp, lindexp, wvcenp, dispsp = generate_patterns(use_tcent, wavedata, self._npix,
                                                             detsrch, lstsrch, pix_tol)
        # Now run pattern recognition assuming pixels anti-correlate with wavelength
        use_tcent, _ = self.get_use_tcent(-1, tcent_ecent)
        dindexm, lindexm, wvcenm, dispsm = generate_patterns(use_tcent, wavedata, self._npix,
                                                             detsrch, lstsrch, pix_tol)
        return (dindexp, lindexp, wvcenp, dispsp,), (dindexm, lindexm, wvcenm, dispsm,)

    def results_kdtree(self, use_tcent, res, dindex, lindex, ordfit=2):
        # Assign wavelengths to each pixel
        nrows = len(res)
        ncols = sum(map(len, res))
        nindx = dindex.shape[1]
        wvdisp = np.zeros(ncols)
        wvcent = np.zeros(ncols)
        dind = np.zeros((ncols, nindx), dtype=int)
        lind = np.zeros((ncols, nindx), dtype=int)
        cnt = 0
        for x in range(nrows):
            for y in range(len(res[x])):
                # dx = use_tcent[dindex[x, -1]] - use_tcent[dindex[x, 0]]
                # dp = self._wvdata[lindex[res[x][y], -1]] - self._wvdata[lindex[res[x][y], 0]]
                # try:
                #     null, cgrad = utils.robust_polyfit(use_tcent[dindex[x, :]], self._wvdata[lindex[res[x][y], :]],
                #                                        1, sigma=2.0, verbose=False)
                #     wvdisp[cnt] = cgrad[1]
                # except:
                #     wvdisp[cnt] = (dp / dx)
                #
                coeff = np.polyfit(use_tcent[dindex[x, :]], self._wvdata[lindex[res[x][y]]], ordfit)
                wvcent[cnt] = np.polyval(coeff, self._npix / 2.0)
                wvdisp[cnt] = abs(np.polyval(coeff, (self._npix+1) / 2.0) - wvcent[cnt])
                dind[cnt, :] = dindex[x, :]
                lind[cnt, :] = lindex[res[x][y], :]
                cnt += 1
        return dind, lind, wvcent, wvdisp

    def solve_slit(self, slit, psols, msols, tcent_ecent, nstore=1, nselw=3, nseld=3):
        """
        Need some docs here. I think this routine creates a 2d histogram
        of the patterns and searches for the most represented wave_cen
        and log10(disp). Then it attempts to fit each value determined
        (default of 1) to try to figure out if it is a reasonable fit.

        Args:
            slit:
            psols:
            msols:
            tcent_ecent: list, [tcent, ecent]
            nstore:
                Number of pattern matches to store and fit
            nselw:
                All solutions around the best central wavelength
                solution within +- nselw are selected to be fit
            nseld:
                All solutions around the best log10(dispersion) solution
                within +- nseld are selected to be fit

        Returns:
            tuple: patt_dict, final_dict
        """

        # Extract the solutions
        dindexp, lindexp, wvcenp, dispsp = psols
        dindexm, lindexm, wvcenm, dispsm = msols

        # Remove any invalid results from correlate
        ww = np.where((self._binw[0] < wvcenp) & (wvcenp < self._binw[-1]) &
                      (10.0 ** self._bind[0] < dispsp) & (dispsp < 10.0 ** self._bind[-1]))
        dindexp = dindexp[ww[0], :]
        lindexp = lindexp[ww[0], :]
        dispsp = dispsp[ww]
        wvcenp = wvcenp[ww]

        # Remove any invalid results from anticorrelate
        ww = np.where((self._binw[0] < wvcenm) & (wvcenm < self._binw[-1]) &
                      (10.0 ** self._bind[0] < dispsm) & (dispsm < 10.0 ** self._bind[-1]))
        dindexm = dindexm[ww[0], :]
        lindexm = lindexm[ww[0], :]
        dispsm = dispsm[ww]
        wvcenm = wvcenm[ww]

        # Construct the histograms
        histimgp, xed, yed = np.histogram2d(wvcenp, np.log10(dispsp), bins=[self._binw, self._bind])
        histimgm, xed, yed = np.histogram2d(wvcenm, np.log10(dispsm), bins=[self._binw, self._bind])
        #histimgp = gaussian_filter(histimgp, 3)
        #histimgm = gaussian_filter(histimgm, 3)
        histimg = histimgp - histimgm
        sm_histimg = scipy.ndimage.gaussian_filter(histimg, [30, 15])

        #histpeaks = patterns.detect_2Dpeaks(np.abs(sm_histimg))
        histpeaks = patterns.detect_2Dpeaks(np.abs(histimg))

        # Find the indices of the nstore largest peaks
        bidx = np.unravel_index(np.argpartition(np.abs(histpeaks*histimg), -nstore, axis=None)[-nstore:], histimg.shape)

        # Get the peak value of central wavelength and dispersion
        allwcen = self._binw[bidx[0]]
        alldisp = self._bind[bidx[1]]
        allhnum = np.abs(histimg[bidx])

        #debug = False
        if self._debug:# or slit==2:
            this_hist = histimg
            plt.clf()
            rect_image = [0.12, 0.05, 0.85, 0.9]
            fx = plt.figure(1, figsize=(12, 8))
            ax_image = fx.add_axes(rect_image)
            extent = [self._binw[0], self._binw[-1], self._bind[0], self._bind[-1]]
            # plt.subplot(221)
            # plt.imshow((np.abs(histimg[:, ::-1].T)), extent=extent, aspect='auto')
            # plt.subplot(222)
            # plt.imshow((np.abs(sm_histimg[:, ::-1].T)), extent=extent, aspect='auto')
            # plt.subplot(223)
            # plt.imshow((np.abs(histimgp[:, ::-1].T)), extent=extent, aspect='auto')
            # plt.subplot(224)
            # plt.imshow((np.abs(histimgm[:, ::-1].T)), extent=extent, aspect='auto')
            #plt.imshow((np.abs(sm_histimg[:, ::-1].T)), extent=extent, aspect='auto')
            cimg = ax_image.imshow(this_hist.T, extent=extent, aspect='auto',vmin=-2.0,vmax=5.0,
                       interpolation='nearest',origin='lower',cmap='Set1')
            nm = histimg.max() - histimg.min()
            ticks = np.arange(this_hist.min(),this_hist.max() + 1,1)
            cbar = fx.colorbar(cimg, ax=ax_image,ticks = ticks,drawedges = True, extend ='both',
                               spacing = 'proportional',orientation ='horizontal')
            cbar.set_ticklabels(ticks)
            cbar.set_label('# of Occurences')
            ax_image.set_xlabel('Central Wavelength (Angstroms)')
            ax_image.set_ylabel('log10(Dispersion/(Ang/pixel))')
            delta_wv = np.median(self._binw - np.roll(self._binw,1))
            delta_disp = np.median(self._bind - np.roll(self._bind,1))
            label = 'Maximum: (lam, log10(disp)) = ({:8.2f}'.format(allwcen[0]) + ',{:5.3f})'.format(alldisp[0])
            ax_image.plot(allwcen + delta_wv/2.0, alldisp + delta_disp/2.0, color='red', marker='+', markersize=10.0, fillstyle='none',
                     linestyle='None', zorder = 10,label=label)
            ax_image.legend()
            #plt.plot(self._binw[bidx[0]], self._bind[bidx[1]], 'r+')
#            ax_image.axvline(allwcen, color='r', linestyle='--')
#            ax_image.axhline(alldisp, color='r', linestyle='--')
            plt.show()
            #if False:
            #    plt.clf()
            #    plt.imshow(histimgp[:, ::-1].T, extent=extent, aspect='auto',vmin=0.0,vmax=1.0)
            #    plt.show()


        # Find all good solutions
        bestlist = []
        for idx in range(nstore):
            # Select all solutions around the best solution within a square of side 2*nsel
            wlo = self._binw[max(0, bidx[0][idx] - nselw)]
            whi = self._binw[min(self._ngridw - 1, bidx[0][idx] + nselw)]
            dlo = 10.0 ** self._bind[max(0, bidx[1][idx] - nseld)]
            dhi = 10.0 ** self._bind[min(self._ngridd - 1, bidx[1][idx] + nseld)]
            if histimgp[bidx][idx] > histimgm[bidx][idx]:
                wgd = np.where((wvcenp > wlo) & (wvcenp < whi) & (dispsp > dlo) & (dispsp < dhi))
                dindex = dindexp[wgd[0], :].flatten()
                lindex = lindexp[wgd[0], :].flatten()
                sign = +1
            else:
                wgd = np.where((wvcenm > wlo) & (wvcenm < whi) & (dispsm > dlo) & (dispsm < dhi))
                dindex = dindexm[wgd[0], :].flatten()
                lindex = lindexm[wgd[0], :].flatten()
                sign = -1
            # Store relevant values in an array to solve for best solution
            bestlist.append([allwcen[idx], alldisp[idx], allhnum[idx], sign, dindex, lindex])

        if self._verbose:
            msgs.info("Fitting the wavelength solution for each slit")
        patt_dict, final_dict = None, None
        for idx in range(nstore):
            # Solve the patterns
            tpatt_dict = self.solve_patterns(slit, bestlist[idx], tcent_ecent)
            if tpatt_dict is None:
                # This pattern wasn't good enough
                continue
            # Fit the full set of lines with the derived patterns
            use_tcent, _ = self.get_use_tcent(tpatt_dict['sign'], tcent_ecent)
            tfinal_dict = wv_fitting.fit_slit(self._spec[:, slit], tpatt_dict, use_tcent, self._line_lists)
            # tfinal_dict = self.fit_slit(slit, tpatt_dict, use_tcent)
            if tfinal_dict is None:
                # This pattern wasn't good enough
                continue
            # Check if this solution is better than the last
            if patt_dict is None:
                # First time a fit is found
                patt_dict, final_dict = tpatt_dict, tfinal_dict
                continue
            elif tfinal_dict['rms'] < wvutils.parse_param(self._par, 'rms_threshold', slit):
                # Has a better fit been identified (i.e. more lines ID)?
                if len(tfinal_dict['pixel_fit']) > len(final_dict['pixel_fit']):
                    patt_dict, final_dict = copy.deepcopy(tpatt_dict), copy.deepcopy(tfinal_dict)
        return patt_dict, final_dict

    def solve_patterns(self, slit, bestlist, tcent_ecent):
        """

        Args:
            slit (int):
               The ID of the slit
            bestlist (list or ndarray):
                A 5 element list, each containing a numpy.ndarray, with the following values required for each index:
                    0: central wavelength of the pattern
                    1: central dispersion of pattern
                    2: sign of the pattern (note, sign = 1 [-1] if pixels correlate [anticorrelate] with wavelength
                    3: index of the full list of patterns that were created from the detected arc lines
                    4: index of the full list of patterns that were created from the line list.
            tcent_ecent (list):
                A list [tcent, ecent] indicating which detection list should be used. Note that if arr_err is set then the weak keyword is ignored.

        Returns:
            patt_dict (dict):
                Dictionary containing information about the best patterns.

        """

        # Obtain a full list of indices that are consistent with the maximum value
        wcen, dcen, sign, dindex, lindex = bestlist[0], bestlist[1], bestlist[3], bestlist[4], bestlist[5]

        # Find the favoured sign and only use those values
        use_tcent, _ = self.get_use_tcent(sign, tcent_ecent)
        if sign == +1:
            signtxt = "correlate"
        else:
            signtxt = "anticorrelate"

        # Initialise the patterns dictionary
        patt_dict = dict(acceptable=False, nmatch=0, ibest=-1, bwv=0.,
                         sigdetect=wvutils.parse_param(self._par, 'sigdetect', slit),
                         mask=np.zeros(use_tcent.size, dtype=bool))
        patterns.solve_triangles(use_tcent, self._wvdata, dindex, lindex, patt_dict)
        # Check if a solution was found
        if not patt_dict['acceptable']:
            return None

        # Fill in the patterns dictionary
        patt_dict['sign'] = sign
        patt_dict['bwv'] = wcen
        patt_dict['bdisp'] = 10.0 ** dcen

        # Check that a solution has been found
        if patt_dict['nmatch'] == 0 and self._verbose:
            msgs.info(msgs.newline() +
                      '---------------------------------------------------' + msgs.newline() +
                      'Initial report:' + msgs.newline() +
                      '  No matches! Try another algorithm' + msgs.newline() +
                      '---------------------------------------------------')
            return None
        elif self._verbose:
            # Report
            msgs.info(msgs.newline() +
                      '---------------------------------------------------' + msgs.newline() +
                      'Initial report:' + msgs.newline() +
                      '  Pixels {:s} with wavelength'.format(signtxt) + msgs.newline() +
                      '  Number of lines recovered    = {:d}'.format(self._all_tcent.size) + msgs.newline() +
                      '  Number of lines analyzed     = {:d}'.format(use_tcent.size) + msgs.newline() +
                      '  Number of acceptable matches = {:d}'.format(patt_dict['nmatch']) + msgs.newline() +
                      '  Best central wavelength      = {:g}A'.format(patt_dict['bwv']) + msgs.newline() +
                      '  Best dispersion              = {:g}A/pix'.format(patt_dict['bdisp']) + msgs.newline() +
                      '  Best wave/disp                = {:g}'.format(patt_dict['bwv']/patt_dict['bdisp']) + msgs.newline() +
                      '---------------------------------------------------')
        return patt_dict

    def finalize_fit(self, detections):
        """
        Once the best IDs have been found for each slit, perform a final
        fit to all slits and save the results
        """

        for slit in range(self._nslit):
            if slit not in self._ok_mask:
                continue
            if self._all_patt_dict[str(slit)] is None:
                continue
            # Save the QA for the best solution
            slittxt = '_Slit{0:03d}'.format(slit+1)
            use_tcent, use_ecent = self.get_use_tcent(self._all_patt_dict[str(slit)]['sign'],detections[str(slit)])
            if self._outroot is not None:
                # Write IDs
                out_dict = dict(pix=use_tcent, IDs=self._all_patt_dict[str(slit)]['IDs'])
                jdict = ltu.jsonify(out_dict)
                ltu.savejson(self._outroot + slittxt + '.json', jdict, easy_to_read=True, overwrite=True)
                msgs.info("Wrote: {:s}".format(self._outroot + slittxt + '.json'))

                # Plot
                tmp_list = np.vstack([self._line_lists, self._unknwns])
                match_qa(self._spec[:, slit], use_tcent, tmp_list,
                            self._all_patt_dict[str(slit)]['IDs'], self._all_patt_dict[str(slit)]['scores'],
                            outfile=self._outroot + slittxt + '.pdf')
                msgs.info("Wrote: {:s}".format(self._outroot + slittxt + '.pdf'))
            # Perform the final fit for the best solution
            best_final_fit = wv_fitting.fit_slit(self._spec[:, slit], self._all_patt_dict[str(slit)], use_tcent,
                                                 self._line_lists, outroot=self._outroot, slittxt=slittxt)
            #best_final_fit = self.fit_slit(slit, self._all_patt_dict[str(slit)], use_tcent, outroot=self._outroot, slittxt=slittxt)
            self._all_final_fit[str(slit)] = copy.deepcopy(best_final_fit)

    def report_prelim(self, slit, best_patt_dict, best_final_fit):

        good_fit = False
        # Report on the best preliminary result
        if best_final_fit is None:
            msgs.warn('---------------------------------------------------' + msgs.newline() +
                      'Preliminary report for slit {0:d}/{1:d}:'.format(slit, self._nslit-1) + msgs.newline() +
                      '  No matches! Attempting to cross match.' + msgs.newline() +
                      '---------------------------------------------------')
            self._all_patt_dict[str(slit)] = None
            self._all_final_fit[str(slit)] = None
        elif best_final_fit['rms'] > wvutils.parse_param(self._par, 'rms_threshold', slit):
            msgs.warn('---------------------------------------------------' + msgs.newline() +
                      'Preliminary report for slit {0:d}/{1:d}:'.format(slit, self._nslit-1) + msgs.newline() +
                      '  Poor RMS ({0:.3f})! Attempting to cross match.'.format(best_final_fit['rms']) + msgs.newline() +
                      '---------------------------------------------------')
            self._all_patt_dict[str(slit)] = None
            self._all_final_fit[str(slit)] = None
        else:
            good_fit = True
            if best_patt_dict['sign'] == +1:
                signtxt = 'correlate'
            else:
                signtxt = 'anitcorrelate'
            # Report
            msgs.info('---------------------------------------------------' + msgs.newline() +
                      'Preliminary report for slit {0:d}/{1:d}:'.format(slit, self._nslit-1) + msgs.newline() +
                      '  Pixels {:s} with wavelength'.format(signtxt) + msgs.newline() +
                      '  Number of weak lines         = {:d}'.format(self._det_weak[str(slit)][0].size) + msgs.newline() +
                      '  Number of strong lines       = {:d}'.format(self._det_stro[str(slit)][0].size) + msgs.newline() +
                      '  Number of lines analyzed     = {:d}'.format(len(best_final_fit['pixel_fit'])) + msgs.newline() +
                      '  Number of pattern matches    = {:d}'.format(best_patt_dict['nmatch']) + msgs.newline() +
                      '  Patt match cen wavelength    = {:g}A'.format(best_patt_dict['bwv']) + msgs.newline() +
                      '  Patt match dispersion        = {:g}A/pix'.format(best_patt_dict['bdisp']) + msgs.newline() +
                      '  Best patt match wave/disp    = {:g}'.format(best_patt_dict['bwv']/best_patt_dict['bdisp']) + msgs.newline() +
                      '  Final RMS of fit             = {:g}'.format(best_final_fit['rms']) + msgs.newline() +
                      '---------------------------------------------------')
            self._all_patt_dict[str(slit)] = copy.deepcopy(best_patt_dict)
            self._all_final_fit[str(slit)] = copy.deepcopy(best_final_fit)
        return good_fit

    def report_final(self):
        """Print out the final report of the wavelength calibration"""
        for slit in range(self._nslit):
            # Prepare a message for bad wavelength solutions
            badmsg = '---------------------------------------------------' + msgs.newline() +\
                     'Final report for slit {0:d}/{1:d}:'.format(slit, self._nslit-1) + msgs.newline() +\
                     '  Wavelength calibration not performed!'
            if slit not in self._ok_mask:
                msgs.warn(badmsg)
                continue
            if self._all_patt_dict[str(slit)] is None:
                msgs.warn(badmsg)
                continue
            st = str(slit)
            if self._all_patt_dict[st]['sign'] == +1:
                signtxt = 'correlate'
            else:
                signtxt = 'anitcorrelate'
            # Report
            centwave = self._all_final_fit[st].pypeitfit.eval(0.5)
            tempwave = self._all_final_fit[st].pypeitfit.eval(0.5 + 1.0/self._npix)
            centdisp = abs(centwave-tempwave)
            msgs.info(msgs.newline() +
                      '---------------------------------------------------' + msgs.newline() +
                      'Final report for slit {0:d}/{1:d}:'.format(slit, self._nslit-1) + msgs.newline() +
                      '  Pixels {:s} with wavelength'.format(signtxt) + msgs.newline() +
                      '  Number of weak lines         = {:d}'.format(self._det_weak[str(slit)][0].size) + msgs.newline() +
                      '  Number of strong lines       = {:d}'.format(self._det_stro[str(slit)][0].size) + msgs.newline() +
                      '  Number of lines analyzed     = {:d}'.format(len(self._all_final_fit[st]['pixel_fit'])) + msgs.newline() +
                      '  Central wavelength           = {:g}A'.format(centwave) + msgs.newline() +
                      '  Central dispersion           = {:g}A/pix'.format(centdisp) + msgs.newline() +
                      '  Central wave/disp             = {:g}'.format(centwave/centdisp) + msgs.newline() +
                      '  Final RMS of fit             = {:g}'.format(self._all_final_fit[st]['rms']))
        return


def results_kdtree_nb(use_tcent, wvdata, res, residx, dindex, lindex, nindx, npix, ordfit=1):
    """ A numba speedup of the results_kdtree function in the General class (see above).
    For all of the acceptable pattern matches, estimate the central wavelength and dispersion,
    and record the index in the linelist and the corresponding indices of the detected lines.

    Parameters
    ----------
    use_tcent : ndarray
        detected lines
    wvdata : ndarray
        the linelist
    res : list
        A flattened list of the results of the ball tree query from the
        KDTree (this contains all acceptable matches) This needs to be a
        flattened list for numba
    residx : list
        This contains the original indices of the unflattened (i.e. nested) 'res' list
    dindex : ndarray
        Indices of the lines in the detected lines for all patterns
    lindex : ndarray
        Indices of the lines in the linelist for all patterns
    nindx : int
        Number of acceptable pattens
    npix : int
        Number of pixels in the spectral direction
    ordfit : int
        Order of the polynomial used to fit the pixel/wavelength IDs

    Returns
    -------
    dind : ndarray
        Indices of the lines in the detected lines that were used for each acceptable pattern
    lind : linelist index of patterns
        Indices of the lines in the linelist that were used for each corresponding pattern
    wvcent : ndarray
        Central wavelength of each pattern
    wvdisp : ndarray
        Central dispersion of each pattern
    """
    # Assign wavelengths to each pixel
    ncols = len(res)
    wvdisp = np.zeros(ncols, dtype=float)
    wvcent = np.zeros(ncols, dtype=float)
    dind = np.zeros((ncols, nindx), dtype=np.uint64)
    lind = np.zeros((ncols, nindx), dtype=np.uint64)
    Xmat = np.ones((nindx, ordfit+1), dtype=float)
    for x in range(ncols):
        for ii in range(ordfit, -1, -1):
            Xmat[:, ii] = np.power(use_tcent[dindex[residx[x], :]], ordfit-ii)
        coeff = np.linalg.lstsq(Xmat, wvdata[lindex[res[x]]])[0]
        sumv = 0.0
        for ii in range(ordfit, -1, -1):
            wvcent[x] += coeff[ordfit-ii]*((npix/2.0)**ii)
            sumv += coeff[ordfit-ii]*(((npix+1)/2.0)**ii)
        wvdisp[x] = abs(sumv - wvcent[x])
        dind[x, :] = dindex[residx[x], :]
        lind[x, :] = lindex[res[x], :]
    return dind, lind, wvcent, wvdisp


