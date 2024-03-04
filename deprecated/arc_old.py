from __future__ import (print_function, absolute_import, division, unicode_literals)

import inspect

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt


import scipy
from astropy.stats import sigma_clipped_stats, sigma_clip

from pypeit import ararclines
from pypeit import debugger
from pypeit import msgs
from pypeit import utils
from pypeit.core import parse
from pypeit.core import pixels
from pypeit.core.wavecal import autoid
from pypeit import debugger
from pypeit.core import pydl
from pypeit.core import qa
from skimage.transform import resize


def fit_double_poly_old(all_wv_order, work2d, thismask, nspec_coeff, norder_coeff):
    """ This perform the actual fit of the 2D wavelength solution of an echelle spectrograph.

    Parameters
    ----------
    all_wv_order: np.array
     wavelength*order of the identified lines
    work2d: np.array
      matrix containing the coefficient of the legendre polinomial
    thismask: boolean
      mask for good lines
    nspec_coeff : np.int
      order of the fitting along the spectral (pixel) direction for each order
    norder_coeff : np.int
      order of the fitting in the order direction

    Returns:
    -------
      coeffs, wv_order_mod

    """
    work2di = np.transpose(work2d * np.outer(np.ones(norder_coeff * nspec_coeff, dtype=np.float64), thismask))
    alpha = work2d.dot(work2di)
    beta = all_wv_order.dot(work2di)
    coeffs = np.linalg.solve(alpha, beta)
    wv_order_mod = coeffs.dot(work2d)
    return coeffs, wv_order_mod


def fit2darc_old(all_wv, all_pix, all_orders, nspec, nspec_coeff=4, norder_coeff=4, sigrej=3.0, debug=False):
    """Routine to obtain the 2D wavelength solution for an echelle spectrograph. This is calculated from the spec direction
    pixelcentroid and the order number of identified arc lines. The fit is a simple least-squares with rejections.
    This is a port of the XIDL code: x_fit2darc.pro

    Parameters
    ----------
    all_wv: np.array
     wavelength of the identified lines
    all_pix: np.array
      y-centroid position of the identified lines
    all_orders: np.array
      order number of the identified lines
    nspec: int
      Size of the image in the spectral direction
    nspec_coeff : np.int
      order of the fitting along the spectral (pixel) direction for each order
    norder_coeff : np.int
      order of the fitting in the order direction
    sigrej: np.float
      sigma level for the rejection
    debug: boolean
      Extra plots to check the status of the procedure

    Returns:
    -------
    """


    # To use the legendre polynomial pixels and orders
    # need to be normalized in the -1,+1 range
    # Normalize pixels
    mnx = 0  # np.min(all_pix)
    mxx = float(nspec - 1)  # np.max(all_pix)
    norm_pixel = np.array([0.5 * (mnx + mxx), mxx - mnx])
    pix_nrm = 2. * (all_pix - norm_pixel[0]) / norm_pixel[1]
    # Normalize orders
    mnx, mxx = np.min(all_orders), np.max(all_orders)
    norm_order = np.array([0.5 * (mnx + mxx), mxx - mnx])
    orders_nrm = 2. * (all_orders - norm_order[0]) / norm_order[1]

    if debug:
        # set some plotting parameters
        utils.pyplot_rcparams()
        plt.figure(figsize=(7, 5))
        msgs.info("Plot identified lines")
        cm = plt.cm.get_cmap('RdYlBu_r')
        sc = plt.scatter(orders_nrm, pix_nrm, c=all_wv / 10000., cmap=cm)
        cbar = plt.colorbar(sc)
        cbar.set_label(r'Wavelength [$\mu$m]', rotation=270,
                       labelpad=20)
        plt.xlabel(r'Normalized Orders')
        plt.ylabel(r'Normalized Pixels')
        plt.title(r'Location of the identified lines')
        plt.show()

    # Setup some things for the fits
    all_wv_order = all_wv * all_orders
    work2d = np.zeros((nspec_coeff * norder_coeff, len(all_wv)), dtype=np.float64)
    worky = pydl.flegendre(pix_nrm, nspec_coeff)
    workt = pydl.flegendre(orders_nrm, norder_coeff)
    for i in range(norder_coeff):
        for j in range(nspec_coeff):
            work2d[j * norder_coeff + i, :] = worky[j, :] * workt[i, :]

    # ToDO add upper lower to inputs
    lower = np.abs(sigrej)
    upper = np.abs(sigrej)
    maxiter = 25
    iIter = 0
    qdone = False
    thismask = np.ones_like(all_wv, dtype=bool)
    while (not qdone) and (iIter < maxiter):
        coeffs, wv_order_mod = fit_double_poly(all_wv_order, work2d, thismask.astype(float), nspec_coeff, norder_coeff)
        thismask, qdone = pydl.djs_reject(all_wv_order, wv_order_mod, outmask=thismask,
                                          lower=np.float64(lower), upper=np.float64(upper), use_mad=True, sticky=True)
        iIter += 1
        if debug:
            utils.pyplot_rcparams()
            plt.figure(figsize=(7, 5))
            plt.axhline(y=np.average(wv_order_mod[thismask] / all_orders[thismask] - all_wv[thismask]), color='r',
                        linestyle='--')
            plt.axhline(y=+np.std(wv_order_mod[thismask] / all_orders[thismask] - all_wv[thismask]), color='r',
                        linestyle=':')
            plt.axhline(y=-np.std(wv_order_mod[thismask] / all_orders[thismask] - all_wv[thismask]), color='r',
                        linestyle=':')
            plt.scatter(all_wv[~thismask] / 10000., wv_order_mod[~thismask] / all_orders[~thismask] - all_wv[~thismask],
                        marker="v", label=r'Rejected values')
            plt.scatter(all_wv[thismask] / 10000., wv_order_mod[thismask] / all_orders[thismask] - all_wv[thismask],
                        marker="v", label=r'Good values')
            plt.text(np.min(all_wv / 10000),
                     np.average(wv_order_mod[thismask] / all_orders[thismask] - all_wv[thismask]),
                     r'Average={0:.1f}$\AA$'.format(
                         np.average(wv_order_mod[thismask] / all_orders[thismask] - all_wv[thismask])),
                     ha="left", va="bottom",
                     bbox=dict(boxstyle="square", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8), alpha=0.7, ))
            plt.text(np.max(all_wv / 10000), np.std(wv_order_mod[thismask] / all_orders[thismask] - all_wv[thismask]),
                     r'Sigma={0:.1f}$\AA$'.format(
                         np.std(wv_order_mod[thismask] / all_orders[thismask] - all_wv[thismask])), ha="right",
                     va="bottom", bbox=dict(boxstyle="square", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8), alpha=0.7, ))
            plt.legend()
            plt.title(r'Residuals after rejection iteration #{:d}'.format(iIter))
            plt.xlabel(r'Wavelength [$\mu$m]')
            plt.ylabel(r'Residuals [$\AA$]')
            plt.show()
    if iIter == maxiter:
        msgs.warn('Maximum number of iterations maxiter={:}'.format(maxiter) + ' reached in robust_polyfit_djs')

    # Final fit
    coeffs, wv_order_mod = fit_double_poly(all_wv_order, work2d,
                                           thismask.astype(float),
                                           nspec_coeff, norder_coeff)

    # Check quality
    resid = (wv_order_mod[thismask] - all_wv_order[thismask])
    fin_rms = np.sqrt(np.mean(resid ** 2))
    msgs.info("RMS: {0:.5f} Ang*Order#".format(fin_rms))

    orders = np.unique(all_orders)
    fit_dict = dict(coeffs=coeffs, orders=orders,
                    nspec_coeff=nspec_coeff, norder_coeff=norder_coeff,
                    pixel_cen=norm_pixel[0], pixel_norm=norm_pixel[1],
                    order_cen=norm_order[0], order_norm=norm_order[1],
                    nspec=nspec, all_pix=all_pix, all_wv=all_wv,
                    all_orders=all_orders, all_mask=thismask)

    if debug:
        fit2darc_global_qa(fit_dict)
        fit2darc_orders_qa(fit_dict)

    return fit_dict


def fit2darc_global_qa_old(fit_dict, outfile=None):
    """ QA on 2D fit of the wavelength solution.

    Parameters
    ----------
    fit_dict: dict
      dict of the 2D arc solution
    outfile:
      parameter for QA

    Returns
    -------
    """

    msgs.info("Creating QA for 2D wavelength solution")

    utils.pyplot_rcparams()

    # Extract info from fit_dict
    nspec = fit_dict['nspec']
    orders = fit_dict['orders']
    pixel_norm = fit_dict['pixel_norm']
    pixel_cen = fit_dict['pixel_cen']
    nspec_coeff = fit_dict['nspec_coeff']
    norder_coeff = fit_dict['norder_coeff']
    all_wv = fit_dict['all_wv']
    all_pix = fit_dict['all_pix']
    all_orders = fit_dict['all_orders']
    thismask = fit_dict['all_mask']
    resid_wl_global = []

    # Define pixels array
    all_pixels = np.arange(nspec)

    # Define figure properties
    plt.figure(figsize=(8 ,5))

    # Variable where to store the max wavelength covered by the
    # spectrum
    mx = 0.

    # Loop over orders
    for ii in orders:

        # define the color
        rr = (i i -np.max(orders) ) /(np.min(orders ) -np.max(orders))
        gg = 0.0
        bb = (i i -np.min(orders) ) /(np.max(orders ) -np.min(orders))

        # evaluate solution
        wv_order_mod = eval2dfit(fit_dict, all_pixels, ii)

        # Plot solution
        plt.plot(wv_order_mo d /ii, all_pixels ,color=(rr ,gg ,bb),
                 linestyle='-', linewidth=2.5)

        # Evaluate residuals at each order
        this_pix = all_pix[all_orders == ii]
        this_wv = all_wv[all_orders == ii]
        this_msk = thismask[all_orders == ii]
        wv_order_mod_resid = eval2dfit(fit_dict, this_pix, ii)
        resid_wl = (wv_order_mod_resid /i i -this_wv)
        resid_wl_global = np.append(resid_wl_global ,resid_wl[this_msk])
        plt.scatter((wv_order_mod_resid[~this_msk ] /ii )+ \
                    100 . *resid_wl[~this_msk], this_pix[~this_msk], \
                    marker='x', color='black', linewidths=2.5, s=16.)
        plt.scatter((wv_order_mod_resid[this_msk ] /ii )+ \
                    100 . *resid_wl[this_msk], this_pix[this_msk], \
                    color=(rr ,gg ,bb), linewidth=2.5, s=16.)
        if np.max(wv_order_mod_resi d /ii) > mx :
            mx = np.max(wv_order_mod_resi d /ii)

    rms_global = np.sqrt(np.mean((resid_wl_global )* *2))

    plt.text(mx ,np.max(all_pixels) ,r'residuals $\times$100', \
             ha="right", va="top")
    plt.title(r'Arc 2D FIT, norder_coeff={:d}, nspec_coeff={:d}, RMS={:5.3f} Ang*Order#'.format(
        norder_coeff, nspec_coeff, rms_global))
    plt.xlabel(r'Wavelength [$\AA$]')
    plt.ylabel(r'Row [pixel]')

    # Finish
    if outfile is not None:
        plt.savefig(outfile, dpi=800)
        plt.close()
    else:
        plt.show()

    # restore default rcparams
    utils.pyplot_rcparams_default()


def fit2darc_orders_qa_old(fit_dict, outfile=None):
    """ QA on 2D fit of the wavelength solution of an Echelle spectrograph.
    Each panel contains a single order with the global fit and the
    residuals.

    Parameters
    ----------
    fit_dict: dict
      dict of the 2D arc solution
    outfile:
      parameter for QA

    Returns
    -------
    """

    msgs.info("Creating QA for 2D wavelength solution")

    utils.pyplot_rcparams()

    # Extract info from fit_dict
    nspec = fit_dict['nspec']
    orders = fit_dict['orders']
    pixel_norm = fit_dict['pixel_norm']
    pixel_cen = fit_dict['pixel_cen']
    nspec_coeff = fit_dict['nspec_coeff']
    norder_coeff = fit_dict['norder_coeff']
    all_wv = fit_dict['all_wv']
    all_pix = fit_dict['all_pix']
    all_orders = fit_dict['all_orders']
    thismask = fit_dict['all_mask']
    resid_wl_global = []

    # Define pixels array
    all_pixels = np.arange(nspec)

    # set the size of the plot
    nrow = np.int(2)
    ncol = np.int(np.ceil(len(orders ) /2.))
    fig = plt.figure(figsize=( 5 *ncol , 6 *nrow))

    outer = gridspec.GridSpec(nrow, ncol, wspace=0.3, hspace=0.2)

    for ii_row in range(nrow):
        for ii_col in range(ncol):
            if (ii_ro w *ncol + ii_col) < len(orders):
                inner = gridspec.GridSpecFromSubplotSpec(2, 1,
                                                         height_ratios=[2 ,1], width_ratios=[1],
                                                         subplot_spec=outer[ii_ro w *ncol + ii_col],
                                                         wspace=0.1, hspace=0.0)
                ax0 = plt.Subplot(fig, inner[0])
                ax1 = plt.Subplot(fig, inner[1], sharex=ax0)
                plt.setp(ax0.get_xticklabels(), visible=False)

                ii = orders[ii_ro w *ncol + ii_col]

                # define the color
                rr = (i i -np.max(orders) ) /(np.min(orders ) -np.max(orders))
                gg = 0.0
                bb = (i i -np.min(orders) ) /(np.max(orders ) -np.min(orders))

                # Evaluate function
                wv_order_mod = eval2dfit(fit_dict, all_pixels, ii)
                # Evaluate delta lambda
                dw l =(wv_order_mod[-1 ] -wv_order_mod[0] ) /i i /(all_pixels[-1 ] -all_pixels[0])

                # Estimate the residuals
                this_pix = all_pix[all_orders == ii]
                this_wv = all_wv[all_orders == ii]
                this_msk = thismask[all_orders == ii]

                wv_order_mod_resid = eval2dfit(fit_dict, this_pix, ii)
                resid_wl = (wv_order_mod_resi d /i i -this_wv)
                resid_wl_global = np.append(resid_wl_global ,resid_wl[this_msk])

                # Plot the fit
                ax0.set_title('Order = {0:0.0f}'.format(ii))
                ax0.plot(all_pixels, wv_order_mo d /i i /10000. ,color=(rr ,gg ,bb), linestyle='-',
                         linewidth=2.5)
                ax0.scatter(this_pix[~this_msk], (wv_order_mod_resid[~this_msk ] /i i /10000. )+ \
                            100 . *resid_wl[~this_msk ] /10000., marker='x', color='black', \
                            linewidth=2.5, s=16.)
                ax0.scatter(this_pix[this_msk], (wv_order_mod_resid[this_msk ] /i i /10000. )+ \
                            100 . *resid_wl[this_msk ] /10000., color=(rr ,gg ,bb), \
                            linewidth=2.5, s=16.)

                ax0.set_ylabel(r'Wavelength [$\mu$m]')

                # Plot the residuals
                ax1.scatter(this_pix[~this_msk] ,(resid_wl[~this_msk ] /dwl) ,marker='x', color='black', \
                            linewidth=2.5, s=16.)
                ax1.scatter(this_pix[this_msk], (resid_wl[this_msk ] /dwl), color=(rr ,gg ,bb), \
                            linewidth=2.5, s=16.)
                ax1.axhline(y=0., color=(rr ,gg ,bb), linestyle=':', linewidth=2.5)
                ax1.get_yaxis().set_label_coords(-0.15 ,0.5)

                rms_order = np.sqrt(np.mean((resid_wl[this_msk] )* *2))

                ax1.set_ylabel(r'Res. [pix]')

                ax0.text(0.1 ,0.9 ,r'RMS={0:.3f} Pixel'.format(rms_orde r /np.abs(dwl)) ,ha="left", va="top",
                         transform = ax0.transAxes)
                ax0.text(0.1 ,0.8 ,r'$\Delta\lambda$={0:.3f} Pixel/$\AA$'.format(np.abs(dwl)) ,ha="left", va="top",
                         transform = ax0.transAxes)
                ax0.get_yaxis().set_label_coords(-0.15 ,0.5)

                fig.add_subplot(ax0)
                fig.add_subplot(ax1)

    rms_global = np.sqrt(np.mean((resid_wl_global )* *2))

    fig.text(0.5, 0.04, r'Row [pixel]', ha='center', size='large')
    fig.suptitle \
        (r'Arc 2D FIT, norder_coeff={:d}, nspec_coeff={:d}, RMS={:5.3f} Ang*Order#, residuals $\times$100'.format
            (norder_coeff,
                                                                                                                          nspec_coeff ,rms_global))

    # Finish
    if outfile is not None:
        plt.savefig(outfile, dpi=800)
        plt.close()
    else:
        plt.show()


def eval2dfit_old(fit_dict, pixels, order):
    """ Evaluate the 2D fit at a given pixel and order.

    Parameters
    ----------
    fit_dict: dict
      dictionary containing the result of the fit
    pixels: np.array
      pixels where you want to evaluate the fit
    order: np.array
      order where you want to evaluate the fit

    Returns
    -------
    wv_order_mod
      wavelength*order evaluated at the given pixel and order
    """

    if pixels.ndim != 1:
        msgs.error('pixels must be a one dimensional array')

    nspec_coeff = fit_dict['nspec_coeff']
    norder_coeff = fit_dict['norder_coeff']
    coeffs = fit_dict['coeffs']
    npix = pixels.size
    # legendre basis for the order direction
    osub = np.ones_like(pixels, dtype=np.float64) * order
    order_nrm = 2. 0 *(osub - fit_dict['order_cen'] ) /fit_dict['order_norm']
    work_order = pydl.flegendre(order_nrm, norder_coeff)
    # legendre basis for the spectral direction
    pix_nrm = 2. 0 *(pixels - fit_dict['pixel_cen'] ) /fit_dict['pixel_norm']
    work_pix = pydl.flegendre(pix_nrm, nspec_coeff)
    # array to hold the
    work2d = np.zeros((nspec_coef f *norder_coeff, npix), dtype=float)
    for i in range(norder_coeff):
        for j in range(nspec_coeff):
            work2d[ j *norder_coeff + i, :] = work_pix[j, : ] *work_order[i, :]
    wv_order_mod = coeffs.dot(work2d)

    return wv_order_mod


# TODO Make det_arxiv an optional input. In the event that the Line IDs don't exist in the arxiv, simply run peak finding and
# interpolate the wavelength solution onto those line locations to the get initial IDs

def reidentify(spec, spec_arxiv_in, wave_soln_arxiv_in, det_arxiv, line_list, nreid_min, detections=None, cc_thresh=0.8,cc_local_thresh = 0.8,
               match_toler=2.0, nlocal_cc=11, nonlinear_counts=1e10,sigdetect=5.0,fwhm=4.0, debug_xcorr=False, debug_reid=False, debug_peaks = False):
    """ Determine  a wavelength solution for a set of spectra based on archival wavelength solutions

    Parameters
    ----------
    spec:  float ndarray shape (nspec)
       Arc spectrum for which wavelength identifications are desired.

    spec_arxiv:  float ndarray shape (nspec, narxiv) or (nspec)
       Collection of archival arc spectra for which wavelength solution and line identifications are known

    wave_soln_arxiv:  float ndarray shape (nspec, narxiv) or (nspec)
       Wavelength solutions for the archival arc spectra spec_arxiv

    det_arxiv:  dict, the dict has narxiv keys which are '0','1', ... up to str(narxiv-1). det_arxiv['0'] points to an
                an ndarray of size determined by the number of lines that were detected.

       Arc line pixel locations in the spec_arxiv spectra that were used in combination with line identifications from the
       line list to determine the wavelength solution wave_soln_arxiv.

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
    (detections, patt_dict)

    detections: ndarray,
       Pixel locations of arc lines detected.
    patt_dict: dict
       Arc lines pattern dictionary with some information about the IDs as well as the cross-correlation values

    Revision History
    ----------------
    November 2018 by J.F. Hennawi. Built from an initial version of cross_match code written by Ryan Cooke.
    """

    # Determine the seed for scipy.optimize.differential_evolution optimizer. Just take the sum of all the elements
    # and round that to an integer
    from IPython import embed

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

    # TODO Implement different binnings between archival and input spectrum
    if spec_arxiv_in.ndim == 1:
        spec_arxiv = spec_arxiv_in.reshape(spec_arxiv_in.size,1)
        wave_soln_arxiv = wave_soln_arxiv_in.reshape(wave_soln_arxiv_in.size,1)
    elif spec_arxiv_in.ndim == 2:
        spec_arxiv = spec_arxiv_in.copy()
        wave_soln_arxiv = wave_soln_arxiv_in.copy()
    else:
        msgs.error('Unrecognized shape for spec_arxiv. It must be either a one dimensional or two dimensional numpy array')

    nspec_arxiv, narxiv = spec_arxiv.shape

    this_soln = wave_soln_arxiv[:,0]
    sign = 1 if (this_soln[this_soln.size // 2] > this_soln[this_soln.size // 2 - 1]) else -1



    xrng = np.arange(nspec_arxiv)
    if nspec_arxiv != nspec:
        msgs.error('Different spectral binning is not supported yet but it will be soon')

    # Search for lines to continuum subtract the spectrum.
    tcent, ecent, cut_tcent, icut, spec_cont_sub = wvutils.arc_lines_from_spec(spec, sigdetect=sigdetect,nonlinear_counts=nonlinear_counts, fwhm = fwhm, debug = debug_peaks)
    # If the detections were not passed in assing them
    if detections is None:
        detections = tcent[icut]

    # For convenience pull out all the spectra from the wv_calib_arxiv archive
    wvc_arxiv = np.zeros(narxiv, dtype=float)
    disp_arxiv = np.zeros(narxiv, dtype=float)
    for iarxiv in range(narxiv):
        wvc_arxiv[iarxiv] = wave_soln_arxiv[nspec_arxiv//2, iarxiv]
        disp_arxiv[iarxiv] = np.median(wave_soln_arxiv[:,iarxiv] - np.roll(wave_soln_arxiv[:,iarxiv], 1))

    marker_tuple = ('o','v','<','>','8','s','p','P','*','X','D','d','x')
    color_tuple = ('black','green','red','cyan','magenta','blue','darkorange','yellow','dodgerblue','purple','lightgreen','cornflowerblue')
    marker = itertools.cycle(marker_tuple)
    colors = itertools.cycle(color_tuple)

    # Cross-correlate with each arxiv spectrum to identify lines
    line_indx = np.array([], dtype=np.int)
    det_indx = np.array([], dtype=np.int)
    line_cc = np.array([], dtype=float)
    line_iarxiv = np.array([], dtype=np.int)
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
            wvutils.xcorr_shift_stretch(spec_cont_sub, spec_arxiv[:, iarxiv], cc_thresh=cc_thresh, fwhm = fwhm, seed = random_state,
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
        prod_smooth = scipy.ndimage.filters.convolve1d(spec_cont_sub*spec_arxiv_ss, window)
        spec2_smooth = scipy.ndimage.filters.convolve1d(spec_cont_sub**2, window)
        arxiv2_smooth = scipy.ndimage.filters.convolve1d(spec_arxiv_ss**2, window)
        denom = np.sqrt(spec2_smooth*arxiv2_smooth)
        corr_local = np.zeros_like(denom)
        corr_local[denom > 0] = prod_smooth[denom > 0]/denom[denom > 0]
        corr_local[denom == 0.0] = -1.0

        # Loop over the current slit line pixel detections and find the nearest arxiv spectrum line
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


# MOVE TO DEPRECATED
def simple_calib_driver(llist, censpec, ok_mask, n_final=5, get_poly=False,
                        sigdetect=10.,
                        IDpixels=None, IDwaves=None, nonlinear_counts=1e10):
    wv_calib = {}
    for slit in ok_mask:
        iwv_calib = simple_calib(llist, censpec[:, slit], n_final=n_final,
                                 get_poly=get_poly, IDpixels=IDpixels, IDwaves=IDwaves,
                                 nonlinear_counts=nonlinear_counts, sigdetect=sigdetect)
        wv_calib[str(slit)] = iwv_calib.copy()
    return wv_calib


def simple_calib(llist, censpec, n_final=5, get_poly=False,
                 IDpixels=None, IDwaves=None, debug=False, sigdetect=10.,
                 nonlinear_counts=1e10):
    """Simple calibration algorithm for longslit wavelengths

    Parameters
    ----------
    llist : `astropy.table.Table`_
    censpec : `numpy.ndarray`_
    get_poly : bool, optional
      Pause to record the polynomial pix = b0 + b1*lambda + b2*lambda**2
    IDpixels : list
    IDwaves : list

    Returns
    -------
    final_fit : dict
      Dict of fit info
    """

    # Extract the arc
    msgs.work("Detecting lines..")
    #tampl, tcent, twid, _, w, yprep, nsig = detect_lines(censpec, nfitpix=nfitpix,
    #                                                     sigdetect=sigdetect,
    #                                                     nonlinear_counts = aparm['nonlinear_counts'])
    tcent, ecent, cut_tcent, icut, spec_cont_sub = wvutils.arc_lines_from_spec(
        censpec, sigdetect=sigdetect, nonlinear_counts=nonlinear_counts)#, debug = debug_peaks)

    # Cut down to the good ones
    tcent = tcent[icut]

    # IDs were input by hand
    # Check that there are at least 4 values
    pixels = np.array(IDpixels) # settings.argflag['arc']['calibrate']['IDpixels'])
    if np.sum(pixels > 0.) < 4:
        msgs.error("Need to give at least 4 pixel values!")
    #
    msgs.info("Using input lines to seed the wavelength solution")
    # Calculate median offset
    mdiff = [np.min(np.abs(tcent-pix)) for pix in pixels]
             #settings.argflag['arc']['calibrate']['IDpixels']]
    med_poff = np.median(np.array(mdiff))
    msgs.info("Will apply a median offset of {:g} pixels".format(med_poff))

    # Match input lines to observed spectrum
    nid = pixels.size
    idx_str = np.ones(nid).astype(int)
    ids = np.zeros(nid)
    idsion = np.array(['     ']*nid)
    gd_str = np.arange(nid).astype(int)
    for jj,pix in enumerate(pixels):
        diff = np.abs(tcent-pix-med_poff)
        if np.min(diff) > 2.:
            msgs.error("No match with input pixel {:g}!".format(pix))
        else:
            imn = np.argmin(diff)
        # Set
        idx_str[jj] = imn
        # Take wavelength from linelist instead of input value
        wdiff = np.abs(llist['wave']-IDwaves[jj]) # settings.argflag['arc']['calibrate']['IDwaves'][jj])
        imnw = np.argmin(wdiff)
        if wdiff[imnw] > 0.015:  # Arbitrary tolerance
            msgs.error("Input IDwaves={:g} is not in the linelist.  Fix".format(
                IDwaves[jj]))
                    #settings.argflag['arc']['calibrate']['IDwaves'][jj]))
        else:
            ids[jj] = llist['wave'][imnw]
            #idsion[jj] = llist['Ion'][imnw]
            msgs.info("Identifying arc line: {:s} {:g}".format(idsion[jj],ids[jj]))

    # Debug
    disp = (ids[-1]-ids[0])/(tcent[idx_str[-1]]-tcent[idx_str[0]])
    final_fit = wv_fitting.iterative_fitting(censpec, tcent, idx_str, ids,
                                          llist, disp, verbose=False, n_final=n_final)
    # Return
    return final_fit


def order_saturation(satmask, ordcen, ordwid):
    """
    .. todo::
        Document this!
    """
    sz_y, sz_x = satmask.shape
    sz_o = ordcen.shape[1]

    xmin = ordcen - ordwid
    xmax = ordcen + ordwid + 1
    xmin[xmin < 0] = 0
    xmax[xmax >= sz_x] = sz_x

    ordsat = np.zeros((sz_y, sz_o), dtype=int)
    for o in range(sz_o):
        for y in range(sz_y):
            ordsat[y,o] = (xmax[y,o] > xmin[y,o]) & np.any(satmask[y,xmin[y,o]:xmax[y,o]] == 1)

    return ordsat


def search_for_saturation_edge(a, x, y, sy, dx, satdown, satlevel, mask):
    sx = dx
    localx = a[x+sx,y+sy]
    while True:
        mask[x+sx,y+sy] = True
        sx += dx
        if x+sx > a.shape[0]-1 or x+sx < 0:
            break
        if a[x+sx,y+sy] >= localx/satdown and a[x+sx,y+sy]<satlevel:
            break
        localx = a[x+sx,y+sy]
    return mask


def determine_saturation_region(a, x, y, sy, dy, satdown, satlevel, mask):
    localy = a[x,y+sy]
    while True:
        mask[x,y+sy] = True
        mask = search_for_saturation_edge(a, x, y, sy, 1, satdown, satlevel, mask)
        mask = search_for_saturation_edge(a, x, y, sy, -1, satdown, satlevel, mask)

        sy += dy
        if y+sy > a.shape[1]-1 or y+sy < 0:
            return mask
        if a[x,y+sy] >= localy/satdown and a[x,y+sy] < satlevel:
            return mask
        localy = a[x,y+sy]


def saturation_mask(a, satlevel):
    """
    ... todo::
        Document this!
    """
    mask = np.zeros(a.shape, dtype=bool)
    a_is_saturated = a >= satlevel
    if not np.any(a_is_saturated):
        return mask.astype(int)

    satdown = 1.001
    sz_x, sz_y = a.shape

    for y in range (0,sz_y):
        for x in range(0,sz_x):
            if a_is_saturated[x,y] and not mask[x,y]:
                mask[x,y] = True
                mask = determine_saturation_region(a, x, y, 0, 1, satdown, satlevel, mask)
                mask = determine_saturation_region(a, x, y, -1, -1, satdown, satlevel, mask)

    return mask.astype(int)

