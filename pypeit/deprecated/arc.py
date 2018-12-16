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
        resid_wl = (wv_order_mod_resi d /i i -this_wv)
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
