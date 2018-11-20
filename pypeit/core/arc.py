from __future__ import (print_function, absolute_import, division, unicode_literals)

import inspect

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt



from scipy.ndimage.filters import gaussian_filter

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


def fit_double_poly(all_wv_order, work2d, thismask, nspec_coeff, norder_coeff):
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


def fit2darc(all_wv,all_pix,all_orders,nspec, nspec_coeff=4,norder_coeff=4,sigrej=3.0,debug=True, skip_QA=False):
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
    mnx = 0 #np.min(all_pix)
    mxx = float(nspec - 1) #np.max(all_pix)
    norm_pixel = np.array([0.5 * (mnx + mxx), mxx - mnx])
    pix_nrm = 2. * (all_pix - norm_pixel[0])/norm_pixel[1]
    # Normalize orders
    mnx, mxx = np.min(all_orders), np.max(all_orders)
    norm_order = np.array([0.5 * (mnx + mxx), mxx - mnx])
    orders_nrm = 2. * (all_orders - norm_order[0])/norm_order[1]

    if debug:
        # set some plotting parameters
        utils.pyplot_rcparams()
        plt.figure(figsize=(7,5))
        msgs.info("Plot identified lines")
        cm = plt.cm.get_cmap('RdYlBu_r')
        sc = plt.scatter(orders_nrm, pix_nrm,c=all_wv/10000., cmap=cm)
        cbar = plt.colorbar(sc)
        cbar.set_label(r'Wavelength [$\mu$m]', rotation=270,
                       labelpad=20)
        plt.xlabel(r'Normalized Orders')
        plt.ylabel(r'Normalized Pixels')
        plt.title(r'Location of the identified lines')
        plt.show()


    # Setup some things for the fits
    all_wv_order = all_wv * all_orders
    work2d = np.zeros((nspec_coeff*norder_coeff, len(all_wv)), dtype=np.float64)
    worky = pydl.flegendre(pix_nrm, nspec_coeff)
    workt = pydl.flegendre(orders_nrm, norder_coeff)
    for i in range(norder_coeff):
        for j in range(nspec_coeff):
            work2d[j*norder_coeff+i,:] = worky[j,:] * workt[i,:]

    #ToDO add upper lower to inputs
    lower = np.abs(sigrej)
    upper = np.abs(sigrej)
    maxiter = 25
    iIter = 0
    qdone = False
    thismask = np.ones_like(all_wv,dtype=bool)
    while (not qdone) and (iIter < maxiter):
        coeffs, wv_order_mod = fit_double_poly(all_wv_order, work2d, thismask.astype(float), nspec_coeff, norder_coeff)
        thismask, qdone = pydl.djs_reject(all_wv_order, wv_order_mod, outmask=thismask,
                                          lower=np.float64(lower),upper=np.float64(upper),use_mad=True,sticky=True)
        iIter += 1
        if debug:
            utils.pyplot_rcparams()
            plt.figure(figsize=(7,5))
            plt.axhline(y=np.average(wv_order_mod[thismask] /all_orders[thismask] - all_wv[thismask]),color='r', linestyle='--')
            plt.axhline(y=+np.std(wv_order_mod[thismask] /all_orders[thismask] - all_wv[thismask]),color='r', linestyle=':')
            plt.axhline(y=-np.std(wv_order_mod[thismask] /all_orders[thismask] - all_wv[thismask]),color='r', linestyle=':')
            plt.scatter(all_wv[~thismask]/10000.,wv_order_mod[~thismask] /all_orders[~thismask] - all_wv[~thismask],marker="v",label=r'Rejected values')
            plt.scatter(all_wv[thismask]/10000.,wv_order_mod[thismask]/all_orders[thismask] - all_wv[thismask],marker="v",label=r'Good values')
            plt.text(np.min(all_wv/10000), np.average(wv_order_mod[thismask]/all_orders[thismask]-all_wv[thismask]),
                     r'Average={0:.1f}$\AA$'.format(np.average(wv_order_mod[thismask]/all_orders[thismask]-all_wv[thismask])),
                     ha="left", va="bottom",bbox=dict(boxstyle="square",ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8),alpha=0.7,))
            plt.text(np.max(all_wv/10000), np.std(wv_order_mod[thismask]/all_orders[thismask]-all_wv[thismask]),
                     r'Sigma={0:.1f}$\AA$'.format(np.std(wv_order_mod[thismask]/all_orders[thismask]-all_wv[thismask])),ha="right",
                     va="bottom",bbox=dict(boxstyle="square",ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8),alpha=0.7,))
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
    resid = (wv_order_mod[thismask]-all_wv_order[thismask])
    fin_rms = np.sqrt(np.mean(resid**2))
    msgs.info("RMS: {0:.5f} Ang*Order#".format(fin_rms))

    orders = np.unique(all_orders)
    fit_dict = dict(coeffs=coeffs, orders=orders,
                    nspec_coeff=nspec_coeff, norder_coeff=norder_coeff,
                    pixel_cen=norm_pixel[0], pixel_norm=norm_pixel[1],
                    order_cen=norm_order[0], order_norm=norm_order[1],
                    nspec=nspec)

    debug = True
    if debug:
        fit2darc_qa(fit_dict, fin_rms, all_wv, all_pix, all_orders, thismask)

    return fit_dict

def fit2darc_qa(fit_dict, fin_rms, all_wv, all_pix, all_orders, thismask, setup=None, out_dir=None):
    """ QA on 2D fit of the wavelength solution.
    
    Parameters
    ----------
    fit_dict: dict
      dict of the 2D arc solution
    fin_rms: float
      rms of the model
    all_wv: np.array
     wavelength of the identified lines
    all_pix: np.array
      y-centroid position of the identified lines
    all_orders: np.array
      order number of the identified lines
    thismas: boolean
      mask containing the rejected lines
    setup:
      parameter for QA
    out_dir:
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

    # Define pixels array
    all_pixels_qa = np.arange(nspec)

    # Global plot
    '''
    outfile_global = qa.set_qa_filename(setup, method='fit2darc_global',
                                        out_dir=out_dir)
    '''

    # Define figure properties
    plt.figure(figsize=(8,5))
    plt.title(r'Arc 2D FIT, norder_coeff={:d}, nspec_coeff={:d}, RMS={:5.3f} Ang*Order#'.format(
              norder_coeff, nspec_coeff, fin_rms))
    plt.xlabel(r'Wavelength [$\AA$]')
    plt.ylabel(r'Row [pixel]')

    # Variable where to store the max wavelength covered by the
    # spectrum
    mx = 0.

    # Loop over orders
    for ii in orders:

        # define the color
        rr = (ii-np.max(orders))/(np.min(orders)-np.max(orders))
        gg = 0.0
        bb = (ii-np.min(orders))/(np.max(orders)-np.min(orders))

        # evaluate solution
        wv_order_mod_qa = eval2dfit(fit_dict, all_pixels_qa, ii)

        # Plot solution
        plt.plot(wv_order_mod_qa/ii, all_pixels_qa,color=(rr,gg,bb),
                 linestyle='-', linewidth=2.5)

        # Evaluate residuals at each order
        this_pix = all_pix[all_orders == ii]
        this_wv = all_wv[all_orders == ii]
        this_msk = thismask[all_orders == ii]
        wv_order_mod_resid_qa = eval2dfit(fit_dict, this_pix, ii)
        resid_qa = (wv_order_mod_resid_qa/ii-this_wv)
        plt.scatter((wv_order_mod_resid_qa[~this_msk]/ii)+ \
                    100.*resid_qa[~this_msk], this_pix[~this_msk], \
                    marker='x', color='black', linewidths=2.5, s=16.)
        plt.scatter((wv_order_mod_resid_qa[this_msk]/ii)+ \
                    100.*resid_qa[this_msk], this_pix[this_msk], \
                    color=(rr,gg,bb), linewidth=2.5, s=16.)

        if np.max(wv_order_mod_resid_qa/ii) > mx :
            mx = np.max(wv_order_mod_resid_qa/ii)

    plt.text(mx,np.max(all_pixels_qa),r'residuals $\times$100', \
             ha="right", va="top")

    # Finish
    # plt.savefig(outfile_global, dpi=800)
    # plt.close()
    plt.show()

    # Individual plots
    '''
    outfile_order = qa.set_qa_filename(setup, method='fit2darc_order',
                                       out_dir=out_dir)
    '''

    # set the size of the plot
    nrow = np.int(2)
    ncol = np.int(np.ceil(len(orders)/2.))

    fig = plt.figure(figsize=(5*ncol,6*nrow))

    outer = gridspec.GridSpec(nrow, ncol, wspace=0.3, hspace=0.2)

    for ii_row in range(nrow):
        for ii_col in range(ncol):
            if (ii_row*ncol + ii_col) < len(orders):
                inner = gridspec.GridSpecFromSubplotSpec(2, 1,
                        height_ratios=[2,1], width_ratios=[1],
                        subplot_spec=outer[ii_row*ncol + ii_col],
                        wspace=0.1, hspace=0.0)
                ax0 = plt.Subplot(fig, inner[0])
                ax1 = plt.Subplot(fig, inner[1], sharex=ax0)
                plt.setp(ax0.get_xticklabels(), visible=False)

                ii = orders[ii_row*ncol + ii_col]

                # define the color
                rr = (ii-np.max(orders))/(np.min(orders)-np.max(orders))
                gg = 0.0
                bb = (ii-np.min(orders))/(np.max(orders)-np.min(orders))

                # Evaluate function
                wv_order_mod_qa = eval2dfit(fit_dict, all_pixels_qa, ii)
                # Evaluate delta lambda
                dwl=(wv_order_mod_qa[-1]-wv_order_mod_qa[0])/ii/(all_pixels_qa[-1]-all_pixels_qa[0])

                # Estimate the residuals
                this_pix = all_pix[all_orders == ii]
                this_wv = all_wv[all_orders == ii]
                this_msk = thismask[all_orders == ii]

                wv_order_mod_resid_qa = eval2dfit(fit_dict, this_pix, ii)
                resid_wl_qa = (wv_order_mod_resid_qa/ii-this_wv)

                # Plot the fit
                ax0.set_title('Order = {0:0.0f}'.format(ii))
                ax0.plot(all_pixels_qa, wv_order_mod_qa/ii/10000.,color=(rr,gg,bb), linestyle='-',
                         linewidth=2.5)
                ax0.scatter(this_pix[~this_msk], (wv_order_mod_resid_qa[~this_msk]/ii/10000.)+ \
                            100.*resid_wl_qa[~this_msk]/10000., marker='x', color='black', \
                            linewidth=2.5, s=16.)
                ax0.scatter(this_pix[this_msk], (wv_order_mod_resid_qa[this_msk]/ii/10000.)+ \
                            100.*resid_wl_qa[this_msk]/10000., color=(rr,gg,bb), \
                            linewidth=2.5, s=16.)

                ax0.set_ylabel(r'Wavelength [$\mu$m]')


                # Plot the residuals
                ax1.scatter(this_pix[~this_msk],(resid_wl_qa[~this_msk]/dwl),marker='x', color='black', \
                            linewidth=2.5, s=16.)
                ax1.scatter(this_pix[this_msk], (resid_wl_qa[this_msk]/dwl), color=(rr,gg,bb), \
                            linewidth=2.5, s=16.)
                ax1.axhline(y=0., color=(rr,gg,bb), linestyle=':', linewidth=2.5)
                ax1.get_yaxis().set_label_coords(-0.15,0.5)

                rms_qa = np.sqrt(np.mean((resid_wl_qa[this_msk])**2))

                ax1.set_ylabel(r'Res. [pix]')

                ax0.text(0.1,0.9,r'RMS={0:.3f} Pixel'.format(rms_qa/np.abs(dwl)),ha="left", va="top",
                         transform = ax0.transAxes)
                ax0.text(0.1,0.8,r'$\Delta\lambda$={0:.3f} Pixel/$\AA$'.format(np.abs(dwl)),ha="left", va="top",
                         transform = ax0.transAxes)
                ax0.get_yaxis().set_label_coords(-0.15,0.5)

                fig.add_subplot(ax0)
                fig.add_subplot(ax1)

    fig.text(0.5, 0.04, r'Row [pixel]', ha='center', size='large')
    fig.suptitle(r'Arc 2D FIT, norder_coeff={:d}, nspec_coeff={:d}, RMS={:5.3f} Ang*Order#, residuals $\times$100'.format(norder_coeff, nspec_coeff,fin_rms))

    # # Finish
    # plt.savefig(outfile_order, dpi=800)
    # plt.close()
    plt.show()

def eval2dfit(fit_dict, pixels, order):
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
    order_nrm = 2.0*(osub - fit_dict['order_cen'])/fit_dict['order_norm']
    work_order = pydl.flegendre(order_nrm, norder_coeff)
    # legendre basis for the spectral direction
    pix_nrm = 2.0*(pixels - fit_dict['pixel_cen'])/fit_dict['pixel_norm']
    work_pix = pydl.flegendre(pix_nrm, nspec_coeff)
    # array to hold the
    work2d = np.zeros((nspec_coeff*norder_coeff, npix), dtype=float)
    for i in range(norder_coeff):
        for j in range(nspec_coeff):
            work2d[j*norder_coeff + i, :] = work_pix[j, :]*work_order[i, :]
    wv_order_mod = coeffs.dot(work2d)

    return wv_order_mod

def get_censpec(slit_left, slit_righ, slitpix, arcimg, inmask = None, box_rad = 3.0, xfrac = 0.5):

    """Fit a non-parametric object profile to an object spectrum, unless the S/N ratio is low (> sn_gauss) in which
    fit a simple Gaussian. Port of IDL LOWREDUX long_gprofile.pro


    Parameters
    ----------
    slit_left:  float ndarray
        Left boundary of slit/order to be extracted (given as floating pt pixels). This a 1-d array with shape (nspec, 1)
        or (nspec)

    slit_righ:  float ndarray
        Left boundary of slit/order to be extracted (given as floating pt pixels). This a 1-d array with shape (nspec, 1)
        or (nspec)


    slitpix:  float  ndarray
        Mask image specifying the pixels which lie on the slit/order to search for objects on. This is created by
        traceslits in the tslits_dict, and has the convention that each slit/order has an integer index starting with one.
        Locations with zeros are not on slits/orders.

    arcimg:  float ndarray
        Image to extract the arc from. This should be an arcimage or perhaps a frame with night sky lines.

    Optional Parameters
    ----------
    inmask: boolean ndararay
         Input mask image with same shape as arcimg. Convention True = good and False = bad. The default is None.

    box_rad: float, [default = 3.0]
        Size of boxcar window in pixels (as a floating point number) in the spatial direction used to extract the arc.

    xfrac: float [default = 0.5]
        Fraction location along the slit to extract the arc. The default is at the midpoint of the slit which is 0.5,
        but this can be adjusted to an off center location.

    Returns
    -------
    :func:`tuple`
         A tuple containing the (arc_spec, maskslit)

            arc_spec: float ndarray with shape (nspec, nslits)
                Array containing the extracted arc spectrum for each slit.

            maskslit: int ndarray with shape (nslits)
               output mask indicating whether a slit is good or bad. 0 = good, 1 = bad.
     """

    if inmask is None:
        inmask = slitpix > -1

    nslits = slit_left.shape[1]
    (nspec, nspat) = arcimg.shape
    maskslit = np.zeros(nslits, dtype=np.int)
    trace = slit_left + xfrac*(slit_righ - slit_left)
    arc_spec = np.zeros((nspec, nslits))
    spat_img = np.outer(np.ones(nspec,dtype=int), np.arange(nspat,dtype=int)) # spatial position everywhere along image


    for islit in range(nslits):
        msgs.info("Extracting an approximate arc spectrum at the centre of slit {:d}".format(islit + 1))
        # Create a mask for the pixels that will contribue to the arc
        trace_img = np.outer(trace[:,islit], np.ones(nspat))  # left slit boundary replicated spatially
        arcmask = (slitpix > -1) & inmask & (spat_img > (trace_img - box_rad)) & (spat_img < (trace_img + box_rad))
        # Trimming the image makes this much faster
        left = np.fmax(spat_img[arcmask].min() - 4,0)
        righ = np.fmin(spat_img[arcmask].max() + 5,nspat)
        this_mean, this_med, this_sig = sigma_clipped_stats(arcimg[:,left:righ], mask=np.invert(arcmask[:,left:righ])
                                                            , sigma=3.0, axis=1)
        arc_spec[:,islit] = this_med.data
        if not np.any(arc_spec[:,islit]):
            maskslit[islit] = 1

    return arc_spec, maskslit



def detect_peaks(x, mph=None, mpd=1, threshold=0, edge='rising',
                 kpsh=False, valley=False, show=False, ax=None):
    """Detect peaks in data based on their amplitude and other features.

    Parameters
    ----------
    x : 1D array_like
        data.
    mph : {None, number}, optional (default = None)
        detect peaks that are greater than minimum peak height (if parameter
        `valley` is False) or peaks that are smaller than maximum peak height
         (if parameter `valley` is True).
    mpd : positive integer, optional (default = 1)
        detect peaks that are at least separated by minimum peak distance (in
        number of data).
    threshold : positive number, optional (default = 0)
        detect peaks (valleys) that are greater (smaller) than `threshold`
        in relation to their immediate neighbors.
    edge : {None, 'rising', 'falling', 'both'}, optional (default = 'rising')
        for a flat peak, keep only the rising edge ('rising'), only the
        falling edge ('falling'), both edges ('both'), or don't detect a
        flat peak (None).
    kpsh : bool, optional (default = False)
        keep peaks with same height even if they are closer than `mpd`.
    valley : bool, optional (default = False)
        if True (1), detect valleys (local minima) instead of peaks.
    show : bool, optional (default = False)
        if True (1), plot data in matplotlib figure.
    ax : a matplotlib.axes.Axes instance, optional (default = None).

    Returns
    -------
    ind : 1D array_like
        indeces of the peaks in `x`.

    Notes
    -----
    The detection of valleys instead of peaks is performed internally by simply
    negating the data: `ind_valleys = detect_peaks(-x)`

    The function can handle NaN's

    See this IPython Notebook [1]_.

    __author__ = "Marcos Duarte, https://github.com/demotu/BMC"
    __version__ = "1.0.5"
    __license__ = "MIT"

    References
    ----------
    .. [1] http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb

    Examples
    --------
    >>> from detect_peaks import detect_peaks
    >>> x = np.random.randn(100)
    >>> x[60:81] = np.nan
    >>> # detect all peaks and plot data
    >>> ind = detect_peaks(x, show=True)
    >>> print(ind)

    >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
    >>> # set minimum peak height = 0 and minimum peak distance = 20
    >>> detect_peaks(x, mph=0, mpd=20, show=True)

    >>> x = [0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0]
    >>> # set minimum peak distance = 2
    >>> detect_peaks(x, mpd=2, show=True)

    >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
    >>> # detection of valleys instead of peaks
    >>> detect_peaks(x, mph=-1.2, mpd=20, valley=True, show=True)

    >>> x = [0, 1, 1, 0, 1, 1, 0]
    >>> # detect both edges
    >>> detect_peaks(x, edge='both', show=True)

    >>> x = [-2, 1, -2, 2, 1, 1, 3, 0]
    >>> # set threshold = 2
    >>> detect_peaks(x, threshold = 2, show=True)

    Version history
    ---------------
    '1.0.5':
        The sign of `mph` is inverted if parameter `valley` is True

    """

    x = np.atleast_1d(x).astype('float64')
    if x.size < 3:
        return np.array([], dtype=int)
    if valley:
        x = -x
        if mph is not None:
            mph = -mph
    # find indices of all peaks
    dx = x[1:] - x[:-1]
    # handle NaN's
    indnan = np.where(np.isnan(x))[0]
    if indnan.size:
        x[indnan] = np.inf
        dx[np.where(np.isnan(dx))[0]] = np.inf
    ine, ire, ife = np.array([[], [], []], dtype=int)
    if not edge:
        ine = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) > 0))[0]
    else:
        if edge.lower() in ['rising', 'both']:
            ire = np.where((np.hstack((dx, 0)) <= 0) & (np.hstack((0, dx)) > 0))[0]
        if edge.lower() in ['falling', 'both']:
            ife = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) >= 0))[0]
    ind = np.unique(np.hstack((ine, ire, ife)))
    # handle NaN's
    if ind.size and indnan.size:
        # NaN's and values close to NaN's cannot be peaks
        ind = ind[np.in1d(ind, np.unique(np.hstack((indnan, indnan - 1, indnan + 1))), invert=True)]
    # first and last values of x cannot be peaks
    if ind.size and ind[0] == 0:
        ind = ind[1:]
    if ind.size and ind[-1] == x.size - 1:
        ind = ind[:-1]
    # remove peaks < minimum peak height
    if ind.size and mph is not None:
        ind = ind[x[ind] >= mph]
    # remove peaks - neighbors < threshold
    if ind.size and threshold > 0:
        dx = np.min(np.vstack([x[ind] - x[ind - 1], x[ind] - x[ind + 1]]), axis=0)
        ind = np.delete(ind, np.where(dx < threshold)[0])
    # detect small peaks closer than minimum peak distance
    if ind.size and mpd > 1:
        ind = ind[np.argsort(x[ind])][::-1]  # sort ind by peak height
        idel = np.zeros(ind.size, dtype=bool)
        for i in range(ind.size):
            if not idel[i]:
                # keep peaks with the same height if kpsh is True
                idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
                       & (x[ind[i]] > x[ind] if kpsh else True)
                idel[i] = 0  # Keep current peak
        # remove the small peaks and sort back the indices by their occurrence
        ind = np.sort(ind[~idel])

    if show:
        if indnan.size:
            x[indnan] = np.nan
        if valley:
            x = -x
            if mph is not None:
                mph = -mph
        _plot(x, mph, mpd, threshold, edge, valley, ax, ind)

    return ind


def _plot(x, mph, mpd, threshold, edge, valley, ax, ind):
    """Plot results of the detect_peaks function, see its help."""


    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(8, 4))

    ax.plot(x, 'b', lw=1)
    if ind.size:
        label = 'valley' if valley else 'peak'
        label = label + 's' if ind.size > 1 else label
        ax.plot(ind, x[ind], '+', mfc=None, mec='r', mew=2, ms=8,
                label='%d %s' % (ind.size, label))
        ax.legend(loc='best', framealpha=.5, numpoints=1)
    ax.set_xlim(-.02 * x.size, x.size * 1.02 - 1)
    ymin, ymax = x[np.isfinite(x)].min(), x[np.isfinite(x)].max()
    yrange = ymax - ymin if ymax > ymin else 1
    ax.set_ylim(ymin - 0.1 * yrange, ymax + 0.1 * yrange)
    ax.set_xlabel('Data #', fontsize=14)
    ax.set_ylabel('Amplitude', fontsize=14)
    mode = 'Valley detection' if valley else 'Peak detection'
    ax.set_title("%s (mph=%s, mpd=%f, threshold=%s, edge='%s')"
                 % (mode, str(mph), mpd, str(threshold), edge))
    # plt.grid()
    plt.show()

# ToDO JFH nfitpix should be chosen based on the spectral sampling of the spectroscopic setup
def detect_lines(censpec, sigdetect = 5.0, input_thresh = None, cont_subtract=True, fwhm = 4.0, fit_frac_fwhm=1.25, mask_frac_fwhm = 1.0, max_frac_fwhm = 2.5, cont_samp = 30,
                 nonlinear_counts=1e10, niter_cont = 3,nfind = None, verbose = False, debug=False, debug_peak_find = False):
    """
    Extract an arc down the center of the chip and identify
    statistically significant lines for analysis.

    Parameters
    ----------
    censpec : ndarray
      A 1D spectrum to be searched for significant detections

    Optional Parameters
    -------------------
    sigdetect: float, default=20.
       Sigma threshold above fluctuations for arc-line detection. Arcs are continuum subtracted and the fluctuations are
       computed after continuum subtraction.

    input_thresh: float, str, default= None
       Optionally the user can specify the threhsold that peaks must be above to be kept. In this case the sigdetect parameter
       will be ignored. This is most useful for example for cases where cont_subtract =False, and the user prefers to determine
       the significance  threhsold outside of this routine, rather than using this routines defaults to determine the continuum level
       and standard deviation of the continuum subtracted spetrum. If a string input of 'None' is set then the code will simply return
        all peaks irrespective of any threshold. This is equivalent to setting the mph parameter to None in the detect_peaks code.

    fwhm:  float, default = 4.0
       Number of pixels per fwhm resolution element.

    fit_frac_fwhm: float, default 0.5
       Number of pixels that are used in the fits for Gaussian arc line centroiding expressed as a fraction of the fwhm parameter

    max_frac_fwhm:  float, default = 2.5
       maximum width allowed for usable arc lines expressed relative to the fwhm.

    mask_frac_fwhm float, default = 1.0
       width used for masking peaks in the spectrum when the continuum is being defined. Expressed as a fraction of the fwhm
       parameter

    cont_subtract: bool, default = True
       If true, the code will continuum subtract the input array by iteratively determining the continuum

    cont_samp: float, default = 30.0
       The number of samples across the spectrum used for continuum subtraction. Continuum subtraction is done via
       median filtering, with a width of ngood/cont_samp, where ngood is the number of good pixels for estimating the continuum
       (i.e. that don't have peaks).

    niter_cont: int, default = 3
       Number of iterations of peak finding, masking, and continuum fitting used to define the continuum.


    nonlinear_counts: float, default = 1e10
       Value above which to mask saturated arc lines. This should be nonlinear_counts= nonlinear*saturation according to pypeit parsets.
       Default is 1e10 which is to not mask.

    nfind: int, default = None
       Return only the nfind highest significance lines. The default is None, which means the code will
       return all the lines above the significance threshold.

    verbose: bool, default = False
       Output more stuff to the screen.

    debug: boolean, default = False
       Make plots showing results of peak finding and final arc lines that are used.

    Returns
    -------
    tampl : ndarray
      The amplitudes of the line detections in the true arc
    tampl_cont : ndarray
      The amplitudes of the line detections in the continuum subtracted arc
    tcent : ndarray
      The centroids of the line detections
    twid : ndarray
      The 1sigma Gaussian widths of the line detections
    centerr : ndarray
      The variance on tcent
    w : ndarray
      An index array indicating which detections are the most reliable.
    arc : ndarray
      The continuum sutracted arc used to find detections.
    nsig : ndarray
      The significance of each line detected relative to the 1sigma variation in the continuum subtracted arc in the
      the line free region. Bad lines are assigned a significance of -1, since they don't have an amplitude fit
    """

    # debug = True

    # Detect the location of the arc lines
    if verbose:
        msgs.info("Detecting lines...isolating the strongest, nonsaturated lines")

    if len(censpec.shape) == 3:
        detns = censpec[:, 0].flatten()
    else:
        detns = censpec.copy()
    detns = detns.astype(np.float)
    xrng = np.arange(detns.size, dtype=np.float)

    cont_mask = np.ones(detns.size, dtype=bool)
    if cont_subtract:
        nspec = detns.size
        spec_vec = np.arange(nspec)
        cont_now = np.arange(nspec)
        mask_sm = np.round(mask_frac_fwhm*fwhm).astype(int)
        mask_odd = mask_sm + 1 if mask_sm % 2 == 0 else mask_sm
        for iter in range(niter_cont):
            arc_now = detns - cont_now
            (mean, med, stddev) = sigma_clipped_stats(arc_now[cont_mask], sigma_lower=3.0, sigma_upper=3.0)
            # be very liberal in determining threshold for continuum determination
            thresh = med + 2.0*stddev
            pixt_now = detect_peaks(arc_now, mph=thresh, mpd=fwhm*0.75)
            # mask out the peaks we find for the next continuum iteration
            cont_mask_fine = np.ones_like(cont_now)
            cont_mask_fine[pixt_now] = 0.0
            # cont_mask is the mask for defining the continuum regions: True is good,  False is bad
            cont_mask = (utils.smooth(cont_mask_fine,mask_odd) > 0.999)
            # If more than half the spectrum is getting masked than short circuit this masking
            frac_mask = np.sum(~cont_mask)/float(nspec)
            if (frac_mask > 0.70):
                msgs.warn('Too many pixels masked in arc continuum definiton: frac_mask = {:5.3f}'.format(frac_mask) + ' . Not masking....')
                cont_mask = np.ones_like(cont_mask)
            ngood = np.sum(cont_mask)
            samp_width = np.ceil(ngood/cont_samp).astype(int)
            cont_med = utils.fast_running_median(detns[cont_mask], samp_width)
            cont_now = np.interp(spec_vec,spec_vec[cont_mask],cont_med)
    else:
        cont_now = np.zeros_like(detns)

    #msgs.info('Detected {:d} lines in the arc spectrum.'.format(len(w[0])))
    # Final peak detection

    arc = detns - cont_now
    if input_thresh is None:
        (mean, med, stddev) = sigma_clipped_stats(arc[cont_mask], sigma_lower=3.0, sigma_upper=3.0)
        thresh = med + sigdetect*stddev
    else:
        med = 0.0
        if isinstance(input_thresh,(float, int)):
            thresh = input_thresh
        elif isinstance(input_thresh, str):
            if input_thresh == 'None':
                thresh = None
        else:
            msgs.error('Unrecognized value for thresh')
        stddev = 1.0
    pixt = detect_peaks(arc, mph=thresh, mpd=fwhm*0.75, show=debug_peak_find)
    nfitpix = np.round(fit_frac_fwhm*fwhm).astype(int)
    fwhm_max = max_frac_fwhm*fwhm
    tampl_fit, tcent, twid, centerr = fit_arcspec(xrng, arc, pixt, nfitpix)
    # This is the amplitude of the lines in the actual detns spectrum not continuum subtracted
    tampl_true = np.interp(pixt, xrng, detns)
    tampl = np.interp(pixt, xrng, arc)
    #         width is fine & width > 0.0 & width < FWHM/2.35 &  center positive  &  center on detector
    #        & amplitude not nonlinear
    good = (np.invert(np.isnan(twid))) & (twid > 0.0) & (twid < fwhm_max/2.35) & (tcent > 0.0) & (tcent < xrng[-1]) & \
           (tampl_true < nonlinear_counts)
    ww = np.where(good)
    # Compute the significance of each line, set the significance of bad lines to be -1
    nsig = (tampl - med)/stddev

    # If the user requested the nfind most significant peaks have been requested, then grab and return only these lines
    if nfind is not None:
        if nfind > len(nsig):
            msgs.warn('Requested nfind = {:}'.format(nfind) + ' peaks but only npeak = {:}'.format(len(tampl)) +
                      ' were found. Returning all the peaks found.')
        else:
            ikeep = (nsig.argsort()[::-1])[0:nfind]
            tampl_true = tampl_true[ikeep]
            tampl = tampl[ikeep]
            tcent = tcent[ikeep]
            twid = twid[ikeep]
            centerr = centerr[ikeep]
            ww = np.where(good[ikeep])
            nsig = nsig[ikeep]
            good = good[ikeep]

    if debug:
        plt.figure(figsize=(14, 6))
        plt.plot(xrng, arc, color='black', drawstyle = 'steps-mid', lw=3, label = 'arc', linewidth = 1.0)
        plt.plot(tcent[~good], tampl[~good],'r+', markersize =6.0, label = 'bad peaks')
        plt.plot(tcent[good], tampl[good],'g+', markersize =6.0, label = 'good peaks')
        if thresh is not None:
            plt.hlines(thresh, xrng.min(), xrng.max(), color='cornflowerblue', linestyle=':', linewidth=2.0,
                       label='threshold', zorder=10)
        if nonlinear_counts < 1e9:
            plt.hlines(nonlinear_counts,xrng.min(), xrng.max(), color='orange', linestyle='--',linewidth=2.0,
                       label='nonlinear', zorder=10)
        plt.title('Good Lines = {:d}'.format(np.sum(good)) + ',  Bad Lines = {:d}'.format(np.sum(~good)))
        plt.ylim(arc.min(), 1.5*arc.max())
        plt.legend()
        plt.show()

    return tampl_true, tampl, tcent, twid, centerr, ww, arc, nsig


def fit_arcspec(xarray, yarray, pixt, fitp):
    """
    Fit an arc spectrum line

    Args:
        xarray:
        yarray:
        pixt:
        fitp: int
          Number of pixels to fit with

    Returns:

    """

    fitp_even = fitp if fitp % 2 == 0 else fitp + 1
    fit_interval = fitp_even//2

    # Setup the arrays with fit parameters
    sz_p = pixt.size
    sz_a = yarray.size
    b      = -1.0*np.ones(sz_p, dtype=np.float)
    ampl   = -1.0*np.ones(sz_p, dtype=np.float)
    cent   = -1.0*np.ones(sz_p, dtype=np.float)
    widt   = -1.0*np.ones(sz_p, dtype=np.float)
    centerr = -1.0*np.ones(sz_p, dtype=np.float)

    for p in range(sz_p):
#        pmin = pixt[p]-(fitp-1)//2
#        pmax = pixt[p]-(fitp-1)//2 + fitp
        pmin = pixt[p]-fit_interval
        pmax = pixt[p]-fit_interval + fitp + 1
        if pmin < 0:
            pmin = 0
        if pmax > sz_a:
            pmax = sz_a
        if pmin == pmax:
            continue
        if (pmax - pmin) < (fit_interval):
            continue # Probably won't be a good solution
#       JFH Removed below
#       if pixt[p]-pmin <= 1 or pmax-pixt[p] <= 1:
#            continue  # Probably won't be a good solution
        # Fit the gaussian
        try:
            popt, pcov = utils.func_fit(xarray[pmin:pmax], yarray[pmin:pmax], "gaussian", 3, return_errors=True)
            ampl[p] = popt[0]
            cent[p] = popt[1]
            widt[p] = popt[2]
            centerr[p] = pcov[1, 1]
            #popt, pcov = utils.func_fit(xarray[pmin:pmax], yarray[pmin:pmax], "gaussian", 4, return_errors=True)
            #b[p]    = popt[0]
            #ampl[p] = popt[1]
            #cent[p] = popt[2]
            #widt[p] = popt[3]
            #centerr[p] = pcov[2, 2]
        except RuntimeError:
            pass
    return ampl, cent, widt, centerr


def simple_calib_driver(msarc, aparm, censpec, ok_mask, nfitpix=5, get_poly=False,
                        IDpixels=None, IDwaves=None):
    wv_calib = {}
    for slit in ok_mask:
        iwv_calib = simple_calib(msarc, aparm, censpec[:, slit], nfitpix=nfitpix,
                                 get_poly=get_poly, IDpixels=IDpixels, IDwaves=IDwaves)
        wv_calib[str(slit)] = iwv_calib.copy()
    return wv_calib


def simple_calib(msarc, aparm, censpec, nfitpix=5, get_poly=False,
                 IDpixels=None, IDwaves=None, debug=False, sigdetect=5.):
    """Simple calibration algorithm for longslit wavelengths

    Uses slf._arcparam to guide the analysis

    Parameters
    ----------
    msarc : ndarray
    aparm : dict
    censpec : ndarray
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
    tampl, tcent, twid, _, w, yprep, nsig = detect_lines(censpec, nfitpix=nfitpix,
                                                         sigdetect=sigdetect,
                                                         nonlinear_counts = aparm['nonlinear_counts'])

    # Cut down to the good ones
    tcent = tcent[w]
    tampl = tampl[w]
    msgs.info('Detected {:d} lines in the arc spectrum.'.format(len(w[0])))

    # Read Arc linelist
    llist = aparm['llist']

    # IDs were input by hand
    if IDpixels is not None:
        # Check that there are at least 5 values
        pixels = np.array(IDpixels) # settings.argflag['arc']['calibrate']['IDpixels'])
        if np.sum(pixels > 0.) < 5:
            msgs.error("Need to give at least 5 pixel values!")
        #
        msgs.info("Using input lines to seed the wavelength solution")
        # Calculate median offset
        mdiff = [np.min(np.abs(tcent-pix)) for pix in pixels]
                 #settings.argflag['arc']['calibrate']['IDpixels']]
        med_poff = np.median(np.array(mdiff))
        msgs.info("Will apply a median offset of {:g} pixels".format(med_poff))

        # Match input lines to observed spectrum
        nid = pixels.size # len(settings.argflag['arc']['calibrate']['IDpixels'])
        idx_str = np.ones(nid).astype(int)
        ids = np.zeros(nid)
        idsion = np.array(['     ']*nid)
        gd_str = np.arange(nid).astype(int)
        for jj,pix in enumerate(pixels): #settings.argflag['arc']['calibrate']['IDpixels']):
            diff = np.abs(tcent-pix-med_poff)
            if np.min(diff) > 2.:
                debugger.set_trace()
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
                idsion[jj] = llist['Ion'][imnw]
                msgs.info("Identifying arc line: {:s} {:g}".format(idsion[jj],ids[jj]))
    else:
        # Generate dpix pairs
        msgs.info("Using pair algorithm for wavelength solution")
        nlist = len(llist)
        dpix_list = np.zeros((nlist,nlist))
        for kk,row in enumerate(llist):
            #dpix_list[kk,:] = (np.array(row['wave'] - llist['wave']))/disp
            dpix_list[kk,:] = msarc.shape[0]*(aparm['b1']*(np.array(row['wave'] - llist['wave'])) + aparm['b2']*np.array(row['wave']**2 - llist['wave']**2) )

        # Lambda pairs for the strongest N lines
        srt = np.argsort(tampl)
        idx_str = srt[-aparm['Nstrong']:]
        idx_str.sort()
        dpix_obs = np.zeros((aparm['Nstrong'], aparm['Nstrong']))
        for kk,idx in enumerate(idx_str):
            dpix_obs[kk,:] = np.array(tcent[idx] - tcent[idx_str])

        # Match up (ugly loops)
        ids = np.zeros(aparm['Nstrong'])
        idsion = np.array(['     ']*aparm['Nstrong'])
        for kk in range(aparm['Nstrong']):
            med_off = np.zeros(nlist)
            for ss in range(nlist):
                dpix = dpix_list[ss]
                min_off = []
                for jj in range(aparm['Nstrong']):
                    min_off.append(np.min(np.abs(dpix_obs[kk,jj]-dpix)))
                med_off[ss] = np.median(min_off)
            # Set by minimum
            idm = np.argmin(med_off)
            ids[kk] = llist['wave'][idm]
            idsion[kk] = llist['Ion'][idm]

        # Calculate disp of the strong lines
        disp_str = np.zeros(aparm['Nstrong'])
        for kk in range(aparm['Nstrong']):
            disp_val = (ids[kk]-ids)/(tcent[idx_str[kk]]-tcent[idx_str])
            isf = np.isfinite(disp_val)
            disp_str[kk] = np.median(disp_val[isf])
        # Consider calculating the RMS with clipping
        gd_str = np.where( np.abs(disp_str-aparm['disp'])/aparm['disp'] < aparm['disp_toler'])[0]
        msgs.info('Found {:d} lines within the dispersion threshold'.format(len(gd_str)))
        if len(gd_str) < 5:
            if debug:
                msgs.warn('You should probably try your best to ID lines now.')
                debugger.set_trace()
                debugger.plot1d(yprep)
            else:
                msgs.error('Insufficient lines to auto-fit.')

    # Debug
    msgs.work('Cross correlate here?')

    # Setup for fitting
    ifit = idx_str[gd_str]
    sv_ifit = list(ifit) # Keep the originals
    all_ids = -999.*np.ones(len(tcent))
    all_idsion = np.array(['12345']*len(tcent))
    all_ids[ifit] = ids[gd_str]
    all_idsion[ifit] = idsion[gd_str]
    # Fit
    n_order = aparm['n_first']
    flg_quit = False
    fmin, fmax = -1., 1.
    msgs.info('Iterative wavelength fitting..')
    #debug=True
    while (n_order <= aparm['n_final']) and (flg_quit is False):
        #msgs.info('n_order={:d}'.format(n_order))
        # Fit with rejection
        xfit, yfit = tcent[ifit], all_ids[ifit]
        mask, fit = utils.robust_polyfit(xfit, yfit, n_order, function=aparm['func'], sigma=aparm['nsig_rej'], minv=fmin, maxv=fmax)
        rms_ang = utils.calc_fit_rms(xfit, yfit, fit, aparm['func'], minv=fmin, maxv=fmax)
        # Reject but keep originals (until final fit)
        ifit = list(ifit[mask == 0]) + sv_ifit
        # Find new points (should we allow removal of the originals?)
        twave = utils.func_val(fit, tcent, aparm['func'], minv=fmin, maxv=fmax)
        for ss,iwave in enumerate(twave):
            mn = np.min(np.abs(iwave-llist['wave']))
            if mn/aparm['disp'] < aparm['match_toler']:
                imn = np.argmin(np.abs(iwave-llist['wave']))
                if debug:
                    print('Adding {:g} at {:g}'.format(llist['wave'][imn],tcent[ss]))
                # Update and append
                all_ids[ss] = llist['wave'][imn]
                all_idsion[ss] = llist['Ion'][imn]
                ifit.append(ss)
        # Keep unique ones
        ifit = np.unique(np.array(ifit,dtype=int))
        if debug:
            debugger.set_trace()
        # Increment order
        if n_order < aparm['n_final']:
            n_order += 1
        else:
            # This does 2 iterations at the final order
            flg_quit = True

    # Final fit (originals can now be rejected)
    fmin, fmax = 0., 1.
    xfit, yfit = tcent[ifit]/(msarc.shape[0]-1), all_ids[ifit]
    mask, fit = utils.robust_polyfit(xfit, yfit, n_order, function=aparm['func'], sigma=aparm['nsig_rej_final'], minv=fmin, maxv=fmax)
    irej = np.where(mask==1)[0]
    if len(irej) > 0:
        xrej = xfit[irej]
        yrej = yfit[irej]
        for imask in irej:
            msgs.info('Rejecting arc line {:g}'.format(yfit[imask]))
    else:
        xrej = []
        yrej = []
    xfit = xfit[mask==0]
    yfit = yfit[mask==0]
    ions = all_idsion[ifit][mask==0]
    #
    '''
    if debug:
        wave = utils.func_val(fit, np.arange(msarc.shape[0])/float(msarc.shape[0]),
            'legendre', minv=fmin, maxv=fmax)
        debugger.set_trace()

        #debugger.xplot(xfit, np.ones(len(xfit)), scatter=True,
        #    xtwo=np.arange(msarc.shape[0]),ytwo=yprep)
        #debugger.xplot(xfit,yfit, scatter=True, xtwo=np.arange(msarc.shape[0]),
        #    ytwo=wave)
        #debugger.set_trace()
        #wave = utils.func_val(fit, np.arange(msarc.shape[0])/float(msarc.shape[0]),
        #    'legendre', min=fmin, max=fmax)
    '''

    # 2nd order Poly fit for archival
    #get_poly=True
    if get_poly:
        poly_fit = utils.func_fit(yfit,xfit, 'polynomial',2, minv=fmin, maxv=fmax)
        print(' Most likely you with to record these values:')
        print(poly_fit)
        debugger.set_trace()
    # Pack up fit
    final_fit = dict(fitc=fit, function=aparm['func'], xfit=xfit, yfit=yfit,
        ions=ions, fmin=fmin, fmax=fmax, xnorm=float(msarc.shape[0]),
        xrej=xrej, yrej=yrej, mask=mask, spec=yprep, nrej=aparm['nsig_rej_final'],
        shift=0., tcent=tcent)
    # QA
    #arc_fit_qa(slf, final_fit, slit)
    # RMS
    rms_ang = utils.calc_fit_rms(xfit, yfit, fit, aparm['func'], minv=fmin, maxv=fmax)
    wave = utils.func_val(fit, np.arange(msarc.shape[0])/float(msarc.shape[0]),
                            aparm['func'], minv=fmin, maxv=fmax)
    rms_pix = rms_ang/np.median(np.abs(wave-np.roll(wave,1)))
    msgs.info("Fit RMS = {} pix".format(rms_pix))
    # Return
    return final_fit

# JFH I think this is all deprecated code as it is in wavecalib.py now
def calib_with_arclines(aparm, spec, ok_mask=None, use_method="general"):
    """Holy grail algorithms for wavelength calibration

    Uses arcparam to guide the analysis

    Parameters
    ----------
    aparm
    spec
    use_method : str, optional

    Returns
    -------
    final_fit : dict
      Dict of fit info
    """
    raise DeprecationWarning("THIS HAS BEEN MOVED INSIDE wavecalib.py")
    assert False

    if ok_mask is None:
        ok_mask = np.arange(spec.shape[1])

    if use_method == "semi-brute":
        final_fit = {}
        for slit in ok_mask:
            best_dict, ifinal_fit = autoid.semi_brute(spec[:, slit], aparm['lamps'], aparm['wv_cen'], aparm['disp'],
                                                      fit_parm=aparm, min_nsig=aparm['min_nsig'], nonlinear_counts= aparm['nonlinear_counts'])
            final_fit[str(slit)] = ifinal_fit.copy()
    elif use_method == "basic":
        final_fit = {}
        for slit in ok_mask:
            status, ngd_match, match_idx, scores, ifinal_fit =\
                autoid.basic(spec[:, slit], aparm['lamps'], aparm['wv_cen'], aparm['disp'], nonlinear_counts = aparm['nonlinear_counts'])
            final_fit[str(slit)] = ifinal_fit.copy()
    else:
        # Now preferred
        arcfitter = autoid.General(spec, aparm['lamps'], ok_mask=ok_mask, fit_parm=aparm, min_nsig=aparm['min_nsig'],
                                   lowest_nsig=aparm['lowest_nsig'], nonlinear_counts = aparm['nonlinear_counts'],
                                   rms_threshold=aparm['rms_threshold'])
        patt_dict, final_fit = arcfitter.get_results()
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



