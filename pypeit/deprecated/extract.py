import time
import copy
import inspect

import numpy as np
import scipy


#from matplotlib import gridspec, font_manager

from astropy import stats

from pypeit import msgs
from pypeit.core import pydl
from pypeit import utils
from pypeit.core import pixels
from pypeit import ginga
from matplotlib import pyplot as plt
from pypeit.core import trace_slits
from pypeit.core import arc
from scipy import interpolate

from sklearn.decomposition import PCA
from pypeit import specobjs
#from pypeit import tracepca
from pypeit.core.pydl import spheregroup

from IPython import embed

def extract_boxcar(image,trace_in, radius_in, ycen = None):
    """ Extract the total flux within a boxcar window at many positions. The ycen position is optional. If it is not provied, it is assumed to be integers
     in the spectral direction (as is typical for traces). Traces are expected to run vertically to be consistent with other
     extract_  routines. Based on idlspec2d/spec2d/extract_boxcar.pro

     Parameters
     ----------
     image :  float ndarray
         Image to extract from. It is a 2-d array with shape (nspec, nspat)

     trace_in :  float ndarray
         Trace for the region to be extracted (given as floating pt pixels). This can either be an 2-d  array with shape
         (nspec, nTrace) array, or a 1-d array with shape (nspec) forthe case of a single trace.

     radius :  float or ndarray
         boxcar radius in floating point pixels. This can be either be in put as a scalar or as an array to perform
         boxcar extraction a varaible radius. If an array is input it must have the same size and shape as trace_in, i.e.
         a 2-d  array with shape (nspec, nTrace) array, or a 1-d array with shape (nspec) for the case of a single trace.


     Optional Parameters
     -------------------
     ycen :  float ndarray
         Y positions corresponding to trace_in (expected as integers). Will be rounded to the nearest integer if floats
         are provided. This needs to have the same shape as trace_in  provided above. In other words,
         either a  2-d  array with shape (nspec, nTrace) array, or a 1-d array with shape (nspec) forthe case of a single trace.


     Returns
     -------
     fextract:   ndarray
         Extracted flux at positions specified by (left<-->right, ycen). The output will have the same shape as
         Left and Right, i.e.  an 2-d  array with shape (nspec, nTrace) array if multiple traces were input, or a 1-d array with shape (nspec) for
         the case of a single trace.

     Revision History
     ----------------
     24-Mar-1999  Written by David Schlegel, Princeton.
     22-Apr-2018  Ported to python by Joe Hennawi
     """


    # Checks on radius
    if (isinstance(radius_in,int) or isinstance(radius_in,float)):
        radius = radius_in
    elif ((np.size(radius_in)==np.size(trace_in)) & (np.shape(radius_in) == np.shape(trace_in))):
        radius = radius_in.T
    else:
        raise ValueError('Boxcar radius must a be either an integer, a floating point number, or an ndarray '
                         'with the same shape and size as trace_in')

    trace = trace_in.T

    dim = trace.shape
    ndim = len(dim)
    if (ndim == 1):
        nTrace = 1
        npix = dim[0]
    else:
        nTrace = dim[0]
        npix = dim[1]

    if ycen is None:
        if ndim == 1:
            ycen_out = np.arange(npix, dtype='int')
        elif ndim == 2:
            ycen_out = np.outer(np.ones(nTrace, dtype=int), np.arange(npix, dtype=int))
        else:
            raise ValueError('trace is not 1 or 2 dimensional')
    else:
        ycen_out = ycen.T
        ycen_out = np.rint(ycen_out).astype(int)

    if ((np.size(trace) != np.size(ycen_out)) | (np.shape(trace) != np.shape(ycen_out))):
        raise ValueError('Number of elements and shape of trace and ycen must be equal')



    left = trace - radius
    right = trace + radius
    fextract = extract_asymbox2(image, left, right, ycen_out)

    return fextract

def extract_asymbox2(image,left_in,right_in, ycen=None, weight_image=None):
    """ Extract the total flux within a variable window at many positions. This routine will accept an asymmetric/variable window
    specified by the left_in and right_in traces.  The ycen position is optional. If it is not provied, it is assumed to be integers
    in the spectral direction (as is typical for traces). Traces are expected to run vertically to be consistent with other
    extract_  routines. Based on idlspec2d/spec2d/extract_asymbox2.pro

    Args:
    image :  float ndarray
        Image to extract from. It is a 2-d array with shape (nspec, nspat)
    left  :  float ndarray
        Left boundary of region to be extracted (given as floating pt pixels). This can either be an 2-d  array with shape
        (nspec, nTrace) array, or a 1-d array with shape (nspec) forthe case of a single trace.

    right  :  float ndarray
        Right boundary of region to be extracted (given as floating pt pixels). This can either be an 2-d  array with shape
        (nspec, nTrace) array, or a 1-d array with shape (nspec) forthe case of a single trace.


    Returns:
    ycen :  float ndarray
        Y positions corresponding to "Left"  and "Right" (expected as integers). Will be cast to an integer if floats
        are provided. This needs to have the same shape as left and right broundarys provided above. In other words,
        either a  2-d  array with shape (nspec, nTrace) array, or a 1-d array with shape (nspec) forthe case of a single trace.

    weight_image: float ndarray
        Weight map to be applied to image before boxcar. It is a 2-d array with shape (nspec, nspat)

    Returns
    -------
    fextract:   ndarray
       Extracted flux at positions specified by (left<-->right, ycen). The output will have the same shape as
       Left and Right, i.e.  an 2-d  array with shape (nspec, nTrace) array if multiple traces were input, or a 1-d array with shape (nspec) for
       the case of a single trace.


    Revision History
    ----------------
    24-Mar-1999  Written by David Schlegel, Princeton.
    17-Feb-2003  Written with slow IDL routine, S. Burles, MIT
    22-Apr-2018  Ported to python by Joe Hennawi
    """

    # ToDO it would be nice to avoid this transposing, but I got confused during the IDL port
    left = left_in.T
    right = right_in.T

    dim = left.shape
    ndim = left.ndim
    if (ndim == 1):
        nTrace = 1
        npix = dim[0]
    else:
        nTrace = dim[0]
        npix = dim[1]

    if ycen is None:
        if ndim == 1:
            ycen_out = np.arange(npix, dtype=int)
        elif ndim == 2:
            ycen_out = np.outer(np.ones(nTrace, dtype=int), np.arange(npix, dtype=int))
        else:
            raise ValueError('trace is not 1 or 2 dimensional')
    else:
        ycen_out = ycen.T
        ycen_out = np.rint(ycen_out).astype(int)

    if ((np.size(left) != np.size(ycen_out)) | (np.shape(left) != np.shape(ycen_out))):
        raise ValueError('Number of elements and left of trace and ycen must be equal')

    idims = image.shape
    nspat = idims[1]
    nspec = idims[0]

    maxwindow = np.max(right - left)
    tempx = np.int(maxwindow + 3.0)

    bigleft = np.outer(left[:], np.ones(tempx))
    bigright = np.outer(right[:], np.ones(tempx))
    spot = np.outer(np.ones(npix * nTrace), np.arange(tempx)) + bigleft - 1
    bigy = np.outer(ycen_out[:], np.ones(tempx, dtype='int'))

    fullspot = np.array(np.fmin(np.fmax(np.round(spot + 1) - 1, 0), nspat - 1), int)
    fracleft = np.fmax(np.fmin(fullspot - bigleft, 0.5), -0.5)
    fracright = np.fmax(np.fmin(bigright - fullspot, 0.5), -0.5)
    del bigleft
    del bigright
    bool_mask1 = (spot >= -0.5) & (spot < (nspat - 0.5))
    bool_mask2 = (bigy >= 0) & (bigy <= (nspec - 1))
    weight = (np.fmin(np.fmax(fracleft + fracright, 0), 1)) * bool_mask1 * bool_mask2
    del spot
    del fracleft
    del fracright
    bigy = np.fmin(np.fmax(bigy, 0), nspec - 1)

    if weight_image is not None:
        temp = np.array([weight_image[x1, y1] * image[x1, y1] for (x1, y1) in zip(bigy.flatten(), fullspot.flatten())])
        temp2 = np.reshape(weight.flatten() * temp, (nTrace, npix, tempx))
        fextract = np.sum(temp2, axis=2)
        temp_wi = np.array([weight_image[x1, y1] for (x1, y1) in zip(bigy.flatten(), fullspot.flatten())])
        temp2_wi = np.reshape(weight.flatten() * temp_wi, (nTrace, npix, tempx))
        f_ivar = np.sum(temp2_wi, axis=2)
        fextract = fextract / (f_ivar + (f_ivar == 0)) * (f_ivar > 0)
    else:
        # Might be a more pythonic way to code this. I needed to switch the flattening order in order to get
        # this to work
        temp = np.array([image[x1, y1] for (x1, y1) in zip(bigy.flatten(), fullspot.flatten())])
        temp2 = np.reshape(weight.flatten() * temp, (nTrace, npix, tempx))
        fextract = np.sum(temp2, axis=2)

    # IDL version model functionality not implemented yet
    # At the moment I'm not reutnring the f_ivar for the weight_image mode. I'm not sure that this functionality is even
    # ever used

    if(nTrace ==1):
        fextract = fextract.reshape(npix)
    return fextract.T


def iter_tracefit(image, xinit_in, ncoeff, inmask = None, trc_inmask = None, fwhm = 3.0, maxdev = 2.0, maxiter = 25,
                  niter=9, gweight=False, show_fits=False, idx = None, verbose=False, xmin= None, xmax = None):
    """ Utility routine for object find to iteratively trace and fit. Used by both objfind and ech_objfind

    Parameters
    ----------
    image: ndaarray, float
        Image of objects to be traced
    xinit_in: ndarray, float
        Initial guesses for spatial direction trace. This can either be an 2-d  array with shape
        (nspec, nTrace) array, or a 1-d array with shape (nspec) for the case of a single trace.
    ncoeff: int
        Order of polynomial fits to trace

    Optional Parameter
    ------------------
    inmask: ndarray, bool
        Input mask for the image
    trc_inmask: ndarray, bool
        Input mask for the trace, i.e. places where you know the trace is going to be bad that you always want to mask in the
        fits. Same size as xinit_in (nspec, nTrace)
    fwhm: float
        fwhm width parameter for the flux or gaussian weighted tracing. For flux weighted trace the code does a third
        of the iterations with window 1.3*fwhm, a third with 1.1*fwhm, and a third with fwhm. For Gaussian weighted tracing
        it uses the fwhm to determine the sigma of the Gausisan which is used for all iterations.
    gweight: bool, default = False
        If gweight is True the routine does Gaussian weighted tracing, if it is False it will do flux weighted tracing.
        Normally the approach is to do a round of flux weighted tracing first, and then refine the traces with Gaussian
        weighted tracing.
    show_fits: bool, default = False
        Will plot the data and the fits.
    idx: ndarray of strings, default = None
        Array of idx IDs for each object. Used only if show_fits is true for the plotting.
    xmin: float, default = None
        Lower reference for robust_polyfit polynomial fitting. Default is to use zero
    xmax: float, defualt = None
        Upper refrence for robust_polyfit polynomial fitting.  Default is to use the image size in nspec direction
    Returns
    -------
    xpos: ndarray, float
       The output has the same size as xinit_in and contains the fit to the spatial direction of trace for each
       object.



    Revision History
    ----------------
    23-June-2018  Written by J. Hennawi
    """


    if inmask is None:
        inmask = np.ones_like(image,dtype=bool)

    # Allow for single vectors as input as well:
    nspec = xinit_in.shape[0]

    if xmin is None:
        xmin = 0.0
    if xmax is None:
        xmax = float(nspec-1)

    # Deal with the possibility of vectors as inputs instead of 2d arrays
    if xinit_in.ndim == 1:
        nobj = 1
        xinit = xinit_in.reshape(nspec,1)
        if trc_inmask is not None:
            trc_inmask_out = trc_inmask.reshape(nspec,1)
        else:
            trc_inmask_out = np.ones_like(xinit,dtype=bool)
    else:
        nobj = xinit_in.shape[1]
        xinit = xinit_in
        if trc_inmask is not None:
            trc_inmask_out = trc_inmask
        else:
            trc_inmask_out = np.ones_like(xinit,dtype=bool)

    spec_vec = np.arange(nspec)

    if verbose:
        msgs.info('Fitting the object traces')

    # Iterate flux weighted centroiding
    fwhm_vec = np.zeros(niter)
    fwhm_vec[0:niter//3] = 1.3*fwhm
    fwhm_vec[niter//3:2*niter//3] = 1.1*fwhm
    fwhm_vec[2*niter//3:] = fwhm

    if gweight:
        title_text = 'Gaussian Weighted'
    else:
        title_text = 'Flux Weighted'

    xfit1 = np.copy(xinit)

    for iiter in range(niter):
        if gweight:
            xpos1, xerr1 = trace_slits.trace_gweight(image*inmask,xfit1, invvar=inmask.astype(float),sigma=fwhm/2.3548)
        else:
            xpos1, xerr1 = trace_slits.trace_fweight(image*inmask,xfit1, invvar = inmask.astype(float), radius = fwhm_vec[iiter])

        # If a trc_inmask was input, always set the masked values to the initial input crutch. The point is that the crutch
        # initially comes from either the standard or the slit boundaries, and if we continually replace it for all iterations
        # we will naturally always extraplate the trace to match the shape of a high S/N ratio fit (i.e. either the standard)
        # or the flat which was used to determine the slit edges.
        xpos1[np.invert(trc_inmask_out)] = xinit[np.invert(trc_inmask_out)]

        # Do not do any kind of masking based on xerr1. Trace fitting is much more robust when masked pixels are simply
        # replaced by the tracing crutch. We thus do not do weighted fits, i.e. uniform weights, but we set the relative
        # weight of the trc_inmask pixels to be lower. This way they still get a say but do not overly influence the fit.
        xinvvar = np.ones_like(xpos1.T)
        xinvvar[np.invert(trc_inmask_out.T)] = 0.1
        pos_set1 = pydl.xy2traceset(np.outer(np.ones(nobj),spec_vec), xpos1.T,
                                    #inmask = trc_inmask_out.T,
                                    ncoeff=ncoeff, maxdev=maxdev,
                                    maxiter=maxiter, invvar=xinvvar, xmin=xmin, xmax =xmax)
        xfit1 = pos_set1.yfit.T
        # bad pixels have errors set to 999 and are returned to lie on the input trace. Use this only for plotting below
        #errmask = (xerr1 > 990.0)  # bad pixels have errors set to 999 and are returned to lie on the input trace
        outmask = pos_set1.outmask.T
        # Plot all the points that were not masked initially
        if(show_fits) & (iiter == niter - 1):
            for iobj in range(nobj):
                # The sum of all these masks adds up to the number of pixels.
                inmask_trc = np.invert(trc_inmask_out[:,iobj]) # masked on the way in
                errmask = xerr1[:,iobj] > 990.0 # masked by fweight or gweight, was set to input trace and still fit
                rejmask = np.invert(outmask[:, iobj]) & np.invert(inmask_trc) # was good on the way in, masked by the poly fit
                nomask = outmask[:, iobj] & np.invert(errmask) # was actually fit and not replaced to input trace
                plt.plot(spec_vec[nomask],xpos1[nomask,iobj],marker='o', c='k', markersize=3.0,linestyle='None',label=title_text + ' Centroid')
                plt.plot(spec_vec,xinit[:,iobj],c='g', zorder = 25, linewidth=2.0,linestyle='--', label='initial guess')
                plt.plot(spec_vec,xfit1[:,iobj],c='red',zorder=30,linewidth = 2.0, label ='fit to trace')
                if np.any(errmask):
                    plt.plot(spec_vec[errmask],xfit1[errmask,iobj], c='blue',marker='+',
                             markersize=5.0,linestyle='None',zorder= 20, label='masked by tracing, set to init guess')
                if np.any(rejmask):
                    plt.plot(spec_vec[rejmask],xpos1[rejmask,iobj], c='cyan',marker='v',
                             markersize=5.0,linestyle='None',zorder= 20, label='masked by polynomial fit')
                if np.any(inmask_trc):
                    plt.plot(spec_vec[inmask_trc],xpos1[inmask_trc,iobj],
                             c='orange',marker='s',markersize=3.0,linestyle='None',zorder= 20, label='input masked points, not fit')
                try:
                    plt.title(title_text + ' Centroid to object {:s}.'.format(idx[iobj]))
                except TypeError:
                    plt.title(title_text + ' Centroid to object {:d}.'.format(iobj))
                plt.ylim((0.995*xfit1[:, iobj].min(), 1.005*xfit1[:, iobj].max()))
                plt.xlabel('Spectral Pixel')
                plt.ylabel('Spatial Pixel')
                plt.legend()
                plt.show()

    # Returns the fit, the actual weighted traces, and the pos_set1 object
    return xfit1, xpos1, xerr1, pos_set1



# TODO: JFH It would be really ideal if we could replace this pca with a weighted PCA!!
def pca_trace(xinit_in, spec_min_max=None, predict = None, npca = None, pca_explained_var=99.0,
              coeff_npoly = None, coeff_weights=None, debug=True, order_vec = None, lower = 3.0,
              upper = 3.0, minv = None,maxv = None, maxrej=1,
              xinit_mean = None):

    """
    Use a PCA model to determine the best object (or slit edge) traces for echelle spectrographs.

    Args:
      xinit:  ndarray, (nspec, norders)
         Array of input traces that one wants to PCA model. For object finding this will be the traces for orders where
         an object was detected. If an object was not detected on some orders (see ech_objfind), the standard star
         (or order boundaries)  will be  assigned to these orders at the correct fractional slit position, and a joint PCA
         fit will be performed to the detected traces and the standard/slit traces.

    spec_min_max: float or int ndarray, (2, norders), default=None.
         This is a 2-d array which defines the minimum and maximum of each order in the
         spectral direction on the detector. This should only be used for echelle spectrographs for which the orders do not
         entirely cover the detector, and each order passed in for xinit_in is a succession of orders on the detector.
         The code will re-map the traces such that they all have the same length, compute the PCA, and then re-map the orders
         back. This improves performanc for echelle spectrographs by removing the nonlinear shrinking of the orders so that
         the linear pca operation can better predict the traces. THIS IS AN EXPERIMENTAL FEATURE. INITIAL TESTS WITH
         XSHOOTER-NIR INDICATED THAT IT DID NOT IMPROVE PERFORMANCE AND SIMPLY LINEAR EXTRAPOLATION OF THE ORDERS INTO THE
         REGIONS THAT ARE NOT ILLUMINATED PERFORMED SIGNIFICANTLY BETTER. DO NOT USE UNTIL FURTHER TESTING IS PERFORMED. IT
         COULD HELP WITH OTHER MORE NONLINEAR SPECTROGRAPHS.
    predict: ndarray, bool (norders,), default = None
         Orders which have True are those that will be predicted by extrapolating the fit of the PCA coefficents for those
         orders which have False set in this array. The default is None, which means that the coefficients of all orders
         will be fit simultaneously and no extrapolation will be performed. For object finding, we use the standard star
         (or slit boundaries) as the input for orders for which a trace is not identified and fit the coefficients of all
         simultaneously. Thus no extrapolation is performed. For tracing slit boundaries it may be useful to perform
          extrapolations.
    npca: int, default = None
         number of PCA components to be kept. The maximum number of possible PCA components would be = norders, which is to say
         that no PCA compression woulud be performed. For the default of None, npca will be automatically determinedy by
         calculating the minimum number of components required to explain 99% (pca_explained_var) of the variance in the different orders.
    pca_explained_var: float, default = 99
         Amount of explained variance cut used to determine where to truncate the PCA, i.e. to determine npca.
    coeff_npoly: int, default = None
         Order of polynomial fits used for PCA coefficients fitting. The defualt is None, which means that coeff_noly
         will be automatically determined by taking the number of orders into account. PCA components that explain
         less variance (and are thus much noiser) are fit with lower order.
    coeff_weights (np.ndarray): shape = (norders,), default=None
         If input these weights will be used for the polynomial fit to the PCA coefficients. Even if you are predicting
         orders and hence only fitting a subset of the orders != norders, the shape of coeff_weights must be norders.
         Just give the orders you don't plan to fit a weight of zero. This option is useful for fitting object
         traces since the weights can be set to (S/N)^2 of each order.
         TODO: Perhaps we should get rid of the predict option and simply allow the user to set the weights of the orders
         they want predicted to be zero. That would be more straightforward, but would require a rework of the code.

    debug: bool, default = False
         Show plots useful for debugging.

    Returns:
    --------
    pca_fit:  ndarray, float (nspec, norders)
        Array with the same size as xinit, which contains the pca fitted orders.
    """

    nspec, norders = xinit_in.shape

    if order_vec is None:
        order_vec = np.arange(norders,dtype=float)

    if predict is None:
        predict = np.zeros(norders,dtype=bool)

    # use_order = True orders used to predict the predict = True bad orders
    use_order = np.invert(predict)
    ngood = np.sum(use_order)

    if ngood < 2:
        msgs.warn('There are no good traces to PCA fit. There is probably a bug somewhere. Exiting and returning input traces.')
        return xinit_in, {}, None, None

    if spec_min_max is not None:
        xinit = remap_orders(xinit_in, spec_min_max)
    else:
        xinit = xinit_in

    # Take out the mean position of each input trace
    if xinit_mean is None:
        xinit_mean = np.mean(xinit, axis=0)

    xpca = xinit - xinit_mean
    xpca_use = xpca[:, use_order].T
    pca_full = PCA()
    pca_full.fit(xpca_use)
    var = np.cumsum(np.round(pca_full.explained_variance_ratio_, decimals=6) * 100)
    npca_full = var.size
    if npca is None:
        if var[0]>=pca_explained_var:
            npca = 1
            msgs.info('The first PCA component contains more than {:5.3f} of the information'.format(pca_explained_var))
        else:
            npca = int(np.ceil(np.interp(pca_explained_var, var,np.arange(npca_full)+1)))
            msgs.info('Truncated PCA to contain {:5.3f}'.format(pca_explained_var) + '% of the total variance. ' +
                      'Number of components to keep is npca = {:d}'.format(npca))
    else:
        npca = int(npca)
        var_trunc = np.interp(float(npca),np.arange(npca_full)+1.0, var)
        msgs.info('Truncated PCA with npca={:d} components contains {:5.3f}'.format(npca, var_trunc) + '% of the total variance.')

    if npca_full < npca:
        msgs.warn('Not enough good traces for a PCA fit of the requested dimensionality. The full (non-compressing) PCA has size: '
                  'npca_full = {:d}'.format(npca_full) + ' is < npca = {:d}'.format(npca))
        msgs.warn('Using the input trace for now. But you should lower npca <= npca_full')
        return xinit_in, {}, None, None

    if coeff_npoly is None:
        coeff_npoly = int(np.fmin(np.fmax(np.floor(3.3*ngood/norders),1.0),3.0))


    # Polynomial coefficient for PCA coefficients
    npoly_vec =np.zeros(npca, dtype=int)
    # Fit first pca dimension (with largest variance) with a higher order npoly depending on number of good orders.
    # Fit all higher dimensions (with lower variance) with a line
    # Cascade down and use lower order polynomial for PCA directions that contain less variance
    for ipoly in range(npca):
        npoly_vec[ipoly] = np.fmax(coeff_npoly - ipoly,1)

        pca = PCA(n_components=npca)
        pca_coeffs_use = pca.fit_transform(xpca_use)
        pca_vectors = pca.components_

    pca_coeffs_new = np.zeros((norders, npca))
    fit_dict = {}
    # Now loop over the dimensionality of the compression and perform a polynomial fit to
    for idim in range(npca):
        # Only fit the use_order orders, then use this to predict the others
        xfit = order_vec[use_order]
        yfit = pca_coeffs_use[:,idim]
        ncoeff = npoly_vec[idim]
        # Apply a 10% relative error to each coefficient. This performs better than use_mad, since larger coefficients
        # will always be considered inliers, if the coefficients vary rapidly with order as they sometimes do.
        sigma = np.fmax(0.1*np.abs(yfit), 0.1)
        invvar = utils.inverse(sigma**2)
        use_weights =  coeff_weights[use_order] if coeff_weights is not None else None
        # TODO Note that we are doing a weighted fit using the coeff_weights, but the rejection is still done
        # usnig the ad-hoc invvar created in the line above. I cannot think of a better way.
        msk_new, poly_out = utils.robust_polyfit_djs(xfit, yfit, ncoeff, invvar = invvar, weights=use_weights,
                                                     function='polynomial', maxiter=25,
                                                     lower=lower, upper=upper,
                                                     maxrej=maxrej,
                                                     sticky=False, use_mad=False, minx = minv, maxx = maxv)
        # ToDO robust_poly_fit needs to return minv and maxv as outputs for the fits to be usable downstream
        pca_coeffs_new[:,idim] = utils.func_val(poly_out, order_vec, 'polynomial')
        fit_dict[str(idim)] = {}
        fit_dict[str(idim)]['coeffs'] = poly_out
        fit_dict[str(idim)]['minv'] = minv
        fit_dict[str(idim)]['maxv'] = maxv
        if debug:
            # Evaluate the fit
            xvec = np.linspace(order_vec.min(),order_vec.max(),num=100)
            robust_mask_new = msk_new == 1
            plt.plot(xfit, yfit, 'ko', mfc='None', markersize=8.0, label='pca coeff')
            plt.plot(xfit[~robust_mask_new], yfit[~robust_mask_new], 'r+', markersize=20.0,label='robust_polyfit_djs rejected')
            plt.plot(xvec, utils.func_val(poly_out, xvec, 'polynomial'),ls='-.', color='steelblue',
                     label='Polynomial fit of order={:d}'.format(ncoeff))
            plt.xlabel('Order Number', fontsize=14)
            plt.ylabel('PCA Coefficient', fontsize=14)
            plt.title('PCA Fit for Dimension #{:d}/{:d}'.format(idim + 1,npca))
            plt.legend()
            plt.show()

    pca_model = np.outer(pca.mean_,np.ones(norders)) + (np.dot(pca_coeffs_new, pca_vectors)).T
#   pca_model_mean = np.mean(pca_model,0)
#   pca_fit = np.outer(np.ones(nspec), xinit_mean) + (pca_model - pca_model_mean)
#   JFH which is correct?
    pca_fit = np.outer(np.ones(nspec), xinit_mean) + (pca_model)

    if spec_min_max is not None:
        pca_out = remap_orders(pca_fit, spec_min_max, inverse=True)
    else:
        pca_out = pca_fit

    return pca_out, fit_dict, pca.mean_, pca_vectors



