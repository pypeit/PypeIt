"""
Module for PypeIt extraction code

.. include:: ../include/links.rst

"""
import copy

import numpy as np
import scipy
from matplotlib import pyplot as plt

from IPython import embed

from astropy import stats

from pypeit import msgs
from pypeit import utils
from pypeit import specobj
from pypeit import specobjs
from pypeit import tracepca
from pypeit import bspline
from pypeit.display import display
from pypeit.core import pydl
from pypeit.core import pixels
from pypeit.core import arc
from pypeit.core import fitting
from pypeit.core import procimg
from pypeit.core.trace import fit_trace
from pypeit.core.moment import moment1d


def extract_optimal(sciimg, ivar, mask, waveimg, skyimg, thismask, oprof,
                    spec, min_frac_use=0.05, base_var=None, count_scale=None, noise_floor=None):

    """
    Calculate the spatial FWHM from an object profile. Utility routine for
    fit_profile

    The specobj object is changed in place with the boxcar and optimal
    dictionaries being filled with the extraction parameters.

    Parameters
    ----------
    sciimg : float `numpy.ndarray`_, shape (nspec, nspat)
        Science frame
    ivar : float `numpy.ndarray`_, shape (nspec, nspat)
        Inverse variance of science frame. Can be a model or deduced from the
        image itself.
    mask : boolean `numpy.ndarray`_, shape (nspec, nspat)
        Good-pixel mask, indicating which pixels should or should not be
        used. Good pixels = True, Bad Pixels = False
    waveimg : float `numpy.ndarray`_, shape (nspec, nspat)
        Wavelength image.
    skyimg : float `numpy.ndarray`_, shape (nspec, nspat)
        Image containing our model of the sky
    thismask : boolean `numpy.ndarray`_, shape (nspec, nspat)
        Image indicating which pixels are on the slit/order in question.
        True=Good.
    oprof : float `numpy.ndarray`_, shape (nspec, nspat)
        Image containing the profile of the object that we are extracting.
    spec : :class:`~pypeit.specobj.SpecObj`
        This is the container that holds object, trace, and extraction
        information for the object in question. This routine operates one object
        at a time.  **This object is altered in place!**
    min_frac_use : :obj:`float`, optional
        If the sum of object profile across the spatial direction are less than
        this value, the optimal extraction of this spectral pixel is masked
        because the majority of the object profile has been masked.
    base_var : `numpy.ndarray`_, shape is (nspec, nspat), optional
        The "base-level" variance in the data set by the detector properties and
        the image processing steps.  See
        :func:`~pypeit.core.procimg.base_variance`.
    count_scale : :obj:`float`, `numpy.ndarray`_, optional
        A scale factor, :math:`s`, that *has already been applied* to the
        provided science image.  For example, if the image has been flat-field
        corrected, this is the inverse of the flat-field counts.  If None, set
        to 1.  If a single float, assumed to be constant across the full image.
        If an array, the shape must match ``base_var``.  The variance will be 0
        wherever :math:`s \leq 0`, modulo the provided ``adderr``.  This is one
        of the components needed to construct the model variance; see
        ``model_noise``.
    noise_floor : :obj:`float`, optional
        A fraction of the counts to add to the variance, which has the effect of
        ensuring that the S/N is never greater than ``1/noise_floor``; see
        :func:`~pypeit.core.procimg.variance_model`.  If None, no noise floor is
        added.
    """
    # Setup
    imgminsky = sciimg - skyimg
    nspat = imgminsky.shape[1]
    nspec = imgminsky.shape[0]

    spec_vec = np.arange(nspec)
    spat_vec = np.arange(nspat)

    # TODO This makes no sense for difference imaging? Not sure we need NIVAR anyway
    var_no = None if base_var is None \
                else procimg.variance_model(base_var, counts=skyimg, count_scale=count_scale,
                                            noise_floor=noise_floor)

    ispec, ispat = np.where(oprof > 0.0)

    # Exit gracefully if we have no positive object profiles, since that means something was wrong with object fitting
    if not np.any(oprof > 0.0):
        msgs.warn('Object profile is zero everywhere. This aperture is junk.')
        return

    mincol = np.min(ispat)
    maxcol = np.max(ispat) + 1
    nsub = maxcol - mincol

    mask_sub = mask[:,mincol:maxcol]
    thismask_sub = thismask[:, mincol:maxcol]
    wave_sub = waveimg[:,mincol:maxcol]
    ivar_sub = np.fmax(ivar[:,mincol:maxcol],0.0) # enforce positivity since these are used as weights
    vno_sub = None if var_no is None else np.fmax(var_no[:,mincol:maxcol],0.0)

    base_sub = None if base_var is None else base_var[:,mincol:maxcol]
    img_sub = imgminsky[:,mincol:maxcol]
    sky_sub = skyimg[:,mincol:maxcol]
    oprof_sub = oprof[:,mincol:maxcol]
    # enforce normalization and positivity of object profiles
    norm = np.nansum(oprof_sub,axis = 1)
    norm_oprof = np.outer(norm, np.ones(nsub))
    oprof_sub = np.fmax(oprof_sub/norm_oprof, 0.0)

    ivar_denom = np.nansum(mask_sub*oprof_sub, axis=1)
    mivar_num = np.nansum(mask_sub*ivar_sub*oprof_sub**2, axis=1)
    mivar_opt = mivar_num/(ivar_denom + (ivar_denom == 0.0))
    flux_opt = np.nansum(mask_sub*ivar_sub*img_sub*oprof_sub, axis=1)/(mivar_num + (mivar_num == 0.0))
    # Optimally extracted noise variance (sky + read noise) only. Since
    # this variance is not the same as that used for the weights, we
    # don't get the usual cancellation. Additional denom factor is the
    # analog of the numerator in Horne's variance formula. Note that we
    # are only weighting by the profile (ivar_sub=1) because
    # otherwise the result depends on the signal (bad).
    nivar_num = np.nansum(mask_sub*oprof_sub**2, axis=1) # Uses unit weights
    if vno_sub is None:
        nivar_opt = None
    else:
        nvar_opt = ivar_denom * np.nansum(mask_sub * vno_sub * oprof_sub**2, axis=1) \
                            / (nivar_num**2 + (nivar_num**2 == 0.0))
        nivar_opt = 1.0/(nvar_opt + (nvar_opt == 0.0))
    # Optimally extract sky and (read noise)**2 in a similar way
    sky_opt = ivar_denom*(np.nansum(mask_sub*sky_sub*oprof_sub**2, axis=1))/(nivar_num**2 + (nivar_num**2 == 0.0))
    if base_var is None:
        base_opt = None
    else:
        base_opt = ivar_denom * np.nansum(mask_sub * base_sub * oprof_sub**2, axis=1) \
                        / (nivar_num**2 + (nivar_num**2 == 0.0))
        base_opt = np.sqrt(base_opt)
        base_opt[np.isnan(base_opt)]=0.0

    tot_weight = np.nansum(mask_sub*ivar_sub*oprof_sub, axis=1)
    prof_norm = np.nansum(oprof_sub, axis=1)
    frac_use = (prof_norm > 0.0)*np.nansum((mask_sub*ivar_sub > 0.0)*oprof_sub, axis=1)/(prof_norm + (prof_norm == 0.0))

    # Use the same weights = oprof^2*mivar for the wavelenghts as the flux.
    # Note that for the flux, one of the oprof factors cancels which does
    # not for the wavelengths.
    wave_opt = np.nansum(mask_sub*ivar_sub*wave_sub*oprof_sub**2, axis=1)/(mivar_num + (mivar_num == 0.0))
    mask_opt = (tot_weight > 0.0) & (frac_use > min_frac_use) & (mivar_num > 0.0) & (ivar_denom > 0.0) & \
               np.isfinite(wave_opt) & (wave_opt > 0.0)

    # Interpolate wavelengths over masked pixels
    badwvs = (mivar_num <= 0) | np.invert(np.isfinite(wave_opt)) | (wave_opt <= 0.0)
    if badwvs.any():
        oprof_smash = np.nansum(thismask_sub*oprof_sub**2, axis=1)
        # Can we use the profile average wavelengths instead?
        oprof_good = badwvs & (oprof_smash > 0.0)
        if oprof_good.any():
            wave_opt[oprof_good] = np.nansum(
                wave_sub[oprof_good,:]*thismask_sub[oprof_good,:]*oprof_sub[oprof_good,:]**2, axis=1)/\
                                   np.nansum(thismask_sub[oprof_good,:]*oprof_sub[oprof_good,:]**2, axis=1)
        oprof_bad = badwvs & ((oprof_smash <= 0.0) | (np.isfinite(oprof_smash) == False) | (wave_opt <= 0.0) | (np.isfinite(wave_opt) == False))
        if oprof_bad.any():
            # For pixels with completely bad profile values, interpolate from trace.
            f_wave = scipy.interpolate.RectBivariateSpline(spec_vec,spat_vec, waveimg*thismask)
            wave_opt[oprof_bad] = f_wave(spec.trace_spec[oprof_bad], spec.TRACE_SPAT[oprof_bad],
                                         grid=False)

    flux_model = np.outer(flux_opt,np.ones(nsub))*oprof_sub
    chi2_num = np.nansum((img_sub - flux_model)**2*ivar_sub*mask_sub,axis=1)
    chi2_denom = np.fmax(np.nansum(ivar_sub*mask_sub > 0.0, axis=1) - 1.0, 1.0)
    chi2 = chi2_num/chi2_denom

    # Fill in the optimally extraction tags
    spec.OPT_WAVE = wave_opt    # Optimally extracted wavelengths
    spec.OPT_COUNTS = flux_opt    # Optimally extracted flux
    spec.OPT_COUNTS_IVAR = mivar_opt   # Inverse variance of optimally extracted flux using modelivar image
    spec.OPT_COUNTS_SIG = np.sqrt(utils.inverse(mivar_opt))
    spec.OPT_COUNTS_NIVAR = nivar_opt  # Optimally extracted noise variance (sky + read noise) only
    spec.OPT_MASK = mask_opt    # Mask for optimally extracted flux
    spec.OPT_COUNTS_SKY = sky_opt      # Optimally extracted sky
    spec.OPT_COUNTS_SIG_DET = base_opt      # Square root of optimally extracted read noise squared
    spec.OPT_FRAC_USE = frac_use    # Fraction of pixels in the object profile subimage used for this extraction
    spec.OPT_CHI2 = chi2            # Reduced chi2 of the model fit for this spectral pixel


def extract_asym_boxcar(sciimg, left_trace, righ_trace, gpm=None, ivar=None):
    """
    Perform assymetric boxcar extraction of the flux between two traces.

    Parameters
    ----------
    sciimg : float `numpy.ndarray`_, shape (nspec, nspat)
        Science frame for the extraction

    left_trace : float `numpy.ndarray`_, shape (nspec, napertures)
        Left trace boundary of the extraction region.

    right_trace : float `numpy.ndarray`_, shape (nspec, napertures)
        Right trace boundary of the extraction region.

    gpm : boolean `numpy.ndarray`_, shape (nspec, nspat)
        Good-pixel mask, indicating which pixels are should or should not be
        used. Good pixels = True, Bad Pixels = False. Optional.
    ivar : float `numpy.ndarray`_, shape (nspec, nspat)
        Inverse variance of science frame. Can be a model or deduced from the
        image itself. If ivar is set then this code will return the
        error on the extraction. Optional

    Returns
    -------
    flux_out : float `numpy.ndarray`_, shape (nspec, napertures)
        Array containing the boxcar extracted flux as a function of spectral position for each aperture.
    gpm_box : boool `numpy.ndarray`_, shape (nspec, napertures)
        Good pixel mask for the boxcar extracted flux
    box_npix :   float `numpy.ndarray`_, shape (nspec, napertures)
        Array containing the number of pixels which contributed to the boxcar sum of the flux for each
        spectral position for each aperture.
    ivar_out : float `numpy.ndarray`_, shape (nspec, napertures)
        Array containing the inverse variance of the boxcar extracted flux as a function of spectral position
        for each aperture. This will only be returned if the input parameter ivar is not None.

    """
    ivar1 = np.ones_like(sciimg) if ivar is None else ivar
    gpm1 = ivar1 > 0.0 if gpm is None else gpm

    flux_box = moment1d(sciimg*gpm1, (left_trace+righ_trace)/2.0, (righ_trace-left_trace))[0]
    #box_denom = moment1d(gpm1, (left_trace+righ_trace)/2.0, (righ_trace-left_trace))[0]

    pixtot = moment1d(sciimg*0 + 1.0, (left_trace+righ_trace)/2.0, (righ_trace-left_trace))[0]
    pixmsk = moment1d(ivar1*gpm1 == 0.0, (left_trace+righ_trace)/2.0, (righ_trace-left_trace))[0]

    # If every pixel is masked then mask the boxcar extraction
    gpm_box = (pixmsk != pixtot)

    varimg = 1.0 / (ivar1 + (ivar1 == 0.0))
    var_box = moment1d(varimg * gpm1, (left_trace+righ_trace)/2.0, (righ_trace-left_trace))[0]

    ivar_box = 1.0/(var_box + (var_box == 0.0))

    flux_out = flux_box*gpm_box
    ivar_out = ivar_box*gpm_box
    box_npix = pixtot - pixmsk

    if ivar is None:
        return flux_out, gpm_box, box_npix
    else:
        return flux_out, gpm_box, box_npix, ivar_out


def extract_boxcar(sciimg, ivar, mask, waveimg, skyimg, spec, base_var=None,
                   count_scale=None, noise_floor=None):
    """
    Perform boxcar extraction for a single SpecObj.  The size of the boxcar must
    be available as an attribute of :class:`~pypeit.specobj.SpecObj`.


    Note that the provided :class:`~pypeit.specobj.SpecObj` (``spec``) is
    modified in place.

    Parameters
    ----------
    sciimg : float `numpy.ndarray`_, shape (nspec, nspat)
        Science frame
    ivar : float `numpy.ndarray`_, shape (nspec, nspat)
        Inverse variance of science frame. Can be a model or deduced from the
        image itself.
    mask : boolean `numpy.ndarray`_, shape (nspec, nspat)
        Good-pixel mask, indicating which pixels are should or should not be
        used. Good pixels = True, Bad Pixels = False
    waveimg : float `numpy.ndarray`_, shape (nspec, nspat)
        Wavelength image.
    skyimg : float `numpy.ndarray`_, shape (nspec, nspat)
        Image containing our model of the sky
    spec : :class:`~pypeit.specobj.SpecObj`
        This is the container that holds object, trace, and extraction
        information for the object in question. This routine operates one object
        at a time.  **This object is altered in place!**
    base_var : `numpy.ndarray`_, shape is (nspec, nspat), optional
        The "base-level" variance in the data set by the detector properties and
        the image processing steps.  See
        :func:`~pypeit.core.procimg.base_variance`.
    count_scale : :obj:`float`, `numpy.ndarray`_, optional
        A scale factor, :math:`s`, that *has already been applied* to the
        provided science image.  For example, if the image has been flat-field
        corrected, this is the inverse of the flat-field counts.  If None, set
        to 1.  If a single float, assumed to be constant across the full image.
        If an array, the shape must match ``base_var``.  The variance will be 0
        wherever :math:`s \leq 0`, modulo the provided ``adderr``.  This is one
        of the components needed to construct the model variance; see
        ``model_noise``.
    noise_floor : :obj:`float`, optional
        A fraction of the counts to add to the variance, which has the effect of
        ensuring that the S/N is never greater than ``1/noise_floor``; see
        :func:`~pypeit.core.procimg.variance_model`.  If None, no noise floor is
        added.
    """
    # Setup
    imgminsky = sciimg - skyimg
    nspat = imgminsky.shape[1]
    nspec = imgminsky.shape[0]

    spec_vec = np.arange(nspec)
    spat_vec = np.arange(nspat)
    if spec.trace_spec is None:
        spec.trace_spec = spec_vec

    # get boxcar_radius
    box_radius = spec.BOX_RADIUS

    # TODO This makes no sense for difference imaging? Not sure we need NIVAR anyway
    var_no = None if base_var is None \
                else procimg.variance_model(base_var, counts=skyimg, count_scale=count_scale,
                                            noise_floor=noise_floor)

    # Fill in the boxcar extraction tags
    flux_box = moment1d(imgminsky*mask, spec.TRACE_SPAT, 2*box_radius, row=spec.trace_spec)[0]
    # Denom is computed in case the trace goes off the edge of the image
    box_denom = moment1d(waveimg*mask > 0.0, spec.TRACE_SPAT, 2*box_radius,
                         row=spec.trace_spec)[0]
    wave_box = moment1d(waveimg*mask, spec.TRACE_SPAT, 2*box_radius,
                        row=spec.trace_spec)[0] / (box_denom + (box_denom == 0.0))
    varimg = 1.0/(ivar + (ivar == 0.0))
    var_box = moment1d(varimg*mask, spec.TRACE_SPAT, 2*box_radius, row=spec.trace_spec)[0]
    nvar_box = None if var_no is None \
                else moment1d(var_no*mask, spec.TRACE_SPAT, 2*box_radius, row=spec.trace_spec)[0]
    sky_box = moment1d(skyimg*mask, spec.TRACE_SPAT, 2*box_radius, row=spec.trace_spec)[0]
    if base_var is None:
        base_box = None
    else:
        _base_box = moment1d(base_var*mask, spec.TRACE_SPAT, 2*box_radius, row=spec.trace_spec)[0]
        base_posind = (_base_box > 0.0)
        base_box = np.zeros(_base_box.shape, dtype=float)
        base_box[base_posind] = np.sqrt(_base_box[base_posind])
    pixtot = moment1d(ivar*0 + 1.0, spec.TRACE_SPAT, 2*box_radius, row=spec.trace_spec)[0]
    pixmsk = moment1d(ivar*mask == 0.0, spec.TRACE_SPAT, 2*box_radius, row=spec.trace_spec)[0]
    # If every pixel is masked then mask the boxcar extraction
    mask_box = (pixmsk != pixtot) & np.isfinite(wave_box) & (wave_box > 0.0)
    bad_box = (wave_box <= 0.0) | np.invert(np.isfinite(wave_box)) | (box_denom == 0.0)
    # interpolate bad wavelengths over masked pixels
    if bad_box.any():
        f_wave = scipy.interpolate.RectBivariateSpline(spec_vec, spat_vec, waveimg)
        wave_box[bad_box] = f_wave(spec.trace_spec[bad_box], spec.TRACE_SPAT[bad_box], grid=False)

    ivar_box = 1.0/(var_box + (var_box == 0.0))
    nivar_box = None if nvar_box is None else 1.0/(nvar_box + (nvar_box == 0.0))

    # Fill em up!
    spec.BOX_WAVE = wave_box
    spec.BOX_COUNTS = flux_box*mask_box
    spec.BOX_COUNTS_IVAR = ivar_box*mask_box
    spec.BOX_COUNTS_SIG = np.sqrt(utils.inverse(ivar_box*mask_box))
    spec.BOX_COUNTS_NIVAR = None if nivar_box is None else nivar_box*mask_box
    spec.BOX_MASK = mask_box
    spec.BOX_COUNTS_SKY = sky_box
    spec.BOX_COUNTS_SIG_DET = base_box
    # TODO - Confirm this should be float, not int
    spec.BOX_NPIX = pixtot-pixmsk


def findfwhm(model, sig_x):
    """ Calculate the spatial FWHM from an object profile. Utitlit routine for fit_profile

    Parameters
    ----------
    model :   numpy float 2-d array [nspec, nspat]
    x :

    Returns
    -------
    peak :  Peak value of the profile model
    peak_x:  sig_x location where the peak value is obtained
    lwhm:   Value of sig_x at the left width at half maximum
    rwhm:   Value of sig_x at the right width at half maximum

    Notes
    -----
    Revision History
        - 11-Mar-2005  Written by J. Hennawi and S. Burles David Schlegel, Princeton.
        - 28-May-2018  Ported to python by J. Hennawi
    """


    peak = (model*(np.abs(sig_x) < 1.)).max()
    peak_x = sig_x[(model*(np.abs(sig_x) < 1.)).argmax()]

    lrev = ((sig_x < peak_x) & (model < 0.5*peak))[::-1]
    lind, = np.where(lrev)
    if(lind.size > 0):
        lh = lind.min()
        lwhm = (sig_x[::-1])[lh]
    else:
        lwhm = -0.5*2.3548

    rind, = np.where((sig_x > peak_x) & (model < 0.5*peak))
    if(rind.size > 0):
        rh = rind.min()
        rwhm = sig_x[rh]
    else:
        rwhm = 0.5 * 2.3548

    return (peak, peak_x, lwhm, rwhm)



def qa_fit_profile(x_tot,y_tot, model_tot, l_limit = None, r_limit = None, ind = None,
                   title =' ', xtrunc = 1e6, xlim = None, ylim = None, qafile = None):

    # Plotting pre-amble
    plt.close("all")
    #plt.clf()
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')
    width = 10.0 # Golden ratio 1.618
    fig, ax = plt.subplots(1, figsize=(width, width/1.618))

    if ind is None:
        indx = np.slice(x_tot.size)
    else:
        if len(ind) == 0:
            indx = np.slice(x_tot.size)
            title = title + ': no good pixels, showing all'
        else:
            indx = ind

    x = x_tot.flat[indx]
    y = y_tot.flat[indx]
    model = model_tot.flat[indx]

    ax.plot(x,y,color='k',marker='o',markersize= 0.3, mfc='k',fillstyle='full',linestyle='None')


    max_model = np.fmin(model.max(), 2.0)
    if xlim is None:
        goodpix = model > 0.001*max_model
        if goodpix.any():
            minx = np.fmax(1.1*(x[goodpix]).min(), -xtrunc)
            maxx = np.fmin(1.1*(x[goodpix]).max(),  xtrunc)
        else:
            minx = -5.0
            maxx = 5.0
    else:
        minx = -xlim
        maxx = xlim
    xlimit = (minx, maxx)

    nsamp = 150
    half_bin = (maxx - minx)/nsamp/2.0
    if ylim is None:
        ymax = np.fmax(1.5*model.max(), 0.1)
        ymin = np.fmin(-0.1*ymax, -0.05)
        ylim = (ymin, ymax)
    plot_mid = (np.arange(nsamp) + 0.5)/nsamp*(maxx - minx) + minx

    y20 = np.zeros(nsamp)
    y80 = np.zeros(nsamp)
    y50 = np.zeros(nsamp)
    model_samp = np.zeros(nsamp)
    nbin = np.zeros(nsamp)

    for i in range(nsamp):
        dist = np.abs(x - plot_mid[i])
        close = dist < half_bin
        yclose = y[close]
        nclose = close.sum()
        nbin[i] = nclose
        if close.any():
            closest = (dist[close]).argmin()
            model_samp[i] = (model[close])[closest]
        if nclose > 3:
            s = yclose.argsort()
            y50[i] = yclose[s[int(np.rint((nclose - 1)*0.5))]]
            y80[i] = yclose[s[int(np.rint((nclose - 1)*0.8))]]
            y20[i] = yclose[s[int(np.rint((nclose - 1)*0.2))]]
            ax.plot([plot_mid[i],plot_mid[i]], [y20[i],y80[i]], linewidth=1.2, color='orange')

    icl = nbin > 3
    if icl.any():
        ax.plot(plot_mid[icl],y50[icl],marker = 'o', color='lime', markersize=2, fillstyle='full', linestyle='None')
    else:
        ax.plot(plot_mid, y50, marker='o', color='lime', markersize=2, fillstyle = 'full', linestyle='None')

    isort = x.argsort()
    ax.plot(x[isort], model[isort], color='red', linewidth=1.0)



    if l_limit is not None:
        ax.axvline(x =l_limit, color='cornflowerblue',linewidth=2.0)
    if r_limit is not None:
        ax.axvline(x=r_limit, color='cornflowerblue',linewidth=2.0)

    ax.set_xlim(xlimit)
    ax.set_ylim(ylim)
    ax.set_title(title)
    ax.set_xlabel(r'$x/\sigma$')
    ax.set_ylabel('Normalized Profile')

    fig.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
    plt.show()
    return

# TODO JFH Split up the profile part and the QA part for cleaner code
def return_gaussian(sigma_x, norm_obj, fwhm, med_sn2, obj_string, show_profile,
                    ind = None, l_limit = None, r_limit=None, xlim = None, xtrunc = 1e6):

    """
    Utility function to return Gaussian object profile in the case of low S/N ratio or too many rejected pixels.

    Args:
        sigma_x: ndarray (nspec, nspat)
            Spatial of gaussian
        norm_obj: ndarray (nspec, nspat)
            Normalized 2-d spectrum.
        fwhm: float
            FWHM parameter for Gaussian profile
        med_sn2: float
            Median (S/N)^2 used only for QA
        obj_string: str
            String identifying object. Used only for QA.
        show_profile: bool
            Is set, qa plot will be shown to screen
        ind: array, int
            Good indices in object profile. Used by QA routine.
        l_limit: float
            Left limit of profile fit where derivative is evaluated for Gaussian apodization. Used by QA routine.
        r_limit: float
            Right limit of profile fit where derivative is evaluated for Gaussian apodization. Used by QA routine.
        xlim: float
            Spatial location to trim object profile for plotting in QA routine.
        xtrunc: float
            Spatial nsigma to truncate object profile for plotting in QA routine.

    Returns:

    """

    profile_model = np.exp(-0.5*sigma_x**2)/np.sqrt(2.0 * np.pi)*(sigma_x ** 2 < 25.)
    info_string = "FWHM=" + "{:6.2f}".format(fwhm) + ", S/N=" + "{:8.3f}".format(np.sqrt(med_sn2))
    title_string = obj_string + ', ' + info_string
    msgs.info(title_string)
    inf = np.isfinite(profile_model) == False
    ninf = np.sum(inf)
    if ninf != 0:
        msgs.warn("Nan pixel values in object profile... setting them to zero")
        profile_model[inf] = 0.0
    if show_profile:
        qa_fit_profile(sigma_x, norm_obj, profile_model, title = title_string, l_limit = l_limit, r_limit = r_limit,
                       ind=ind, xlim = xlim, xtrunc=xtrunc)
    return profile_model


def fit_profile(image, ivar, waveimg, thismask, spat_img, trace_in, wave, flux, fluxivar,
                inmask=None, thisfwhm=4.0, max_trace_corr=2.0, sn_gauss=4.0, percentile_sn2=70.0,
                prof_nsigma=None, no_deriv=False, gauss=False, obj_string='',
                show_profile=False):

    """
    Fit a non-parametric object profile to an object spectrum, unless
    the S/N ratio is low (> sn_gauss) in which fit a simple Gaussian.
    Port of IDL LOWREDUX long_gprofile.pro

    Parameters
    ----------
    image : numpy float 2-d array (nspec, nspat)
        sky-subtracted image
    ivar : numpy float 2-d array (nspec, nspat)
        inverse variance of sky-subtracted image
    waveimg numpy float 2-d array (nspec, nspat)
        2-d wavelength map
    spat_img: float ndarray, shape (nspec, nspat)
        Image containing the spatial location of pixels. If not input,
        it will be computed via spat_img = np.outer(np.ones(nspec), np.arange(nspat))
    trace_in : numpy 1-d array (nspec,)
        object trace
    wave : numpy 1-d array (nspec,)
        extracted wavelength of spectrum
    flux : numpy 1-d array (nspec,)
        extracted flux of spectrum
    fluxivar : numpy 1-d array (nspec,)
        inverse variance of extracted flux spectrum
    thisfwhm : float, optional
        fwhm of the object trace
    max_trace_corr : float [default = 2.0], optional
        maximum trace correction to apply
    sn_gauss : float [default = 3.0], optional
        S/N ratio below which code just uses a Gaussian
    percentile_sn2: float [default = 70.0], optional
        Estimates the S/N of an object from pixels in the upper percentile_sn2 percentile of wavelength values.
        For example if percentile_sn2 = 70.0 then the upper 30% of spectrals are used.
        This ensures the code can still fit line only objects and/or high redshift quasars which might only have
        signal in reddest part of a spectrum.
    maskwidth : float [default = None], optional
        object maskwidth determined from object finding algorithm. If = None,
        code defaults to use 3.0*(np.max(thisfwhm) + 1.0)
        THIS PARAMETER IS NOT USED IN THIS METHOD
    prof_nsigma : float [default = None], optional
        Number of sigma to include in the profile fitting. This option is only needed for bright objects that are not
        point sources, which allows the profile fitting to fit the high S/N wings (rather than the default behavior
        which truncates exponentially). This allows for extracting all the flux and results in better sky-subtraction
        for bright extended objects.
    no_deriv : boolean [default = False]
        disables determination of derivatives and exponential apodization

    Returns
    -------
    sset: object
        bspline object
    outmask: : :class:`numpy.ndarray`
        output mask which the same size as xdata
    yfit  : :class:`numpy.ndarray`
        result of the bspline fit (same size as xdata)
    reduced_chi: float
        value of the reduced chi^2
    """

    if inmask is None:
        inmask = (ivar > 0.0) & thismask

    totmask = inmask & (ivar > 0.0) & thismask

    if prof_nsigma is not None:
        no_deriv = True


    thisfwhm = np.fmax(thisfwhm,1.0) # require the FWHM to be greater than 1 pixel

    nspat = image.shape[1]
    nspec = image.shape[0]

    # dspat is the spatial position along the image centered on the object trace
    dspat = spat_img - np.outer(trace_in, np.ones(nspat))
    # create some images we will need
    sn2_img = np.zeros((nspec,nspat))
    spline_img = np.zeros((nspec,nspat))

    flux_sm = scipy.ndimage.filters.median_filter(flux, size=5, mode = 'reflect')
    fluxivar_sm0 =  scipy.ndimage.filters.median_filter(fluxivar, size = 5, mode = 'reflect')
    fluxivar_sm0 = fluxivar_sm0*(fluxivar > 0.0)
    wave_min = waveimg[thismask].min()
    wave_max = waveimg[thismask].max()

    # sigma_x represents the profile argument, i.e. (x-x0)/sigma
    sigma = np.full(nspec, thisfwhm/2.3548)
    fwhmfit = sigma*2.3548
    trace_corr = np.zeros(nspec)
    sigma_x = dspat/np.outer(sigma, np.ones(nspat)) - np.outer(trace_corr, np.ones(nspat))

    # This adds an error floor to the fluxivar_sm, preventing too much rejection at high-S/N (i.e. standard stars)
    # TODO implement the ivar_cap function here in utils.
    sn_cap = 100.0
    fluxivar_sm = utils.clip_ivar(flux_sm, fluxivar_sm0, sn_cap)
    indsp = (wave >= wave_min) & (wave <= wave_max) & np.isfinite(flux_sm) & (flux_sm > -1000.0) & (fluxivar_sm > 0.0)
    eligible_pixels = np.sum((wave >= wave_min) & (wave <= wave_max))
    good_pix_frac = 0.05
    if np.sum(indsp) < good_pix_frac*eligible_pixels:
        msgs.warn('There are no pixels eligible to be fit for the object profile.' + msgs.newline() +
                  'There is likely an issue in local_skysub_extract. Returning a Gassuain with fwhm={:5.3f}'.format(thisfwhm))
        profile_model = return_gaussian(sigma_x, None, thisfwhm, 0.0, obj_string, False)
        return profile_model, trace_in, fwhmfit, 0.0

    b_answer, bmask   = fitting.iterfit(wave[indsp], flux_sm[indsp], invvar = fluxivar_sm[indsp],
                                        kwargs_bspline={'everyn': 1.5}, kwargs_reject={'groupbadpix':True,'maxrej':1})
    b_answer, bmask2  = fitting.iterfit(wave[indsp], flux_sm[indsp], invvar = fluxivar_sm[indsp]*bmask,
                                     kwargs_bspline={'everyn': 1.5}, kwargs_reject={'groupbadpix':True,'maxrej':1})
    c_answer, cmask   = fitting.iterfit(wave[indsp], flux_sm[indsp], invvar = fluxivar_sm[indsp]*bmask2,
                                     kwargs_bspline={'everyn': 30}, kwargs_reject={'groupbadpix':True,'maxrej':1})
    spline_flux, _ = b_answer.value(wave[indsp])
    try:
        cont_flux, _ = c_answer.value(wave[indsp])
    except:
        msgs.warn('Problem estimating S/N ratio of spectrum' + msgs.newline() +
                  'There is likely an issue in local_skysub_extract. Returning a Gassuain with fwhm={:5.3f}'.format(thisfwhm))
        profile_model = return_gaussian(sigma_x, None, thisfwhm, 0.0, obj_string, False)
        return profile_model, trace_in, fwhmfit, 0.0

    sn2 = (np.fmax(spline_flux*(np.sqrt(np.fmax(fluxivar_sm[indsp], 0))*bmask2),0))**2
    ind_nonzero = (sn2 > 0)
    nonzero = np.sum(ind_nonzero)
    if(nonzero >0):
        ## Select the top 30% data for estimating the med_sn2. This ensures the code still fit line only object and/or
        ## high redshift quasars which might only have signal in part of the spectrum.
        sn2_percentile = np.percentile(sn2,percentile_sn2)
        mean, med_sn2, stddev = stats.sigma_clipped_stats(sn2[sn2>sn2_percentile],sigma_lower=3.0,sigma_upper=5.0)
    else:
        med_sn2 = 0.0

    #min_wave = np.min(wave[indsp])
    #max_wave = np.max(wave[indsp])
    # TODO -- JFH document this
    spline_flux1 = np.zeros(nspec)
    cont_flux1 = np.zeros(nspec)
    sn2_1 = np.zeros(nspec)
    ispline = (wave >= wave_min) & (wave <= wave_max)
    spline_tmp, _ = b_answer.value(wave[ispline])
    spline_flux1[ispline] = spline_tmp
    cont_tmp, _ = c_answer.value(wave[ispline])
    cont_flux1[ispline] = cont_tmp
    isrt = np.argsort(wave[indsp])
    s2_1_interp = scipy.interpolate.interp1d(wave[indsp][isrt], sn2[isrt],assume_sorted=False, bounds_error=False,fill_value = 0.0)
    sn2_1[ispline] = s2_1_interp(wave[ispline])
    bmask = np.zeros(nspec,dtype='bool')
    bmask[indsp] = bmask2
    spline_flux1 = pydl.djs_maskinterp(spline_flux1,(bmask == False))
    cmask2 = np.zeros(nspec,dtype='bool')
    cmask2[indsp] = cmask
    cont_flux1 = pydl.djs_maskinterp(cont_flux1,(cmask2 == False))
    (_, _, sigma1) = stats.sigma_clipped_stats(flux[indsp],sigma_lower=3.0,sigma_upper=5.0)

    sn2_med_filt = scipy.ndimage.filters.median_filter(sn2, size=9, mode='reflect')
    if np.any(totmask):
        sn2_interp = scipy.interpolate.interp1d(wave[indsp][isrt],sn2_med_filt[isrt],assume_sorted=False,
                                                bounds_error=False,fill_value = 'extrapolate')
        sn2_img[totmask] = sn2_interp(waveimg[totmask])
    else:
        msgs.warn('All pixels are masked')

    msgs.info('sqrt(med(S/N)^2) = ' + "{:5.2f}".format(np.sqrt(med_sn2)))

    # TODO -- JFH document this
    if(med_sn2 <= 2.0):
        spline_img[totmask]= np.fmax(sigma1,0)
    else:
        if((med_sn2 <=5.0) and (med_sn2 > 2.0)):
            spline_flux1 = cont_flux1
        # Interp over points <= 0 in boxcar flux or masked points using cont model
        badpix = (spline_flux1 <= 0.5) | (bmask == False)
        goodval = (cont_flux1 > 0.0) & (cont_flux1 < 5e5)
        indbad1 = badpix & goodval
        nbad1 = np.sum(indbad1)
        if(nbad1 > 0):
            spline_flux1[indbad1] = cont_flux1[indbad1]
        indbad2 = badpix & np.invert(goodval)
        nbad2 = np.sum(indbad2)
        ngood0 = np.sum(np.invert(badpix))
        if((nbad2 > 0) or (ngood0 > 0)):
            spline_flux1[indbad2] = np.median(spline_flux1[np.invert(badpix)])
        # take a 5-pixel median to filter out some hot pixels
        spline_flux1 = scipy.ndimage.filters.median_filter(spline_flux1,size=5,mode ='reflect')

        # Create the normalized object image
        if np.any(totmask):
            igd = (wave >= wave_min) & (wave <= wave_max)
            isrt1 = np.argsort(wave[igd])
            #plt.plot(wave[igd][isrt1], spline_flux1[igd][isrt1])
            #plt.show()
            spline_img_interp = scipy.interpolate.interp1d(wave[igd][isrt1],spline_flux1[igd][isrt1],assume_sorted=False,
                                                           bounds_error=False,fill_value = 'extrapolate')
            spline_img[totmask] = spline_img_interp(waveimg[totmask])
        else:
            spline_img[totmask] = np.fmax(sigma1, 0)

    # TODO -- JFH document this
    norm_obj = (spline_img != 0.0)*image/(spline_img + (spline_img == 0.0))
    norm_ivar = ivar*spline_img**2

    # Cap very large inverse variances
    ivar_mask = (norm_obj > -0.2) & (norm_obj < 0.7) & totmask & np.isfinite(norm_obj) & np.isfinite(norm_ivar)
    norm_ivar = norm_ivar*ivar_mask
    good = (norm_ivar.flatten() > 0.0)
    ngood = np.sum(good)

    xtemp = (np.cumsum(np.outer(4.0 + np.sqrt(np.fmax(sn2_1, 0.0)),np.ones(nspat)))).reshape((nspec,nspat))
    xtemp = xtemp/xtemp.max()

    msgs.info("Gaussian vs b-spline of width " + "{:6.2f}".format(thisfwhm) + " pixels")
    area = 1.0

    # If we have too few pixels to fit a profile or S/N is too low, just use a Gaussian profile
    if((ngood < 10) or (med_sn2 < sn_gauss**2) or (gauss is True)):
        msgs.info("Too few good pixels or S/N <" + "{:5.1f}".format(sn_gauss) + " or gauss flag set")
        msgs.info("Returning Gaussian profile")
        profile_model = return_gaussian(sigma_x, norm_obj, thisfwhm, med_sn2, obj_string,show_profile,ind=good,xtrunc=7.0)
        return profile_model, trace_in, fwhmfit, med_sn2

    mask = np.full(nspec*nspat, False, dtype=bool)

    # The following lines set the limits for the b-spline fit
    limit = scipy.special.erfcinv(0.1/np.sqrt(med_sn2))*np.sqrt(2.0)
    if(prof_nsigma is None):
        sinh_space = 0.25*np.log10(np.fmax((1000./np.sqrt(med_sn2)),10.))
        abs_sigma = np.fmin((np.abs(sigma_x.flat[good])).max(),2.0*limit)
        min_sigma = np.fmax(sigma_x.flat[good].min(), (-abs_sigma))
        max_sigma = np.fmin(sigma_x.flat[good].max(), (abs_sigma))
        nb = (np.arcsinh(abs_sigma)/sinh_space).astype(int) + 1
    else:
        msgs.info("Using prof_nsigma= " + "{:6.2f}".format(prof_nsigma) + " for extended/bright objects")
        nb = np.round(prof_nsigma > 10)
        max_sigma = prof_nsigma
        min_sigma = -1*prof_nsigma
        sinh_space = np.arcsinh(prof_nsigma)/nb

    rb = np.sinh((np.arange(nb) + 0.5) * sinh_space)
    bkpt = np.concatenate([(-rb)[::-1], rb])
    keep = ((bkpt >= min_sigma) & (bkpt <= max_sigma))
    bkpt = bkpt[keep]

    # Attempt B-spline first
    GOOD_PIX = (sn2_img > sn_gauss**2) & (norm_ivar > 0)
    IN_PIX   = (sigma_x >= min_sigma) & (sigma_x <= max_sigma) & (norm_ivar > 0)
    ngoodpix = np.sum(GOOD_PIX)
    ninpix     = np.sum(IN_PIX)

    if (ngoodpix >= 0.2*ninpix):
        inside,  = np.where((GOOD_PIX & IN_PIX).flatten())
    else:
        inside, = np.where(IN_PIX.flatten())


    si = inside[np.argsort(sigma_x.flat[inside])]
    sr = si[::-1]

    bset, bmask = fitting.iterfit(sigma_x.flat[si],norm_obj.flat[si], invvar = norm_ivar.flat[si],
                                      nord = 4, bkpt = bkpt, maxiter = 15, upper = 1, lower = 1)
    mode_fit, _ = bset.value(sigma_x.flat[si])
    median_fit = np.median(norm_obj[norm_ivar > 0.0])

    # TODO I don't follow the logic behind this statement but I'm leaving it for now. If the median is large it is used, otherwise we  user zero???
    if (np.abs(median_fit) > 0.01):
        msgs.info("Median flux level in profile is not zero: median = " + "{:7.4f}".format(median_fit))
    else:
        median_fit = 0.0

    # Find the peak and FWHM based this profile fit
    (peak, peak_x, lwhm, rwhm) = findfwhm(mode_fit - median_fit, sigma_x.flat[si])
    trace_corr = np.full(nspec, peak_x)
    min_level = peak*np.exp(-0.5*limit**2)

    bspline_fwhm = (rwhm - lwhm)*thisfwhm/2.3548
    msgs.info("Bspline FWHM: " + "{:7.4f}".format(bspline_fwhm) + ", compared to initial object finding FWHM: " + "{:7.4f}".format(thisfwhm) )
    sigma = sigma * (rwhm-lwhm)/2.3548

    limit = limit * (rwhm-lwhm)/2.3548

    rev_fit = mode_fit[::-1]
    lind, = np.where(((rev_fit < (min_level+median_fit)) & (sigma_x.flat[sr] < peak_x)) | (sigma_x.flat[sr] < (peak_x-limit)))
    if (lind.size > 0):
        lp = lind.min()
        l_limit = sigma_x.flat[sr[lp]]
    else:
        l_limit = min_sigma

    rind, = np.where(((mode_fit < (min_level+median_fit)) & (sigma_x.flat[si] > peak_x)) | (sigma_x.flat[si] > (peak_x+limit)))
    if (rind.size > 0):
        rp = rind.min()
        r_limit = sigma_x.flat[si[rp]]
    else:
        r_limit = max_sigma

    msgs.info("Trace limits: limit = " + "{:7.4f}".format(limit) + ", min_level = " + "{:7.4f}".format(min_level) +
              ", l_limit = " + "{:7.4f}".format(l_limit) + ", r_limit = " + "{:7.4f}".format(r_limit))

    # Just grab the data points within the limits
    mask[si]=((norm_ivar.flat[si] > 0) & (np.abs(norm_obj.flat[si] - mode_fit) < 0.1))
    inside, = np.where((sigma_x.flat[si] > l_limit) & (sigma_x.flat[si] < r_limit) & mask[si])
    ninside = inside.size

    # If we have too few pixels after this step, then again just use a Gaussian profile and return.
    if(ninside < 10):
        msgs.info("Too few pixels inside l_limit and r_limit")
        msgs.info("Returning Gaussian profile")
        profile_model = return_gaussian(sigma_x, norm_obj, bspline_fwhm, med_sn2, obj_string,show_profile,
                                                          ind=good, l_limit=l_limit, r_limit=r_limit, xlim=7.0)
        return (profile_model, trace_in, fwhmfit, med_sn2)

    sigma_iter = 3
    isort = (xtemp.flat[si[inside]]).argsort()
    inside = si[inside[isort]]
    pb = np.ones(inside.size)

    for iiter in range(1,sigma_iter + 1):
        mode_zero, _ = bset.value(sigma_x.flat[inside])
        mode_zero = mode_zero*pb

        mode_min05, _ = bset.value(sigma_x.flat[inside]-0.5)
        mode_plu05, _ = bset.value(sigma_x.flat[inside]+0.5)
        mode_shift = (mode_min05  - mode_plu05)*pb*((sigma_x.flat[inside] > (l_limit + 0.5)) &
                                                (sigma_x.flat[inside] < (r_limit - 0.5)))

        mode_by13, _ = bset.value(sigma_x.flat[inside]/1.3)
        mode_stretch = mode_by13*pb/1.3 - mode_zero

        nbkpts = (np.log10(np.fmax(med_sn2, 11.0))).astype(int)

        xx = np.sum(xtemp, 1)/nspat
        profile_basis = np.column_stack((mode_zero,mode_shift))

        mode_shift_out = fitting.bspline_profile(xtemp.flat[inside], norm_obj.flat[inside],
                                               norm_ivar.flat[inside], profile_basis,
                                               maxiter=1, kwargs_bspline={'nbkpts':nbkpts})
        # Check to see if the mode fit failed, if so punt and return a Gaussian
        if not np.any(mode_shift_out[1]):
            msgs.info('B-spline fit to trace correction failed for fit to ninside = {:}'.format(ninside) + ' pixels')
            msgs.info("Returning Gaussian profile")
            profile_model = return_gaussian(sigma_x, norm_obj, bspline_fwhm, med_sn2, obj_string,
                                            show_profile, ind=good, l_limit=l_limit, r_limit=r_limit, xlim=7.0)
            return (profile_model, trace_in, fwhmfit, med_sn2)


        mode_shift_set = mode_shift_out[0]
        temp_set = bspline.bspline(None, fullbkpt=mode_shift_set.breakpoints,
                                   nord=mode_shift_set.nord)
        temp_set.coeff = mode_shift_set.coeff[0, :]
        h0, _ = temp_set.value(xx)
        temp_set.coeff = mode_shift_set.coeff[1, :]
        h1, _ = temp_set.value(xx)
        ratio_10 = (h1/(h0 + (h0 == 0.0)))
        delta_trace_corr = ratio_10/(1.0 + np.abs(ratio_10)/0.1)
        trace_corr = trace_corr + delta_trace_corr

        profile_basis = np.column_stack((mode_zero,mode_stretch))
        mode_stretch_out = fitting.bspline_profile(xtemp.flat[inside], norm_obj.flat[inside],
                                                 norm_ivar.flat[inside], profile_basis, maxiter=1,
                                                 fullbkpt=mode_shift_set.breakpoints)
        if not np.any(mode_stretch_out[1]):
            msgs.info('B-spline fit to width correction failed for fit to ninside = {:}'.format(ninside) + ' pixels')
            msgs.info("Returning Gaussian profile")
            profile_model  = return_gaussian(sigma_x, norm_obj, bspline_fwhm, med_sn2, obj_string,
                                                              show_profile,ind=good, l_limit=l_limit, r_limit=r_limit, xlim=7.0)
            return (profile_model, trace_in, fwhmfit, med_sn2)

        mode_stretch_set = mode_stretch_out[0]
        temp_set = bspline.bspline(None, fullbkpt=mode_stretch_set.breakpoints,
                                   nord=mode_stretch_set.nord)
        temp_set.coeff = mode_stretch_set.coeff[0, :]
        h0, _ = temp_set.value(xx)
        temp_set.coeff = mode_stretch_set.coeff[1, :]
        h2, _ = temp_set.value(xx)
        h0 = np.fmax(h0 + h2*mode_stretch.sum()/mode_zero.sum(),0.1)
        ratio_20 = (h2 / (h0 + (h0 == 0.0)))
        sigma_factor = 0.3 * ratio_20 / (1.0 + np.abs(ratio_20))

        msgs.info("Iteration# " + "{:3d}".format(iiter))
        msgs.info("Median abs value of trace correction = " + "{:8.3f}".format(np.median(np.abs(delta_trace_corr))))
        msgs.info("Median abs value of width correction = " + "{:8.3f}".format(np.median(np.abs(sigma_factor))))

        sigma = sigma*(1.0 + sigma_factor)
        area = area * h0/(1.0 + sigma_factor)

        sigma_x = dspat/np.outer(sigma, np.ones(nspat)) - np.outer(trace_corr, np.ones(nspat))

        # Update the profile B-spline fit for the next iteration
        if iiter < sigma_iter-1:
            ss = sigma_x.flat[inside].argsort()
            pb = (np.outer(area, np.ones(nspat,dtype=float))).flat[inside]
            keep = (bkpt >= sigma_x.flat[inside].min()) & (bkpt <= sigma_x.flat[inside].max())
            if keep.sum() == 0:
                keep = np.ones(bkpt.size, dtype=bool)
            bset_out = fitting.bspline_profile(sigma_x.flat[inside[ss]], norm_obj.flat[inside[ss]],
                                             norm_ivar.flat[inside[ss]],pb[ss], nord=4,
                                             bkpt=bkpt[keep], maxiter=2)
            if not np.any(bset_out[1]):
                msgs.info('B-spline to profile in trace and width correction loop failed for fit to ninside = {:}'.format(ninside) + ' pixels')
                msgs.info("Returning Gaussian profile")
                profile_model = return_gaussian(sigma_x, norm_obj, bspline_fwhm, med_sn2, obj_string,
                                                                  show_profile, ind=good, l_limit=l_limit,r_limit=r_limit, xlim=7.0)
                return (profile_model, trace_in, fwhmfit, med_sn2)

            bset = bset_out[0] # This updated bset used for the next set of trace corrections

    # Apply trace corrections only if they are small (added by JFH)
    if np.median(np.abs(trace_corr*sigma)) < max_trace_corr:
        xnew = trace_corr * sigma + trace_in
    else:
        xnew = trace_in

    fwhmfit = sigma*2.3548
    ss=sigma_x.flatten().argsort()
    inside, = np.where((sigma_x.flat[ss] >= min_sigma) &
                       (sigma_x.flat[ss] <= max_sigma) &
                       mask[ss] &
                       np.isfinite(norm_obj.flat[ss]) &
                       np.isfinite(norm_ivar.flat[ss]))
    pb = (np.outer(area, np.ones(nspat,dtype=float)))
    bset_out = fitting.bspline_profile(sigma_x.flat[ss[inside]], norm_obj.flat[ss[inside]],
                                     norm_ivar.flat[ss[inside]], pb.flat[ss[inside]], nord=4,
                                     bkpt=bkpt, upper=10, lower=10)
    bset = bset_out[0]
    outmask = bset_out[1]

    # igood = False for pixels within (min_sigma, max_sigma), True outside
    igood = (sigma_x.flatten() > min_sigma) & (sigma_x.flatten() < max_sigma)
    full_bsp = np.zeros(nspec*nspat, dtype=float)
    sigma_x_igood = sigma_x.flat[igood]
    yfit_out, _  = bset.value(sigma_x_igood)
    full_bsp[igood] = yfit_out
    isrt2 = sigma_x_igood.argsort()
    (peak, peak_x, lwhm, rwhm) = findfwhm(yfit_out[isrt2] - median_fit, sigma_x_igood[isrt2])


    left_bool = (((full_bsp[ss] < (min_level+median_fit)) & (sigma_x.flat[ss] < peak_x)) | (sigma_x.flat[ss] < (peak_x-limit)))[::-1]
    ind_left, = np.where(left_bool)
    if ind_left.size == 0:
        lp = 0
    else:
        lp = np.fmax(ind_left.min(), 0)
    righ_bool = ((full_bsp[ss] < (min_level+median_fit)) & (sigma_x.flat[ss] > peak_x))  | (sigma_x.flat[ss] > (peak_x+limit))
    ind_righ, = np.where(righ_bool)
    if ind_righ.size == 0:
        rp = 0
    else:
        rp = np.fmax(ind_righ.min(), 0)

    l_limit = ((sigma_x.flat[ss])[::-1])[lp] - 0.1
    r_limit = sigma_x.flat[ss[rp]] + 0.1

    # Determine the left and right locations (l_limit and r_limit) where the profile logarithmic derivative crosses 1.0 for apodization
    # of the object profiles. If the profile slope is never actually one, then just find the maximum value in the interval between
    # (l_limit, -1.0) or (1.0, r_limit)
    l_lim_vec = np.arange(l_limit+0.1,-1.0, 0.1)
    l_lim_vec = np.asarray([-1.1]) if len(l_lim_vec) == 0 else l_lim_vec
    l_fit1, _ = bset.value(l_lim_vec)
    l_fit2, _ = bset.value(l_lim_vec*0.9)
    l_deriv_vec = (np.log(l_fit2) - np.log(l_fit1))/(0.1*l_lim_vec)
    l_deriv_max = np.fmax(l_deriv_vec.min(), -1.0)

    r_lim_vec = np.arange(r_limit-0.1,1.0, -0.1)
    r_lim_vec = np.asarray([1.1]) if len(r_lim_vec) == 0 else r_lim_vec

    r_fit1, _ = bset.value(r_lim_vec)
    r_fit2, _ = bset.value(r_lim_vec*0.9)
    r_deriv_vec = (np.log(r_fit2) - np.log(r_fit1))/(0.1*r_lim_vec)
    r_deriv_max = np.fmin(r_deriv_vec.max(), 1.0)


    while True:
        l_limit += 0.1
        l_fit, _ = bset.value(np.asarray([l_limit]))
        l2, _ = bset.value(np.asarray([l_limit])* 0.9)
        l_deriv = (np.log(l2[0]) - np.log(l_fit[0]))/(0.1*l_limit)
        if (l_deriv <= l_deriv_max) | (l_limit >= -1.0):
            break

    while True:
        r_limit -= 0.1
        r_fit, _ = bset.value(np.asarray([r_limit]))
        r2, _ = bset.value(np.asarray([r_limit])* 0.9)
        r_deriv = (np.log(r2[0]) - np.log(r_fit[0]))/(0.1*r_limit)
        if (r_deriv >= r_deriv_max) | (r_limit <= 1.0):
            break

    # JXP kludge
    if prof_nsigma is not None:
       #By setting them to zero we ensure QA won't plot them in the profile QA.
       l_limit = 0.0
       r_limit = 0.0

    # Apodization of object profiles with exponential, ensuring continuity of first derivative
    if (l_deriv < 0) and (r_deriv > 0) and no_deriv is False:
        left = sigma_x.flatten() < l_limit
        full_bsp[left] =  np.exp(-(sigma_x.flat[left]-l_limit)*l_deriv) * l_fit[0]
        right = sigma_x.flatten() > r_limit
        full_bsp[right] = np.exp(-(sigma_x.flat[right] - r_limit) * r_deriv) * r_fit[0]

    # Final object profile
    full_bsp = full_bsp.reshape(nspec,nspat)
    profile_model = full_bsp*pb
    res_mode = (norm_obj.flat[ss[inside]] - profile_model.flat[ss[inside]])*np.sqrt(norm_ivar.flat[ss[inside]])
    chi_good = (outmask == True) & (norm_ivar.flat[ss[inside]] > 0)
    chi_med = np.median(res_mode[chi_good]**2)
    chi_zero = np.median(norm_obj.flat[ss[inside]]**2*norm_ivar.flat[ss[inside]])

    msgs.info("--------------------  Results of Profile Fit --------------------")
    msgs.info(" min(fwhmfit)={:5.2f}".format(fwhmfit.min()) +
              " max(fwhmfit)={:5.2f}".format(fwhmfit.max()) + " median(chi^2)={:5.2f}".format(chi_med) +
              " nbkpts={:2d}".format(bkpt.size))
    msgs.info("-----------------------------------------------------------------")

    nxinf = np.sum(np.isfinite(xnew) == False)
    if (nxinf != 0):
        msgs.warn("Nan pixel values in trace correction")
        msgs.warn("Returning original trace....")
        xnew = trace_in
    inf = np.isfinite(profile_model) == False
    ninf = np.sum(inf)
    if (ninf != 0):
        msgs.warn("Nan pixel values in object profile... setting them to zero")
        profile_model[inf] = 0.0
    # Normalize profile
    norm = np.outer(np.sum(profile_model, 1), np.ones(nspat))
    if (np.sum(norm) > 0.0):
        profile_model = (norm > 0.0)*profile_model/(norm + (norm == 0.0))

    info_string = "FWHM range:" + "{:5.2f}".format(fwhmfit.min()) + " - {:5.2f}".format(fwhmfit.max()) \
                  + ", S/N=" + "{:8.3f}".format(np.sqrt(med_sn2)) + ", median(chi^2)={:8.3f}".format(chi_med)
    title_string = obj_string + ' ' + info_string
    if(show_profile):
        qa_fit_profile(sigma_x, norm_obj/(pb + (pb == 0.0)), full_bsp,
                       l_limit = l_limit, r_limit = r_limit, ind = ss[inside], xlim = prof_nsigma, title = title_string)

    return (profile_model, xnew, fwhmfit, med_sn2)



def remap_orders(xinit, spec_min_max, inverse=False):

    """
    This code remaps echelle orders to all extend over the same numer of pixels. It is meant to deal with cases
    where the echelle orders do not completely span the detector. Initial tests with PCA using this remapping
    indicate (for vlt_xshooter_nir) that the remapping does not work as well as simply linear extroplating the traces
    in the PCA. The traces appear to compress better onto a PCA when they are not remapped. So this functionality
    is experimemtnal and not currently used.


    Args:
        xinit: ndarray, (nspec, norders)
           Array of input traces that one wants to remap.
        spec_min_max: ndarray, (2, norders)
           Array containing the spec_min and spec_max defined for each order
        inverse: bool, default = False,
           If True, perform the inverse re-mapping rather than the re-mapping.

    Returns:

    """

    nspec, norders = xinit.shape
    spec_vec = np.arange(nspec)
    spec_vec_norm = spec_vec/float(nspec-1)
    xinit_remap = np.zeros_like(xinit)
    for iord in range(norders):
        igood = (spec_vec >= spec_min_max[0,iord]) & (spec_vec <= spec_min_max[1,iord])
        spec_norm_iord = (spec_vec - spec_vec[igood].min())/(spec_vec[igood].max()  - spec_vec[igood].min())
        if inverse:
            xinit_remap[:, iord] = scipy.interpolate.interp1d(spec_vec_norm, xinit[:, iord], kind='linear',
                                                              bounds_error=False, fill_value='extrapolate')(spec_norm_iord)
        else:
            xinit_remap[:,iord] = scipy.interpolate.interp1d(spec_norm_iord[igood], xinit[igood,iord], kind='linear',
                                                             bounds_error=False, fill_value='extrapolate')(spec_vec_norm)
    return xinit_remap

