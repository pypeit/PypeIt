""" Module for finding patterns in arc line spectra
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)


import numpy as np

from pypeit import utils
from pypeit.core.wavecal import qa
from pypeit import msgs

from pypeit import debugger

def iterative_fitting(spec, tcent, ifit, IDs, llist, disp,
                      match_toler = 3.0, func = 'legendre', n_first = 2, sigrej_first = 2.0, n_final = 4, sigrej_final = 3.0,
                      weights=None, plot_fil=None, verbose = False):

    """ Routine for iteratively fitting wavelength solutions.

    Parameters
    ----------
    spec : ndarray, shape = (nspec,)
      arcline spectrum
    tcent : ndarray
      Centroids in pixels of lines identified in spec
    ifit : ndarray
      Indices of the lines that will be fit
    IDs: ndarray
      wavelength IDs of the lines that will be fit (I think?)
    llist: dict
      Linelist dictionary
    disp: float
      dispersion

    Optional Parameters
    -------------------
    match_toler: float, default = 3.0
      Matching tolerance when searching for new lines. This is the difference in pixels between the wavlength assigned to
      an arc line by an iteration of the wavelength solution to the wavelength in the line list.
    func: str, default = 'legendre'
      Name of function used for the wavelength solution
    n_first: int, default = 2
      Order of first guess to the wavelength solution.
    sigrej_first: float, default = 2.0
      Number of sigma for rejection for the first guess to the wavelength solution.
    n_final: int, default = 4
      Order of the final wavelength solution fit
    sigrej_final: float, default = 3.0
      Number of sigma for rejection for the final fit to the wavelength solution.
    weights: ndarray
      Weights to be used?
    verbose : bool
      If True, print out more information.
    plot_fil:
      Filename for plotting some QA?

    Returns
    -------
    final_fit: dict
      Dictionaty containing the full fitting results and the final best guess of the line IDs
    """

    from IPython import embed
    embed()
    if weights is None:
        weights = np.ones(tcent.size)

    npix = spec.size
    # Setup for fitting
    sv_ifit = list(ifit)  # Keep the originals
    all_ids = -999.*np.ones(len(tcent))
    all_idsion = np.array(['UNKNWN']*len(tcent))
    all_ids[ifit] = IDs

    # Fit
    n_order = n_first
    flg_quit = False
    fmin, fmax = -1., 1.
    while (n_order <= n_final) and (flg_quit is False):
        # Fit with rejection
        xfit, yfit, wfit = tcent[ifit], all_ids[ifit], weights[ifit]
        mask, fit = utils.robust_polyfit(xfit, yfit, n_order, function=func, sigma=sigrej_first,
                                         minv=fmin, maxv=fmax, verbose=verbose, weights=wfit)

        rms_ang = utils.calc_fit_rms(xfit[mask == 0], yfit[mask == 0], fit, func, minv=fmin, maxv=fmax,
                                     weights=wfit[mask == 0])
        rms_pix = rms_ang/disp
        if verbose:
            msgs.info("RMS = {:g}".format(rms_pix))

        # Reject but keep originals (until final fit)
        ifit = list(ifit[mask == 0]) + sv_ifit
        # Find new points (should we allow removal of the originals?)
        twave = utils.func_val(fit, tcent, func, minv=fmin, maxv=fmax)
        for ss, iwave in enumerate(twave):
            mn = np.min(np.abs(iwave-llist['wave']))
            if mn/disp < match_toler:
                imn = np.argmin(np.abs(iwave-llist['wave']))
                #if verbose:
                #    print('Adding {:g} at {:g}'.format(llist['wave'][imn],tcent[ss]))
                # Update and append
                all_ids[ss] = llist['wave'][imn]
                all_idsion[ss] = llist['ion'][imn]
                ifit.append(ss)
        # Keep unique ones
        ifit = np.unique(np.array(ifit, dtype=int))
        # Increment order
        if n_order < (n_final+2):
            n_order += 1
        else:
            # This does 2 iterations at the final order
            flg_quit = True

    # Final fit (originals can now be rejected)
    fmin, fmax = 0., 1.
    xfit, yfit, wfit = tcent[ifit]/(npix-1), all_ids[ifit], weights[ifit]
    mask, fit = utils.robust_polyfit(xfit, yfit, n_order, function=func, sigma=sigrej_final,
                                     minv=fmin, maxv=fmax, verbose=verbose, weights=wfit)#, debug=True)
    irej = np.where(mask == 1)[0]
    if len(irej) > 0:
        xrej = xfit[irej]
        yrej = yfit[irej]
        if verbose:
            for kk, imask in enumerate(irej):
                wave = utils.func_val(fit, xrej[kk], func, minv=fmin, maxv=fmax)
                msgs.info('Rejecting arc line {:g}; {:g}'.format(yfit[imask], wave))
    else:
        xrej = []
        yrej = []
    xfit = xfit[mask == 0]
    yfit = yfit[mask == 0]
    wfit = wfit[mask == 0]
    ions = all_idsion[ifit][mask == 0]
    # Final RMS
    rms_ang = utils.calc_fit_rms(xfit, yfit, fit, func,
                                 minv=fmin, maxv=fmax, weights=wfit)
    rms_pix = rms_ang/disp

    # Pack up fit
    final_fit = dict(fitc=fit, function=func, xfit=xfit, yfit=yfit, weights=wfit,
                     ions=ions, fmin=fmin, fmax=fmax, xnorm=float(npix),
                     xrej=xrej, yrej=yrej, mask=mask, spec=spec, nrej=sigrej_final,
                     shift=0., tcent=tcent, rms=rms_pix)

    # If set to True, this will output a file that can then be included in the tests
    saveit = False
    if saveit:
        from linetools import utils as ltu
        jdict = ltu.jsonify(final_fit)
        if plot_fil is None:
            outname = "temp"
            print("You should have set the plot_fil directory to save wavelength fits... using 'temp' as a filename")
        else:
            outname = plot_fil
        ltu.savejson(outname + '.json', jdict, easy_to_read=True, overwrite=True)
        print(" Wrote: {:s}".format(outname + '.json'))

    # QA
    if plot_fil is not None:
        qa.arc_fit_qa(final_fit, plot_fil)
    # Return
    return final_fit
