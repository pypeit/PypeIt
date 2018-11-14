""" Module for finding patterns in arc line spectra
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)


import numpy as np

from pypeit import utils
from pypeit.core.wavecal import qa
from pypeit import msgs


def fit_slit(spec, patt_dict, tcent, line_lists, outroot=None, slittxt="Slit", thar=False,match_toler=3.0,
             func='legendre', n_first=2,sigrej_first=2.0,n_final=4,sigrej_final=3.0,verbose=False):

    """ Perform a fit to the wavelength solution. Wrapper for iterative fitting code.

    Parameters
    ----------
    spec : ndarray
      arc spectrum
    patt_dict : dict
      dictionary of patterns
    tcent: ndarray
      List of the detections in this slit to be fit using the patt_dict
    outroot : str
      root directory to save QA

    Optional Parameters
    -------------------
    outroot: str
      Path for QA file.
    slittxt : str
      Label used for QA
    thar: bool, default = False
      True if this is a ThAr fit
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
    verbose : bool
      If True, print out more information.
    plot_fil:
      Filename for plotting some QA?

    Returns
    -------
    final_fit : dict
      A dictionary containing all of the information about the fit
    """

    # Check that patt_dict and tcent refer to each other
    if patt_dict['mask'].shape != tcent.shape:
        msgs.error('patt_dict and tcent do not refer to each other. Something is very wrong')

    # Perform final fit to the line IDs
    if thar:
        NIST_lines = (line_lists['NIST'] > 0) & (np.char.find(line_lists['Source'].data, 'MURPHY') >= 0)
    else:
        NIST_lines = line_lists['NIST'] > 0
    ifit = np.where(patt_dict['mask'])[0]

    if outroot is not None:
        plot_fil = outroot + slittxt + '_fit.pdf'
    else:
        plot_fil = None

    # TODO Profx maybe you can add a comment on what this is doing. Why do we have use_unknowns=True only to purge them later??
    # Purge UNKNOWNS from ifit
    imsk = np.ones(len(ifit), dtype=np.bool)
    for kk, idwv in enumerate(np.array(patt_dict['IDs'])[ifit]):
        if np.min(np.abs(line_lists['wave'][NIST_lines] - idwv)) > 0.01:
            imsk[kk] = False
    ifit = ifit[imsk]
    # Fit
    try:
        final_fit = iterative_fitting(spec, tcent, ifit,np.array(patt_dict['IDs'])[ifit], line_lists[NIST_lines],
                                      patt_dict['bdisp'],match_toler=match_toler, func=func, n_first=n_first,
                                      sigrej_first=sigrej_first,n_final=n_final, sigrej_final=sigrej_final,
                                      plot_fil=plot_fil, verbose=verbose)
    except TypeError:
        # A poor fitting result, this can be ignored.
        return None

    if plot_fil is not None:
        print("Wrote: {:s}".format(plot_fil))

    # Return
    return final_fit


def iterative_fitting(spec, tcent, ifit, IDs, llist, disp,
                      match_toler = 2.0, func = 'legendre', n_first = 2, sigrej_first = 2.0, n_final = 4, sigrej_final = 3.0,
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
      Dictionary containing the full fitting results and the final best guess of the line IDs
    """

    if weights is None:
        weights = np.ones(tcent.size)

    nspec = spec.size
    # Setup for fitting
    sv_ifit = list(ifit)  # Keep the originals
    all_ids = -999.*np.ones(len(tcent))
    all_idsion = np.array(['UNKNWN']*len(tcent))
    all_ids[ifit] = IDs

    # Fit
    n_order = n_first
    flg_quit = False
    #fmin , fmax = -1.0, 1.0
    fmin, fmax = 0.0, float(nspec-1)
    while (n_order <= n_final) and (flg_quit is False):
        # Fit with rejection
        xfit, yfit, wfit = tcent[ifit], all_ids[ifit], weights[ifit]
        mask, fit = utils.robust_polyfit(xfit, yfit, n_order, function=func, sigma=sigrej_first,
                                         minv=fmin, maxv=fmax, verbose=verbose, weights=wfit)

        rms_ang = utils.calc_fit_rms(xfit[mask == 0], yfit[mask == 0], fit, func, minv=fmin, maxv=fmax,
                                     weights=wfit[mask == 0])
        rms_pix = rms_ang/disp
        if verbose:
            msgs.info('n_order = {:d}'.format(n_order) + ': RMS = {:g}'.format(rms_pix))

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
    #fmin, fmax = 0., 1.
    #xfit, yfit, wfit = tcent[ifit]/(nspec-1), all_ids[ifit], weights[ifit]
    xfit, yfit, wfit = tcent[ifit], all_ids[ifit], weights[ifit]
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
    cen_wave = utils.func_val(fit, float(nspec)/2, func, minv=fmin, maxv=fmax)
    cen_wave_min1 = utils.func_val(fit, float(nspec)/2 - 1.0, func, minv=fmin, maxv=fmax)
    cen_disp = cen_wave - cen_wave_min1

    final_fit = dict(fitc=fit, function=func, pixels_fit=xfit, wave_fit=yfit, weights=wfit, ions=ions,
                     fmin=fmin, fmax=fmax, nspec=nspec, cen_wave = cen_wave, cen_disp = cen_disp,
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
