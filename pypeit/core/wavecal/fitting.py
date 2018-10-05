""" Module for finding patterns in arc line spectra
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)


import numpy as np

from pypeit import utils
from pypeit.core.wavecal import qa

from pypeit import debugger

def iterative_fitting(spec, tcent, ifit, IDs, llist, disp, plot_fil=None,
                      verbose=False, aparm=None, weights=None):

    if aparm is None:
        aparm = dict(llist='',
                    disp=disp,           # Ang/unbinned pixel
                    match_toler=3.,      # Matcing tolerance (pixels)
                    func='legendre',     # Function for fitting
                    n_first=2,           # Order of polynomial for first fit
                    n_final=4,           # Order of polynomial for final fit
                    nsig_rej=2.,         # Number of sigma for rejection
                    nsig_rej_final=3.0)  # Number of sigma for rejection (final fit)
            #  disp_toler=0.1,      # 10% tolerance  JFH disp_toler is never used.

    if weights is None:
        weights = np.ones(tcent.size)

    npix = spec.size
    aparm['disp'] = disp

    # Setup for fitting
    sv_ifit = list(ifit)  # Keep the originals
    all_ids = -999.*np.ones(len(tcent))
    all_idsion = np.array(['UNKNWN']*len(tcent))
    all_ids[ifit] = IDs

    # Fit
    n_order = aparm['n_first']
    flg_quit = False
    fmin, fmax = -1., 1.
    while (n_order <= aparm['n_final']) and (flg_quit is False):
        # Fit with rejection
        xfit, yfit, wfit = tcent[ifit], all_ids[ifit], weights[ifit]
        mask, fit = utils.robust_polyfit(xfit, yfit, n_order, function=aparm['func'], sigma=aparm['nsig_rej'],
                                         minv=fmin, maxv=fmax, verbose=verbose, weights=wfit)

        rms_ang = utils.calc_fit_rms(xfit[mask == 0], yfit[mask == 0], fit,
                                     aparm['func'], minv=fmin, maxv=fmax, weights=wfit[mask == 0])
        rms_pix = rms_ang/disp
        if verbose:
            print("RMS = {:g}".format(rms_pix))

        # Reject but keep originals (until final fit)
        ifit = list(ifit[mask == 0]) + sv_ifit
        # Find new points (should we allow removal of the originals?)
        twave = utils.func_val(fit, tcent, aparm['func'], minv=fmin, maxv=fmax)
        for ss, iwave in enumerate(twave):
            mn = np.min(np.abs(iwave-llist['wave']))
            if mn/aparm['disp'] < aparm['match_toler']:
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
        if n_order < (aparm['n_final']+2):
            n_order += 1
        else:
            # This does 2 iterations at the final order
            flg_quit = True

    # Final fit (originals can now be rejected)
    fmin, fmax = 0., 1.
    xfit, yfit, wfit = tcent[ifit]/(npix-1), all_ids[ifit], weights[ifit]
    mask, fit = utils.robust_polyfit(xfit, yfit, n_order, function=aparm['func'], sigma=aparm['nsig_rej_final'],
                                     minv=fmin, maxv=fmax, verbose=verbose, weights=wfit)#, debug=True)
    irej = np.where(mask == 1)[0]
    if len(irej) > 0:
        xrej = xfit[irej]
        yrej = yfit[irej]
        if verbose:
            for kk, imask in enumerate(irej):
                wave = utils.func_val(fit, xrej[kk], aparm['func'], minv=fmin, maxv=fmax)
                print('Rejecting arc line {:g}; {:g}'.format(yfit[imask], wave))
    else:
        xrej = []
        yrej = []
    xfit = xfit[mask == 0]
    yfit = yfit[mask == 0]
    wfit = wfit[mask == 0]
    ions = all_idsion[ifit][mask == 0]
    # Final RMS
    rms_ang = utils.calc_fit_rms(xfit, yfit, fit, aparm['func'],
                                 minv=fmin, maxv=fmax, weights=wfit)
    rms_pix = rms_ang/disp

    # Pack up fit
    final_fit = dict(fitc=fit, function=aparm['func'], xfit=xfit, yfit=yfit, weights=wfit,
                     ions=ions, fmin=fmin, fmax=fmax, xnorm=float(npix),
                     xrej=xrej, yrej=yrej, mask=mask, spec=spec, nrej=aparm['nsig_rej_final'],
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
