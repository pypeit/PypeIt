""" Module for finding patterns in arc line spectra
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

from scipy.ndimage.filters import gaussian_filter
import numpy as np
import pdb

from pypeit.core.wavecal import waveio
from pypeit.core.wavecal import patterns
from pypeit.core.wavecal import fitting
from pypeit.core.wavecal import utils
from pypeit.core.wavecal import qa

from pypeit import msgs


def basic(spec, lines, wv_cen, disp, min_ampl=300.,
          swv_uncertainty=350., pix_tol=2, plot_fil=None, min_nmatch=5,
          **kwargs):
    """ Basic algorithm to wavelength calibrate spectroscopic data

    Parameters
    ----------
    spec : ndarray
      Extracted 1D Arc Spectrum
    lines : list
      List of arc lamps on
    wv_cen : float
      Guess at central wavelength
    disp : float
      Dispersion A/pix
    min_ampl : float
      Minimum amplitude of the arc lines that will be used in the fit
    swv_uncertainty : float

    pix_tol : float
      Tolerance in units of pixels to match to
    plot_fil : str, optional
      Name of output file
    min_nmatch : int
      Minimum number of acceptable matches before a solution is deemed to be found
    Returns
    -------
    status : int
      If successful, status=1

    """

    # Init line-lists and wavelength 'guess'
    npix = spec.size
    wave = wv_cen + (np.arange(npix) - npix/2.)*disp

    line_lists = waveio.load_line_lists(lines, unknown=True)
    wvdata = line_lists['wave'].data  # NIST + Extra
    isrt = np.argsort(wvdata)
    wvdata = wvdata[isrt]

    # Find peaks
    all_tcent, cut_tcent, icut = utils.arc_lines_from_spec(spec, min_ampl=min_ampl)

    # Matching
    match_idx, scores = patterns.run_quad_match(cut_tcent, wave, wvdata,
                                                disp, swv_uncertainty=swv_uncertainty,
                                                pix_tol=pix_tol)

    # Check quadrants
    xquad = npix//4 + 1
    msgs.info("================================================================" + msgs.newline() +
              "Checking quadrants:" + msgs.newline() +
              "----------------------------------------------------------------")
    for jj in range(4):
        tc_in_q = (cut_tcent >= jj*xquad) & (cut_tcent < (jj+1)*xquad)
        cstat = '  quad {:d}: ndet={:d}'.format(jj, np.sum(tc_in_q))
        # Stats
        for key in ['Perf', 'Good', 'OK', 'Amb']:
            in_stat = scores[tc_in_q] == key
            cstat += ' {:s}={:d}'.format(key, np.sum(in_stat))
        # Print
            msgs.indent(cstat)
    msgs.indent("----------------------------------------------------------------")

    # Go for it!?
    mask = np.array([False]*len(all_tcent))
    IDs = []
    for kk,score in enumerate(scores):
        if score in ['Perf', 'Good', 'Ok']:
            mask[icut[kk]] = True
            uni, counts = np.unique(match_idx[kk]['matches'], return_counts=True)
            imx = np.argmax(counts)
            IDs.append(wvdata[uni[imx]])
    ngd_match = np.sum(mask)
    if ngd_match < min_nmatch:
        msgs.warn("Insufficient matches to continue")
        status = -1
        return status, ngd_match, match_idx, scores, None

    # Fit
    NIST_lines = line_lists['NIST'] > 0
    ifit = np.where(mask)[0]
    final_fit = fitting.iterative_fitting(spec, all_tcent, ifit,
                                          IDs, line_lists[NIST_lines], disp, plot_fil=plot_fil)
    # Return
    status = 1
    return status, ngd_match, match_idx, scores, final_fit


def semi_brute(spec, lines, wv_cen, disp, min_ampl=300.,
               outroot=None, debug=False, do_fit=True, verbose=False,
               fit_parm=None, min_nmatch=3, lowest_ampl=200.):
    """
    Parameters
    ----------
    spec
    lines
    wv_cen
    disp
    siglev
    min_ampl
    outroot
    debug
    do_fit
    verbose
    fit_parm
    min_nmatch
    lowest_ampl

    Returns
    -------
    best_dict : dict
    final_fit : dict

    """
    # imports
    from astropy.table import vstack
    from linetools import utils as ltu

    # Load line lists
    line_lists = waveio.load_line_lists(lines)
    unknwns = waveio.load_unknown_list(lines)

    npix = spec.size

    # Lines
    all_tcent, cut_tcent, icut = utils.arc_lines_from_spec(spec, min_ampl=min_ampl)

    # Best
    best_dict = dict(nmatch=0, ibest=-1, bwv=0., min_ampl=min_ampl, unknown=False,
                     pix_tol=1, ampl=min_ampl)

    # 3 things to fiddle:
    #  pix_tol -- higher for fewer lines  1/2
    #  unknowns -- on for fewer lines  off/on
    #  scoring -- weaken for more lines ??

    # Loop on unknowns
    #for unknown in [False, True]:
    for unknown in [True]:
        if unknown:
            tot_list = vstack([line_lists,unknwns])
        else:
            tot_list = line_lists
        wvdata = np.array(tot_list['wave'].data) # Removes mask if any
        wvdata.sort()
        sav_nmatch = best_dict['nmatch']

        # Loop on pix_tol
        for pix_tol in [1., 2.]:
            # Scan on wavelengths
            patterns.scan_for_matches(wv_cen, disp, npix, cut_tcent, wvdata,
                                      best_dict=best_dict, pix_tol=pix_tol)
            # Lower minimum amplitude
            ampl = min_ampl
            #pdb.set_trace()
            while(best_dict['nmatch'] < min_nmatch):
                ampl /= 2.
                if ampl < lowest_ampl:
                    break
                all_tcent, cut_tcent, icut = utils.arc_lines_from_spec(spec, min_ampl=ampl)
                patterns.scan_for_matches(wv_cen, disp, npix, cut_tcent, wvdata,
                                          best_dict=best_dict, pix_tol=pix_tol, ampl=ampl)

        #if debug:
        #    pdb.set_trace()
        # Save linelist?
        if best_dict['nmatch'] > sav_nmatch:
            best_dict['line_list'] = tot_list.copy()
            best_dict['unknown'] = unknown
            best_dict['ampl'] = ampl
            best_dict['pix_tol'] = pix_tol

    # Try to pick up some extras by turning off/on unknowns
    if best_dict['unknown']:
        tot_list = line_lists
    else:
        tot_list = vstack([line_lists,unknwns])
    wvdata = np.array(tot_list['wave'].data) # Removes mask if any
    wvdata.sort()
    tmp_dict = best_dict.copy()
    tmp_dict['nmatch'] = 0
    patterns.scan_for_matches(best_dict['bwv'], disp, npix, cut_tcent, wvdata,
                              best_dict=tmp_dict, pix_tol=best_dict['pix_tol'],
                              ampl=best_dict['ampl'], wvoff=1.)
    for kk,ID in enumerate(tmp_dict['IDs']):
        if (ID > 0.) and (best_dict['IDs'][kk] == 0.):
            best_dict['IDs'][kk] = ID
            best_dict['scores'][kk] = tmp_dict['scores'][kk]
            best_dict['mask'][kk] = True
            best_dict['midx'][kk] = tmp_dict['midx'][kk]
            best_dict['nmatch'] += 1
    #pdb.set_trace()

    if best_dict['nmatch'] == 0:
        msgs.info('---------------------------------------------------' + msgs.newline() +
                  'Report:' + msgs.newline() +
                  '  No matches!  Could be you input a bad wvcen or disp value' + msgs.newline() +
                  '---------------------------------------------------')
        return

    # Report
    msgs.info('---------------------------------------------------' + msgs.newline() +
              'Report:' + msgs.newline() +
              '  Number of lines recovered    = {:d}'.format(all_tcent.size) + msgs.newline() +
              '  Number of lines analyzed     = {:d}'.format(cut_tcent.size) + msgs.newline() +
              '  Number of acceptable matches = {:d}'.format(best_dict['nmatch']) + msgs.newline() +
              '  Best central wavelength      = {:g}A'.format(best_dict['bwv']) + msgs.newline() +
              '  Best solution used pix_tol   = {}'.format(best_dict['pix_tol']) + msgs.newline() +
              '  Best solution had unknown    = {}'.format(best_dict['unknown']) + msgs.newline())

    if debug:
        match_idx = best_dict['midx']
        for kk in match_idx.keys():
            uni, counts = np.unique(match_idx[kk]['matches'], return_counts=True)
            msgs.info('kk={}, {}, {}, {}'.format(kk, uni, counts, np.sum(counts)))

    # Write scores
    #out_dict = best_dict['scores']
    #jdict = ltu.jsonify(out_dict)
    #ltu.savejson(pargs.outroot+'.scores', jdict, easy_to_read=True, overwrite=True)

    # Write IDs
    if outroot is not None:
        out_dict = dict(pix=cut_tcent, IDs=best_dict['IDs'])
        jdict = ltu.jsonify(out_dict)
        ltu.savejson(outroot+'.json', jdict, easy_to_read=True, overwrite=True)
        msgs.info("Wrote: {:s}".format(outroot+'.json'))

    # Plot
    if outroot is not None:
        tmp_list = vstack([line_lists,unknwns])
        qa.match_qa(spec, cut_tcent, tmp_list, best_dict['IDs'], best_dict['scores'], outroot+'.pdf')
        msgs.info("Wrote: {:s}".format(outroot+'.pdf'))

    # Fit
    final_fit = None
    if do_fit:
        '''
        # Read in Full NIST Tables
        full_NIST = waveio.load_line_lists(lines, NIST=True)
        # KLUDGE!!!!!
        keep = full_NIST['wave'] > 8800.
        pdb.set_trace()
        line_lists = vstack([line_lists, full_NIST[keep]])
        '''
        #
        NIST_lines = line_lists['NIST'] > 0
        ifit = np.where(best_dict['mask'])[0]
        if outroot is not None:
            plot_fil = outroot+'_fit.pdf'
        else:
            plot_fil = None
        # Purge UNKNOWNS from ifit
        imsk = np.array([True]*len(ifit))
        for kk, idwv in enumerate(np.array(best_dict['IDs'])[ifit]):
            if np.min(np.abs(line_lists['wave'][NIST_lines]-idwv)) > 0.01:
                imsk[kk] = False
        ifit = ifit[imsk]
        # Allow for weaker lines in the fit
        all_tcent, weak_cut_tcent, icut = utils.arc_lines_from_spec(spec, min_ampl=lowest_ampl)
        add_weak = []
        for weak in weak_cut_tcent:
            if np.min(np.abs(cut_tcent-weak)) > 5.:
                add_weak += [weak]
        if len(add_weak) > 0:
            cut_tcent = np.concatenate([cut_tcent, np.array(add_weak)])
        # Fit
        final_fit = fitting.iterative_fitting(spec, cut_tcent, ifit,
                                              np.array(best_dict['IDs'])[ifit], line_lists[NIST_lines],
                                              disp, plot_fil=plot_fil, verbose=verbose, aparm=fit_parm)
        if plot_fil is not None:
            print("Wrote: {:s}".format(plot_fil))

    # Return
    return best_dict, final_fit


def general(spec, lines, ok_mask=None, min_ampl=1000., islinelist=False,
            outroot=None, debug=False, do_fit=True, verbose=False,
            fit_parm=None, lowest_ampl=200.,
            binw=None, bind=None, nstore=10, use_unknowns=True):
    """ General algorithm to wavelength calibrate spectroscopic data

    Parameters
    ----------
    spec : ndarray
      Extracted 1D Arc Spectrum
    lines : list
      List of arc lamps on
    ok_mask : ndarray

    min_ampl : float
      Minimum amplitude of the arc lines that will be used in the fit
    islinelist : bool
      Is lines a linelist (True), or a list of ions (False)
    outroot : str, optional
      Name of output file
    debug : bool
      Used to debug the algorithm
    do_fit : bool
      If True, a fit and iterative identification of arc lines will be performed.
      If False, the final fit will not be computed, and only the initial IDs will
      be returned (as well as a blank list of empty dicts for the final fit).
    verbose : bool
      If True, the final fit will print out more detail as the RMS is refined,
      and lines are rejected. This is mostly helpful for developing the algorithm.
    fit_parm : dict
      Fitting parameter dictionary (see fitting.iterative_fitting)
    lowest_ampl : float
    binw : ndarray, optional
      Set the wavelength grid when identifying the best solution
    bind : ndarray, optional
      Set the dispersion grid when identifying the best solution
    nstore : int
      The number of "best" initial solutions to consider
    use_unknowns : bool
      If True, arc lines that are known to be present in the spectra, but
      have not been attributed to an element+ion, will be included in the fit.

    Returns
    -------
    all_patt_dict : list of dicts
      A list of dictionaries, which contain the results from the preliminary
      pattern matching algorithm providing the first guess at the ID lines
    all_final_fit : list of dicts
      A list of dictionaries, which contain the full fitting results and
      final best guess of the line IDs
    """
    from astropy.table import vstack
    from linetools import utils as ltu

    # Import the triangles algorithm
    from pypeit.core.wavecal.patterns import triangles

    npix, nslit = spec.shape

    if ok_mask is None:
        ok_mask = np.arange(nslit)

    # Load line lists
    if islinelist:
        line_lists = lines
        unknwns = lines[:0].copy()
    else:
        if 'ThAr' in lines:
            line_lists_all = waveio.load_line_lists(lines)
            line_lists = line_lists_all[np.where(line_lists_all['ion'] != 'UNKNWN')]
            unknwns = line_lists_all[np.where(line_lists_all['ion'] == 'UNKNWN')]
        else:
            line_lists = waveio.load_line_lists(lines)
            unknwns = waveio.load_unknown_list(lines)

    if use_unknowns:
        tot_list = vstack([line_lists, unknwns])
    else:
        tot_list = line_lists
    wvdata = np.array(tot_list['wave'].data)  # Removes mask if any
    wvdata.sort()

    # Setup grid parameters
    #ngridw, ngridd = 1000, 1000  # Longslit
    #ngridw, ngridd = 100000, 100  # Echelle

    nselw, nseld = 5, 25  # Longslit
    # nselw, nseld = 1, 1  # Echelle

    # The wavelength grid (i.e. the binw size) should depend on the dispersion.

    # Set the wavelength grid
    if binw is None:
        # Ideally, you want binw to roughly sample the A/pix of the spectrograph
        ngridw = 1000
        binw = np.linspace(np.min(wvdata), np.max(wvdata), ngridw)
    else:
        ngridw = binw.size
    # Set the dispersion grid
    if bind is None:
        ngridd = 1000
        bind = np.linspace(-3.0, 1.0, ngridd)
    else:
        ngridd = bind.size

    bestlist = []
    allwcen, alldisp, allhnum = np.array([]), np.array([]), np.array([])
    slit_tcent = []
    for cnt, slit in enumerate(ok_mask):
        bestlist.append([])
        # Lines
        all_tcent, cut_tcent, icut = utils.arc_lines_from_spec(spec[:, slit], min_ampl=min_ampl)

        # Decide which tcent to use (either all_tcent or cut_tcent)
        use_tcent = all_tcent.copy()
        slit_tcent.append(use_tcent.copy())

        if use_tcent.size == 0:
            msgs.warn("No lines to identify in slit {0:d}!".format(slit))
            bestlist[cnt].append([None]*6)
            continue

        # Loop on pix_tol
        # TODO: Allow for different pixel tolerance?
        for pix_tol in [1.0]:
            # First run pattern recognition assuming pixels correlate with wavelength

            # Triangle pattern matching
            dindexp, lindexp, wvcenp, dispsp = triangles(use_tcent, wvdata, npix, 5, 10, pix_tol)
            # dindexp, lindexp, wvcenp, dispsp = triangles(use_tcent, wvdata, npix, 3, 6, pix_tol)
            # Remove any invalid results
            ww = np.where((binw[0] < wvcenp) & (wvcenp < binw[-1]) & (10.0**bind[0] < dispsp) & (dispsp < 10.0**bind[-1]))
            dindexp = dindexp[ww[0], :]
            lindexp = lindexp[ww[0], :]
            dispsp = dispsp[ww]
            wvcenp = wvcenp[ww]

            # Now run pattern recognition assuming pixels correlate with wavelength
            use_tcent = (npix - 1.0) - all_tcent.copy()[::-1]
            # Triangle pattern matching
            dindexm, lindexm, wvcenm, dispsm = triangles(use_tcent, wvdata, npix, 5, 10, pix_tol)
            #dindexm, lindexm, wvcenm, dispsm = triangles(use_tcent, wvdata, npix, 3, 6, pix_tol)
            # Remove any invalid results
            ww = np.where((binw[0] < wvcenm) & (wvcenm < binw[-1]) & (10.0**bind[0] < dispsm) & (dispsm < 10.0**bind[-1]))
            dindexm = dindexm[ww[0], :]
            lindexm = lindexm[ww[0], :]
            dispsm = dispsm[ww]
            wvcenm = wvcenm[ww]
            # Construct the histograms
            histimgp, xed, yed = np.histogram2d(wvcenp, np.log10(dispsp), bins=[binw, bind])
            histimgm, xed, yed = np.histogram2d(wvcenm, np.log10(dispsm), bins=[binw, bind])
            #histimgp = gaussian_filter(histimgp, 3)
            #histimgm = gaussian_filter(histimgm, 3)
            histimg = histimgp - histimgm
            #histimg = gaussian_filter(histimg, 6)

            histpeaks = patterns.detect_peaks(np.abs(histimg))

            # Find the indices of the nstore largest peaks
            bidx = np.unravel_index(np.argpartition(np.abs(histpeaks*histimg), -nstore, axis=None)[-nstore:], histimg.shape)

            debug = False
            if debug:
                from matplotlib import pyplot as plt
                plt.clf()
                plt.imshow(np.log10(np.abs(histimg[:, ::-1].T)), extent=[binw[0], binw[-1], bind[0], bind[-1]], aspect='auto')
                #plt.imshow(histimg[:, ::-1].T, extent=[binw[0], binw[-1], bind[0], bind[-1]], aspect='auto')
                plt.plot(binw[bidx[0]], bind[bidx[1]], 'ro')
                #plt.axvline(binw[bidx[0]], color='r', linestyle='--')
                #plt.axhline(bind[bidx[1]], color='r', linestyle='--')
                plt.show()
                pdb.set_trace()
                plt.imshow(histimgp[:, ::-1].T, extent=[binw[0], binw[-1], bind[0], bind[-1]], aspect='auto')

            # Get the peak value of central wavelength and dispersion
            wcenval = binw[bidx[0]]
            dispval = bind[bidx[1]]
            histnum = np.abs(histimg[bidx])

            # Find all good solutions
            for idx in range(nstore):
                # Select all solutions around the best solution within a square of side 2*nsel
                wlo = binw[max(0, bidx[0][idx] - nselw)]
                whi = binw[min(ngridw - 1, bidx[0][idx] + nselw)]
                dlo = 10.0 ** bind[max(0, bidx[1][idx] - nseld)]
                dhi = 10.0 ** bind[min(ngridd - 1, bidx[1][idx] + nseld)]
                if histimgp[bidx][idx] > histimgm[bidx][idx]:
                    wgd = np.where((wvcenp > wlo) & (wvcenp < whi) & (dispsp > dlo) & (dispsp < dhi))
                    dindex = dindexp[wgd[0], :].flatten()
                    lindex = lindexp[wgd[0], :].flatten()
                    sign = +1.0
                else:
                    wgd = np.where((wvcenm > wlo) & (wvcenm < whi) & (dispsm > dlo) & (dispsm < dhi))
                    dindex = dindexm[wgd[0], :].flatten()
                    lindex = lindexm[wgd[0], :].flatten()
                    sign = -1.0
                # Store relevant values in an array to solve for best solution
                bestlist[cnt].append([wcenval[idx], dispval[idx], histnum[idx], sign, dindex, lindex])
            allwcen = np.append(allwcen, wcenval)
            alldisp = np.append(alldisp, dispval)
            allhnum = np.append(allhnum, histnum)

    # Using the results from all slits, decide on the best solutions (assume all slits have the same dispersion)
    dhist, dedge = np.histogram(alldisp, bins=bind, weights=allhnum)
    dhmax = np.argmax(dhist)
    if debug:
        from matplotlib import pyplot as plt
        null = plt.hist(alldisp, bins=bind, weights=allhnum)
        plt.show()
    msgs.info("Best initial guess for spectrograph dispersion: {0:.4f}A/pixel".format(10.0**np.mean(dedge[dhmax:dhmax+2])))
    msgs.info("Fitting the wavelength solution for each slit")

    # Fit the wavelength solution for each slit
    all_patt_dict, all_final_fit = {}, {}
    for cnt, slit in enumerate(ok_mask):
        # patt_dict
        patt_dict = dict(nmatch=0, ibest=-1, bwv=0., min_ampl=min_ampl)

        # Check there are lines in this slit
        if slit_tcent[cnt].size == 0:
            msgs.warn("No lines to identify in slit {0:d}!".format(slit))
            all_patt_dict[str(slit)] = None
            all_final_fit[str(slit)] = None
            continue

        # Obtain a full list of indices that are consistent with the maximum value
        dindex, lindex, allsgn = np.array([]), np.array([]), np.array([])
        dcen, wcen = np.array([]), np.array([])
        for ss in range(len(bestlist[cnt])):
            if dedge[dhmax-nseld] <= bestlist[cnt][ss][1] <= dedge[dhmax+1+nseld]:
                wcen = np.append(wcen, bestlist[cnt][ss][0])
                dcen = np.append(dcen, bestlist[cnt][ss][1])
                allsgn = np.append(allsgn, bestlist[cnt][ss][3]*np.ones(bestlist[cnt][ss][4].size))
                dindex = np.append(dindex, bestlist[cnt][ss][4])
                lindex = np.append(lindex, bestlist[cnt][ss][5])
        # Find the favoured sign and only use those values
        if np.sum(allsgn) > 0.0:
            use_tcent = slit_tcent[cnt].copy()
            sign = +1.0
            signtxt = "correlate"
        else:
            use_tcent = (npix - 1.0) - slit_tcent[cnt].copy()[::-1]
            sign = -1.0
            signtxt = "anticorrelate"
        dindex = dindex[np.where(allsgn == sign)]
        lindex = lindex[np.where(allsgn == sign)]
        patterns.solve_triangles(use_tcent, wvdata, dindex, lindex, patt_dict)

        # Fill in the patterns dictionary
        patt_dict['bwv'] = np.mean(wcen)
        patt_dict['bdisp'] = 10.0**np.mean(dcen)

        # Check that a solution has been found
        if patt_dict['nmatch'] == 0:
            msgs.info('---------------------------------------------------' + msgs.newline() +
                      'Initial report for slit {0:d}/{1:d}:'.format(slit+1, nslit) + msgs.newline() +
                      '  No matches! Try another algorithm' + msgs.newline() +
                      '---------------------------------------------------')
            all_patt_dict[str(slit)] = None
            all_final_fit[str(slit)] = None
            continue

        # Report
        msgs.info('---------------------------------------------------' + msgs.newline() +
                  'Initial report for slit {0:d}/{1:d}:'.format(slit+1, nslit) + msgs.newline() +
                  '  Pixels {:s} with wavelength'.format(signtxt) + msgs.newline() +
                  '  Number of lines recovered    = {:d}'.format(all_tcent.size) + msgs.newline() +
                  '  Number of lines analyzed     = {:d}'.format(use_tcent.size) + msgs.newline() +
                  '  Number of acceptable matches = {:d}'.format(patt_dict['nmatch']) + msgs.newline() +
                  '  Best central wavelength      = {:g}A'.format(patt_dict['bwv']) + msgs.newline() +
                  '  Best dispersion              = {:g}A/pix'.format(patt_dict['bdisp']) + msgs.newline() +
                  '  Best solution had unknown    = {}'.format(use_unknowns) + msgs.newline() +
                  '---------------------------------------------------')

        slittxt = '_Slit{0:03d}'.format(slit)
        if outroot is not None:
            # Write IDs
            out_dict = dict(pix=use_tcent, IDs=patt_dict['IDs'])
            jdict = ltu.jsonify(out_dict)
            ltu.savejson(outroot + slittxt + '.json', jdict, easy_to_read=True, overwrite=True)
            msgs.info("Wrote: {:s}".format(outroot + slittxt + '.json'))

            # Plot
            tmp_list = vstack([line_lists, unknwns])
            qa.match_qa(spec, use_tcent, tmp_list,
                        patt_dict['IDs'], patt_dict['scores'], outroot + slittxt + '.pdf')
            msgs.info("Wrote: {:s}".format(outroot + slittxt + '.pdf'))

        # Perform final fit to the line IDs
        final_fit = dict()
        if do_fit:
            NIST_lines = line_lists['NIST'] > 0
            ifit = np.where(patt_dict['mask'])[0]
            if outroot is not None:
                plot_fil = outroot + slittxt + '_fit.pdf'
            else:
                plot_fil = None
            # Purge UNKNOWNS from ifit
            imsk = np.ones(len(ifit), dtype=np.bool)
            for kk, idwv in enumerate(np.array(patt_dict['IDs'])[ifit]):
                if np.min(np.abs(line_lists['wave'][NIST_lines]-idwv)) > 0.01:
                    imsk[kk] = False
            ifit = ifit[imsk]
            # Allow for weaker lines in the fit
            all_tcent, weak_cut_tcent, icut = utils.arc_lines_from_spec(spec[:, slit], min_ampl=lowest_ampl)
            use_weak_tcent = all_tcent.copy()
            add_weak = []
            for weak in use_weak_tcent:
                if np.min(np.abs(all_tcent-weak)) > 5.:
                    add_weak += [weak]
            if len(add_weak) > 0:
                if sign == +1.0:
                    use_weak = np.array(add_weak)
                else:
                    use_weak = (npix - 1.0) - np.array(add_weak)[::-1]
                use_tcent = np.concatenate([use_tcent, use_weak])
            # Fit
            final_fit = fitting.iterative_fitting(spec, use_tcent, ifit,
                                                  np.array(patt_dict['IDs'])[ifit], line_lists[NIST_lines],
                                                  patt_dict['bdisp'], plot_fil=plot_fil, verbose=verbose,
                                                  aparm=fit_parm)

            if plot_fil is not None:
                print("Wrote: {:s}".format(plot_fil))

        # Append the results to the full list
        all_patt_dict[str(slit)] = patt_dict.copy()
        all_final_fit[str(slit)] = final_fit.copy()

    # Return
    return all_patt_dict, all_final_fit
