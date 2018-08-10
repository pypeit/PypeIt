""" Module for finding patterns in arc line spectra
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

from scipy.ndimage.filters import gaussian_filter
import numpy as np
import pdb

from arclines import io as arcl_io
from arclines.holy import patterns as arch_patt
from arclines.holy import fitting as arch_fit
from arclines.holy import utils as arch_utils


def basic(spec, lines, wv_cen, disp, siglev=20., min_ampl=300.,
          swv_uncertainty=350., pix_tol=2, plot_fil=None, min_match=5,
          **kwargs):
    """ Basic holy grail algorithm

    Parameters
    ----------
    spec : spectrum
    lines : list
      List of arc lamps on
    wv_cen : float
      Guess at central wavelength
    disp : float
      Dispersion A/pix
    siglev
    min_ampl
    swv_uncertainty
    pix_tol
    plot_fil

    Returns
    -------
    status : int

    """

    # Init line-lists and wavelength 'guess'
    npix = spec.size
    wave = wv_cen + (np.arange(npix) - npix/2.)*disp

    line_lists = arcl_io.load_line_lists(lines, unknown=True)
    wvdata = line_lists['wave'].data  # NIST + Extra
    isrt = np.argsort(wvdata)
    wvdata = wvdata[isrt]

    # Find peaks
    all_tcent, cut_tcent, icut = arch_utils.arc_lines_from_spec(spec, min_ampl=min_ampl)

    # Matching
    match_idx, scores = arch_patt.run_quad_match(cut_tcent, wave, wvdata,
                                                 disp, swv_uncertainty=swv_uncertainty,
                                                 pix_tol=pix_tol)

    # Check quadrants
    xquad = npix//4 + 1
    print("================================================================")
    print("Checking quadrants:")
    print("----------------------------------------------------------------")
    for jj in range(4):
        tc_in_q = (cut_tcent >= jj*xquad) & (cut_tcent < (jj+1)*xquad)
        cstat = 'quad {:d}: ndet={:d}'.format(jj, np.sum(tc_in_q))
        # Stats
        for key in ['Perf', 'Good', 'OK', 'Amb']:
            in_stat = scores[tc_in_q] == key
            cstat += ' {:s}={:d}'.format(key, np.sum(in_stat))
        # Print
        print(cstat)
    print("----------------------------------------------------------------")

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
    if ngd_match < min_match:
        print("Insufficient matches to continue")
        status = -1
        return status, ngd_match, match_idx, scores, None

    # Fit
    NIST_lines = line_lists['NIST'] > 0
    ifit = np.where(mask)[0]
    final_fit = arch_fit.iterative_fitting(spec, all_tcent, ifit,
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
    from arclines import plots as arcl_plots
    # Load line lists
    line_lists = arcl_io.load_line_lists(lines)
    unknwns = arcl_io.load_unknown_list(lines)

    npix = spec.size

    # Lines
    all_tcent, cut_tcent, icut = arch_utils.arc_lines_from_spec(spec, min_ampl=min_ampl)

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
        for pix_tol in [1.,2.]:
            # Scan on wavelengths
            arch_patt.scan_for_matches(wv_cen, disp, npix, cut_tcent, wvdata,
                                      best_dict=best_dict, pix_tol=pix_tol)
            # Lower minimum amplitude
            ampl = min_ampl
            #pdb.set_trace()
            while(best_dict['nmatch'] < min_nmatch):
                ampl /= 2.
                if ampl < lowest_ampl:
                    break
                all_tcent, cut_tcent, icut = arch_utils.arc_lines_from_spec(spec, min_ampl=ampl)
                arch_patt.scan_for_matches(wv_cen, disp, npix, cut_tcent, wvdata,
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
    arch_patt.scan_for_matches(best_dict['bwv'], disp, npix, cut_tcent, wvdata,
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
        print('---------------------------------------------------')
        print('Report:')
        print('::   No matches!  Could be you input a bad wvcen or disp value')
        print('---------------------------------------------------')
        return

    # Report
    print('---------------------------------------------------')
    print('Report:')
    print('::   Number of lines recovered = {:d}'.format(all_tcent.size))
    print('::   Number of lines analyzed = {:d}'.format(cut_tcent.size))
    print('::   Number of Perf/Good/Ok matches = {:d}'.format(best_dict['nmatch']))
    print('::   Best central wavelength = {:g}A'.format(best_dict['bwv']))
    print('::   Best solution used pix_tol = {}'.format(best_dict['pix_tol']))
    print('::   Best solution had unknown = {}'.format(best_dict['unknown']))
    print('---------------------------------------------------')

    if debug:
        match_idx = best_dict['midx']
        for kk in match_idx.keys():
            uni, counts = np.unique(match_idx[kk]['matches'], return_counts=True)
            print('kk={}, {}, {}, {}'.format(kk, uni, counts, np.sum(counts)))

    # Write scores
    #out_dict = best_dict['scores']
    #jdict = ltu.jsonify(out_dict)
    #ltu.savejson(pargs.outroot+'.scores', jdict, easy_to_read=True, overwrite=True)

    # Write IDs
    if outroot is not None:
        out_dict = dict(pix=cut_tcent, IDs=best_dict['IDs'])
        jdict = ltu.jsonify(out_dict)
        ltu.savejson(outroot+'.json', jdict, easy_to_read=True, overwrite=True)
        print("Wrote: {:s}".format(outroot+'.json'))

    # Plot
    if outroot is not None:
        tmp_list = vstack([line_lists,unknwns])
        arcl_plots.match_qa(spec, cut_tcent, tmp_list,
                            best_dict['IDs'], best_dict['scores'], outroot+'.pdf')
        print("Wrote: {:s}".format(outroot+'.pdf'))

    # Fit
    final_fit = None
    if do_fit:
        '''
        # Read in Full NIST Tables
        full_NIST = arcl_io.load_line_lists(lines, NIST=True)
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
        all_tcent, weak_cut_tcent, icut = arch_utils.arc_lines_from_spec(spec, min_ampl=lowest_ampl)
        add_weak = []
        for weak in weak_cut_tcent:
            if np.min(np.abs(cut_tcent-weak)) > 5.:
                add_weak += [weak]
        if len(add_weak) > 0:
            cut_tcent = np.concatenate([cut_tcent, np.array(add_weak)])
        # Fit
        final_fit = arch_fit.iterative_fitting(spec, cut_tcent, ifit,
                                               np.array(best_dict['IDs'])[ifit], line_lists[NIST_lines],
                                               disp, plot_fil=plot_fil, verbose=verbose, aparm=fit_parm)
        if plot_fil is not None:
            print("Wrote: {:s}".format(plot_fil))

    # Return
    return best_dict, final_fit


def general(spec, lines, min_ampl=300.,
            outroot=None, debug=False, do_fit=True, verbose=False,
            fit_parm=None, lowest_ampl=200.):
    """
    Parameters
    ----------
    spec
    lines
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
    from arclines import plots as arcl_plots

    # Import the triangles algorithm
    from arclines.holy.patterns import triangles

    # Load line lists
    line_lists = arcl_io.load_line_lists(lines)
    unknwns = arcl_io.load_unknown_list(lines)

    npix = spec.size

    # Lines
    all_tcent, cut_tcent, icut = arch_utils.arc_lines_from_spec(spec, min_ampl=min_ampl)
    use_tcent = all_tcent.copy()
    #use_tcent = cut_tcent.copy()  # min_ampl is having not effect at present

    # Best
    best_dict = dict(nmatch=0, ibest=-1, bwv=0., min_ampl=min_ampl)

    ngrid = 1000

    # Loop on unknowns
    for unknown in [False, True]:
        if unknown:
            tot_list = vstack([line_lists,unknwns])
        else:
            tot_list = line_lists
        wvdata = np.array(tot_list['wave'].data)  # Removes mask if any
        wvdata.sort()

        sav_nmatch = best_dict['nmatch']

        # Loop on pix_tol
        for pix_tol in [1.]:#, 2.]:
            # Triangle pattern matching
            dindex, lindex, wvcen, disps = triangles(use_tcent, wvdata, npix, 5, 10, pix_tol)

            # Remove any invalid results
            ww = np.where((wvcen > 0.0) & (disps > 0.0))
            dindex = dindex[ww[0], :]
            lindex = lindex[ww[0], :]
            disps = disps[ww]
            wvcen = wvcen[ww]

            # Setup the grids and histogram
            binw = np.linspace(max(np.min(wvcen), np.min(wvdata)), min(np.max(wvcen), np.max(wvdata)), ngrid)
            bind = np.linspace(np.min(np.log10(disps)), np.max(np.log10(disps)), ngrid)
            histimg, xed, yed = np.histogram2d(wvcen, np.log10(disps), bins=[binw, bind])
            histimg = gaussian_filter(histimg, 3)

            # Find the best combination of central wavelength and dispersion
            bidx = np.unravel_index(np.argmax(histimg), histimg.shape)

            debug = False
            if debug:
                from matplotlib import pyplot as plt
                plt.clf()
                plt.imshow(histimg[:, ::-1].T, extent=[binw[0], binw[-1], bind[0], bind[-1]], aspect='auto')
                plt.axvline(binw[bidx[0]], color='r', linestyle='--')
                plt.axhline(bind[bidx[1]], color='r', linestyle='--')
                plt.show()
                print(histimg[bidx], binw[bidx[0]], 10.0**bind[bidx[1]])
                pdb.set_trace()

            # Find all good solutions
            nsel = 5  # Select all solutions around the best solution within a square of side 2*nsel
            wlo = binw[bidx[0] - nsel]
            whi = binw[bidx[0] + nsel]
            dlo = 10.0 ** bind[bidx[1] - 5*nsel]
            dhi = 10.0 ** bind[bidx[1] + 5*nsel]
            wgd = np.where((wvcen > wlo) & (wvcen < whi) & (disps > dlo) & (disps < dhi))
            dindex = dindex[wgd[0], :].flatten()
            lindex = lindex[wgd[0], :].flatten()

            # Given this solution, fit for all detlines
            arch_patt.solve_triangles(use_tcent, wvdata, dindex, lindex, best_dict)
            if best_dict['nmatch'] > sav_nmatch:
                best_dict['pix_tol'] = pix_tol

        # Save linelist?
        if best_dict['nmatch'] > sav_nmatch:
            best_dict['bwv'] = binw[bidx[0]]
            best_dict['bdisp'] = 10.0**bind[bidx[1]]
            best_dict['line_list'] = tot_list.copy()
            best_dict['unknown'] = unknown
            best_dict['ampl'] = unknown

    # Try to pick up some extras by turning off/on unknowns
    if best_dict['unknown']:
        tot_list = line_lists
    else:
        tot_list = vstack([line_lists,unknwns])

    # Retrieve the wavelengths of the linelist and sort
    wvdata = np.array(tot_list['wave'].data)  # Removes mask if any
    wvdata.sort()

    if best_dict['nmatch'] == 0:
        print('---------------------------------------------------')
        print('Report:')
        print('::   No matches! Try another algorithm')
        print('---------------------------------------------------')
        return

    # Report
    print('---------------------------------------------------')
    print('Report:')
    print('::   Number of lines recovered    = {:d}'.format(all_tcent.size))
    print('::   Number of lines analyzed     = {:d}'.format(use_tcent.size))
    print('::   Number of acceptable matches = {:d}'.format(best_dict['nmatch']))
    print('::   Best central wavelength      = {:g}A'.format(best_dict['bwv']))
    print('::   Best dispersion              = {:g}A/pix'.format(best_dict['bdisp']))
    print('::   Best solution used pix_tol   = {}'.format(best_dict['pix_tol']))
    print('::   Best solution had unknown    = {}'.format(best_dict['unknown']))
    print('---------------------------------------------------')

    # Write IDs
    if outroot is not None:
        out_dict = dict(pix=use_tcent, IDs=best_dict['IDs'])
        jdict = ltu.jsonify(out_dict)
        ltu.savejson(outroot+'.json', jdict, easy_to_read=True, overwrite=True)
        print("Wrote: {:s}".format(outroot+'.json'))

    # Plot
    if outroot is not None:
        tmp_list = vstack([line_lists, unknwns])
        arcl_plots.match_qa(spec, use_tcent, tmp_list,
                            best_dict['IDs'], best_dict['scores'], outroot+'.pdf')
        print("Wrote: {:s}".format(outroot+'.pdf'))

    # Fit
    final_fit = None
    if do_fit:
        # Good lines = NIST or OH
        good_lines = np.any([line_lists['NIST'] > 0, line_lists['ion'] == 'OH'], axis=0)
        #
        ifit = np.where(best_dict['mask'])[0]
        if outroot is not None:
            plot_fil = outroot+'_fit.pdf'
        else:
            plot_fil = None
        # Purge UNKNOWNS from ifit
        imsk = np.array([True]*len(ifit))
        for kk, idwv in enumerate(np.array(best_dict['IDs'])[ifit]):
            if np.min(np.abs(line_lists['wave'][good_lines]-idwv)) > 0.01:
                imsk[kk] = False
        ifit = ifit[imsk]
        # Allow for weaker lines in the fit
        all_tcent, weak_cut_tcent, icut = arch_utils.arc_lines_from_spec(spec, min_ampl=lowest_ampl)
        use_weak_tcent = all_tcent.copy()
        add_weak = []
        for weak in use_weak_tcent:
            if np.min(np.abs(use_tcent-weak)) > 5.:
                add_weak += [weak]
        if len(add_weak) > 0:
            use_tcent = np.concatenate([use_tcent, np.array(add_weak)])
        # Fit
        final_fit = arch_fit.iterative_fitting(spec, use_tcent, ifit,
                                               np.array(best_dict['IDs'])[ifit], line_lists[good_lines],
                                               best_dict['bdisp'], plot_fil=plot_fil, verbose=verbose,
                                               aparm=fit_parm)
        if plot_fil is not None:
            print("Wrote: {:s}".format(plot_fil))

    # Return
    return best_dict, final_fit
