""" Module for finding patterns in arc line spectra
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

from scipy.ndimage.filters import gaussian_filter
from scipy.spatial import cKDTree
from scipy import signal
from scipy import stats
import scipy
from linetools import utils as ltu
from astropy.table import vstack
import copy
import numpy as np
import pdb

from pypeit.core.wavecal import kdtree_generator
from pypeit.core.wavecal import waveio
from pypeit.core.wavecal import patterns
from pypeit.core.wavecal import fitting
from pypeit.core.wavecal import wvutils
from pypeit.core.wavecal import qa
from pypeit import utils

from pypeit import msgs
from pypeit import debugger


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
    all_tcent, cut_tcent, icut = wvutils.arc_lines_from_spec(spec, min_ampl=min_ampl)

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
    all_tcent, cut_tcent, icut = wvutils.arc_lines_from_spec(spec, min_ampl=min_ampl)

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
                all_tcent, cut_tcent, icut = wvutils.arc_lines_from_spec(spec, min_ampl=ampl)
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
    tmp_dict = copy.deepcopy(best_dict)
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
        all_tcent, weak_cut_tcent, icut = wvutils.arc_lines_from_spec(spec, min_ampl=lowest_ampl)
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


class General:
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
    verbose : bool
      If True, the final fit will print out more detail as the RMS is refined,
      and lines are rejected. This is mostly helpful for developing the algorithm.
    fit_parm : dict
      Fitting parameter dictionary (see fitting.iterative_fitting)
    lowest_ampl : float
      Lowest amplitude of an arc line that will be used in teh final fit
    rms_threshold : float
      Maximum RMS dispersion that is considered acceptable for a good solution
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

    def __init__(self, spec, lines, ok_mask=None, min_ampl=1000., islinelist=False,
              outroot=None, debug=False, verbose=False,
              fit_parm=None, lowest_ampl=200., rms_threshold=0.1,
              binw=None, bind=None, nstore=1, use_unknowns=True):

        # Set some default parameters
        self._spec = spec
        self._lines = lines
        self._npix, self._nslit = spec.shape
        self._nstore = nstore
        self._binw = binw
        self._bind = bind

        if ok_mask is None:
            self._ok_mask = np.arange(self._nslit)
        else:
            self._ok_mask = ok_mask

        self._min_ampl = min_ampl
        self._lowest_ampl = lowest_ampl

        self._use_unknowns = use_unknowns
        self._islinelist = islinelist

        self._outroot = outroot
        self._fit_parm = fit_parm
        self._rms_threshold = rms_threshold

        self._debug = debug
        self._verbose = verbose

        # Load the linelist to be used for pattern matching
        self.load_linelist()

        # Set up the grids to be used for pattern matching
        self.set_grids()

        # Find the wavelength solution!
        # KD Tree algorithm only works for ThAr - check first that this is what is being used
        if 'ThAr' in lines and len(lines) == 1:
            msgs.info("Using KD Tree pattern matching algorithm to wavelength calibrate")
            self.run_kdtree()
        else:
            msgs.info("Using brute force pattern matching algorithm to wavelength calibrate")
            self.run_brute()

    def get_results(self):
        return copy.deepcopy(self._all_patt_dict), copy.deepcopy(self._all_final_fit)

    def load_linelist(self):
        # Load line lists
        if self._islinelist:
            self._line_lists = self._lines
            self._unknwns = self._lines[:0].copy()
        else:
            if 'ThAr' in self._lines:
                line_lists_all = waveio.load_line_lists(self._lines)
                self._line_lists = line_lists_all[np.where(line_lists_all['ion'] != 'UNKNWN')]
                self._unknwns = line_lists_all[np.where(line_lists_all['ion'] == 'UNKNWN')]
            else:
                self._line_lists = waveio.load_line_lists(self._lines)
                self._unknwns = waveio.load_unknown_list(self._lines)

        if self._use_unknowns:
            self._tot_list = vstack([self._line_lists, self._unknwns])
        else:
            self._tot_list = self._line_lists

        # Generate the final linelist and sort
        self._wvdata = np.array(self._tot_list['wave'].data)  # Removes mask if any
        self._wvdata.sort()
        return

    def set_grids(self):
        # Set the wavelength grid
        if self._binw is None:
            # Ideally, you want binw to roughly sample the A/pix of the spectrograph
            self._ngridw = 200
            self._binw = np.linspace(np.min(self._wvdata), np.max(self._wvdata), self._ngridw)
        else:
            self._ngridw = self._binw.size
        # Set the dispersion grid
        if self._bind is None:
            self._ngridd = 2000
            self._bind = np.linspace(-3.0, 1.0, self._ngridd)
        else:
            self._ngridd = self._bind.size
        return

    def run_brute(self):
        """Run through the parameter space and determine the best solution
        """

        # Set the parameter space that gets searched
        rng_poly = [3, 4]            # Range of algorithms to check (only trigons+tetragons are supported)
        rng_list = range(3, 10)      # Number of lines to search over for the linelist
        rng_detn = range(3, 10)      # Number of lines to search over for the detected lines
        rng_pixt = [1.0]             # Pixel tolerance

        self._all_patt_dict = {}
        self._all_final_fit = {}
        good_fit = np.zeros(self._nslit, dtype=np.bool)
        self._detections = {}
        for slit in range(self._nslit):
            if slit not in self._ok_mask:
                continue
            # Detect lines, and decide which tcent to use
            self._all_tcent, self._cut_tcent, self._icut =\
                wvutils.arc_lines_from_spec(self._spec[:, slit].copy(), min_ampl=self._min_ampl)
            self._all_tcent_weak, self._cut_tcent_weak, self._icut_weak =\
                wvutils.arc_lines_from_spec(self._spec[:, slit].copy(), min_ampl=self._lowest_ampl)
            if self._all_tcent.size == 0:
                msgs.warn("No lines to identify in slit {0:d}!".format(slit))
                continue
            self._detections[str(slit)] = self._all_tcent_weak.copy()
            best_patt_dict, best_final_fit = None, None
            # Loop through parameter space
            for poly in rng_poly:
                for lstsrch in rng_list:
                    for detsrch in rng_detn:
                        for pix_tol in rng_pixt:
                            psols, msols = self.results_brute(poly=poly, pix_tol=pix_tol,
                                                              detsrch=detsrch, lstsrch=lstsrch)
                            patt_dict, final_fit = self.solve_slit(slit, psols, msols)
                            if final_fit is None:
                                # This is not a good solution
                                continue
                            # Test if this solution is better than the currently favoured solution
                            if best_patt_dict is None:
                                # First time a fit is found
                                best_patt_dict, best_final_fit = copy.deepcopy(patt_dict), copy.deepcopy(final_fit)
                                continue
                            elif final_fit['rms'] < self._rms_threshold:
                                # Has a better fit been identified (i.e. more lines identified)?
                                if len(final_fit['xfit']) > len(best_final_fit['xfit']):
                                    best_patt_dict, best_final_fit = copy.deepcopy(patt_dict), copy.deepcopy(final_fit)
            # Report on the best preliminary result
            if best_final_fit is None:
                msgs.warn('---------------------------------------------------' + msgs.newline() +
                          'Preliminary report for slit {0:d}/{1:d}:'.format(slit+1, self._nslit) + msgs.newline() +
                          '  No matches! Attempting to cross match.' + msgs.newline() +
                          '---------------------------------------------------')
                self._all_patt_dict[str(slit)] = None
                self._all_final_fit[str(slit)] = None
            elif best_final_fit['rms'] > self._rms_threshold:
                msgs.warn('---------------------------------------------------' + msgs.newline() +
                          'Preliminary report for slit {0:d}/{1:d}:'.format(slit + 1, self._nslit) + msgs.newline() +
                          '  Poor RMS ({0:.3f})! Attempting to cross match.'.format(best_final_fit['rms']) + msgs.newline() +
                          '---------------------------------------------------')
                self._all_patt_dict[str(slit)] = None
                self._all_final_fit[str(slit)] = None
            else:
                good_fit[slit] = True
                if best_patt_dict['sign'] == +1:
                    signtxt = 'correlate'
                else:
                    signtxt = 'anitcorrelate'
                # Report
                msgs.info('---------------------------------------------------' + msgs.newline() +
                          'Preliminary report for slit {0:d}/{1:d}:'.format(slit+1, self._nslit) + msgs.newline() +
                          '  Pixels {:s} with wavelength'.format(signtxt) + msgs.newline() +
                          '  Number of lines recovered    = {:d}'.format(self._all_tcent.size) + msgs.newline() +
                          '  Number of lines analyzed     = {:d}'.format(len(best_final_fit['xfit'])) + msgs.newline() +
                          '  Number of pattern matches    = {:d}'.format(best_patt_dict['nmatch']) + msgs.newline() +
                          '  Best central wavelength      = {:g}A'.format(best_patt_dict['bwv']) + msgs.newline() +
                          '  Best dispersion              = {:g}A/pix'.format(best_patt_dict['bdisp']) + msgs.newline() +
                          '  Final RMS of fit             = {:g}'.format(best_final_fit['rms']) + msgs.newline() +
                          '---------------------------------------------------')
                self._all_patt_dict[str(slit)] = copy.deepcopy(best_patt_dict)
                self._all_final_fit[str(slit)] = copy.deepcopy(best_final_fit)

        # Now that all slits have been inspected, cross match to generate a
        # master list of all lines in every slit, and refit all spectra
        if self._nslit > 1:
            obad_slits = self.cross_match(good_fit)
            cntr = 0  # Introduce a counter to stop the while loop, just in case the loop gets stuck
            while obad_slits.size > 0:
                good_fit = np.ones(self._nslit, dtype=np.bool)
                good_fit[obad_slits] = False
                bad_slits = self.cross_match(good_fit)
                if np.array_equal(bad_slits, obad_slits):
                    break
                obad_slits = bad_slits.copy()
                cntr += 1
                if cntr > 100:
                    msgs.warn("Breaking while loop before convergence. Check the wavelength solution!")
                    break

        # With the updates to the fits of each slit, determine the final fit, and save the QA
        for slit in range(self._nslit):
            if slit not in self._ok_mask:
                continue
            if self._all_patt_dict[str(slit)] is None:
                continue
            # Save the QA for the best solution
            slittxt = '_Slit{0:03d}'.format(slit+1)
            use_tcent = self.get_use_tcent(self._all_patt_dict[str(slit)]['sign'],
                                           arr=self._detections[str(slit)], weak=True)
            if self._outroot is not None:
                # Write IDs
                out_dict = dict(pix=use_tcent, IDs=self._all_patt_dict[str(slit)]['IDs'])
                jdict = ltu.jsonify(out_dict)
                ltu.savejson(self._outroot + slittxt + '.json', jdict, easy_to_read=True, overwrite=True)
                msgs.info("Wrote: {:s}".format(self._outroot + slittxt + '.json'))

                # Plot
                tmp_list = vstack([self._line_lists, self._unknwns])
                qa.match_qa(self._spec[:, slit], use_tcent, tmp_list,
                            self._all_patt_dict[str(slit)]['IDs'], self._all_patt_dict[str(slit)]['scores'],
                            self._outroot + slittxt + '.pdf')
                msgs.info("Wrote: {:s}".format(self._outroot + slittxt + '.pdf'))
            # Perform the final fit for the best solution
            best_final_fit = self.fit_slit(slit, self._all_patt_dict[str(slit)], tcent=use_tcent,
                                           outroot=self._outroot, slittxt=slittxt)
            self._all_final_fit[str(slit)] = copy.deepcopy(best_final_fit)

        # Print the final report of all lines
        self.final_report()
        return

    def run_kdtree(self, polygon=4, detsrch=4):
        """ KD Tree algorithm to wavelength calibrate spectroscopic data.
        Currently, this is only designed for ThAr lamp spectra. See the
        'run_brute' function if you want to calibrate longslit spectra.
        """

        # Load the linelist KD Tree
        lsttree, lindex = waveio.load_tree(polygon=polygon)

        # Set the search error to be 1 pixel
        err = 1.0 / self._npix

        self._all_patt_dict = {}
        self._all_final_fit = {}
        good_fit = np.zeros(self._nslit, dtype=np.bool)
        for slit in range(self._nslit):
            if slit not in self._ok_mask:
                continue
            # Detect lines, and decide which tcent to use
            self._all_tcent, self._cut_tcent, self._icut =\
                wvutils.arc_lines_from_spec(self._spec[:, slit], min_ampl=self._min_ampl)
            self._all_tcent_weak, self._cut_tcent_weak, self._icut_weak =\
                wvutils.arc_lines_from_spec(self._spec[:, slit], min_ampl=self._lowest_ampl)
            if self._all_tcent.size == 0:
                msgs.warn("No lines to identify in slit {0:d}!".format(slit))
                continue
            best_patt_dict, best_final_fit = None, None

            use_tcentp = self.get_use_tcent(1)
            use_tcentm = self.get_use_tcent(-1)
            if use_tcentp.size < detsrch:
                if self._verbose:
                    msgs.info("Not enough lines to test this solution, will attempt another.")
                return None, None

            # Create a detlines KD Tree
            maxlinear = 0.25*self._npix
            if polygon == 3:
                msgs.info("Generating patterns for a trigon")
                patternp, indexp = kdtree_generator.trigon(use_tcentp, detsrch, maxlinear)
                patternm, indexm = kdtree_generator.trigon(use_tcentm, detsrch, maxlinear)
            elif polygon == 4:
                msgs.info("Generating patterns for a tetragon")
                patternp, indexp = kdtree_generator.tetragon(use_tcentp, detsrch, maxlinear)
                patternm, indexm = kdtree_generator.tetragon(use_tcentm, detsrch, maxlinear)
            elif polygon == 5:
                msgs.info("Generating patterns for a pentagon")
                patternp, indexp = kdtree_generator.pentagon(use_tcentp, detsrch, maxlinear)
                patternm, indexm = kdtree_generator.pentagon(use_tcentm, detsrch, maxlinear)
            elif polygon == 6:
                msgs.info("Generating patterns for a hexagon")
                patternp, indexp = kdtree_generator.hexagon(use_tcentp, detsrch, maxlinear)
                patternm, indexm = kdtree_generator.hexagon(use_tcentm, detsrch, maxlinear)
            else:
                msgs.warn("Patterns can only be generated with 3 <= polygon <= 6")
                return None

            dettreep = cKDTree(patternp, leafsize=30)
            dettreem = cKDTree(patternm, leafsize=30)

            # Query the detections tree
            msgs.info("Querying KD tree patterns (slit {0:d}/{1:d})".format(slit+1, self._nslit))
            resultp = dettreep.query_ball_tree(lsttree, r=err)
            resultm = dettreem.query_ball_tree(lsttree, r=err)

            msgs.info("Identifying wavelengths")
            psols = self.results_kdtree(use_tcentp, resultp, indexp, lindex)
            msols = self.results_kdtree(use_tcentm, resultm, indexm, lindex)
            patt_dict, final_fit = self.solve_slit(slit, psols, msols)

            pdb.set_trace()

        # Return
        return all_patt_dict, all_final_fit

    def cross_match(self, good_fit):
        """Cross-correlate the spectra across all slits to ID all of the lines.
        good_fit : ndarray (bool)
          Indicates which slits are deemed to be a good fit (although, they need not necessarily be a good fit).
        """
        # Steps:
        # Check that all of the "good" slits are indeed good
        # For all of the bad slits, cross-correlate against all of the good slits to label each line
        # For all newly labeled lines, create a patt_dict of these labeled lines
        # Perform a final fit on these lines

        # First, sort spectra according to increasing central wavelength
        ngd = good_fit.sum()
        idx_gd = np.zeros(ngd, dtype=np.int)
        wvc_gd = np.zeros(ngd, dtype=np.float)
        dsp_gd = np.zeros(ngd, dtype=np.float)
        cntr = 0
        for slit in range(self._nslit):
            if good_fit[slit]:
                idx_gd[cntr] = slit
                wvc_gd[cntr] = self._all_patt_dict[str(slit)]["bwv"]
                dsp_gd[cntr] = self._all_patt_dict[str(slit)]["bdisp"]
                cntr += 1
        srt = np.argsort(wvc_gd)
        sort_idx = idx_gd[srt]
        sort_wvc = wvc_gd[srt]
        sort_dsp = dsp_gd[srt]

        # Cross correlate all good spectra with each other, in order of wavelength
        ncrco = np.arange(sort_idx.size).sum()
        #ccor_val = np.zeros(ncrco)
        dwvc_val = np.zeros(ncrco)
        slit_ids = np.zeros((ncrco, 2), dtype=np.int)
        cntr = 0
        for gd in range(0, sort_idx.size-1):
            for gc in range(gd+1, sort_idx.size):
                corr = scipy.signal.correlate(self._spec[:, sort_idx[gd]], self._spec[:, sort_idx[gc]], mode='same')
                amax = np.argmax(corr)
                #ccor_val[cntr] = (amax - self._spec.shape[0] // 2)
                dwvc_val[cntr] = (sort_wvc[gc]-sort_wvc[gd]) / (0.5*(sort_dsp[gc]+sort_dsp[gd])) - (amax - self._spec.shape[0] // 2)
                slit_ids[cntr, 0] = gd
                slit_ids[cntr, 1] = gc
                cntr += 1

        # Identify the good orders
        sigrej = 3.0
        mad = 1.4826 * np.median(np.abs(dwvc_val))
        gdmsk = np.where(np.abs(dwvc_val) < sigrej * mad)[0]
        for ii in range(100):  # Limit to 100 iterations - this will likely never be reached...
            ogdmsk = gdmsk.copy()
            mad = 1.4826 * np.median(np.abs(dwvc_val[gdmsk]))
            gdmsk = np.where(np.abs(dwvc_val) < sigrej*mad)[0]
            if np.array_equal(gdmsk, ogdmsk):
                break

        debug = False
        if debug:
            from matplotlib import pyplot as plt
            xplt = np.arange(dwvc_val.size)
            plt.plot(xplt, dwvc_val, 'rx')
            plt.plot(xplt[gdmsk], dwvc_val[gdmsk], 'bo')
            plt.plot([0.0,1000.0],[0.0,0.0], 'b-')
            plt.show()

        # Catalogue the good and bad slits
        good_slits = np.sort(sort_idx[np.unique(slit_ids[gdmsk, :].flatten())])
        bad_slits = np.setdiff1d(np.arange(self._nslit), good_slits, assume_unique=True)
        # Get the sign (i.e. if pixels correlate/anticorrelate with wavelength)
        # and dispersion (A/pix). Assume these are the same for all slits
        sign = self._all_patt_dict[str(good_slits[0])]['sign']
        disp = self._all_patt_dict[str(good_slits[0])]['bdisp']

        # For all of the bad slits, estimate some line wavelengths
        new_bad_slits = np.array([], dtype=np.int)
        for bs in bad_slits:
            if bs not in self._ok_mask:
                continue
            bsdet = self.get_use_tcent(sign, arr=self._detections[str(bs)], weak=True)
            lindex = np.array([], dtype=np.int)
            dindex = np.array([], dtype=np.int)
            wcen = np.zeros(good_slits.size)
            for cntr, gs in enumerate(good_slits):
                # Match the peaks between the two spectra.
                # spec_gs_adj is the stretched spectrum
                stretch, shift = wvutils.match_peaks(self._spec[:, bs], self._spec[:, gs])
                if stretch is None:
                    continue
                # Estimate wcen for this slit
                wcen[cntr] = self._all_patt_dict[str(gs)]['bwv'] - shift*disp
                # For each peak in the gs spectrum, identify the corresponding peaks in the bs spectrum
                strfact = (self._npix + stretch - 1)/(self._npix - 1)
                gsdet = self.get_use_tcent(sign, arr=self._detections[str(gs)], weak=True)
                gsdet_ss = shift + gsdet * strfact
                debug = False
                if debug:
                    from matplotlib import pyplot as plt
                    xplt = np.arange(self._npix)
                    plt.plot(xplt, self._spec[:, bs], 'k-', drawstyle='steps')
                    plt.plot(bsdet, np.zeros(bsdet.size), 'ro')
                    plt.plot(gsdet_ss, 0.01*np.max(self._spec[:, bs])*np.ones(gsdet_ss.size), 'bo')
                    plt.show()
                    pdb.set_trace()
                # Calculate wavelengths for the gsdet detections
                fitc = self._all_final_fit[str(gs)]['fitc']
                xfit = gsdet/(self._npix - 1)
                fitfunc = self._all_final_fit[str(gs)]['function']
                fmin, fmax = 0.0, 1.0
                wvval = utils.func_val(fitc, xfit, fitfunc, minv=fmin, maxv=fmax)
                for dd in range(bsdet.size):
                    pdiff = np.abs(bsdet[dd]-gsdet_ss)
                    bstpx = np.argmin(pdiff)
                    # If a match is found within 2 pixels, consider this a successful match
                    if pdiff[bstpx] < 2.0:
                        bstwv = np.abs(self._wvdata - wvval[bstpx])
                        if bstwv[np.argmin(bstwv)] > 2.0*disp:
                            # This is probably not a good match
                            continue
                        lindex = np.append(lindex, np.argmin(bstwv))
                        dindex = np.append(dindex, dd)
            # Finalize the best guess of each line
            # Initialise the patterns dictionary
            patt_dict = dict(acceptable=False, nmatch=0, ibest=-1, bwv=0., min_ampl=self._min_ampl,
                             mask=np.zeros(bsdet.size, dtype=np.bool))
            patt_dict['sign'] = sign
            patt_dict['bwv'] = np.median(wcen[wcen != 0.0])
            patt_dict['bdisp'] = disp
            patterns.solve_triangles(bsdet, self._wvdata, dindex, lindex, patt_dict)
            # Check if a solution was found
            if not patt_dict['acceptable']:
                new_bad_slits = np.append(new_bad_slits, bs)
                continue
            final_fit = self.fit_slit(bs, patt_dict, tcent=bsdet)
            if final_fit is None:
                # This pattern wasn't good enough
                new_bad_slits = np.append(new_bad_slits, bs)
                continue
            if final_fit['rms'] > self._rms_threshold:
                msgs.warn('---------------------------------------------------' + msgs.newline() +
                          'Cross-match report for slit {0:d}/{1:d}:'.format(bs + 1, self._nslit) + msgs.newline() +
                          '  Poor RMS ({0:.3f})! Will try cross matching iteratively'.format(final_fit['rms']) + msgs.newline() +
                          '---------------------------------------------------')
                # Store this result in new_bad_slits, so the iteration can be performed,
                # but make sure to store the result, as this might be the best possible.
                new_bad_slits = np.append(new_bad_slits, bs)
            self._all_patt_dict[str(bs)] = copy.deepcopy(patt_dict)
            self._all_final_fit[str(bs)] = copy.deepcopy(final_fit)
            debug = False
            if debug:
                from matplotlib import pyplot as plt
                xplt = np.linspace(0.0, 1.0, self._npix)
                yplt = utils.func_val(final_fit['fitc'], xplt, 'legendre', minv=0.0, maxv=1.0)
                plt.plot(final_fit['xfit'], final_fit['yfit'], 'bx')
                plt.plot(xplt, yplt, 'r-')
                plt.show()
                pdb.set_trace()
        return new_bad_slits

    def get_use_tcent(self, corr, arr=None, weak=False):
        """Set if pixels correlate with wavelength (corr==1) or anticorrelate (corr=-1)
        """
        # Decide which array to use
        if arr is None:
            if weak:
                arr = self._all_tcent_weak.copy()
            else:
                arr = self._all_tcent.copy()
        # Return the appropriate tcent
        if corr == 1:
            return arr
        else:
            return (self._npix - 1.0) - arr[::-1]

    def results_brute(self, poly=3, pix_tol=0.5, detsrch=5, lstsrch=5):

        # Import the pattern matching algorithms
        if poly == 3:
            from pypeit.core.wavecal.patterns import triangles as generate_patterns
        elif poly == 4:
            from pypeit.core.wavecal.patterns import quadrangles as generate_patterns
        else:
            msgs.warn("Pattern matching is only available for trigons and tetragons.")
            return None, None

        # Test if there are enough lines to generate a solution
        use_tcent = self.get_use_tcent(1)
        if use_tcent.size < lstsrch or use_tcent.size < detsrch:
            if self._verbose:
                msgs.info("Not enough lines to test this solution, will attempt another.")
            return None, None

        if self._verbose:
            msgs.info("Begin pattern matching")
        # First run pattern recognition assuming pixels correlate with wavelength
        dindexp, lindexp, wvcenp, dispsp = generate_patterns(use_tcent, self._wvdata, self._npix,
                                                             detsrch, lstsrch, pix_tol)
        # Now run pattern recognition assuming pixels correlate with wavelength
        use_tcent = self.get_use_tcent(-1)
        dindexm, lindexm, wvcenm, dispsm = generate_patterns(use_tcent, self._wvdata, self._npix,
                                                             detsrch, lstsrch, pix_tol)
        return (dindexp, lindexp, wvcenp, dispsp,), (dindexm, lindexm, wvcenm, dispsm,)

    def results_kdtree(self, use_tcent, res, dindex, lindex, ordfit=2):
        # Assign wavelengths to each pixel
        waveid = [np.array([]) for xx in use_tcent]
        nrows = len(res)
        ncols = sum(map(len, res))
        nindx = dindex.shape[1]
        wvdisp = np.zeros(ncols)
        wvcent = np.zeros(ncols)
        dind = np.zeros((ncols, nindx), dtype=np.int)
        lind = np.zeros((ncols, nindx), dtype=np.int)
        cnt = 0
        for x in range(nrows):
            for y in range(len(res[x])):
                dx = use_tcent[dindex[x, -1]] - use_tcent[dindex[x, 0]]
                dp = self._wvdata[lindex[res[x][y], -1]] - self._wvdata[lindex[res[x][y], 0]]
                try:
                    null, cgrad = utils.robust_polyfit(use_tcent[dindex[x, :]], self._wvdata[lindex[res[x][y], :]],
                                                       1, sigma=2.0)
                    wvdisp[cnt] = cgrad[1]
                except:
                    wvdisp[cnt] = (dp / dx)

                coeff = np.polyfit(use_tcent[dindex[x, :]], self._wvdata[lindex[res[x][y]]], ordfit)
                wvcent[cnt] = np.polyval(coeff, self._npix / 2.0)
                dind[cnt, :] = dindex[x, :]
                lind[cnt, :] = lindex[res[x][y], :]
                cnt += 1
        return dind, lind, wvcent, wvdisp

    def solve_slit(self, slit, psols, msols, nstore=1, nselw=3, nseld=3):

        # Extract the solutions
        dindexp, lindexp, wvcenp, dispsp = psols
        dindexm, lindexm, wvcenm, dispsm = msols

        # Remove any invalid results from correlate
        ww = np.where((self._binw[0] < wvcenp) & (wvcenp < self._binw[-1]) &
                      (10.0 ** self._bind[0] < dispsp) & (dispsp < 10.0 ** self._bind[-1]))
        dindexp = dindexp[ww[0], :]
        lindexp = lindexp[ww[0], :]
        dispsp = dispsp[ww]
        wvcenp = wvcenp[ww]

        # Remove any invalid results from anticorrelate
        ww = np.where((self._binw[0] < wvcenm) & (wvcenm < self._binw[-1]) &
                      (10.0 ** self._bind[0] < dispsm) & (dispsm < 10.0 ** self._bind[-1]))
        dindexm = dindexm[ww[0], :]
        lindexm = lindexm[ww[0], :]
        dispsm = dispsm[ww]
        wvcenm = wvcenm[ww]

        # Construct the histograms
        histimgp, xed, yed = np.histogram2d(wvcenp, np.log10(dispsp), bins=[self._binw, self._bind])
        histimgm, xed, yed = np.histogram2d(wvcenm, np.log10(dispsm), bins=[self._binw, self._bind])
        #histimgp = gaussian_filter(histimgp, 3)
        #histimgm = gaussian_filter(histimgm, 3)
        histimg = histimgp - histimgm
        sm_histimg = gaussian_filter(histimg, [3, 15])

        #histpeaks = patterns.detect_2Dpeaks(np.abs(sm_histimg))
        histpeaks = patterns.detect_2Dpeaks(np.abs(histimg))

        # Find the indices of the nstore largest peaks
        bidx = np.unravel_index(np.argpartition(np.abs(histpeaks*histimg), -nstore, axis=None)[-nstore:], histimg.shape)

        debug = False
        if debug:
            from matplotlib import pyplot as plt
            plt.clf()
            extent = [self._binw[0], self._binw[-1], self._bind[0], self._bind[-1]]
            plt.imshow((np.abs(sm_histimg[:, ::-1].T)), extent=extent, aspect='auto')
            #plt.imshow(histimg[:, ::-1].T, extent=extent, aspect='auto')
            plt.plot(self._binw[bidx[0]], self._bind[bidx[1]], 'r+')
            #plt.axvline(self._binw[self._bidx[0]], color='r', linestyle='--')
            #plt.axhline(self._bind[self._bidx[1]], color='r', linestyle='--')
            plt.show()
            if False:
                pdb.set_trace()
                plt.clf()
                plt.imshow(histimgp[:, ::-1].T, extent=extent, aspect='auto')
                plt.show()

        # Get the peak value of central wavelength and dispersion
        allwcen = self._binw[bidx[0]]
        alldisp = self._bind[bidx[1]]
        allhnum = np.abs(histimg[bidx])

        # Find all good solutions
        bestlist = []
        for idx in range(nstore):
            # Select all solutions around the best solution within a square of side 2*nsel
            wlo = self._binw[max(0, bidx[0][idx] - nselw)]
            whi = self._binw[min(self._ngridw - 1, bidx[0][idx] + nselw)]
            dlo = 10.0 ** self._bind[max(0, bidx[1][idx] - nseld)]
            dhi = 10.0 ** self._bind[min(self._ngridd - 1, bidx[1][idx] + nseld)]
            if histimgp[bidx][idx] > histimgm[bidx][idx]:
                wgd = np.where((wvcenp > wlo) & (wvcenp < whi) & (dispsp > dlo) & (dispsp < dhi))
                dindex = dindexp[wgd[0], :].flatten()
                lindex = lindexp[wgd[0], :].flatten()
                sign = +1
            else:
                wgd = np.where((wvcenm > wlo) & (wvcenm < whi) & (dispsm > dlo) & (dispsm < dhi))
                dindex = dindexm[wgd[0], :].flatten()
                lindex = lindexm[wgd[0], :].flatten()
                sign = -1
            # Store relevant values in an array to solve for best solution
            bestlist.append([allwcen[idx], alldisp[idx], allhnum[idx], sign, dindex, lindex])

        if self._verbose:
            msgs.info("Fitting the wavelength solution for each slit")
        patt_dict, final_dict = None, None
        for idx in range(nstore):
            # Solve the patterns
            tpatt_dict = self.solve_patterns(bestlist[idx])
            if tpatt_dict is None:
                # This pattern wasn't good enough
                continue
            # Fit the full set of lines with the derived patterns
            tfinal_dict = self.fit_slit(slit, tpatt_dict)
            if tfinal_dict is None:
                # This pattern wasn't good enough
                continue
            # Check if this solution is better than the last
            if patt_dict is None:
                # First time a fit is found
                patt_dict, final_dict = tpatt_dict, tfinal_dict
                continue
            elif tfinal_fit['rms'] < self._rms_threshold:
                # Has a better fit been identified (i.e. more lines ID)?
                if len(tfinal_fit['xfit']) > len(final_fit['xfit']):
                    patt_dict, final_dict = copy.deepcopy(tpatt_dict), copy.deepcopy(tfinal_dict)
        return patt_dict, final_dict

    def solve_patterns(self, bestlist):

        # Obtain a full list of indices that are consistent with the maximum value
        wcen, dcen, sign, dindex, lindex = bestlist[0], bestlist[1], bestlist[3], bestlist[4], bestlist[5]

        # Find the favoured sign and only use those values
        use_tcent = self.get_use_tcent(sign)
        if sign == +1:
            signtxt = "correlate"
        else:
            signtxt = "anticorrelate"

        # Initialise the patterns dictionary
        patt_dict = dict(acceptable=False, nmatch=0, ibest=-1, bwv=0., min_ampl=self._min_ampl,
                         mask=np.zeros(use_tcent.size, dtype=np.bool))
        patterns.solve_triangles(use_tcent, self._wvdata, dindex, lindex, patt_dict)
        # Check if a solution was found
        if not patt_dict['acceptable']:
            return None

        # Fill in the patterns dictionary
        patt_dict['sign'] = sign
        patt_dict['bwv'] = wcen
        patt_dict['bdisp'] = 10.0 ** dcen

        # Check that a solution has been found
        if patt_dict['nmatch'] == 0 and self._verbose:
            msgs.info('---------------------------------------------------' + msgs.newline() +
                      'Initial report:' + msgs.newline() +
                      '  No matches! Try another algorithm' + msgs.newline() +
                      '---------------------------------------------------')
            return None
        elif self._verbose:
            # Report
            msgs.info('---------------------------------------------------' + msgs.newline() +
                      'Initial report:' + msgs.newline() +
                      '  Pixels {:s} with wavelength'.format(signtxt) + msgs.newline() +
                      '  Number of lines recovered    = {:d}'.format(self._all_tcent.size) + msgs.newline() +
                      '  Number of lines analyzed     = {:d}'.format(use_tcent.size) + msgs.newline() +
                      '  Number of acceptable matches = {:d}'.format(patt_dict['nmatch']) + msgs.newline() +
                      '  Best central wavelength      = {:g}A'.format(patt_dict['bwv']) + msgs.newline() +
                      '  Best dispersion              = {:g}A/pix'.format(patt_dict['bdisp']) + msgs.newline() +
                      '---------------------------------------------------')
        return patt_dict

    def fit_slit(self, slit, patt_dict, outroot=None, slittxt="Slit", tcent=None):
        # Perform final fit to the line IDs
        NIST_lines = self._line_lists['NIST'] > 0
        ifit = np.where(patt_dict['mask'])[0]

        if outroot is not None:
            plot_fil = outroot + slittxt + '_fit.pdf'
        else:
            plot_fil = None
        # Purge UNKNOWNS from ifit
        imsk = np.ones(len(ifit), dtype=np.bool)
        for kk, idwv in enumerate(np.array(patt_dict['IDs'])[ifit]):
            if np.min(np.abs(self._line_lists['wave'][NIST_lines]-idwv)) > 0.01:
                imsk[kk] = False
        ifit = ifit[imsk]
        # Allow for weaker lines in the fit
        if tcent is None:
            tcent = self.get_use_tcent(patt_dict['sign'], weak=True)
        # Fit
        try:
            final_fit = fitting.iterative_fitting(self._spec[:, slit], tcent, ifit,
                                                  np.array(patt_dict['IDs'])[ifit], self._line_lists[NIST_lines],
                                                  patt_dict['bdisp'], plot_fil=plot_fil, verbose=self._verbose,
                                                  aparm=self._fit_parm)
        except TypeError:
            # A poor fitting result, this can be ignored.
            return None

        if plot_fil is not None:
            print("Wrote: {:s}".format(plot_fil))

        # Return
        return final_fit

    def final_report(self):
        """Print out the final report of the wavelength calibration"""
        msgs.info('###################################################')
        for slit in range(self._nslit):
            # Prepare a message for bad wavelength solutions
            badmsg = '---------------------------------------------------' + msgs.newline() +\
                     'Final report for slit {0:d}/{1:d}:'.format(slit + 1, self._nslit) + msgs.newline() +\
                     '  Wavelength calibration not performed!'
            if slit not in self._ok_mask:
                msgs.warn(badmsg)
                continue
            if self._all_patt_dict[str(slit)] is None:
                msgs.warn(badmsg)
                continue
            st = str(slit)
            if self._all_patt_dict[st]['sign'] == +1:
                signtxt = 'correlate'
            else:
                signtxt = 'anitcorrelate'
            # Report
            centwave = utils.func_val(self._all_final_fit[st]['fitc'], 0.5,
                                      self._all_final_fit[st]['function'], minv=0.0, maxv=1.0)
            tempwave = utils.func_val(self._all_final_fit[st]['fitc'], 0.5 + 1.0/self._npix,
                                      self._all_final_fit[st]['function'], minv=0.0, maxv=1.0)
            centdisp = abs(centwave-tempwave)
            msgs.info('---------------------------------------------------' + msgs.newline() +
                      'Final report for slit {0:d}/{1:d}:'.format(slit + 1, self._nslit) + msgs.newline() +
                      '  Pixels {:s} with wavelength'.format(signtxt) + msgs.newline() +
                      '  Number of lines recovered    = {:d}'.format(self._detections[st].size) + msgs.newline() +
                      '  Number of lines analyzed     = {:d}'.format(len(self._all_final_fit[st]['xfit'])) + msgs.newline() +
                      '  Central wavelength           = {:g}A'.format(centwave) + msgs.newline() +
                      '  Central dispersion           = {:g}A/pix'.format(centdisp) + msgs.newline() +
                      '  Final RMS of fit             = {:g}'.format(self._all_final_fit[st]['rms']))
        msgs.info('###################################################')
        return
