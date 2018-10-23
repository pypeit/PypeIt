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
import numba as nb
import numpy as np
import pdb

from pypeit.par import pypeitpar
from pypeit.core.wavecal import kdtree_generator
from pypeit.core.wavecal import waveio
from pypeit.core.wavecal import patterns
from pypeit.core.wavecal import fitting
from pypeit.core.wavecal import wvutils
from pypeit.core.wavecal import qa
from pypeit.core import pca
from pypeit import utils

from pypeit import msgs
from pypeit import debugger
from matplotlib import pyplot as plt


def basic(spec, lines, wv_cen, disp, min_nsig=20.,nonlinear_counts = 1e10,
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
    min_nsig : float
      Minimum significance of the arc lines that will be used in the fit
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
    all_tcent, cut_tcent, icut = wvutils.arc_lines_from_spec(spec, min_nsig=min_nsig, nonlinear_counts = nonlinear_counts)

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


def semi_brute(spec, lines, wv_cen, disp, min_nsig=30., nonlinear_counts = 1e10,
               outroot=None, debug=False, do_fit=True, verbose=False,
               min_nmatch=3, lowest_nsig=20.,
               match_toler=3.0, func='legendre', n_first=2, sigrej_first=2.0, n_final=4, sigrej_final=3.0):
    """
    Parameters
    ----------
    spec
    lines
    wv_cen
    disp
    siglev
    min_nsig
    outroot
    debug
    do_fit
    verbose
    min_nmatch
    lowest_nsig

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
    all_tcent, cut_tcent, icut, _ = wvutils.arc_lines_from_spec(spec, min_nsig=min_nsig, nonlinear_counts = nonlinear_counts)

    # Best
    best_dict = dict(nmatch=0, ibest=-1, bwv=0., min_nsig=min_nsig, unknown=False,
                     pix_tol=1, nsig=min_nsig)

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
            # Lower minimum significance
            nsig = min_nsig
            #pdb.set_trace()
            while(best_dict['nmatch'] < min_nmatch):
                nsig /= 2.
                if nsig < lowest_nsig:
                    break
                all_tcent, cut_tcent, icut, _ = wvutils.arc_lines_from_spec(spec, min_nsig=nsig, nonlinear_counts = nonlinear_counts)
                patterns.scan_for_matches(wv_cen, disp, npix, cut_tcent, wvdata,
                                          best_dict=best_dict, pix_tol=pix_tol)#, nsig=nsig)

        #if debug:
        #    pdb.set_trace()
        # Save linelist?
        if best_dict['nmatch'] > sav_nmatch:
            best_dict['line_list'] = tot_list.copy()
            best_dict['unknown'] = unknown
            best_dict['nsig'] = nsig
            best_dict['pix_tol'] = pix_tol

    # Try to pick up some extras by turning off/on unknowns
    if best_dict['unknown']:
        tot_list = line_listsarc_lines_from_spec
    else:
        tot_list = vstack([line_lists,unknwns])
    wvdata = np.array(tot_list['wave'].data) # Removes mask if any
    wvdata.sort()
    tmp_dict = copy.deepcopy(best_dict)
    tmp_dict['nmatch'] = 0
    patterns.scan_for_matches(best_dict['bwv'], disp, npix, cut_tcent, wvdata,
                              best_dict=tmp_dict, pix_tol=best_dict['pix_tol'], wvoff=1.)
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
        qa.match_qa(spec, cut_tcent, tmp_list, best_dict['IDs'], best_dict['scores'], outfile = outroot+'.pdf')
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
        all_tcent, weak_cut_tcent, icut = wvutils.arc_lines_from_spec(spec, min_nsig=lowest_nsig, nonlinear_counts = nonlinear_counts)
        add_weak = []
        for weak in weak_cut_tcent:
            if np.min(np.abs(cut_tcent-weak)) > 5.:
                add_weak += [weak]
        if len(add_weak) > 0:
            cut_tcent = np.concatenate([cut_tcent, np.array(add_weak)])
        # Fit
        final_fit = fitting.iterative_fitting(spec, cut_tcent, ifit,
                                              np.array(best_dict['IDs'])[ifit], line_lists[NIST_lines],
                                              disp, plot_fil=plot_fil, verbose=verbose,
                                              match_toler=match_toler, func=func, n_first=n_first, sigrej_first=sigrej_first,
                                              n_final=n_final,sigrej_final=sigrej_final)
        if plot_fil is not None:
            print("Wrote: {:s}".format(plot_fil))

    # Return
    return best_dict, final_fit

class General:
    """ General algorithm to wavelength calibrate spectroscopic data

    Parameters
    ----------
    spec : ndarray
      2D array of arcline spectra (nspec,nslit)
    lines : list
      List of arc lamps on
    par:
    ok_mask : ndarray
      Array of good slits
    islinelist : bool
      Is lines a linelist (True), or a list of ions (False)
    outroot : str, optional
      Name of output file
    debug : bool
      Used to debug the algorithm
    verbose : bool
      If True, the final fit will print out more detail as the RMS is refined,
      and lines are rejected. This is mostly helpful for developing the algorithm.
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

    def __init__(self, spec, par = None, ok_mask=None, islinelist=False, outroot=None, debug = False, verbose=False,
                 binw=None, bind=None, nstore=1, use_unknowns=True):

        # Set some default parameters
        self._spec = spec
        self._par = pypeitpar.WavelengthSolutionPar() if par is None else par
        self._lines = self._par['lamps']
        self._npix, self._nslit = spec.shape
        self._nstore = nstore
        self._binw = binw
        self._bind = bind


        # Mask info
        if ok_mask is None:
            self._ok_mask = np.arange(self._nslit)
        else:
            self._ok_mask = ok_mask
        self._bad_slits = []  # List of bad slits

        # Set the input parameters
        self._nonlinear_counts = self._par['nonlinear_counts']
        self._min_nsig = self._par['min_nsig']
        self._lowest_nsig = self._par['lowest_nsig']
        self._rms_threshold = self._par['rms_threshold']
        self._match_toler = self._par['match_toler']
        self._func = self._par['func']
        self._n_first= self._par['n_first']
        self._sigrej_first= self._par['sigrej_first']
        self._n_final= self._par['n_final']
        self._sigrej_final= self._par['sigrej_final']

        self._use_unknowns = use_unknowns
        self._islinelist = islinelist

        self._outroot = outroot

        self._debug = debug
        self._verbose = verbose

        # Load the linelist to be used for pattern matching
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

        # Find the wavelength solution!
        # KD Tree algorithm only works for ThAr - check first that this is what is being used
        self._thar = False
        if 'ThAr' in self._lines and len(self._lines) == 1:
            self._thar = True
            # Set up the grids to be used for pattern matching
            self.set_grids(ngridw=5000, ngridd=1000)
            msgs.info("Using KD Tree pattern matching algorithm to wavelength calibrate")
            self.run_kdtree()
        else:
            # Set up the grids to be used for pattern matching
            self.set_grids()
            msgs.info("Using brute force pattern matching algorithm to wavelength calibrate")
            self.run_brute()

    def get_results(self):
        return copy.deepcopy(self._all_patt_dict), copy.deepcopy(self._all_final_fit)

    def set_grids(self, ngridw = 200, ngridd=2000): #ngridw = 200, ngridd=2000):
        # Set the wavelength grid
        if self._binw is None:
            # Ideally, you want binw to roughly sample the A/pix of the spectrograph
            self._ngridw = ngridw
            self._binw = np.linspace(np.min(self._wvdata), np.max(self._wvdata), self._ngridw)
        else:
            self._ngridw = self._binw.size
        # Set the dispersion grid
        if self._bind is None:
            self._ngridd = ngridd
            #self._bind = np.linspace(-3.0, 1.0, self._ngridd)
            # JFH I have no idea why this goes down as low as -3.0. 3000/10^(-3.0) would be R ~ 3e6. No spectrograph
            # has a dispersion that high. I'm changing this to be -1.5 which would be R ~ 100,000 at 3000A. In this regime
            # one would anyway use the ThAr routine. I'm rasing
            # the upper limit to be 2.0 to handle low-resolution data (i.e in the near-IR 2.5e4/100 = R ~ 250

            self._bind = np.linspace(-1.5, 2.0, self._ngridd)
        else:
            self._ngridd = self._bind.size
        return

    def run_brute_loop(self, slit, arrerr=None, wavedata=None):
        # Set the parameter space that gets searched
        rng_poly = [3, 4]            # Range of algorithms to check (only trigons+tetragons are supported)
        rng_list = range(3, 6)       # Number of lines to search over for the linelist
        rng_detn = range(3, 6)       # Number of lines to search over for the detected lines
        rng_pixt = [1.0]             # Pixel tolerance
        idthresh = 0.5               # Criteria for early return (at least this fraction of lines must have
                                     # an ID on either side of the spectrum)

        best_patt_dict, best_final_fit = None, None
        # Loop through parameter space
        for poly in rng_poly:
            for detsrch in rng_detn:
                for lstsrch in rng_list:
                    for pix_tol in rng_pixt:
                        psols, msols = self.results_brute(poly=poly, pix_tol=pix_tol, detsrch=detsrch, lstsrch=lstsrch,
                                                          wavedata=wavedata, arrerr=arrerr)
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
                            # Decide if an early return is acceptable
                            nlft = np.sum(best_final_fit['tcent'] < best_final_fit['xnorm']/2.0)
                            nrgt = best_final_fit['tcent'].size-nlft
                            if np.sum(best_final_fit['xfit'] < 0.5)/nlft > idthresh and\
                                np.sum(best_final_fit['xfit'] >= 0.5) / nrgt > idthresh:
                                # At least half of the lines on either side of the spectrum have been identified
                                return best_patt_dict, best_final_fit

        return best_patt_dict, best_final_fit

    def run_brute(self):
        """Run through the parameter space and determine the best solution
        """

        self._all_patt_dict = {}
        self._all_final_fit = {}
        good_fit = np.zeros(self._nslit, dtype=np.bool)
        self._detections = {}
        for slit in range(self._nslit):
            if slit not in self._ok_mask:
                continue
            # Detect lines, and decide which tcent to use
            self._all_tcent, self._all_ecent, self._cut_tcent, self._icut =\
                wvutils.arc_lines_from_spec(self._spec[:, slit].copy(), min_nsig=self._min_nsig, nonlinear_counts = self._nonlinear_counts)
            self._all_tcent_weak, self._all_ecent_weak, self._cut_tcent_weak, self._icut_weak =\
                wvutils.arc_lines_from_spec(self._spec[:, slit].copy(), min_nsig=self._lowest_nsig, nonlinear_counts = self._nonlinear_counts)
            if self._all_tcent.size == 0:
                msgs.warn("No lines to identify in slit {0:d}!".format(slit))
                self._detections[str(slit)] = [None,None]
                continue
            # Setup
            #self._detections[str(slit)] = [self._all_tcent_weak.copy(), self._all_ecent_weak.copy()]
            self._detections[str(slit)] = [self._all_tcent_weak[self._icut_weak].copy(),
                                           self._all_ecent_weak[self._icut_weak].copy()]
            # Run it
            best_patt_dict, best_final_fit = self.run_brute_loop(slit)

            # Print preliminary report
            good_fit[slit] = self.report_prelim(slit, best_patt_dict, best_final_fit)

        # Now that all slits have been inspected, cross match to generate a
        # master list of all lines in every slit, and refit all spectra
        if self._nslit > 1:
            msgs.info('Checking wavelength solution by cross-correlating with all slits')

            msgs.info('Cross-correlation iteration #1')
            obad_slits = self.cross_match(good_fit)
            cntr = 2
            while obad_slits.size > 0:
                msgs.info('Cross-correlation iteration #{:d}'.format(cntr))
                good_fit = np.ones(self._nslit, dtype=np.bool)
                good_fit[obad_slits] = False
                bad_slits = self.cross_match(good_fit)
                if np.array_equal(bad_slits, obad_slits):
                    break
                obad_slits = bad_slits.copy()
                cntr += 1
                if cntr > 10:
                    msgs.warn("Breaking while loop before convergence. Check the wavelength solution!")
                    break

        # With the updates to the fits of each slit, determine the final fit, and save the QA
        self.finalize_fit()

        # Print the final report of all lines
        self.report_final()
        return

    def run_kdtree(self, polygon=4, detsrch=7, lstsrch=10, pixtol=5):
        """ KD Tree algorithm to wavelength calibrate spectroscopic data.
        Currently, this is only designed for ThAr lamp spectra. See the
        'run_brute' function if you want to calibrate longslit spectra.

        Parameters
        ----------
        polygon : int
          Number of sides to the polygon used in pattern matching:
            polygon=3  -->  trigon (two anchor lines and one floating line)
            polygon=4  -->  tetragon (two anchor lines and two floating lines)
            polygon=5  -->  pentagon (two anchor lines and three floating lines)
            ...
        detsrch : int
          Number of consecutive detected lines used to generate a pattern. For
          example, if detsrch is 4, then for a trigon, the following patterns will
          be generated (assuming line #1 is the left anchor):
          1 2 3  (in this case line #3 is the right anchor)
          1 2 4  (in this case line #4 is the right anchor)
          1 3 4  (in this case line #4 is the right anchor)
        lstsrch : int
          Number of consecutive lines in the linelist used to generate a pattern.
          See example above for detsrch
        pixtol : float
          Tolerance used to find good patterns. An acceptable match if
          the closest distance to a pattern is < pixtol/npix, where npix
          is the number of pixels in the spectral direction. Ideally, this
          should depend on the pattern...

        Internals
        ---------

        detections : list of lists
        """

        # Load the linelist KD Tree
        lsttree, lindex = waveio.load_tree(polygon=polygon, numsearch=lstsrch)

        # Set the search error to be 5 pixels
        err = pixtol / self._npix

        self._all_patt_dict = {}
        self._all_final_fit = {}
        good_fit = np.zeros(self._nslit, dtype=np.bool)
        self._detections = {}
        for slit in range(self._nslit):
            if slit not in self._ok_mask:
                continue
            # Detect lines, and decide which tcent to use
            self._all_tcent, self._all_ecent, self._cut_tcent, self._icut =\
                wvutils.arc_lines_from_spec(self._spec[:, slit], min_nsig=self._min_nsig, nonlinear_counts = self._nonlinear_counts)
            self._all_tcent_weak, self._all_ecent_weak, self._cut_tcent_weak, self._icut_weak =\
                wvutils.arc_lines_from_spec(self._spec[:, slit], min_nsig=self._lowest_nsig, nonlinear_counts = self._nonlinear_counts)
            if self._all_tcent.size == 0:
                msgs.warn("No lines to identify in slit {0:d}!".format(slit+ 1))
                continue

            # Save the detections
            self._detections[str(slit)] = [self._all_tcent_weak.copy(), self._all_ecent_weak.copy()]

            use_tcentp, use_ecentp = self.get_use_tcent(1)
            use_tcentm, use_ecentm = self.get_use_tcent(-1)
            if use_tcentp.size < detsrch:
                if self._verbose:
                    msgs.info("Not enough lines to test this solution, will attempt another.")
                return None, None

            # Create a detlines KD Tree
            maxlinear = 0.5*self._npix
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

            msgs.info("Identifying wavelengths for each pattern")
            # First flatten the KD Tree query results so numba can handle the input array
            flatresp = [item for sublist in resultp for item in sublist]
            flatresm = [item for sublist in resultm for item in sublist]
            flatidxp = [ii for ii, sublist in enumerate(resultp) for item in sublist]
            flatidxm = [ii for ii, sublist in enumerate(resultm) for item in sublist]
            # Obtain the correlate and anti-correlate solutions
            psols = results_kdtree_nb(use_tcentp, self._wvdata, flatresp, flatidxp, indexp,
                                      lindex, indexp.shape[1], self._npix)
            msols = results_kdtree_nb(use_tcentm, self._wvdata, flatresm, flatidxm, indexm,
                                      lindex, indexm.shape[1], self._npix)

            msgs.info("Identifying the best solution")
            patt_dict, final_fit = self.solve_slit(slit, psols, msols, nselw=1, nseld=2)

            # Print preliminary report
            good_fit[slit] = self.report_prelim(slit, patt_dict, final_fit)

        # Now that all slits have been inspected, attempt to find a better
        # solution for orders that were not fit well, by estimating the
        # wavelength coverage of that slit
        #self.cross_match_order(good_fit)

        # With the updates to the fits of each slit, determine the final fit, and save the QA
        self.finalize_fit()

        # Print the final report of all lines
        self.report_final()
        return

    def cross_match(self, good_fit):
        """Cross-correlate the spectra across all slits to ID all of the lines.

        Parameters
        ----------
        good_fit : ndarray (bool)
          Indicates which slits are deemed to be a good fit (although, sometimes a bad fit can be
          labelled as a good fit). To remedy this, the true good fits are determined in this routine.
        """
        # Steps:
        # Check that all of the "good" slits are indeed good
        # For all of the bad slits, cross-correlate against all of the good slits to label each line
        # For all newly labeled lines, create a patt_dict of these labeled lines
        # Perform a final fit on these lines

        #self._debug = True
        # First, sort spectra according to increasing central wavelength
        ngd = good_fit.sum()
        idx_gd = np.zeros(ngd, dtype=np.int)
        wvc_gd = np.zeros(ngd, dtype=np.float)
        dsp_gd = np.zeros(ngd, dtype=np.float)
        wvc_gd_jfh = np.zeros(ngd, dtype=np.float)
        dsp_gd_jfh = np.zeros(ngd, dtype=np.float)
        xrng = np.arange(self._npix)
        cntr = 0
        for slit in range(self._nslit):
            if good_fit[slit]:
                idx_gd[cntr] = slit
                wvc_gd[cntr] = self._all_patt_dict[str(slit)]["bwv"]
                dsp_gd[cntr] = self._all_patt_dict[str(slit)]["bdisp"]
                # JFH stuff
                fitc = self._all_final_fit[str(slit)]['fitc']
                xfit = xrng/(self._npix - 1)
                fitfunc = self._all_final_fit[str(slit)]['function']
                fmin, fmax = 0.0, 1.0
                wave_soln = utils.func_val(fitc, xfit, fitfunc, minv=fmin, maxv=fmax)
                wvc_gd_jfh[cntr] = wave_soln[self._npix//2]
                dsp_gd_jfh[cntr]= np.median(wave_soln - np.roll(wave_soln,1))
                # JFH end of JFH stuff
                cntr += 1
        srt = np.argsort(wvc_gd)
        sort_idx = idx_gd[srt]
        sort_wvc = wvc_gd[srt]
        sort_dsp = dsp_gd[srt]
        sort_wvc_jfh = wvc_gd_jfh[srt]
        sort_dsp_jfh = dsp_gd_jfh[srt]

        # Cross correlate all good spectra with each other, in order of wavelength
        ncrco = ngd*(ngd-1)//2
        ccorr_val = np.zeros(ncrco)
        shift_val = np.zeros(ncrco)
        dwvc_val = np.zeros(ncrco)
        slit_ids = np.zeros((ncrco, 2), dtype=np.int)
        cntr = 0
        # JFH Consider adding something in here that takes advantage of the
        for gd in range(0, sort_idx.size-1):
            for gc in range(gd+1, sort_idx.size):
                #corr = scipy.signal.correlate(self._spec[:, sort_idx[gd]], self._spec[:, sort_idx[gc]], mode='same')
                #amax = np.argmax(corr)
                # dwvc_val[cntr] = (sort_wvc[gc]-sort_wvc[gd]) / (0.5*(sort_dsp[gc]+sort_dsp[gd])) - (amax - self._spec.shape[0] // 2)
                # JFH replaced with more robust xcorr
                shift_val[cntr], ccorr_val[cntr]= wvutils.xcorr_shift(self._spec[:, sort_idx[gd]],self._spec[:, sort_idx[gc]],
                                                                      smooth=5.0, percent_ceil=90.0)
                #dwvc_val[cntr] = (sort_wvc[gc]-sort_wvc[gd]) / (0.5*(sort_dsp[gc]+sort_dsp[gd])) - shift
                # JFH TESTING
                dwvc_val[cntr] = (sort_wvc_jfh[gc]-sort_wvc_jfh[gd]) / (0.5*(sort_dsp_jfh[gc]+sort_dsp_jfh[gd])) - shift_val[cntr]
                slit_ids[cntr, 0] = gd
                slit_ids[cntr, 1] = gc
                cntr += 1


        # TODO Replace this code below with code based on either sigma_clipped_stats or djs_reject
        # Identify the good slits as those for which the cross-correlation is consistent with the mad of all the slits.
        # Bad slits are then the outliers.
        sigrej = 3.0
        mad = 1.4826 * np.median(np.abs(dwvc_val))
        gdmsk = np.abs(dwvc_val) < sigrej * mad
        for ii in range(100):  # Limit to 100 iterations - this will likely never be reached...
            ogdmsk = gdmsk.copy()
            mad = 1.4826 * np.median(np.abs(dwvc_val[gdmsk]))
            gdmsk = np.abs(dwvc_val) < sigrej*mad
            if np.array_equal(gdmsk, ogdmsk):
                break

        if self._debug:
            # TODO Add something here indicating slit indices? Like plot the cc pairs next to the points?
            xplt = np.arange(dwvc_val.size)
            plt.plot(xplt[~gdmsk], dwvc_val[~gdmsk], 'rx',label ='bad slit')
            plt.plot(xplt[gdmsk], dwvc_val[gdmsk], 'bo',label = 'good slit')
            plt.hlines(0,xplt.min(),xplt.max(), color='black',linestyle='--')
            plt.xticks(xplt)
            plt.legend()
            plt.show()

        # Catalogue the good and bad slits.

        # ToDO Basically a slit needs to have a bad cross-correlation with every other slit in order
        # to be classified as a bad slit here. Is this the behavior we want?? Maybe we should be more
        # conservative and call a bad any slit which results in an outlier here?
        good_slits = np.sort(sort_idx[np.unique(slit_ids[gdmsk, :].flatten())])
        bad_slits = np.setdiff1d(np.arange(self._nslit), good_slits, assume_unique=True)
        nbad = bad_slits.size
        if nbad > 0:
            msgs.info('Working on {:d}'.format(nbad) + ' bad slits: {:}'.format(bad_slits + 1))

        # Get the sign (i.e. if pixels correlate/anticorrelate with wavelength)
        # and dispersion (A/pix). Assume these are the same for all slits

        # JFH Changed this to take the median which is more robust. Could even reject outliers
        disp_good = np.zeros(good_slits.size,dtype=float)
        sign_good = np.zeros(good_slits.size,dtype=int)
        wvc_good  = np.zeros(good_slits.size,dtype=float)
        for islit in range(good_slits.size):
            sign_good[islit] =  self._all_patt_dict[str(good_slits[islit])]['sign']
            # JFH stuff
            fitc = self._all_final_fit[str(good_slits[islit])]['fitc']
            xfit = xrng / (self._npix - 1)
            fitfunc = self._all_final_fit[str(good_slits[islit])]['function']
            fmin, fmax = 0.0, 1.0
            wave_soln = utils.func_val(fitc, xfit, fitfunc, minv=fmin, maxv=fmax)
            wvc_good[islit] = wave_soln[self._npix // 2]
            disp_good[islit] = np.median(wave_soln - np.roll(wave_soln, 1))


        disp_med = np.median(disp_good)
        sign = np.median(sign_good)
        #disp = self._all_patt_dict[str(good_slits[0])]['bdisp']
        #sign = self._all_patt_dict[str(good_slits[0])]['sign']

        # For all of the bad slits, estimate some line wavelengths
        new_bad_slits = np.array([], dtype=np.int)
        for bs in bad_slits:
            if bs not in self._ok_mask:
                continue
            if self._detections[str(bs)][0] is None:  # No detections at all; slit is hopeless
                msgs.warn("Slit {} has no arc line detections.  Likely this slit is junk!")
                self._bad_slits.append(bs)
                continue
            bsdet, _ = self.get_use_tcent(sign, arrerr=self._detections[str(bs)], weak=True)
            lindex = np.array([], dtype=np.int)
            dindex = np.array([], dtype=np.int)
            wcen = np.zeros(good_slits.size)
            disp = np.zeros(good_slits.size)
            shift_vec = np.zeros(good_slits.size)
            stretch_vec = np.zeros(good_slits.size)
            ccorr_vec = np.zeros(good_slits.size)
            for cntr, gs in enumerate(good_slits):
                msgs.info('Cross-correlating bad slit # {:d}'.format(bs + 1) + ' with good slit # {:d}'.format(gs + 1))
                # Match the peaks between the two spectra.
                # spec_gs_adj is the stretched spectrum
                success, shift_vec[cntr], stretch_vec[cntr], ccorr_vec[cntr], _, _ =  \
                    wvutils.xcorr_shift_stretch(self._spec[:, bs],self._spec[:, gs], debug = self._debug)
                if not success:
                    continue
                # JFH Put in a cut on the cross-correlation value here in this logic
                # so that we only consider slits that are sufficiently similar?

                # Estimate wcen and disp for this bad slit based on its shift/stretch relative to the good slit
                disp[cntr] = disp_good[cntr]/stretch_vec[cntr]
                wcen[cntr] = wvc_good[cntr] - shift_vec[cntr]*disp[cntr]

                # For each peak in the gs spectrum, identify the corresponding peaks in the bs spectrum. Do this by
                # transform these good slit line pixel locations into the (shifted and stretched) bs frame
                gsdet, _ = self.get_use_tcent(sign, arrerr=self._detections[str(gs)], weak=True)
                gsdet_ss = gsdet*stretch_vec[cntr] + shift_vec[cntr]
                if self._debug:
                    plt.figure(figsize=(14, 6))
                    xrng = np.arange(self._npix)
                    tampl_bs = np.interp(bsdet, xrng, self._spec[:, bs])
                    plt.plot(xrng, self._spec[:, bs], color='red', drawstyle='steps-mid', label='bad slit arc', linewidth=1.0, zorder= 10)
                    plt.plot(bsdet, tampl_bs, 'r.', markersize=10.0, label='bad slit lines', zorder= 10)
                    tampl_gs = np.interp(gsdet, xrng, self._spec[:, gs])
                    plt.plot(xrng, self._spec[:, gs], color='black', drawstyle='steps-mid', linestyle=':',
                             label='good slit arc', linewidth=0.5)
                    plt.plot(gsdet, tampl_gs, 'k+', markersize=8.0, label='good slit lines')
                    gdarc_ss = wvutils.shift_and_stretch(self._spec[:, gs], shift_vec[cntr], stretch_vec[cntr])
                    #tampl_ss = np.interp(gsdet_ss, xrng, gdarc_ss)
                    for iline in range(gsdet_ss.size):
                        plt.plot([gsdet[iline],gsdet_ss[iline]],[tampl_gs[iline], tampl_gs[iline]], color='cornflowerblue', linewidth = 1.0)
                    plt.plot(xrng, gdarc_ss, color='black', drawstyle='steps-mid', label='good slit arc shift/stretch', linewidth=1.0)
                    plt.plot(gsdet_ss, tampl_gs, 'k.', markersize=10.0, label='predicted good slit lines')
                    plt.title('Cross-correlation of bad slit # {:d}'.format(bs+1) + ' and good slit # {:d}'.format(gs+1) +
                              ': ccor = {:5.3f}'.format(ccorr_vec[cntr]) +
                              ', shift = {:6.1f}'.format(shift_vec[cntr]) +
                              ', stretch = {:5.4f}'.format(stretch_vec[cntr]) +
                              ', wv_cen = {:7.1f}'.format(wcen[cntr]) +
                              ', disp = {:5.3f}'.format(disp[cntr]))
                    plt.ylim(-5.0, 1.5*self._spec[:, bs].max())
                    plt.legend()
                    plt.show()

                # Calculate wavelengths for all of the gsdet detections
                fitc = self._all_final_fit[str(gs)]['fitc']
                xfit = gsdet/(self._npix - 1)
                fitfunc = self._all_final_fit[str(gs)]['function']
                fmin, fmax = 0.0, 1.0
                wvval = utils.func_val(fitc, xfit, fitfunc, minv=fmin, maxv=fmax)
                # Loop over the bad slit line pixel detections and find the nearest good slit line
                for dd in range(bsdet.size):
                    pdiff = np.abs(bsdet[dd]-gsdet_ss)
                    bstpx = np.argmin(pdiff)
                    # If a match is found within 2 pixels, consider this a successful match
                    if pdiff[bstpx] < 2.0:
                        # Using the good slit wavelength solution, search for the nearest line in the line list
                        bstwv = np.abs(self._wvdata - wvval[bstpx])
                        # This is probably not a good match
                        if bstwv[np.argmin(bstwv)] > 2.0*disp_med:
                            continue
                        lindex = np.append(lindex, np.argmin(bstwv)) # index in the line list self._wvdata
                        dindex = np.append(dindex, dd)               # index in the array of pixel detections bsdet
            # Finalize the best guess of each line
            # Initialise the patterns dictionary
            patt_dict = dict(acceptable=False, nmatch=0, ibest=-1, bwv=0., min_nsig=self._min_nsig,
                             mask=np.zeros(bsdet.size, dtype=np.bool))
            patt_dict['sign'] = sign
            patt_dict['bwv'] = np.median(wcen[wcen != 0.0])
            patt_dict['bdisp'] = np.median(disp[disp != 0.0])
            patterns.solve_triangles(bsdet, self._wvdata, dindex, lindex, patt_dict = patt_dict)

            if self._debug:
                tmp_list = vstack([self._line_lists, self._unknwns])
                qa.match_qa(self._spec[:, bs], bsdet, tmp_list,patt_dict['IDs'], patt_dict['scores'])

            # Use only the perfect IDs
            iperfect = np.array(patt_dict['scores']) != 'Perfect'
            patt_dict['mask'][iperfect] = False
            patt_dict['nmatch'] = np.sum(patt_dict['mask'])
            if patt_dict['nmatch'] < 3:
                patt_dict['acceptable'] = False

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
            if self._debug:
                xrng = np.arange(self._npix)
                xplt = np.linspace(0.0, 1.0, self._npix)
                yplt = utils.func_val(final_fit['fitc'], xplt, 'legendre', minv=0.0, maxv=1.0)
                plt.plot(final_fit['xfit'], final_fit['yfit'], 'bx')
                plt.plot(xplt, yplt, 'r-')
                plt.show()
                #pdb.set_trace()
        return new_bad_slits

    def cross_match_order(self, good_fit):
        """Using the solutions of all orders, identify the good solutions, and refit the bad ones!

        TODO: This function needs work... The first few lines of code successfully pick up the good orders,
        but we need a new routine that (based on an estimated central wavelength and dispersion) can successfully
        ID all of the lines.
        """

        # First determine the central wavelength and dispersion of every slit, using the known good solutions
        xplt = np.arange(self._nslit)
        yplt, dplt = np.zeros(self._nslit), np.zeros(self._nslit)
        imsk = np.ones(self._nslit, dtype=np.int)
        for slit in range(self._nslit):
            if good_fit[slit]:
                yplt[slit] = self._all_patt_dict[str(slit)]['bwv']
                dplt[slit] = self._all_patt_dict[str(slit)]['bdisp']
                imsk[slit] = 0

        mask, fit = utils.robust_polyfit(xplt, yplt, 2, function='polynomial', sigma=2,
                                         initialmask=imsk, forceimask=True)
        good_fit[mask == 1] = False
        wavemodel = utils.func_val(fit, xplt, 'polynomial')
        disp = np.median(dplt[good_fit])

        # TODO: maybe rethink the model at this point? Using the derived
        # central wavelength and dispersion identify liens in all orders?

        if self._debug:
            plt.subplot(211)
            plt.plot(xplt, wavemodel, 'r-')
            ww = np.where(mask==0)
            plt.plot(xplt[ww], yplt[ww], 'bx')
            ww = np.where(mask==1)
            plt.plot(xplt[ww], yplt[ww], 'rx')
            plt.subplot(212)
            plt.plot(xplt, dplt, 'bx')
            plt.show()
            pdb.set_trace()

        fact_nl = 1.2  # Non linear factor
        new_good_fit = np.zeros(self._nslit, dtype=np.bool)
        for slit in range(self._nslit):
            wmin = wavemodel[slit] - fact_nl*disp*self._npix/2
            wmax = wavemodel[slit] + fact_nl*disp*self._npix/2
            ww = np.where((self._wvdata > wmin) & (self._wvdata < wmax))
            wavedata = self._wvdata[ww]
            msgs.info('Brute force ID for slit {0:d}/{1:d}'.format(slit+1, self._nslit))
            best_patt_dict, best_final_fit =\
                self.run_brute_loop(slit, arrerr=self._detections[str(slit)], wavedata=wavedata)

            self._all_patt_dict[str(slit)] = copy.deepcopy(best_patt_dict)
            self._all_final_fit[str(slit)] = copy.deepcopy(best_final_fit)
            new_good_fit[slit] = self.report_prelim(slit, best_patt_dict, best_final_fit)
        return new_good_fit


        # Set some fitting parameters
        if self._n_final is None:
            order = 4
        else:
            order = self._n_final

        ofit = [5, 3, 1, 0]
        lnpc = len(ofit) - 1

        # Prepare the fitting coefficients
        xv = np.arange(self._npix)/(self._npix-1)
        ords = np.arange(self._nslit)
        xcen = xv[:, np.newaxis].repeat(self._nslit, axis=1)
        extrapord = ~good_fit
        maskord = np.where(extrapord)[0]

        coeffs = None
        waves = np.zeros(xcen.shape, dtype=np.float)
        for slit in range(self._nslit):
            if good_fit[slit]:
                func = self._all_final_fit[str(slit)]['function']
                fmin = self._all_final_fit[str(slit)]['fmin']
                fmax = self._all_final_fit[str(slit)]['fmax']
                fitc = self._all_final_fit[str(slit)]['fitc']
                if coeffs is None:
                    coeffs = np.zeros((fitc.size, self._nslit))
                coeffs[:, slit] = fitc.copy()
                waves[:, slit] = utils.func_val(fitc, xv, func, minv=fmin, maxv=fmax)

        msgs.info("Performing a PCA on the order wavelength solutions")
        pdb.set_trace()
        pca_wave, outpar = pca.basis(xcen, waves, coeffs, lnpc, ofit, x0in=ords, mask=maskord, skipx0=False, function=func)

        # Report the QA
        # TODO: fix setup passing
        setup = "BLAH"
        pca.pca_plot(setup, outpar, ofit, "wave_cross_match", pcadesc="Wavelength calibration PCA")


        # Extrapolate the remaining orders requested
        #extrap_wave, outpar = pca.extrapolate(outpar, ords)

        # Determine if pixels correlate and anti-correlate with wavelength
        signs = np.zeros(self._nslit, dtype=np.int)
        for slit in range(self._nslit):
            wvval = pca_wave[:, slit]
            if wvval[wvval.size//2] > wvval[wvval.size//2-1]:
                signs[slit] = 1
            else:
                signs[slit] = -1
        sign = 1
        if np.sum(signs) < 0:
            sign = -1

        new_bad_slits = np.array([], dtype=np.int)
        # Using the first guesses at the wavelength solution, identify lines
        for slit in range(self._nslit):
            # Get the detections
            dets, _ = self.get_use_tcent(sign, arrerr=self._detections[str(slit)], weak=True)
            lindex = np.array([], dtype=np.int)
            dindex = np.array([], dtype=np.int)
            # Calculate wavelengths for the gsdet detections
            wvval = pca_wave[:, slit]
            wvcen = wvval[wvval.size//2]
            disp = abs(wvval[wvval.size//2] - wvval[wvval.size//2-1])
            for dd in range(dets.size):
                pdiff = np.abs(dets[dd] - xv)
                bstpx = np.argmin(pdiff)
                bstwv = np.abs(self._wvdata - wvval[bstpx])
                if bstwv[np.argmin(bstwv)] > 10.0 * disp:
                    # This is probably not a good match
                    continue
                lindex = np.append(lindex, np.argmin(bstwv))
                dindex = np.append(dindex, dd)

            # Finalize the best guess of each line
            # Initialise the patterns dictionary
            patt_dict = dict(acceptable=False, nmatch=0, ibest=-1, bwv=0., min_nsig=self._min_nsig,
                             mask=np.zeros(dets.size, dtype=np.bool))
            patt_dict['sign'] = sign
            patt_dict['bwv'] = wvcen
            patt_dict['bdisp'] = disp

            patterns.solve_triangles(dets, self._wvdata, dindex, lindex, patt_dict)
            # Check if a solution was found
            if not patt_dict['acceptable']:
                new_bad_slits = np.append(new_bad_slits, slit)
                msgs.warn('---------------------------------------------------' + msgs.newline() +
                          'Cross-match report for slit {0:d}/{1:d}:'.format(slit + 1, self._nslit) + msgs.newline() +
                          '  Lines could not be identified! Will try cross matching iteratively' + msgs.newline() +
                          '---------------------------------------------------')
                continue
            final_fit = self.fit_slit(slit, patt_dict, tcent=dets)
            if final_fit is None:
                # This pattern wasn't good enough
                new_bad_slits = np.append(new_bad_slits, slit)
                msgs.warn('---------------------------------------------------' + msgs.newline() +
                          'Cross-match report for slit {0:d}/{1:d}:'.format(slit + 1, self._nslit) + msgs.newline() +
                          '  Fit was not good enough! Will try cross matching iteratively' + msgs.newline() +
                          '---------------------------------------------------')
                continue
            if final_fit['rms'] > self._rms_threshold:
                msgs.warn('---------------------------------------------------' + msgs.newline() +
                          'Cross-match report for slit {0:d}/{1:d}:'.format(slit + 1, self._nslit) + msgs.newline() +
                          '  Poor RMS ({0:.3f})! Will try cross matching iteratively'.format(final_fit['rms']) + msgs.newline() +
                          '---------------------------------------------------')
                # Store this result in new_bad_slits, so the iteration can be performed,
                # but make sure to store the result, as this might be the best possible.
                new_bad_slits = np.append(new_bad_slits, slit)
            self._all_patt_dict[str(slit)] = copy.deepcopy(patt_dict)
            self._all_final_fit[str(slit)] = copy.deepcopy(final_fit)
            if self._debug:
                xplt = np.linspace(0.0, 1.0, self._npix)
                yplt = utils.func_val(final_fit['fitc'], xplt, 'legendre', minv=0.0, maxv=1.0)
                plt.plot(final_fit['xfit'], final_fit['yfit'], 'bx')
                plt.plot(xplt, yplt, 'r-')
                plt.show()
                pdb.set_trace()

        # debugging
        if self._debug:
            # First determine the central wavelength and dispersion of every slit, using the known good solutions
            xplt = np.arange(self._nslit)
            yplt, dplt = np.zeros(self._nslit), np.zeros(self._nslit)
            imsk = np.ones(self._nslit, dtype=np.int)
            for slit in range(self._nslit):
                if good_fit[slit]:
                    yplt[slit] = self._all_patt_dict[str(slit)]['bwv']
                    dplt[slit] = self._all_patt_dict[str(slit)]['bdisp']
                    imsk[slit] = 0

            mask, fit = utils.robust_polyfit(xplt, yplt, 2, function='polynomial', sigma=2,
                                             initialmask=imsk, forceimask=True)

            ymodel = utils.func_val(fit, xplt, 'polynomial')
            plt.subplot(211)
            plt.plot(xplt, ymodel, 'r-')
            ww = np.where(mask==0)
            plt.plot(xplt[ww], yplt[ww], 'bx')
            ww = np.where(mask==1)
            plt.plot(xplt[ww], yplt[ww], 'rx')
            plt.subplot(212)
            plt.plot(xplt, dplt, 'bx')
            plt.show()
            pdb.set_trace()

        return new_bad_slits

    def get_use_tcent(self, corr, cut=True, arrerr=None, weak=False):
        """ Grab the lines to use
            Args:
                corr:  int
                  Set if pixels correlate with wavelength (corr==1) or anticorrelate (corr=-1)
                arrerr:
                weak: bool, optional
                   If True, return the weak lines
                cut: bool, optional
                   Cut on the lines according to significance

            Returns:
                arr: ndarray
                err: ndarray

            """
        # Decide which array to use
        if arrerr is None:
            if weak:
                if cut:
                    arr = self._all_tcent_weak.copy()[self._icut_weak]
                    err = self._all_ecent_weak.copy()[self._icut_weak]
                else:
                    debugger.set_trace()
            else:
                if cut:
                    arr = self._all_tcent.copy()[self._icut]
                    err = self._all_ecent.copy()[self._icut]
                else:
                    debugger.set_trace()
        else:
            arr, err = arrerr[0], arrerr[1]
        # Return the appropriate tcent
        if corr == 1:
            return arr, err
        else:
            return (self._npix - 1.0) - arr[::-1], err[::-1]

    def results_brute(self, poly=3, pix_tol=0.5, detsrch=5, lstsrch=5, wavedata=None, arrerr=None):
        """
        Need some docs here. I think this routine generates the patterns, either triangles are quadrangles.

        :param poly:     algorithms to use for pattern matching. Only triangles (3) and quadrangles (4) are supported
        :param pix_tol:  tolerance that is used to determine if a pattern match is successful (in units of pixels)
        :param detsrch:  Number of lines to search over for the detected lines
        :param lstsrch:  Number of lines to search over for the detected lines
        :param wavedata:
        :param arrerr:
        :return:
        """
        # Import the pattern matching algorithms
        if poly == 3:
            from pypeit.core.wavecal.patterns import triangles as generate_patterns
        elif poly == 4:
            from pypeit.core.wavecal.patterns import quadrangles as generate_patterns
        else:
            msgs.warn("Pattern matching is only available for trigons and tetragons.")
            return None, None

        if wavedata is None:
            wavedata = self._wvdata

        # Test if there are enough lines to generate a solution
        use_tcent, _ = self.get_use_tcent(1, arrerr=arrerr)
        if use_tcent.size < lstsrch or use_tcent.size < detsrch:
            if self._verbose:
                msgs.info("Not enough lines to test this solution, will attempt another.")
            return None, None

        if self._verbose:
            msgs.info("Begin pattern matching")
        # First run pattern recognition assuming pixels correlate with wavelength
        dindexp, lindexp, wvcenp, dispsp = generate_patterns(use_tcent, wavedata, self._npix,
                                                             detsrch, lstsrch, pix_tol)
        # Now run pattern recognition assuming pixels correlate with wavelength
        use_tcent, _ = self.get_use_tcent(-1, arrerr=arrerr)
        dindexm, lindexm, wvcenm, dispsm = generate_patterns(use_tcent, wavedata, self._npix,
                                                             detsrch, lstsrch, pix_tol)
        return (dindexp, lindexp, wvcenp, dispsp,), (dindexm, lindexm, wvcenm, dispsm,)

    def results_kdtree(self, use_tcent, res, dindex, lindex, ordfit=2):
        # Assign wavelengths to each pixel
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
                # dx = use_tcent[dindex[x, -1]] - use_tcent[dindex[x, 0]]
                # dp = self._wvdata[lindex[res[x][y], -1]] - self._wvdata[lindex[res[x][y], 0]]
                # try:
                #     null, cgrad = utils.robust_polyfit(use_tcent[dindex[x, :]], self._wvdata[lindex[res[x][y], :]],
                #                                        1, sigma=2.0, verbose=False)
                #     wvdisp[cnt] = cgrad[1]
                # except:
                #     wvdisp[cnt] = (dp / dx)
                #
                coeff = np.polyfit(use_tcent[dindex[x, :]], self._wvdata[lindex[res[x][y]]], ordfit)
                wvcent[cnt] = np.polyval(coeff, self._npix / 2.0)
                wvdisp[cnt] = abs(np.polyval(coeff, (self._npix+1) / 2.0) - wvcent[cnt])
                dind[cnt, :] = dindex[x, :]
                lind[cnt, :] = lindex[res[x][y], :]
                cnt += 1
        return dind, lind, wvcent, wvdisp

    def solve_slit(self, slit, psols, msols, nstore=1, nselw=3, nseld=3):
        """
        Need some docs here. I think this routine creates a 2d histogram of the patterns and searches for the most
        represented wave_cen and log10(disp). Then it attempts to fit each value determined (default of 1) to
        try to figure out if it is a reasonable fit.

        :param slit:
        :param psols:
        :param msols:
        :param nstore: Number of pattern matches to store and fit
        :param nselw:  All solutions around the best central wavelength solution within +- nselw are selected to be fit
        :param nseld:  All solutions around the best log10(dispersion) solution within +- nseld are selected to be fit
        :return:  patt_dict, final_dict
                patt_dict = ??
                final_dict = ??
        """

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
        sm_histimg = gaussian_filter(histimg, [30, 15])

        #histpeaks = patterns.detect_2Dpeaks(np.abs(sm_histimg))
        histpeaks = patterns.detect_2Dpeaks(np.abs(histimg))

        # Find the indices of the nstore largest peaks
        bidx = np.unravel_index(np.argpartition(np.abs(histpeaks*histimg), -nstore, axis=None)[-nstore:], histimg.shape)

        # Get the peak value of central wavelength and dispersion
        allwcen = self._binw[bidx[0]]
        alldisp = self._bind[bidx[1]]
        allhnum = np.abs(histimg[bidx])

        #debug = False
        if self._debug:# or slit==2:
            this_hist = histimg
            plt.clf()
            rect_image = [0.12, 0.05, 0.85, 0.9]
            fx = plt.figure(1, figsize=(12, 8))
            ax_image = fx.add_axes(rect_image)
            extent = [self._binw[0], self._binw[-1], self._bind[0], self._bind[-1]]
            # plt.subplot(221)
            # plt.imshow((np.abs(histimg[:, ::-1].T)), extent=extent, aspect='auto')
            # plt.subplot(222)
            # plt.imshow((np.abs(sm_histimg[:, ::-1].T)), extent=extent, aspect='auto')
            # plt.subplot(223)
            # plt.imshow((np.abs(histimgp[:, ::-1].T)), extent=extent, aspect='auto')
            # plt.subplot(224)
            # plt.imshow((np.abs(histimgm[:, ::-1].T)), extent=extent, aspect='auto')
            #plt.imshow((np.abs(sm_histimg[:, ::-1].T)), extent=extent, aspect='auto')
            cimg = ax_image.imshow(this_hist.T, extent=extent, aspect='auto',vmin=-2.0,vmax=5.0,
                       interpolation='nearest',origin='lower',cmap='Set1')
            nm = histimg.max() - histimg.min()
            ticks = np.arange(this_hist.min(),this_hist.max() + 1,1)
            cbar = fx.colorbar(cimg, ax=ax_image,ticks = ticks,drawedges = True, extend ='both',
                               spacing = 'proporational',orientation ='horizontal')
            cbar.set_ticklabels(ticks)
            cbar.set_label('# of Occurences')
            ax_image.set_xlabel('Central Wavelength (Angstroms)')
            ax_image.set_ylabel('log10(Dispersion/(Ang/pixel))')
            delta_wv = np.median(self._binw - np.roll(self._binw,1))
            delta_disp = np.median(self._bind - np.roll(self._bind,1))
            label = 'Maximum: (lam, log10(disp)) = ({:8.2f}'.format(allwcen[0]) + ',{:5.3f})'.format(alldisp[0])
            ax_image.plot(allwcen + delta_wv/2.0, alldisp + delta_disp/2.0, color='red', marker='+', markersize=10.0, fillstyle='none',
                     linestyle='None', zorder = 10,label=label)
            ax_image.legend()
            #plt.plot(self._binw[bidx[0]], self._bind[bidx[1]], 'r+')
#            ax_image.axvline(allwcen, color='r', linestyle='--')
#            ax_image.axhline(alldisp, color='r', linestyle='--')
            plt.show()
            #if False:
            #    plt.clf()
            #    plt.imshow(histimgp[:, ::-1].T, extent=extent, aspect='auto',vmin=0.0,vmax=1.0)
            #    plt.show()


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
            elif tfinal_dict['rms'] < self._rms_threshold:
                # Has a better fit been identified (i.e. more lines ID)?
                if len(tfinal_dict['xfit']) > len(final_dict['xfit']):
                    patt_dict, final_dict = copy.deepcopy(tpatt_dict), copy.deepcopy(tfinal_dict)
        return patt_dict, final_dict

    def solve_patterns(self, bestlist):

        # Obtain a full list of indices that are consistent with the maximum value
        wcen, dcen, sign, dindex, lindex = bestlist[0], bestlist[1], bestlist[3], bestlist[4], bestlist[5]

        # Find the favoured sign and only use those values
        use_tcent, _ = self.get_use_tcent(sign)
        if sign == +1:
            signtxt = "correlate"
        else:
            signtxt = "anticorrelate"

        # Initialise the patterns dictionary
        patt_dict = dict(acceptable=False, nmatch=0, ibest=-1, bwv=0., min_nsig=self._min_nsig,
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
                      '  Best wave/disp                = {:g}'.format(patt_dict['bwv']/patt_dict['bdisp']) + msgs.newline() +
                      '---------------------------------------------------')
        return patt_dict

    def fit_slit(self, slit, patt_dict, outroot=None, slittxt="Slit", tcent=None, ecent=None):
        """ Perform a fit to the wavelength solution

        Parameters
        ----------
        slit : int
          slit number
        patt_dict : dict
          dictionary of patterns
        outroot : str
          root directory to save QA
        slittxt : str
          Label used for QA
        tcent : ndarray
          List of the detections in this slit...
        ecent : ndarray
          ... and their corresponding errors

        Returns
        -------
        final_fit : dict
          A dictionary containing all of the information about the fit
        """
        # Perform final fit to the line IDs
        if self._thar:
            NIST_lines = (self._line_lists['NIST'] > 0) & (np.char.find(self._line_lists['Source'].data, 'MURPHY') >= 0)
        else:
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
            tcent, ecent = self.get_use_tcent(patt_dict['sign'], weak=True)
            #weights = 1.0/ecent
            weights = np.ones(tcent.size)
        else:
            if ecent is None:
                weights = np.ones(tcent.size)
            else:
                #weights = 1.0/ecent
                weights = np.ones(tcent.size)
        # Fit
        try:
            final_fit = fitting.iterative_fitting(self._spec[:, slit], tcent, ifit,
                                                  np.array(patt_dict['IDs'])[ifit], self._line_lists[NIST_lines],
                                                  patt_dict['bdisp'], weights=weights,
                                                  match_toler=self._match_toler, func=self._func, n_first=self._n_first,
                                                  sigrej_first=self._sigrej_first,
                                                  n_final=self._n_final, sigrej_final=self._sigrej_final,
                                                  plot_fil = plot_fil, verbose = self._verbose)
        except TypeError:
            # A poor fitting result, this can be ignored.
            return None

        if plot_fil is not None:
            print("Wrote: {:s}".format(plot_fil))

        # Return
        return final_fit

    def finalize_fit(self):
        """ Once the best IDs have been found for each slit, perform a final fit to all slits and save the results
        """
        for slit in range(self._nslit):
            if slit not in self._ok_mask:
                continue
            if self._all_patt_dict[str(slit)] is None:
                continue
            # Save the QA for the best solution
            slittxt = '_Slit{0:03d}'.format(slit+1)
            use_tcent, use_ecent = self.get_use_tcent(self._all_patt_dict[str(slit)]['sign'],
                                                      arrerr=self._detections[str(slit)], weak=True)
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
                            outfile=self._outroot + slittxt + '.pdf')
                msgs.info("Wrote: {:s}".format(self._outroot + slittxt + '.pdf'))
            # Perform the final fit for the best solution
            best_final_fit = self.fit_slit(slit, self._all_patt_dict[str(slit)], tcent=use_tcent, ecent=use_ecent,
                                           outroot=self._outroot, slittxt=slittxt)
            self._all_final_fit[str(slit)] = copy.deepcopy(best_final_fit)

    def report_prelim(self, slit, best_patt_dict, best_final_fit):

        good_fit = False
        # Report on the best preliminary result
        if best_final_fit is None:
            msgs.warn('---------------------------------------------------' + msgs.newline() +
                      'Preliminary report for slit {0:d}/{1:d}:'.format(slit + 1, self._nslit) + msgs.newline() +
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
            good_fit = True
            if best_patt_dict['sign'] == +1:
                signtxt = 'correlate'
            else:
                signtxt = 'anitcorrelate'
            # Report
            msgs.info('---------------------------------------------------' + msgs.newline() +
                      'Preliminary report for slit {0:d}/{1:d}:'.format(slit + 1, self._nslit) + msgs.newline() +
                      '  Pixels {:s} with wavelength'.format(signtxt) + msgs.newline() +
                      '  Number of lines recovered    = {:d}'.format(self._all_tcent.size) + msgs.newline() +
                      '  Number of lines analyzed     = {:d}'.format(len(best_final_fit['xfit'])) + msgs.newline() +
                      '  Number of pattern matches    = {:d}'.format(best_patt_dict['nmatch']) + msgs.newline() +
                      '  Best central wavelength      = {:g}A'.format(best_patt_dict['bwv']) + msgs.newline() +
                      '  Best dispersion              = {:g}A/pix'.format(best_patt_dict['bdisp']) + msgs.newline() +
                      '  Best wave/disp               = {:g}'.format(best_patt_dict['bwv']/best_patt_dict['bdisp']) + msgs.newline() +
                      '  Final RMS of fit             = {:g}'.format(best_final_fit['rms']) + msgs.newline() +
                      '---------------------------------------------------')
            self._all_patt_dict[str(slit)] = copy.deepcopy(best_patt_dict)
            self._all_final_fit[str(slit)] = copy.deepcopy(best_final_fit)
        return good_fit

    def report_final(self):
        """Print out the final report of the wavelength calibration"""
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
            msgs.info(msgs.newline() +
                      '---------------------------------------------------' + msgs.newline() +
                      'Final report for slit {0:d}/{1:d}:'.format(slit + 1, self._nslit) + msgs.newline() +
                      '  Pixels {:s} with wavelength'.format(signtxt) + msgs.newline() +
                      '  Number of lines recovered    = {:d}'.format(self._detections[st][0].size) + msgs.newline() +
                      '  Number of lines analyzed     = {:d}'.format(len(self._all_final_fit[st]['xfit'])) + msgs.newline() +
                      '  Central wavelength           = {:g}A'.format(centwave) + msgs.newline() +
                      '  Central dispersion           = {:g}A/pix'.format(centdisp) + msgs.newline() +
                      '  Central wave/disp             = {:g}'.format(centwave/centdisp) + msgs.newline() +
                      '  Final RMS of fit             = {:g}'.format(self._all_final_fit[st]['rms']))
        return


@nb.jit(nopython=True, cache=True)
def results_kdtree_nb(use_tcent, wvdata, res, residx, dindex, lindex, nindx, npix, ordfit=1):
    """ A numba speedup of the results_kdtree function in the General class (see above).
    For all of the acceptable pattern matches, estimate the central wavelength and dispersion,
    and record the index in the linelist and the corresponding indices of the detected lines.

    Parameters
    ----------
    use_tcent : ndarray
      detected lines
    wvdata : ndarray
      the linelist
    res : list
      A flattened list of the results of the ball tree query from the KDTree (this contains all acceptable matches)
      This needs to be a flattened list for numba
    residx : list
      This contains the original indices of the unflattened (i.e. nested) 'res' list
    dindex : ndarray
      Indices of the lines in the detected lines for all patterns
    lindex : ndarray
      Indices of the lines in the linelist for all patterns
    nindx : int
      Number of acceptable pattens
    npix : int
      Number of pixels in the spectral direction
    ordfit : int
      Order of the polynomial used to fit the pixel/wavelength IDs

    Returns
    -------
    dind : ndarray
      Indices of the lines in the detected lines that were used for each acceptable pattern
    lind : linelist index of patterns
      Indices of the lines in the linelist that were used for each corresponding pattern
    wvcent : ndarray
      Central wavelength of each pattern
    wvdisp : ndarray
      Central dispersion of each pattern
    """
    # Assign wavelengths to each pixel
    ncols = len(res)
    wvdisp = np.zeros(ncols, dtype=nb.types.float64)
    wvcent = np.zeros(ncols, dtype=nb.types.float64)
    dind = np.zeros((ncols, nindx), dtype=nb.types.uint64)
    lind = np.zeros((ncols, nindx), dtype=nb.types.uint64)
    Xmat = np.ones((nindx, ordfit+1), dtype=nb.types.float64)
    for x in range(ncols):
        for ii in range(ordfit, -1, -1):
            Xmat[:, ii] = np.power(use_tcent[dindex[residx[x], :]], ordfit-ii)
        coeff = np.linalg.lstsq(Xmat, wvdata[lindex[res[x]]])[0]
        sumv = 0.0
        for ii in range(ordfit, -1, -1):
            wvcent[x] += coeff[ordfit-ii]*((npix/2.0)**ii)
            sumv += coeff[ordfit-ii]*(((npix+1)/2.0)**ii)
        wvdisp[x] = abs(sumv - wvcent[x])
        dind[x, :] = dindex[residx[x], :]
        lind[x, :] = lindex[res[x], :]
    return dind, lind, wvcent, wvdisp
