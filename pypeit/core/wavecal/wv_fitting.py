""" Module for finding patterns in arc line spectra
"""
import numpy as np
import inspect

from pypeit.core.wavecal import autoid
from pypeit.core.wavecal import defs
from pypeit.core import fitting
from pypeit import msgs

from pypeit import datamodel

from IPython import embed


class WaveFit(datamodel.DataContainer):
    """
    DataContainer for the output from BuildWaveCalib

    All of the items in the datamodel are required for instantiation,
      although they can be None (but shouldn't be)

    """
    minimum_version = '1.0.0'
    version = '1.0.0'

    # I/O
    output_to_disk = None  # This writes all items that are not None

    datamodel = {
        'pypeitfit': dict(otype=fitting.PypeItFit, desc='Fit to 1D wavelength solutions'),
        'pixel_fit': dict(otype=np.ndarray, atype=np.floating, desc='Pixel values of arc lines'),
        'wave_fit': dict(otype=np.ndarray, atype=np.floating, desc='Wavelength IDs assigned'),
        'xnorm': dict(otype=float, desc='Normalization for fit'),
        'ion_bits': dict(otype=np.ndarray, atype=np.integer, desc='Ion bit values for the Ion names'),
        'cen_wave': dict(otype=float, desc='Central wavelength'),
        'cen_disp': dict(otype=float, desc='Approximate wavelength dispersion'),
        'spec': dict(otype=np.ndarray, atype=np.floating, desc='Arc spectrum'),
        'wave_soln': dict(otype=np.ndarray, atype=np.floating, desc='Evaluated wavelengths at pixel_fit'),
        'sigrej': dict(otype=float, desc='Final sigma rejection applied'),
        'shift': dict(otype=float, desc='Shift applied'),
        'tcent': dict(otype=np.ndarray, atype=np.floating, desc='Pixel centroids of all arc lines found'),
        'rms': dict(otype=float, desc='RMS of the solution'),
    }

    bitmask = defs.LinesBitMask()

    def __init__(self, pypeitfit=None, pixel_fit=None, wave_fit=None, ion_bits=None,
                 cen_wave=None, cen_disp=None, spec=None, wave_soln=None,
                 sigrej=None, shift=None, tcent=None, rms=None, xnorm=None):
        # Parse
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=d)

    def _init_internals(self):
        # Needs to be here so we can set it in WaveCalib
        self.hdu_prefix = None

    def _bundle(self, **kwargs):
        """
        Over-ride DataContainer._bundle() to deal with DETECTOR

        Args:
            kwargs:
                Passed to DataContainer._bundle()

        Returns:
            list:

        """
        _d = super(WaveFit, self)._bundle(ext='WAVEFIT', **kwargs)
        # Pull the fit out
        if _d[0]['WAVEFIT']['pypeitfit'] is not None:
            _d.append(dict(PYPEITFIT=_d[0]['WAVEFIT'].pop('pypeitfit')))
        # Return
        return _d


    # TODO -- This snippet shows up in 3 places now.  Can we genrealize somehow KW??
    def to_hdu(self, hdr=None, add_primary=False, primary_hdr=None,
               limit_hdus=None, force_to_bintbl=True):
        """
        Over-ride :func:`pypeit.datamodel.DataContainer.to_hdu` to force to
        a BinTableHDU

        See that func for Args and Returns
        """
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = dict([(k,values[k]) for k in args[1:]])
        # Force
        _d['force_to_bintbl'] = True
        # Do it
        return super(WaveFit, self).to_hdu(**_d)

    @property
    def ions(self):
        """

        Returns:
            `numpy.ndarray`_:  Array of the ion label for each line as recorded in ion_bits

        """
        ionlist = []
        for ionbit in self.ion_bits:
            ionlist += self.bitmask.flagged_bits(ionbit)
        # Return
        return np.asarray(ionlist)



def fit_slit(spec, patt_dict, tcent, line_lists, vel_tol = 1.0, outroot=None, slittxt="Slit", thar=False,match_toler=3.0,
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
    line_lists: astropy Table
      Table containing the line list
    Optional Parameters
    -------------------
    vel_tol: float, default = 1.0
      Tolerance in km/s for matching lines in the IDs to lines in the NIST database. The default is 1.0 km/s
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
        if (np.min(np.abs(line_lists['wave'][NIST_lines] - idwv)))/idwv*3.0e5 > vel_tol:
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
        msgs.warn('Fit failed for this slit')
        return None

    if plot_fil is not None:
        print("Wrote: {:s}".format(plot_fil))

    # Return
    return final_fit


def iterative_fitting(spec, tcent, ifit, IDs, llist, disp,
                      match_toler = 2.0, func = 'legendre', n_first=2, sigrej_first=2.0,
                      n_final=4, sigrej_final=3.0, input_only=False,
                      weights=None, plot_fil=None, verbose=False):

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
    input_only: bool
      If True, the routine will only perform a robust polyfit to the input IDs.
      If False, the routine will fit the input IDs, and then include additional
      lines in the linelist that are a satisfactory fit.
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

    #TODO JFH add error checking here to ensure that IDs and ifit have the same size!

    if weights is None:
        weights = np.ones(tcent.size)

    nspec = spec.size
    xnspecmin1 = float(nspec-1)
    # Setup for fitting
    sv_ifit = list(ifit)  # Keep the originals
    all_ids = -999.*np.ones(len(tcent))
    all_idsion = np.array(['UNKNWN']*len(tcent))
    all_ids[ifit] = IDs

    # Fit
    n_order = n_first
    flg_continue = True
    flg_penultimate = False
    fmin, fmax = 0.0, 1.0
    # Note the number of parameters is actually n_order and not n_order+1
    while flg_continue:
        if flg_penultimate:
            flg_continue = False
        # Fit with rejection
        xfit, yfit, wfit = tcent[ifit], all_ids[ifit], weights[ifit]
        #mask, fit = utils.robust_polyfit(xfit/xnspecmin1, yfit, n_order, function=func, sigma=sigrej_first,
        #                                 minx=fmin, maxx=fmax, verbose=verbose, weights=wfit)
        pypeitFit = fitting.robust_fit(xfit/xnspecmin1, yfit, n_order, function=func,
                                       lower=sigrej_first, upper=sigrej_first,
                                       minx=fmin, maxx=fmax, weights=wfit)

        rms_ang = pypeitFit.calc_fit_rms()#xfit[pypeitFit.gpm == 0]/xnspecmin1,
                                     #yfit[pypeitFit.gpm == 0],
                                     #weights=wfit[pypeitFit.gpm == 0])
        rms_pix = rms_ang/disp
        if verbose:
            msgs.info('n_order = {:d}'.format(n_order) + ': RMS = {:g}'.format(rms_pix))

        # Reject but keep originals (until final fit)
        ifit = list(ifit[pypeitFit.gpm == 0]) + sv_ifit
        if not input_only:
            # Find new points from the linelist (should we allow removal of the originals?)
            twave = pypeitFit.val(tcent/xnspecmin1)#, func, minx=fmin, maxx=fmax)
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
        # Increment order?
        if n_order < n_final:
            n_order += 1
        else:
            flg_penultimate = True

    # Final fit (originals can now be rejected)
    #fmin, fmax = 0., 1.
    #xfit, yfit, wfit = tcent[ifit]/(nspec-1), all_ids[ifit], weights[ifit]
    xfit, yfit, wfit = tcent[ifit], all_ids[ifit], weights[ifit]
    #mask, fit = utils.robust_polyfit(xfit/xnspecmin1, yfit, n_order, function=func, sigma=sigrej_final,
    pypeitFit = fitting.robust_fit(xfit/xnspecmin1, yfit, n_order, function=func,
                                   lower=sigrej_final, upper=sigrej_final,
                                   minx=fmin, maxx=fmax, weights=wfit)#, debug=True)
    irej = np.where(pypeitFit.gpm == 1)[0]
    if len(irej) > 0:
        xrej = xfit[irej]
        yrej = yfit[irej]
        if verbose:
            for kk, imask in enumerate(irej):
                wave = pypeitFit.val(xrej[kk]/xnspecmin1)#, func, minx=fmin, maxx=fmax)
                msgs.info('Rejecting arc line {:g}; {:g}'.format(yfit[imask], wave))
    else:
        xrej = []
        yrej = []

    #xfit = xfit[mask == 0]
    #yfit = yfit[mask == 0]
    #wfit = wfit[mask == 0]
    ions = all_idsion[ifit]
#    ions = all_idsion[ifit][mask == 0]
    # Final RMS
    rms_ang = pypeitFit.calc_fit_rms()#xfit[pypeitFit.gpm == 0] / xnspecmin1,
                                     #yfit[pypeitFit.gpm == 0],
                                     #weights=wfit[pypeitFit.gpm == 0])
#    rms_ang = utils.calc_fit_rms(xfit, yfit, fit, func,
#                                 minx=fmin, maxx=fmax, weights=wfit)
    rms_pix = rms_ang/disp

    # Pack up fit
    spec_vec = np.arange(nspec)
    #wave_soln = utils.func_val(fit,spec_vec/xnspecmin1, func, minx=fmin, maxx=fmax)
    wave_soln = pypeitFit.val(spec_vec/xnspecmin1)#, func, minx=fmin, maxx=fmax)
    #cen_wave = utils.func_val(fit, float(nspec)/2/xnspecmin1, func, minx=fmin, maxx=fmax)
    cen_wave = pypeitFit.val(float(nspec)/2/xnspecmin1)#, func, minx=fmin, maxx=fmax)
    #cen_wave_min1 = utils.func_val(fit, (float(nspec)/2 - 1.0)/xnspecmin1, func, minx=fmin, maxx=fmax)
    cen_wave_min1 = pypeitFit.val((float(nspec)/2 - 1.0)/xnspecmin1)#, func, minx=fmin, maxx=fmax)
    cen_disp = cen_wave - cen_wave_min1

    # Ions bit
    ion_bits = np.zeros(len(ions), dtype=WaveFit.bitmask.minimum_dtype())
    for kk,ion in enumerate(ions):
        ion_bits[kk] = WaveFit.bitmask.turn_on(ion_bits[kk], ion)

    # DataContainer time
    final_fit = WaveFit(pypeitfit=pypeitFit, pixel_fit=xfit, wave_fit=yfit, #weights=wfit,
                        ion_bits=ion_bits, xnorm=xnspecmin1, #nspec=nspec,
                        cen_wave=cen_wave, cen_disp=cen_disp,
                       #xrej=xrej, yrej=yrej, mask=(mask == 0),
                        spec=spec, wave_soln = wave_soln, sigrej=sigrej_final,
                        shift=0., tcent=tcent, rms=rms_pix)

    # If set to True, this will output a file that can then be included in the tests
    #saveit = False
    #if saveit:
        #from linetools import utils as ltu
        #jdict = ltu.jsonify(final_fit)
        #if plot_fil is None:
            #outname = "temp"
            #print("You should have set the plot_fil directory to save wavelength fits... using 'temp' as a filename")
        #else:
            #outname = plot_fil
        #ltu.savejson(outname + '.json', jdict, easy_to_read=True, overwrite=True)
        #print(" Wrote: {:s}".format(outname + '.json'))

    # QA
    if plot_fil is not None:
        autoid.arc_fit_qa(final_fit, plot_fil)
    # Return
    return final_fit
