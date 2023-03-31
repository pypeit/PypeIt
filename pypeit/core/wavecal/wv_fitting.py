""" Module for finding patterns in arc line spectra

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import numpy as np
import inspect

from astropy.io import fits

from pypeit.core.wavecal import autoid
from pypeit.core.wavecal import defs
from pypeit.core import fitting
from pypeit import msgs

from pypeit import datamodel

from IPython import embed


class WaveFit(datamodel.DataContainer):
    """
    DataContainer for the output from BuildWaveCalib

    All of the items in the datamodel are required for instantiation, although
    they can be None (but shouldn't be).

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_wavefit.rst

    When written to an output-file HDU, all `numpy.ndarray`_ elements are
    bundled into an `astropy.io.fits.BinTableHDU`_, the ``pypeitfit`` attribute
    is written to a separate extension (see
    :class:`~pypeit.core.fitting.PypeItFit`), and the other elements are written
    as header keywords.  Any datamodel elements that are None are *not* included
    in the output.  The two HDU extensions are given names according to their
    spatial ID; see :func:`hduext_prefix_from_spatid`.

    """
    version = '1.1.0'

    datamodel = {'spat_id': dict(otype=(int,np.integer), descr='Spatial position of slit/order for this fit. Required for I/O'),
                 'pypeitfit': dict(otype=fitting.PypeItFit,
                                   descr='Fit to 1D wavelength solutions'),
                 'pixel_fit': dict(otype=np.ndarray, atype=np.floating,
                                   descr='Pixel values of arc lines'),
                 'wave_fit': dict(otype=np.ndarray, atype=np.floating,
                                  descr='Wavelength IDs assigned'),
                 'xnorm': dict(otype=float, descr='Normalization for fit'),
                 'fwhm': dict(otype=float, descr='Estimate FWHM of arc lines in binned pixels of the input arc frame'),
                 'ion_bits': dict(otype=np.ndarray, atype=np.integer,
                                  descr='Ion bit values for the Ion names'),
                 'cen_wave': dict(otype=float, descr='Central wavelength'),
                 'cen_disp': dict(otype=float, descr='Approximate wavelength dispersion'),
                 'spec': dict(otype=np.ndarray, atype=np.floating, descr='Arc spectrum'),
                 'wave_soln': dict(otype=np.ndarray, atype=np.floating,
                                   descr='Evaluated wavelengths at pixel_fit'),
                 'sigrej': dict(otype=float, descr='Final sigma rejection applied'),
                 'shift': dict(otype=float, descr='Shift applied'),
                 'tcent': dict(otype=np.ndarray, atype=np.floating,
                               descr='Pixel centroids of all arc lines found'),
                 'rms': dict(otype=float, descr='RMS of the solution')}

    bitmask = defs.LinesBitMask()

    @staticmethod
    def hduext_prefix_from_spatid(spat_id):
        """ Naming for HDU extensions"""
        return 'SPAT_ID-{}_'.format(spat_id)

    def __init__(self, spat_id, pypeitfit=None, pixel_fit=None, wave_fit=None, ion_bits=None,
                 cen_wave=None, cen_disp=None, spec=None, wave_soln=None,
                 sigrej=None, shift=None, tcent=None, rms=None, xnorm=None,
                 fwhm=None):
        # Parse
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=d)

    def _bundle(self, **kwargs):
        """
        Over-ride DataContainer._bundle() to deal with PYPEITFIT

        Args:
            kwargs:
                Passed to DataContainer._bundle()

        Returns:
            list:

        """
        # Extension prefix (for being unique with slits)
        hdue_pref = self.hduext_prefix_from_spatid(self.spat_id)
        # Without PypeItFit
        _d = super(WaveFit, self)._bundle(
            ext=hdue_pref+'WAVEFIT', **kwargs)
        # Deal with PypeItFit
        if _d[0][hdue_pref+'WAVEFIT']['pypeitfit'] is not None:
            _d.append({hdue_pref+'PYPEITFIT': _d[0][hdue_pref + 'WAVEFIT'].pop('pypeitfit')})
        # Return
        return _d

    def to_hdu(self, **kwargs):
        """
        Over-ride base-class function to force ``force_to_bintbl`` to always be
        true.

        See base-class :class:`~pypeit.datamodel.DataContainer.to_hdu` for arguments.

        Args:
            kwargs:
                Passed directly to the base-class
                :class:`~pypeit.datamodel.DataContainer.to_hdu`.  If 
                ``force_to_bintbl`` is present, it is forced to be True.

        Returns:
            :obj:`list`, `astropy.io.fits.HDUList`_: A list of HDUs,
            where the type depends on the value of ``add_primary``.
        """
        if 'force_to_bintbl' in kwargs:
            if not kwargs['force_to_bintbl']:
                msgs.warn(f'{self.__class__.__name__} objects must always use '
                          'force_to_bintbl = True!')
            kwargs.pop('force_to_bintbl')
        return super().to_hdu(force_to_bintbl=True, **kwargs)

    @classmethod
    def from_hdu(cls, hdu, **kwargs):
        """
        Parse the data from the provided HDU.

        See the base-class :func:`~pypeit.datamodel.DataContainer.from_hdu` for
        the argument descriptions.

        Args:
            kwargs:
                Passed directly to the base-class
                :class:`~pypeit.datamodel.DataContainer.from_hdu`.  Should not
                include ``hdu_prefix`` because this is set directly by the
                class, read from the SPAT_ID card in a relevant header.

        Returns:
            :class:`WaveFit`: Object instantiated from data in the provided HDU.
        """
        if 'hdu_prefix' in kwargs:
            kwargs.pop('hdu_prefix')
        # Set hdu_prefix
        spat_id = hdu[1].header['SPAT_ID'] if isinstance(hdu, fits.HDUList)\
                    else hdu.header['SPAT_ID']
        hdu_prefix = cls.hduext_prefix_from_spatid(spat_id)
        # Run the default parser to get the data
        return super(WaveFit, cls).from_hdu(hdu, hdu_prefix=hdu_prefix, **kwargs)

    @property
    def ions(self):
        """
        Returns an array of ion labels

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
    imsk = np.ones(len(ifit), dtype=bool)
    for kk, idwv in enumerate(np.array(patt_dict['IDs'])[ifit]):
        if (np.min(np.abs(line_lists['wave'][NIST_lines] - idwv)))/idwv*3.0e5 > vel_tol:
            imsk[kk] = False
    ifit = ifit[imsk]
    # Fit
    final_fit = iterative_fitting(spec, tcent, ifit, np.array(patt_dict['IDs'])[ifit], line_lists[NIST_lines],
                                  patt_dict['bdisp'],match_toler=match_toler, func=func, n_first=n_first,
                                  sigrej_first=sigrej_first,n_final=n_final, sigrej_final=sigrej_final,
                                  plot_fil=plot_fil, verbose=verbose)
    if plot_fil is not None and final_fit is not None:
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
    final_fit: :class:`pypeit.core.wavecal.wv_fitting.WaveFit`
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
        maxiter = xfit.size - n_order - 2
        #
        if xfit.size == 0:
            msgs.warn("All points rejected !!")
            return None
        # Fit
        pypeitFit = fitting.robust_fit(xfit/xnspecmin1, yfit, n_order, function=func, maxiter=maxiter,
                                       lower=sigrej_first, upper=sigrej_first, maxrej=1, sticky=True,
                                       minx=fmin, maxx=fmax, weights=wfit)
        # Junk fit?
        if pypeitFit is None:
            msgs.warn("Bad fit!!")
            return None

        rms_ang = pypeitFit.calc_fit_rms(apply_mask=True)
        rms_pix = rms_ang/disp
        if verbose:
            msgs.info('n_order = {:d}'.format(n_order) + ': RMS = {:g}'.format(rms_pix))

        # Reject but keep originals (until final fit)
        ifit = list(ifit[pypeitFit.gpm == 1]) + sv_ifit
        if not input_only:
            # Find new points from the linelist (should we allow removal of the originals?)
            twave = pypeitFit.eval(tcent/xnspecmin1)#, func, minx=fmin, maxx=fmax)
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
    xfit, yfit, wfit = tcent[ifit], all_ids[ifit], weights[ifit]
    pypeitFit = fitting.robust_fit(xfit/xnspecmin1, yfit, n_order, function=func,
                                   lower=sigrej_final, upper=sigrej_final, maxrej=1, sticky=True,
                                   minx=fmin, maxx=fmax, weights=wfit)#, debug=True)
    irej = np.where(np.logical_not(pypeitFit.bool_gpm))[0]
    if len(irej) > 0:
        xrej = xfit[irej]
        yrej = yfit[irej]
        if verbose:
            for kk, imask in enumerate(irej):
                wave = pypeitFit.eval(xrej[kk]/xnspecmin1)#, func, minx=fmin, maxx=fmax)
                msgs.info('Rejecting arc line {:g}; {:g}'.format(yfit[imask], wave))
    else:
        xrej = []
        yrej = []

    ions = all_idsion[ifit]
    # Final RMS
    rms_ang = pypeitFit.calc_fit_rms(apply_mask=True)
    rms_pix = rms_ang/disp

    # Pack up fit
    spec_vec = np.arange(nspec)
    wave_soln = pypeitFit.eval(spec_vec/xnspecmin1)
    cen_wave = pypeitFit.eval(float(nspec)/2/xnspecmin1)
    cen_wave_min1 = pypeitFit.eval((float(nspec)/2 - 1.0)/xnspecmin1)
    cen_disp = cen_wave - cen_wave_min1

    # Ions bit
    ion_bits = np.zeros(len(ions), dtype=WaveFit.bitmask.minimum_dtype())
    for kk,ion in enumerate(ions):
        ion_bits[kk] = WaveFit.bitmask.turn_on(ion_bits[kk], ion.replace(' ', ''))

    # DataContainer time
    # spat_id is set to an arbitrary -1 here and is updated in wavecalib.py
    final_fit = WaveFit(-1, pypeitfit=pypeitFit, pixel_fit=xfit, wave_fit=yfit,
                        ion_bits=ion_bits, xnorm=xnspecmin1,
                        cen_wave=cen_wave, cen_disp=cen_disp,
                        spec=spec, wave_soln = wave_soln, sigrej=sigrej_final,
                        shift=0., tcent=tcent, rms=rms_pix)

    # QA
    if plot_fil is not None:
        autoid.arc_fit_qa(final_fit, plot_fil)
    # Return
    return final_fit
