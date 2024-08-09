
## the following massive function is deprecated
def generate_sensfunc_old(wave, counts, counts_ivar, airmass, exptime, spectrograph, telluric=False, star_type=None,
                      star_mag=None, ra=None, dec=None, std_file = None, BALM_MASK_WID=5., norder=4, nresln=None,debug=False):
    """ Function to generate the sensitivity function.
    This can work in different regimes:
    - If telluric=False and RA=None and Dec=None
      the code creates a sintetic standard star spectrum using the Kurucz models,
      and from this it generates a sens func using nresln=20.0 and masking out
      telluric regions.
    - If telluric=False and RA and Dec are assigned
      the standard star spectrum is extracted from the archive, and a sens func
      is generated using nresln=20.0 and masking out telluric regions.
    - If telluric=True
      the code creates a sintetic standard star spectrum using the Kurucz models,
      the sens func is created setting nresln=1.5 it contains the correction for
      telluric lines.

    Parameters:
    ----------
    wave : array
      Wavelength of the star [no longer with units]
    counts : array
      Flux (in counts) of the star
    counts_ivar : array
      Inverse variance of the star
    airmass : float
      Airmass
    exptime : float
      Exposure time in seconds
    spectrograph : dict
      Instrument specific dict
      Used for extinction correction
    telluric : bool
      if True performs a telluric correction
    star_type : str
      Spectral type of the telluric star (used if telluric=True)
    star_mag : float
      Apparent magnitude of telluric star (used if telluric=True)
    RA : float
      deg, RA of the telluric star
      if assigned, the standard star spectrum will be extracted from
      the archive
    DEC : float
      deg, DEC of the telluric star
      if assigned, the standard star spectrum will be extracted from
      the archive
    BALM_MASK_WID : float
      Mask parameter for Balmer absorption. A region equal to
      BALM_MASK_WID*resln is masked wher resln is the estimate
      for the spectral resolution.
    nresln : float
      Number of resolution elements for break-point placement.
      If assigned, overwrites the settings imposed by the code.
    norder: int
      Order number of polynomial fit.

    Returns:
    -------
    sens_dict : dict
      sensitivity function described by a dict
    """

    # Create copy of the arrays to avoid modification and convert to
    # electrons / s
    wave_star = wave.copy()
    flux_star = counts.copy() / exptime
    ivar_star = counts_ivar.copy() * exptime ** 2

    # Units
    if not isinstance(wave_star, units.Quantity):
        wave_star = wave_star * units.AA

    # ToDo
    # This should be changed. At the moment the extinction correction procedure
    # requires the spectra to be in the optical. For the NIR is probably enough
    # to extend the tables to longer wavelength setting the extinction to 0.0mag.
    msgs.warn("Extinction correction applyed only if the spectra covers <10000Ang.")
    # Apply Extinction if optical bands
    if np.max(wave_star) < 10000. * units.AA:
        msgs.info("Applying extinction correction")
        extinct = load_extinction_data(spectrograph.telescope['longitude'],
                                       spectrograph.telescope['latitude'])
        ext_corr = extinction_correction(wave_star, airmass, extinct)
        # Correct for extinction
        flux_star = flux_star * ext_corr
        ivar_star = ivar_star / ext_corr ** 2
    else:
        msgs.info("Extinction correction not applied")

    # Create star model
    if (ra is not None) and (dec is not None) and (star_mag is None) and (star_type is None):
        # Pull star spectral model from archive
        msgs.info("Get standard model")
        # Grab closest standard within a tolerance
        std_dict = find_standard_file(ra, dec)
        if std_dict is not None:
            # Load standard
            load_standard_file(std_dict)
            # Interpolate onto observed wavelengths
            #std_xspec = XSpectrum1D.from_tuple((std_dict['wave'], std_dict['flux']))
            debugger.set_trace()
            xspec = std_xspec.rebin(wave_star)  # Conserves flambda
            flux_true = xspec.flux.value
        else:
            msgs.error('No spectrum found in our database for your standard star. Please use another standard star \
                       or consider add it into out database.')
    elif (star_mag is not None) and (star_type is not None):
        # Create star spectral model
        msgs.info("Creating standard model")
        # Create star model
        star_loglam, star_flux, std_dict = telluric_sed(star_mag, star_type)
        star_lam = 10 ** star_loglam
        # Generate a dict matching the output of find_standard_file
        std_dict = dict(cal_file='KuruczTelluricModel', name=star_type, fmt=1,
                        std_ra=None, std_dec=None)
        std_dict['wave'] = star_lam * units.AA
        std_dict['flux'] = 1e17 * star_flux * units.erg / units.s / units.cm ** 2 / units.AA
        # ToDO If the Kuruck model is used, rebin create weird features
        # I using scipy interpolate to avoid this
        flux_true = scipy.interpolate.interp1d(std_dict['wave'], std_dict['flux'],
                                               bounds_error=False,
                                               fill_value='extrapolate')(wave_star)
    else:
        debugger.set_trace()
        msgs.error('Insufficient information provided for fluxing. '
                   'Either the coordinates of the standard or a stellar type and magnitude are needed.')


    if np.min(flux_true) <= 0.:
        msgs.warn('Your spectrum extends beyond calibrated standard star, extrapolating the spectra with polynomial.')
        # ToDo: should we extrapolate it using graybody model?
        mask_model = flux_true<=0
        msk_poly, poly_coeff = utils.robust_polyfit_djs(std_dict['wave'].value, std_dict['flux'].value,8,function='polynomial',
                                                    invvar=None, guesses=None, maxiter=50, inmask=None, sigma=None, \
                                                    lower=3.0, upper=3.0, maxdev=None, maxrej=3, groupdim=None,
                                                    groupsize=None,groupbadpix=False, grow=0, sticky=True, use_mad=True)
        star_poly = utils.func_val(poly_coeff, wave_star.value, 'polynomial')
        #flux_true[mask_model] = star_poly[mask_model]
        flux_true = star_poly.copy()
        if debug:
            plt.plot(std_dict['wave'], std_dict['flux'],'bo',label='Raw Star Model')
            plt.plot(std_dict['wave'],  utils.func_val(poly_coeff, std_dict['wave'].value, 'polynomial'), 'k-',label='robust_poly_fit')
            plt.plot(wave_star,flux_true,'r-',label='Your Final Star Model used for sensfunc')
            plt.show()

    # Set nresln
    if nresln is None:
        if telluric:
            nresln = 1.5
            msgs.info("Set nresln to 1.5")
        else:
            nresln = 20.0
            msgs.info("Set nresln to 20.0")

    # ToDo
    # Compute an effective resolution for the standard. This could be improved
    # to setup an array of breakpoints based on the resolution. At the
    # moment we are using only one number
    msgs.work("Should pull resolution from arc line analysis")
    msgs.work("At the moment the resolution is taken as the PixelScale")
    msgs.work("This needs to be changed!")
    std_pix = np.median(np.abs(wave_star - np.roll(wave_star, 1)))
    std_res = std_pix
    resln = std_res
    if (nresln * std_res) < std_pix:
        msgs.warn("Bspline breakpoints spacing shoud be larger than 1pixel")
        msgs.warn("Changing input nresln to fix this")
        nresln = std_res / std_pix

    # Mask bad pixels, edges, and Balmer, Paschen, Brackett, and Pfund lines
    # Mask (True = good pixels)
    msgs.info("Masking spectral regions:")
    msk_star = np.ones_like(flux_star).astype(bool)

    # Mask bad pixels
    msgs.info(" Masking bad pixels")
    msk_star[ivar_star <= 0.] = False
    msk_star[flux_star <= 0.] = False

    # Mask edges
    msgs.info(" Masking edges")
    msk_star[:1] = False
    msk_star[-1:] = False

    # Mask Balmer
    msgs.info(" Masking Balmer")
    lines_balm = np.array([3836.4, 3969.6, 3890.1, 4102.8, 4102.8, 4341.6, 4862.7, 5407.0,
                           6564.6, 8224.8, 8239.2]) * units.AA
    for line_balm in lines_balm:
        ibalm = np.abs(wave_star - line_balm) <= BALM_MASK_WID * resln
        msk_star[ibalm] = False

    # Mask Paschen
    msgs.info(" Masking Paschen")
    # air wavelengths from:
    # https://www.subarutelescope.org/Science/Resources/lines/hi.html
    lines_pasc = np.array([8203.6, 9229.0, 9546.0, 10049.4, 10938.1,
                           12818.1, 18751.0]) * units.AA
    for line_pasc in lines_pasc:
        ipasc = np.abs(wave_star - line_pasc) <= BALM_MASK_WID * resln
        msk_star[ipasc] = False

    # Mask Brackett
    msgs.info(" Masking Brackett")
    # air wavelengths from:
    # https://www.subarutelescope.org/Science/Resources/lines/hi.html
    lines_brac = np.array([14584.0, 18174.0, 19446.0, 21655.0,
                           26252.0, 40512.0]) * units.AA
    for line_brac in lines_brac:
        ibrac = np.abs(wave_star - line_brac) <= BALM_MASK_WID * resln
        msk_star[ibrac] = False

    # Mask Pfund
    msgs.info(" Masking Pfund")
    # air wavelengths from:
    # https://www.subarutelescope.org/Science/Resources/lines/hi.html
    lines_pfund = np.array([22788.0, 32961.0, 37395.0, 46525.0,
                            74578.0]) * units.AA
    for line_pfund in lines_pfund:
        ipfund = np.abs(wave_star - line_pfund) <= BALM_MASK_WID * resln
        msk_star[ipfund] = False

    # Mask Atm. cutoff
    msgs.info(" Masking Below the atmospheric cutoff")
    atms_cutoff = wave_star <= 3000.0 * units.AA
    msk_star[atms_cutoff] = False

    #if ~telluric: #Feige:  This is a bug
    if not telluric:
        # Mask telluric absorption
        msgs.info("Masking Telluric")
        tell = np.any([((wave_star >= 7580.00 * units.AA) & (wave_star <= 7750.00 * units.AA)),
                       ((wave_star >= 7160.00 * units.AA) & (wave_star <= 7340.00 * units.AA)),
                       ((wave_star >= 6860.00 * units.AA) & (wave_star <= 6930.00 * units.AA)),
                       ((wave_star >= 9310.00 * units.AA) & (wave_star <= 9665.00 * units.AA)),
                       ((wave_star >= 11120.0 * units.AA) & (wave_star <= 11615.0 * units.AA)),
                       ((wave_star >= 12610.0 * units.AA) & (wave_star <= 12720.0 * units.AA)),
                       ((wave_star >= 13160.0 * units.AA) & (wave_star <= 15065.0 * units.AA)),
                       ((wave_star >= 15700.0 * units.AA) & (wave_star <= 15770.0 * units.AA)),
                       ((wave_star >= 16000.0 * units.AA) & (wave_star <= 16100.0 * units.AA)),
                       ((wave_star >= 16420.0 * units.AA) & (wave_star <= 16580.0 * units.AA)),
                       ((wave_star >= 17310.0 * units.AA) & (wave_star <= 20775.0 * units.AA)),
                       (wave_star >= 22680.0 * units.AA)], axis=0)
        msk_star[tell] = False

    # Apply mask
    ivar_star[~msk_star] = 0.0

    # Fit in magnitudes
    kwargs_bspline = {'bkspace': resln.value * nresln}
    kwargs_reject = {'maxrej': 5}
    sensfunc, sensfit = bspline_magfit(wave_star.value, flux_star, ivar_star, flux_true, inmask=msk_star,
                              kwargs_bspline=kwargs_bspline, kwargs_reject=kwargs_reject,debug=debug)

    #Cleaning sensfunc
    ## ToDo: currently I'm fitting the sensfunc in the masked region with a polynomial, should we change the algorithm to
    ##   fit polynomial first and then bsline the poly-subtracted flux ???
    ## keep tell free region for poly fit. tell2 is different from tell since tell2 include more small trunk of telluric free
    ## regions. tell2 might be not suitable for the bspline fitting. We need to select a more robust telluric region for both purpose.
    tell2 = np.any([((wave_star >= 7580.00 * units.AA) & (wave_star <= 7750.00 * units.AA)),
                   ((wave_star >= 7160.00 * units.AA) & (wave_star <= 7340.00 * units.AA)),
                   ((wave_star >= 6860.00 * units.AA) & (wave_star <= 6930.00 * units.AA)),
                   ((wave_star >= 9310.00 * units.AA) & (wave_star <= 9665.00 * units.AA)),
                   ((wave_star >= 11120.0 * units.AA) & (wave_star <= 11545.0 * units.AA)),
                   ((wave_star >= 12610.0 * units.AA) & (wave_star <= 12720.0 * units.AA)),
                   ((wave_star >= 13400.0 * units.AA) & (wave_star <= 14830.0 * units.AA)),
                   ((wave_star >= 15700.0 * units.AA) & (wave_star <= 15770.0 * units.AA)),
                   ((wave_star >= 16000.0 * units.AA) & (wave_star <= 16100.0 * units.AA)),
                   ((wave_star >= 16420.0 * units.AA) & (wave_star <= 16580.0 * units.AA)),
                   ((wave_star >= 17630.0 * units.AA) & (wave_star <= 19690.0 * units.AA)),
                   ((wave_star >= 19790.0 * units.AA) & (wave_star <= 19810.0 * units.AA)),
                   ((wave_star >= 19950.0 * units.AA) & (wave_star <= 20310.0 * units.AA)),
                   ((wave_star >= 20450.0 * units.AA) & (wave_star <= 20920.0 * units.AA)),
                   ((wave_star >= 24000.0 * units.AA) & (wave_star <= 24280.0 * units.AA)),
                   ((wave_star >= 24320.0 * units.AA) & (wave_star <= 24375.0 * units.AA)),
                   (wave_star >= 24450.0 * units.AA)], axis=0)
    msk_all = msk_star.copy() # mask for polynomial fitting
    msk_sens = msk_star.copy() # mask for sensfunc
    med, mad = utils.robust_meanstd(sensfunc)

    msk_crazy = (sensfunc<=0) | (sensfunc>1e3*med)
    msk_all[tell2] = False
    msk_all[msk_crazy] = False
    msk_sens[msk_crazy] = False

    if (len(wave_star.value[msk_all]) < norder+1) or (len(wave_star.value[msk_all]) < 0.1*len(wave_star.value)):
        msgs.warn('It seems this order/spectrum well within the telluric region. No polynomial fit will be performed.')
    else:
        #polyfit the sensfunc
        msk_poly, poly_coeff = utils.robust_polyfit_djs(wave_star.value[msk_all],np.log10(sensfunc[msk_all]), norder, function='polynomial',
                                               invvar=None,guesses = None, maxiter = 50, inmask = None, sigma = None,\
                                               lower = 3.0, upper = 3.0,maxdev=None,maxrej=3,groupdim=None,groupsize=None,\
                                               groupbadpix=False, grow=0,sticky=True,use_mad=True)
        sensfunc_poly = 10**(utils.func_val(poly_coeff, wave_star.value, 'polynomial'))
        sensfunc[~msk_sens] =  sensfunc_poly[~msk_sens]
        if debug:
            plt.rcdefaults()
            plt.rcParams['font.family'] = 'times new roman'
            plt.plot(wave_star.value[~msk_sens], sensfunc[~msk_sens], 'bo')
            plt.plot(wave_star.value, sensfunc_poly, 'r-',label='Polyfit')
            plt.plot(wave_star.value, sensfunc, 'k-',label='bspline fitting')
            plt.ylim(0.0, 100.0)
            plt.legend()
            plt.xlabel('Wavelength [ang]')
            plt.ylabel('Sensfunc')
            plt.show()
            plt.close()

            plt.figure(figsize=(10, 6))
            plt.clf()
            plt.plot(wave_star.value,flux_star*sensfunc, label='Calibrated Spectrum')
            plt.plot(wave_star.value,flux_true, label='Model')
            plt.plot(wave_star.value,np.sqrt(1/ivar_star))
            plt.legend()
            plt.xlabel('Wavelength [ang]')
            plt.ylabel('Flux [erg/s/cm2/Ang.]')
            plt.ylim(0,np.median(flux_true)*2.5)
            plt.title('Final corrected spectrum')
            plt.show()
            plt.close()


    # JFH Left off here.
    # Creating the dict
    #msgs.work("Is min, max and wave_min, wave_max a duplicate?")
    #sens_dict = dict(wave=wave_sens, sensfunc=sensfunc, min=None, max=None, std=std_dict)

    # Add in wavemin,wavemax
    sens_dict = {}
    sens_dict['wave'] = wave_star
    sens_dict['sensfunc'] = sensfunc
    sens_dict['wave_min'] = np.min(wave_star)
    sens_dict['wave_max'] = np.max(wave_star)
    sens_dict['exptime']= exptime
    sens_dict['airmass']= airmass
    sens_dict['std_file']= std_file
    # Get other keys from standard dict
    sens_dict['std_ra'] = std_dict['std_ra']
    sens_dict['std_dec'] = std_dict['std_dec']
    sens_dict['std_name'] = std_dict['name']
    sens_dict['cal_file'] = std_dict['cal_file']
    sens_dict['flux_true'] = flux_true
    #sens_dict['std_dict'] = std_dict
    #sens_dict['msk_star'] = msk_star
    #sens_dict['mag_set'] = mag_set

    return sens_dict


## bspline_magfit is deprecated at this moment.
def bspline_magfit(wave, flux, ivar, flux_std, inmask=None, maxiter=35, upper=2, lower=2,
                   kwargs_bspline={}, kwargs_reject={}, debug=False, show_QA=False):
    """
    Perform a bspline fit to the flux ratio of standard to
    observed counts. Used to generate a sensitivity function.

    Parameters
    ----------
    wave : ndarray
      wavelength as observed
    flux : ndarray
      counts/s as observed
    ivar : ndarray
      inverse variance
    flux_std : Quantity array
      standard star true flux (erg/s/cm^2/A)
    inmask : ndarray
      bspline mask
    maxiter : integer
      maximum number of iterations for bspline_iterfit
    upper : integer
      number of sigma for rejection in bspline_iterfit
    lower : integer
      number of sigma for rejection in bspline_iterfit
    kwargs_bspline : dict, optional
      keywords for bspline_iterfit
    kwargs_reject : dict, optional
      keywords for bspline_iterfit
    debug : bool
      if True shows some dubugging plots

    Returns
    -------
    bset_log1
    """
    # Create copy of the arrays to avoid modification
    wave_obs = wave.copy()
    flux_obs = flux.copy()
    ivar_obs = ivar.copy()

    # preparing arrays to run in bspline_iterfit

    if np.all(~np.isfinite(ivar_obs)):
        msgs.warn("NaN are present in the inverse variance")

    # Preparing arrays to run in bspline_iterfit

    if np.all(~np.isfinite(ivar_obs)):
        msgs.warn("NaN are present in the inverse variance")

    # Removing outliers

    # Calculate log of flux_obs setting a floor at TINY
    logflux_obs = 2.5 * np.log10(np.maximum(flux_obs, TINY))
    # Set a fix value for the variance of logflux
    logivar_obs = np.ones_like(logflux_obs) * (10.0 ** 2)

    # Calculate log of flux_std model setting a floor at TINY
    logflux_std = 2.5 * np.log10(np.maximum(flux_std, TINY))

    # Calculate ratio setting a floor at MAGFUNC_MIN and a ceiling at
    # MAGFUNC_MAX
    magfunc = logflux_std - logflux_obs
    magfunc = np.maximum(np.minimum(magfunc, MAGFUNC_MAX), MAGFUNC_MIN)
    magfunc_mask = (magfunc < 0.99 * MAGFUNC_MAX) & (magfunc > 0.99 * MAGFUNC_MIN)

    # Mask outliners
    # masktot=True means good pixel
    if inmask is None:
        masktot = (ivar_obs > 0.0) & np.isfinite(logflux_obs) & np.isfinite(ivar_obs) & \
                  np.isfinite(logflux_std) & magfunc_mask
    else:
        masktot = inmask & (ivar_obs > 0.0) & np.isfinite(logflux_obs) & np.isfinite(ivar_obs) & \
                  np.isfinite(logflux_std) & magfunc_mask
    logivar_obs[~masktot] = 0.

    # Calculate sensfunc
    sensfunc = 10.0 ** (0.4 * magfunc)

    msgs.info("Initialize bspline for flux calibration")

    init_bspline = pydl.bspline(wave_obs, bkspace=kwargs_bspline['bkspace'])
    fullbkpt = init_bspline.breakpoints

    # TESTING turning off masking for now
    # remove masked regions from breakpoints
    msk_obs = np.ones_like(wave_obs).astype(bool)
    msk_obs[~masktot] = False
    msk_bkpt = scipy.interpolate.interp1d(wave_obs, msk_obs, kind='nearest', fill_value='extrapolate')(fullbkpt)
    init_breakpoints = fullbkpt[msk_bkpt > 0.999]

    # init_breakpoints = fullbkpt

    #  First round of the fit:
    msgs.info("Bspline fit: step 1")
    bset1, bmask = pydl.iterfit(wave_obs, magfunc, invvar=logivar_obs, inmask=masktot, upper=upper, lower=lower,
                                fullbkpt=init_breakpoints, maxiter=maxiter, kwargs_bspline=kwargs_bspline,
                                kwargs_reject=kwargs_reject)
    logfit1, _ = bset1.value(wave_obs)
    logfit_bkpt, _ = bset1.value(init_breakpoints)

    if debug:
        # Check for calibration
        plt.figure(1)
        plt.plot(wave_obs, magfunc, drawstyle='steps-mid', color='black', label='magfunc')
        plt.plot(wave_obs, logfit1, color='cornflowerblue', label='logfit1')
        plt.plot(wave_obs[~masktot], magfunc[~masktot], '+', color='red', markersize=5.0, label='masked magfunc')
        plt.plot(wave_obs[~masktot], logfit1[~masktot], '+', color='red', markersize=5.0, label='masked logfit1')
        plt.plot(init_breakpoints, logfit_bkpt, '.', color='green', markersize=4.0, label='breakpoints')
        plt.plot(init_breakpoints, np.interp(init_breakpoints, wave_obs, magfunc), '.', color='green', markersize=4.0,
                 label='breakpoints')
        plt.plot(wave_obs, 1.0 / np.sqrt(logivar_obs), color='orange', label='sigma')
        plt.legend()
        plt.xlabel('Wavelength [ang]')
        plt.ylim(0.0, 1.2 * MAGFUNC_MAX)
        plt.title('1st Bspline fit')
        plt.show()

    modelfit1 = np.power(10.0, 0.4 * np.maximum(np.minimum(logfit1, MAGFUNC_MAX), MAGFUNC_MIN))
    residual = sensfunc / (modelfit1 + (modelfit1 == 0)) - 1.
    # new_mask = masktot & (sensfunc > 0)

    # residual_ivar = (modelfit1 * flux_obs / (sensfunc + (sensfunc == 0.0))) ** 2 * ivar_obs
    residual_ivar = np.ones_like(residual) / (0.1 ** 2)
    residual_ivar = residual_ivar * masktot

    (mean, med, stddev) = sigma_clipped_stats(residual[masktot], sigma_lower=3.0, sigma_upper=3.0)

    if np.median(stddev > 0.01):
        #  Second round of the fit:
        msgs.info("Bspline fit: step 2")
        #  Now do one more fit to the ratio of data/model - 1.
        bset_residual, bmask2 = pydl.iterfit(wave_obs, residual, invvar=residual_ivar, inmask=masktot, upper=upper,
                                             lower=lower, maxiter=maxiter, fullbkpt=bset1.breakpoints,
                                             kwargs_bspline=kwargs_bspline, kwargs_reject=kwargs_reject)
        bset_log1 = bset1.copy()
        bset_log1.coeff = bset_log1.coeff + bset_residual.coeff
        if debug:
            # Check for calibration
            resid_fit, _ = bset_residual.value(wave_obs)
            logfit2, _ = bset_log1.value(wave_obs)
            logfit2_bkpt, _ = bset_log1.value(bset1.breakpoints)
            plt.figure(1)
            plt.plot(wave_obs, residual, drawstyle='steps-mid', color='black', label='residual')
            plt.plot(wave_obs, resid_fit, color='cornflowerblue', label='resid_fit')
            plt.plot(wave_obs[~masktot], residual[~masktot], '+', color='red', markersize=5.0, label='masked residual')
            plt.plot(wave_obs[~masktot], resid_fit[~masktot], '+', color='red', markersize=5.0, label='masked resid_fit')
            plt.plot(init_breakpoints, logfit2_bkpt, '.', color='green', markersize=4.0, label='breakpoints')
            plt.plot(wave_obs, 1.0 / np.sqrt(residual_ivar), color='orange', label='sigma')
            plt.legend()
            plt.xlabel('Wavelength [ang]')
            plt.ylim(-0.1, 0.1)
            plt.title('2nd Bspline fit')
            plt.show()
    else:
        bset_log1 = bset1.copy()

    # ToDo JFH I think we should move towards writing this out as a vector in a fits table
    # rather than the b-spline.

    # Create sensitivity function
    newlogfit, _ = bset_log1.value(wave_obs)
    sensfit = np.power(10.0, 0.4 * np.maximum(np.minimum(newlogfit, MAGFUNC_MAX), MAGFUNC_MIN))
    sensfit[~magfunc_mask] = 0.0

    if debug:

        # Check for calibration
        plt.figure(1)
        plt.plot(wave_obs, sensfunc, drawstyle='steps-mid', color='black', label='sensfunc')
        plt.plot(wave_obs, sensfit, color='cornflowerblue', label='sensfunc fit')
        plt.plot(wave_obs[~masktot], sensfunc[~masktot], '+', color='red', markersize=5.0, label='masked sensfunc')
        plt.plot(wave_obs[~masktot], sensfit[~masktot], '+', color='red', markersize=5.0, label='masked sensfuncfit')
        plt.legend()
        plt.xlabel('Wavelength [ang]')
        plt.ylim(0.0, 100.0)
        plt.show()

    # Check quality of the fit
    absdev = np.median(np.abs(sensfit / modelfit1 - 1))
    msgs.info('Difference between fits is {:g}'.format(absdev))

    # Check for residual of the fit
    if debug:
        # scale = np.power(10.0, 0.4 * sensfit)
        flux_cal = flux_obs * sensfit
        ivar_cal = ivar_obs / sensfit ** 2.

        plt.rcdefaults()
        plt.rcParams['font.family']= 'times new roman'
        plt.figure(figsize=(11, 8.5))
        plt.clf()
        plt.plot(wave_obs,flux_cal, label='Calibrated Spectrum')
        plt.plot(wave_obs,flux_std, label='Model')
        plt.plot(wave_obs,np.sqrt(1/ivar_cal))
        plt.legend()
        plt.xlabel('Wavelength [ang]')
        plt.ylabel('Flux [erg/s/cm2/Ang.]')
        plt.ylim(0,np.median(flux_std)*2.5)
        plt.show()
        plt.close()

    # QA
    msgs.work("Add QA for sensitivity function")
    if show_QA:
        qa_bspline_magfit(wave_obs, bset_log1, magfunc, masktot)

    return sensfunc,sensfit

def qa_bspline_magfit(wave, bset, magfunc, mask):
    plt.close("all")
    plt.rcParams['savefig.dpi'] = 600
    plt.rcParams['xtick.top'] = True
    plt.rcParams['ytick.right'] = True
    plt.rcParams['xtick.minor.visible'] = True
    plt.rcParams['ytick.minor.visible'] = True
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['xtick.major.size'] = 6
    plt.rcParams['ytick.major.size'] = 6
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['axes.linewidth'] = 1
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['legend.frameon'] = False
    plt.rcParams['legend.handletextpad'] = 1.0
    final_fit, _ = bset.value(wave)
    final_fit_bkpt, _ = bset.value(bset.breakpoints)

    plt.figure(1)
    plt.plot(bset.breakpoints, final_fit_bkpt, '.', color='green', markersize=4.0, label='breakpoints')
    plt.plot(wave, magfunc, drawstyle='steps-mid', color='black', label='magfunc')
    plt.plot(wave, final_fit, color='cornflowerblue', label='bspline fit')
    plt.plot(wave[~mask], magfunc[~mask], '+', color='red', markersize=5.0, label='masked points')
    plt.legend()
    plt.xlabel('Wavelength [ang]')
    plt.title('Final Result of the Bspline fit')

    plt.show()
    return
