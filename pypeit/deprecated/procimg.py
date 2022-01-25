
def lacosmic(sciframe, saturation, nonlinear, varframe=None, maxiter=1, grow=1.5,
             remove_compact_obj=True, sigclip=5.0, sigfrac=0.3, objlim=5.0):
    """
    Identify cosmic rays using the L.A.Cosmic algorithm
    U{http://www.astro.yale.edu/dokkum/lacosmic/}
    (article : U{http://arxiv.org/abs/astro-ph/0108003})
    This routine is mostly courtesy of Malte Tewes

    Args:
        sciframe:
        saturation:
        nonlinear:
        varframe:
        maxiter:
        grow:
        remove_compact_obj:
        sigclip (float):
            Threshold for identifying a CR
        sigfrac:
        objlim:

    Returns:
        ndarray: mask of cosmic rays (0=no CR, 1=CR)

    """
    msgs.info("Detecting cosmic rays with the L.A.Cosmic algorithm")
#    msgs.work("Include these parameters in the settings files to be adjusted by the user")
    # Set the settings
    scicopy = sciframe.copy()
    crmask = np.cast['bool'](np.zeros(sciframe.shape))
    sigcliplow = sigclip * sigfrac

    # Determine if there are saturated pixels
    satpix = np.zeros_like(sciframe)
#    satlev = settings_det['saturation']*settings_det['nonlinear']
    satlev = saturation*nonlinear
    wsat = np.where(sciframe >= satlev)
    if wsat[0].size == 0: satpix = None
    else:
        satpix[wsat] = 1.0
        satpix = np.cast['bool'](satpix)

    # Define the kernels
    laplkernel = np.array([[0.0, -1.0, 0.0], [-1.0, 4.0, -1.0], [0.0, -1.0, 0.0]])  # Laplacian kernal
    growkernel = np.ones((3,3))
    for i in range(1, maxiter+1):
        msgs.info("Convolving image with Laplacian kernel")
        # Subsample, convolve, clip negative values, and rebin to original size
        subsam = utils.subsample(scicopy)
        conved = signal.convolve2d(subsam, laplkernel, mode="same", boundary="symm")
        cliped = conved.clip(min=0.0)
        lplus = utils.rebin_evlist(cliped, np.array(cliped.shape)/2.0)

        msgs.info("Creating noise model")
        # Build a custom noise map, and compare  this to the laplacian
        m5 = ndimage.filters.median_filter(scicopy, size=5, mode='mirror')
        if varframe is None:
            noise = np.sqrt(np.abs(m5))
        else:
            noise = np.sqrt(varframe)
        msgs.info("Calculating Laplacian signal to noise ratio")

        # Laplacian S/N
        s = lplus / (2.0 * noise)  # Note that the 2.0 is from the 2x2 subsampling

        # Remove the large structures
        sp = s - ndimage.filters.median_filter(s, size=5, mode='mirror')

        msgs.info("Selecting candidate cosmic rays")
        # Candidate cosmic rays (this will include HII regions)
        candidates = sp > sigclip
        nbcandidates = np.sum(candidates)

        msgs.info("{0:5d} candidate pixels".format(nbcandidates))

        # At this stage we use the saturated stars to mask the candidates, if available :
        if satpix is not None:
            msgs.info("Masking saturated pixels")
            candidates = np.logical_and(np.logical_not(satpix), candidates)
            nbcandidates = np.sum(candidates)

            msgs.info("{0:5d} candidate pixels not part of saturated stars".format(nbcandidates))

        msgs.info("Building fine structure image")

        # We build the fine structure image :
        m3 = ndimage.filters.median_filter(scicopy, size=3, mode='mirror')
        m37 = ndimage.filters.median_filter(m3, size=7, mode='mirror')
        f = m3 - m37
        f /= noise
        f = f.clip(min=0.01)

        msgs.info("Removing suspected compact bright objects")

        # Now we have our better selection of cosmics :

        if remove_compact_obj:
            cosmics = np.logical_and(candidates, sp/f > objlim)
        else:
            cosmics = candidates
        nbcosmics = np.sum(cosmics)

        msgs.info("{0:5d} remaining candidate pixels".format(nbcosmics))

        # What follows is a special treatment for neighbors, with more relaxed constains.

        msgs.info("Finding neighboring pixels affected by cosmic rays")

        # We grow these cosmics a first time to determine the immediate neighborhod  :
        growcosmics = np.cast['bool'](signal.convolve2d(np.cast['float32'](cosmics), growkernel, mode="same", boundary="symm"))

        # From this grown set, we keep those that have sp > sigmalim
        # so obviously not requiring sp/f > objlim, otherwise it would be pointless
        growcosmics = np.logical_and(sp > sigclip, growcosmics)

        # Now we repeat this procedure, but lower the detection limit to sigmalimlow :

        finalsel = np.cast['bool'](signal.convolve2d(np.cast['float32'](growcosmics), growkernel, mode="same", boundary="symm"))
        finalsel = np.logical_and(sp > sigcliplow, finalsel)

        # Unmask saturated pixels:
        if satpix is not None:
            msgs.info("Masking saturated stars")
            finalsel = np.logical_and(np.logical_not(satpix), finalsel)

        ncrp = np.sum(finalsel)

        msgs.info("{0:5d} pixels detected as cosmics".format(ncrp))

        # We find how many cosmics are not yet known :
        newmask = np.logical_and(np.logical_not(crmask), finalsel)
        nnew = np.sum(newmask)

        # We update the mask with the cosmics we have found :
        crmask = np.logical_or(crmask, finalsel)

        msgs.info("Iteration {0:d} -- {1:d} pixels identified as cosmic rays ({2:d} new)".format(i, ncrp, nnew))
        if ncrp == 0:
            break

    # Additional algorithms (not traditionally implemented by LA cosmic) to
    # remove some false positives.
    msgs.work("The following algorithm would be better on the rectified, tilts-corrected image")
    filt  = ndimage.sobel(sciframe, axis=1, mode='constant')
    filty = ndimage.sobel(filt/np.sqrt(np.abs(sciframe)), axis=0, mode='constant')
    filty[np.where(np.isnan(filty))]=0.0

    sigimg = cr_screen(filty)

    sigsmth = ndimage.filters.gaussian_filter(sigimg,1.5)
    sigsmth[np.where(np.isnan(sigsmth))]=0.0
    sigmask = np.cast['bool'](np.zeros(sciframe.shape))
    sigmask[np.where(sigsmth>sigclip)] = True
    crmask = np.logical_and(crmask, sigmask)
    msgs.info("Growing cosmic ray mask by 1 pixel")
    crmask = grow_mask(crmask, grow)

    return crmask.astype(bool)



# TODO: Add sigdev to the high-level parameter set so that it can be
# changed by the user?
def find_bad_pixels(bias, numamplifiers, datasec, sigdev=10.0, trim=True):
    """
    Identify bad pixels in the datasection of the bias frame based on
    their robust deviation from the median.

    Args:
        bias (:obj:`numpy.ndarray`):
            Bias frame
        numamplifiers (int):
            Number of amplifiers
        datasec (list):
            List of slices, one per amplifier, that contain the data in
            the raw frame.  The slices and be lists of slice ojects or
            strings.  If they are strings, :func:`parse.sec2slice` is
            used to convert them for use in the function.
        sigdev (:obj:`float`, optional):
            Number of robust standard deviations beyond which to flag
            pixels as bad.
        trim (:obj:`bool`, optional):
            Flag to trim image to the data section.

    Returns:
        :obj:`numpy.ndarray`: An integer array with bad pixels set to 1.

    Raises:
        ValueError:
            Raised if the number of data sections does not match the
            number of amplifiers.
    """
    # Check input
    if len(datasec) != numamplifiers:
        raise ValueError('Number of amplifiers does not match provided data sections.')

    # If the input image sections are strings, convert them
    if isinstance(datasec[0], str):
        _datasec = datasec.copy()
        for i in range(numamplifiers):
            _datasec[i] = parse.sec2slice(datasec[i], require_dim=2)
    else:
        _datasec = datasec

    # Find the bad pixels and mask them
    mask = np.zeros_like(bias, dtype=np.int8)
    is_data = np.zeros_like(bias, dtype=bool)
    for i in range(numamplifiers):
        is_data[datasec[i]] = True
        temp = np.abs(np.median(bias[datasec[i]])-bias[datasec[i]])
        sigval = 1.4826*max(np.median(temp), 1)
        mask[datasec[i]][temp > sigdev*sigval] = 1

    msgs.info("Identified {0:d} bad pixels".format(int(np.sum(mask))))
    return trim_frame(mask, np.invert(is_data)) if trim else mask

#def badpix(frame, numamplifiers, datasec, sigdev=10.0):
#    """
#    frame is a master bias frame
#    numamplifiers : int
#    datasec : list
#    sigdev is the number of standard deviations away from the median that a pixel needs to be in order to be classified as a bad pixel
#    """
#    bpix = np.zeros_like(frame, dtype=np.int)
#    subfr, tframe, temp = None, None, None
#    #for i in range(settings.spect[dnum]['numamplifiers']):
#    for i in range(numamplifiers):
#        #datasec = "datasec{0:02d}".format(i+1)
#        x0, x1 = datasec[i][0][0], datasec[i][0][1]
#        y0, y1 = datasec[i][1][0], datasec[i][1][1]
#        #x0, x1 = settings.spect[dnum][datasec][0][0], settings.spect[dnum][datasec][0][1]
#        #y0, y1 = settings.spect[dnum][datasec][1][0], settings.spect[dnum][datasec][1][1]
#        xv = np.arange(x0, x1)
#        yv = np.arange(y0, y1)
#        # Construct an array with the rows and columns to be extracted
#        w = np.ix_(xv,yv)
#        tframe = frame[w]
#        temp = np.abs(np.median(tframe)-tframe)
#        sigval = max(np.median(temp)*1.4826, 1.4826)
#        ws = np.where(temp > sigdev*sigval)
#        subfr = np.zeros(tframe.shape, dtype=np.int)
#        subfr[ws] = 1
#        bpix[w] = subfr
#    del subfr, tframe, temp
#    # Finally, trim the bad pixel frame
#    bpix = trim(bpix, numamplifiers, datasec)
#    msgs.info("Identified {0:d} bad pixels".format(int(np.sum(bpix))))
#    return bpix


#def bias_subtract(rawframe, msbias, numamplifiers=None, datasec=None, oscansec=None):
#    """ Core routine for bias subtraction
#    Calls sub_overscan if msbias == 'overscan'
#
#    Parameters
#    ----------
#    rawframe : ndarray
#    msbias : ndarray or str
#
#    Returns
#    -------
#    newframe : ndarray
#      Bias subtracted frame
#
#    """
#    if isinstance(msbias, np.ndarray):
#        msgs.info("Subtracting bias image from raw frame")
#        # Subtract the master bias frame
#        return rawframe-msbias
#    elif isinstance(msbias, basestring) and msbias == 'overscan':
#        return subtract_overscan(rawframe, numamplifiers, datasec, oscansec, settings=None)
#
#    msgs.error('Could not subtract bias level with the input bias approach.')


'''
def error_frame_postext(sciframe, idx, fitsdict, settings_spect):
    # Dark Current noise
    dnoise = settings.spect['det']['darkcurr'] * float(fitsdict["exptime"][idx])/3600.0
    # The effective read noise
    rnoise = settings.spect['det']['ronoise']**2 + (0.5*settings.spect['det']['gain'])**2
    errframe = np.zeros_like(sciframe)
    w = np.where(sciframe != -999999.9)
    errframe[w] = np.sqrt(sciframe[w] + rnoise + dnoise)
    w = np.where(sciframe == -999999.9)
    errframe[w] = 999999.9
    return errframe
'''

'''
def get_datasec_trimmed(spectrograph, scifile, det, settings_det,
                        naxis0=None, naxis1=None):
    """
    Primarily a wrapper with calls to get_datasec and pix_to_amp()

    Parameters
    ----------
    spectrograph : str
    scifile : str
    numamplifiers : int
    det : int
    settings_det : dict
    naxis0 : int, optional
    naxis1 : int, optional

    Returns
    -------
    datasec_img : ndarray
    naxis0 : int
    naxis1 : int
    """
    # Instrument specific bits
    # TODO -- Remove instrument specific items in a method like this
    if spectrograph in ['keck_lris_blue', 'keck_lris_red', 'keck_deimos']:
        # Grab
        datasec, oscansec, naxis0, naxis1 = get_datasec(spectrograph, scifile,
                                                        numamplifiers=settings_det['numamplifiers'], det=det)
        # Fill (for backwards compatability)
        for kk in range(settings_det['numamplifiers']):
            sdatasec = "datasec{0:02d}".format(kk+1)
            settings_det[sdatasec] = datasec[kk]
            soscansec = "oscansec{0:02d}".format(kk+1)
            settings_det[soscansec] = oscansec[kk]
        #fitstbl['naxis0'][scidx] = naxis0
        #fitstbl['naxis1'][scidx] = naxis1

    # Build the datasec lists for pix_to_amp
    datasec = []
    for i in range(settings_det['numamplifiers']):
        sdatasec = "datasec{0:02d}".format(i+1)
        datasec.append(settings_det[sdatasec])
    # Call
    #naxis0, naxis1 = int(fitstbl['naxis0'][scidx]), int(fitstbl['naxis1'][scidx])
    datasec_img = arpixels.pix_to_amp(naxis0, naxis1, datasec, settings_det['numamplifiers'])
    return datasec_img, naxis0, naxis1
'''


'''
def sn_frame(slf, sciframe, idx):
    # Dark Current noise
    dnoise = settings.spect['det']['darkcurr'] * float(slf._fitsdict["exptime"][idx])/3600.0
    # The effective read noise
    rnoise = np.sqrt(settings.spect['det']['ronoise']**2 + (0.5*settings.spect['det']['gain'])**2)
    errframe = np.abs(sciframe) + rnoise + dnoise
    # If there are negative pixels, mask them as bad pixels
    w = np.where(errframe <= 0.0)
    if w[0].size != 0:
        msgs.warn("The error frame is negative for {0:d} pixels".format(w[0].size)+msgs.newline()+"Are you sure the bias frame is correct?")
        msgs.info("Masking these {0:d} pixels".format(w[0].size))
        errframe[w]  = 0.0
        slf._bpix[w] = 1.0
    w = np.where(errframe > 0.0)
    snframe = np.zeros_like(sciframe)
    snframe[w] = sciframe[w]/np.sqrt(errframe[w])
    return snframe
'''



#def sub_overscan(rawframe, numamplifiers, datasec, oscansec, method='savgol', params=[5, 65]):
#    """
#    Subtract overscan
#
#    Args:
#        frame (:obj:`numpy.ndarray`):
#            Frame from which to subtract overscan
#        numamplifiers (int):
#            Number of amplifiers for this detector.
#        datasec (list):
#            Specifies the data sections, one sub-list per amplifier
#        oscansec (list):
#            Specifies the overscan sections, one sub-list per amplifier
#        method (:obj:`str`, optional):
#            The method used to fit the overscan region.  Options are
#            polynomial, savgol, median.
#        params (:obj:`list`, optional):
#            Parameters for the overscan subtraction.  For
#            method=polynomial, set params = order, number of pixels,
#            number of repeats ; for method=savgol, set params = order,
#            window size ; for method=median, params are ignored.
#    Returns:
#        :obj:`numpy.ndarray`: The input frame with the overscan region
#        subtracted
#    """
#    for i in range(numamplifiers):
#        # Determine the section of the chip that contains data
#        dx0, dx1 = datasec[i][0][0], datasec[i][0][1]
#        dy0, dy1 = datasec[i][1][0], datasec[i][1][1]
#        if dx0 < 0: dx0 += rawframe.shape[0]
#        if dx1 <= 0: dx1 += rawframe.shape[0]
#        if dy0 < 0: dy0 += rawframe.shape[1]
#        if dy1 <= 0: dy1 += rawframe.shape[1]
#        xds = np.arange(dx0, dx1)
#        yds = np.arange(dy0, dy1)
#
#        # Determine the section of the chip that contains the overscan
#        # region
#        ox0, ox1 = oscansec[i][0][0], oscansec[i][0][1]
#        oy0, oy1 = oscansec[i][1][0], oscansec[i][1][1]
#        if ox0 < 0: ox0 += rawframe.shape[0]
#        if ox1 <= 0: ox1 += min(rawframe.shape[0], dx1)  # Truncate to datasec
#        if oy0 < 0: oy0 += rawframe.shape[1]
#        if oy1 <= 0: oy1 += min(rawframe.shape[1], dy1)  # Truncate to datasec
#        xos = np.arange(ox0, ox1)
#        yos = np.arange(oy0, oy1)
#        w = np.ix_(xos, yos)
#        oscan = rawframe[w]
#
#        # Make sure the overscan section has at least one side consistent with datasec
#        if dx1-dx0 == ox1-ox0:
#            osfit = np.median(oscan, axis=1)  # Mean was hit by CRs
#        elif dy1-dy0 == oy1-oy0:
#            osfit = np.median(oscan, axis=0)
#        elif method.lower() == 'median':
#            osfit = np.median(oscan)
#        else:
#            msgs.error('Overscan sections do not match amplifier sections for'
#                       'amplifier {0}'.format(i+1))
#
#        # Fit/Model the overscan region
#        if method.lower() == 'polynomial':
#            c = np.polyfit(np.arange(osfit.size), osfit, params[0])
#            ossub = np.polyval(c, np.arange(osfit.size))#.reshape(osfit.size,1)
#        elif method.lower() == 'savgol':
#            ossub = signal.savgol_filter(osfit, params[1], params[0])
#        elif method.lower() == 'median':  # One simple value
#            ossub = osfit * np.ones(1)
#        else:
#            # TODO: Should we raise an exception instead?
#            msgs.warn('Unknown overscan subtraction method: {0}'.format(method))
#            msgs.info('Using a linear fit to the overscan region')
#            c = np.polyfit(np.arange(osfit.size), osfit, 1)
#            ossub = np.polyval(c, np.arange(osfit.size))#.reshape(osfit.size,1)
#
#        # Determine the section of the chip that contains data for this amplifier
#        if i==0:
#            frame = rawframe.copy()
#        wd = np.ix_(xds, yds)
#        ossub = ossub.reshape(osfit.size, 1)
#        if wd[0].shape[0] == ossub.shape[0]:
#            frame[wd] -= ossub
#        elif wd[1].shape[1] == ossub.shape[0]:
#            frame[wd] -= ossub.T
#        elif method.lower() == 'median':
#            frame[wd] -= osfit
#        else:
#            msgs.error("Could not subtract bias from overscan region --"
#                       + msgs.newline() + "size of extracted regions does not match")
#
#    # Return
#    del xds, yds, xos, yos, oscan
#    return frame


#def trim(frame, numamplifiers, datasec):
#    """ Core method to trim an input image
#
#    Parameters
#    ----------
#    frame : ndarray
#    numamplifiers : int
#    datasec : list of datasecs
#      One per amplifier
#
#    Returns
#    -------
#    frame : ndarray
#      Trimmed
#    """
#    for i in range(numamplifiers):
#        #datasec = "datasec{0:02d}".format(i+1)
#        #x0, x1 = settings.spect[dnum][datasec][0][0], settings.spect[dnum][datasec][0][1]
#        #y0, y1 = settings.spect[dnum][datasec][1][0], settings.spect[dnum][datasec][1][1]
#        x0, x1 = datasec[i][0][0], datasec[i][0][1]
#        y0, y1 = datasec[i][1][0], datasec[i][1][1]
#        # Fuss with edges
#        if x0 < 0:
#            x0 += frame.shape[0]
#        if x1 <= 0:
#            x1 += frame.shape[0]
#        if y0 < 0:
#            y0 += frame.shape[1]
#        if y1 <= 0:
#            y1 += frame.shape[1]
#        if i == 0:
#            xv = np.arange(x0, x1)
#            yv = np.arange(y0, y1)
#        else:
#            xv = np.unique(np.append(xv, np.arange(x0, x1)))
#            yv = np.unique(np.append(yv, np.arange(y0, y1)))
#    # Construct and array with the rows and columns to be extracted
#    w = np.ix_(xv, yv)
##	if len(file.shape) == 2:
##		trimfile = file[w]
##	elif len(file.shape) == 3:
##		trimfile = np.zeros((w[0].shape[0],w[1].shape[1],file.shape[2]))
##		for f in range(file.shape[2]):
##			trimfile[:,:,f] = file[:,:,f][w]
##	else:
##		msgs.error("Cannot trim {0:d}D frame".format(int(len(file.shape))))
#    try:
#        return frame[w]
#    except:
#        msgs.bug("Odds are datasec is set wrong. Maybe due to transpose")
#        debugger.set_trace()
#        msgs.error("Cannot trim file")
