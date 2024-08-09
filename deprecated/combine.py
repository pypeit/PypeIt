
# TODO JFH This code is a giant piece of crapp that should be rewritten from scratch
def comb_frames(frames_arr, saturation=None,
                     maskvalue=1048577, method='weightmean', satpix='reject', cosmics=None,
                     n_lohi=[0,0], sig_lohi=[3.,3.], replace='maxnonsat'):
    """
    Combine several frames

    .. todo::
        - Make better use of np.ma.MaskedArray objects throughout?
        - More testing of replacement code necessary?
        - Improve docstring...

    Parameters
    ----------
    frames_arr : ndarray (3D)
      Array of frames to be combined
    weights : str, or None (optional)
      How should the frame combination by weighted (not currently
      implemented)
    maskvalue : int (optional)
      What should the masked values be set to (should be greater than
      the detector's saturation value -- Default = 1 + 2**20)
    reject : dict, optional
      Set the rejection parameters:  cosmics, lowhigh, level, replace
      Perhaps these should be called out separately
    satpix : str, optional
      Method for handling saturated pixels
    saturation : float, optional
      Saturation value;  only required for some choices of reject['replace']

    Returns
    -------
    comb_frame : ndarray
    """
    ###########
    # FIRST DO SOME CHECKS ON THE INPUT
    ###########
    # Was printtype specified
    if frames_arr is None:
        msgs.error("No frames were given to comb_frames to combine")
    (sz_x, sz_y, num_frames) = np.shape(frames_arr)
    if num_frames == 1:
        msgs.info("Only one frame to combine!")
        msgs.info("Returning input frame")
        return frames_arr[:, :, 0]
    else:
        msgs.info("Combining {0:d} frames".format(num_frames))

    # Check if the user has allowed the combination of long and short
    # frames (e.g. different exposure times)
    msgs.work("lscomb feature has not been included here yet...")
    # Check the user hasn't requested to reject more frames than available
    if n_lohi[0] > 0 and n_lohi[1] > 0 and n_lohi[0] + n_lohi[1] >= num_frames:
        msgs.error('You cannot reject more frames than are available with \'n_lohi\'.'
                   + msgs.newline() + 'There are {0:d} frames '.format(num_frames)
                   + 'and n_lohi will reject {0:d} low and {1:d} high values.'.format(
                                                                n_lohi[0], n_lohi[1]))

    # Calculate the values to be used if all frames are rejected in some pixels
    if replace == 'min':
        allrej_arr = np.amin(frames_arr, axis=2)
    elif replace == 'max':
        allrej_arr = np.amax(frames_arr, axis=2)
    elif replace == 'mean':
        allrej_arr = np.mean(frames_arr, axis=2)
    elif replace == 'median':
        allrej_arr = np.median(frames_arr, axis=2)
    elif replace == 'weightmean':
        msgs.work("No weights are implemented yet")
        allrej_arr = frames_arr.copy()
        allrej_arr = masked_weightmean(allrej_arr, maskvalue)
    elif replace == 'maxnonsat':
        allrej_arr = frames_arr.copy()
        allrej_arr = maxnonsat(allrej_arr, saturation)
    else:
        msgs.error("You must specify what to do in case all pixels are rejected")

    ################
    # Saturated Pixels
    msgs.info("Finding saturated and non-linear pixels")
    if satpix == 'force':
        # If a saturated pixel is in one of the frames, force them to
        # all have saturated pixels
#		satw = np.zeros_like(frames_arr)
#		satw[np.where(frames_arr > settings.spect['det']['saturation']*settings.spect['det']['nonlinear'])] = 1.0
#		satw = np.any(satw,axis=2)
#		del satw
        setsat = np.zeros_like(frames_arr)
        setsat[frames_arr > saturation] = 1
    elif satpix == 'reject':
        # Ignore saturated pixels in frames if possible
        frames_arr[frames_arr > saturation] = maskvalue
    elif satpix == 'nothing':
        # Don't do anything special for saturated pixels (Hopefully the
        # user has specified how to deal with them below!)
        pass
    else:
        msgs.error('Option \'{0}\' '.format(satpix)
                   + 'for dealing with saturated pixels was not recognised.')

    ################
    # Cosmic Rays
    if cosmics > 0.0:
        msgs.info("Rejecting cosmic rays")  # Use a robust statistic
        masked_fa = np.ma.MaskedArray(frames_arr, mask=frames_arr==maskvalue)
        medarr = np.ma.median(masked_fa, axis=2)
        stdarr = 1.4826*np.ma.median(np.ma.absolute(masked_fa - medarr[:,:,None]), axis=2)
        indx = (frames_arr != maskvalue) \
                    & (frames_arr > (medarr.data + cosmics * stdarr.data)[:,:,None])
        frames_arr[indx] = maskvalue
        # Delete unecessary arrays
        del medarr, stdarr
    else:
        msgs.info("Not rejecting cosmic rays")

    ################
    # Low and High pixel rejection --- Masks *additional* pixels
    rejlo, rejhi = n_lohi
    if n_lohi[0] > 0 or n_lohi[1] > 0:

        # First reject low pixels
        frames_arr = np.sort(frames_arr, axis=2)
        if n_lohi[0] > 0:
            msgs.info("Rejecting {0:d} deviant low pixels".format(n_lohi[0]))
            while rejlo > 0:
                xi, yi = np.indices(sz_x, sz_y)
                frames_arr[xi, yi, np.argmin(frames_arr, axis=2)] = maskvalue
                del xi, yi
                rejlo -= 1

        # Now reject high pixels
        if n_lohi[1] > 0:
            msgs.info("Rejecting {0:d} deviant high pixels".format(n_lohi[1]))
            frames_arr[np.where(frames_arr == maskvalue)] *= -1
            while rejhi > 0:
                xi, yi = np.indices(sz_x, sz_y)
                frames_arr[xi, yi, np.argmax(frames_arr, axis=2)] = -maskvalue
                del xi, yi
                rejhi -= 1
            frames_arr[np.where(frames_arr) == -maskvalue] *= -1

# TODO: Do we need this?
# The following is an example of *not* masking additional pixels
#		if reject['lowhigh'][1] > 0:
#			msgs.info("Rejecting {0:d} deviant high pixels".format(reject['lowhigh'][1]))
#			masktemp[:,:,-reject['lowhigh'][0]:] = True
    else:
        msgs.info("Not rejecting any low/high pixels")

    ################
    # Deviant Pixels
    # TODO: sig_lohi (what was level) is not actually used, instead this
    # just selects if cosmics should be used.  Is this intentional?  Why
    # not just do: `if cosmics > 0:`?
    if sig_lohi[0] > 0.0 or sig_lohi[1] > 0.0:
        msgs.info("Rejecting deviant pixels")  # Use a robust statistic

        masked_fa = np.ma.MaskedArray(frames_arr, mask=frames_arr==maskvalue)
        medarr = np.ma.median(masked_fa, axis=2)
        stdarr = 1.4826*np.ma.median(np.ma.absolute(masked_fa - medarr[:,:,None]), axis=2)
        indx = (frames_arr != maskvalue) \
                    & ( (frames_arr > (medarr.data + cosmics*stdarr.data)[:,:,None])
                        | (frames_arr < (medarr.data - cosmics*stdarr.data)[:,:,None]))
        frames_arr[indx] = maskvalue

        # Delete unecessary arrays
        del medarr, stdarr
    else:
        msgs.info("Not rejecting deviant pixels")

    ##############
    # Combine the arrays
    msgs.info("Combining frames with a {0:s} operation".format(method))
    if method == 'mean':
        comb_frame = np.ma.mean(np.ma.MaskedArray(frames_arr, mask=frames_arr==maskvalue), axis=2)
    elif method == 'median':
        comb_frame = np.ma.median(np.ma.MaskedArray(frames_arr, mask=frames_arr==maskvalue), axis=2)
    elif method == 'weightmean':
        comb_frame = frames_arr.copy()
        comb_frame = masked_weightmean(comb_frame, maskvalue)
    else:
        msgs.error("Combination type '{0:s}' is unknown".format(method))

    ##############
    # If any pixels are completely masked, apply user-specified function
    msgs.info("Replacing completely masked pixels with the {0:s} value of the input frames".format(replace))
    indx = comb_frame == maskvalue
    comb_frame[indx] = allrej_arr[indx]
    # Delete unecessary arrays
    del allrej_arr

    ##############
    # Apply the saturated pixels:
    if satpix == 'force':
        msgs.info("Applying saturated pixels to final combined image")
        comb_frame[setsat] = saturation # settings.spect[dnum]['saturation']

    ##############
    # And return a 2D numpy array
    msgs.info("{0:d} frames combined successfully!".format(num_frames))
    # Make sure the returned array is the correct type
    comb_frame = np.array(comb_frame, dtype=np.float)
    return comb_frame


def masked_weightmean(a, maskvalue):
    """
    .. todo::
        Document this!
    """
    num = np.ma.MaskedArray(a.copy(), mask=(a==maskvalue))
    num[np.invert(num.mask) & (num <= 1.0)] = 0.0
    num = np.ma.sum(np.ma.sqrt(num)*num, axis=2)
    den = np.ma.MaskedArray(a.copy(), mask=(a==maskvalue))
    den[np.invert(den.mask) & (den <= 1.0)] = 1.0
    den = np.ma.sum(np.sqrt(den), axis=2)
    return np.ma.divide(num, den).filled(maskvalue)


def maxnonsat(array, saturated):
    """
    .. todo::
        Document this!
    """
    minimum = np.amin(np.clip(array, None, saturated), axis=2)
    _array = np.ma.MaskedArray(array, mask=np.invert((array > 0.0) & (array<saturated)))
    maximum = np.ma.amax(_array, axis=2)
    maximum[maximum.mask] = minimum[maximum.mask]
    return maximum.data


