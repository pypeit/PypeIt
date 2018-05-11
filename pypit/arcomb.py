from __future__ import (print_function, absolute_import, division, unicode_literals)

import time

import numpy as np

from pypit import msgs
from pypit import arparse as settings


def comb_frames(frames_arr, det, frametype, **kwargs):
    """ This method has been reduced to a simple wrapper to the core method.
    It will be deprecated in a future refactor

    Parameters
    ----------
    frames_arr : ndarray (3D)
      Array of frames to be combined
    frames_arr
    det : int
      Detector index
    frametype : str, optional
      What is the type of frame being combining?

    Returns
    -------
    comb_frame : ndarray

    """
    dnum = settings.get_dnum(det)
    reject = settings.argflag[frametype]['combine']['reject']
    method = settings.argflag[frametype]['combine']['method']
    satpix = settings.argflag[frametype]['combine']['satpix']
    saturation = settings.spect[dnum]['saturation']*settings.spect[dnum]['nonlinear']
    return core_comb_frames(frames_arr, frametype=frametype,
                method=method, reject=reject, satpix=satpix, saturation=saturation, **kwargs)

def core_comb_frames(frames_arr, maskvalue=1048577, printtype=None, frametype='Unknown',
                method='weightmean', reject=None, satpix='reject', saturation=None):
    """ Combine several frames

    .. todo::
        - I've just replaced the arcycomb calls, but this function
          should probably be rewritten so that it makes better use of
          np.ma.MaskedArray objects throughout.

        - Some of the replacement code still needs to be tested!

    Parameters
    ----------
    frames_arr : ndarray (3D)
      Array of frames to be combined
    weights : str, or None (optional)
      How should the frame combination by weighted (not currently
      implemented)
    frametype : str, optional
      What is the type of frame being combining?
    maskvalue : int (optional)
      What should the masked values be set to (should be greater than
      the detector's saturation value -- Default = 1 + 2**20)
    printtype : str (optional)
      The frame type string that should be printed by armsgs. If None,
      frametype will be used
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
    if reject is None:
        reject = {'cosmics': 20., 'lowhigh': [0,0], 'level': [3.,3.], 'replace': 'maxnonsat'}

    ###########
    # FIRST DO SOME CHECKS ON THE INPUT
    ###########
    # Was printtype specified
    if printtype is None:
        printtype = frametype
    # Check the number of frames
    if frames_arr is None:
        msgs.error("No '{0:s}' frames were given to comb_frames to combine".format(printtype))
    (sz_x, sz_y, num_frames) = np.shape(frames_arr)
    if num_frames == 1:
        msgs.info("Only one frame to combine!")
        msgs.info("Returning input frame")
        return frames_arr[:, :, 0]
    else:
        msgs.info("Combining {0:d} {1:s} frames".format(num_frames, printtype))

    # Check if the user has allowed the combination of long and short
    # frames (e.g. different exposure times)
    msgs.work("lscomb feature has not been included here yet...")
    # Check the user hasn't requested to reject more frames than available
    if reject['lowhigh'][0] > 0 and reject['lowhigh'][1] > 0 and \
        reject['lowhigh'][0] + reject['lowhigh'][1] >= num_frames:
        msgs.error('You cannot reject more frames than is available with \'reject lowhigh\'.'
                   + msgs.newline() + 'There are {0:d} frames '.format(num_frames)
                   + 'and reject lowhigh will reject {0:d} low '.format(reject['lowhigh'][0])
                   + 'and {0:d} high'.format(reject['lowhigh'][1]))


    # Check that some information on the frames was supplied
    #if settings.spect is None:
    #    msgs.error('When combining the {0:s} frames, spectrograph information'.format(printtype)
    #               + msgs.newline() + 'was not provided.')
    # Calculate the values to be used if all frames are rejected in some pixels
    if reject['replace'] == 'min':
#        allrej_arr = arcycomb.minmax(frames_arr, 0)
        allrej_arr = np.amin(frames_arr, axis=2)
    elif reject['replace'] == 'max':
#        allrej_arr = arcycomb.minmax(frames_arr, 1)
        allrej_arr = np.amax(frames_arr, axis=2)
    elif reject['replace'] == 'mean':
#        allrej_arr = arcycomb.mean(frames_arr)
        allrej_arr = np.mean(frames_arr, axis=2)
    elif reject['replace'] == 'median':
#        allrej_arr = arcycomb.median(frames_arr)
        allrej_arr = np.median(frames_arr, axis=2)
    elif reject['replace'] == 'weightmean':
        msgs.work("No weights are implemented yet")
#        print('calling masked_weightmean')
#        _allrej_arr = frames_arr.copy()
#        t = time.clock()
#        _allrej_arr = arcycomb.masked_weightmean(_allrej_arr, maskvalue)
#        print('Old masked_weightmean: {0} seconds'.format(time.clock() - t))
        __allrej_arr = frames_arr.copy()
#        t = time.clock()
        __allrej_arr = new_masked_weightmean(__allrej_arr, maskvalue)
#        print('New masked_weightmean: {0} seconds'.format(time.clock() - t))
#        print(__allrej_arr.shape)
#        assert np.sum(__allrej_arr != _allrej_arr) == 0, \
#                    'Difference between old and new masked_weightmean'
        allrej_arr = __allrej_arr
##        allrej_arr = arcycomb.masked_weightmean(frames_arr, maskvalue)
#        allrej_arr = new_masked_weightmean(frames_arr, maskvalue)
    elif reject['replace'] == 'maxnonsat':
#        print('calling maxnonsat')
#        _allrej_arr = frames_arr.copy()
#        t = time.clock()
#        _allrej_arr = arcycomb.maxnonsat(_allrej_arr, saturation)
#        print('Old maxnonsat: {0} seconds'.format(time.clock() - t))
        __allrej_arr = frames_arr.copy()
#        t = time.clock()
        __allrej_arr = new_maxnonsat(__allrej_arr, saturation)
#        print('New maxnonsat: {0} seconds'.format(time.clock() - t))
#       Bug in arcycomb.maxnonsat ?
#        assert np.sum(__allrej_arr != _allrej_arr) == 0, 'Difference between old and new maxnonsat'
        allrej_arr = __allrej_arr
##        allrej_arr = arcycomb.maxnonsat(frames_arr, settings.spect[dnum]['saturation']*settings.spect[dnum]['nonlinear'])
#        allrej_arr = new_maxnonsat(frames_arr, saturation)
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
#        setsat = arcycomb.masked_limitget(frames_arr, settings.spect[dnum]['saturation']*settings.spect[dnum]['nonlinear'], 2)
        setsat = np.zeros_like(frames_arr)
        setsat[frames_arr > saturation] = 1
    elif satpix == 'reject':
        # Ignore saturated pixels in frames if possible
#        frames_arr = arcycomb.masked_limitset(frames_arr, settings.spect[dnum]['saturation']*settings.spect[dnum]['nonlinear'], 2, maskvalue)
        frames_arr[frames_arr > saturation] = maskvalue
    elif satpix == 'nothing':
        # Don't do anything special for saturated pixels (Hopefully the
        # user has specified how to deal with them below!)
        pass
    else:
        msgs.error('Option \'{0}\' '.format(satpix)
                   + 'for dealing with saturated pixels was not recognised.')
    # Delete unecessary arrays
    # None!
    ################
    # Cosmic Rays
    if reject['cosmics'] > 0.0:
        msgs.info("Rejecting cosmic rays")  # Use a robust statistic
        masked_fa = np.ma.MaskedArray(frames_arr, mask=frames_arr==maskvalue)
        medarr = np.ma.median(masked_fa, axis=2)
        stdarr = 1.4826*np.ma.median(np.ma.absolute(masked_fa - medarr[:,:,None]), axis=2)
        indx = (frames_arr != maskvalue) \
                    & (frames_arr > (medarr.data + reject['cosmics']*stdarr.data)[:,:,None])
        frames_arr[indx] = maskvalue
#        medarr = arcycomb.masked_median(frames_arr, maskvalue)
#        stdarr = 1.4826*arcycomb.masked_median(np.abs(frames_arr-medarr[:, :, np.newaxis]), maskvalue)
#        frames_arr = arcycomb.masked_limitsetarr(frames_arr, (medarr + reject['cosmics']*stdarr), 2, maskvalue)
        # Delete unecessary arrays
        del medarr, stdarr
    else:
        msgs.info("Not rejecting cosmic rays")
    ################
    # Low and High pixel rejection --- Masks *additional* pixels
    rejlo, rejhi = reject['lowhigh']
    if reject['lowhigh'][0] > 0 or reject['lowhigh'][1] > 0:
        # First reject low pixels
        frames_arr = np.sort(frames_arr, axis=2)
        if reject['lowhigh'][0] > 0:
            msgs.info("Rejecting {0:d} deviant low pixels".format(reject['lowhigh'][0]))
            while rejlo > 0:
                xi, yi = np.indices(sz_x, sz_y)
                frames_arr[xi, yi, np.argmin(frames_arr, axis=2)] = maskvalue
                del xi, yi
                rejlo -= 1
        # Now reject high pixels
        if reject['lowhigh'][1] > 0:
            msgs.info("Rejecting {0:d} deviant high pixels".format(reject['lowhigh'][1]))
            frames_arr[np.where(frames_arr == maskvalue)] *= -1
            while rejhi > 0:
                xi, yi = np.indices(sz_x, sz_y)
                frames_arr[xi, yi, np.argmax(frames_arr, axis=2)] = -maskvalue
                del xi, yi
                rejhi -= 1
            frames_arr[np.where(frames_arr) == -maskvalue] *= -1
# The following is an example of *not* masking additional pixels
#		if reject['lowhigh'][1] > 0:
#			msgs.info("Rejecting {0:d} deviant high pixels".format(reject['lowhigh'][1]))
#			masktemp[:,:,-reject['lowhigh'][0]:] = True
    else:
        msgs.info("Not rejecting any low/high pixels")
    ################
    # Deviant Pixels
    if reject['level'][0] > 0.0 or reject['level'][1] > 0.0:
        msgs.info("Rejecting deviant pixels")  # Use a robust statistic

        masked_fa = np.ma.MaskedArray(frames_arr, mask=frames_arr==maskvalue)
        medarr = np.ma.median(masked_fa, axis=2)
        stdarr = 1.4826*np.ma.median(np.ma.absolute(masked_fa - medarr[:,:,None]), axis=2)
        indx = (frames_arr != maskvalue) \
                    & ( (frames_arr > (medarr.data + reject['cosmics']*stdarr.data)[:,:,None])
                        | (frames_arr < (medarr.data - reject['cosmics']*stdarr.data)[:,:,None]))
        frames_arr[indx] = maskvalue

#        medarr = arcycomb.masked_median(frames_arr, maskvalue)
#        stdarr = 1.4826*arcycomb.masked_median(np.abs(frames_arr-medarr[:, :, np.newaxis]), maskvalue)
#        frames_arr = arcycomb.masked_limitsetarr(frames_arr, (medarr - reject['level'][0]*stdarr), -2, maskvalue)
#        frames_arr = arcycomb.masked_limitsetarr(frames_arr, (medarr + reject['level'][1]*stdarr), 2, maskvalue)
        # Delete unecessary arrays
        del medarr, stdarr
    else:
        msgs.info("Not rejecting deviant pixels")

    ##############
    # Combine the arrays
    msgs.info("Combining frames with a {0:s} operation".format(method))
    if method == 'mean':
#        frames_arr = arcycomb.masked_mean(frames_arr, maskvalue)
        comb_frame = np.ma.mean(np.ma.MaskedArray(frames_arr, mask=frames_arr==maskvalue), axis=2)
    elif method == 'median':
#        frames_arr = arcycomb.masked_median(frames_arr, maskvalue)
        comb_frame = np.ma.median(np.ma.MaskedArray(frames_arr, mask=frames_arr==maskvalue), axis=2)
    elif method == 'weightmean':
#        print('calling masked_weightmean')
#        _frames_arr = frames_arr.copy()
#        t = time.clock()
#        _frames_arr = arcycomb.masked_weightmean(_frames_arr, maskvalue)
#        print('Old masked_weightmean: {0} seconds'.format(time.clock() - t))
        __frames_arr = frames_arr.copy()
#        t = time.clock()
        __frames_arr = new_masked_weightmean(__frames_arr, maskvalue)
#        print('New masked_weightmean: {0} seconds'.format(time.clock() - t))
#        print(__frames_arr.shape)
#
#        if np.sum(np.absolute(__frames_arr-_frames_arr) > 1e-10) != 0:
#            print(np.sum(np.absolute(__frames_arr-_frames_arr) > 1e-10))
#            plt.imshow(_frames_arr, origin='lower', interpolation='nearest', aspect='auto')
#            plt.show()
#            plt.imshow(__frames_arr, origin='lower', interpolation='nearest', aspect='auto')
#            plt.show()
#            plt.imshow(_frames_arr - __frames_arr, origin='lower', interpolation='nearest', aspect='auto')
#            plt.show()
#            plt.imshow(np.ma.divide(_frames_arr,__frames_arr) - 1, origin='lower', interpolation='nearest', aspect='auto')
#            plt.colorbar()
#            plt.show()
#                    
#        assert np.sum( np.absolute(__frames_arr-_frames_arr) > 1e-10 ) == 0, \
#                    'Difference between old and new masked_weightmean'
        comb_frame = __frames_arr
    else:
        msgs.error("Combination type '{0:s}' is unknown".format(method))
    ##############
    # If any pixels are completely masked, apply user-specified function
    msgs.info("Replacing completely masked pixels with the {0:s} value of the input frames".format(reject['replace']))
#    frames_arr = arcycomb.masked_replace(frames_arr, allrej_arr, maskvalue)
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
    msgs.info("{0:d} {1:s} frames combined successfully!".format(num_frames, printtype))
    # Make sure the returned array is the correct type
    comb_frame = np.array(comb_frame, dtype=np.float)
    return comb_frame


def new_masked_weightmean(a, maskvalue):
    num = np.ma.MaskedArray(a.copy(), mask=(a==maskvalue))
    num[np.invert(num.mask) & (num <= 1.0)] = 0.0
    num = np.ma.sum(np.ma.sqrt(num)*num, axis=2)
    den = np.ma.MaskedArray(a.copy(), mask=(a==maskvalue))
    den[np.invert(den.mask) & (den <= 1.0)] = 1.0
    den = np.ma.sum(np.sqrt(den), axis=2)
    return np.ma.divide(num, den).filled(maskvalue)


def new_maxnonsat(array, saturated):
    """
    sz_x, sz_y, nfr = array.shape
    mmarr = np.zeros((sz_x,sz_y), dtype=float)
    d = 0
    for x in range(sz_x):
        for y in range(sz_y):
            # Sort the array
            temp = 0.0
            minv = saturated
            for n in range(nfr):
                if array[x,y,n] > temp and (BUG: temp) < saturated:
                    temp = array[x,y,n]
                if array[x,y,n] < minv:
                    minv = array[x,y,n]
            if temp == 0.0:
                mmarr[x,y] = minv
            else:
                mmarr[x,y] = temp
    return mmarr
    """
    minimum = np.amin(np.clip(array, None, saturated), axis=2)
    _array = np.ma.MaskedArray(array, mask=np.invert((array > 0.0) & (array<saturated)))
    maximum = np.ma.amax(_array, axis=2)
    maximum[maximum.mask] = minimum[maximum.mask]
    return maximum.data

