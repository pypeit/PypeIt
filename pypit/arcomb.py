from __future__ import absolute_import, division, print_function

import numpy as np
from pypit import armsgs
from pypit import arparse as settings

# Logging
msgs = armsgs.get_logger()


def comb_frames(frames_arr, det, frametype, weights=None, maskvalue=1048577, printtype=None):
    """ Combine several frames

    Parameters
    ----------
    frames_arr : ndarray (3D)
      Array of frames to be combined
    det : int
      Detector index
    frametype : str
      What is the type of frame being combining? (only used for screen printout)
    weights : str, or None (optional)
      How should the frame combination by weighted (not currently implemented)
    maskvalue : int (optional)
      What should the masked values be set to (should be greater than the detector's saturation value -- Default = 1 + 2**20)
    printtype : str (optional)
      The frame type string that should be printed by armsgs. If None, frametype will be used
    """

    from pypit import arcycomb
    dnum = settings.get_dnum(det)
    reject = settings.argflag[frametype]['combine']['reject']
    method = settings.argflag[frametype]['combine']['method']
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
    # Check if the user has allowed the combination of long and short frames (e.g. different exposure times)
    msgs.work("lscomb feature has not been included here yet...")
    # Check the user hasn't requested to reject more frames than available
    if reject['lowhigh'][0] > 0 and reject['lowhigh'][1] > 0 and \
        reject['lowhigh'][0] + reject['lowhigh'][1] >= num_frames:
        msgs.error("You cannot reject more frames than is available with 'reject lowhigh'." + msgs.newline() +
                   "There are {0:d} frames and reject lowhigh will reject {1:d} low and {2:d} high".format(num_frames, reject['lowhigh'][0], reject['lowhigh'][1]))
    # Check that some information on the frames was supplied
    if settings.spect is None:
        msgs.error("When combining the {0:s} frames, spectrograph information".format(printtype)+msgs.newline()+"was not provided")
    # Calculate the values to be used if all frames are rejected in some pixels
    if reject['replace'] == 'min':
        allrej_arr = arcycomb.minmax(frames_arr, 0)
    elif reject['replace'] == 'max':
        allrej_arr = arcycomb.minmax(frames_arr, 1)
    elif reject['replace'] == 'mean':
        allrej_arr = arcycomb.mean(frames_arr)
    elif reject['replace'] == 'median':
        allrej_arr = arcycomb.median(frames_arr)
    elif reject['replace'] == 'weightmean':
        msgs.work("No weights are implemented yet")
        allrej_arr = arcycomb.masked_weightmean(frames_arr, maskvalue)
    elif reject['replace'] == 'maxnonsat':
        allrej_arr = arcycomb.maxnonsat(frames_arr, settings.spect[dnum]['saturation']*settings.spect[dnum]['nonlinear'])
    else:
        msgs.error("You must specify what to do in case all pixels are rejected")
    ################
    # Saturated Pixels
    msgs.info("Finding saturated and non-linear pixels")
    if settings.argflag[frametype]['combine']['satpix'] == 'force':
        # If a saturated pixel is in one of the frames, force them to all have saturated pixels
#		satw = np.zeros_like(frames_arr)
#		satw[np.where(frames_arr > settings.spect['det']['saturation']*settings.spect['det']['nonlinear'])] = 1.0
#		satw = np.any(satw,axis=2)
        setsat = arcycomb.masked_limitget(frames_arr, settings.spect[dnum]['saturation']*settings.spect[dnum]['nonlinear'], 2)
#		del satw
    elif settings.argflag[frametype]['combine']['satpix'] == 'reject':
        # Ignore saturated pixels in frames if possible
        frames_arr = arcycomb.masked_limitset(frames_arr, settings.spect[dnum]['saturation']*settings.spect[dnum]['nonlinear'], 2, maskvalue)
    elif settings.argflag[frametype]['combine']['satpix'] == 'nothing':
        # Don't do anything special for saturated pixels (Hopefully the user has specified how to deal with them below!)
        pass
    else:
        msgs.error("Option '{0:s}' for dealing with saturated pixels was not recognised".format(settings.argflag[frametype]['combine']['satpix']))
    # Delete unecessary arrays
    # None!
    ################
    # Cosmic Rays
    if reject['cosmics'] > 0.0:
        msgs.info("Rejecting cosmic rays")  # Use a robust statistic
        medarr = arcycomb.masked_median(frames_arr, maskvalue)
        stdarr = 1.4826*arcycomb.masked_median(np.abs(frames_arr-medarr[:, :, np.newaxis]), maskvalue)
        frames_arr = arcycomb.masked_limitsetarr(frames_arr, (medarr + reject['cosmics']*stdarr), 2, maskvalue)
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
        medarr = arcycomb.masked_median(frames_arr, maskvalue)
        stdarr = 1.4826*arcycomb.masked_median(np.abs(frames_arr-medarr[:, :, np.newaxis]), maskvalue)
        frames_arr = arcycomb.masked_limitsetarr(frames_arr, (medarr - reject['level'][0]*stdarr), -2, maskvalue)
        frames_arr = arcycomb.masked_limitsetarr(frames_arr, (medarr + reject['level'][1]*stdarr), 2, maskvalue)
        # Delete unecessary arrays
        del medarr, stdarr
    else:
        msgs.info("Not rejecting deviant pixels")
    ##############
    # Combine the arrays
    msgs.info("Combining frames with a {0:s} operation".format(method))
    if method == 'mean':
        frames_arr = arcycomb.masked_mean(frames_arr, maskvalue)
    elif method == 'median':
        frames_arr = arcycomb.masked_median(frames_arr, maskvalue)
    elif method == 'weightmean':
        frames_arr = arcycomb.masked_weightmean(frames_arr, maskvalue)
    else:
        msgs.error("Combination type '{0:s}' is unknown".format(method))
    ##############
    # If any pixels are completely masked, apply user-specified function
    msgs.info("Replacing completely masked pixels with the {0:s} value of the input frames".format(reject['replace']))
    frames_arr = arcycomb.masked_replace(frames_arr, allrej_arr, maskvalue)
    # Delete unecessary arrays
    del allrej_arr
    ##############
    # Apply the saturated pixels:
    if settings.argflag[frametype]['combine']['satpix'] == 'force':
        msgs.info("Applying saturated pixels to final combined image")
        frames_arr[setsat] = settings.spect[dnum]['saturation']
    ##############
    # And return a 2D numpy array
    msgs.info("{0:d} {1:s} frames combined successfully!".format(num_frames, printtype))
    # Make sure the returned array is the correct type
    frames_arr = np.array(frames_arr, dtype=np.float)
    return frames_arr
