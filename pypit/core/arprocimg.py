""" Module for image processing core methods
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import time

import numpy as np
import os

from scipy import signal, ndimage

from pypit import msgs

from pypit.core import arlris
from pypit.core import ardeimos
from pypit import arutils
from pypit import arparse
from pypit import arpixels

from pypit import ardebug as debugger

try:
    basestring
except NameError:
    basestring = str

def badpix(frame, numamplifiers, datasec, sigdev=10.0):
    """
    frame is a master bias frame
    numamplifiers : int
    datasec : list
    sigdev is the number of standard deviations away from the median that a pixel needs to be in order to be classified as a bad pixel
    """
    bpix = np.zeros_like(frame, dtype=np.int)
    subfr, tframe, temp = None, None, None
    #for i in range(settings.spect[dnum]['numamplifiers']):
    for i in range(numamplifiers):
        #datasec = "datasec{0:02d}".format(i+1)
        x0, x1 = datasec[i][0][0], datasec[i][0][1]
        y0, y1 = datasec[i][1][0], datasec[i][1][1]
        #x0, x1 = settings.spect[dnum][datasec][0][0], settings.spect[dnum][datasec][0][1]
        #y0, y1 = settings.spect[dnum][datasec][1][0], settings.spect[dnum][datasec][1][1]
        xv = np.arange(x0, x1)
        yv = np.arange(y0, y1)
        # Construct an array with the rows and columns to be extracted
        w = np.ix_(xv,yv)
        tframe = frame[w]
        temp = np.abs(np.median(tframe)-tframe)
        sigval = max(np.median(temp)*1.4826, 1.4826)
        ws = np.where(temp > sigdev*sigval)
        subfr = np.zeros(tframe.shape, dtype=np.int)
        subfr[ws] = 1
        bpix[w] = subfr
    del subfr, tframe, temp
    # Finally, trim the bad pixel frame
    bpix = trim(bpix, numamplifiers, datasec)
    msgs.info("Identified {0:d} bad pixels".format(int(np.sum(bpix))))
    return bpix


def bias_subtract(rawframe, msbias, numamplifiers=None, datasec=None, oscansec=None):
    """ Core routine for bias subtraction
    Calls sub_overscan if msbias == 'overscan'

    Parameters
    ----------
    rawframe : ndarray
    msbias : ndarray or str

    Returns
    -------
    newframe : ndarray
      Bias subtracted frame

    """
    if type(msbias) is np.ndarray:
        msgs.info("Subtracting bias image from raw frame")
        newframe = rawframe-msbias  # Subtract the master bias frame
    elif isinstance(msbias, basestring):
        if msbias == "overscan":
            newframe = sub_overscan(rawframe, numamplifiers, datasec, oscansec, settings=None)
        else:
            msgs.error("Could not subtract bias level with the input bias approach") #when loading {0:s} frames".format(frametype))
    return newframe


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


def get_datasec_trimmed(spectrograph, scifile, det, settings_det,
                        naxis0=None, naxis1=None):
    """
    Primarily a wrapper with calls to get_datasec and pix_to_amp()

    Parameters
    ----------
    fitstbl
    det : int
    scidx : int

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


def get_datasec(spectrograph, scifile, numamplifiers=None, det=None):
    """  Determine the data and overscan sections of an image

    Currently only used for LRIS and DEIMOS (with their multiple detectors
    packed in funny ways).  Should consider another approach.

    Parameters
    ----------
    spectrograph : str
    scifile : str
    numamplifiers : int (optional)
    det : int (optional)
      Detector number, starts at 1

    Returns
    -------
    datasec : list
    oscansec : list
    naxis0 : int
    naxis1 : int
    """
    # Get naxis0, naxis1, datasec, oscansec, ampsec for specific instruments
    datasec, oscansec, naxis0, naxis1 = [], [], 0, 0
    # TODO -- Remove instrument specific items in a method like this
    if spectrograph in ['keck_lris_blue', 'keck_lris_red']:
        msgs.info("Parsing datasec and oscansec from headers")
        temp, head0, secs = arlris.read_lris(scifile, det)
        for kk in range(numamplifiers):
            #datasec = "datasec{0:02d}".format(kk+1)
            #settings.spect[dnum][datasec] = settings.load_sections(secs[0][kk], fmt_iraf=False)
            datasec.append(arparse.load_sections(secs[0][kk], fmt_iraf=False))
            #oscansec = "oscansec{0:02d}".format(kk+1)
            #settings.spect[dnum][oscansec] = settings.load_sections(secs[1][kk], fmt_iraf=False)
            oscansec.append(arparse.load_sections(secs[1][kk], fmt_iraf=False))
    elif spectrograph in ['keck_deimos']:
        msgs.info("Parsing datasec and oscansec from headers")
        temp, head0, secs = ardeimos.read_deimos(scifile, det=det)
        datasec.append(arparse.load_sections(secs[0][0], fmt_iraf=False))
        oscansec.append(arparse.load_sections(secs[1][0], fmt_iraf=False))
    else:  # Other instruments are set in their settings file
        msgs.warn("Should not have called get_datasec!")
        return datasec, oscansec, naxis0, naxis1

    naxis0 = temp.shape[0]
    naxis1 = temp.shape[1]
    # Return
    return datasec, oscansec, naxis0, naxis1



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


def lacosmic(det, sciframe, settings_det, maxiter=1, grow=1.5,
             varframe=None, remove_compact_obj=True):
    """
    settings_det : settings.spect[dnum]
      Detector info

    Identify cosmic rays using the L.A.Cosmic algorithm
    U{http://www.astro.yale.edu/dokkum/lacosmic/}
    (article : U{http://arxiv.org/abs/astro-ph/0108003})
    This routine is mostly courtesy of Malte Tewes

    :param grow: Once CRs are identified, grow each CR detection by all pixels within this radius
    :return: mask of cosmic rays (0=no CR, 1=CR)
    """
    dnum = arparse.get_dnum(det)

    msgs.info("Detecting cosmic rays with the L.A.Cosmic algorithm")
    msgs.work("Include these parameters in the settings files to be adjusted by the user")
    sigclip = 5.0
    sigfrac = 0.3
    objlim  = 5.0
    # Set the settings
    scicopy = sciframe.copy()
    crmask = np.cast['bool'](np.zeros(sciframe.shape))
    sigcliplow = sigclip * sigfrac

    # Determine if there are saturated pixels
    satpix = np.zeros_like(sciframe)
    satlev = settings_det['saturation']*settings_det['nonlinear']
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
        #set_trace()
        subsam = arutils.subsample(scicopy)
        conved = signal.convolve2d(subsam, laplkernel, mode="same", boundary="symm")
        cliped = conved.clip(min=0.0)
        lplus = arutils.rebin(cliped, np.array(cliped.shape)/2.0)

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
        #debugger.set_trace()

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
        if ncrp == 0: break
    # Additional algorithms (not traditionally implemented by LA cosmic) to remove some false positives.
    msgs.work("The following algorithm would be better on the rectified, tilts-corrected image")
    filt  = ndimage.sobel(sciframe, axis=1, mode='constant')
    filty = ndimage.sobel(filt/np.sqrt(np.abs(sciframe)), axis=0, mode='constant')
    filty[np.where(np.isnan(filty))]=0.0

    # Old cr_screen can yield nan pixels, new one does not, meaning that
    # there are differences between the returned arrays.  For all
    # non-nan pixels, the two algorithms are identical.
#    print('calling cr_screen')
#    t = time.clock()
#    _sigimg  = arcyproc.cr_screen(filty,0.0)
#    print('Old cr_screen: {0} seconds'.format(time.clock() - t))
#    print(np.sum(np.invert(np.isfinite(_sigimg))))
#    t = time.clock()
    sigimg  = new_cr_screen(filty)
#    print('New cr_screen: {0} seconds'.format(time.clock() - t))
#    print(np.sum(np.invert(np.isfinite(sigimg))))
#    if np.sum(_sigimg != sigimg) != 0:
#        plt.imshow(_sigimg, origin='lower', interpolation='nearest', aspect='auto')
#        plt.show()
#        plt.imshow(sigimg, origin='lower', interpolation='nearest', aspect='auto')
#        plt.show()
#        plt.imshow(_sigimg-sigimg, origin='lower', interpolation='nearest', aspect='auto')
#        plt.colorbar()
#        plt.show()
#        r = np.ma.divide(_sigimg,sigimg)-1
#        print(np.ma.sum(r))
#        plt.imshow(np.ma.divide(_sigimg,sigimg)-1, origin='lower', interpolation='nearest', aspect='auto')
#        plt.colorbar()
#        plt.show()
#    assert np.sum(_sigimg != sigimg) == 0, 'Difference between old and new cr_screen'

#    print(sigimg.shape)
#    print(new_sigimg.shape)
#    print(np.ma.mean(np.ma.log10(sigimg)))
#    print(np.ma.mean(np.ma.log10(new_sigimg)))
#    print(np.ma.mean(np.ma.log10(sigimg)-np.ma.log10(new_sigimg)))
#
#    plt.imshow(np.ma.log10(sigimg), origin='lower', interpolation='nearest', aspect='auto')
#    plt.colorbar()
#    plt.show()
#
#    plt.imshow(np.ma.log10(new_sigimg), origin='lower', interpolation='nearest', aspect='auto')
#    plt.colorbar()
#    plt.show()
#
#    plt.imshow(np.ma.log10(new_sigimg)-np.ma.log10(sigimg), origin='lower', interpolation='nearest', aspect='auto')
#    plt.colorbar()
#    plt.show()

    sigsmth = ndimage.filters.gaussian_filter(sigimg,1.5)
    sigsmth[np.where(np.isnan(sigsmth))]=0.0
    sigmask = np.cast['bool'](np.zeros(sciframe.shape))
    sigmask[np.where(sigsmth>sigclip)] = True
    crmask = np.logical_and(crmask, sigmask)
    msgs.info("Growing cosmic ray mask by 1 pixel")
#    print('calling grow_masked')
#    t = time.clock()
#    _crmask = arcyutils.grow_masked(crmask.astype(np.float), grow, 1.0)
#    print('Old grow_masked: {0} seconds'.format(time.clock() - t))
#    t = time.clock()
    crmask = new_grow_masked(crmask.astype(np.float), grow, 1.0)
#    print('New grow_masked: {0} seconds'.format(time.clock() - t))
#    assert np.sum(_crmask != crmask) == 0, 'Difference between old and new grow_masked'

    return crmask


def new_cr_screen(a, mask_value=0.0, spatial_axis=1):
    r"""
    Calculate the significance of pixel deviations from the median along
    the spatial direction.

    No type checking is performed of the input array; however, the
    function assumes floating point values.

    Args:
        a (numpy.ndarray): Input 2D array
        mask_value (float): (**Optional**) Values to ignore during the
            calculation of the median.  Default is 0.0.
        spatial_axis (int): (**Optional**) Axis along which to calculate
            the median.  Default is 1.

    Returns:
        numpy.ndarray: Returns a map of :math:`|\Delta_{i,j}|/\sigma_j`,
        where :math:`\Delta_{i,j}` is the difference between the pixel
        value and the median along axis :math:`i` and :math:`\sigma_j`
        is robustly determined using the median absolute deviation,
        :math:`sigma_j = 1.4826 MAD`.
    """
    # Check input
    if len(a.shape) != 2:
        msgs.error('Input array must be two-dimensional.')
    if spatial_axis not in [0,1]:
        msgs.error('Spatial axis must be 0 or 1.')

    # Mask the pixels equal to mask value: should use np.isclose()
    _a = np.ma.MaskedArray(a, mask=(a==mask_value))
    # Get the median along the spatial axis
    meda = np.ma.median(_a, axis=spatial_axis)
    # Get a robust measure of the standard deviation using the median
    # absolute deviation; 1.4826 factor is the ratio of sigma/MAD
    d = np.absolute(_a - meda[:,None])
    mada = 1.4826*np.ma.median(d, axis=spatial_axis)
    # Return the ratio of the difference to the standard deviation
    return np.ma.divide(d, mada[:,None]).filled(mask_value)


def new_grow_masked(img, grow, growval):

    if not np.any(img == growval):
        return img

    _img = img.copy()
    sz_x, sz_y = img.shape
    d = int(1+grow)
    rsqr = grow*grow

    # Grow any masked values by the specified amount
    for x in range(sz_x):
        for y in range(sz_y):
            if img[x,y] != growval:
                continue

            mnx = 0 if x-d < 0 else x-d
            mxx = x+d+1 if x+d+1 < sz_x else sz_x
            mny = 0 if y-d < 0 else y-d
            mxy = y+d+1 if y+d+1 < sz_y else sz_y

            for i in range(mnx,mxx):
                for j in range(mny, mxy):
                    if (i-x)*(i-x)+(j-y)*(j-y) <= rsqr:
                        _img[i,j] = growval
    return _img


def gain_frame(datasec_img, namp, gain_list):
    """ Generate a gain image

    Parameters
    ----------
    datasec_img : ndarray
    namp : int
    gain_list : list

    Returns
    -------
    gain_img : ndarray

    """
    #namp = settings.spect[dnum]['numamplifiers'])
    #gains = settings.spect[dnum]['gain'][amp - 1]
    msgs.warn("Should probably be measuring the gain across the amplifier boundary")

    # Loop on amplifiers
    gain_img = np.zeros_like(datasec_img)
    for ii in range(namp):
        amp = ii+1
        amppix = datasec_img == amp
        gain_img[amppix] = gain_list[ii]
    # Return
    return gain_img


def rn_frame(det, datasec_img, settings_det):
    """ Generate a RN image

    Parameters
    ----------
    det
    settings_det : settings.spect[dnum]

    Returns
    -------
    rn_img : ndarray
      Read noise *variance* image (i.e. RN**2)
    """
    dnum = arparse.get_dnum(det)

    # Loop on amplifiers
    rnimg = np.zeros_like(datasec_img)
    for ii in range(settings_det['numamplifiers']):
        amp = ii+1
        amppix = datasec_img == amp
        rnimg[amppix] = (settings_det['ronoise'][ii]**2 +
                         (0.5*settings_det['gain'][ii])**2)
    # Return
    return rnimg


def sub_overscan(rawframe, numamplifiers, datasec, oscansec, settings=None):
    """
    Subtract overscan

    Parameters
    ----------
    frame : ndarray
      frame which should have the overscan region subtracted
    numamplifiers : int
    datasec : list
      Specifies the data sections, one sub-list per amplifier
    oscansec : list
      Specifies the overscan sections, one sub-list per amplifier
    settings : dict, optional
      Describes the internal options for the overscan subtraction
      Perhaps we would prefer all of these be on the method call eventually??

    Returns
    -------
    frame : ndarray
      The input frame with the overscan region subtracted
    """
    #dnum = settings.get_dnum(det)
    if settings is None:
        settings = dict(reduce={'overscan': {'method': 'savgol', 'params': [5,65]}})


    for i in range(numamplifiers):
        # Determine the section of the chip that contains data
        #datasec = "datasec{0:02d}".format(i+1)
        #dx0, dx1 = settings.spect[dnum][datasec][0][0], settings.spect[dnum][datasec][0][1]
        #dy0, dy1 = settings.spect[dnum][datasec][1][0], settings.spect[dnum][datasec][1][1]
        dx0, dx1 = datasec[i][0][0], datasec[i][0][1]
        dy0, dy1 = datasec[i][1][0], datasec[i][1][1]
        if dx0 < 0: dx0 += rawframe.shape[0]
        if dx1 <= 0: dx1 += rawframe.shape[0]
        if dy0 < 0: dy0 += rawframe.shape[1]
        if dy1 <= 0: dy1 += rawframe.shape[1]
        xds = np.arange(dx0, dx1)
        yds = np.arange(dy0, dy1)
        # Determine the section of the chip that contains the overscan region
        #oscansec = "oscansec{0:02d}".format(i+1)
        #ox0, ox1 = settings.spect[dnum][oscansec][0][0], settings.spect[dnum][oscansec][0][1]
        #oy0, oy1 = settings.spect[dnum][oscansec][1][0], settings.spect[dnum][oscansec][1][1]
        ox0, ox1 = oscansec[i][0][0], oscansec[i][0][1]
        oy0, oy1 = oscansec[i][1][0], oscansec[i][1][1]
        if ox0 < 0: ox0 += rawframe.shape[0]
        if ox1 <= 0: ox1 += min(rawframe.shape[0], dx1)  # Truncate to datasec
        if oy0 < 0: oy0 += rawframe.shape[1]
        if oy1 <= 0: oy1 += min(rawframe.shape[1], dy1)  # Truncate to datasec
        xos = np.arange(ox0, ox1)
        yos = np.arange(oy0, oy1)
        w = np.ix_(xos, yos)
        oscan = rawframe[w]
        # Make sure the overscan section has at least one side consistent with datasec
        if dx1-dx0 == ox1-ox0:
            osfit = np.median(oscan, axis=1)  # Mean was hit by CRs
        elif dy1-dy0 == oy1-oy0:
            osfit = np.median(oscan, axis=0)
        elif settings['reduce']['overscan']['method'].lower() == "median":
            osfit = np.median(oscan)
        else:
            msgs.error("Overscan sections do not match amplifier sections for amplifier {0:d}".format(i+1))
        # Fit/Model the overscan region
        if settings['reduce']['overscan']['method'].lower() == "polynomial":
            c = np.polyfit(np.arange(osfit.size), osfit, settings['reduce']['overscan']['params'][0])
            ossub = np.polyval(c, np.arange(osfit.size))#.reshape(osfit.size,1)
        elif settings['reduce']['overscan']['method'].lower() == "savgol":
            ossub = signal.savgol_filter(osfit,
                                         settings['reduce']['overscan']['params'][1],
                                         settings['reduce']['overscan']['params'][0])
        elif settings['reduce']['overscan']['method'].lower() == "median":  # One simple value
            ossub = osfit * np.ones(1)
        else:
            msgs.warn("Overscan subtraction method {0:s} is not implemented".format(settings['reduce']['overscan']['method']))
            msgs.info("Using a linear fit to the overscan region")
            c = np.polyfit(np.arange(osfit.size), osfit, 1)
            ossub = np.polyval(c, np.arange(osfit.size))#.reshape(osfit.size,1)
        # Determine the section of the chip that contains data for this amplifier
        if i==0:
            frame = rawframe.copy()
        wd = np.ix_(xds, yds)
        ossub = ossub.reshape(osfit.size, 1)
        if wd[0].shape[0] == ossub.shape[0]:
            frame[wd] -= ossub
        elif wd[1].shape[1] == ossub.shape[0]:
            frame[wd] -= ossub.T
        elif settings['reduce']['overscan']['method'].lower() == "median":
            frame[wd] -= osfit
        else:
            msgs.error("Could not subtract bias from overscan region --"+msgs.newline()+"size of extracted regions does not match")
    # Return
    del xds, yds, xos, yos, oscan
    return frame


def replace_columns(img, bad_cols, replace_with='mean'):
    """ Replace bad columns with values from the neighbors

    Parameters
    ----------
    img : ndarray
    bad_cols: ndarray (bool, 1D, shape[1] of img)
      True = bad column
      False = ok column
    replace_with : str, optional
      Option for replacement
       mean -- Use the mean of the closest left/right columns

    Returns
    -------
    img2 : ndarray
      Copy of the input image with the bad columns replaced
    """
    # Prep
    img2 = img.copy()
    # Find the starting/ends of the bad column sets
    tmp = np.zeros(img.shape[1], dtype=int)
    tmp[bad_cols] = 1
    tmp2 = tmp - np.roll(tmp,1)
    # Deal with first column
    if bad_cols[0]:
        tmp2[0]=1
    ledges = np.where(tmp2 == 1)[0]
    redges = np.where(tmp2 == -1)[0]
    # Last column
    if tmp2[-1] == 1:
        redges = np.concatenate([redges, np.array([bad_cols.size-1])])
    # Loop on em
    for kk, ledge in enumerate(ledges):
        lval = img[:,ledge-1]
        rval = img[:,redges[kk]]
        # Replace
        if replace_with == 'mean':
            mval = (lval+rval)/2.
            for ii in range(ledge, redges[kk]+1):
                img2[:,ii] = mval
        else:
            msgs.error("Bad option to replace_columns")
    # Return
    return img2


def trim(frame, numamplifiers, datasec):
    """ Core method to trim an input image

    Parameters
    ----------
    frame : ndarray
    numamplifiers : int
    datasec : list of datasecs
      One per amplifier

    Returns
    -------
    frame : ndarray
      Trimmed
    """
    for i in range(numamplifiers):
        #datasec = "datasec{0:02d}".format(i+1)
        #x0, x1 = settings.spect[dnum][datasec][0][0], settings.spect[dnum][datasec][0][1]
        #y0, y1 = settings.spect[dnum][datasec][1][0], settings.spect[dnum][datasec][1][1]
        x0, x1 = datasec[i][0][0], datasec[i][0][1]
        y0, y1 = datasec[i][1][0], datasec[i][1][1]
        # Fuss with edges
        if x0 < 0:
            x0 += frame.shape[0]
        if x1 <= 0:
            x1 += frame.shape[0]
        if y0 < 0:
            y0 += frame.shape[1]
        if y1 <= 0:
            y1 += frame.shape[1]
        if i == 0:
            xv = np.arange(x0, x1)
            yv = np.arange(y0, y1)
        else:
            xv = np.unique(np.append(xv, np.arange(x0, x1)))
            yv = np.unique(np.append(yv, np.arange(y0, y1)))
    # Construct and array with the rows and columns to be extracted
    w = np.ix_(xv, yv)
#	if len(file.shape) == 2:
#		trimfile = file[w]
#	elif len(file.shape) == 3:
#		trimfile = np.zeros((w[0].shape[0],w[1].shape[1],file.shape[2]))
#		for f in range(file.shape[2]):
#			trimfile[:,:,f] = file[:,:,f][w]
#	else:
#		msgs.error("Cannot trim {0:d}D frame".format(int(len(file.shape))))
    try:
        return frame[w]
    except:
        msgs.bug("Odds are datasec is set wrong. Maybe due to transpose")
        debugger.set_trace()
        msgs.error("Cannot trim file")


def variance_frame(datasec_img, det, sciframe, settings_det=None,
                   fitsdict=None, skyframe=None, objframe=None,
                   idx=None, dnoise=None):
    """ Calculate the variance image including detector noise
    Parameters
    ----------
    datasec_img : ndarray
    det
    sciframe
    settings_det : settings.spect[dnum]
      Detector info
    fitsdict : dict, optional
      Contains relevant information from fits header files
    idx : int, optional
    objframe : ndarray, optional
      Model of object counts
    Returns
    -------
    variance image : ndarray
    """
    # The effective read noise (variance image)
    rnoise = rn_frame(det, datasec_img, settings_det)
    if skyframe is not None:
        if objframe is None:
            objframe = np.zeros_like(skyframe)
        varframe = np.abs(skyframe + objframe - np.sqrt(2)*np.sqrt(rnoise)) + rnoise
        return varframe
    else:
        scicopy = sciframe.copy()
        # Dark Current noise
        if dnoise is None:
            dnoise = (settings_det['darkcurr'] * float(fitsdict["exptime"][idx])/3600.0)
        # Return
        return np.abs(scicopy) + rnoise + dnoise


