import numpy as np
import scipy
import armsgs as msgs
import arcyutils
from matplotlib import pyplot as plt
import pdb

try:
    from xastropy.xutils import xdebug as xdb
except:
    pass

def boxcar(wvimg, sciframe, varframe, skyframe, crmask, scitrace, maskval=-999999.9, weighted=True):
    """

    :param wvimg: wavelength image
    :param sciframe: science frame
    :param varframe: variance image
    :param bgframe: sky background frame
    :param crmask: mask of cosmic ray hits
    :param scitrace: object and background trace images
    :param maskval: value to be used to mask bad pixels
    :param weighted: weight by the object flux
    :return:
    """
    bgfitord = 1  # Polynomial order used to fit the background
    nobj = scitrace['traces'].shape[1]
    mask = 1.0-crmask
    bgfit = np.linspace(0.0, 1.0, sciframe.shape[1])
    for o in range(nobj):
        #pdb.set_trace()
        msgs.info("Performing boxcar extraction on object {0:d}/{1:d}".format(o+1,nobj))
        # Fit the background
        msgs.info("   Fitting the background")
        bgframe = arcyutils.func2d_fit_val(bgfit, sciframe, scitrace['background'][:,:,o]*mask, bgfitord)
        if weighted: weight = np.abs(scitrace['object'][:,:,o]*mask*(sciframe-bgframe))
        else: weight = scitrace['object'][:,:,o]*mask
        sumweight = np.sum(weight, axis=1)
        # Generate wavelength array
        wvsum = np.sum(wvimg*weight, axis=1)
        wvsum /= sumweight
        # Generate sky spectrum (flux per pixel)
        skysum = np.sum(skyframe*weight, axis=1)
        skysum /= sumweight
        # Total the object flux
        msgs.info("   Summing object counts")
        scisum = np.sum((sciframe-bgframe)*weight, axis=1)
        if weighted:
            scisum /= sumweight
        # Total the variance array
        msgs.info("   Summing variance array")
        varsum = np.sum(varframe*weight, axis=1)
        if weighted:
            varsum /= sumweight
        # Mask zero weights
        w = sumweight <= 0.0
        if np.sum(w) > 0:
            scisum[w] = maskval
            varsum[w] = 0. #abs(maskval)
            # Wavelength -- Need to fill these in
            ival = np.arange(wvsum.size)
            fwv = scipy.interpolate.InterpolatedUnivariateSpline(ival[~w], wvsum[~w], k=2)
            wvsum[w] = fwv(ival[w]) # Includes extrapolation
            skysum[w] = 0. #abs(maskval)
        '''
        pltv = np.arange(scisum.size)
        plt.clf()
        plt.plot(pltv, scisum, 'k-', drawstyle='steps')
        plt.plot(pltv, np.sqrt(varsum), 'r-')
        #plt.plot(pltv, -np.sqrt(varsum), 'r-')
        plt.show()
        plt.close()
        '''
        #xdb.xplot(scisum/np.sqrt(varsum)) # S/N
        #xdb.set_trace()
    return wvsum, scisum, varsum, skysum
