import numpy as np
import armsgs as msgs
import arcyutils
from matplotlib import pyplot as plt
import pdb

def boxcar(sciframe, varframe, crmask, scitrace, maskval=-999999.9, weighted=True):
    """

    :param sciframe: science frame
    :param varframe: variance image
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
        pdb.set_trace()
        msgs.info("Performing boxcar extraction on object {0:d}/{1:d}".format(o+1,nobj))
        # Fit the background
        msgs.info("   Fitting the background")
        bgframe = arcyutils.func2d_fit_val(bgfit, sciframe, scitrace['background'][:,:,o]*mask, bgfitord)
        if weighted: weight = np.abs(scitrace['object'][:,:,o]*mask*(sciframe-bgframe))
        else: weight = scitrace['object'][:,:,o]*mask
        sumweight = np.sum(weight, axis=1)
        # Total the object flux
        msgs.info("   Summing object counts")
        scisum = np.sum((sciframe-bgframe)*weight, axis=1)
        scisum /= sumweight
        # Total the variance array
        msgs.info("   Summing variance array")
        varsum = np.sum(varframe*weight, axis=1)
        varsum /= sumweight
        # Mask zero weights
        w = np.where(sumweight == 0.0)
        if w[0].size != 0:
            scisum[w] = maskval
            varsum[w] = abs(maskval)
        pltv = np.arange(scisum.size)
        plt.clf()
        plt.plot(pltv, scisum, 'k-', drawstyle='steps')
        plt.plot(pltv, np.sqrt(varsum), 'r-')
        plt.plot(pltv, -np.sqrt(varsum), 'r-')
        plt.show()
    return None
