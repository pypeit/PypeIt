import numpy as np
from matplotlib import pyplot as plt

MAX_REJECT = 0.5
MIN_NPIXELS = 5
GOOD_PIXEL = 0
BAD_PIXEL = 1
KREJ = 2.5
MAX_ITERATIONS = 5


def get_dimen(x, maxp=25):
    """
    Assign the plotting dimensions to be the "most square"
    x     : An integer that equals the number of panels to be plot
    maxp  : The maximum number of panels to plot on a single page
    """
    pages, npp = [], []
    xr = x
    while xr > 0:
        if xr > maxp: xt = maxp
        else: xt = xr
        ypg = int(np.sqrt(np.float(xt)))
        if int(xt) % ypg == 0: xpg = int(xt)/ypg
        else: xpg = 1 + int(xt)/ypg
        pages.append([xpg, ypg])
        npp.append(xt)
        xr -= xpg*ypg
    return pages, npp


def zscale(image, nsamples=1000, contrast=0.25, bpmask=None, zmask=None):
    """
    Implement IRAF zscale algorithm
    nsamples=1000 and contrast=0.25 are the IRAF display task defaults
    bpmask and zmask not implemented yet
    image is a 2-d np array
    returns (z1, z2)
    """

    # Sample the image
    samples = zsc_sample(image, nsamples, bpmask, zmask)
    npix = len(samples)
    samples.sort()
    zmin = samples[0]
    zmax = samples[-1]
    # For a zero-indexed array
    center_pixel = (npix - 1) // 2
    if npix % 2 == 1:
        median = samples[center_pixel]
    else:
        median = 0.5 * (samples[center_pixel] + samples[center_pixel + 1])

    #
    # Fit a line to the sorted array of samples
    minpix = np.max([MIN_NPIXELS, int(npix * MAX_REJECT)])
    ngrow = np.max([1, int (npix * 0.01)])
    ngoodpix, zstart, zslope = zsc_fit_line(samples, npix, KREJ, ngrow, MAX_ITERATIONS)

    if ngoodpix < minpix:
        z1 = zmin
        z2 = zmax
    else:
        if contrast > 0: zslope /= contrast
        z1 = np.max([zmin, median - (center_pixel - 1) * zslope])
        z2 = np.min([zmax, median + (npix - center_pixel) * zslope])
    return z1, z2


def zsc_sample(image, maxpix, bpmask=None, zmask=None):
    # Figure out which pixels to use for the zscale algorithm
    # Returns the 1-d array samples
    # Don't worry about the bad pixel mask or zmask for the moment
    # Sample in a square grid, and return the first maxpix in the sample
    nc = image.shape[0]
    nl = image.shape[1]
    stride = np.max([1.0, np.sqrt((nc - 1) * (nl - 1) / float(maxpix))])
    stride = int(stride)
    samples = image[::stride,::stride].flatten()
    return samples[:maxpix]


def zsc_fit_line(samples, npix, krej, ngrow, maxiter):
    # First re-map indices from -1.0 to 1.0
    xscale = 2.0 / (npix - 1)
    xnorm = np.arange(npix)
    xnorm = xnorm * xscale - 1.0

    ngoodpix = npix
    minpix = np.max([MIN_NPIXELS, int (npix*MAX_REJECT)])
    last_ngoodpix = npix + 1

    # This is the mask used in k-sigma clipping.0 is good, 1 is bad
    badpix = np.zeros(npix, dtype="int32")

    #Iterate
    for niter in xrange(maxiter):

        if (ngoodpix >= last_ngoodpix) or (ngoodpix < minpix):
            break

        # Accumulate sums to calculate straight line fit
        goodpixels = np.where(badpix == GOOD_PIXEL)
        sumx = xnorm[goodpixels].sum()
        sumxx = (xnorm[goodpixels]*xnorm[goodpixels]).sum()
        sumxy = (xnorm[goodpixels]*samples[goodpixels]).sum()
        sumy = samples[goodpixels].sum()
        sum = len(goodpixels[0])

        delta = sum * sumxx - sumx * sumx
        # Slope and intercept
        intercept = (sumxx * sumy - sumx * sumxy) / delta
        slope = (sum * sumxy - sumx * sumy) / delta

        # Subtract fitted line from the data array
        fitted = xnorm*slope + intercept
        flat = samples - fitted

        # Compute the k-sigma rejection threshold
        ngoodpix, mean, sigma = zsc_compute_sigma(flat, badpix)

        threshold = sigma * krej

        # Detect and reject pixels further than k*sigma from the fitted line
        lcut = -threshold
        hcut = threshold
        below = np.where(flat < lcut)
        above = np.where(flat > hcut)

        badpix[below] = BAD_PIXEL
        badpix[above] = BAD_PIXEL

        # Convolve with a kernel of length ngrow
        kernel = np.ones(ngrow, dtype="int32")
        badpix = np.convolve(badpix, kernel, mode='same')

        ngoodpix = len(np.where(badpix == GOOD_PIXEL)[0])

        niter += 1

    # Transform the line coefficients back to the X range [0:npix-1]
    zstart = intercept - slope
    zslope = slope * xscale

    return ngoodpix, zstart, zslope


def zsc_compute_sigma(flat, badpix):
    """
    Compute the rms deviation from the mean of a flattened array.
    Ignore rejected pixels
    """

    # Accumulate sum and sum of squares
    goodpixels = np.where(badpix == GOOD_PIXEL)
    sumz = flat[goodpixels].sum()
    sumsq = (flat[goodpixels]*flat[goodpixels]).sum()
    ngoodpix = len(goodpixels[0])
    if ngoodpix == 0:
        mean = None
        sigma = None
    elif ngoodpix == 1:
        mean = sumz
        sigma = None
    else:
        mean = sumz / ngoodpix
        temp = sumsq / (ngoodpix - 1) - sumz*sumz / (ngoodpix * (ngoodpix - 1))
        if temp < 0:
            sigma = 0.0
        else:
            sigma = np.sqrt (temp)

    return ngoodpix, mean, sigma


def plot_orderfits(slf, model, ydata, xdata=None, xmodl=None, textplt="Order", maxp=4, desc="", maskval=-999999.9):
    """
    Given data and a model, the model+fit of each order is saved to a pdf
    """
    npix, nord = ydata.shape
    pages, npp = get_dimen(nord, maxp=maxp)
    if xdata is None: xdata = np.arange(npix).reshape((npix, 1)).repeat(nord, axis=1)
    if xmodl is None: xmodl = np.arange(model.shape[0])
    # Loop through all pages and plot the results
    ndone = 0
    axesIdx = True
    for i in xrange(len(pages)):
        f, axes = plt.subplots(pages[i][1], pages[i][0])
        ipx, ipy = 0, 0
        for j in xrange(npp[i]):
            if pages[i][0] == 1 and pages[i][1] == 1: axesIdx = False
            elif pages[i][1] == 1: ind = (ipx,)
            elif pages[i][0] == 1: ind = (ipy,)
            else: ind = (ipy, ipx)
            if axesIdx:
                axes[ind].plot(xdata[:,ndone+j], ydata[:,ndone+j], 'bx', drawstyle='steps')
                axes[ind].plot(xmodl, model[:,ndone+j], 'r-')
            else:
                axes.plot(xdata[:,ndone+j], ydata[:,ndone+j], 'bx', drawstyle='steps')
                axes.plot(xmodl, model[:,ndone+j], 'r-')
            ytmp = ydata[:,ndone+j]
            ytmp = ytmp[np.where(ytmp != maskval)]
            if ytmp.size != 0: amn = min(np.min(ytmp), np.min(model[:,ndone+j]))
            else: amn = np.min(model[:,ndone+j])
            if ytmp.size != 0: amx = max(np.max(ytmp), np.max(model[:,ndone+j]))
            else: amx = np.max(model[:,ndone+j])
            xtmp = xdata[:,ndone+j]
            xtmp = xtmp[np.where(xtmp != maskval)]
            if xtmp.size == 0:
                xmn = np.min(xmodl)
                xmx = np.max(xmodl)
            else:
                xmn = min(np.min(xtmp), np.min(xmodl))
                xmx = max(np.max(xtmp), np.max(xmodl))
            if axesIdx:
                axes[ind].axis([xmn, xmx, amn, amx])
                axes[ind].set_title("{0:s} {1:d}".format(textplt, ndone+j+1))
            else:
                axes.axis([xmn, xmx, amn, amx])
                axes.set_title("{0:s} {1:d}".format(textplt, ndone+j+1))
            ipx += 1
            if ipx == pages[i][0]:
                ipx = 0
                ipy += 1
        # Delete the unnecessary axes
        if axesIdx:
            for j in xrange(npp[i], axes.size):
                if pages[i][1] == 1: ind = (ipx,)
                elif pages[i][0] == 1: ind = (ipy,)
                else: ind = (ipy, ipx)
                f.delaxes(axes[ind])
                if ipx == pages[i][0]:
                    ipx = 0
                    ipy += 1
        ndone += npp[i]
        # Save the figure
        if axesIdx: axsz = axes.size
        else: axsz = 1.0
        if pages[i][1] == 1 or pages[i][0] == 1: ypngsiz = 11.0/axsz
        else: ypngsiz = 11.0*axes.shape[0]/axes.shape[1]
        f.set_size_inches(11.0, ypngsiz)
        if desc != "":
            pgtxt = ""
            if len(pages) != 1:
                pgtxt = ", page {0:d}/{1:d}".format(i+1, len(pages))
            f.suptitle(desc + pgtxt, y=1.02, size=16)
        f.tight_layout()
        slf._qa.savefig(dpi=200, orientation='landscape', bbox_inches='tight')
        plt.close()
        f.clf()
    del f
    return
