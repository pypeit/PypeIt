"""
Module with slit - mask design matching routines.

Routines are primarily used for matching the traced slit edges to the predicted trace
from the mask design/optical model.

TODO: These routines are specific for DEIMOS. Can they be generalized?

These routines are taken from the DEEP2 IDL-based pipeline.

"""
# from IPython import embed

import numpy
from scipy.optimize import curve_fit

from matplotlib import pyplot as plt

from astropy.stats import sigma_clipped_stats

# from pypeit import msgs
# from pypeit import utils


def linefit(x, a, b):
    return a + x*b


def slit_coeff_eval(x, y, coeff):
    """
    Taken from DEEP2/spec2d/pro/deimos_slit_match.pro
    """
    return coeff[0] + coeff[1]*x + coeff[2]*y



def best_offset(x_det, x_model, step=1, xlag_range=None):
    """
    Taken from discrete_correlate.pro in DEEP2/spec2d/pro/
    x_det==x1, x_model==x2
    """
    nbest = int(x_det.size * .85)
    if xlag_range is not None:
        xlag = numpy.arange(xlag_range[0], xlag_range[1]+step, step)
        min_x_det, max_x_det = numpy.min(x_det), numpy.max(x_det)
        wkeep =(x_model > min_x_det+xlag_range[0]) & (x_model < max_x_det+xlag_range[1])
        if x_model[wkeep].size<2:
            print('Working between {} and {}'.format(min_x_det+xlag_range[0], max_x_det+xlag_range[1]))
            print('Not enough lines to run!!!')
            sdev=1e10
            return 0.
        x_model_trim=x_model[wkeep]
    else:
        min_x_model, max_x_model = numpy.min(x_model), numpy.max(x_model)
        max_x_det = numpy.max(x_det)
        xlag = numpy.arange((max_x_model-min_x_model+max_x_det)/step)*step+min_x_model-max_x_det
        x_model_trim = x_model
    # -------- Store results in sdev
    sdev = numpy.zeros(xlag.size)

    # -------- Loop over lags
    for j in range(xlag.size):
        # for each element of x_det, get closest value of x_model
        x_det_lag = x_det+xlag[j]
        join = numpy.concatenate([x_det_lag, x_model_trim])
        sind = numpy.argsort(join)
        nj = sind.size
        w1 = numpy.where(sind<x_det.size)

        # [IDL-version comment] the following code is incorrect if the first element or last element in the
        # joined array is an element of x_det
        if x_det.size > 10:
            offs1 = (join[sind[((w1[0] + 1)>(nj - 1)).choose((w1[0] + 1),(nj - 1))]] - x_det_lag)
            offs2 = (x_det_lag - join[sind[((w1[0] - 1)<0).choose(w1[0] - 1,0)]])
            offs = (offs1 > offs2).choose(offs1, offs2)
        else:
            # [IDL-version comment] so added this brute-force version - still not quite right,
            # as assumes don't match 2 x_det's to the same x_model
            offs=numpy.zeros((x_det.size))
            for xelement in range(x_det.size):
                offs[xelement]=numpy.min(numpy.abs(x_det_lag[xelement]-x_model_trim))

            # use only nbest best matches
            soffs = numpy.argsort(numpy.abs(offs))
            nbest2= nbest if nbest<x_det.size else x_det.size
            offs = offs[soffs[0:nbest2]]

        sdev[j] = (offs**2).sum()          # big if match is bad

    best_sdev = int(numpy.mean(numpy.where(sdev==sdev.min())[0]))      # average ind
    #sdev = sdev[best_sdev]

    return xlag[best_sdev]


def discrete_correlate_match(x_det, x_model, step=1, xlag_range=[-50, 50]):
    """
    Taken from discrete_correlate_match.pro in DEEP2/spec2d/pro/
    x_det==x1, x_model==x2_in
    """
    best_off = best_offset(x_det, x_model, step=step, xlag_range=xlag_range)
    x_model_new = x_model - best_off

    ind = numpy.zeros(x_det.size, dtype=int)
    for i in range(ind.size):
        # ind[i] = numpy.where(numpy.min(numpy.abs(x_det[i] - x_model_new)) == numpy.abs(x_det[i] - x_model_new))[0]
        ind[i] = numpy.argmin(numpy.abs(x_det[i] - x_model_new))

    dx = x_det - x_model_new[ind]
    coeff, cov = curve_fit(linefit, x_det, dx, p0=[0, 1])
    scale = 1 + coeff[1] if x_det.size > 4 else 1
    x_model_new *= scale

    new_best_off = best_offset(x_det, x_model_new, step=step, xlag_range=xlag_range)

    x_model_new -= new_best_off

    for i in range(ind.size):
        # ind[i] = numpy.where(numpy.min(numpy.abs(x_det[i] - x_model_new)) == numpy.abs(x_det[i] - x_model_new))[0]
        ind[i] = numpy.argmin(numpy.abs(x_det[i] - x_model_new))

    dx = x_det - x_model_new[ind]
    x_model_new += numpy.median(dx)

    for i in range(ind.size):
        # ind[i] = numpy.where(numpy.min(numpy.abs(x_det[i] - x_model_new)) == numpy.abs(x_det[i] - x_model_new))[0]
        ind[i] = numpy.argmin(numpy.abs(x_det[i] - x_model_new))

    return ind


def slit_match(x_det, x_model, step=1, xlag_range=[-50,50], sigrej=3, print_matches=False, edge=None):
    """
    Taken from DEEP2/spec2d/pro/deimos_slit_match.pro
    """

    ind=discrete_correlate_match(x_det, x_model, step=step, xlag_range=xlag_range)

    residual = (x_det-x_model[ind]) - numpy.median(x_det-x_model[ind])
    weights = numpy.zeros(residual.size, dtype=int)
    weights[numpy.abs(residual) < 100.] = 1
    if weights.sum() == 0:
        weights=numpy.ones(residual.size, dtype=int)
    coeff, cov= curve_fit(linefit, x_model[ind]*weights, x_det, p0=[0,1])
    coeff=numpy.append(coeff, 0)
    yfit = slit_coeff_eval(x_model[ind], x_model[ind]*0., coeff)

    # compute residual
    res = yfit - x_det
    sigres = sigma_clipped_stats(res, sigma = sigrej)[2]   # RMS residuals
    cut = 5 if res.size<5 else sigrej*sigres
    out = numpy.abs(res)> cut

    # check for duplicate indices
    bad = numpy.ones(ind.size, dtype=int)
    bad[numpy.unique(ind, return_index=True)[1]] = 0
    wbad = numpy.where(bad)[0]
    if wbad.size > 0:
        print('Trimming {} bad matches in deimos_slit_match'.format(wbad.size))
        for i in range(wbad.size):
            badind = ind[wbad[i]]
            w = ind == badind
            bad[w] = 1    # set all with (ind eq badind) bad for a moment
            wdif = numpy.argmin(numpy.abs(res[w]))
            bad[w[wdif]] = 0  # set this one not bad
        bad = bad|out
        good = bad == 0
        ind = ind[good]
        x_det=x_det[good]
    if print_matches is True:
        if edge is not None:
            print('-----------------------------------------------')
            print('             {} slit edge               '.format(edge))
        print('-----------------------------------------------')
        print('Index      xblupix       found_x               ')
        print('-----------------------------------------------')
        for i in range(ind.size):
            print('{}  {}  {}'.format(ind[i], x_model[ind][i], x_det[i]))
    return x_det, ind, coeff, sigres


def plot_matches(edgetrace, ind, x_model, x_det, yref, slit_index, edge=None, shape=None):
    yref_xdet=numpy.tile(yref, x_det.size)
    yref_x_model=numpy.tile(yref, x_model.size)

    _shape = (4096, 2048) if shape is None else shape
    buffer = 20
    dist = _shape[0] - yref

    plt.rc('xtick', direction='in')
    plt.rc('ytick', direction='in')
    fig=plt.figure(figsize=(10, 4.5))
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)
    if edge is not None:
        plt.title('{} slit edges cross-matching'.format(edge))
    for x in edgetrace.T:
        plt.plot(x, numpy.arange(x.size), color='k', lw=0.5, zorder=0)
    plt.scatter(x_det, yref_xdet, marker='D', s=10, lw=0, color='m', zorder=1,
                                      label='Image trace midpoint (x_det)')
    plt.scatter(x_model, yref_x_model, marker='o', s=10, lw=0, color='b', zorder=1,
                                      label='Optical model trace BEFORE x-corr (x_model)')
    plt.scatter(x_model[ind], yref_x_model[ind], marker='o', s=40, facecolors='none', edgecolors='g', zorder=1,
                                      label='Optical model trace AFTER x-corr (x_model[ind])')
    for i in range(x_model[ind].size):
        plt.text(x_model[ind][i], yref_x_model[ind][i]+0.05*dist, slit_index[ind][i], rotation=45, color='g',
                                     fontsize=8, horizontalalignment='center')
    for i in range(x_model.size):
        if (x_model[i] >= buffer) and (x_model[i] <= _shape[1] + buffer):
            plt.text(x_model[i], yref_x_model[i]-0.15*dist, slit_index[i], rotation=45, color='b',
                                     fontsize=8, horizontalalignment='center')
    plt.xlabel('Spatial pixels')
    plt.ylabel('Spectral pixels')
    plt.xlim(buffer, _shape[1]+buffer)
    plt.ylim(0, _shape[0])
    plt.legend(loc=1)


def slit_match_fix(need, ind):
    """
    Taken from deimos_slit_match_fix in DEEP2/spec2d/pro/deimos_slit_match.pro
    """

    # synth = numpy.zeros(ind.size, dtype=int)
    needadd = need.copy()
    needadd[ind] = False
    needind = numpy.where(needadd)[0]  # slits we are missing
    # if needind.size > 0:
    #     ind = numpy.append(ind, needind)
    #     synth = numpy.append(synth, numpy.ones(needind.shape, dtype=int))
    # sind = numpy.argsort(ind)
    # ind = ind[sind]
    # synth = synth[sind]

    return needind     # ind, synth
