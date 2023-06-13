"""
Module with slit - mask design matching routines.

Routines are primarily used for matching the traced slit edges to the predicted trace
from the mask design/optical model.

TODO: These routines are specific for DEIMOS. Can they be generalized?

These routines are taken from the DEEP2 IDL-based pipeline.

.. include:: ../include/links.rst

"""
from IPython import embed

import numpy as np
from matplotlib import pyplot as plt

from astropy.stats import sigma_clipped_stats
from pypeit.core import fitting

from pypeit import msgs


def best_offset(x_det, x_model, step=1, xlag_range=None):
    """
    Script to determine the best offset between the slit edge predicted
    by the optical model and the one found in the image. This is used iteratively.

    Taken from DEEP2/spec2d/pro/ discrete_correlate.pro
    x_det==x1, x_model==x2

    Args:
        x_det (`numpy.ndarray`_):
            1D array of slit edge spatial positions found from image
        x_model (`numpy.ndarray`_):
            1D array of slit edge spatial positions predicted by the optical model
        step (:obj:`int`):
            step size in pixels used to generate a list of possible offsets within the `offsets_range`
        xlag_range (:obj:`list`, optional):
            range of offsets in pixels allowed between the slit positions predicted by
            the mask design and the traced slit positions.

    Returns:
        :obj:`float`: best offset between the slit edge predicted by the optical model
        and the one found in the image

    """
    # This select the number of best matches used later fo statistics
    nbest = int(x_det.size * .85)
    # Genarate an array of offsets
    if xlag_range is not None:
        xlag = np.arange(xlag_range[0], xlag_range[1]+step, step)
        min_x_det, max_x_det = np.min(x_det), np.max(x_det)
        # we keep only the x_model values that are in the current detector
        wkeep =(x_model > min_x_det+xlag_range[0]) & (x_model < max_x_det+xlag_range[1])
        if x_model[wkeep].size<2:
            msgs.warn('Working between {} and {}'.format(min_x_det+xlag_range[0], max_x_det+xlag_range[1]))
            msgs.warn('Not enough lines to run!!!')
            sdev = 1e10
            return 0.
        x_model_trim = x_model[wkeep]
    else:
        min_x_model, max_x_model = np.ma.min(x_model), np.ma.max(x_model)
        max_x_det = np.max(x_det)
        xlag = np.arange(min_x_model-max_x_det, max_x_model, step)
        x_model_trim = x_model
    # The results will be stored in sdev
    sdev = np.zeros(xlag.size)

    # Loop over the array of offsets
    # so that for each element of x_det we get the closest value of x_model
    for j in range(xlag.size):
        x_det_lag = x_det+xlag[j]
        join = np.ma.concatenate([x_det_lag, x_model_trim])
        sind = np.argsort(join)
        nj = sind.size
        w1 = np.where(sind < x_det.size)

        # [IDL-version comment] the following code is incorrect if the first element or last element in the
        # joined array is an element of x_det
        if x_det.size > 10:
            offs1 = (join[sind[((w1[0] + 1) > (nj - 1)).choose((w1[0] + 1), (nj - 1))]] - x_det_lag)
            offs2 = (x_det_lag - join[sind[((w1[0] - 1) < 0).choose(w1[0] - 1, 0)]])
            offs = (offs1 > offs2).choose(offs1, offs2)  # set all values of offs1 > offs2 equal to offs2
        else:
            # [IDL-version comment] so added this brute-force version - still not quite right,
            # as assumes don't match 2 x_det's to the same x_model
            offs = np.amin(np.absolute(x_det_lag[:, None] - x_model_trim[None, :]), axis=1)

            # use only nbest best matches
            soffs = np.argsort(np.abs(offs))
            nbest2 = nbest if nbest<x_det.size else x_det.size
            offs = offs[soffs[0:nbest2]]

        # record the offset due to `xlag[j]` in `sdev` as a sum over `x_det` size.
        sdev[j] = (offs**2).sum()          # big if match is bad

    # The best offset will be the one with the smallest `sdev`
    best_sdev = int(np.mean(np.argmin(sdev)))      # average ind

    return xlag[best_sdev]


def discrete_correlate_match(x_det, x_model, step=1, xlag_range=[-50, 50]):
    """
    Script to find the the x_model values that match the traced edges.

    This method uses :func:`best_offset` to determine the best offset between
    slit edge predicted by the optical model and the one found in the image, given a range of
    offsets. This is used iteratively.

    Taken from in DEEP2/spec2d/pro/discrete_correlate_match.pro
    x_det==x1, x_model==x2_in

    Args:
        x_det (`numpy.ndarray`_):
            1D array of slit edge spatial positions found from image
        x_model (`numpy.ndarray`_):
            1D array of slit edge spatial positions predicted by the optical model
        step (:obj:`int`):
            step size in pixels used to generate a list of possible offsets within the `offsets_range`
        xlag_range (:obj:`list`, optional):
            range of offsets in pixels allowed between the slit positions predicted by
            the mask design and the traced slit positions.

    Returns:
        `numpy.ndarray`_: array of indices for x_model, which defines the matches to x_det,
        i.e., x_det matches x_model[ind]

    """
    # -------- PASS 1: get offset between x1 and x2

    # Determine the offset between x_det and x_model
    best_off = best_offset(x_det, x_model, step=step, xlag_range=xlag_range)
    # apply the offset to x_model
    x_model_new = x_model - best_off

    # for each traced edge (`x_det`) determine the value of x_model that gives the smallest offset
    ind = np.ma.argmin(np.ma.absolute(x_det[:, None] - x_model_new[None, :]), axis=1)

    # -------- PASS 2: remove linear trend (i.e. adjust scale)
    # fit the offsets to `x_det` to find the scale and apply it to x_model
    dx = np.ma.compressed(x_det - x_model_new[ind])
    pypeitFit = fitting.robust_fit(x_det, dx, 1, maxiter=100, lower=2, upper=2)
    coeff = pypeitFit.fitc
    scale = 1 + coeff[1] if x_det.size > 4 else 1
    x_model_new *= scale

    # Find again the best offset and apply it to x_model
    new_best_off = best_offset(x_det, x_model_new, step=step, xlag_range=xlag_range)

    x_model_new -= new_best_off

    # find again `ind`
    ind = np.ma.argmin(np.ma.absolute(x_det[:,None] - x_model_new[None,:]), axis=1)

    # -------- PASS 3: tweak offset
    dx = x_det - x_model_new[ind]
    x_model_new += np.ma.median(dx)

    # find again `ind`
    ind = np.ma.argmin(np.ma.absolute(x_det[:,None] - x_model_new[None,:]), axis=1)

    return ind


def slit_match(x_det, x_model, step=1, xlag_range=[-50,50], sigrej=3, print_matches=False,
               edge=None):
    """
    Script that perform the slit edges matching.

    This method uses :func:`discrete_correlate_match` to find the
    indices of x_model that match x_det.

    Taken from DEEP2/spec2d/pro/deimos_slit_match.pro

    Parameters
    ----------
    x_det: `numpy.ndarray`_
        1D array of slit edge spatial positions found from image.
    x_model: `numpy.ndarray`_
        1D array of slit edge spatial positions predicted by the
        optical model.
    step: :obj:`int`, optional
        Step size in pixels used to generate a list of possible
        offsets within the `offsets_range`.
    xlag_range: :obj:`list`, optional
        Range of offsets in pixels allowed between the slit
        positions predicted by the mask design and the traced
        slit positions.
    sigrej: :obj:`float`, optional
        Reject slit matches larger than this number of sigma in
        the match residuals.
    print_matches: :obj:`bool`, optional
        Print the result of the matching.
    edge: :obj:`str`, optional
        String that indicates which edges are being plotted,
        i.e., left of right. Ignored if ``print_matches`` is
        False.

    Returns
    -------
    ind: `numpy.ndarray`_
        1D array of indices for `x_model`, which defines the matches
        to `x_det`, i.e., `x_det` matches `x_model[ind]`
    dupl: `numpy.ndarray`_
        1D array of `bool` that flags which `ind` are duplicates.
    coeff: `numpy.ndarray`_
        pypeitFit coefficients of the fitted relation between `x_det`
        and `x_model[ind]`
    sigres: :obj:`float`
        RMS residual for the fitted relation between `x_det` and
        `x_model[ind]`

    """
    # Determine the indices of `x_model` that match `x_det`
    ind = discrete_correlate_match(x_det, np.ma.masked_equal(x_model, -1), step=step, xlag_range=xlag_range)

    # Define the weights for the fitting
    residual = (x_det-x_model[ind]) - np.median(x_det-x_model[ind])
    weights = np.zeros(residual.size, dtype=int)
    weights[np.abs(residual) < 100.] = 1
    if weights.sum() == 0:
        weights = np.ones(residual.size, dtype=int)
    # Fit between `x_det` and `x_model[ind]`
    pypeitFit = fitting.robust_fit(x_model[ind], x_det, 1, maxiter=100, weights=weights, lower=3, upper=3)
    coeff = pypeitFit.fitc
    yfit = pypeitFit.eval(x_model[ind])

    # compute residuals
    res = yfit - x_det
    sigres = sigma_clipped_stats(res, sigma=sigrej)[2]   # RMS residuals
    # flag the matches that have residuals > `sigrej` times the RMS, or if res>5
    cut = 5 if res.size < 5 else sigrej*sigres
    out = np.abs(res) > cut

    # check for duplicate indices
    dupl = np.ones(ind.size, dtype=bool)
    # If there are duplicates of `ind`, for now we keep only the first one. We don't remove the others yet
    dupl[np.unique(ind, return_index=True)[1]] = False
    wdupl = np.where(dupl)[0]
    # Iterate over the duplicates flagged as bad
    if wdupl.size > 0:
        for i in range(wdupl.size):
            duplind = ind[wdupl[i]]
            # Where are the other duplicates of this `ind`?
            w = np.where(ind == duplind)[0]
            # set those to be bad (for the moment)
            dupl[w] = True
            # Among the duplicates of this particular `ind`, which one has the smallest residual?
            wdif = np.argmin(np.abs(res[w]))
            # The one with the smallest residuals, is then set to not bad
            dupl[w[wdif]] = False
        # Both duplicates and matches with high RMS are considered bad
        dupl = dupl | out
        if edge is not None:
            msgs.warn('{} duplicate match(es) for {} edges'.format(dupl[dupl == 1].size, edge))
        else:
            msgs.warn('{} duplicate match(es)'.format(dupl[dupl == 1].size))
        # I commented the 3 lines below because I don't really need to trim the duplicate matches. I just
        # propagate the flag.
        # good = dupl == 0
        # ind = ind[good]
        # x_det=x_det[good]
    if print_matches:
        if edge is not None:
            msgs.info('-----------------------------------------------')
            msgs.info('             {} slit edges               '.format(edge))
        msgs.info('-----------------------------------------------')
        msgs.info('Index      omodel_edge       spat_edge               ')
        msgs.info('-----------------------------------------------')
        for i in range(ind.size):
            msgs.info('{}  {}  {}'.format(ind[i], x_model[ind][i], x_det[i]))
        msgs.info('-----------------------------------------------')
    return ind, dupl, coeff, sigres


def plot_matches(edgetrace, ind, x_model, yref, slit_index, nspat=2048, duplicates=None, missing=None, edge=None):
    r"""
    Plot the slit mask matching results.

    Args:
        edgetrace (`numpy.ndarray`_):
            2D array with the location of the slit edges for each
            spectral pixel as measured from the trace image. Shape is
            :math:`(N_{\rm spec},N_{\rm trace})`.
        ind (`numpy.ndarray`_):
            1D array of indices for `x_model`, which defines the
            matches to `x_det`.
        x_model (`numpy.ndarray`_):
            1D array of slit edge spatial positions predicted by the
            optical model.
        yref (:obj:`float`):
            Reference pixel in the `spec` direction.
        slit_index (`numpy.ndarray`_):
            1D array of slit-mask design indices.
        nspat (:obj:`int`, optional):
            Spatial dimension of the detector, for plotting purpose.
        duplicates (`numpy.ndarray`_, optional):
            1D array of `bool` that flags which `ind` are duplicates.
        missing (`numpy.ndarray`_, optional):
            1D array of indices for `x_model`, which defines the
            missing traces, if any.
        edge (:obj:`str`, optional):
            String that indicates which edges are being plotted,
            i.e., left of right.
    """

    # Slit edge spatial positions found from image at yref
    x_det = edgetrace[yref, :]

    yref_xdet = np.tile(yref, x_det.size)
    yref_x_model = np.tile(yref, x_model.size)

    buffer = 200
    dist = edgetrace.shape[0] - yref

    plt.rc('xtick', direction='in')
    plt.rc('ytick', direction='in')

    fig = plt.figure(figsize=(10, 4.5))
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)

    if edge is not None:
        plt.title('{} slit edges cross-matching'.format(edge))

    # Plot the traced edges
    for x in range(edgetrace.shape[1]):
        if duplicates is not None and duplicates[x]:
            plt.plot(edgetrace[:, x], np.arange(edgetrace[:, x].size), color='orange', lw=0.5, zorder=0)
        else:
            plt.plot(edgetrace[:, x], np.arange(edgetrace[:, x].size), color='k', lw=0.5, zorder=0)

    # Plot `x_det`, `x_model`, and `x_model[ind]` at a reference pixel in the `spec` direction
    plt.scatter(x_det, yref_xdet, marker='D', s=30, lw=1.2, facecolors='none', edgecolors='m', zorder=1,
                label='Image trace midpoint')
    plt.scatter(x_model[x_model != -1], yref_x_model[x_model != -1], marker='o', s=10, lw=0, color='b', zorder=1,
                label='Predicted optical model trace')
    plt.scatter(x_model[ind], yref_x_model[ind], marker='o', s=150, facecolors='none', edgecolors='g', zorder=1,
                label='Optical model trace matched to the image trace')
    if missing is not None:
        plt.scatter(x_model[missing], yref_x_model[missing], marker='x', s=40, color='r', zorder=1,
                    label='Optical model trace missing in image trace')

    # Print in the plot the values of `slintindx` for the matched edges
    for i in range(x_det.size):
        plt.text(x_det[i]+0.01*nspat, yref_xdet[i]+0.05*dist, slit_index[ind][i], rotation=45,
                 color='g', fontsize=8, horizontalalignment='center')
    for i in range(x_model.size):
        if x_model[i] != -1:
            plt.text(x_model[i]+0.01*nspat, yref_x_model[i]-0.15*dist, slit_index[i], rotation=45, color='b',
                     fontsize=8, horizontalalignment='center')
    if missing is not None:
        for i in range(x_model[missing].size):
            plt.text(x_model[missing][i]+0.01*nspat, yref_x_model[missing][i]+0.05*dist, slit_index[missing][i],
                     rotation=45, color='r', fontsize=8, horizontalalignment='center')

    plt.xlabel('Spatial pixels')
    plt.ylabel('Spectral pixels')
    plt.xlim(-buffer, nspat+buffer)
    plt.ylim(0, edgetrace.shape[0])
    plt.legend(loc=1)
    plt.show()
