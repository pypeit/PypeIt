"""
Module with slit - mask design matching routines.

Routines are primarily used for matching the traced slit edges to the predicted trace
from the mask design/optical model.

TODO: These routines are specific for DEIMOS. Can they be generalized?

These routines are taken from the DEEP2 IDL-based pipeline.

.. include:: ../include/links.rst

"""
# from IPython import embed

import numpy
from matplotlib import pyplot as plt

from astropy.stats import sigma_clipped_stats
from pypeit.core import fitting

# from pypeit import msgs


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
        xlag (:obj:`float`):
            best offset between the slit edge predicted by the optical model and the one found in the image

    """
    # This select the number of best matches used later fo statistics
    nbest = int(x_det.size * .85)
    # Genarate an array of offsets
    if xlag_range is not None:
        xlag = numpy.arange(xlag_range[0], xlag_range[1]+step, step)
        min_x_det, max_x_det = numpy.min(x_det), numpy.max(x_det)
        # we keep only the x_model values that are in the current detector
        wkeep =(x_model > min_x_det+xlag_range[0]) & (x_model < max_x_det+xlag_range[1])
        if x_model[wkeep].size<2:
            print('Working between {} and {}'.format(min_x_det+xlag_range[0], max_x_det+xlag_range[1]))
            print('Not enough lines to run!!!')
            sdev = 1e10
            return 0.
        x_model_trim = x_model[wkeep]
    else:
        min_x_model, max_x_model = numpy.min(x_model), numpy.max(x_model)
        max_x_det = numpy.max(x_det)
        xlag = numpy.arange(min_x_model-max_x_det, max_x_model, step)
        x_model_trim = x_model
    # The results will be stored in sdev
    sdev = numpy.zeros(xlag.size)

    # Loop over the array of offsets
    # so that for each element of x_det we get the closest value of x_model
    for j in range(xlag.size):
        x_det_lag = x_det+xlag[j]
        join = numpy.concatenate([x_det_lag, x_model_trim])
        sind = numpy.argsort(join)
        nj = sind.size
        w1 = numpy.where(sind < x_det.size)

        # [IDL-version comment] the following code is incorrect if the first element or last element in the
        # joined array is an element of x_det
        if x_det.size > 10:
            offs1 = (join[sind[((w1[0] + 1) > (nj - 1)).choose((w1[0] + 1), (nj - 1))]] - x_det_lag)
            offs2 = (x_det_lag - join[sind[((w1[0] - 1) < 0).choose(w1[0] - 1, 0)]])
            offs = (offs1 > offs2).choose(offs1, offs2)  # set all values of offs1 > offs2 equal to offs2
        else:
            # [IDL-version comment] so added this brute-force version - still not quite right,
            # as assumes don't match 2 x_det's to the same x_model
            offs = numpy.zeros(x_det.size)
            for xelement in range(x_det.size):
                offs[xelement] = numpy.min(numpy.abs(x_det_lag[xelement]-x_model_trim))

            # use only nbest best matches
            soffs = numpy.argsort(numpy.abs(offs))
            nbest2 = nbest if nbest<x_det.size else x_det.size
            offs = offs[soffs[0:nbest2]]

        # record the offset due to `xlag[j]` in `sdev` as a sum over `x_det` size.
        sdev[j] = (offs**2).sum()          # big if match is bad

    # The best offset will be the one with the smallest `sdev`
    best_sdev = int(numpy.mean(numpy.argmin(sdev)))      # average ind

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
        ind (`numpy.ndarray`_):
            array of indices for x_model, which defines the matches to x_det,
            i.e., x_det matches x_model[ind]

    """
    # -------- PASS 1: get offset between x1 and x2

    # Determine the offset between x_det and x_model
    best_off = best_offset(x_det, x_model, step=step, xlag_range=xlag_range)
    # apply the offset to x_model
    x_model_new = x_model - best_off

    # for each traced edge (`x_det`) determine the value of x_model that gives the smallest offset
    ind = numpy.zeros(x_det.size, dtype=int)
    for i in range(ind.size):
        ind[i] = numpy.argmin(numpy.abs(x_det[i] - x_model_new))

    # -------- PASS 2: remove linear trend (i.e. adjust scale)
    # fit the offsets to `x_det` to find the scale and apply it to x_model
    dx = x_det - x_model_new[ind]
    pypeitFit = fitting.robust_fit(x_det, dx, 1, maxiter=100, lower=3, upper=3)
    coeff = pypeitFit.fitc
    scale = 1 + coeff[1] if x_det.size > 4 else 1
    x_model_new *= scale

    # Find again the best offset and apply it to x_model
    new_best_off = best_offset(x_det, x_model_new, step=step, xlag_range=xlag_range)

    x_model_new -= new_best_off

    # find again `ind`
    for i in range(ind.size):
        ind[i] = numpy.argmin(numpy.abs(x_det[i] - x_model_new))

    # -------- PASS 3: tweak offset
    dx = x_det - x_model_new[ind]
    x_model_new += numpy.median(dx)

    # find again `ind`
    for i in range(ind.size):
        ind[i] = numpy.argmin(numpy.abs(x_det[i] - x_model_new))

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
    x_det: `numpy.ndarray`_
        1D array of slit edge spatial positions found from image
        trimmed if duplicated matches exist.
    ind: `numpy.ndarray`_
        1D array of indices for `x_model`, which defines the matches
        to `x_det`, i.e., `x_det` matches `x_model[ind]`
    coeff: `numpy.ndarray`_
        pypeitFit coefficients of the fitted relation between `x_det`
        and `x_model[ind]`
    sigres: :obj:`float`
        RMS residual for the fitted relation between `x_det` and
        `x_model[ind]`

    """
    # Determine the indices of `x_model` that match `x_det`
    ind = discrete_correlate_match(x_det, x_model, step=step, xlag_range=xlag_range)

    # Define the weights for the fitting
    residual = (x_det-x_model[ind]) - numpy.median(x_det-x_model[ind])
    weights = numpy.zeros(residual.size, dtype=int)
    weights[numpy.abs(residual) < 100.] = 1
    if weights.sum() == 0:
        weights = numpy.ones(residual.size, dtype=int)
    # Fit between `x_det` and `x_model[ind]`
    pypeitFit = fitting.robust_fit(x_model[ind], x_det, 1, maxiter=100, weights=weights, lower=3, upper=3)
    coeff = pypeitFit.fitc
    yfit = pypeitFit.eval(x_model[ind])

    # compute residuals
    res = yfit - x_det
    sigres = sigma_clipped_stats(res, sigma=sigrej)[2]   # RMS residuals
    # flag the matches that have residuals > `sigrej` times the RMS, or if res>5
    cut = 5 if res.size < 5 else sigrej*sigres
    out = numpy.abs(res) > cut

    # check for duplicate indices
    bad = numpy.ones(ind.size, dtype=int)
    # If there are duplicates of `ind`, for now we keep only the first one. We don't remove the others yet
    bad[numpy.unique(ind, return_index=True)[1]] = 0
    wbad = numpy.where(bad)[0]
    # Iterate over the duplicates flagged as bad
    if wbad.size > 0:
        for i in range(wbad.size):
            badind = ind[wbad[i]]
            # Where are the other duplicates of this `ind`?
            w = numpy.where(ind == badind)[0]
            # set those to be bad (for the moment)
            bad[w] = 1
            # Among the duplicates of this particular `ind`, which one has the smallest residual?
            wdif = numpy.argmin(numpy.abs(res[w]))
            # The one with the smallest residuals, is then set to not bad
            bad[w[wdif]] = 0
        # Both duplicates and matches with high RMS are trimmed
        bad = bad|out
        print('Trimming {} bad match(es) in deimos_slit_match'.format(bad[bad == 1].size))
        good = bad == 0
        ind = ind[good]
        x_det=x_det[good]
    if print_matches is True:
        if edge is not None:
            print('-----------------------------------------------')
            print('             {} slit edges               '.format(edge))
        print('-----------------------------------------------')
        print('Index      omodel_edge       spat_edge               ')
        print('-----------------------------------------------')
        for i in range(ind.size):
            print('{}  {}  {}'.format(ind[i], x_model[ind][i], x_det[i]))
    return x_det, ind, coeff, sigres


def plot_matches(edgetrace, ind, x_model, x_det, yref, slit_index, trimmed=None, edge=None, shape=None):
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
        x_det (`numpy.ndarray`_):
            1D array of slit edge spatial positions found from image.
        yref (:obj:`float`):
            Reference pixel in the `spec` direction.
        slit_index (`numpy.ndarray`_):
            1D array of slit-mask design indices.
        trimmed (:obj:`bool`, optional):
            True if the duplicated matches have been trimmed during
            :func:`slit_match`.
        edge (:obj:`str`, optional):
            String that indicates which edges are being plotted,
            i.e., left of right.
        shape (:obj:`tuple`, optional):
            Shape of the detector, for plotting purpose. It should be
            :math:`(N_{\rm spec}, N_{\rm spat})`.
    """
    yref_xdet = numpy.tile(yref, x_det.size)
    yref_x_model = numpy.tile(yref, x_model.size)

    _shape = (4096, 2048) if shape is None else shape
    buffer = 20
    dist = _shape[0] - yref

    plt.rc('xtick', direction='in')
    plt.rc('ytick', direction='in')

    fig = plt.figure(figsize=(10, 4.5))
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)

    if edge is not None:
        plt.title('{} slit edges cross-matching'.format(edge))

    # Plot the traced edges
    for x in range(edgetrace.T.shape[0]):
        if trimmed is not None and x in trimmed:
            plt.plot(edgetrace.T[x,:], numpy.arange(edgetrace.T[x,:].size), color='orange', lw=0.5, zorder=0)
        else:
            plt.plot(edgetrace.T[x, :], numpy.arange(edgetrace.T[x, :].size), color='k', lw=0.5, zorder=0)

    # Plot `x_det`, `x_model`, and `x_model[ind]` at a reference pixel in the `spec` direction
    plt.scatter(x_det, yref_xdet, marker='D', s=10, lw=0, color='m', zorder=1, label='Image trace midpoint (x_det)')
    plt.scatter(x_model, yref_x_model, marker='o', s=10, lw=0, color='b', zorder=1,
                label='Optical model trace BEFORE x-corr (x_model)')
    plt.scatter(x_model[ind], yref_x_model[ind], marker='o', s=40, facecolors='none', edgecolors='g', zorder=1,
                label='Optical model trace AFTER x-corr (x_model[ind])')

    # Print in the plot the values of `slintindx` for the matched edges
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
