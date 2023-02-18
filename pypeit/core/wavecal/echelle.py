"""
Module for echelle-specific wavelength calibration functions.

.. include:: ../include/links.rst
"""
from IPython import embed

import numpy as np
from scipy import interpolate

from astropy.io import fits
from astropy.table import Table

from pypeit import msgs
from pypeit.core import fitting
from pypeit.core.wavecal import wvutils
from pypeit import data


def predict_ech_order_coverage(angle_fits_params, xd_angle_coeffs, xdisp, xdangle, norders, pad=0):
    """
    Predict the coverage of orders in the echelle spectrum using the disperser dependent
    fits of the reddest order as a function of xdangle.

    Args:
        angle_fits_params (astropy.table.Table):
            Table holding the arxiv parameters
        xd_angle_coeffs
            Table holding the arxiv data
        xdisp (str):
            Corss disperser. For HIRES this is either 'UV' or 'RED'
        xdangle (float):
            Cross-disperser angle.
        norders (int):
            Number of orders identified on the detector
        pad (int):
            Number of orders to pad the coverage by on the blue and red side.

    Returns:
        `numpy.ndarray`_: Array of order numbers for the predicted coverage.
    """

    # Evaluate the fits for reddest order vs xdanalge which using hte values stored in the angle_fits_params
    xd_min, xd_max = angle_fits_params['xd_xmin'], angle_fits_params['xd_xmax']
    idisp = angle_fits_params['xdisp_vec'] == xdisp
    reddest_order_fit = int(np.round(
        fitting.evaluate_fit(xd_angle_coeffs[idisp, :].flatten(), angle_fits_params['xd_func'], xdangle,
                             minx=xd_min, maxx=xd_max)))
    order_vec = reddest_order_fit + (np.arange(norders + 2*pad) - pad)[::-1]

    return order_vec

def predict_ech_wave_soln(angle_fits_params, ech_angle_coeffs, ech_angle, order_vec, nspec):
    """
    Predict an echelle spectrum wavelength solution for each order by evluating the polynomial fits of
    wavelength solution coefficients vs echelle angle at the given echelle angle.

    Args:
        angle_fits_params (astropy.table.Table):
            Table holding the parameters governing the echelle angle fits
        ech_angle_coeffs (numpy.ndarray):
            Array holding the polynomial coefficients for the fits of the wavelength solution polynomial coefficients
            vs echelle angle.
        ech_angle (float):
            Echelle angle
        order_vec (numpy.ndarray):
            Array of order numbers for the deisred predicted spectrum. Shape = (norders,)
        nspec (int):
            Number of spectral pixels in the echelle spectrum

    Returns:
        `numpy.ndarray`_: Array containing the predicted echelle spectrum. Shape
        is (nspec, norders)
    """

    norders = order_vec.size
    wave_soln_guess = np.zeros((nspec, norders))

    xnspecmin1 = float(nspec - 1)
    xnspec = np.arange(nspec)/xnspecmin1

    for iord, order in enumerate(order_vec):
        # Index of the order in the total order vector used cataloguing the fits in the coeff arxiv
        indx = order - angle_fits_params['order_min']
        coeff_predict = np.zeros(angle_fits_params['ech_n_final'] + 1)
        # Evaluate the coefficients for this order and the current ech_angle
        for ic in range(angle_fits_params['ech_n_final'] + 1):
            coeff_predict[ic] = fitting.evaluate_fit(
                ech_angle_coeffs[indx, ic, :], angle_fits_params['ech_func'],
                ech_angle, minx=angle_fits_params['ech_xmin'], maxx=angle_fits_params['ech_xmax'])

        wave_soln_guess[:, iord] = fitting.evaluate_fit(coeff_predict, angle_fits_params['wave_func'], xnspec,
        minx=angle_fits_params['wave_xmin'], maxx=angle_fits_params['wave_xmax'])


    return wave_soln_guess


def predict_ech_arcspec(angle_fits_file, composite_arc_file, echangle, xdangle, xdisp, nspec, norders, pad=3):
    """
    Predict the echelle arc spectrum using the fits to wavelength solution vs echangle and xdangle  and the archived
    composite arcs.

    Parameters
    ----------
    angle_fits_file : str
        File containing the fits to wavelength solution vs echangle and xdangle
    composite_arc_file : str
        File containing the archived composite arcs for each order.
    echangle : float
        Echelle angle
    xdangle : float
        Cross-disperser angle
    xdisp : str
        Cross disperser. E.g. for Keck HIRES this is either 'UV' or 'RED'
    nspec : int
        Number of spectral pixels in the echelle spectrum
    norders : int
        Number of orders in the echelle spectrum
    pad : int
        Number of orders to pad the coverage by on the blue and red side.

    Returns
    -------
    order_vec_guess : `numpy.ndarray`_
        Vector of order numbers for the predicted echelle spectrum. Shape = (norders,)
    wave_soln_guess :  `numpy.ndarray`_
        Predicted wavelength solution. Shape = (nspec, norders)
    arcspec_guess :  `numpy.ndarray`_
        Predicted echelle arc spectrum. Shape = (nspec, norders)

    """

    # Read in the echelle angle fits
    angle_fits_file, _ = data.get_reid_arxiv_filepath(angle_fits_file)
    hdu = fits.open(angle_fits_file)
    angle_fits_params = Table(hdu[1].data)[0]
    ech_angle_coeffs = hdu[2].data
    xd_angle_coeffs = hdu[3].data

    # Read in the composite arc spectrum
    composite_arc_file, _ = data.get_reid_arxiv_filepath(composite_arc_file)
    hdu = fits.open(composite_arc_file)
    composite_arc_params = Table(hdu[1].data)[0]
    wave_composite = hdu[2].data
    arc_composite = hdu[3].data
    gpm_composite = (hdu[4].data).astype(bool)

    order_vec_guess = predict_ech_order_coverage(angle_fits_params, xd_angle_coeffs, xdisp, xdangle, norders, pad=pad)
    norders_guess = order_vec_guess.size
    wave_soln_guess = predict_ech_wave_soln(angle_fits_params, ech_angle_coeffs, echangle, order_vec_guess, nspec)


    order_min, order_max = angle_fits_params['order_min'], angle_fits_params['order_max']

    arcspec_guess = np.zeros_like(wave_soln_guess)
    # Interpolate the composite arc spectrum onto the predicted wavelength solution
    for iord, order in enumerate(order_vec_guess):
        indx = order - order_min
        igood = gpm_composite[:, indx]
        arcspec_guess[:, iord] = interpolate.interp1d(wave_composite[igood, indx], arc_composite[igood, indx],
                                                      kind='cubic', bounds_error=False,
                                                      fill_value=-1e10)(wave_soln_guess[:, iord])


    return order_vec_guess, wave_soln_guess, arcspec_guess

def identify_ech_orders(arcspec, echangle, xdangle, dispname, angle_fits_file, 
                        composite_arc_file, debug=False, pad=3):
    """
    Identify the orders in the echelle spectrum via cross correlation with the best guess predicted arc based
    on echangle, xdangle, and cross-disperser

    Parameters
    ----------
    arcspec : `numpy.ndarray`_
        Extracted arc spectrum, shape = (nspec, norders)
    echangle : float
        Echelle angle
    xdangle : float
        Cross-disperser angle
    dispname : str
        Cross-disperser. E.g. for Keck HIRES this is either 'UV' or 'RED'
    angle_fits_file : str
        File containing the fits to wavelength solution vs echangle and xdangle
    composite_arc_file : str
        File containing the archived composite arcs for each order.
    pad : int, optional
        Number of orders to pad the coverage by on the blue and red side.
    debug : bool, optional
        Passed to xcorr_shift

    Returns
    -------
    order_vec : `numpy.ndarray`_
        Array of order numbers corresponding to the input arcspec, shape = (norders,)
    wave_soln_guess : `numpy.ndarray`_
        Array containing the predicted wavelength solution, shape = (nspec, norders)
    arcspec_guess : `numpy.ndarray`_
        Array containing the predicted arc spectrum, shape = (nspec, norders)

    """

    nspec, norders = arcspec.shape

    # Predict the echelle order coverage and wavelength solution
    order_vec_guess, wave_soln_guess_pad, arcspec_guess_pad = predict_ech_arcspec(
        angle_fits_file, composite_arc_file, echangle, xdangle, dispname, nspec, norders, pad=pad)
    norders_guess = order_vec_guess.size

    # Since we padded the guess we need to pad the data to the same size
    arccen_pad = np.zeros((nspec, norders_guess))
    arccen_pad[:nspec, :norders] = arcspec

    # Cross correlate the data with the predicted arc spectrum
    # TODO Does it make sense for xcorr_shift to continuum subtract here?
    shift_cc, corr_cc = wvutils.xcorr_shift(
        arccen_pad.flatten('F'), arcspec_guess_pad.flatten('F'), 
        percent_ceil=50.0, sigdetect=5.0, sig_ceil=10.0, fwhm=4.0, debug=debug)
    
    # Finish
    x_ordr_shift = shift_cc / nspec
    ordr_shift = int(np.round(shift_cc / nspec))
    spec_shift = int(np.round(shift_cc - ordr_shift * nspec))
    msgs.info('Shift in order number between prediction and reddest order: {:.3f}'.format(ordr_shift + pad))
    msgs.info('Shift in spectral pixels between prediction and data: {:.3f}'.format(spec_shift))

    order_vec = order_vec_guess[-1] - ordr_shift + np.arange(norders)[::-1]
    ind = np.isin(order_vec_guess, order_vec, assume_unique=True)


    return order_vec, wave_soln_guess_pad[:, ind], arcspec_guess_pad[:, ind]


