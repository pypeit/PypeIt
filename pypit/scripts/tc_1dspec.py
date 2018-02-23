#!/usr/bin/env python

import sys
import numpy as np
from pypit import arutils
from pypit import arload, arflux
from linetools.spectra.xspectrum1d import XSpectrum1D
from pkg_resources import resource_filename
#from pdb as debugger

#msgs = pyputils.get_dummy_logger()

def parser(options=None):
    import argparse

    description = (
                  'Script to telluric correct a spec!D file. '
                  'Currently only works for LRIS 800/10000 grating.'
                  )
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("infile", type=str, help="Input file (fits)")
    parser.add_argument("atm_tran", type=str, help="Atmospheric transmission spectrum")
    #parser.add_argument("--debug", default=False, i
    #                    action='store_true', help="Turn debugging on")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args

def get_transmission(atm_file, data):
    """ Convolve the atmospheric absorption model to the
    resolution and wavelengthg gird of the data.
    Parameters
    ----------
    atm_file: string
        filename passed through script call
    data: XSpectrum1d object
    Returns
    -------
    model: 2D ndarray
    """

    fname = resource_filename('pypit', 'data/extinction/{}'.format(atm_file))
    tmp = np.loadtxt(fname)

    # Convert transmission spectrum to angstroms
    # and clip it to wavelength range of data
    i = ((tmp[:,0]*10 >= min(data.wavelength.value)) &
         (tmp[:,0]*10 <= max(data.wavelength.value)))
    trans = XSpectrum1D.from_tuple((tmp[:,0][i]*10, tmp[:,1][i]))

    # Convolve the atm. transsissiom spectrum to LRIS_R resolution
    res         = arflux.lrisr_800_10000_r(trans.wavelength.value)
    fwhm_pix    = arutils.get_fwhm_pix(trans, inval=20000., outval=res)

    smooth_spec = np.zeros((len(trans.wavelength.value)))
    for i, fp in enumerate(fwhm_pix):
        tmp = trans.gauss_smooth(fp)
        smooth_spec[i] = tmp.flux.value[i]

    # Get model transmission on same wavelength grid
    # as data
    model = np.zeros((len(data.wavelength.value), 2))
    model[:,0] = data.wavelength.value
    model[:,1] = np.interp(data.wavelength.value,
                           trans.wavelength.value, smooth_spec)

    return model

def get_fscale(data, tran):
    """ Determine the best scale factor by minimizing the
    differences between the data and the transmission template
    Parameters
    ----------
    data: XSpectrum1D object
    tran: 2D ndarray
    Returns
    -------
    fscale: float
        The scale factor of bestfit template
    """
    from numpy.polynomial.chebyshev import chebfit, chebval
    from scipy.optimize import minimize

    i = ((data.wavelength.value >= 9250) & (data.wavelength.value <= 9650))

    coef = chebfit(data.wavelength.value[i], data.flux.value[i], 4)
    poly = chebval(data.wavelength.value[i], coef)
    tmp_data = np.zeros((len(data.wavelength.value[i]), 3))
    tmp_data[:,0] = data.wavelength.value[i]
    tmp_data[:,1] = data.flux.value[i]/poly
    tmp_data[:,2] = data.sig.value[i]

    coef = chebfit(tran[:,0][i], tran[:,1][i], 4)
    poly = chebval(tran[:,0][i], coef)
    tmp_tran = np.zeros((len(tran[:,0][i]), 2))
    tmp_tran[:,0] = tran[:,0][i]
    tmp_tran[:,1] = tran[:,1][i]/poly

    soln = minimize(arutils.opposite_lnlike_tf, 1,
                    args=(tmp_data, tmp_tran))

    return soln['x'][0]

# This function will eventuall write out a fits/hdf5 file
# and QA plots
def tcorrect_data(fscale, data, tran):
    from copy import deepcopy
    import matplotlib.pyplot as plt

    template = arutils.get_atm_template(fscale, tran[:,1])

    tc_data      = np.zeros((len(data.wavelength.value), 2))
    tc_data[:,0] = data.wavelength.value
    tc_data[:,1] = data.flux.value/template

    fig = plt.figure(figsize=(10,4))
    ax1 = plt.subplot(1,2,1)
    ax2 = plt.subplot(1,2,2)

    i = ((data.wavelength.value >= 8050) & (data.wavelength.value <= 8400))
    ax1.plot(data.wavelength.value[i], data.flux.value[i]/np.median(data.flux.value[i]) + 0.1,
             color='k', label='Data')
    ax1.plot(data.wavelength.value[i], template[i]/np.median(template[i]),
             color='r', label='Atm Model')
    ax1.plot(tc_data[:,0][i], tc_data[:,1][i]/np.median(tc_data[:,1][i]) - 0.12,
             color='b', label='TC Data')

    i = ((data.wavelength.value >= 8900) & (data.wavelength.value <= 9200))
    ax2.plot(data.wavelength.value[i], data.flux.value[i]/np.median(data.flux.value[i]) + 0.1,
             color='k', label='Data')
    ax2.plot(data.wavelength.value[i], template[i]/np.median(template[i]),
             color='r', label='Atm Model')
    ax2.plot(tc_data[:,0][i], tc_data[:,1][i]/np.median(tc_data[:,1][i]) - 0.12,
             color='b', label='TC Data')

    ax1.legend()

    plt.tight_layout()
    plt.show()

def main(args, unit_test=False, path=''):
    import glob
    import yaml
    import os.path

    # Load the input file
    data = arload.load_1dspec(args.infile, exten=1)
    tran = get_transmission(args.atm_tran, data)
    fscale = get_fscale(data, tran)
    tcorrect_data(fscale, data, tran)

