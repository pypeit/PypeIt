#!/usr/bin/env python

import sys
import numpy as np
from astropy.io import fits
from pypit import arutils
from pypit import arload, arflux
from linetools.spectra.xspectrum1d import XSpectrum1D
#from pdb as debugger

#msgs = pyputils.get_dummy_logger()

def parser(options=None):
    import argparse

    description = (
                  'Script to telluric correct a '
                  'spec1D file. Currently only works '
                  'for LRIS 600/10000 grating.'
                  )
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("infile", type=str,
                        help="Input file (YAML)")
    #parser.add_argument("--debug", default=False, i
    #                    action='store_true', help="Turn debugging on")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args

def get_fscale(data, tran, l1=9250, l2=9650):
    """ Determine the best scale factor by minimizing the
    differences between the data and the transmission template
    Parameters
    ----------
    data: XSpectrum1D object
    tran: 2D ndarray
    l1: int
        Lower limit of range over which to do the normalization
    l2: int
        Upper limit of range over which to do the normalization
    Returns
    -------
    fscale: float
        The scale factor of bestfit template
    """
    from numpy.polynomial.chebyshev import chebfit, chebval
    from scipy.optimize import minimize

    i = ((data.wavelength.value >= l1) &
         (data.wavelength.value <= l2))

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

def tcorrect_data(fscale, data, tran, original):
    """ Correct the full spectrum
    Parameters
    ----------
    fscale: float
    data: XSpectrum1d object
    tran: 2D ndarray
    original: fits
        The original fits file
    Returns
    -------
    """
    #from pkg_resources import resource_filename
    import matplotlib.pyplot as plt

    template = arutils.get_atm_template(fscale, tran[:,1])

    tc_data      = np.zeros((len(data.wavelength.value), 2))
    tc_data[:,0] = data.wavelength.value
    tc_data[:,1] = data.flux.value/template

    # Need to copy original header in here
    with fits.open(original) as hdu:
        header = hdu[0].header
    primary_hdu = fits.PrimaryHDU(header=header)

    c1 = fits.Column(name='wave', array=tc_data[:,0], format='K')
    c2 = fits.Column(name='o_flux', array=data.flux.value, format='K')
    c3 = fits.Column(name='tc_flux', array=tc_data[:,1], format='K')
    c4 = fits.Column(name='error', array=data.sig.value, format='K')
    t  = fits.BinTableHDU.from_columns([c1, c2, c3, c4])

    # Save to directory that script is being run in
    t.writeto(original.strip('.fits') + '_tc.fits')


def main(args, unit_test=False, path=''):
    import glob
    import yaml
    import os.path

    with open(args.infile, 'r') as infile:
        tc_dict = yaml.load(infile)
        print(tc_dict)
        try:
            assert (tc_dict['grating'] == '600/10000')
        except AssertionError:
            error_message = ('Telluric correction only '
                             'valid for LRIS 600/10000 '
                             'grating currently. Please '
                             'open an issue on GitHub for '
                             'other gratings. \n '
                            )
            print(error_message)
            raise

        atm_tran = tc_dict['transmission']
        files    = tc_dict['filenames']
        if 'region' in tc_dict.keys():
            if tc_dict['region'] is not None:
                l1 = tc_dict['region'][0]
                l2 = tc_dict['region'][1]

    for fname in files.keys():
        exten  = tc_dict['filenames'][fname]
        data   = arload.load_1dspec(fname, exten=exten)
        tran   = arflux.get_transmission(atm_tran, data)
        try:
            assert (l1 in locals() and l2 in locals())
            fscale = get_fscale(data, tran, l1=l1, l2=l2)
        except:
            fscale = get_fscale(data, tran)
        tcorrect_data(fscale, data, tran, fname)

