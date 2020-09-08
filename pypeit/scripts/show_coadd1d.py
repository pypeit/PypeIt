#!/usr/bin/env python
"""
Wrapper to the linetools XSpecGUI
"""
import argparse
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import gridspec

from pypeit import coadd1d
from pypeit.core import coadd
import numpy as np
from pypeit import utils

def rebinspec(wave, flux, sig, gpm=None, nbin=3, wave_method='linear', wave_grid_min=None, wave_grid_max=None,
              dwave=None, dv=None, dloglam=None):

    if gpm is None:
        gpm = sig>0.
    # Get the new wave_grid
    wave_grid, wave_grid_mid, dsamp = coadd.get_wave_grid(wave, masks=gpm, wave_method=wave_method,
                                                          wave_grid_min=wave_grid_min, wave_grid_max=wave_grid_max,
                                                          dwave=dwave, dv=dv, dloglam=dloglam, samp_fact=1.0/nbin)

    # Using the ivar weights
    ivar = utils.calc_ivar(sig**2)
    weight = np.copy(ivar)

    # Do the rebin with histogram algorithm
    wave_binned, flux_binned, ivar_binned, gpm_binned, nused = coadd.compute_stack(wave_grid, wave, flux, ivar, gpm,
                                                                                   weight)
    sig_binned = np.sqrt(utils.calc_ivar(ivar_binned))
    indsort = np.argsort(wave_binned)
    wave_binned, flux_binned, ivar_binned, gpm_binned, nused = wave_binned[indsort], flux_binned[indsort],\
                                                               ivar_binned[indsort], gpm_binned[indsort], nused[indsort]

    gpm_wave = (wave_binned>np.min(wave_grid)) & (wave_binned<np.max(wave_grid))
    return wave_binned[gpm_wave], flux_binned[gpm_wave], sig_binned[gpm_wave], gpm_binned[gpm_wave]

def parser(options=None):
    parser = argparse.ArgumentParser(description='Parse')
    parser.add_argument("file", type=str, help="Spectral file")
    parser.add_argument("-b", "--rebin", type=int, default=0, help="rebin the spectrum by how many pixels")
    parser.add_argument("-o", "--outfile", type=str, default=None, help="output file name of the plot")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args

def main(args):
    """ Runs the XSpecGui on an input file
    """

    onespec = coadd1d.OneSpec.from_file(args.file)

    # Unpack
    wave = onespec.wave
    flux = onespec.flux
    telluric = onespec.telluric
    model = onespec.obj_model #* telluric
    ivar = onespec.ivar
    gpm = onespec.mask.astype(bool)
    #spect_dict = onespec.spect_meta
    head = onespec.head0
    sig = np.sqrt(utils.inverse(ivar))

    if (args.rebin!=0) and (telluric is not None):
        tck_model = interpolate.splrep(wave, model, s=0)
        tck_telluric = interpolate.splrep(wave, telluric, s=0)

        if 'Echelle' in head['PYPELINE']:
            wave, flux, sig, gpm = rebinspec(wave, flux, sig, gpm=gpm, nbin=args.rebin, wave_method='velocity')
        else:
            wave, flux, sig, gpm = rebinspec(wave, flux, sig, gpm=gpm, nbin=args.rebin, wave_method='linear')

        model = interpolate.splev(wave, tck_model, der=0)
        telluric = interpolate.splev(wave, tck_telluric, der=0)


    ymax = 1.5*np.percentile(flux[gpm],95)
    ymin = np.fmin(0.,np.fmax(np.percentile(flux[gpm],5),-0.1*ymax))

    plt.rcParams["xtick.top"] = True
    plt.rcParams["ytick.right"] = True
    plt.rcParams["xtick.minor.visible"] = True
    plt.rcParams["ytick.minor.visible"] = True
    plt.rcParams["ytick.direction"] = 'in'
    plt.rcParams["xtick.direction"] = 'in'
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.default"] = "regular"
    plt.rcParams["text.usetex"] = True

    if telluric is not None:
        f = plt.figure(figsize=(10, 5))
        f.subplots_adjust(left=0.09, right=0.97, bottom=0.12, top=0.95, wspace=0, hspace=0)
        gs = gridspec.GridSpec(2, 1, hspace=0, wspace=0.2)
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
        ax1.plot(wave,telluric*0.9*ymax,linestyle='-',color='0.7',lw=1.0)
        ax1.plot(wave,flux,'k-',lw=0.5)
        ax1.plot(wave,sig,'-',color='steelblue',lw=0.5)
        ax1.set_xlim(wave[gpm].min(),wave[gpm].max())
        ax1.set_ylim(ymin,ymax)

        ax2.plot(wave, flux*telluric,'k-',lw=0.5)
        ax2.plot(wave, model*telluric,'-',color='darkorange',lw=1)
        ax2.set_xlim(wave[gpm].min(),wave[gpm].max())
        ax2.set_ylim(ymin,ymax)
        ax2.set_xlabel(r'Wavelength ($\rm \AA$)', fontsize=14)
        ax2.set_ylabel(r'$f_{\rm \lambda} ~ {\rm (10^{-17}~erg~s^{-1}~cm^{-2}~\AA^{-1})}$', fontsize=14)
        ax2.yaxis.set_label_coords(-0.06, 1.0)

    else:
        f = plt.figure(figsize=(10, 5))
        f.subplots_adjust(left=0.09, right=0.97, bottom=0.12, top=0.95, wspace=0, hspace=0)
        plt.plot(wave,flux,'k-',lw=0.5)
        plt.plot(wave,sig,':',color='lightskyblue',lw=0.5)
        plt.xlim(wave[gpm].min(),wave[gpm].max())
        plt.ylim(ymin,ymax)
        plt.xlabel(r'Wavelength ($\rm \AA$)', fontsize=14)
        plt.ylabel(r'$f_{\rm \lambda} ~ {\rm (10^{-17}~erg~s^{-1}~cm^{-2}~\AA^{-1})}$', fontsize=14)

    if args.outfile is not None:
        if args.outfile[-3:] not in ['png','pdf','jpg']:
            outfile = args.outfile+'.png'
        else:
            outfile = args.outfile
        plt.savefig(outfile, dpi=300)

    plt.show()

