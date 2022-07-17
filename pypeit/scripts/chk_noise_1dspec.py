"""
This script explores the noise in a slit or spectrum
for one of the outputs of PypeIt

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import numpy as np

import os
from matplotlib import pyplot as plt
from astropy.stats import sigma_clip, mad_std
from astropy.modeling.models import Gaussian1D
from astropy.io import fits

from pypeit.scripts import scriptbase
from pypeit import utils
from pypeit import msgs
from pypeit import specobjs
from pypeit.onespec import OneSpec

from IPython import embed


def plot(args, line_wav_z:np.ndarray, line_names:np.ndarray,
         flux:np.ndarray, err:np.ndarray, mask:np.ndarray, 
         input_mask:np.ndarray, ratio:np.ndarray, 
         lbda:np.ndarray, basename:str, folder:str=None, z:float=None, plot_shaded=True):
    """ 
    Plot the 1D spectrum

    Args:
        args (`argparse.ArgumentParser`_):
            Command-line argument parser.
        line_wav_z (`numpy.ndarray`_):
            Array of wavelength of spectral features to be plotted.
        line_names (`numpy.ndarray`_):
            Array of names of spectral features to be plotted.
        flux (`numpy.ndarray`_):
            Flux of the spectrum.
        err (`numpy.ndarray`_):
            Sig of the spectrum.
        lbda (`numpy.ndarray`_):
            Wavelength values of the spectrum.
        mask (`numpy.ndarray`_):
        input_mask (`numpy.ndarray`_):
        ratio (`numpy.ndarray`_):
            Chi values
        basename (:obj:`str`):
            Basename
        folder (:obj:`str`, optional):
            Output folder, if saving. Defaults to None.
        z (:obj:`float`, optional):
            Redshift of the object that we want to plot. Defaults to None.
        plot_shaded (:obj:`bool`, optional):
            Do you want to plot the shaded area that shows the
            wavelength range used in the analysis? Defaults to True.
    """

    _folder = 'chk_noise_1dspec' if folder is None else folder

    fig = plt.figure(figsize=(23, 6.))
    ax = plt.subplot2grid((1, 4), (0, 0), rowspan=1, colspan=3)
    ax.minorticks_on()
    drawstyle = 'steps-mid' if args.step else 'default'

    sec=mask
    if lbda[sec].size == 0:
        sec = lbda > 0

    ax.plot(lbda[sec], flux[sec], lw=1, zorder=1, drawstyle=drawstyle, label='{}'.format(basename))
    ax.plot(lbda[sec], np.zeros(lbda[sec].size), ls='--', color='Gray', zorder=-1)
    
    if args.ploterr: 
        ax.plot(lbda[sec], err[sec], drawstyle=drawstyle, lw=1, color='#ff6666', zorder=2, label='noise')
    if z is not None:
        line_num = 0
        for i in range(line_wav_z.shape[0]):
            if lbda[sec].min() < line_wav_z[i] < lbda[sec].max():
                line_num += 1
                yannot = 0.94 if (line_num % 2) == 0 else 0.90
                ax.axvline(line_wav_z[i], color='Gray', ls='dotted', zorder=-1)
                ax.annotate('{}'.format(line_names[i]), xy=(line_wav_z[i],1), xytext=(line_wav_z[i], yannot),
                            xycoords=('data', 'axes fraction'), arrowprops=dict(facecolor='None', edgecolor='None',
                                                                                headwidth=0., headlength=0, width=0,
                                                                                shrink=0.),
                            annotation_clip=True, horizontalalignment='center', color='k', fontsize=10)
    if plot_shaded:
        ax.axvspan(lbda[input_mask].min(), lbda[input_mask].max(), color='tab:green', alpha=0.2)
    ax.set_xlim(np.min(lbda[sec]), np.max(lbda[sec]))
    ymax = sigma_clip(flux[sec], sigma=10, return_bounds=True)[2]*1.3
    ymin = sigma_clip(flux[sec], sigma=5, return_bounds=True)[1]*1.05
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('Wavelength  (Angstrom)')
    ax.set_ylabel(r'Counts')
    plt.legend(loc=3)
    ax2 = plt.subplot2grid((1, 4), (0, 3), rowspan=1, colspan=1)
    ax2.minorticks_on()
    bin_step = (np.nanmax(ratio) - np.nanmin(ratio))*0.01
    bins = np.arange(np.nanmin(ratio), np.nanmax(ratio), bin_step)
    hist_n, hist_bins, _ = ax2.hist(ratio, bins=bins, histtype='stepfilled')
    mod_mods = Gaussian1D(amplitude=hist_n.max(), mean=np.median(ratio), stddev=1.)
    ax2.plot(bins, mod_mods(bins), label=r"$\sigma=1$")
    ax2.axvline(0, ls='dotted', color='Gray')
    ax2.set_xlim(-6.5, 6.5)
    ax2.set_ylim(-0.02, hist_n.max()*1.5)
    ax2.set_xlabel(r'Flux/Noise')
    ax2.set_ylabel(r'#')
    err_over_flux = (np.median(err[input_mask])/mad_std(flux[input_mask]))
    ax2.text(0.97, 0.93, r'Median Noise= {:.1f} - Flux RMS= {:.1f} --> {:.2f}x'.format(np.median(err[input_mask]),
                                                                                       mad_std(flux[input_mask]),
                                                                                       err_over_flux),
             color='k', fontsize=9, horizontalalignment='right', transform=ax2.transAxes)
    ax2.text(0.97, 0.87, r'Chi:  Median = {:.2f}, Std = {:.2f}'.format(np.median(ratio), mad_std(ratio)),
             color='k', fontsize=12, horizontalalignment='right', transform=ax2.transAxes, weight='bold')
    ax2.legend(loc=2)
    plt.tight_layout()

    # Finish
    if args.plot_or_save == 'plot': 
        plt.show()
    elif args.plot_or_save == 'save': 
        plt.savefig('{}/noisecheck_{}.png'.format(_folder, basename), bbox_inches='tight', dpi=200)
    plt.close()


def get_basename(header, extract_type, pypeit_name=None, maskdef_objname=None, maskdef_extract=False):

    if maskdef_extract:
        return '{}_{}obj{}_{}_{}'.format(header.get('DECKER'), extract_type, maskdef_objname, pypeit_name, 'maskdef_extract')
    elif maskdef_objname is not None:
        return '{}_{}obj{}_{}'.format(header.get('DECKER'), extract_type, maskdef_objname, pypeit_name)
    elif pypeit_name is not None:
        return '{}_{}_{}'.format(header.get('DECKER'), extract_type, pypeit_name)
    else:
        return '{}_{}'.format(header.get('DECKER'), extract_type)


class ChkNoise1D(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Examine the noise in a PypeIt spectrum',
                                    width=width)
        parser.add_argument('files', type=str, nargs='*', help='PypeIt spec1d file(s)')
        parser.add_argument('--fileformat', default='spec1d', type=str, help='Is this coadd1d or spec1d?')
        parser.add_argument('--extraction', default='opt', type=str, help='If spec1d, which extraction? opt or box')
        parser.add_argument('--ploterr', default=False, help='Plot noise spectrum', action='store_true')
        parser.add_argument('--step', default=False, help='Use `steps-mid` as linestyle', action='store_true')
        parser.add_argument('--z', default=None, type=float, nargs='*', help='Object redshift')
        parser.add_argument('--maskdef_objname', default=None, type=str, help='MASKDEF_OBJNAME of the '
                                                                              'target that you want to plot. '
                                                                              'If maskdef_objname is not provided, '
                                                                              'nor a pypeit_name, all the 1D spectra '
                                                                              'in the file(s) will be plotted.')
        parser.add_argument('--pypeit_name', default=None, type=str, help='PypeIt name of the target that '
                                                                          'you want to plot. '
                                                                          'If pypeit_name is not provided, '
                                                                          'nor a maskdef_objname, all the 1D spectra '
                                                                          'in the file(s) will be plotted.')
        parser.add_argument('--wavemin', default=None, type=float, help='Wavelength min. This is for selecting '
                                                                        'a region of the spectrum to analyze.')
        parser.add_argument('--wavemax', default=None, type=float, help='Wavelength max.This is for selecting '
                                                                        'a region of the spectrum to analyze.')
        parser.add_argument('--plot_or_save', default='plot', type=str, help='Do you want to save to disk or open '
                                                                             'a plot in a mpl window. If you choose '
                                                                             'save, a folder called spec1d*_noisecheck'
                                                                             ' will be created and all the relevant '
                                                                             'plot will be placed there.')

        return parser

    @staticmethod
    def main(args):
        # Load em
        line_names, line_wav = utils.list_of_spectral_lines()
            
        files = np.array(args.files)

        if args.z is not None:
            zs = np.array(args.z) 

        # Loop on the files
        for i in range(files.size):    
            # reinitialize lines wave
            line_wav_z = line_wav.copy()

            # Deal with redshifts
            if args.z is not None:
                z = zs[i] if zs.size == files.size else zs[0]
                line_wav_z *= (1+z)   # redshifted linelist
            else:
                z = None

            # this file and header
            file = files[i]
            head = fits.getheader(file)

            # I/O spec object
            specObjs = [OneSpec.from_file(file)] if args.fileformat == 'coadd1d' else \
                            specobjs.SpecObjs.from_fitsfile(file, chk_version=False)

            # loop on the spectra
            for spec in specObjs:
                if args.fileformat == 'coadd1d':
                    lbda = spec.wave
                    flux = spec.flux
                    err = spec.sig
                    mask = spec.mask == 1
                    extract_type = spec.ext_mode
                    maskdef_objname = None
                    pypeit_name = None
                    maskdef_extract = False
                else:
                    maskdef_objname = spec.MASKDEF_OBJNAME
                    pypeit_name = spec.NAME
                    maskdef_extract = spec.MASKDEF_EXTRACT
                    if args.extraction == 'opt' and spec['OPT_COUNTS'] is not None:
                        lbda = spec['OPT_WAVE']
                        flux = spec['OPT_COUNTS']
                        err = spec['OPT_COUNTS_SIG']
                        mask = spec['OPT_MASK']
                        extract_type = 'OPT'
                    else:
                        lbda = spec['BOX_WAVE']
                        flux = spec['BOX_COUNTS']
                        err = spec['BOX_COUNTS_SIG']
                        mask = spec['BOX_MASK']
                        extract_type = 'BOX'

                # plot specific spectrum or all of them?
                plot_this = True
                if args.maskdef_objname is not None:
                    plot_this = True if maskdef_objname == args.maskdef_objname else False
                if args.pypeit_name is not None:
                    plot_this = True if pypeit_name == args.pypeit_name else False

                # OK plot
                if plot_this:
                    basename = get_basename(head, extract_type, pypeit_name=pypeit_name,
                                            maskdef_objname=maskdef_objname, maskdef_extract=maskdef_extract)
                    # Cut down
                    input_mask = mask.copy()
                    if args.wavemin is not None:
                        input_mask &= lbda > args.wavemin
                    if args.wavemax is not None:
                        input_mask &= lbda < args.wavemax

                    if lbda[input_mask].size < 10:
                        msgs.warn("The spectrum was cut down to <10 pixels.  Skipping")
                        continue

                    # determine if plotting the shaded area in the plot that shows the
                    # wavelength range used to compute the chi
                    plot_shaded = False if args.wavemin is None and args.wavemax is None else True

                    # Save?
                    if args.plot_or_save == 'save':
                        folder = '{}_noisecheck'.format(file.split('.fits')[0])
                        if not os.path.exists(folder):
                            os.makedirs(folder)
                    else:
                        folder = None

                    # Plot
                    ratio = flux[input_mask]/err[input_mask]
                    plot(args, line_wav_z, line_names, flux, err, mask, input_mask, ratio, lbda, basename,
                         folder=folder, z=z, plot_shaded=plot_shaded)

