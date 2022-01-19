"""
This script explores the noise in a slit or spectrum
for one of the outputs of PypeIt

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os

import numpy as np

import os
from matplotlib import pyplot as plt
from astropy.stats import sigma_clip, mad_std
from astropy.modeling.models import Gaussian1D
from astropy.io import fits
from astropy.table import Table

from pypeit.scripts import scriptbase
from pypeit import utils

from IPython import embed


def plot(args, line_wav_z, line_names, flux, err, mask, input_mask,
         ratio, lbda, filename:str, folder:str=None, z:float=None):

    ########### PLOT ######
    fig=plt.figure(figsize=(23,6.))
    ax=plt.subplot2grid((1, 4), (0, 0), rowspan=1, colspan=3)
    ax.minorticks_on()
    if args.step:
        drawstyle='steps-mid'
    else:
        drawstyle='default'
    sec=mask

    if lbda[sec].size ==0:
        sec = lbda >0
    ax.plot(lbda[sec], flux[sec], lw=1, zorder=1, drawstyle=drawstyle, label='{}'.format(filename))
    ax.plot(lbda[sec], np.zeros(lbda[sec].size), ls='--', color='Gray', zorder=-1)
    
    if args.ploterr: 
        ax.plot(lbda[sec], err[sec], drawstyle=drawstyle, lw=1, color='#ff6666', zorder=0, label='noise')
    if z is not None:
        for i in range(line_wav_z.shape[0]):
            if (line_wav_z[i]>lbda[sec].min())&(line_wav_z[i]<lbda[sec].max()):
                ax.axvline(line_wav_z[i], color='Gray', ls='dotted', zorder=-1)
                ax.annotate('{}'.format(line_names[i]), xy=(line_wav_z[i],1),xytext=(line_wav_z[i],0.95), xycoords=('data', 'axes fraction'), arrowprops=dict(facecolor='None', edgecolor='None', headwidth=0., headlength=0, width=0, shrink=0.), annotation_clip=True, horizontalalignment='center', color='k', fontsize=10)

    ax.axvspan(lbda[input_mask].min(), lbda[input_mask].max(), color='tab:green', alpha=0.2)
    ax.set_xlim(np.min(lbda[sec]), np.max(lbda[sec]))
    ymax = sigma_clip(flux[sec], sigma=10, return_bounds=True)[2]*1.3
    ymin = sigma_clip(flux[sec], sigma=5, return_bounds=True)[1]*1.05
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('Wavelength  (Angstrom)')
    ax.set_ylabel(r'Counts')
    plt.legend(loc=2)
    ax2=plt.subplot2grid((1, 4), (0, 3), rowspan=1, colspan=1)
    ax2.minorticks_on()
    bins=np.arange(ratio.min(), ratio.max(), 0.1)
    hist_n, hist_bins, _ = ax2.hist(ratio, bins=bins, histtype='stepfilled')
    mod_mods=Gaussian1D(amplitude=hist_n.max(), mean=np.median(ratio), 
                        stddev=1.)
    ax2.plot(bins, mod_mods(bins), label=r"Gaussian ($\sigma=1$)")
    ax2.axvline(0, ls='dotted', color='Gray')
    ax2.set_xlim(hist_bins[:-1][hist_n > 10].min()*2, 
                 hist_bins[:-1][hist_n > 10].max()*2)
    ax2.set_ylim(-0.02, hist_n.max()*1.5)
    ax2.set_xlabel(r'Flux/Noise')
    ax2.set_ylabel(r'#')
    err_over_flux = (np.median(err[input_mask])/mad_std(flux[input_mask]))
    ax2.text(0.99, 0.95, r'Median Noise= {:.1f} - Flux RMS= {:.1f} --> {:.2f}x'.format(np.median(err[input_mask]), mad_std(flux[input_mask]), err_over_flux), color='k', fontsize=9, horizontalalignment='right', transform=ax2.transAxes)
    ax2.text(0.99, 0.90, r'Chi:  Median = {:.2f}, Std = {:.2f}'.format(np.median(ratio), mad_std(ratio)), color='k', fontsize=12, horizontalalignment='right', transform=ax2.transAxes, weight='bold')
    ax2.legend(loc=2)
    plt.tight_layout()
    if args.plot_or_save == 'plot': 
        plt.show()
    if args.plot_or_save == 'save': 
        plt.savefig('{}/noisecheck_{}.png'.format(folder, filename), 
                    bbox_inches='tight', dpi=400)
    plt.close()


class ChkNoise1D(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Examine the noise in a PypeIt spectrum',
                                    width=width)
        parser.add_argument('files', type = str, nargs='*', help = 'PypeIt spec1d file(s)')
        parser.add_argument('--fileformat', default='spec1d', type=str, help='Is this coadd1d or spec1d?')
        parser.add_argument('--extraction', default='opt', type=str, help='If spec1d, which extraction? OPT or BOX')
        parser.add_argument('--ploterr', default=False, help='Plot noise spectrum',action='store_true')
        parser.add_argument('--step', default=False, help='Use `steps-mid` as linestyle',action='store_true')
        parser.add_argument('--z', default=None, type=float, nargs='*', help='Object redshift')
        parser.add_argument('--maskdef_objname', default=None, type=str, help='MASKDEF_OBJNAME of the target that you want to plot')
        parser.add_argument('--pypeit_name', default=None, type=str, help='PypeIt name of the target that you want to plot')
        parser.add_argument('--wavemin', default=None, type=float, help='Wavelength min. This is for selecting a region of the spectrum to analyze.')
        parser.add_argument('--wavemax', default=None, type=float, help='Wavelength max.This is for selecting a region of the spectrum to analyze.')
        parser.add_argument('--plot_or_save', default='plot', type=str, help='Do you want to save to disk or open a plot in a mpl window. If you choose save, a folder called spec1d*_noisecheck will be created and all the relevant plot will be placed there.')
        args = parser.parse_args()

        return parser

    @staticmethod
    def main(args):

        # Load em
        line_names, line_wav = utils.list_of_spectral_lines()
            
        files=np.array(args.files)

        if args.z is not None:
            zs = np.array(args.z) 

        # Loop on the files
        for i in range(files.size):    
            # reinitialize lines wave
            line_wav_z=line_wav.copy()

            file = files[i]

            # Deal with redshifts
            if args.z is not None:
                z = zs[i] if zs.size == files.size else zs[0]
                line_wav_z *= (1+z)   #redshift linelist
            else:
                z = None

            # Save?
            if args.fileformat == 'spec1d':
                if args.plot_or_save == 'save':
                    folder = '{}_noisecheck'.format(file.split('.fits')[0])
                    if not os.path.exists(folder): os.makedirs(folder)

            # I/O
            hdu= fits.open(file)
            hdr = hdu[0].header

            for i in range(len(hdu)-1):
                hdr_spec1d = hdu[i+1].header
                data = Table(hdu[i+1].data)
                if args.fileformat == 'coadd1d':
                    lbda = data['wave'].data
                    flux = data['flux'].data     # counts
                    err = 1./np.sqrt(data['ivar'])
                    mask = data['mask'] == 1
                    filename = file.split('.fits')[0]
                    maskdef_objname = None
                    pypeit_name = None
                if (args.fileformat == 'spec1d') and (hdr_spec1d.get('EXTNAME') is not None):
                    if hdr_spec1d.get('EXTNAME')[:4] != 'SPAT':
                        pass
                    else:
                        if args.extraction == 'box':
                            lbda = data['BOX_WAVE'].data
                            flux = data['BOX_COUNTS'].data     # counts
                            err = data['BOX_COUNTS_SIG']
                            mask = data['BOX_MASK']
                            extraction = 'BOX'
                        else:
                            try:
                                lbda = data['OPT_WAVE'].data
                                flux = data['OPT_COUNTS'].data     # counts
                                err = data['OPT_COUNTS_SIG']
                                mask = data['OPT_MASK']
                                extraction = 'OPT'
                            except KeyError:
                                lbda = data['BOX_WAVE'].data
                                flux = data['BOX_COUNTS'].data     # counts
                                err = data['BOX_COUNTS_SIG']
                                mask = data['BOX_MASK']
                                extraction = 'BOX'

                        maskdef_objname = hdr_spec1d.get('HIERARCH MASKDEF_OBJNAME')
                        pypeit_name = hdr_spec1d.get('NAME')
                    if maskdef_objname is not None:
                        if hdr_spec1d.get('MASKDEF_EXTRACT') is not False:
                            filename = '{}_{}obj{}_{}_{}'.format(hdr.get('DECKER'), extraction, maskdef_objname, pypeit_name, 'maskdef_extract')
                        else:
                            filename = '{}_{}obj{}_{}'.format(hdr.get('DECKER'),extraction, maskdef_objname, pypeit_name)
                    else:
                        filename = '{}_{}_{}'.format(hdr.get('DECKER'), extraction, pypeit_name)

                if (maskdef_objname is not None and maskdef_objname==args.maskdef_objname) or (pypeit_name is not None and pypeit_name==args.pypeit_name) or (args.maskdef_objname is None and args.pypeit_name is None):
                    input_mask = mask.copy()
                    if args.wavemin is not None:
                        input_mask &= lbda > args.wavemin
                    if args.wavemax is not None:
                        input_mask &= lbda < args.wavemax

                    if lbda[input_mask].size < 10:
                        pass
                    else:
                        ratio = flux[input_mask]/err[input_mask]
                        plot(args, line_wav_z, line_names, flux, err, mask, 
                             input_mask, ratio, lbda, filename, 
                             folder=None, z=None)