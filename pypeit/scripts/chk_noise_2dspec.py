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
from astropy.table import Table
import argparse

from pypeit import spec2dobj
from pypeit import msgs
import pypeit
from pypeit.scripts import scriptbase
from pypeit.images.detector_container import DetectorContainer

from IPython import embed

def grab_lines():
    ## spectral features
    Lyalphanam, Lyalphawav='Lyalpha', 1215.7
    OIInam, OIIwav='[OII]', 3727.5 #average between 3726,3729
    OIIInam1, OIIIwav1='[OIII]', 5007.
    OIIInam2, OIIIwav2='[OIII]', 4959.
    OIIInam3, OIIIwav3='[OIII]', 4363.
    Halphanam, Halphawav='Halpha', 6563.
    Hbetanam, Hbetawav='Hbeta', 4861.
    Hdeltanam, Hdeltawav='Hdelta', 4101.
    Hgammanam, Hgammawav='Hgamma', 4341.

    NeIIInam, NeIIIwav = '[NeIII]', 3869.
    NeVnam, NeVwav = '[NeV]', 3426.
    SIInam, SIIwav = '[SII]', 6724. #average between 6717,6731


    ##absorption 
    H13nam, H13wav = 'H13', 3734.
    H12nam, H12wav = 'H12', 3750.
    H11nam, H11wav = 'H11', 3771.
    H10nam, H10wav = 'H10', 3798.
    H9nam, H9wav = 'H9', 3835.
    H8nam, H8wav = 'H8', 3889.
    HeInam, HeIwav = 'HeI', 3889.

    CAII_Knam, CaII_Kwav = 'CaK', 3934. 
    CAII_Hnam, CaII_Hwav = 'CaH', 3968.

    Gbandnam, Gbandwav = 'Gband', 4305.

    line_names=np.array([Lyalphanam, OIInam, OIIInam1, OIIInam2, OIIInam3, Halphanam, Hbetanam, Hdeltanam, Hgammanam, NeIIInam, NeVnam,SIInam, H13nam, H12nam, H11nam, H10nam, 
                            H9nam, HeInam, CAII_Knam, CAII_Hnam, Gbandnam,])

    line_wav=np.array([Lyalphawav, OIIwav, OIIIwav1, OIIIwav2, OIIIwav3, Halphawav, Hbetawav,Hdeltawav, Hgammawav, NeIIIwav, NeVwav, SIIwav, H13wav,H12wav, H11wav, H10wav, 
                        H9wav, HeIwav,CaII_Kwav, CaII_Hwav, Gbandwav])
    return line_names, line_wav


def plot(image:np.ndarray, line_wav:list, line_names:list, 
         lbda:np.ndarray, lbda_min:float, lbda_max:float, aspect_ratio, 
         chi_select, flux_select, err_select, filename:str):
    fig=plt.figure(figsize=(23,4.))
    ax=plt.subplot2grid((1, 4), (0, 0), rowspan=1, colspan=3)
    ax.minorticks_on()
    zmax = sigma_clip(image[image!=0], sigma=2, return_bounds=True)[2]*1.3
    zmin = sigma_clip(image[image!=0], sigma=2, return_bounds=True)[1]*1.3
    ax.imshow(image.T, origin ='lower', interpolation='nearest', aspect=aspect_ratio, vmin=zmin, vmax=zmax, cmap=plt.get_cmap('gist_gray'))
    ax.text(0.005, 1.1, '{}'.format(filename), color='k', fontsize=13, horizontalalignment='left', transform=ax.transAxes, bbox=dict(edgecolor='black', facecolor='white', linewidth=1))
    if line_wav.size>0:
        for i in range(line_wav.size):
            ax.axvline(line_wav[i], color='black', ls='dotted', zorder=2)
            ax.annotate('{}'.format(line_names[i]), 
                        xy=(line_wav[i],1),
                        xytext=(line_wav[i],1.03), 
                        xycoords=('data', 'axes fraction'), 
                        arrowprops=dict(facecolor='None', edgecolor='None', headwidth=0., headlength=0, width=0, shrink=0.), annotation_clip=True, horizontalalignment='center', color='k', fontsize=10)
    #ax.axvspan(lbda.searchsorted(lbda_min), lbda.searchsorted(lbda_max), 
    #           color='tab:green', alpha=0.2, zorder=2)
    ax.set_xticks([])
    ax.set_yticks([])

    # Guassian stats
    ax2=plt.subplot2grid((1, 4), (0, 3), rowspan=1, colspan=1)
    ax2.minorticks_on()

    bins=np.arange(chi_select[chi_select!=0].min(), 
                    chi_select[chi_select!=0].max(), 0.1)
    hist_n, hist_bins, _ = ax2.hist(chi_select[chi_select!=0], bins=bins, histtype='stepfilled')
    mod_mods=Gaussian1D(amplitude=hist_n.max(), mean=np.median(chi_select[chi_select!=0]), stddev=1.)
    ax2.plot(bins, mod_mods(bins), label=r"Gaussian ($\sigma=1$)")
    ax2.axvline(0, ls='dotted', color='Gray')
    ax2.set_xlim(hist_bins[:-1][hist_n > 10].min(), hist_bins[:-1][hist_n > 10].max())
    ax2.set_ylim(-0.02, hist_n.max()*1.5)
    ax2.set_xlabel(r'(sciimg - skymodel) * sqrt(ivarmodel) * (bpmmask == 0)')
    ax2.set_ylabel(r'#')
    err_over_flux = (np.median(err_select[flux_select!=0])/mad_std(flux_select[flux_select!=0]))
    ax2.text(0.99, 0.95, r'Median Noise= {:.1f} - Flux RMS= {:.1f} --> {:.2f}x'.format(np.median(err_select[flux_select!=0]), mad_std(flux_select[flux_select!=0]), err_over_flux), color='k', fontsize=9, horizontalalignment='right', transform=ax2.transAxes)
    ax2.text(0.99, 0.90, r'Chi:  Median = {:.2f}, Std = {:.2f}'.format(
        np.median(chi_select[chi_select!=0]), 
        mad_std(chi_select[chi_select!=0])), color='k', fontsize=12, horizontalalignment='right', transform=ax2.transAxes, weight='bold')
    ax2.legend(loc=2)
    plt.tight_layout()


class ChkNoise2D(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Examine the noise in a PypeIt slit/order',
                                    width=width)

        parser.add_argument('files', type = str, nargs='*', help = 'PypeIt spec2d file(s)')
        parser.add_argument('--det', default='1', type=str,
                            help='Detector name or number.  If a number, the name is constructed '
                                 'assuming the reduction is for a single detector.  If a string, '
                                 'it must match the name of the detector object (e.g., DET01 for '
                                 'a detector, MSC01 for a mosaic).')
        parser.add_argument('--z', default=None, type=float, nargs='*', help='Object redshift')
        parser.add_argument('--maskdef_id', default=None, type=int, help='MASKDEF_ID of the slit that you want to plot')
        parser.add_argument('--pypeit_id', default=None, type=int, help='PypeIt ID of the slit that you want to plot')
        parser.add_argument('--pad', default=-5, type=int, help='[spec2d only] Padding for the selected slit. Negative value will trim. [default: -5]')
        parser.add_argument('--aspect_ratio', default=3, type=int, help='Aspect ratio when plotting the spec2d')
        parser.add_argument('--wavemin', default=None, type=float, help='Wavelength min. This is for selecting a region of the spectrum to analyze.')
        parser.add_argument('--wavemax', default=None, type=float, help='Wavelength max.This is for selecting a region of the spectrum to analyze.')
        parser.add_argument('--mode', default='plot', type=str, help='Do you want to save to disk or open a plot in a mpl window. If you choose save, a folder called spec2d*_noisecheck will be created and all the relevant plot will be placed there.')
        #parser.add_argument('--det', default=1, type=int, help='Detector number')
        return parser



    @staticmethod
    def main(args):

        # Parse the detector name
        try:
            det = int(args.det)
        except:
            detname = args.det
        else:
            detname = DetectorContainer.get_name(det)

        # Load em
        line_names, line_wav = grab_lines()
            
        files=np.array(args.files)

        if args.z is not None:
            zs = np.array(args.z) 

        # Loop on the files
        for i in range(files.size):    
            # reinitialize lines wave
            line_wav_z=line_wav.copy()

            # Load 2D object
            file = files[i]
            spec2DObj = spec2dobj.Spec2DObj.from_file(file, detname, chk_version=False)

            # Deal with redshifts
            if args.z is not None:
                z = zs[i] if zs.size == files.size else zs[0]
                line_wav_z *= (1+z)   #redshift linelist
            else:
                z = None

            # Save?
            if args.mode == 'save':
                folder = '{}_noisecheck'.format(file.split('.fits')[0])
                if not os.path.exists(folder): os.makedirs(folder)
            elif args.mode == 'print':
                # Generate a Table for pretty printing
                tbl = Table()
                tbl['Slit'] = spec2DObj.slits.slitord_id
                tbl['med_chis'] = spec2DObj.med_chis
                tbl['std_chis'] = spec2DObj.std_chis
                print('')
                print(tbl)
                print('-----------------------------------------------------')
                return

            # Generate chi image
            chi = (spec2DObj.sciimg - spec2DObj.skymodel) * np.sqrt(spec2DObj.ivarmodel) * (spec2DObj.bpmmask == 0)

            # Find the slit of interest
            all_maskdef_ids = spec2DObj.slits.maskdef_id
            all_pypeit_ids = spec2DObj.slits.slitord_id
            if args.maskdef_id is not None and all_maskdef_ids is None:
                msgs.error('This spec2d does not have maskdef_id. Choose a pypeit_id insteed.')

            # Build the mask
            input_mask = spec2DObj.bpmmask == 0
            if args.wavemin is not None:
                input_mask *= spec2DObj.waveimg > args.wavemin
            if args.wavemax is not None:
                input_mask *= spec2DObj.waveimg < args.wavemax

            # Decide on slits to show
            if args.pypeit_id is not None or args.maskdef_id is not None:
                if args.maskdef_id is not None and args.maskdef_id in all_maskdef_ids:
                    pypeit_id = all_pypeit_ids[all_maskdef_ids==args.maskdef_id][0]
                    slitidx = np.where(all_maskdef_ids==args.maskdef_id)[0][0]
                elif args.pypeit_id is not None and args.pypeit_id in all_pypeit_ids:
                    pypeit_id = args.pypeit_id
                    slitidx = np.where(all_pypeit_ids==args.pypeit_id)[0][0]
                show_slits = range(slitidx, slitidx+1)
            else:
                show_slits = range(all_pypeit_ids.size)

                '''
                # Chi
                chi_slit, _, _ = spec2DObj.calc_chi_slit(slitidx, pad=args.pad)

                # Cut down on bad pixels and wavelengths (optional)
                chi_select = chi_slit * input_mask
                if np.all(chi_select == 0):
                    continue

                # Flux in slit
                flux_select = (spec2DObj.sciimg - spec2DObj.skymodel) * input_mask
                err_select = 1/np.sqrt(spec2DObj.ivarmodel)* input_mask

                left, right, _ = spec2DObj.slits.select_edges()
                spat_start = int(left[:, slitidx].min())
                spat_end = int(right[:, slitidx].max())

                lbda = spec2DObj.waveimg[:,pypeit_id]
                if lbda[lbda!=0].size == 0:
                    continue

                line_wav_plt = np.array([])
                line_names_plt = np.array([])
                if z is not None:
                    for i in range(line_wav_z.shape[0]):
                        if (line_wav_z[i]>lbda[lbda!=0].min())&(line_wav_z[i]<lbda[lbda!=0].max()):
                            line_wav_plt = np.append(line_wav_plt, lbda.searchsorted(line_wav_z[i]))
                            line_names_plt = np.append(line_names_plt, line_names[i])

                lbda_min = args.wavemin if args.wavemin is not None else lbda[lbda!=0].min()
                lbda_max = args.wavemax if args.wavemax is not None else lbda[lbda!=0].max()

                # Plot!
                plot(chi_slit[:, spat_start:spat_end], 
                          line_wav_plt, line_names_plt, lbda,
                          lbda_min, lbda_max, args.aspect_ratio, 
                          chi_select, flux_select, err_select, filename)
                if args.plot_or_save == 'plot': plt.show()
                if args.plot_or_save == 'save': plt.savefig('{}/noisecheck_{}.png'.format(folder, filename), bbox_inches='tight', dpi=400)
                plt.close()
                '''


            for i in show_slits:
                pypeit_id = all_pypeit_ids[i]
                if all_maskdef_ids is not None:
                    filename = '{}_DET{}_maskdefID{}_pypeitID{}'.format(
                        spec2DObj.head0['DECKER'], args.det, all_maskdef_ids[i], pypeit_id)
                else:
                    filename = '{}_DET{}_pypeitID{}'.format(
                        spec2DObj.head0['DECKER'], args.det, pypeit_id)

                # Chi
                chi_slit, _, _ = spec2DObj.calc_chi_slit(i, pad=args.pad)

                # Cut down
                chi_select = chi_slit * input_mask
                if np.all(chi_select == 0):
                    msgs.warn(f"All of the chi values are masked in {filename}!")
                    continue

                # Flux to show
                flux_select = spec2DObj.sciimg - spec2DObj.skymodel
                if spec2DObj.objmodel is not None:
                    flux_select -= spec2DObj.objmodel
                flux_select *= input_mask
                # Error
                err_select = 1/np.sqrt(spec2DObj.ivarmodel)* input_mask

                # Wavelengths
                left, right, _ = spec2DObj.slits.select_edges()
                spat_start = int(left[:, i].min())
                spat_end = int(right[:, i].max())

                slit_select = spec2DObj.slits.slit_img(pad=args.pad, slitidx=slitidx)
                in_slit = slit_select == spec2DObj.slits.spat_id[slitidx]

                lbda = spec2DObj.waveimg[in_slit]
                if lbda[lbda!=0].size == 0:
                    msgs.warn(f"None of the wavelength values work for {filename}!")
                    continue

                line_wav_plt = np.array([])
                line_names_plt = np.array([])
                if z is not None:
                    for i in range(line_wav_z.shape[0]):
                        if (line_wav_z[i]>lbda[lbda!=0].min())&(line_wav_z[i]<lbda[lbda!=0].max()):
                            line_wav_plt = np.append(line_wav_plt, lbda.searchsorted(line_wav_z[i]))
                            line_names_plt = np.append(line_names_plt, line_names[i])

                lbda_min = args.wavemin if args.wavemin is not None else lbda[lbda!=0].min()
                lbda_max = args.wavemax if args.wavemax is not None else lbda[lbda!=0].max()
                plot(chi_slit[:, spat_start:spat_end], line_wav_plt, line_names_plt, 
                    lbda, lbda_min, lbda_max, 
                    args.aspect_ratio, chi_select, flux_select, 
                    err_select, filename)
                if args.mode == 'plot': plt.show()
                if args.mode == 'save': plt.savefig('{}/noisecheck_{}.png'.format(folder, filename), bbox_inches='tight', dpi=400)
                plt.close()
