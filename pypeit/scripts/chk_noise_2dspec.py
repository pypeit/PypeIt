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
from astropy.table import Table

from pypeit import spec2dobj
from pypeit import msgs
from pypeit import io
from pypeit import utils
from pypeit.scripts import scriptbase
from pypeit.images.detector_container import DetectorContainer
from pypeit.utils import list_of_spectral_lines

from IPython import embed


def plot(image:np.ndarray, chi_select:np.ndarray, flux_select:np.ndarray,
         err_select:np.ndarray, basename:str, line_wav, line_names,
         lbda_1d, lbda_min=None, lbda_max=None, aspect_ratio=1):
    """ Generate the plot

    Args:
        image (np.ndarray):
            Image of the slit to plot.
        chi_select (`numpy.ndarray`_):
            2D array of the chi value for a selected part of the slit.
        line_wav (`numpy.ndarray`_):
            Array of wavelength of spectral features to be plotted.
        line_names (`numpy.ndarray`_):
            Array of names of spectral features to be plotted.
        lbda_1d (np.ndarray):
            1D array of the wavelength
        lbda_min (float):
            Minimum wavelength of the select pat of the slit
        lbda_max (float):
            Maximum wavelength of the select pat of the slit
        aspect_ratio ([type]):
            Aspect ratio for plotting the spec2d image.
        flux_select (`numpy.ndarray`_):
            Flux of the 2D spectrum
        err_select (`numpy.ndarray`_):
            Sig of the 2D spectrum
        basename (:obj:`str`):
            Basename
    """
    fig = plt.figure(figsize=(23,4.))
    ax = plt.subplot2grid((1, 4), (0, 0), rowspan=1, colspan=3)
    ax.minorticks_on()
    zmax = sigma_clip(image[image!=0], sigma=2, return_bounds=True)[2]*1.3
    zmin = sigma_clip(image[image!=0], sigma=2, return_bounds=True)[1]*1.3
    ax.imshow(image.T, origin ='lower', interpolation='nearest', aspect=aspect_ratio,
              vmin=zmin, vmax=zmax, cmap=plt.get_cmap('gist_gray'))
    ax.text(0.005, -0.15, '{}'.format(basename), color='k', fontsize=13, horizontalalignment='left',
            transform=ax.transAxes, bbox=dict(edgecolor='black', facecolor='white', linewidth=1))
    if line_wav.size>0:
        line_num=0
        for i in range(line_wav.size):
            line_num += 1
            yannot = 1.03 if (i % 2) == 0 else 1.13
            ax.axvline(line_wav[i], color='black', ls='dotted', zorder=2)
            ax.annotate('{}'.format(line_names[i]), xy=(line_wav[i],1), xytext=(line_wav[i], yannot),
                        xycoords=('data', 'axes fraction'),arrowprops=dict(facecolor='None', edgecolor='None',
                                                                           headwidth=0., headlength=0, width=0,
                                                                           shrink=0.),
                        annotation_clip=True, horizontalalignment='center', color='k', fontsize=10)

    if lbda_min is not None or lbda_max is not None:
        if lbda_min is None:
            lbda_min = lbda_1d.min()
        if lbda_max is None:
            lbda_max = lbda_1d.max()
        ax.axvspan(lbda_1d.searchsorted(lbda_min), lbda_1d.searchsorted(lbda_max), color='tab:green', alpha=0.2, zorder=2)
    ax.set_xticks([])
    ax.set_yticks([])

    # Guassian stats
    ax2=plt.subplot2grid((1, 4), (0, 3), rowspan=1, colspan=1)
    ax2.minorticks_on()

    bins=np.arange(chi_select[chi_select!=0].min(), chi_select[chi_select!=0].max(), 0.1)
    hist_n, hist_bins, _ = ax2.hist(chi_select[chi_select!=0], bins=bins, histtype='stepfilled')
    mod_mods=Gaussian1D(amplitude=hist_n.max(), mean=np.median(chi_select[chi_select!=0]), stddev=1.)
    ax2.plot(bins, mod_mods(bins), label=r"$\sigma=1$")
    ax2.axvline(0, ls='dotted', color='Gray')
    ax2.set_xlim(-6.5,6.5)
    ax2.set_ylim(-0.02, hist_n.max()*1.5)
    ax2.set_xlabel(r'(sciimg - skymodel - objmodel) * sqrt(ivarmodel) * (bpmmask == 0)')
    ax2.set_ylabel(r'#')
    err_over_flux = (np.median(err_select[flux_select!=0])/mad_std(flux_select[flux_select!=0]))
    ax2.text(0.97, 0.93, r'Median Noise= {:.1f} - Flux RMS= {:.1f} --> {:.2f}x'.format(
        np.median(err_select[flux_select!=0]), mad_std(flux_select[flux_select!=0]), err_over_flux),
             color='k', fontsize=9, horizontalalignment='right', transform=ax2.transAxes)
    ax2.text(0.97, 0.87, r'Chi:  Median = {:.2f}, Std = {:.2f}'.format(
        np.median(chi_select[chi_select!=0]), mad_std(chi_select[chi_select!=0])),
             color='k', fontsize=12, horizontalalignment='right', transform=ax2.transAxes, weight='bold')
    ax2.legend(loc=2)
    plt.tight_layout()


def get_flux_slit(spec2DObj, slitidx, pad=0):
    """
    Returns the flux and error of a specific slit.
    The flux would be sky subtracted and object removed.

    Args:
        spec2DObj (:class:`~pypeit.spec2dobj.Spec2DObj`): 2D spectra object
        slitidx (int): Given slit/order
        pad (int, optional):  Ignore pixels within pad of edges.

    Returns:
        :obj:`tuple`: tuple of `numpy.ndarray`_ with flux and error of the 2D spectrum

    """
    slit_select = spec2DObj.slits.slit_img(pad=pad, slitidx=slitidx)

    flux = spec2DObj.sciimg - spec2DObj.skymodel
    if spec2DObj.objmodel is not None:
        flux -= spec2DObj.objmodel
    flux_slit = flux * (slit_select == spec2DObj.slits.spat_id[slitidx])
    # Error
    err_slit = np.sqrt(utils.inverse(spec2DObj.ivarmodel)) * (slit_select == spec2DObj.slits.spat_id[slitidx])
    return flux_slit, err_slit


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
        parser.add_argument('--maskdef_id', default=None, type=int, help='MASKDEF_ID of the slit that '
                                                                         'you want to plot. If maskdef_id is '
                                                                         'not provided, nor a pypeit_id, all the'
                                                                         ' 2D spectra in the file(s) will be plotted.')
        parser.add_argument('--pypeit_id', default=None, type=int, help='PypeIt ID of the slit that '
                                                                        'you want to plot. If pypeit_id is not '
                                                                        'provided, nor a maskdef_id, all the '
                                                                        '2D spectra in the file(s) will be plotted.')
        parser.add_argument('--pad', default=-5, type=int, help='Padding for the selected slit. '
                                                                'Negative value will trim.')
        parser.add_argument('--aspect_ratio', default=3, type=int, help='Aspect ratio when plotting the spec2d')
        parser.add_argument('--wavemin', default=None, type=float, help='Wavelength min. This is for selecting a '
                                                                        'region of the spectrum to analyze.')
        parser.add_argument('--wavemax', default=None, type=float, help='Wavelength max.This is for selecting a '
                                                                        'region of the spectrum to analyze.')
        parser.add_argument('--mode', default='plot', type=str, help='Options are: plot, save, print'
                                                                     'Do you want to save to disk or open a plot '
                                                                     'in a mpl window. If you choose save, a '
                                                                     'folder called spec2d*_noisecheck will be '
                                                                     'created and all the relevant plot will be '
                                                                     'placed there. If you choose print, check noise '
                                                                     'value are printed in the terminal')
        parser.add_argument('--list', default=False, help='List the extensions only?',
                            action='store_true')
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
        line_names, line_wav = list_of_spectral_lines()
            
        files=np.array(args.files)

        if args.z is not None:
            zs = np.array(args.z) 

        # Loop on the files
        for i in range(files.size):    
            # reinitialize lines wave
            line_wav_z = line_wav.copy()

            # Load 2D object
            file = files[i]
            # List only?
            if args.list:
                io.fits_open(file).info()
                continue
            spec2DObj = spec2dobj.Spec2DObj.from_file(file, detname, chk_version=False)

            # Deal with redshifts
            if args.z is not None:
                z = zs[i] if zs.size == files.size else zs[0]
                line_wav_z *= (1+z)   #redshift linelist
            else:
                z = None

            # Save?
            folder = None
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

            # Find the slit of interest
            all_maskdef_ids = spec2DObj.slits.maskdef_id
            all_pypeit_ids = spec2DObj.slits.slitord_id
            if args.maskdef_id is not None and all_maskdef_ids is None:
                msgs.error('This spec2d does not have maskdef_id. Choose a pypeit_id insteed.')

            # Build the mask
            input_mask = spec2DObj.bpmmask.mask == 0
            if args.wavemin is not None:
                input_mask *= spec2DObj.waveimg > args.wavemin
            if args.wavemax is not None:
                input_mask *= spec2DObj.waveimg < args.wavemax

            # Decide on slits to show
            show_slits = range(all_pypeit_ids.size)
            if args.pypeit_id is not None or args.maskdef_id is not None:
                if args.maskdef_id is not None and args.maskdef_id in all_maskdef_ids:
                    slitidx = np.where(all_maskdef_ids == args.maskdef_id)[0][0]
                elif args.pypeit_id is not None and args.pypeit_id in all_pypeit_ids:
                    slitidx = np.where(all_pypeit_ids == args.pypeit_id)[0][0]
                show_slits = range(slitidx, slitidx+1)

            # loop on em
            for i in show_slits:
                pypeit_id = all_pypeit_ids[i]
                if all_maskdef_ids is not None:
                    basename = '{}_{}_maskdefID{}_pypeitID{}'.format(spec2DObj.head0['DECKER'],
                                                                        detname, all_maskdef_ids[i], pypeit_id)
                else:
                    basename = '{}_{}_pypeitID{}'.format(spec2DObj.head0['DECKER'], detname, pypeit_id)

                # Chi
                chi_slit, _, _ = spec2DObj.calc_chi_slit(i, pad=args.pad)

                if chi_slit is None:
                    continue

                # Cut down
                chi_select = chi_slit * input_mask
                if np.all(chi_select == 0):
                    msgs.warn(f"All of the chi values are masked in slit {pypeit_id} of {basename}!")
                    continue

                # Flux to show
                # get flux and err from in this slit
                flux_slit, err_slit = get_flux_slit(spec2DObj, i, pad=args.pad)
                # flux in the wavelength range
                flux_select = flux_slit * input_mask
                # Error in the wavelength range
                err_select = err_slit * input_mask

                # get edges of the slit to plot
                left, right, _ = spec2DObj.slits.select_edges()
                spat_start = int(left[:, i].min())
                spat_end = int(right[:, i].max())
                mid_spat = int((spat_end + spat_start)/2.)

                # Wavelengths
                if spec2DObj.waveimg[input_mask].size == 0:
                    msgs.warn(f"None of the wavelength values work in slit {pypeit_id} of {basename}!")
                    continue
                lbda_1darray = spec2DObj.waveimg[:, mid_spat]

                line_wav_plt = np.array([])
                line_names_plt = np.array([])
                if z is not None:
                    for i in range(line_wav_z.shape[0]):
                        if lbda_1darray[lbda_1darray != 0].min() < line_wav_z[i] < lbda_1darray[lbda_1darray != 0].max():
                            line_wav_plt = np.append(line_wav_plt, lbda_1darray.searchsorted(line_wav_z[i]))
                            line_names_plt = np.append(line_names_plt, line_names[i])

                plot(chi_slit[:, spat_start:spat_end], chi_select, flux_select, err_select, basename,
                     line_wav_plt, line_names_plt, lbda_1darray, lbda_min=args.wavemin,
                     lbda_max=args.wavemax, aspect_ratio=args.aspect_ratio)
                if args.mode == 'plot':
                    plt.show()
                if args.mode == 'save':
                    plt.savefig('{}/noisecheck_{}.png'.format(folder, basename), bbox_inches='tight', dpi=400)
                plt.close()
