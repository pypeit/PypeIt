""" Routine for Echelle coaddition
"""
import numpy as np
from astropy import stats
from astropy.io import fits
from astropy import units
import matplotlib.pyplot as plt

from pypeit.core import coadd
from pypeit.core import load
from pypeit import msgs

from linetools.spectra.utils import collate
from linetools.spectra.xspectrum1d import XSpectrum1D
from pkg_resources import resource_filename


# setting plot parameters
plt.rcdefaults()
plt.rcParams['font.family'] = 'times new roman'
plt.rcParams["xtick.top"] = True
plt.rcParams["ytick.right"] = True
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams["ytick.direction"] = 'in'
plt.rcParams["xtick.direction"] = 'in'
plt.rcParams["xtick.labelsize"] = 17
plt.rcParams["ytick.labelsize"] = 17
plt.rcParams["axes.labelsize"] = 17

def spec_from_array(wave,flux,sig,**kwargs):
    """
    Make an XSpectrum1D from numpy arrays of wave, flux and sig
    Parameters
    ----------
        If wave is unitless, Angstroms are assumed
        If flux is unitless, it is made dimensionless
        The units for sig and co are taken from flux.
    Return spectrum from arrays of wave, flux and sigma
    """

    # Get rid of 0 wavelength
    good_wave = (wave>1.0*units.AA)
    wave,flux,sig = wave[good_wave],flux[good_wave],sig[good_wave]
    ituple = (wave, flux, sig)
    spectrum = XSpectrum1D.from_tuple(ituple, **kwargs)
    # Polish a bit -- Deal with NAN, inf, and *very* large values that will exceed
    #   the floating point precision of float32 for var which is sig**2 (i.e. 1e38)
    bad_flux = np.any([np.isnan(spectrum.flux), np.isinf(spectrum.flux),
                       np.abs(spectrum.flux) > 1e30,
                       spectrum.sig ** 2 > 1e10,
                       ], axis=0)
    if np.sum(bad_flux):
        msgs.warn("There are some bad flux values in this spectrum.  Will zero them out and mask them (not ideal)")
        spectrum.data['flux'][spectrum.select][bad_flux] = 0.
        spectrum.data['sig'][spectrum.select][bad_flux] = 0.
    return spectrum

def order_phot_scale(spectra, phot_scale_dicts, nsig=3.0, niter=5, debug=False):
    '''
    Scale coadded spectra with photometric data.
    Parameters:
      spectra: XSpectrum1D spectra (longslit) or spectra list (echelle)
      phot_scale_dicts: A dict contains photometric information of each orders (if echelle).
        An example is given below.
        phot_scale_dicts = {0: {'filter': None, 'mag': None, 'mag_type': None, 'masks': None},
                            1: {'filter': 'UKIRT-Y', 'mag': 20.33, 'mag_type': 'AB', 'masks': None},
                            2: {'filter': 'UKIRT-J', 'mag': 20.19, 'mag_type': 'AB', 'masks': None},
                            3: {'filter': 'UKIRT-H', 'mag': 20.02, 'mag_type': 'AB', 'masks': None},
                            4: {'filter': 'UKIRT-K', 'mag': 19.92, 'mag_type': 'AB', 'masks': None}}
      Show QA plot if debug=True
    Return a new scaled XSpectrum1D spectra
    '''

    from pypeit.core.flux_calib import scale_in_filter

    norder = spectra.nspec

    # scaling spectrum order by order.
    spectra_list_new = []
    for iord in range(norder):
        phot_scale_dict = phot_scale_dicts[iord]
        if (phot_scale_dict['filter'] is not None) & (phot_scale_dict['mag'] is not None):
            speci = scale_in_filter(spectra[iord], phot_scale_dict)
        else:
            #ToDo: Think a better way to do the following
            try:
                spec0 = scale_in_filter(spectra[iord-1], phot_scale_dicts[iord-1])
                speci = spectra[iord]
                med_flux = spec0.data['flux'] / speci.data['flux']
                mn_scale, med_scale, std_scale = stats.sigma_clipped_stats(med_flux, sigma=nsig, iters=niter)
                med_scale = np.minimum(med_scale, 5.0)
                spectra.data['flux'] *= med_scale
                spectra.data['sig'] *= med_scale
                msgs.warn("Not enough photometric information given. Scaled order {:d} to order {:d}".format(iord, iord-1))
            except KeyError:
                msgs.warn("Not enough photometric information given. Scale order {:d} to order {:d} failed".format(iord, iord-1))
                try:
                    spec0 = scale_in_filter(spectra[iord + 1], phot_scale_dicts[iord + 1])
                    speci = spectra[iord]
                    med_flux = spec0.data['flux'] / speci.data['flux']
                    mn_scale, med_scale, std_scale = stats.sigma_clipped_stats(med_flux, sigma=nsig, iters=niter)
                    med_scale = np.minimum(med_scale, 5.0)
                    speci.data['flux'] *= med_scale
                    speci.data['sig'] *= med_scale
                    msgs.warn("Not enough photometric information given. Scaled order {:d} to order {:d}".format(iord, iord+1))
                except:
                    msgs.warn("Not enough photometric information given. No scaling on order {:d}".format(iord))
                    speci = spectra[iord]
        spectra_list_new.append(speci)

        if debug:
            gdp = speci.sig>0
            plt.plot(spectra[iord].wavelength[gdp], spectra[iord].flux[gdp], 'k-', label='raw spectrum')
            plt.plot(speci.wavelength[gdp], speci.flux[gdp], 'b-',
                     label='scaled spectrum')
            mny, medy, stdy = stats.sigma_clipped_stats(speci.flux[gdp], sigma=3, iters=5)
            plt.ylim([0.1 * medy, 4.0 * medy])
            plt.legend()
            plt.xlabel('wavelength')
            plt.ylabel('Flux')
            plt.show()

    return collate(spectra_list_new)

def order_median_scale(spectra, nsig=3.0, niter=5, overlapfrac=0.03, num_min_pixels=50, SN_MIN_MEDSCALE=0.5, debug=False):
    '''
    Scale different orders using the median of overlap regions. It starts from the reddest order, i.e. scale H to K,
      and then scale J to H+K, etc.
    Parameters:
      spectra: XSpectrum1D spectra
      nsig: float
        sigma used for sigma_clipping median
      niter: int
        number of iterations for sigma_clipping median
      overlapfrac: float
        minmum overlap fraction (number of overlapped pixels devided by number of pixels of the whole spectrum) between orders.
      num_min_pixels: int
        minum required good pixels. The code only scale orders when the overlapped
        pixels > max(num_min_pixels,overlapfrac*len(wave))
      SN_MIN_MEDSCALE: float
        Maximum RMS S/N allowed to automatically apply median scaling
      Show QA plot if debug=True
    Return:
        No return, but the spectra is already scaled after executing this function.
    '''
    norder = spectra.nspec
    fluxes, sigs, wave = coadd.unpack_spec(spectra, all_wave=False)
    fluxes_raw = fluxes.copy()

    # scaling spectrum order by order. We use the reddest order as the reference since slit loss in redder is smaller
    for i in range(norder - 1):
        iord = norder - i - 1
        sn_iord_iref = fluxes[iord] * (1. / sigs[iord])
        sn_iord_scale = fluxes[iord - 1] * (1. / sigs[iord - 1])
        allok = (sigs[iord - 1, :] > 0) & (sigs[iord, :] > 0) & (sn_iord_iref > SN_MIN_MEDSCALE) & (
        sn_iord_scale > SN_MIN_MEDSCALE)
        if sum(allok) > np.maximum(num_min_pixels, len(wave) * overlapfrac):
            # Ratio
            med_flux = spectra.data['flux'][iord, allok] / spectra.data['flux'][iord - 1, allok]
            # Clip
            mn_scale, med_scale, std_scale = stats.sigma_clipped_stats(med_flux, sigma=nsig, iters=niter)
            med_scale = np.minimum(med_scale, 5.0)
            spectra.data['flux'][iord - 1, :] *= med_scale
            spectra.data['sig'][iord - 1, :] *= med_scale
            msgs.info('Scaled %s order by a factor of %s'%(iord,str(med_scale)))

            if debug:
                plt.plot(wave, spectra.data['flux'][iord], 'r-', label='reference spectrum')
                plt.plot(wave, fluxes_raw[iord - 1], 'k-', label='raw spectrum')
                plt.plot(spectra.data['wave'][iord - 1, :], spectra.data['flux'][iord - 1, :], 'b-',
                         label='scaled spectrum')
                mny, medy, stdy = stats.sigma_clipped_stats(fluxes[iord, allok], sigma=nsig, iters=niter)
                plt.ylim([0.1 * medy, 4.0 * medy])
                plt.xlim([np.min(wave[sigs[iord - 1, :] > 0]), np.max(wave[sigs[iord, :] > 0])])
                plt.legend()
                plt.xlabel('wavelength')
                plt.ylabel('Flux')
                plt.show()
        else:
            msgs.warn('Not enough overlap region for sticking different orders.')

def ech_coadd(files,objids=None,extract='OPT',flux=True,giantcoadd=False,orderscale='median',mergeorder=True,
              wave_grid_method='velocity', niter=5,wave_grid_min=None, wave_grid_max=None,v_pix=None,
              scale_method='auto', do_offset=False, sigrej_final=3.,do_var_corr=False,
              SN_MIN_MEDSCALE = 0.5, overlapfrac = 0.01, num_min_pixels=10,phot_scale_dicts=None,
              qafile=None, outfile=None,do_cr=True, debug=False,**kwargs):
    """
    routines for coadding spectra observed with echelle spectrograph.
    parameters:
        files (list): file names
        objids (str): objid
        extract (str): 'OPT' or 'BOX'
        flux (bool): fluxed or not
        giantcoadd (bool): coadding order by order or do it at once?
        wave_grid_method (str): default velocity
        niter (int): number of iteration for rejections
        wave_grid_min (float): min wavelength, None means it will find the min value from your spectra
        wave_grid_max (float): max wavelength, None means it will find the max value from your spectra
        v_pix (float): delta velocity, see coadd.py
        scale_method (str): see coadd.py
        do_offset (str): see coadd.py, not implemented yet.
        sigrej_final (float): see coadd.py
        do_var_corr (bool): see coadd.py, default False. It seems True will results in a large error
        SN_MIN_MEDSCALE (float): minimum SNR for scaling different orders
        overlapfrac (float): minimum overlap fraction for scaling different orders.
        qafile (str): name of qafile
        outfile (str): name of coadded spectrum
        do_cr (bool): remove cosmic rays?
        debug (bool): show debug plots?
        kwargs: see coadd.py
    returns:
        spec1d: coadded XSpectrum1D
    """

    nfile = len(files)
    if nfile <=1:
        msgs.info('Only one spectrum exits coadding...')
        return

    fname = files[0]
    ext_final = fits.getheader(fname, -1)
    norder = ext_final['ECHORDER'] + 1
    msgs.info('spectrum {:s} has {:d} orders'.format(fname, norder))
    if norder <= 1:
        msgs.error('The number of orders have to be greater than one for echelle. Longslit data?')

    if giantcoadd:
        msgs.info('Coadding all orders and exposures at once')
        spectra = load.ech_load_spec(files, objid=objids,order=None, extract=extract, flux=flux)
        wave_grid = np.zeros((2,spectra.nspec))
        for i in range(spectra.nspec):
            wave_grid[0, i] = spectra[i].wvmin.value
            wave_grid[1, i] = spectra[i].wvmax.value
        ech_kwargs = {'echelle': True, 'wave_grid_min': np.min(wave_grid), 'wave_grid_max': np.max(wave_grid),
                      'v_pix': v_pix}
        kwargs.update(ech_kwargs)
        # Coadding
        spec1d = coadd.coadd_spectra(spectra, wave_grid_method=wave_grid_method, niter=niter,
                                          scale_method=scale_method, do_offset=do_offset, sigrej_final=sigrej_final,
                                          do_var_corr=do_var_corr, qafile=qafile, outfile=outfile,
                                          do_cr=do_cr, debug=debug,**kwargs)
    else:
        msgs.info('Coadding individual orders first and then merge order')
        spectra_list = []
        # Keywords for Table
        rsp_kwargs = {}
        rsp_kwargs['wave_tag'] = '{:s}_WAVE'.format(extract)
        rsp_kwargs['flux_tag'] = '{:s}_FLAM'.format(extract)
        rsp_kwargs['sig_tag'] = '{:s}_FLAM_SIG'.format(extract)
        #wave_grid = np.zeros((2,norder))
        for iord in range(norder):
            spectra = load.ech_load_spec(files, objid=objids, order=iord, extract=extract, flux=flux)
            ech_kwargs = {'echelle': False, 'wave_grid_min': spectra.wvmin.value, 'wave_grid_max': spectra.wvmax.value, 'v_pix': v_pix}
            #wave_grid[0,iord] = spectra.wvmin.value
            #wave_grid[1,iord] = spectra.wvmax.value
            kwargs.update(ech_kwargs)
            # Coadding the individual orders
            if qafile is not None:
                qafile_iord = qafile+'_%s'%str(iord)
            else:
                qafile_iord =  None
            spec1d_iord = coadd.coadd_spectra(spectra, wave_grid_method=wave_grid_method, niter=niter,
                                       scale_method=scale_method, do_offset=do_offset, sigrej_final=sigrej_final,
                                       do_var_corr=do_var_corr, qafile=qafile_iord, outfile=None,
                                       do_cr=do_cr, debug=debug, **kwargs)
            spectrum = spec_from_array(spec1d_iord.wavelength, spec1d_iord.flux, spec1d_iord.sig,**rsp_kwargs)
            spectra_list.append(spectrum)

        spectra_coadd = collate(spectra_list)

        # Rebin the spectra
        # ToDo: we should read in JFH's wavelength grid here.
        # Join into one XSpectrum1D object
        # Final wavelength array
        kwargs['wave_grid_min'] = np.min(spectra_coadd.data['wave'][spectra_coadd.data['wave'] > 0])
        kwargs['wave_grid_max'] = np.max(spectra_coadd.data['wave'][spectra_coadd.data['wave'] > 0])
        wave_final = coadd.new_wave_grid(spectra_coadd.data['wave'], wave_method=wave_grid_method, **kwargs)
        # The rebin function in linetools can not work on collated spectra (i.e. filled 0).
        # Thus I have to rebin the spectra first and then collate again.
        spectra_list_new = []
        for i in range(spectra_coadd.nspec):
            speci = spectra_list[i].rebin(wave_final * units.AA, all=True, do_sig=True, grow_bad_sig=True,
                                          masking='none')
            spectra_list_new.append(speci)
        spectra_coadd_rebin = collate(spectra_list_new)

        ## Note
        if orderscale == 'photometry':
            # Only tested on NIRES.
            if phot_scale_dicts is not None:
                spectra_coadd_rebin = order_phot_scale(spectra_coadd_rebin, phot_scale_dicts, debug=debug)
            else:
                msgs.warn('No photometric information is provided. Will use median scale.')
                orderscale = 'median'
        elif orderscale == 'median':
            #rmask = spectra_coadd_rebin.data['sig'].filled(0.) > 0.
            #sn2, weights = coadd.sn_weights(fluxes, sigs, rmask, wave)
            ## scaling different orders
            order_median_scale(spectra_coadd_rebin, nsig=sigrej_final, niter=niter, overlapfrac=overlapfrac,
                               num_min_pixels=num_min_pixels, SN_MIN_MEDSCALE=SN_MIN_MEDSCALE, debug=debug)
        else:
            msgs.warn('No any scaling is performed between different orders.')


        if mergeorder:
            fluxes, sigs, wave = coadd.unpack_spec(spectra_coadd_rebin, all_wave=False)
            ## Megering orders
            msgs.info('Merging different orders')
            ## ToDo: Joe claimed not to use pixel depedent weighting.
            weights = 1.0 / sigs**2
            weights[~np.isfinite(weights)] = 0.0
            weight_combine = np.sum(weights, axis=0)
            weight_norm = weights / weight_combine
            weight_norm[np.isnan(weight_norm)] = 1.0
            flux_final = np.sum(fluxes * weight_norm, axis=0)
            sig_final = np.sqrt(np.sum((weight_norm * sigs) ** 2, axis=0))
            spec1d_final = spec_from_array(wave_final * units.AA,flux_final,sig_final,**rsp_kwargs)

            if outfile is not None:
                msgs.info('Saving the final calibrated spectrum as {:s}'.format(outfile))
                coadd.write_to_disk(spec1d_final, outfile)

            if (qafile is not None) or (debug):
                # plot and save qa
                plt.figure(figsize=(12, 6))
                ax1 = plt.axes([0.07, 0.13, 0.9, 0.4])
                ax2 = plt.axes([0.07, 0.55, 0.9, 0.4])
                plt.setp(ax2.get_xticklabels(), visible=False)

                medf = np.median(spec1d_final.flux)
                ylim = (np.sort([0. - 0.3 * medf, 5 * medf]))
                cmap = plt.get_cmap('RdYlBu_r')
                for idx in range(spectra_coadd_rebin.nspec):
                    spectra_coadd_rebin.select = idx
                    color = cmap(float(idx) / spectra_coadd_rebin.nspec)
                    ind_good = spectra_coadd_rebin.sig > 0
                    ax1.plot(spectra_coadd_rebin.wavelength[ind_good], spectra_coadd_rebin.flux[ind_good], color=color)

                if (np.max(spec1d_final.wavelength) > (9000.0 * units.AA)):
                    skytrans_file = resource_filename('pypeit', '/data/skisim/atm_transmission_secz1.5_1.6mm.dat')
                    skycat = np.genfromtxt(skytrans_file, dtype='float')
                    scale = 0.85 * ylim[1]
                    ax2.plot(skycat[:, 0] * 1e4, skycat[:, 1] * scale, 'm-', alpha=0.5)

                ax2.plot(spec1d_final.wavelength, spec1d_final.sig, ls='steps-', color='0.7')
                ax2.plot(spec1d_final.wavelength, spec1d_final.flux, ls='steps-', color='b')

                ax1.set_xlim([np.min(spec1d_final.wavelength.value), np.max(spec1d_final.wavelength.value)])
                ax2.set_xlim([np.min(spec1d_final.wavelength.value), np.max(spec1d_final.wavelength.value)])
                ax1.set_ylim(ylim)
                ax2.set_ylim(ylim)
                ax1.set_xlabel('Wavelength (Angstrom)')
                ax1.set_ylabel('Flux')
                ax2.set_ylabel('Flux')

                plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.2)

                if len(qafile.split('.')) == 1:
                    msgs.info("No fomat given for the qafile, save to PDF format.")
                    qafile = qafile + '.pdf'
                if qafile:
                    plt.savefig(qafile)
                    msgs.info("Wrote coadd QA: {:s}".format(qafile))
                if debug:
                    plt.show()
                plt.close()

            ### Do NOT remove this part althoug it is deprecated.
            # we may need back to using this pieces of code after fixing the coadd.coadd_spectra problem on first order.
            #kwargs['echelle'] = True
            #kwargs['wave_grid_min'] = np.min(wave_grid)
            #kwargs['wave_grid_max'] = np.max(wave_grid)
            #spec1d_final = coadd.coadd_spectra(spectra_coadd_rebin, wave_grid_method=wave_grid_method, niter=niter,
            #                                  scale_method=scale_method, do_offset=do_offset, sigrej_final=sigrej_final,
            #                                  do_var_corr=do_var_corr, qafile=qafile, outfile=outfile,
            #                                  do_cr=do_cr, debug=debug, **kwargs)
            return spec1d_final
        else:
            msgs.warn('Skipped merging orders')
            if outfile is not None:
                for iord in range(len(spectra_list)):
                    outfile_iord = outfile.replace('.fits','_ORDER{:04d}.fits'.format(iord))
                    msgs.info('Saving the final calibrated spectrum of order {:d} as {:s}'.format(iord,outfile))
                    spectra_list[iord].write_to_fits(outfile_iord)
            return spectra_list