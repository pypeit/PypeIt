""" Routine for Echelle coaddition
"""
import numpy as np
from astropy.io import fits
from astropy import units
import matplotlib.pyplot as plt

from pypeit.core import coadd
from pypeit.core import load
from pypeit import msgs

from linetools.spectra.utils import collate
from linetools.spectra.xspectrum1d import XSpectrum1D

## ToDo: change it to a CLASS and modify coadd_1dspec.py

def spec_from_array(wave,flux,sig,**kwargs):
    """
    return spectrum from arrays of wave, flux and sigma
    """

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


def ech_coadd(files,objids=None,extract='OPT',flux=True,giantcoadd=False,
              wave_grid_method='velocity', niter=5,wave_grid_min=None, wave_grid_max=None,v_pix=None,
              scale_method='auto', do_offset=False, sigrej_final=3.,do_var_corr=False,
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
        wave_grid = np.zeros((2,norder))
        for iord in range(norder):
            spectra = load.ech_load_spec(files, objid=objids, order=iord, extract=extract, flux=flux)
            ech_kwargs = {'echelle': False, 'wave_grid_min': spectra.wvmin.value, 'wave_grid_max': spectra.wvmax.value, 'v_pix': v_pix}
            wave_grid[0,iord] = spectra.wvmin.value
            wave_grid[1,iord] = spectra.wvmax.value
            kwargs.update(ech_kwargs)
            # Coadding the individual orders
            if qafile is not None:
                qafile_iord = qafile+'_%s'%str(iord)
            else:
                qafile_iord =  None
            spec1d_iord = coadd.coadd_spectra(spectra, wave_grid_method=wave_grid_method, niter=niter,
                                       scale_method=scale_method, do_offset=do_offset, sigrej_final=sigrej_final,
                                       do_var_corr=do_var_corr, qafile=qafile_iord, outfile=outfile,
                                       do_cr=do_cr, debug=debug, **kwargs)
            spectrum = spec_from_array(spec1d_iord.wavelength, spec1d_iord.flux, spec1d_iord.sig,**rsp_kwargs)
            spectra_list.append(spectrum)
        # Join into one XSpectrum1D object
        spectra_coadd = collate(spectra_list)
        kwargs['echelle'] = True
        kwargs['wave_grid_min'] = np.min(wave_grid)
        kwargs['wave_grid_max'] = np.max(wave_grid)
        spec1d = coadd.coadd_spectra(spectra_coadd, wave_grid_method=wave_grid_method, niter=niter,
                                          scale_method=scale_method, do_offset=do_offset, sigrej_final=sigrej_final,
                                          do_var_corr=do_var_corr, qafile=qafile, outfile=outfile,
                                          do_cr=do_cr, debug=debug, **kwargs)

    return spec1d
