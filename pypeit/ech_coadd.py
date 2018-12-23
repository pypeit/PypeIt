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
def load_spec_order(fname,objid=None,order=None,extract='OPT',flux=True):
    """Loading single order spectrum from a PypeIt 1D specctrum fits file.
        it will be called by ech_load_spec
    Parameters:
        fname (str) : The file name of your spec1d file
        objid (str) : The id of the object you want to load. (default is the first object)
        order (int) : which order you want to load (default is None, loading all orders)
        extract (str) : 'OPT' or 'BOX'
        flux (bool) : default is True, loading fluxed spectra
    returns:
        spectrum_out : XSpectrum1D
    """
    if objid is None:
        objid = 0
    if order is None:
        msgs.error('Please specify which order you want to load')

    # read extension name into a list
    primary_header = fits.getheader(fname, 0)
    nspec = primary_header['NSPEC']
    extnames = [primary_header['EXT0001']] * nspec
    for kk in range(nspec):
        extnames[kk] = primary_header['EXT' + '{0:04}'.format(kk + 1)]
    extnameroot = extnames[0]

    # Figure out which extension is the required data
    ordername = '{0:04}'.format(order)
    extname = extnameroot.replace('OBJ0000', objid)
    extname = extname.replace('ORDER0000', 'ORDER' + ordername)
    try:
        exten = extnames.index(extname) + 1
        msgs.info("Loading extension {:s} of spectrum {:s}".format(extname, fname))
    except:
        msgs.error("Spectrum {:s} does not contain {:s} extension".format(fname, extname))

    spectrum = load.load_1dspec(fname, exten=exten, extract=extract, flux=flux)
    # Polish a bit -- Deal with NAN, inf, and *very* large values that will exceed
    #   the floating point precision of float32 for var which is sig**2 (i.e. 1e38)
    bad_flux = np.any([np.isnan(spectrum.flux), np.isinf(spectrum.flux),
                       np.abs(spectrum.flux) > 1e30,
                       spectrum.sig ** 2 > 1e10,
                       ], axis=0)
    # Sometimes Echelle spectra have zero wavelength
    bad_wave = spectrum.wavelength < 1000.0*units.AA
    bad_all = bad_flux + bad_wave
    ## trim bad part
    wave_out,flux_out,sig_out = spectrum.wavelength[~bad_all],spectrum.flux[~bad_all],spectrum.sig[~bad_all]
    spectrum_out = XSpectrum1D.from_tuple((wave_out,flux_out,sig_out), verbose=False)
    #if np.sum(bad_flux):
    #    msgs.warn("There are some bad flux values in this spectrum.  Will zero them out and mask them (not ideal)")
    #    spectrum.data['flux'][spectrum.select][bad_flux] = 0.
    #    spectrum.data['sig'][spectrum.select][bad_flux] = 0.

    return spectrum_out

def ech_load_spec(files,objid=None,order=None,extract='OPT',flux=True):
    """Loading Echelle spectra from a list of PypeIt 1D spectrum fits files
    Parameters:
        files (str) : The list of file names of your spec1d file
        objid (str) : The id (one per fits file) of the object you want to load. (default is the first object)
        order (int) : which order you want to load (default is None, loading all orders)
        extract (str) : 'OPT' or 'BOX'
        flux (bool) : default is True, loading fluxed spectra
    returns:
        spectrum_out : XSpectrum1D
    """

    nfiles = len(files)
    if objid is None:
        objid = ['OBJ0000'] * nfiles
    elif len(objid) == 1:
        objid = objid * nfiles
    elif len(objid) != nfiles:
        msgs.error('The length of objid should be either 1 or equal to the number of spectra files.')

    fname = files[0]
    ext_final = fits.getheader(fname, -1)
    norder = ext_final['ORDER'] + 1
    msgs.info('spectrum {:s} has {:d} orders'.format(fname, norder))
    if norder <= 1:
        msgs.error('The number of orders have to be greater than one for echelle. Longslit data?')

    # Load spectra
    spectra_list = []
    for ii, fname in enumerate(files):

        if order is None:
            msgs.info('Loading all orders into a gaint spectra')
            for iord in range(norder):
                spectrum = load_spec_order(fname,objid=objid[ii],order=iord,extract=extract,flux=flux)
                # Append
                spectra_list.append(spectrum)
        elif order >= norder:
            msgs.error('order number cannot greater than the total number of orders')
        else:
            spectrum = load_spec_order(fname, objid=objid[ii], order=order, extract=extract, flux=flux)
            # Append
            spectra_list.append(spectrum)
    # Join into one XSpectrum1D object
    spectra = collate(spectra_list)
    # Return
    return spectra

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
    norder = ext_final['ORDER'] + 1
    msgs.info('spectrum {:s} has {:d} orders'.format(fname, norder))
    if norder <= 1:
        msgs.error('The number of orders have to be greater than one for echelle. Longslit data?')

    if giantcoadd:
        msgs.info('Coadding all orders and exposures at once')
        spectra = ech_load_spec(files, objid=objids,order=None, extract=extract, flux=flux)
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
            spectra = ech_load_spec(files, objid=objids, order=iord, extract=extract, flux=flux)
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
        # ToDo: Currently I'm not using the first order due to some problem. Need to add it back after fix the problem.
        spec1d = coadd.coadd_spectra(spectra_coadd[1:], wave_grid_method=wave_grid_method, niter=niter,
                                          scale_method=scale_method, do_offset=do_offset, sigrej_final=sigrej_final,
                                          do_var_corr=do_var_corr, qafile=qafile, outfile=outfile,
                                          do_cr=do_cr, debug=debug, **kwargs)

    return spec1d
