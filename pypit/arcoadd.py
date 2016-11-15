""" Class for coaddition
"""
from __future__ import absolute_import, division, print_function

import numpy as np
import astropy.stats
import scipy.interpolate
from  scipy.signal import medfilt

from astropy import units as u
from astropy.io import fits
from linetools.spectra import xspectrum1d

from pypit import arload
from pypit import armsgs
from pypit import arqa

try:
    from xastropy.xutils import xdebug as debugger
except ImportError:
    import pdb as debugger

# Logging
msgs = armsgs.get_logger()

# TODO
    # Shift spectra
    # Scale by poly
    # Better rejection
    # Grow mask in final_rej?


def new_wave_grid(waves, method='iref', iref=0, A_pix=None, pix_size=None, **kwargs):
    """ Create a new wavelength grid for the
    spectra to be rebinned and coadded on

    Parameters
    ----------
    waves : masked ndarray
        Set of N original wavelength arrays
        Nspec, Npix
    method : str, optional
        Desired method for creating new wavelength grid.
        'iref' -- Use the first wavelength array (default)
        'velocity' -- Constant velocity
        'pixel' -- Constant pixel grid
        'concatenate' -- Meld the input wavelength arrays
    iref : int, optional
      Reference spectrum
    A_pix : float
      Pixel size in same units as input wavelength array (e.g. Angstroms)

    Returns
    -------
    wave_grid : ndarray
        New wavelength grid, not masked
    """
    # Eventually add/change this to also take in slf, which has
    # slf._argflag['reduce']['pixelsize'] = 2.5?? This won't work
    # if running coadding outside of PYPIT, which we'd like as an
    # option!
    from numpy.ma.core import MaskedArray
    if not isinstance(waves, MaskedArray):
        waves = np.ma.array(waves)

    if method == 'velocity': # Constant km/s
        # Loop over spectra and save wavelength arrays to find min, max of
        # wavelength grid
        wave_grid_min = np.min(waves)
        wave_grid_max = np.max(waves)

        wave_grid = [wave_grid_min]
        count = 0

        while max(wave_grid) <= wave_grid_max:
            # How do we determine a reasonable constant velocity? (the 100. here is arbitrary)
            step = wave_grid[count] * (100. / 299792.458)
            wave_grid.append(wave_grid[count] + step)
            count += 1

        wave_grid = np.asarray(wave_grid)

    elif method == 'pixel': # Constant Angstrom
        if pix_size is None:
            msgs.error("Need to provide pixel size with this method")
        #
        wave_grid_min = np.min(waves)
        wave_grid_max = np.max(waves)

        wave_grid = np.arange(wave_grid_min, wave_grid_max + A_pix, A_pix)

    elif method == 'concatenate':  # Concatenate
        # Setup
        loglam = np.log10(waves)
        nspec = waves.shape[0]
        newloglam = loglam[iref, :].compressed()  # Deals with mask
        # Loop
        for j in range(nspec):
            if j == iref:
                continue
            #
            iloglam = loglam[j,:].compressed()
            dloglam_0 = (newloglam[1]-newloglam[0])
            dloglam_n =  (newloglam[-1] - newloglam[-2]) # Assumes sorted
            if (newloglam[0] - iloglam[0]) > dloglam_0:
                kmin = np.argmin(np.abs(iloglam - newloglam[0] - dloglam_0))
                newloglam = np.concatenate([iloglam[:kmin], newloglam])
            #
            if (iloglam[-1] - newloglam[-1]) > dloglam_n:
                kmin = np.argmin(np.abs(iloglam - newloglam[-1] - dloglam_n))
                newloglam = np.concatenate([newloglam, iloglam[kmin:]])
        # Finish
        wave_grid = 10**newloglam
    elif method == 'iref':  # Concatenate
        wave_grid = waves[iref, :].compressed()
    else:
        msgs.error("Bad method for scaling: {:s}".format(method))
    # Concatenate of any wavelengths in other indices that may extend beyond that of wavelengths[0]?

    return wave_grid


def gauss1(x, parameters):
    """ Simple Gaussian
    Parameters
    ----------
    x : ndarray
    parameters : ??

    Returns
    -------

    """
    sz = x.shape[0]
    
    if sz+1 == 5:
        smax = float(26)
    else:
        smax = 13.
    
    if len(parameters) >= 3:
        norm = parameters[2]
    else:
        norm = 1.
        
    u = ( (x - parameters[0]) / max(np.abs(parameters[1]), 1e-20) )**2.
    
    x_mask = np.where(u < smax**2)[0]
    norm = norm / (np.sqrt(2. * np.pi)*parameters[1])
                   
    return norm * x_mask * np.exp(-0.5 * u * x_mask)


def sn_weight(spectra, debug=False):
    """ Calculate the S/N of each input spectrum and
    create an array of weights by which to weight the
    spectra by in coadding.

    Parameters
    ----------
    spectra : XSpectrum1D
        New wavelength grid

    Returns
    -------
    sn2 : array
        Mean S/N^2 value for each input spectra
    weights : ndarray
        Weights to be applied to the spectra
    """
    # Setup
    fluxes = spectra.data['flux']
    sigs = spectra.data['sig']
    wave = spectra.data['wave'][0,:]

    # Calculate
    sn_val = fluxes * (1./sigs)  # Taking flux**2 biases negative values
    sn_sigclip = astropy.stats.sigma_clip(sn_val, sigma=3, iters=5)
    sn = np.mean(sn_sigclip, axis=1).compressed()
    sn2 = sn**2 #S/N^2 value for each spectrum

    rms_sn = np.sqrt(np.mean(sn2)) # Root Mean S/N**2 value for all spectra

    if rms_sn <= 4.0:
        msgs.info("Using constant weights for coadding, RMS S/N = {:g}".format(rms_sn))
        weights = np.ma.outer(np.asarray(sn2), np.ones(fluxes.shape[1]))
    else:
        msgs.info("Using wavelength dependent weights for coadding")
        msgs.warn("If your spectra have very different dispersion, this is *not* accurate")
        sn_med1 = np.ones_like(fluxes) #((fluxes.shape[0], fluxes.shape[1]))
        weights = np.ones_like(fluxes) #((fluxes.shape[0], fluxes.shape[1]))

        bkspace = (10000.0/3.0e5) / (np.log(10.0))
        med_width = wave.shape[0] / ((np.max(wave) - np.min(wave)) / bkspace)
        sig_res = max(med_width, 3)
        nhalf = int(sig_res) * 4L
        xkern = np.arange(0, 2*nhalf+2, dtype='float64')-nhalf

        for spec in xrange(fluxes.shape[0]):
            sn_med1[spec] = medfilt(sn_val[spec], kernel_size = 3)
        
        yvals = gauss1(xkern, [0.0, sig_res, 1, 0])

        for spec in xrange(fluxes.shape[0]):
            weights[spec] = scipy.ndimage.filters.convolve(sn_med1[spec], yvals)**2

    # Give weights the same mask (important later)
    weights.mask = fluxes.mask

    # Finish
    return sn2, weights


def grow_mask(initial_mask, n_grow=1):
    """ Grows sigma-clipped mask by n_grow pixels on each side

    Parameters
    ----------
    initial_mask : ndarray
        Initial mask for the flux + variance arrays
    n_grow : int, optional
        Number of pixels to grow the initial mask by
        on each side. Defaults to 1 pixel

    Returns
    -------
    grow_mask : ndarray
        Final mask for the flux + variance arrays
    """
    if not isinstance(n_grow, int):
        msgs.error("n_grow must be an integer")
    #
    bad_pix_spec = np.where(initial_mask == True)[0]
    bad_pix_loc = np.where(initial_mask == True)[1]
    
    grow_mask = np.ma.copy(initial_mask)
    npix = grow_mask.shape[1]
    
    if len(bad_pix_spec) > 0:
        for i in range(0, len(bad_pix_spec)):
            if initial_mask[bad_pix_spec[i]][bad_pix_loc[i]]:
                msk_p = bad_pix_loc[i] + np.arange(-1*n_grow, n_grow+1)
                gdp = (msk_p >= 0) & (msk_p < npix)
                grow_mask[bad_pix_spec[i]][msk_p[gdp]] = True
    # Return
    return grow_mask


def median_flux(spec, mask=None, nsig=3., niter=5, **kwargs):
    """ Calculate the characteristic, median flux of a spectrum

    Parameters
    ----------
    spec : XSpectrum1D
    mask : ndarray, optional
      Additional input mask with True = masked
      This needs to have the same size as the masked spectrum
    nsig : float, optional
      Clip sigma
    niter : int, optional
      Number of clipping iterations
    **kwargs : optional
      Passed to each call of sigma_clipped_stats

    Returns
    -------
    med_spec, std_spec
    """
    from astropy.stats import sigma_clipped_stats
    #goodpix = WHERE(refivar GT 0.0 AND finite(refflux) AND finite(refivar) $
    #            AND refmask EQ 1 AND refivar LT 1.0d8)
    mean_spec, med_spec, std_spec = sigma_clipped_stats(spec.flux, sigma=nsig, iters=niter, **kwargs)
    # Clip a bit
    #badpix = np.any([spec.flux.value < 0.5*np.abs(med_spec)])
    badpix = spec.flux.value < 0.5*np.abs(med_spec)
    mean_spec, med_spec, std_spec = sigma_clipped_stats(spec.flux.value, mask=badpix,
                                                        sigma=nsig, iters=niter, **kwargs)
    # Return
    return med_spec, std_spec


def scale_spectra(spectra, sn2, iref=0, method='auto', hand_scale=None,
                  SN_MAX_MEDSCALE=2., SN_MIN_MEDSCALE=0.5):
    """
    Parameters
    ----------
    spectra : XSpecrum1D
      Rebinned spectra
      These should be registered, i.e. pixel 0 has the same wavelength for all
    sn2 : ndarray
      S/N**2 estimates for each spectrum
    iref : int, optional
      Index of reference spectrum
    method : str, optional
      Method for scaling
       'auto' -- Use automatic method based on RMS of S/N
       'hand' -- Use input scale factors
       'median' -- Use calcualted median value
    SN_MIN_MEDSCALE : float, optional
      Maximum RMS S/N allowed to automatically apply median scaling
    SN_MAX_MEDSCALE : float, optional
      Maximum RMS S/N allowed to automatically apply median scaling

    Returns
    -------
    scales : list of float or ndarray
      Scale value (or arrays) applied to the data
    omethod : str
      Method applied (mainly useful if auto was adopted)
       'hand'
       'median'
       'none_SN'
    """
    # Init
    med_ref = None
    rms_sn = np.sqrt(np.mean(sn2)) # Root Mean S/N**2 value for all spectra
    # Check for wavelength registration
    gdp = np.all(~spectra.data['flux'].mask, axis=0)
    gidx = np.where(gdp)[0]
    if not np.isclose(spectra.data['wave'][0,gidx[0]],spectra.data['wave'][1,gidx[0]]):
        msgs.error("Input spectra are not registered!")
    # Loop on exposures
    scales = []
    for qq in xrange(spectra.nspec):
        if method == 'hand':
            omethod = 'hand'
            # Input?
            if hand_scale is None:
                msgs.error("Need to provide hand_scale parameter, one value per spectrum")
            spectra.data['flux'][qq,:] *= hand_scale[qq]
            spectra.data['sig'][qq,:] /= hand_scale[qq]
            #arrsky[*, j] = HAND_SCALE[j]*sclsky[*, j]
            scales.append(hand_scale[qq])
            #
        elif ((rms_sn <= SN_MAX_MEDSCALE) and (rms_sn > SN_MIN_MEDSCALE)) or method=='median':
            omethod = 'median'
            # Reference
            if med_ref is None:
                spectra.select = iref
                med_ref, std_ref = median_flux(spectra)
            # Calc
            spectra.select = qq
            med_spec, std_spec = median_flux(spectra)
            # Apply
            med_scale= np.minimum(med_ref/med_spec, 10.0)
            spectra.data['flux'][qq,:] *= med_scale
            spectra.data['sig'][qq,:] /= med_scale
            #
            scales.append(med_scale)
        elif rms_sn <= SN_MIN_MEDSCALE:
            omethod = 'none_SN'
        elif (rms_sn > SN_MAX_MEDSCALE) or method=='poly':
            msgs.work("Should be using poly here, not median")
            omethod = 'median'
            # Reference
            if med_ref is None:
                spectra.select = iref
                med_ref, std_ref = median_flux(spectra)
            # Calc
            spectra.select = qq
            med_spec, std_spec = median_flux(spectra)
            # Apply
            med_scale= np.minimum(med_ref/med_spec, 10.0)
            spectra.data['flux'][qq,:] *= med_scale
            spectra.data['sig'][qq,:] /= med_scale
            #
            scales.append(med_scale)
        else:
            msgs.error('uh oh')
    # Finish
    return scales, omethod


def clean_cr(spectra, n_grow_mask=1, nsig=5.):
    """ Sigma-clips the flux arrays to remove obvious CR

    Parameters
    ----------
    spectra :
    n_grow_mask : int, optional
        Number of pixels to grow the initial mask by
        on each side. Defaults to 1 pixel
    nsig : float, optional
      Number of sigma for rejection

    Returns
    -------
    final_mask : ndarray
        Final mask for the flux + variance arrays
    """
    # This mask may include masked pixels (including padded ones)
    #   We should *not* grow those
    first_mask = spectra.data['flux'].mask.copy()
    # New mask
    new_mask = first_mask.copy()
    new_mask[:] = False

    # Median of the masked arrays
    refflux = np.ma.median(spectra.data['flux'],axis=0)
    diff = spectra.data['flux']-refflux

    # Loop on spectra
    for ispec in xrange(spectra.nspec):
        spectra.select = ispec
        ivar = spectra.ivar
        chi2 = (diff[ispec].compressed())**2 * ivar
        badchi = (ivar > 0.0) & (chi2 > nsig**2)
        nbad = np.sum(badchi)
        if nbad > 0:
            spectra.add_to_mask(badchi, compressed=True)
            msgs.info("Rejecting {:d} CRs in exposure {:d}".format(nbad,ispec))

    # Grow new mask
    if n_grow_mask > 0:
        new_mask = grow_mask(new_mask, n_grow=n_grow_mask)

    # Final mask
    final_mask = first_mask & new_mask

    return final_mask


def one_d_coadd(spectra, weights):
    """ Performs a weighted coadd of the spectra in 1D.

    Parameters
    ----------
    spectra : XSpectrum1D
    weights : ndarray
      Should be masked

    Returns
    -------
    coadd : XSpectrum1D

    """
    from linetools.spectra.xspectrum1d import XSpectrum1D
    # Sum weights
    sum_weights = np.ma.sum(weights, axis=0)

    # Setup
    fluxes = spectra.data['flux']
    variances = spectra.data['sig']**2
    inv_variances = 1./variances

    # Coadd
    new_flux = np.ma.sum(weights*fluxes, axis=0) / (sum_weights + (sum_weights == 0.0).astype(int))
    var = (variances != 0.0).astype(float) / (inv_variances + (inv_variances == 0.0).astype(float))
    new_var = np.ma.sum((weights**2.)*var, axis=0) / ((sum_weights + (sum_weights == 0.0).astype(int))**2.)

    # Replace masked values with zeros
    new_flux = new_flux.filled(0.)
    new_sig = np.sqrt(new_var.filled(0.))

    # New obj
    wv = np.array(spectra.data['wave'][0,:])
    new_spec = XSpectrum1D.from_tuple((wv, new_flux, new_sig), masking='none')

    return new_spec

'''
# Move this to arqa
def qa_plots(wavelengths, fluxes, variances, new_wave, new_flux, new_var):
    """  QA plot
    Parameters
    ----------
    wavelengths
    fluxes
    variances
    new_wave
    new_flux
    new_var

    Returns
    -------

    """
    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    
    qa_plots = PdfPages(target + instru + '.pdf')
    
    plt.figure()

    dev_sig = (np.ma.getdata(fluxes) - new_flux) / (np.sqrt(np.ma.getdata(variances) + new_var))
    dev_sig_clip = astropy.stats.sigma_clip(dev_sig, sigma=4, iters=2)
    std_dev_devsig = np.std(dev_sig_clip)
    flat_dev_sig = dev_sig_clip.flatten()

    xmin = -10
    xmax = 10
    n_bins = 100

    hist, edges = np.histogram(flat_dev_sig, range=(xmin, xmax), bins=n_bins)
    area = len(flat_dev_sig)*((xmax-xmin)/float(n_bins))
    xppf = np.linspace(scipy.stats.norm.ppf(0.0001), scipy.stats.norm.ppf(0.9999), 100)
    plt.plot(xppf, area*scipy.stats.norm.pdf(xppf), color='black', linewidth=2.0)
    plt.gca().bar(edges[:-1], hist, width=((xmax-xmin)/float(n_bins)), alpha=0.5)
    plt.title(std_dev_devsig)
    plt.savefig(qa_plots, format='pdf')


    plt.figure()
    plt.subplots(figsize=(15,8))

    for spec in range(len(wavelengths)):
        if spec == 0:
            line_color = 'blue'
        elif spec == 1:
            line_color = 'green'
        elif spec == 2:
            line_color = 'red'

        plt.plot(wavelengths[spec], fluxes[spec], color=line_color, alpha=0.5, label='individual exposure')

    plt.plot(new_wave, new_flux, color='black', label='coadded spectrum')
    plt.legend()
    plt.title('Coadded + Original Spectra')
    plt.savefig(qa_plots, format='pdf')
    

    qa_plots.close()
    
    return
'''

'''
# What follows should be in XSpectrum1D or pypit.arsave
def save_coadd(new_wave, new_flux, new_var, outfil):
    """ Saves the coadded spectrum as a .fits file

    Parameters
    ----------
    new_wave : array
        New wavelength grid
    new_flux : array
        Coadded flux array
    new_var : array
        Variance of coadded spectrum
    outfil : str, optional
        Name of the coadded spectrum
    """
    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)
    
    col1 = fits.Column(array=new_wave, name='box_wave', format='f8')
    col2 = fits.Column(array=new_flux, name='box_flam_coadd', format='f8')
    col3 = fits.Column(array=new_var, name='box_flam_coadd_var', format='f8')
    cols = fits.ColDefs([col1, col2, col3])
    tbhdu = fits.BinTableHDU.from_columns(cols)

    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(outfil + '.fits', clobber=True)
    
    return
'''


def load_spec(files, iextensions=None, extract='opt'):
    """ Load a list of spectra into one XSpectrum1D object

    Parameters
    ----------
    files : list
      List of filenames
    iextensions : int or list, optional
      List of extensions, 1 per filename
      or an int which is the extension in each file
    extract : str, optional
      Extraction method ('opt', 'box')

    Returns
    -------
    spectra : XSpectrum1D
      -- All spectra are collated in this one object
    """
    from linetools.spectra.utils import collate
    # Extensions
    if iextensions is None:
        msgs.warn("Extensions not provided.  Assuming first extension for all")
        extensions = np.ones(len(files), dtype='int8')
    elif isinstance(iextensions, int):
        extensions = np.ones(len(files), dtype='int8') * iextensions
    else:
        extensions = np.array(iextensions)

    # Load spectra
    spectra_list = []
    for ii,fname in enumerate(files):
        msgs.info("Loading extension {:d} of spectrum {:s}".format(extensions[ii], fname))
        spectrum = arload.load_1dspec(fname, exten=extensions[ii], extract=extract)
        spectra_list.append(spectrum)
    # Join into one XSpectrum1D object
    spectra = collate(spectra_list)
    # Return
    return spectra


def get_std_dev(irspec, ispec1d):
    # Only calculate on regions with 2 or more spectra
    msk = ~irspec.data['flux'].mask
    sum_msk = np.sum(msk, axis=0)
    gdp = sum_msk > 1
    # Here we go [note that dev_sig is still a masked array so we compress it after]
    dev_sig = (irspec.data['flux'][:,gdp] - ispec1d.flux[gdp]) / (irspec.data['sig'][:,gdp]**2 + ispec1d.sig[gdp]**2)
    std_dev = np.std(astropy.stats.sigma_clip(dev_sig.compressed(), sigma=5, iters=2))
    return std_dev, dev_sig.compressed()


def coadd_spectra(spectra, wave_grid_method='concatenate', niter=5,
                  scale_method='auto', sig_clip=False,
                  do_offset=False, sigrej_final=3.,
                  do_var_corr=True, qafile=None, outfile=None,
                  do_cr=True, **kwargs):
    """
    Parameters
    ----------
    spectra : XSpectrum1D
    wave_grid_method :
    sig_clip

    Returns
    -------
    spec1d : XSpectrum1D

    """
    # Init
    if niter <= 0:
        msgs.error('Not prepared for no iterations')
    # Final wavelength array
    new_wave = new_wave_grid(spectra.data['wave'], method=wave_grid_method, **kwargs)

    # Rebin
    rspec = spectra.rebin(new_wave*u.AA, all=True, do_sig=True, masking='none')
    pre_mask = rspec.data['flux'].mask.copy()

    # Clean bad CR
    # -- Do this before sn2 and weights
    # -- Or be sure to reset the weights mask accordingly
    if do_cr:
        clean_cr(rspec)

    # S/N**2, weights
    sn2, weights = sn_weight(rspec)

    # Scale (modifies rspec)
    scales, omethod = scale_spectra(rspec, sn2, method=scale_method, **kwargs)

    # Initial coadd
    spec1d = one_d_coadd(rspec, weights)

    std_dev, _ = get_std_dev(rspec, spec1d)
    msgs.info("Initial std_dev = {:g}".format(std_dev))

    iters = 0
    std_dev = 0.
    var_corr = 1.

    while np.absolute(std_dev - 1.) >= 0.1 and iters < niter:
        iters += 1
        msgs.info("Iterating on coadding... iter={:d}".format(iters))

        # Setup (strip out masks)
        tspec = spec1d.copy()
        tspec.unmask()
        newvar = tspec.data['sig'][0,:].compressed()**2  # JFH Interpolates over bad values?
        newflux = tspec.data['flux'][0,:].compressed()
        newflux_now = newflux  # JFH interpolates
        # Convenient for coadding
        uspec = rspec.copy()
        uspec.unmask()

        # Loop on images to updated noise model for rejection
        for qq in xrange(rspec.nspec):

            # Grab full spectrum (unmasked)
            flux = uspec.data['flux'][qq,:].compressed()
            sig = uspec.data['sig'][qq,:].compressed()
            ivar = np.zeros_like(sig)
            mask = rspec.data['flux'].mask[qq,:]
            gd = sig > 0.
            ivar[gd] = 1./sig[gd]**2

            # var_tot
            var_tot = newvar + (ivar > 0.0)/(ivar + (ivar == 0.0))
            ivar_real = (var_tot > 0.0)/(var_tot + (var_tot == 0.0))
            # smooth out possible outliers in noise
            var_med = medfilt(var_tot, 5)
            var_smooth = medfilt(var_tot, 99)#, boundary = 'reflect')
            # conservatively always take the largest variance
            var_final = np.maximum(var_med, var_smooth)
            ivar_final =  (var_final > 0.0)/ (var_final + (var_final == 0.0))
            # Cap S/N ratio at SN_MAX to prevent overly aggressive rejection
            SN_MAX = 20.0
            ivar_cap = np.minimum(ivar_final,
                                  (SN_MAX/newflux_now + (newflux_now <= 0.0))**2)
            #; adjust rejection to reflect the statistics of the distribtuion
            #; of errors. This fixes cases where for not totally understood
            #; reasons the noise model is not quite right and
            #; many pixels are rejected.

            #; Is the model offset relative to the data? If so take it out
            if do_offset:
                diff1 = flux-newflux_now
                #idum = np.where(arrmask[*, j] EQ 0, nnotmask)
                nnotmask = np.sum(~mask)
                nmed_diff = np.maximum(nnotmask//20, 10)
                #; take out the smoothly varying piece
                #; JXP -- This isnt going to work well if the data has a bunch of
                #; null values in it
                w = np.ones(5, 'd')
                diff_sm = np.convolve(w/w.sum(),
                                      medfilt(diff1*(~mask), nmed_diff), mode='same')
                chi2 = (diff1-diff_sm)**2*ivar_real
                debugger.set_trace()
                goodchi = (~mask) & (ivar_real > 0.0) & (chi2 <= 36.0) # AND masklam, ngd)
                if np.sum(goodchi) == 0:
                    goodchi = np.array([True]*flux.size)
                debugger.set_trace()  # Port next line to Python to use this
                #djs_iterstat, (arrflux[goodchi, j]-newflux_now[goodchi]) $
                #   , invvar = ivar_real[goodchi], mean = offset_mean $
                #   , median = offset $
            else:
                offset = 0.
            chi2 = (flux-newflux_now - offset)**2*ivar_real
            goodchi = (~mask) & (ivar_real > 0.0) & (chi2 <= 36.0) # AND masklam, ngd)
            ngd = np.sum(goodchi)
            if ngd == 0:
                goodchi = np.array([True]*flux.size)
            #; evalute statistics of chi2 for good pixels and excluding
            #; extreme 6-sigma outliers
            chi2_good = chi2[goodchi]
            chi2_srt = chi2_good.copy()
            chi2_srt.sort()
            #; evaluate at 1-sigma and then scale
            gauss_prob = 1.0 - 2.0*(1.-scipy.stats.norm.cdf(1.)) #gaussint(-double(1.0d))
            sigind = int(np.round(gauss_prob*ngd))
            chi2_sigrej = chi2_srt[sigind]
            one_sigma = np.minimum(np.maximum(np.sqrt(chi2_sigrej),1.0),5.0)
            sigrej_eff = sigrej_final*one_sigma
            chi2_cap = (flux-newflux_now - offset)**2*ivar_cap
            # Grow??
            chi_mask = (chi2_cap > sigrej_eff**2) & (~mask)
            nrej = np.sum(chi_mask)
            # Apply
            if nrej > 0:
                msgs.info("Rejecting {:d} pixels in exposure {:d}".format(nrej,qq))
                debugger.set_trace()
                rspec.add_to_mask(chi_mask)
            #outmask[*, j] = (arrmask[*, j] EQ 1) OR (chi2_cap GT sigrej_eff^2)

        # Incorporate saving of each dev/sig panel onto one page? Currently only saves last fit
        #qa_plots(wavelengths, masked_fluxes, masked_vars, new_wave, new_flux, new_var)

        # Coadd anew
        spec1d = one_d_coadd(rspec, weights)
        # Calculate std_dev
        std_dev, _ = get_std_dev(rspec, spec1d)
        #var_corr = var_corr * std_dev
        msgs.info("Desired variance correction: {:g}".format(var_corr))
        msgs.info("New standard deviation: {:g}".format(std_dev))
        if do_var_corr:
            msgs.info("Correcting variance")
            for ispec in xrange(rspec.nspec):
                rspec.data['sig'][ispec] *= np.sqrt(std_dev)
            spec1d = one_d_coadd(rspec, weights)

    if iters == 0:
        msgs.warn("No iterations on coadding done")
        #qa_plots(wavelengths, masked_fluxes, masked_vars, new_wave, new_flux, new_var)
    else: #if iters > 0:
        msgs.info("Final correction to initial variances: {:g}".format(var_corr))

    # QA
    if qafile is not None:
        arqa.coaddspec_qa(spectra, rspec, spec1d, qafile=qafile)

    # Write to disk?
    if outfile is not None:
        msgs.work("Need to include header info")
        if '.hdf5' in outfile:
            spec1d.write_to_hdf5(outfile)
        elif '.fits' in outfile:
            spec1d.write_to_fits(outfile)

    return

def old_sigma_clip(fluxes, variances, sn2, n_grow_mask=1):
    """ Sigma-clips the flux arrays.

   Parameters
   ----------
   initial_mask : mask
      Initial mask for the flux + variance arrays
   n_grow_mask : int
      Number of pixels to grow the initial mask by
      on each side. Defaults to 1 pixel

   Returns
   -------
   grow_mask : mask
      Final mask for the flux + variance arrays
   """
    from functools import reduce

    first_mask = np.ma.getmaskarray(fluxes)
    highest_sn_idx = np.argmax(sn2)

    base_sharp_chi = (fluxes - fluxes[highest_sn_idx]) / (np.sqrt(variances + variances[highest_sn_idx]))

    bad_pix = []

    for row in range(0, base_sharp_chi.shape[0]):
        bad_pix.append(np.where(np.abs(base_sharp_chi[row]) > 3*np.std(base_sharp_chi, axis=1)[row])[0])

    all_bad_pix = reduce(np.union1d, (np.asarray(bad_pix)))

    for idx in range(len(all_bad_pix)):
        spec_to_mask = np.argmax(np.abs(fluxes[:, all_bad_pix[idx]]))
        #print("Masking pixel", all_bad_pix[idx], "in exposure", spec_to_mask+1
        first_mask[spec_to_mask][all_bad_pix[idx]] = True

    final_mask = grow_mask(first_mask, n_grow=n_grow_mask)

    return final_mask


def old_one_d_coadd(wavelengths, fluxes, variances, sig_clip=False, wave_grid_method=None):
    """ Performs a coadding of the spectra in 1D.

    Parameters
    ----------
    wavelengths : nd masked array
        Wavelength arrays of the input spectra
    fluxes : nd masked array
        Flux arrays of the input spectra
    variances : nd masked array
        Variances of the input spectra
    sig_clip : optional
        Perform sigma-clipping of arrays. Defaults to
        no sigma-clipping

    Returns
    -------
    fluxes : nd masked array
        Original flux arrays of the input spectra
    variances : nd masked array
        Original variances of the input spectra
    new_wave : array
        New wavelength grid
    new_flux : array
        Coadded flux array
    new_var : array
        Variance of coadded spectrum
    """
    # Generate the final wavelength grid
    new_wave = new_wave_grid(wavelengths, method=wave_grid_method)

    # Calculate S/N for each spectrum (for weighting)
    sn2, weights = sn_weight(new_wave, fluxes, variances)

    inv_variances = 1./variances

    # Rebin
    for spec in range(fluxes.shape[0]):
        obj = xspectrum1d.XSpectrum1D.from_tuple((np.ma.getdata(wavelengths[spec]),
                                                  np.ma.getdata(fluxes[spec]),
                                                  np.sqrt(np.ma.getdata(variances[spec]))))
        obj = obj.rebin(new_wave*u.AA, do_sig=True)
        fluxes[spec] = np.ma.array(obj.flux)
        variances[spec] = np.ma.array(obj.sig)**2.

    # Clip?
    if sig_clip:
        final_mask = sigma_clip(fluxes, variances, sn2)

        weights = np.ma.array(weights, mask=final_mask)
        fluxes = np.ma.array(fluxes, mask=final_mask)
        variances = np.ma.array(variances, mask=final_mask)

    else:
        final_mask = np.ma.getmaskarray(fluxes)
        weights = np.ma.array(weights, mask=final_mask)

    sum_weights = np.ma.sum(weights, axis=0)

    # Coadd
    new_flux = np.ma.sum(weights*fluxes, axis=0) / (sum_weights + (sum_weights == 0.0).astype(int))
    var = (inv_variances != 0.0).astype(float) / (inv_variances + (inv_variances == 0.0).astype(float))
    new_var = np.ma.sum((weights**2.)*var, axis=0) / ((sum_weights + (sum_weights == 0.0).astype(int))**2.)

    return fluxes, variances, new_wave, new_flux, new_var
