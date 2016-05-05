import numpy as np
from scipy import ndimage
from scipy import interpolate

from linetools.spectra import xspectrum1d
from astropy.io import fits
from astropy import units as u

import ararc as aarc
import armsgs
import arutils

# Logging
msgs = armsgs.get_logger()

try:
    from xastropy.xutils import xdebug as debugger
except:
    import pdb as debugger

def flexure(slf, det):
    """Correct wavelengths for flexure

    Parameters:
    ----------
    slf :
    det : int
      detector number

    Returns:
    ----------
    flex_dict: dict
      dict on flexure results

    """
    #Search for appropriate archived sky spectrum based on latitude, longitude
    #****loop over latitude and longitudes in column indices 1, 2 and find appropriate archived sky spectrum
    #****save name of archived sky spectrum (column index 3) as a string and open it with fits.open
    #****column index 0 is location name maybe? observatory name?
 #   latitude = slf._spect['mosaic']['latitude']
 #   longitude = slf._spect['mosaic']['longitude']
 #   sky_files = Table.read(~/PYPIT/data/sky_files/sky_spectra_table, format='ascii')
 #   archive_wave =
 #   archive_flux =

    # Read archived sky
    root = slf._argflag['run']['pypitdir']
    #root = '/Users/xavier/local/Python/PYPIT'
    #skyspec_fil = 'sky_LRISr_600_7500_5460_7950.fits'
    skyspec_fil = 'paranal_sky.fits'
    msgs.info("Using {:s} file for Sky spectrum".format(skyspec_fil))
    hdu = fits.open(root+'/data/sky_spec/'+skyspec_fil)
    archive_wave = hdu[0].data
    archive_flux = hdu[1].data
    arx_sky = xspectrum1d.XSpectrum1D.from_tuple((archive_wave, archive_flux))

    # Loop on objects
    flex_dict = dict(polyfit=[], roots=[], shift=[], subpix=[], corr=[],
                     spec_file=skyspec_fil, smooth=[])

    for specobj in slf._specobjs[det-1]:  # for convenience

        # Using boxcar
        if slf._argflag['reduce']['flexure'] == 'boxcar':
            sky_wave = specobj.boxcar['wave'].to('AA').value
            sky_flux = specobj.boxcar['sky']
        else:
            msgs.error("Not ready for this flexure method: {}".format(
                    slf._argflag['reduce']['flexure']))

        # Generate 1D spectrum for object
        obj_sky = xspectrum1d.XSpectrum1D.from_tuple((sky_wave, sky_flux))

        #Determine the brightest emission lines
        msgs.warn("If we use Paranal, cut down on wavelength early on")
        arx_amp, arx_cent, arx_wid, arx_w, arx_satsnd, arx_yprep = aarc.detect_lines(slf, det, msarc=None, censpec=arx_sky.flux.value, MK_SATMASK=False)
        obj_amp, obj_cent, obj_wid, obj_w, obj_satsnd, obj_yprep = aarc.detect_lines(slf, det, msarc=None, censpec=obj_sky.flux.value, MK_SATMASK=False)

        #Keep only 5 brightest amplitude lines (xxx_keep is array of indices within arx_w of the 5 brightest)
        arx_keep = np.argsort(arx_amp[arx_w])[-5:]
        obj_keep = np.argsort(obj_amp[obj_w])[-5:]

        #Calculate wavelength (Angstrom per pixel)
        arx_disp = (np.amax(arx_sky.wavelength.value)-np.amin(arx_sky.wavelength.value))/arx_sky.wavelength.size
        obj_disp = (np.amax(obj_sky.wavelength.value)-np.amin(obj_sky.wavelength.value))/obj_sky.wavelength.size

        #Calculate resolution (lambda/delta lambda_FWHM)..maybe don't need this? can just use sigmas
        arx_res = (arx_sky.wavelength.value[0]+(arx_disp*arx_cent[arx_w][arx_keep]))/(
            arx_disp*(2*np.sqrt(2*np.log(2)))*arx_wid[arx_w][arx_keep])
        obj_res = (obj_sky.wavelength.value[0]+(obj_disp*obj_cent[obj_w][obj_keep]))/(
            obj_disp*(2*np.sqrt(2*np.log(2)))*obj_wid[obj_w][obj_keep])
        msgs.info("Resolution of Archive={:g} and Observation={:g}".format(
            np.median(arx_res), np.median(obj_res)))

        #Determine sigma of gaussian for smoothing
        arx_sig = (arx_disp*arx_wid[arx_w][arx_keep])**2.
        obj_sig = (obj_disp*obj_wid[obj_w][obj_keep])**2.

        arx_med_sig = np.median(arx_sig)
        obj_med_sig = np.median(obj_sig)

        if obj_med_sig >= arx_med_sig:
            smooth_sig = np.sqrt(obj_med_sig-arx_med_sig)
        else:
            msgs.warn("Prefer archival sky spectrum to have higher resolution")
            smooth_sig = 0.
            #smooth_sig = np.sqrt(arx_med_sig-obj_med_sig)

        #Determine region of wavelength overlap
        min_wave = max(np.amin(arx_sky.wavelength.value), np.amin(obj_sky.wavelength.value))
        max_wave = min(np.amax(arx_sky.wavelength.value), np.amax(obj_sky.wavelength.value))

        #Smooth higher resolution spectrum by smooth_sig (flux is conserved!)
        if np.median(obj_res) >= np.median(arx_res):
            msgs.warn("New Sky has higher resolution than Archive.  Not smoothing")
            #obj_sky_newflux = ndimage.gaussian_filter(obj_sky.flux, smooth_sig)
        else:
            arx_sky.flux = ndimage.gaussian_filter(arx_sky.flux, smooth_sig)

        # Define wavelengths of overlapping spectra
        keep_idx = np.where((obj_sky.wavelength.value>=min_wave) &
                             (obj_sky.wavelength.value<=max_wave))[0]
        #keep_wave = [i for i in obj_sky.wavelength.value if i>=min_wave if i<=max_wave]

        #Rebin both spectra onto overlapped wavelength range
        if len(keep_idx) <= 50:
            msgs.error("Not enough overlap between sky spectra")

        else: #rebin onto object ALWAYS
            keep_wave = obj_sky.wavelength[keep_idx]
            arx_sky = arx_sky.rebin(keep_wave)
            obj_sky = obj_sky.rebin(keep_wave)

        '''
        if msgs._debug['flexure']:
            debugger.xplot(arx_sky.wavelength, arx_sky.flux, xtwo=obj_sky.wavelength, ytwo=obj_sky.flux)
            debugger.xplot(arx_sky.wavelength, arx_sky.flux, xtwo=np.roll(obj_sky.wavelength.value,9), ytwo=obj_sky.flux*100)
            debugger.set_trace()
        '''
        #deal with bad pixels
        msgs.work("Need to mask bad pixels")

        #deal with underlying continuum
        msgs.work("Need to deal with underlying continuum")

        #Cross correlation of spectra
        corr = np.correlate(arx_sky.flux, obj_sky.flux, "same")

        #Create array around the max of the correlation function for fitting for subpixel max
        max_corr = np.argmax(corr)
        subpix_grid = np.linspace(max_corr-3., max_corr+3., 7.)

        #Fit a 2-degree polynomial to peak of correlation function
        fit = np.polynomial.polynomial.polyfit(subpix_grid, corr[subpix_grid.astype(np.int)], 2)
        roots = np.roots([fit[2], fit[1], fit[0]])
        max_fit = (roots[0]+roots[1])/2.
        #   fit, other = arutils.gauss_lsqfit(subpix_grid, corr[subpix_grid.astype(np.int)], max_corr)

        #Calculate and apply shift in wavelength
        shift = max_fit-(corr.size/2)
        msgs.info("Flexure correction of {:g} pixels".format(shift))
        #model = (fit[2]*(subpix_grid**2.))+(fit[1]*subpix_grid)+fit[0]

        # Simple interpolation
        npix = len(sky_wave)
        x = np.linspace(0., 1., npix)
        # Apply
        for attr in ['boxcar', 'optimal']:
            if not hasattr(specobj,attr):
                continue
            if 'wave' in getattr(specobj, attr).keys():
                msgs.info("Applying flexure correction to {:s} extraction for object {:s}".format(
                    attr, str(specobj)))
                f = interpolate.interp1d(x, specobj.boxcar['wave'], bounds_error=False, fill_value="extrapolate")
                getattr(specobj, attr)['wave'] = f(x+shift/(npix-1))

        # Update dict
        flex_dict['polyfit'].append(fit)
        flex_dict['roots'].append(roots)
        flex_dict['shift'].append(shift)
        flex_dict['subpix'].append(subpix_grid)
        flex_dict['corr'].append(corr[subpix_grid.astype(np.int)])
        flex_dict['smooth'].append(smooth_sig)

    return flex_dict

def airtovac(wave):
    """Convert air-based wavelengths to vacuum

    Parameters:
    ----------
    wave: Quantity array
      Wavelengths 

    Returns:
    ----------
    wave: Quantity array
      Wavelength array corrected to vacuum wavelengths
    """
    # Convert to AA
    wave = wave.to(u.AA)
    wavelength = wave.value

    # Standard conversion format
    sigma_sq = (1.e4/wavelength)**2. #wavenumber squared
    factor = 1 + (5.792105e-2/(238.0185-sigma_sq)) + (1.67918e-3/(57.362-sigma_sq))
    factor = factor*(wavelength>=2000.) + 1.*(wavelength<2000.) #only modify above 2000A

    # Convert
    wavelength = wavelength*factor
    # Units
    new_wave = wavelength*u.AA
    new_wave.to(wave.unit)

    return new_wave

def vactoair(wave):
    """Convert to air-based wavelengths from vacuum

    Parameters:
    ----------
    wave: Quantity array
      Wavelengths 

    Returns:
    ----------
    wave: Quantity array
      Wavelength array corrected to air
    """
    # Convert to AA
    wave = wave.to(u.AA)
    wavelength = wave.value

    # Standard conversion format
    sigma_sq = (1.e4/wavelength)**2. #wavenumber squared
    factor = 1 + (5.792105e-2/(238.0185-sigma_sq)) + (1.67918e-3/(57.362-sigma_sq))
    factor = factor*(wavelength>=2000.) + 1.*(wavelength<2000.) #only modify above 2000A

    # Convert
    wavelength = wavelength/factor
    new_wave = wavelength*u.AA
    new_wave.to(wave.unit)

    return new_wave