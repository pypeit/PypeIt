import numpy as np
import pylab as plt
from linetools.spectra import xspectrum1d
from astropy.io import fits
from astropy import units as u

import arutils


def flexure(slf, det, obj_wave, obj_flux):
    '''Correct for flexure

    Parameters:
    ----------
    obj_wave: Quantity array
    obj_flux: Quantity array
      Wavelength and flux

    Returns:
    ----------
    shift: Float
      Pixel shift to be applied to object

    '''

    #Search for appropriate archived sky spectrum based on latitude, longitude
    #****loop over latitude and longitudes in column indices 1, 2 and find appropriate archived sky spectrum
    #****save name of archived sky spectrum (column index 3) as a string and open it with fits.open
    #****column index 0 is location name maybe? observatory name?
 #   latitude = slf._spect['mosaic']['latitude']
 #   longitude = slf._spect['mosaic']['longitude']
 #   sky_files = Table.read(~/PYPIT/data/sky_files/sky_spectra_table, format='ascii')
 #   archive_wave =
 #   archive_flux =

    #****if archived sky file doesn't exist, send out error
    #****find out where lat+long are saved
 #   archive_sky = fits.open("/Users/tiffanyhsyu/Dropbox/XMP/Kast_Exmpl/kast_sky_blue_600.fits")
 #   archive_wave = archive_sky[0].data
 #   archive_flux = archive_sky[1].data

    #using object as fake archived with known shift..see if shift returns the correct value
    archive_sky=fits.open("/Users/tiffanyhsyu/Dropbox/DropboxData/05192015_Lick/blue_side/combined_frames/extract/kast_sky_blue_fluxed.fits")
    archive_wave1=np.arange(1.5, 2049.5, 1.)
    archive_wave=3428.12+(1.0196*archive_wave1)
    archive_flux=archive_sky[0].data

    #****eventually will need this instead since input is Quantity array
 #   obj_wave = obj_wave.to(u.AA)
 #   arx_sky = xspectrum1d.XSpectrum1D.from_tuple((archive_wave.value, archive_flux.value))
 #   obj_sky = xspectrum1d.XSpectrum1D.from_tuple((obj_wave.value, obj_flux.value))

    arx_sky = xspectrum1d.XSpectrum1D.from_tuple((archive_wave, archive_flux))
    obj_sky = xspectrum1d.XSpectrum1D.from_tuple((obj_wave, obj_flux))

    #Determine the brightest emission lines
    arx_amp, arx_cent, arx_wid, arx_w, arx_satsnd, arx_yprep = ararc.detect_lines(slf, det, msarc=None, censpec=arx_sky.flux.value, MK_SATMASK=False)
    obj_amp, obj_cent, obj_wid, obj_w, obj_satsnd, obj_yprep = ararc.detect_lines(slf, det, msarc=None, censpec=obj_sky.flux.value, MK_SATMASK=False)

    #Keep only 5 brightest amplitude lines (xxx_keep is array of indices within arx_w of the 5 brightest)
    arx_keep = np.argsort(arx_amp[arx_w])[-5:]
    obj_keep = np.argsort(obj_amp[obj_w])[-5:]

    #Calculate dispersion (Angstrom per pixel)
    arx_disp = (np.amax(arx_sky.dispersion.value)-np.amin(arx_sky.dispersion.value))/arx_sky.dispersion.size
    obj_disp = (np.amax(obj_sky.dispersion.value)-np.amin(obj_sky.dispersion.value))/obj_sky.dispersion.size

    #Calculate resolution (lambda/delta lambda_FWHM)..maybe don't need this? can just use sigmas
    arx_res = (arx_sky.dispersion.value[0]+(arx_disp*arx_cent[arx_w][arx_keep]))/(arx_disp*(2*np.sqrt(2*np.log(2)))*arx_wid[arx_w][arx_keep])
    obj_res = (obj_sky.dispersion.value[0]+(obj_disp*obj_cent[obj_w][obj_keep]))/(obj_disp*(2*np.sqrt(2*np.log(2)))*obj_wid[obj_w][obj_keep])

    #Determine sigma of gaussian for smoothing <-- Do I need to convert these to wavelengths?
    arx_sig = (arx_disp*arx_wid[arx_w][arx_keep])**2.
    obj_sig = (obj_disp*obj_wid[obj_w][obj_keep])**2.

    arx_med_sig = np.median(arx_sig)
    obj_med_sig = np.median(obj_sig)

    if arx_med_sig >= obj_med_sig:
        smooth_sig = np.sqrt(arx_med_sig-obj_med_sig)

    else:
        smooth_sig = np.sqrt(obj_med_sig-arx_med_sig)

    #Determine region of wavelength overlap
    min_wave = max(np.amin(arx_sky.dispersion.value), np.amin(obj_sky.dispersion.value))
    max_wave = min(np.amax(arx_sky.dispersion.value), np.amax(obj_sky.dispersion.value))

    #Smooth higher resolution spectrum by smooth_sig (flux is conserved!)
    if np.median(obj_res) >= np.median(arx_res):
        obj_sky_newflux = ndimage.gaussian_filter(obj_sky.flux, smooth_sig)
#       keep_wave = [i for i in arx_sky.dispersion.value if i>=min_wave if i<=max_wave]

    else:
        arx_sky.flux = ndimage.gaussian_filter(arx_sky.flux, smooth_sig)

    keep_wave = [i for i in obj_sky.dispersion.value if i>=min_wave if i<=max_wave]

#    xdb.set_trace() can plot smoothed spectrum over unsmoothed here with:
#    xdb.xplot(obj_sky.flux, xtwo=np.arange(0.,2048.,1.), ytwo=obj_sky_newflux)

    #Rebin both spectra onto overlapped wavelength range
    if len(keep_wave) <= 50:
        msgs.error("Not enough overlap between sky spectra")

    else: #rebin onto object ALWAYS
        arx_sky = arx_sky.rebin(keep_wave*u.AA)
        obj_sky = obj_sky.rebin(keep_wave*u.AA)

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
    model = (fit[2]*(subpix_grid**2.))+(fit[1]*subpix_grid)+fit[0]

#    finer_subpix_grid = np.linspace(max_corr-4,max_corr+4,90.)
#    model2 = (fit[2]*(finer_subpix_grid**2.))+(fit[1]*finer_subpix_grid)+fit[0]
#    model2 = fit[0]*(np.exp((-(finer_subpix_grid-fit[1])**2)/(2*fit[2]**2)))

    #QA plot for cross correlation
 #   axarr[2].plot(subpix_grid, corr[subpix_grid.astype(np.int)])
 #   axarr[2].plot(subpix_grid, model)
 #   axarr[2].set_title(Cross Correlation)
 #   fig.savefig(flexure_correction.pdf, format='pdf')
 #   fig.savefig(os.path.join(path,flexure_lines.pdf), format='pdf')
    debugger.set_trace()

    return shift

def airtovac(wave):
    '''Convert air-based wavelengths to vacuum

    Parameters:
    ----------
    wave: Quantity array
      Wavelengths 

    Returns:
    ----------
    wave: Quantity array
      Wavelength array corrected to vacuum wavelengths
    '''
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