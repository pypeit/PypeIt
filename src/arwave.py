import numpy as np
import pylab as plt
from scipy.optimize import newton
from linetools.spectra import xspectrum1d
from astropy.io import fits
#from xastropy.utils import xdebug as xdebug
from astropy import units as u

import armsgs
import arutils

#testing on my data
'''
object_sky=fits.open("/Users/tiffanyhsyu/Dropbox/DropboxData/05192015_Lick/blue_side/combined_frames/extract/kast_sky_blue_fluxed.fits")
object_wave1=np.arange(0.,2048.,1.)
object_wave2=3428.12+(1.0196*object_wave1)
object_flux=object_sky[0].data
'''

def flexure(obj_wave, obj_flux):
    #look for appropriate archived sky file based on latitude and longitude
    #if archived sky file doesn't exist, send out error
        # --> find out where lat+long are saved
    archive_sky=fits.open("/Users/tiffanyhsyu/Dropbox/XMP/Kast_Exmpl/kast_sky_blue_600.fits")
    archive_wave=archive_sky[0].data
    archive_flux=archive_sky[1].data

    arx_sky=xspectrum1d.XSpectrum1D.from_tuple((archive_wave,archive_flux))
    obj_sky=xspectrum1d.XSpectrum1D.from_tuple((obj_wave,obj_flux))

    #find out resolution of both instruments; assign high and low resolution

    #find brightest sky lines in each of the spectrum and fit gaussians to them
    #find sigma_high^2, sigma_low^2, for those lines, and get a median of sigma
    #then find sigma^2=sigma_low^2-sigma_high^2 and this will be the sigma of the
    #gaussian that we want to smooth the higher resolution spectra (onto lower res.)

    


    #make sure that the two spectra overlap a good amount on wavelength range
    #send error+quit if overlap is less than 10 pixels? what is a good amount to
    #quit at?
    #rebin onto overlapped wavelength range
    min_wave=max(np.amin(arx_sky.dispersion.value),np.amin(obj_sky.dispersion.value))
    max_wave=min(np.amax(arx_sky.dispersion.value),np.amax(obj_sky.dispersion.value))
    print min_wave, max_wave

    arx_res=(np.amax(arx_sky.dispersion.value)-np.amin(arx_sky.dispersion.value))/len(arx_sky.dispersion)
    obj_res=(np.amax(obj_sky.dispersion.value)-np.amin(obj_sky.dispersion.value))/len(obj_sky.dispersion)
    print arx_res, obj_res

    if obj_res <= arx_res: #object has lower wavelength per pixel-->higher resolution
        keep_wave=[i for i in arx_sky.dispersion.value if i>=min_wave if i<=max_wave]
        obj_sort=np.argsort(obj_sky.flux.value)
        print obj_sort

    else:
        keep_wave=[i for i in obj_sky.dispersion.value if i>=min_wave if i<=max_wave]
    #keep_wave=[i for i in obj_sky.dispersion.value if i>=min_wave if i<=max_wave]
    #why use obj_sky.dispersion.value here? over arx_sky.dispersion.value? maybe
    #should check which one has finer resolution over the overlapped range

    if len(keep_wave) <= 10:
        msgs.error("Not enough overlap between sky spectra")

    else:
        arx_sky=arx_sky.rebin(keep_wave*u.AA)
        obj_sky=obj_sky.rebin(keep_wave*u.AA)

    #deal with bad pixels

    #deal with underlying continuum


    #rebin the spectra onto the same pixel-wavelength scale
#    n_arx=len(arx_sky.dispersion)
#    n_obj=len(obj_sky.dispersion)
#
#    if n_arx > n_obj:
#        obj_sky=obj_sky.rebin(arx_sky.dispersion)
#        print "n_archive > n_obj, rebinning object to", len(obj_sky.dispersion)
#
#    elif n_arx < n_obj:
#        print arx_sky.dispersion
#        arx_sky=arx_sky.rebin(obj_sky.dispersion)
#        print "n_obj > n_archive, rebinning archive to", len(arx_sky.dispersion)

    #else:
        # tell code to continue when archived and object wavelengths are the
        # same length? do I need to?


    #cross correlate the two spectra, which are now on the same wavelength
    #range and have the same array sizes
    corr=np.correlate(arx_sky.flux,obj_sky.flux,"same")
    plt.plot(corr)

    #create an array around the max of the correlation function
    max_corr=np.argmax(corr)
    subpix_grid=np.linspace(max_corr-4,max_corr+4,9.)

    #fit a gaussian to where the correlation function peaks
    fit,other=arutils.gauss_lsqfit(subpix_grid,corr[subpix_grid.astype(np.int)],max_corr)

    #calculate the shift and apply it to the wavelength
    shift=fit[1]-0.
    #model=fit[0]*(np.exp((-(subpix_grid-fit[1])**2)/(2*fit[2]**2)))
    #plt.plot(subpix_grid,model)

    finer_subpix_grid=np.linspace(max_corr-4,max_corr+4,90.)
    model=fit[0]*(np.exp((-(finer_subpix_grid-fit[1])**2)/(2*fit[2]**2)))
    plt.plot(finer_subpix_grid,model)
    #plt.show()

    return

#flexure(object_wave2,object_flux)

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
    '''Convert to air-based wavelengths from vacuum

    Parameters:
    ----------
    wave: Quantity array
      Wavelengths 

    Returns:
    ----------
    wave: Quantity array
      Wavelength array corrected to air
    '''
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