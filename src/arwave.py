import numpy as np
from linetools.spectra import xspectrum1d
from xastropy.utils import xdebug as xdebug
from astropy.io import fits
from astropy import units as u

def flexure():
	archived_sky=fits.open("") #going to eventually search for appropriate sky file
	object_sky=fits.open("")

	archived_sky_wave=archived_sky[0].data
	archived_sky_flux=archived_sky[1].data

	#
	object_sky_wave=object_sky[0].data
	object_sky_flux=object_sky[1].data

