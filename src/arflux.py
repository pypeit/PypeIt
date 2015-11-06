# Module for fluxing routines
import numpy as np
import scipy
import glob

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u
from astropy import coordinates as coords

import armsgs as msgs
import arcyextract
import arcyutils
import arcyproc
import arload
import artrace
import arutils
import arplot

try:
    from xastropy.xutils import xdebug as xdb
except:
    pass

def extinction_correction(wave, AM, extinct):
    '''Derive extinction correction
    Based on algorithm in LowRedux (long_extinct)

    Parameters:
    ----------
    wave: Quantity array
      Wavelengths for interpolation
    AM: float
      Airmass
    extinct: Table
      Table of extinction values

    Returns:
    ----------
    flux_corr: ndarray
      Flux corrections at the input wavelengths
    '''
    # Interpolate
    f_mag_ext = scipy.interpolate.interp1d(extinct['wave'],
        extinct['mag_ext'], bounds_error=False, fill_value=0.)
    mag_ext = f_mag_ext(wave.to('AA').value)

    # Deal with outside wavelengths 
    gdv = np.where(mag_ext > 0.)[0]
    if gdv[0] != 0.: # Low wavelengths
        mag_ext[0:gdv[0]] = mag_ext[gdv[0]]
        msgs.warn("Extrapolating at low wavelengths using last valid value")
    if gdv[-1] != 0.: # High wavelengths
        mag_ext[gdv[-1]+1:] = mag_ext[gdv[-1]]
        msgs.warn("Extrapolating at high wavelengths using last valid value")
    # Evaluate
    flux_corr = 10.0**(0.4*mag_ext*AM)
    # Return
    return flux_corr

def find_standard_file(slf, radec, toler=20.*u.arcmin):
    '''Find a match for the input file to one of the archived 
    standard star files (hopefully).  Priority is by order of search.

    Parameters:
    ----------
    radec: tuple
      ra, dec in string format ('05:06:36.6','52:52:01.0')

    Returns:
    --------
    sdict: dict
      'file': str -- Filename
      'fmt': int -- Format flag
           1=Calspec style FITS binary table
      'name': str -- Star name
      'ra': str -- RA(2000)
      'dec': str -- DEC(2000)
    '''
    # Priority
    std_sets = [load_calspec]
    std_file_fmt = [1] # 1=Calspec style FITS binary table

    # SkyCoord
    obj_coord = SkyCoord(radec[0], radec[1], unit=(u.hourangle, u.deg))

    # Loop on standard sets
    for qq,sset in enumerate(std_sets):
        # Stars
        path, star_tbl = sset(slf)
        star_coords = SkyCoord(star_tbl['RA_2000'], star_tbl['DEC_2000'], 
            unit=(u.hourangle, u.deg))
        # Match
        idx, d2d, d3d = coords.match_coordinates_sky(obj_coord, star_coords, nthneighbor=1)
        if d2d < toler:
            # Generate a dict
            std_dict = dict(file=path+star_tbl[int(idx)]['File'], 
                name=star_tbl[int(idx)]['Name'], fmt=std_file_fmt[qq],
                ra=star_tbl[int(idx)]['RA_2000'], 
                dec=star_tbl[int(idx)]['DEC_2000'])
            # Return
            return std_dict

def load_calspec(slf):
    ''' Load the list of calspec standards

    Parameters:
    ----------

    Returns:
    --------
    calspec_path: str
      Path from pypitdir to calspec standard star files
    calspec_stds: Table
      astropy Table of the calspec standard stars (file, Name, RA, DEC)
    '''
    # Read
    calspec_path = '/data/standards/calspec/'
    calspec_file = slf._argflag['run']['pypitdir'] + calspec_path + 'calspec_info.txt'
    calspec_stds = Table.read(calspec_file,comment='#',format='ascii')
    # Return
    return calspec_path, calspec_stds

def load_extinction_data(slf, toler=1.*u.deg):
    '''Find the best extinction file to use, based on longitude and latitude
    Loads it and returns a Table

    Parameters:
    ----------
    slf: 
      Includes mosaic lon/lat
    toler: Angle
      Tolerance for matching detector to site (1 deg)

    Returns:
    ----------
    ext_file: Table
      astropy Table containing the 'wavelength', 'extinct' data for AM=1.
    '''
    # Mosaic coord
    mosaic_coord = SkyCoord(slf._spect['mosaic']['longitude'],
        slf._spect['mosaic']['latitude'], frame='gcrs', unit=u.deg)
    # Read list
    extinct_path = slf._argflag['run']['pypitdir']+'/data/extinction/'
    extinct_summ = extinct_path+'README'
    extinct_files = Table.read(extinct_summ,comment='#',format='ascii')
    # Coords
    ext_coord = SkyCoord(extinct_files['Lon'], extinct_files['Lat'], frame='gcrs', unit=u.deg)
    # Match
    idx, d2d, d3d = coords.match_coordinates_sky(mosaic_coord, ext_coord, nthneighbor=1)
    if d2d < toler:
        extinct_file = extinct_files[int(idx)]['File']
        msgs.info("Using {:s} for extinction corrections.".format(extinct_file))
    else:
        extinct_file = 'atm_trans_am1.0.dat'
        msgs.warn("Using {:s} for extinction corrections.".format(extinct_file))
        msgs.warn("You may wish to generate a site-specific file")
    # Read
    extinct = Table.read(extinct_path+extinct_file,comment='#',format='ascii', names=('wave','mag_ext'))
    # Return
    return extinct


def load_standard_file(slf, std_dict):
    '''Load standard star data

    Parameters:
    ----------
    std_dict: dict
      Info on standard star indcluding filename in 'file'
      May be compressed

    Returns:
    --------
    std_wave: Quantity array
      Wavelengths of standard star array
    std_flux: Quantity array
      Flux of standard star
    '''
    fil = glob.glob(slf._argflag['run']['pypitdir']+
            std_dict['file']+'*')
    if len(fil) == 0:
        msgs.error("No standard star file: {:s}".format(fil))
    else:
        fil = fil[0]
        msgs.info("Loading standard star file: {:s}".format(fil))

    if std_dict['fmt'] == 1:
        std_spec = fits.open(fil)[1].data
        # Return
        return std_spec['WAVELENGTH']*u.AA, 1e17*std_spec['FLUX']*u.erg/u.s/u.cm**2/u.AA
    else:
        msgs.error("Bad Standard Star Format")


