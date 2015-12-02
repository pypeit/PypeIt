# Module for fluxing routines
import numpy as np
import scipy
import glob

from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u
from astropy import coordinates as coords

import armsgs
import arutils

try:
    from linetools.spectra.xspectrum1d import XSpectrum1D
except:
    pass


try:
    from xastropy.xutils import xdebug as xdb
except:
    pass

# Logging
msgs = armsgs.get_logger()

def apply_sensfunc(slf, sc, MAX_EXTRAP=0.05):
    """
    Apply the sensitivity function to the data
    We also correct for extinction.

    Parameters
    ----------
    MAX_EXTRAP : float, optional [0.05]
      Fractional amount to extrapolate sensitivity function
    """

    # Load extinction data
    extinct = load_extinction_data(slf)
    airmass = slf._fitsdict['airmass'][slf._scidx]
    # Loop on objects
    for spobj in slf._specobjs:
        # Loop on extraction modes
        for extract_type in ['boxcar']:
            try:
                extract = getattr(spobj,extract_type)
            except AttributeError:
                continue
            msgs.info("Fluxing {:s} extraction".format(extract_type))
            wave = extract['wave'] # for convenience
            scale = np.zeros(wave.size)
            # Allow for some extrapolation 
            dwv = slf._sensfunc['wave_max']-slf._sensfunc['wave_min']
            inds = (wave >= slf._sensfunc['wave_min']-dwv*MAX_EXTRAP) & (wave <= slf._sensfunc['wave_max']+dwv*MAX_EXTRAP)
            mag_func = arutils.func_val(slf._sensfunc['c'], wave[inds], slf._sensfunc['func'])
            sens = 10.0**(0.4*mag_func)
            # Extinction
            ext_corr = extinction_correction(wave[inds],airmass,extinct)
            scale[inds] = sens*ext_corr
            # Fill
            extract['flam'] = extract['counts']*scale/slf._fitsdict['exptime'][slf._scidx]
            extract['flam_var'] = extract['var']*(scale/slf._fitsdict['exptime'][slf._scidx])**2


def bspline_magfit(wave, flux, var, flux_std, nointerp=False, **kwargs):
    """
    Perform a bspline fit to the flux ratio of standard to
    observed counts.  Used to generate a sensitivity function.

    Parameters
    ----------
    wave : ndarray
    flux : ndarray
      counts/s as observed
    invvar : ndarray
      inverse variance
    flux_std : Quantity array
      standard star true flux (erg/s/cm^2/A)
    nointer : bool, optional [False]
      Skip interpolation over bad points (not recommended)?
    **kwargs : keywords for robust_polyfit

    Returns
    -------
    """
    invvar = (var > 0.)/(var + (var <= 0.))
    nx = wave.size
    pos_error = 1./np.sqrt(np.maximum(invvar,0.) + (invvar == 0))
    pos_mask = (flux > pos_error/10.0) & (invvar > 0) & (flux_std > 0.0)
    #pos = pos_mask==1 npos)

    fluxlog = 2.5*np.log10(np.maximum(flux,pos_error/10))
    logivar = invvar * flux**2 * pos_mask*1.08574

    # cap the magfunc so that sensfunc < 1.0e10
    magfunc = 2.5*np.log10(np.maximum(flux_std,1.0e-2)) - fluxlog
    magfunc = np.minimum(magfunc,25.0)
    sensfunc = 10.0**(0.4*magfunc)*pos_mask

    # Interpolate over masked pixels
    if not nointerp:
        bad = logivar <= 0. 
        if np.sum(bad) > 0:
            f = scipy.interpolate.InterpolatedUnivariateSpline(wave[~bad], magfunc[~bad], k=2)
            magfunc[bad] = f(wave[bad])
            fi = scipy.interpolate.InterpolatedUnivariateSpline(wave[~bad], logivar[~bad], k=2)
            logivar[bad] = fi(wave[bad])

    #  First iteration
    mask, tck = arutils.robust_polyfit(wave, magfunc, 3, function='bspline', weights=np.sqrt(logivar), **kwargs)
    logfit1 = arutils.func_val(tck,wave,'bspline')
    modelfit1 = 10.0**(0.4*(logfit1))

    residual = sensfunc/(modelfit1 + (modelfit1 == 0)) - 1.
    new_mask = pos_mask & (sensfunc > 0)
    residual_ivar = (modelfit1*flux/(sensfunc + (sensfunc == 0.0)))**2*invvar
    residual_ivar = residual_ivar*new_mask

    # Interpolate over masked pixels
    if not nointerp:
        if np.sum(bad) > 0:
            f = scipy.interpolate.InterpolatedUnivariateSpline(wave[~bad], residual[~bad], k=2)
            residual[bad] = f(wave[bad])
            fi = scipy.interpolate.InterpolatedUnivariateSpline(wave[~bad], residual_ivar[~bad], k=2)
            residual_ivar[bad] = fi(wave[bad])

    #  Now do one more fit to the ratio of data/model - 1.
    # Fuss with the knots first ()
    inner_knots = arutils.bspline_inner_knots(tck[0])
    #
    mask, tck_residual = arutils.robust_polyfit(wave, residual, 3, function='bspline', weights=np.sqrt(residual_ivar), knots=inner_knots, **kwargs)
    if tck_residual[1].size != tck[1].size:
        msgs.error('Problem with bspline knots in bspline_magfit')
    #bset_residual = bspline_iterfit(wave, residual, weights=np.sqrt(residual_ivar), knots = tck[0])

    tck_log1 = list(tck)
    tck_log1[1] = tck[1] + tck_residual[1]

    sensfit = 10.0**(0.4*(arutils.func_val(tck_log1,wave, 'bspline')))

    absdev = np.median(np.abs(sensfit/modelfit1-1))
    msgs.info('Difference between fits is {:g}'.format(absdev))

    # QA
    msgs.work("Add QA for sensitivity function")

    return tck_log1

def extinction_correction(wave, airmass, extinct):
    """
    Derive extinction correction
    Based on algorithm in LowRedux (long_extinct)

    Parameters
    ----------
    wave : Quantity array
      Wavelengths for interpolation. Should be sorted
    airmass : float
      Airmass
    extinct : Table
      Table of extinction values

    Returns
    -------
    flux_corr : ndarray
      Flux corrections at the input wavelengths
    """
    # Checks
    if airmass < 1.:
        msgs.error("Bad airmass value in extinction_correction")
    # Interpolate
    f_mag_ext = scipy.interpolate.interp1d(extinct['wave'],
        extinct['mag_ext'], bounds_error=False, fill_value=0.)
    mag_ext = f_mag_ext(wave.to('AA').value)

    # Deal with outside wavelengths 
    gdv = np.where(mag_ext > 0.)[0]
    if len(gdv) == 0:
        msgs.error("None of the input wavelengths are in the extinction correction range.  Presumably something was input wrong.")
    if gdv[0] != 0: # Low wavelengths
        mag_ext[0:gdv[0]] = mag_ext[gdv[0]]
        msgs.warn("Extrapolating at low wavelengths using last valid value")
    if gdv[-1] != (mag_ext.size-1): # High wavelengths
        mag_ext[gdv[-1]+1:] = mag_ext[gdv[-1]]
        msgs.warn("Extrapolating at high wavelengths using last valid value")
    # Evaluate
    flux_corr = 10.0**(0.4*mag_ext*airmass)
    # Return
    return flux_corr


def find_standard_file(argflag, radec, toler=20.*u.arcmin, check=False):
    """
    Find a match for the input file to one of the archived
    standard star files (hopefully).  Priority is by order of search.

    Parameters
    ----------
    argflag : dict
      Arguments and flags used for reduction
    radec : tuple
      ra, dec in string format ('05:06:36.6','52:52:01.0')
    toler : Angle
      Tolerance on matching archived standards to input
    check : bool
      If True, the routine will only check to see if a
      standard star exists within the input ra, dec, and toler range.

    Returns
    -------
    sdict : dict
      'file': str -- Filename
      'fmt': int -- Format flag
           1=Calspec style FITS binary table
      'name': str -- Star name
      'ra': str -- RA(2000)
      'dec': str -- DEC(2000)
    """
    # Priority
    std_sets = [load_calspec]
    std_file_fmt = [1]  # 1=Calspec style FITS binary table

    # SkyCoord
    obj_coord = SkyCoord(radec[0], radec[1], unit=(u.hourangle, u.deg))

    # Loop on standard sets
    closest = dict(sep=999*u.deg)
    for qq,sset in enumerate(std_sets):
        # Stars
        path, star_tbl = sset(argflag)
        star_coords = SkyCoord(star_tbl['RA_2000'], star_tbl['DEC_2000'], 
            unit=(u.hourangle, u.deg))
        # Match
        idx, d2d, d3d = coords.match_coordinates_sky(obj_coord, star_coords, nthneighbor=1)
        if d2d < toler:
            if check: return True
            else:
                # Generate a dict
                std_dict = dict(file=path+star_tbl[int(idx)]['File'],
                    name=star_tbl[int(idx)]['Name'], fmt=std_file_fmt[qq],
                    ra=star_tbl[int(idx)]['RA_2000'],
                    dec=star_tbl[int(idx)]['DEC_2000'])
                # Return
                msgs.info("Using standard star {:s}".format(std_dict['name']))
                return std_dict
        else: # Save closest, if it is
            imind2d = np.argmin(d2d)
            mind2d = d2d[imind2d]
            if mind2d < closest['sep']:
                closest['sep'] = mind2d
                closest.update(dict(name=star_tbl[int(idx)]['Name'], 
                    ra=star_tbl[int(idx)]['RA_2000'], 
                    dec=star_tbl[int(idx)]['DEC_2000']))
    # Standard star not found
    if check: return False
    msgs.warn("No standard star was found within a tolerance of {:g}".format(toler))
    msgs.info("Closest standard was {:s} at separation {:g}".format(closest['name'],closest['sep'].to('arcmin')))
    msgs.warn("Flux calibration will not be performed")
    return None


def load_calspec(argflag):
    """
    Load the list of calspec standards

    Parameters
    ----------

    Returns
    -------
    calspec_path : str
      Path from pypitdir to calspec standard star files
    calspec_stds : Table
      astropy Table of the calspec standard stars (file, Name, RA, DEC)
    """
    # Read
    calspec_path = '/data/standards/calspec/'
    calspec_file = argflag['run']['pypitdir'] + calspec_path + 'calspec_info.txt'
    calspec_stds = Table.read(calspec_file, comment='#', format='ascii')
    # Return
    return calspec_path, calspec_stds


def load_extinction_data(slf, toler=1.*u.deg):
    """
    Find the best extinction file to use, based on longitude and latitude
    Loads it and returns a Table

    Parameters
    ----------
    slf : class
      Includes mosaic lon/lat
    toler : Angle
      Tolerance for matching detector to site (1 deg)

    Returns
    -------
    ext_file : Table
      astropy Table containing the 'wavelength', 'extinct' data for AM=1.
    """
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
        msgs.warn("No file found for extinction corrections.  Applying none") 
        msgs.warn("You should generate a site-specific file")
        return None
    # Read
    extinct = Table.read(extinct_path+extinct_file,comment='#',format='ascii', names=('iwave','mag_ext'))
    wave = Column(np.array(extinct['iwave'])*u.AA, name='wave')
    extinct.add_column(wave)
    # Return
    return extinct[['wave','mag_ext']]


def load_standard_file(slf, std_dict):
    """
    Load standard star data

    Parameters
    ----------
    std_dict : dict
      Info on standard star indcluding filename in 'file'
      May be compressed

    Returns
    -------
    std_wave : Quantity array
      Wavelengths of standard star array
    std_flux : Quantity array
      Flux of standard star
    """
    fil = glob.glob(slf._argflag['run']['pypitdir']+
            std_dict['file']+'*')
    if len(fil) == 0:
        msgs.error("No standard star file: {:s}".format(fil))
    else:
        fil = fil[0]
        msgs.info("Loading standard star file: {:s}".format(fil))
        msgs.info("Fluxes are flambda, normalized to 1e-17")

    if std_dict['fmt'] == 1:
        std_spec = fits.open(fil)[1].data
        # Load
        std_dict['wave'] = std_spec['WAVELENGTH']*u.AA 
        std_dict['flux'] = 1e17*std_spec['FLUX']*u.erg/u.s/u.cm**2/u.AA
    else:
        msgs.error("Bad Standard Star Format")
    return

def generate_sensfunc(slf, sc, BALM_MASK_WID=5., nresln=20):
    """
    Generate sensitivity function from current standard star
    Currently, we are using a bspline generated by bspline_magfit

    Parameters
    ----------
    sc : int
      index for standard  (may not be necessary)
    BALM_MASK_WID : float
      Mask parameter for Balmer absorption.  A region equal to 
      BALM_MASK_WID*resln is masked wher resln is the estimate 
      for the spectral resolution.
    nresln : int  
      Number of resolution elements for break-point placement

    Returns
    -------
    sens_dict : dict
      sensitivity function described by a dict
    """
    # Find brightest object in the exposure
    medfx = []
    for spobj in slf._specobjs:
        medfx.append(np.median(spobj.boxcar['counts']))
    std_obj = slf._specobjs[np.argmax(np.array(medfx))]
    wave = std_obj.boxcar['wave']
    # Apply Extinction
    extinct = load_extinction_data(slf)
    ext_corr = extinction_correction(wave,
        slf._fitsdict['airmass'][slf._scidx], extinct)
    flux_corr = std_obj.boxcar['counts']*ext_corr
    var_corr = std_obj.boxcar['var']*ext_corr**2
    # Convert to electrons / s
    flux_corr /= slf._fitsdict['exptime'][slf._scidx] 
    var_corr /= slf._fitsdict['exptime'][slf._scidx]**2

    # Grab closest standard within a tolerance
    std_dict = find_standard_file(slf._argflag, (slf._fitsdict['ra'][slf._scidx],slf._fitsdict['dec'][slf._scidx]))
    # Load standard
    load_standard_file(slf, std_dict)
    # Interpolate onto observed wavelengths
    std_xspec = XSpectrum1D.from_tuple((std_dict['wave'],std_dict['flux']))
    xspec = std_xspec.rebin(wave) # Conserves flambda
    #flux_interp = scipy.interpolate.interp1d(std_dict['wave'],
    #    std_dict['flux'], bounds_error=False, fill_value=0.)
    #flux_true = flux_interp(wave.to('AA').value)
    flux_true = xspec.flux.value
    if np.min(flux_true) == 0.:
        msgs.warn('Your spectrum extends beyond calibrated standard star.')

    # Mask (True = good pixels)
    msk = np.ones_like(flux_true).astype(bool)
    # Mask bad pixels
    msk[var_corr <= 0.] = False

    # Mask edges 
    msk[flux_true <= 0.] = False
    msk[:10] = False
    msk[-10:] = False
    msgs.info("Masking edges")

    # Mask Balmer
    # Compute an effective resolution for the standard. This could be improved
    # to setup an array of breakpoints based on the resolution. At the 
    # moment we are using only one number
    msgs.warn("Should pull resolution from arc line analysis")
    std_res = 2.0*np.median(np.abs(wave - np.roll(wave, 1)))
    resln = std_res

    msgs.info("Masking Balmer")
    lines_balm = np.array([3836.4, 3969.6, 3890.1, 4102.8, 4102.8, 4341.6, 4862.7, 5407.0, 6564.6, 8224.8, 8239.2])*u.AA
    for line_balm in lines_balm: 
      ibalm = np.abs(wave-line_balm) <= BALM_MASK_WID*resln
      msk[ibalm] = False

    #; Mask telluric absorption
    msgs.info("Masking Telluric")
    tell = np.any([((wave >= 7580.0*u.AA) & (wave <= 7750.0*u.AA)), ((wave >= 7160.0*u.AA) & (wave <= 7340.0*u.AA)),((wave >= 6860.0*u.AA)  & (wave <= 6930.0*u.AA))],axis=0)
    msk[tell] = False

    # Fit in magntiudes
    mag_tck = bspline_magfit(wave.value, flux_corr, var_corr, flux_true, bkspace=resln.value*nresln)
    sens_dict = dict(c=mag_tck, func='bspline',min=None,max=None, std=std_dict)
    # Add in wavemin,wavemax
    sens_dict['wave_min'] = np.min(wave) 
    sens_dict['wave_max'] = np.max(wave) 
    return sens_dict
