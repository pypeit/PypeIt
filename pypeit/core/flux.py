""" Module for fluxing routines
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import glob

import numpy as np
import scipy

from pkg_resources import resource_filename

from astropy import units
from astropy import constants
from astropy import coordinates
from astropy.table import Table, Column
from astropy.io import ascii
from astropy.io import fits

try:
    from linetools.spectra.xspectrum1d import XSpectrum1D
except ImportError:
    pass

from pypeit import msgs
from pypeit import utils
from pypeit import debugger


def apply_sensfunc(spec_obj, sensfunc, airmass, exptime, extinction_data, MAX_EXTRAP=0.05):
    """ Apply the sensitivity function to the data
    We also correct for extinction.

    Parameters
    ----------
    MAX_EXTRAP : float, optional [0.05]
      Fractional amount to extrapolate sensitivity function
    """
    # Load extinction data
    #extinct = load_extinction_data(settings_spec)
    #airmass = fitsdict['airmass'][scidx]

    # Loop on extraction modes
    for extract_type in ['boxcar', 'optimal']:
        extract = getattr(spec_obj, extract_type)
        if len(extract) == 0:
            continue
        msgs.info("Fluxing {:s} extraction for:".format(extract_type) + msgs.newline() +
                  "{}".format(spec_obj))
        wave = extract['wave']  # for convenience
        scale = np.zeros(wave.size)
        # Allow for some extrapolation
        dwv = sensfunc['wave_max']-sensfunc['wave_min']
        inds = ((wave >= sensfunc['wave_min']-dwv*MAX_EXTRAP)
            & (wave <= sensfunc['wave_max']+dwv*MAX_EXTRAP))
        mag_func = utils.func_val(sensfunc['c'], wave[inds],
                                    sensfunc['func'])
        sens = 10.0**(0.4*mag_func)
        # Extinction
        ext_corr = extinction_correction(wave[inds], airmass, extinction_data)
        scale[inds] = sens*ext_corr
        # Fill
        extract['flam'] = extract['counts']*scale/exptime
        extract['flam_var'] = (extract['var']*(scale/exptime)**2)

'''
def apply_sensfunc(slf, det, scidx, fitsdict, MAX_EXTRAP=0.05, standard=False):
    """ Apply the sensitivity function to the data
    We also correct for extinction.

    Parameters
    ----------
    MAX_EXTRAP : float, optional [0.05]
      Fractional amount to extrapolate sensitivity function
    """
    # Load extinction data
    extinct = load_extinction_data()
    airmass = fitsdict['airmass'][scidx]
    # Allow application to standard
    if standard:
        specobjs = slf._msstd[det-1]['spobjs']
    else:
        specobjs = slf._specobjs[det-1]
    if specobjs is None:
        msgs.warn("No objects extracted.  Nothing to apply flux standard to.")
        return
    # Loop on slits
    for sl in range(len(specobjs)):
        # Loop on objects
        for spobj in specobjs[sl]:
            # Loop on extraction modes
            for extract_type in ['boxcar', 'optimal']:
                extract = getattr(spobj, extract_type)
                if len(extract) == 0:
                    continue
#                try:
#                    extract = getattr(spobj, extract_type)
#                except AttributeError:
#                    continue
                msgs.info("Fluxing {:s} extraction for:".format(extract_type) + msgs.newline() +
                          "{}".format(spobj))
                wave = extract['wave']  # for convenience
                scale = np.zeros(wave.size)
                # Allow for some extrapolation
                dwv = slf._sensfunc['wave_max']-slf._sensfunc['wave_min']
                inds = ((wave >= slf._sensfunc['wave_min']-dwv*MAX_EXTRAP)
                        & (wave <= slf._sensfunc['wave_max']+dwv*MAX_EXTRAP))
                mag_func = utils.func_val(slf._sensfunc['c'], wave[inds],
                                            slf._sensfunc['func'])
                sens = 10.0**(0.4*mag_func)
                # Extinction
                ext_corr = extinction_correction(wave[inds], airmass, extinct)
                scale[inds] = sens*ext_corr
                # Fill
                extract['flam'] = extract['counts']*scale/fitsdict['exptime'][scidx]
                extract['flam_var'] = (extract['var']*
                                       (scale/fitsdict['exptime'][scidx])**2)
'''

def bspline_magfit(
        wave,
        flux,
        var,
        flux_std,
        maxiter=10,
        upper=2,
        lower=2,
        kwargs_bspline={},
        kwargs_reject={}):
    """
    Perform a bspline fit to the flux ratio of standard to
    observed counts. Used to generate a sensitivity function.

    Parameters
    ----------
    wave : ndarray
      wavelength as observed
    flux : ndarray
      counts/s as observed
    var : ndarray
      variance
    flux_std : Quantity array
      standard star true flux (erg/s/cm^2/A)
    maxiter : integer
      maximum number of iterations for bspline_iterfit
    upper : integer
      number of sigma for rejection in bspline_iterfit
    lower : integer
      number of sigma for rejection in bspline_iterfit
    kwargs_bspline : dict, optional
      keywords for bspline_iterfit
    kwargs_reject : dict, optional
      keywords for bspline_iterfit

    Returns
    -------
    bset_log1 : 
    """

    # Create copy of the arrays to avoid modification
    wave_obs = wave.copy()
    flux_obs = flux.copy()
    var_obs = var.copy()

    from pypit.core.pydl import bspline
    from pypit.core.pydl import iterfit as bspline_iterfit

    """
    # OLD
    from pydl.pydlutils.bspline import bspline
    from pydl.pydlutils.bspline import iterfit as bspline_iterfit
    """

    # preparing arrays to run in bspline_iterfit

    """
    EPF: this line somehow was not working 
     invvar = (var_obs > 0.) / (np.max(var_obs, 0))
    changed to the less elegant, but effective:
     invvar = 1/var_obs
    """

    from pypit.arutils import calc_ivar
    invvar = calc_ivar(var_obs)
    if (np.all(~np.isfinite(invvar))):
        msgs.warn("NaN are present in the inverse variance")

    """
    # maskregions
    invvar[np.where(var_obs ==  0.0)] = 0.0
    invvar[np.where(var_obs == -1.0)] = 0.0
    """

    # Removing 10sigma outliners
    nx = wave_obs.size
    pos_error = 1. / np.sqrt(np.maximum(invvar, 0.) + (invvar == 0))
    pos_mask = (flux_obs > pos_error / 10.0) & (invvar > 0) & (flux_std > 0.0)

    fluxlog = 2.5 * np.log10(np.maximum(flux_obs, pos_error / 10))
    
    """
    Ema: Old logivar
    logivar = invvar * flux_obs ** 2 * pos_mask * 1.08574
    """
    
    """
    1.08574 converts Log_10 into ln
    """
    logivar = invvar * np.power(flux_obs,2.) * pos_mask * np.power(1.08574,-2.)



    # Exclude extreme values of magfunc
    flux_stdlog = 2.5 * np.log10(np.maximum(flux_std, 1.0e-20))
    magfunc = flux_stdlog - fluxlog
    magfunc = np.minimum(magfunc, 25.)
    sensfunc = 10.0 ** (0.4 * magfunc) * pos_mask

    msgs.info("Initialize bspline for flux calibration")

    init_bspline = bspline(wave_obs, bkspace=kwargs_bspline['bkspace'])
    fullbkpt = init_bspline.breakpoints
    # remove masked regions
    msk_obs = np.ones_like(wave_obs).astype(bool)
    msk_obs[var_obs <= 0.] = False
    import scipy.interpolate as interpolate
    msk_bkpt = interpolate.interp1d(wave_obs, msk_obs, kind='nearest', fill_value='extrapolate')(fullbkpt)

    init_breakpoints = fullbkpt[msk_bkpt == 1.]

    msgs.info("Bspline fit: step 1")

    # Check for magfunc
    
    import matplotlib.pyplot as plt
    plt.figure(1)
    plt.ylim(np.min(logivar),np.max(logivar))
    plt.plot(wave_obs, logivar, label='logivar')
    plt.legend()
    plt.xlabel('Wavelength [ang]')
    plt.show()
    plt.close()



    #  First round of the fit:
    bset1, bmask = bspline_iterfit(wave_obs, magfunc, invvar=logivar,
                                   upper=upper, lower=lower,
                                   maxiter=maxiter,
                                   fullbkpt=init_breakpoints,
                                   kwargs_bspline=kwargs_bspline,
                                   kwargs_reject=kwargs_reject)

    """
    # OLD
    bset1, bmask = bspline_iterfit(wave_obs, magfunc, invvar=logivar,
                                   fullbkpt=init_breakpoints,
                                   maxiter=maxiter,
                                   upper=upper, lower=lower, kwargs_bspline=kwargs_bspline,
                                   kwargs_reject=kwargs_reject)
    """

    # Calculate residuals
    logfit1, _ = bset1.value(wave_obs)
    modelfit1 = 10.0 ** (0.4 * logfit1)
    residual = sensfunc / (modelfit1 + (modelfit1 == 0)) - 1.
    new_mask = pos_mask & (sensfunc > 0)

    residual_ivar = (modelfit1 * flux_obs / (sensfunc + (sensfunc == 0.0))) ** 2 * invvar
    residual_ivar = residual_ivar * new_mask

    """
    # Check for magfunc
    import matplotlib.pyplot as plt
    plt.figure(1)
    plt.ylim(np.min(magfunc),np.max(magfunc))
    plt.plot(wave_obs, magfunc, label='magfunc')
    plt.plot(wave_obs, logfit1, label='logfit1')
    plt.legend()
    plt.xlabel('Wavelength [ang]')
    plt.show()
    plt.close()

    # Check for calibration
    import matplotlib.pyplot as plt
    plt.figure(1)
    plt.plot(wave_obs, sensfunc*flux_obs, label='scaled')
    plt.plot(wave_obs, flux_std, label='model')
    plt.legend()
    plt.xlabel('Wavelength [ang]')
    plt.show()
    plt.close()
    """

    msgs.info("Bspline fit: step 2")

    #  Now do one more fit to the ratio of data/model - 1.
    bset_residual, bmask2 = bspline_iterfit(wave_obs, residual, invvar=residual_ivar,
                                            fullbkpt=bset1.breakpoints, maxiter=maxiter, upper=upper,
                                            lower=lower, kwargs_bspline=kwargs_bspline, kwargs_reject=kwargs_reject)
    """
    OLD
    bset_residual, bmask2 = bspline_iterfit(wave_obs, residual, invvar=residual_ivar,
                                            fullbkpt=bset1.breakpoints,
                                            maxiter=maxiter,
                                            upper=upper, lower=lower, kwargs_bspline=kwargs_bspline,
                                            kwargs_reject=kwargs_reject)
    """

    # Create sensitivity function
    bset_log1 = bset1
    bset_log1.coeff = bset_log1.coeff + bset_residual.coeff
    newlogfit, _ = bset_log1.value(wave_obs)
    sensfit = np.power(10.0, 0.4 * newlogfit)

    ## print(bset_log1)

    bspline_dict={}

    print(bspline_dict)

    # Write the sens_dict to a json file
    msgs.info("Writing bspline_dict into .json file")
    with open('bspline_dict.json', 'w') as fp:
        json.dump(bspline_dict, fp, sort_keys=True, indent=4)


    # bspline_func = bspline_fromdict(wave_obs, from_dict=bspline_dict)


    """
    # Check for calibration
    import matplotlib.pyplot as plt
    plt.figure(1)
    plt.plot(wave_obs, sensfit*flux_obs, label='scaled')
    plt.plot(wave_obs, flux_std, label='model')
    plt.legend()
    plt.xlabel('Wavelength [ang]')
    plt.show()
    plt.close()
    """

    # Check quality of the fit
    absdev = np.median(np.abs(sensfit / modelfit1 - 1))
    msgs.info('Difference between fits is {:g}'.format(absdev))

    """
    # Check for residual of the fit
    import matplotlib.pyplot as plt
    plt.figure(1)
    plt.plot(wave_obs, sensfit / modelfit1 - 1, label='residual')
    plt.legend()
    plt.xlabel('Wavelength [ang]')
    plt.show()
    plt.close()
    """

    # QA
    msgs.work("Add QA for sensitivity function")
    """
    bspline_magfit_new_qa(wave_obs, magfunc, logfit1,
                          newlogfit, bset1.breakpoints,
                          outfile=None, title=None)
    """

    return bset_log1


'''
# Old version
def bspline_magfit(wave, flux, var, flux_std, bspline_par=None):
    """
    Perform a bspline fit to the flux ratio of standard to
    observed counts.  Used to generate a sensitivity function.

    Parameters
    ----------
    wave : ndarray
    flux : ndarray
      counts/s as observed
    var : ndarray
      variance
    flux_std : Quantity array
      standard star true flux (erg/s/cm^2/A)
    bspline_par : dict, optional
      keywords for robust_polyfit

    Returns
    -------
    """
    invvar = (var > 0.)/(np.max(var,0))
    invvar[np.where(var==0.0)] = 0.0
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

    """ Ported from LowRedux but this is a bad idea (and wasn't used there anyhow)
    # Interpolate over masked pixels
    if not nointerp:
        bad = logivar <= 0.
        if np.sum(bad) > 0:
            f = scipy.interpolate.InterpolatedUnivariateSpline(wave[~bad], magfunc[~bad], k=2)
            magfunc[bad] = f(wave[bad])
            fi = scipy.interpolate.InterpolatedUnivariateSpline(wave[~bad], logivar[~bad], k=2)
            logivar[bad] = fi(wave[bad])
    """

    #  First iteration
    mask, tck = utils.robust_polyfit(wave, magfunc, 3, function='bspline', weights=np.sqrt(logivar),
                                       bspline_par=bspline_par)
    logfit1 = utils.func_val(tck,wave,'bspline')
    modelfit1 = 10.0**(0.4*(logfit1))

    residual = sensfunc/(modelfit1 + (modelfit1 == 0)) - 1.
    new_mask = pos_mask & (sensfunc > 0)
    residual_ivar = (modelfit1*flux/(sensfunc + (sensfunc == 0.0)))**2*invvar
    residual_ivar = residual_ivar*new_mask

    """ Ported from LowRedux but this is a bad idea (and wasn't used there anyhow)
    # Interpolate over masked pixels
    if not nointerp:
        if np.sum(bad) > 0:
            f = scipy.interpolate.InterpolatedUnivariateSpline(wave[~bad], residual[~bad], k=2)
            residual[bad] = f(wave[bad])
            fi = scipy.interpolate.InterpolatedUnivariateSpline(wave[~bad], residual_ivar[~bad], k=2)
            residual_ivar[bad] = fi(wave[bad])
    """

    #  Now do one more fit to the ratio of data/model - 1.
    # Fuss with the knots first ()
    inner_knots = utils.bspline_inner_knots(tck[0])
    #
    if bspline_par is None:
        bspline_par = {}
    bspline_par['knots'] = inner_knots # This over-rides everyn
    mask, tck_residual = utils.robust_polyfit(wave, residual, 3, function='bspline',
                                                weights=np.sqrt(residual_ivar), bspline_par=bspline_par) #**kwargs)
    if tck_residual[1].size != tck[1].size:
        msgs.error('Problem with bspline knots in bspline_magfit')
    #bset_residual = bspline_iterfit(wave, residual, weights=np.sqrt(residual_ivar), knots = tck[0])

    tck_log1 = list(tck)
    tck_log1[1] = tck[1] + tck_residual[1]

    sensfit = 10.0**(0.4*(utils.func_val(tck_log1,wave, 'bspline')))

    absdev = np.median(np.abs(sensfit/modelfit1-1))
    msgs.info('Difference between fits is {:g}'.format(absdev))

    # QA
    msgs.work("Add QA for sensitivity function")

    return tck_log1
'''

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


def find_standard_file(radec, toler=20.*units.arcmin, check=False):
    """
    Find a match for the input file to one of the archived
    standard star files (hopefully).  Priority is by order of search.

    Parameters
    ----------
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
    obj_coord = coordinates.SkyCoord(radec[0], radec[1], unit=(units.hourangle, units.deg))
    # Loop on standard sets
    closest = dict(sep=999*units.deg)
    for qq,sset in enumerate(std_sets):
        # Stars
        path, star_tbl = sset()
        star_coords = coordinates.SkyCoord(star_tbl['RA_2000'], star_tbl['DEC_2000'],
                                           unit=(units.hourangle, units.deg))
        # Match
        idx, d2d, d3d = coordinates.match_coordinates_sky(obj_coord, star_coords, nthneighbor=1)
        if d2d < toler:
            if check:
                return True
            else:
                # Generate a dict
                std_dict = dict(file=path+star_tbl[int(idx)]['File'], name=star_tbl[int(idx)]['Name'], fmt=std_file_fmt[qq], ra=star_tbl[int(idx)]['RA_2000'], dec=star_tbl[int(idx)]['DEC_2000'])
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


def load_calspec():
    """
    Load the list of calspec standards

    Parameters
    ----------

    Returns
    -------
    calspec_path : str
      Path from pypeitdir to calspec standard star files
    calspec_stds : Table
      astropy Table of the calspec standard stars (file, Name, RA, DEC)
    """
    # Read
    calspec_path = '/data/standards/calspec/'
    calspec_file = resource_filename('pypeit', calspec_path+'calspec_info.txt')
    calspec_stds = Table.read(calspec_file, comment='#', format='ascii')
    # Return
    return calspec_path, calspec_stds

def load_extinction_data(longitude, latitude, toler=5.*units.deg):
    """
    Find the best extinction file to use, based on longitude and latitude
    Loads it and returns a Table

    Parameters
    ----------
    toler : Angle, optional
      Tolerance for matching detector to site (5 deg)

    Returns
    -------
    ext_file : Table
      astropy Table containing the 'wavelength', 'extinct' data for AM=1.
    """
    # Mosaic coord
    mosaic_coord = coordinates.SkyCoord(longitude, latitude, frame='gcrs', unit=units.deg)
    # Read list
    extinct_path = resource_filename('pypeit', '/data/extinction/')
    extinct_summ = extinct_path+'README'
    extinct_files = Table.read(extinct_summ,comment='#',format='ascii')
    # Coords
    ext_coord = coordinates.SkyCoord(extinct_files['Lon'], extinct_files['Lat'], frame='gcrs',
                                     unit=units.deg)
    # Match
    idx, d2d, d3d = coordinates.match_coordinates_sky(mosaic_coord, ext_coord, nthneighbor=1)
    if d2d < toler:
        extinct_file = extinct_files[int(idx)]['File']
        msgs.info("Using {:s} for extinction corrections.".format(extinct_file))
    else:
        msgs.warn("No file found for extinction corrections.  Applying none")
        msgs.warn("You should generate a site-specific file")
        return None
    # Read
    extinct = Table.read(extinct_path+extinct_file,comment='#',format='ascii',
                         names=('iwave','mag_ext'))
    wave = Column(np.array(extinct['iwave'])*units.AA, name='wave')
    extinct.add_column(wave)
    # Return
    return extinct[['wave','mag_ext']]


def load_standard_file(std_dict):
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
    root = resource_filename('pypeit', std_dict['file']+'*')
    fil = glob.glob(root)
    if len(fil) == 0:
        msgs.error("No standard star file: {:s}".format(fil))
    else:
        fil = fil[0]
        msgs.info("Loading standard star file: {:s}".format(fil))
        msgs.info("Fluxes are flambda, normalized to 1e-17")

    if std_dict['fmt'] == 1:
        std_spec = fits.open(fil)[1].data
        # Load
        std_dict['wave'] = std_spec['WAVELENGTH']*units.AA
        std_dict['flux'] = 1e17*std_spec['FLUX']*units.erg/units.s/units.cm**2/units.AA
    else:
        msgs.error("Bad Standard Star Format")
    return


def find_standard(specobjs):
    """
    Take the median boxcar and then the max object as the standard

    Parameters
    ----------
    specobjs : list

    Returns
    -------

    """
    # Repackage as necessary (some backwards compatability)
    all_specobj = utils.unravel_specobjs(specobjs)
    # Do it
    medfx = []
    for indx, spobj in enumerate(all_specobj):
        if spobj is None:
            medfx.append(0.)
        else:
            medfx.append(np.median(spobj.boxcar['counts']))
    try:
        mxix = np.argmax(np.array(medfx))
    except:
        debugger.set_trace()
    msgs.info("Putative standard star {} has a median boxcar count of {}".format(all_specobj[mxix],
                                                                                 np.max(medfx)))
    # Return
    return mxix


def generate_sensfunc(
        wave,
        flux,
        var,
        airmass,
        exptime,
        spectrograph,
        telluric=False,
        star_type=None,
        star_mag=None,
        RA=None,
        DEC=None,
        BALM_MASK_WID=5.,
        nresln=None):
    """ Function to generate the sensitivity function.
    This can work in different regimes:
    - If telluric=False and RA=None and Dec=None
      the code creates a sintetic standard star spectrum using the Kurucz models,
      and from this it generates a sens func using nresln=20.0 and masking out
      telluric regions.
    - If telluric=False and RA and Dec are assigned
      the standard star spectrum is extracted from the archive, and a sens func 
      is generated using nresln=20.0 and masking out telluric regions.
    - If telluric=True
      the code creates a sintetic standard star spectrum using the Kurucz models,
      the sens func is created setting nresln=1.5 it contains the correction for
      telluric lines.

    Parameters:
    ----------
    wave : array
      Wavelength of the star with units
    flux : array
      Flux (in counts) of the star
    var : array
      Variance of the star
    airmass : float
      Airmass
    exptime : float
      Exposure time in seconds
    spectrograph : dict
      Instrument specific dict
      Used for extinction correction
    telluric : bool
      if True performs a telluric correction
    star_type : str
      Spectral type of the telluric star (used if telluric=True)
    star_mag : float
      Apparent magnitude of telluric star (used if telluric=True)
    RA : float
      deg, RA of the telluric star
      if assigned, the standard star spectrum will be extracted from
      the archive
    DEC : float
      deg, DEC of the telluric star
      if assigned, the standard star spectrum will be extracted from
      the archive
    BALM_MASK_WID : float
      Mask parameter for Balmer absorption. A region equal to
      BALM_MASK_WID*resln is masked wher resln is the estimate
      for the spectral resolution.
    nresln : float
      Number of resolution elements for break-point placement.
      If assigned, overwrites the settings imposed by the code.

    Returns:
    -------
    sens_dict : dict
      sensitivity function described by a dict
    """

    # Create copy of the arrays to avoid modification
    wave_star = wave.copy()
    flux_star = flux.copy()
    var_star = var.copy()

    # ToDo
    # This should be changed. At the moment the extinction correction procedure
    # requires the spectra to be in the optical. For the NIR is probably enough
    # to extend the tables to longer wavelength setting the extinction to 0.0mag.
    msgs.warn("Extinction correction applyed only if the spectra covers <10000Ang.")
    # Apply Extinction if optical bands
    if np.max(wave_star) < 10000. * units.AA:
        msgs.info("Applying extinction correction")
        extinct = flux.load_extinction_data(spectrograph.telescope['longitude'],
                                            spectrograph.telescope['latitude'])
        ext_corr = extinction_correction(wave_star, airmass, extinct)
        # Correct for extinction and convert to electrons / s
        flux_corr = flux_star * ext_corr / exptime
        var_corr = var_star * ext_corr ** 2 / exptime ** 2
    else:
        msgs.info("Extinction correction not applied")
        # Convert to electrons / s
        flux_corr = flux_star / exptime
        var_corr = var_star / exptime ** 2

    # Create star model
    if (RA is not None) and (DEC is not None):
        # Pull star spectral model from archive
        msgs.info("Get standard model")
        # Grab closest standard within a tolerance
        std_dict = find_standard_file((RA, DEC))
        # Load standard
        load_standard_file(std_dict)
    else:
        # Create star spectral model
        msgs.info("Creating standard model")
        # Create star model
        star_loglam, star_flux, std_dict = telluric_sed(star_mag, star_type)
        star_lam = 10 ** star_loglam
        # Generate a dict matching the output of find_standard_file
        std_dict = dict(file='KuruczTelluricModel', name=star_type, fmt=1,
                        ra='00:00:00.0', dec='00:00:00.0')
        std_dict['wave'] = star_lam*units.AA
        std_dict['flux'] = 1e17*star_flux*units.erg/units.s/units.cm**2/units.AA

    # Interpolate onto observed wavelengths
    std_xspec = XSpectrum1D.from_tuple((std_dict['wave'], std_dict['flux']))
    xspec = std_xspec.rebin(wave_star)  # Conserves flambda
    flux_true = xspec.flux.value
    if np.min(flux_true) == 0.:
        msgs.warn('Your spectrum extends beyond calibrated standard star.')

    # Set nresln
    if nresln == None:
        if telluric:
            nresln = 1.5
            msgs.info("Set nresln to 1.5")
        else:
            nresln = 20.0
            msgs.info("Set nresln to 20.0")

    # ToDo
    # Compute an effective resolution for the standard. This could be improved
    # to setup an array of breakpoints based on the resolution. At the
    # moment we are using only one number
    msgs.work("Should pull resolution from arc line analysis")
    msgs.work("At the moment the resolution is taken as 4 x PixelScale")
    msgs.work("This needs to be changed!")
    std_pix = np.median(np.abs(wave_star - np.roll(wave_star, 1)))
    std_res = 4.0 * std_pix
    resln = std_res
    if (nresln * std_res) < std_pix:
        msgs.warn("Bspline breakpoints spacing shoud be larger than 1pixel")
        msgs.warn("Changing input nresln to fix this")
        nresln = std_res / std_pix

    # Mask bad pixels, edges, and Balmer, Paschen, Brackett, and Pfund lines
    # Mask (True = good pixels)
    msgs.info("Masking spectral regions:")
    msk_corr = np.ones_like(flux_corr).astype(bool)

    # Mask bad pixels
    msgs.info(" Masking bad pixels")
    msk_corr[var_corr <= 0.] = False
    msk_corr[flux_corr <= 0.] = False

    # Mask edges
    msgs.info(" Masking edges")
    msk_corr[:10] = False
    msk_corr[-10:] = False

    # Mask Balmer
    msgs.info(" Masking Balmer")
    lines_balm = np.array([3836.4, 3969.6, 3890.1, 4102.8, 4102.8, 4341.6, 4862.7, 5407.0,
                           6564.6, 8224.8, 8239.2]) * units.AA
    for line_balm in lines_balm:
        ibalm = np.abs(wave_star - line_balm) <= BALM_MASK_WID * resln
        msk_corr[ibalm] = False

    # Mask Paschen
    msgs.info(" Masking Paschen")
    # air wavelengths from:
    # https://www.subarutelescope.org/Science/Resources/lines/hi.html
    lines_pasc = np.array([8203.6, 9229.0, 9546.0, 10049.4, 10938.1,
                           12818.1, 18751.0]) * units.AA
    for line_pasc in lines_pasc:
        ipasc = np.abs(wave_star - line_pasc) <= BALM_MASK_WID * resln
        msk_corr[ipasc] = False

    # Mask Brackett
    msgs.info(" Masking Brackett")
    # air wavelengths from:
    # https://www.subarutelescope.org/Science/Resources/lines/hi.html
    lines_brac = np.array([14584.0, 18174.0, 19446.0, 21655.0,
                           26252.0, 40512.0]) * units.AA
    for line_brac in lines_brac:
        ibrac = np.abs(wave_star - line_brac) <= BALM_MASK_WID * resln
        msk_corr[ibrac] = False

    # Mask Pfund
    msgs.info(" Masking Pfund")
    # air wavelengths from:
    # https://www.subarutelescope.org/Science/Resources/lines/hi.html
    lines_pfund = np.array([22788.0, 32961.0, 37395.0, 46525.0,
                            74578.0]) * units.AA
    for line_pfund in lines_pfund:
        ipfund = np.abs(wave_star - line_pfund) <= BALM_MASK_WID * resln
        msk_corr[ipfund] = False

    # Mask Atm. cutoff
    msgs.info(" Masking Below the atmospheric cutoff")
    atms_cutoff = wave_star <= 3000.0 * units.AA
    msk_corr[atms_cutoff] = False

    if telluric == False:
        # Mask telluric absorption
        msgs.info("Masking Telluric")
        tell = np.any([((wave >= 7580.00*units.AA) & (wave <= 7750.00*units.AA)),
                       ((wave >= 7160.00*units.AA) & (wave <= 7340.00*units.AA)),
                       ((wave >= 6860.00*units.AA) & (wave <= 6930.00*units.AA)),
                       ((wave >= 9310.00*units.AA) & (wave <= 9665.00*units.AA)),
                       ((wave >= 11120.0*units.AA) & (wave <= 11615.0*units.AA)),
                       ((wave >= 12610.0*units.AA) & (wave <= 12720.0*units.AA)),
                       ((wave >= 13160.0*units.AA) & (wave <= 15065.0*units.AA)),
                       ((wave >= 15700.0*units.AA) & (wave <= 15770.0*units.AA)),
                       ((wave >= 16000.0*units.AA) & (wave <= 16100.0*units.AA)),
                       ((wave >= 16420.0*units.AA) & (wave <= 16580.0*units.AA)),
                       ((wave >= 17310.0*units.AA) & (wave <= 20775.0*units.AA)),
                       (wave >= 22680.0*units.AA)], axis=0)
        msk_corr[tell] = False

    # Apply mask
    var_corr[msk_corr == False] = -1.

    # Fit in magnitudes
    kwargs_bspline = {'bkspace': resln.value * nresln}
    kwargs_reject = {'maxrej': 5}
    mag_set = bspline_magfit(wave_star.value, flux_corr, var_corr,
                             flux_true, kwargs_bspline=kwargs_bspline,
                             kwargs_reject=kwargs_reject)

    # Creating the dict
    msgs.work("Is min, max and wave_min, wave_max a duplicate?")
    sens_dict = dict(bspline=mag_set, func='bspline', min=None, max=None, std=std_dict)

    # Add in wavemin,wavemax
    sens_dict['wave_min'] = np.min(wave_star)
    sens_dict['wave_max'] = np.max(wave_star)


    """
    # Write the sens_dict to a json file
    msgs.info("Writing sens_dict into .json file")
    with open('sens_dict.json', 'w') as fp:
        json.dump(sens_dict, fp, sort_keys=True, indent=4)
    """

    return sens_dict




'''
def generate_sensfunc(std_obj, RA, DEC, exptime, extinction, BALM_MASK_WID=5., nresln=20):
    """
    Generate sensitivity function from current standard star
    Currently, we are using a bspline generated by bspline_magfit

    Parameters
    ----------
    std_obj : list
      List of spectra
    RA : float
      Deg
    DEC : float
      Deg
    airmass : float
    exptime : float
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
    wave = std_obj.boxcar['wave']
    # Apply Extinction
#    extinct = load_extinction_data(settings_spect) # Observatory specific
#    ext_corr = extinction_correction(wave, airmass, extinction_data)
    flux_corr = std_obj.boxcar['counts'] * extinction
    var_corr = std_obj.boxcar['var'] * np.square(extinction)
    # Convert to electrons / s
    flux_corr /= exptime
    var_corr /= exptime**2

    # Grab closest standard within a tolerance
    std_dict = find_standard_file((RA, DEC))
    # Load standard
    load_standard_file(std_dict)
    # Interpolate onto observed wavelengths
    std_xspec = XSpectrum1D.from_tuple((std_dict['wave'], std_dict['flux']))
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
    lines_balm = np.array([3836.4, 3969.6, 3890.1, 4102.8, 4102.8, 4341.6, 4862.7, 5407.0,
                           6564.6, 8224.8, 8239.2])*units.AA
    for line_balm in lines_balm:
      ibalm = np.abs(wave-line_balm) <= BALM_MASK_WID*resln
      msk[ibalm] = False

    #; Mask telluric absorption
    msgs.info("Masking Telluric")
    tell = np.any([((wave >= 7580.0*units.AA) & (wave <= 7750.0*units.AA)),
                   ((wave >= 7160.0*units.AA) & (wave <= 7340.0*units.AA)),
                   ((wave >= 6860.0*units.AA)  & (wave <= 6930.0*units.AA))],axis=0)
    msk[tell] = False

    # Mask
    msgs.info("Masking Below the atmospheric cutoff")
    atms_cutoff = wave <= 3000.0*units.AA
    msk[atms_cutoff] = False

    # Fit in magnitudes
    var_corr[msk == False] = -1.
    bspline_par = dict(bkspace=resln.value*nresln)
    mag_tck = bspline_magfit(wave.value, flux_corr, var_corr, flux_true, bspline_par=bspline_par) #bkspace=resln.value*nresln)
    sens_dict = dict(c=mag_tck, func='bspline',min=None,max=None, std=std_dict)
    # Add in wavemin,wavemax
    sens_dict['wave_min'] = np.min(wave)
    sens_dict['wave_max'] = np.max(wave)
    return sens_dict
'''


def telluric_params(sptype):
    """Compute physical parameters for a given stellar type.
    This is used by telluric_sed(V, sptype) to create the
    model spectrum.

    Parameters:
    ----------
    sptype: str
      Spectral type of telluric star

    Returns:
    ----------
    tell_param: dict
      Star parameters
    """

    # log(g) of the Sun
    logg_sol = np.log10(6.67259e-8) + np.log10(1.989e33) - 2.0 * np.log10(6.96e10)

    # Load Schmidt-Kaler (1982) table
    sk82_file = resource_filename('pypeit', 'data/standards/kurucz93/schmidt-kaler_table.txt')
    sk82_tab = ascii.read(sk82_file, names=('Sp', 'logTeff', 'Teff', '(B-V)_0', 'M_V', 'B.C.', 'M_bol', 'L/L_sol'))

    # Match input type
    mti = np.where(sptype == sk82_tab['Sp'])[0]
    if len(mti) != 1:
        raise ValueError('Not ready to interpolate yet.')

    # Calculate final quantities
    # relation between radius, temp, and bolometric luminosity
    logR = 0.2 * (42.26 - sk82_tab['M_bol'][mti[0]] - 10.0 * sk82_tab['logTeff'][mti[0]])

    # mass-bolometric luminosity relation 
    # from schimdt-kaler p28 valid for M_bol < 7.5
    logM = 0.46 - 0.10 * sk82_tab['M_bol'][mti[0]]
    logg = logM - 2.0 * logR + logg_sol
    M_V = sk82_tab['M_V'][mti[0]]
    tell_param = dict(logR=logR, logM=logM, logg=logg, M_V=M_V,
                      T=sk82_tab['Teff'][mti[0]])

    # Return
    return tell_param

def telluric_sed(V, sptype):
    """Parse Kurucz SED given T and g
    Also convert absolute/apparent magnitudes

    Parameters:
    ----------
    V: float
      Apparent magnitude of telluric star
    sptype: str
      Spectral type of telluric star

    Returns:
    ----------
    loglam: ndarray
      log wavelengths
    flux: ndarray
      SED f_lambda (cgs units, I think, probably per Ang)
    """

    # Grab Telluric star parameters
    tell_param = telluric_params(sptype)

    # Flux factor (absolute/apparent V mag)
    # Constants
    parsec = constants.pc.cgs  # 3.086e18
    R_sol = constants.R_sun.cgs  # 6.96e10
    # distance modulus
    logd = 0.2 * (V - tell_param['M_V']) + 1.0
    D = parsec * 10. ** logd
    R = R_sol * 10. ** tell_param['logR']
    # factor converts the kurucz surface flux densities to flux observed on Earth
    flux_factor = (R / D.value) ** 2

    # Grab closest T in Kurucz SEDs
    T1 = 3000. + np.arange(28) * 250
    T2 = 10000. + np.arange(6) * 500
    T3 = 13000. + np.arange(22) * 1000
    T4 = 35000. + np.arange(7) * 2500

    Tk = np.concatenate([T1, T2, T3, T4])
    indT = np.argmin(np.abs(Tk - tell_param['T']))

    # Grab closest g in Kurucz SEDs
    loggk = np.arange(11) * 0.5
    indg = np.argmin(np.abs(loggk - tell_param['logg']))

    # Grab Kurucz filename
    std_file = resource_filename('pypeit', '/data/standards/kurucz93/kp00/kp00_{:d}.fits.gz'.format(int(Tk[indT])))
    std = Table.read(std_file)

    # Grab specific spectrum
    loglam = np.array(np.log10(std['WAVELENGTH']))
    gdict = {0: 'g00', 1: 'g05', 2: 'g10', 3: 'g15', 4: 'g20',
             5: 'g25', 6: 'g30', 7: 'g35', 8: 'g40', 9: 'g45', 10: 'g50'}
    flux = std[gdict[indg]]

    # Generate the standard star dict
    std_dict = dict(stellar_type=sptype, Vmag=V)

    # Return
    return loglam, flux.data * flux_factor, std_dict
