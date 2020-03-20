""" Routines related to flexure, air2vac, etc. """
import inspect

import numpy as np

from matplotlib import pyplot as plt
from matplotlib import gridspec


from astropy import units
from astropy.coordinates import solar_system, ICRS
from astropy.coordinates import UnitSphericalRepresentation, CartesianRepresentation
from astropy.time import Time


from pypeit import msgs
from pypeit.core import qa
from pypeit import utils
from IPython import embed


def geomotion_calculate(radec, time, longitude, latitude, elevation, refframe):
    """
    Correct the wavelength calibration solution to the desired reference frame
    """

    # Time
    loc = (longitude * units.deg, latitude * units.deg, elevation * units.m,)
    obstime = Time(time.value, format=time.format, scale='utc', location=loc)
    return geomotion_velocity(obstime, radec, frame=refframe)


def geomotion_correct(specObjs, radec, time, maskslits, longitude, latitude,
                      elevation, refframe):
    """
    Correct the wavelength of every pixel to a barycentric/heliocentric frame.

    Args:
        specObjs (SpecObjs object):
        radec (astropy.coordiantes.SkyCoord):
        time (:obj:`astropy.time.Time`):
        maskslits
        fitstbl : Table/PypeItMetaData
            Containing the properties of every fits file
        longitude (float): deg
        latitude (float): deg
        elevation (float): m
        refframe (str):

    Returns:
        tuple: Two objects are returned:

            - float: The velocity correction that should be applied to
              the wavelength array.
            - float: The relativistic velocity correction that should be
              multiplied by the wavelength array to convert each
              wavelength into the user-specified reference frame.

    """
    # Calculate
    vel = geomotion_calculate(radec, time, longitude, latitude, elevation, refframe)
    vel_corr = np.sqrt((1. + vel/299792.458) / (1. - vel/299792.458))

    gdslits = np.where(np.invert(maskslits))[0]
    # Loop on slits to apply
    for slit in gdslits:
        indx = specObjs.slitorder_indices(slit)
        this_specobjs = specObjs[indx]
        # Loop on objects
        for specobj in this_specobjs:
            if specobj is None:
                continue
            specobj.apply_helio(vel_corr, refframe)
    # Return
    return vel, vel_corr  # Mainly for debugging


def geomotion_velocity(time, skycoord, frame="heliocentric"):
    """ Perform a barycentric/heliocentric velocity correction.

    For the correciton, this routine uses the ephemeris:  astropy.coordinates.solar_system_ephemeris.set
    For more information see `~astropy.coordinates.solar_system_ephemeris`.

    Parameters
    ----------
    time : astropy.time.Time
        The time of observation, including the location.
    skycoord: astropy.coordinates.SkyCoord
        The RA and DEC of the pointing, as a SkyCoord quantity.
    frame : str
        The reference frame that should be used for the calculation.

    Returns
    -------
    vcorr : float
        The velocity correction that should be added to the original velocity.
    """

    # Check that the RA/DEC of the object is ICRS compatible
    if not skycoord.is_transformable_to(ICRS()):
        msgs.error("Cannot transform RA/DEC of object to the ICRS")

    # Calculate ICRS position and velocity of Earth's geocenter
    ep, ev = solar_system.get_body_barycentric_posvel('earth', time)
    # Calculate GCRS position and velocity of observatory
    op, ov = time.location.get_gcrs_posvel(time)
    # ICRS and GCRS are axes-aligned. Can add the velocities
    velocity = ev + ov
    if frame == "heliocentric":
        # ICRS position and velocity of the Sun
        sp, sv = solar_system.get_body_barycentric_posvel('sun', time)
        velocity += sv

    # Get unit ICRS vector in direction of SkyCoord
    sc_cartesian = skycoord.icrs.represent_as(UnitSphericalRepresentation).represent_as(CartesianRepresentation)
    return sc_cartesian.dot(velocity).to(units.km / units.s).value


def airtovac(wave):
    """ Convert air-based wavelengths to vacuum

    Parameters
    ----------
    wave: Quantity array
        Wavelengths 

    Returns
    -------
    wave: Quantity array
        Wavelength array corrected to vacuum wavelengths
    """
    # Convert to AA
    wave = wave.to(units.AA)
    wavelength = wave.value

    # Standard conversion format
    sigma_sq = (1.e4/wavelength)**2. #wavenumber squared
    factor = 1 + (5.792105e-2/(238.0185-sigma_sq)) + (1.67918e-3/(57.362-sigma_sq))
    factor = factor*(wavelength>=2000.) + 1.*(wavelength<2000.) #only modify above 2000A

    # Convert
    wavelength = wavelength*factor
    # Units
    new_wave = wavelength*units.AA
    new_wave.to(wave.unit)

    return new_wave


def vactoair(wave):
    """Convert to air-based wavelengths from vacuum

    Parameters
    ----------
    wave: Quantity array
        Wavelengths 

    Returns
    -------
    wave: Quantity array
        Wavelength array corrected to air

    """
    # Convert to AA
    wave = wave.to(units.AA)
    wavelength = wave.value

    # Standard conversion format
    sigma_sq = (1.e4/wavelength)**2. #wavenumber squared
    factor = 1 + (5.792105e-2/(238.0185-sigma_sq)) + (1.67918e-3/(57.362-sigma_sq))
    factor = factor*(wavelength>=2000.) + 1.*(wavelength<2000.) #only modify above 2000A

    # Convert
    wavelength = wavelength/factor
    new_wave = wavelength*units.AA
    new_wave.to(wave.unit)

    return new_wave

# TODO I don't see why maskslits is needed in these routine, since if the slits are masked in arms, they won't be extracted
#  AND THIS IS WHY THE CODE IS CRASHING
def flexure_qa(specobjs, maskslits, basename, det, flex_list,
               slit_cen=False, out_dir=None):
    """

    Args:
        specobjs:
        maskslits (np.ndarray):
        basename (str):
        det (int):
        flex_list (list):
        slit_cen:
        out_dir:

    """
    plt.rcdefaults()
    plt.rcParams['font.family']= 'times new roman'

    # Grab the named of the method
    method = inspect.stack()[0][3]
    #
    gdslits = np.where(np.invert(maskslits))[0]

    # Loop over slits, and then over objects here
    for slit in gdslits:
        indx = specobjs.slitorder_indices(slit)
        this_specobjs = specobjs[indx]
        this_flex_dict = flex_list[slit]

        # Setup
        if slit_cen:
            nobj = 1
            ncol = 1
        else:
            nobj = np.sum(indx)
            ncol = min(3, nobj)
        #
        if nobj == 0:
            continue
        nrow = nobj // ncol + ((nobj % ncol) > 0)
        # Outfile, one QA file per slit
        outfile = qa.set_qa_filename(basename, method + '_corr', det=det,slit=(slit + 1), out_dir=out_dir)
        plt.figure(figsize=(8, 5.0))
        plt.clf()
        gs = gridspec.GridSpec(nrow, ncol)
        for iobj, specobj in enumerate(this_specobjs):
            if specobj is None or (specobj.BOX_WAVE is None and specobj.OPT_WAVE is None):
                continue
            # Correlation QA
            ax = plt.subplot(gs[iobj//ncol, iobj % ncol])
            # Fit
            fit = this_flex_dict['polyfit'][iobj]
            xval = np.linspace(-10., 10, 100) + this_flex_dict['corr_cen'][iobj] #+ flex_dict['shift'][o]
            #model = (fit[2]*(xval**2.))+(fit[1]*xval)+fit[0]
            model = utils.func_val(fit, xval, 'polynomial')
            mxmod = np.max(model)
            ylim_min = np.min(model/mxmod) if np.isfinite(np.min(model/mxmod)) else 0.0
            ylim = [ylim_min, 1.3]
            ax.plot(xval-this_flex_dict['corr_cen'][iobj], model/mxmod, 'k-')
            # Measurements
            ax.scatter(this_flex_dict['subpix'][iobj]-this_flex_dict['corr_cen'][iobj],
                       this_flex_dict['corr'][iobj]/mxmod, marker='o')
            # Final shift
            ax.plot([this_flex_dict['shift'][iobj]]*2, ylim, 'g:')
            # Label
            if slit_cen:
                ax.text(0.5, 0.25, 'Slit Center', transform=ax.transAxes, size='large', ha='center')
            else:
                ax.text(0.5, 0.25, '{:s}'.format(specobj.NAME), transform=ax.transAxes, size='large', ha='center')
            ax.text(0.5, 0.15, 'flex_shift = {:g}'.format(this_flex_dict['shift'][iobj]),
                    transform=ax.transAxes, size='large', ha='center')#, bbox={'facecolor':'white'})
            # Axes
            ax.set_ylim(ylim)
            ax.set_xlabel('Lag')
        # Finish
        plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
        plt.savefig(outfile, dpi=400)
        plt.close()

        # Sky line QA (just one object)
        if slit_cen:
            iobj = 0
        else:
            iobj = 0
            specobj = this_specobjs[iobj]

        if len(this_flex_dict['shift']) == 0:
            return

        # Repackage
        sky_spec = this_flex_dict['sky_spec'][iobj]
        arx_spec = this_flex_dict['arx_spec'][iobj]

        # Sky lines
        sky_lines = np.array([3370.0, 3914.0, 4046.56, 4358.34, 5577.338, 6300.304,
                              7340.885, 7993.332, 8430.174, 8919.610, 9439.660,
                              10013.99, 10372.88])*units.AA
        dwv = 20.*units.AA
        gdsky = np.where((sky_lines > sky_spec.wvmin) & (sky_lines < sky_spec.wvmax))[0]
        if len(gdsky) == 0:
            msgs.warn("No sky lines for Flexure QA")
            return
        if len(gdsky) > 6:
            idx = np.array([0, 1, len(gdsky)//2, len(gdsky)//2+1, -2, -1])
            gdsky = gdsky[idx]

        # Outfile
        outfile = qa.set_qa_filename(basename, method+'_sky', det=det,slit=(slit + 1), out_dir=out_dir)
        # Figure
        plt.figure(figsize=(8, 5.0))
        plt.clf()
        nrow, ncol = 2, 3
        gs = gridspec.GridSpec(nrow, ncol)
        if slit_cen:
            plt.suptitle('Sky Comparison for Slit Center', y=1.05)
        else:
            plt.suptitle('Sky Comparison for {:s}'.format(specobj.NAME), y=1.05)

        for ii, igdsky in enumerate(gdsky):
            skyline = sky_lines[igdsky]
            ax = plt.subplot(gs[ii//ncol, ii % ncol])
            # Norm
            pix = np.where(np.abs(sky_spec.wavelength-skyline) < dwv)[0]
            f1 = np.sum(sky_spec.flux[pix])
            f2 = np.sum(arx_spec.flux[pix])
            norm = f1/f2
            # Plot
            ax.plot(sky_spec.wavelength[pix], sky_spec.flux[pix], 'k-', label='Obj',
                    drawstyle='steps-mid')
            pix2 = np.where(np.abs(arx_spec.wavelength-skyline) < dwv)[0]
            ax.plot(arx_spec.wavelength[pix2], arx_spec.flux[pix2]*norm, 'r-', label='Arx',
                    drawstyle='steps-mid')
            # Axes
            ax.xaxis.set_major_locator(plt.MultipleLocator(dwv.value))
            ax.set_xlabel('Wavelength')
            ax.set_ylabel('Counts')

        # Legend
        plt.legend(loc='upper left', scatterpoints=1, borderpad=0.3,
                   handletextpad=0.3, fontsize='small', numpoints=1)

        # Finish
        plt.savefig(outfile, dpi=400)
        plt.close()
        #plt.close()

    plt.rcdefaults()

    return


def flexure_qa_oldbuggyversion(specobjs, maskslits, basename, det, flex_list, slit_cen=False):
    """ QA on flexure measurement

    Parameters
    ----------
    det
    flex_list : list
      list of dict containing flexure results
    slit_cen : bool, optional
      QA on slit center instead of objects

    Returns
    -------

    """
    plt.rcdefaults()
    plt.rcParams['font.family']= 'times new roman'

    # Grab the named of the method
    method = inspect.stack()[0][3]
    #
    gdslits = np.where(~maskslits)[0]
    for sl in range(len(specobjs)):
        if sl not in gdslits:
            continue
        if specobjs[sl][0] is None:
            continue
        # Setup
        if slit_cen:
            nobj = 1
            ncol = 1
        else:
            nobj = len(specobjs[sl])
            ncol = min(3, nobj)
        #
        if nobj==0:
            continue
        nrow = nobj // ncol + ((nobj % ncol) > 0)

        # Get the flexure dictionary
        flex_dict = flex_list[sl]

        # Outfile
        outfile = qa.set_qa_filename(basename, method+'_corr', det=det,
                                       slit=specobjs[sl][0].SLITID)

        plt.figure(figsize=(8, 5.0))
        plt.clf()
        gs = gridspec.GridSpec(nrow, ncol)

        # Correlation QA
        for o in range(nobj):
            ax = plt.subplot(gs[o//ncol, o % ncol])
            # Fit
            fit = flex_dict['polyfit'][o]
            xval = np.linspace(-10., 10, 100) + flex_dict['corr_cen'][o] #+ flex_dict['shift'][o]
            #model = (fit[2]*(xval**2.))+(fit[1]*xval)+fit[0]
            model = utils.func_val(fit, xval, 'polynomial')
            mxmod = np.max(model)
            ylim = [np.min(model/mxmod), 1.3]
            ax.plot(xval-flex_dict['corr_cen'][o], model/mxmod, 'k-')
            # Measurements
            ax.scatter(flex_dict['subpix'][o]-flex_dict['corr_cen'][o],
                       flex_dict['corr'][o]/mxmod, marker='o')
            # Final shift
            ax.plot([flex_dict['shift'][o]]*2, ylim, 'g:')
            # Label
            if slit_cen:
                ax.text(0.5, 0.25, 'Slit Center', transform=ax.transAxes, size='large', ha='center')
            else:
                ax.text(0.5, 0.25, '{:s}'.format(specobjs[sl][o].NAME), transform=ax.transAxes, size='large', ha='center')
            ax.text(0.5, 0.15, 'flex_shift = {:g}'.format(flex_dict['shift'][o]),
                    transform=ax.transAxes, size='large', ha='center')#, bbox={'facecolor':'white'})
            # Axes
            ax.set_ylim(ylim)
            ax.set_xlabel('Lag')

        # Finish
        plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
        plt.savefig(outfile, dpi=400)
        plt.close()

        # Sky line QA (just one object)
        if slit_cen:
            o = 0
        else:
            o = 0
            specobj = specobjs[sl][o]
        sky_spec = flex_dict['sky_spec'][o]
        arx_spec = flex_dict['arx_spec'][o]

        # Sky lines
        sky_lines = np.array([3370.0, 3914.0, 4046.56, 4358.34, 5577.338, 6300.304,
                              7340.885, 7993.332, 8430.174, 8919.610, 9439.660,
                              10013.99, 10372.88])*units.AA
        dwv = 20.*units.AA
        gdsky = np.where((sky_lines > sky_spec.wvmin) & (sky_lines < sky_spec.wvmax))[0]
        if len(gdsky) == 0:
            msgs.warn("No sky lines for Flexure QA")
            return
        if len(gdsky) > 6:
            idx = np.array([0, 1, len(gdsky)//2, len(gdsky)//2+1, -2, -1])
            gdsky = gdsky[idx]

        # Outfile
        outfile = qa.set_qa_filename(basename, method+'_sky', det=det,
                                       slit=specobjs[sl][0].SLITID)
        # Figure
        plt.figure(figsize=(8, 5.0))
        plt.clf()
        nrow, ncol = 2, 3
        gs = gridspec.GridSpec(nrow, ncol)
        if slit_cen:
            plt.suptitle('Sky Comparison for Slit Center', y=1.05)
        else:
            plt.suptitle('Sky Comparison for {:s}'.format(specobj.NAME), y=1.05)

        for ii, igdsky in enumerate(gdsky):
            skyline = sky_lines[igdsky]
            ax = plt.subplot(gs[ii//ncol, ii % ncol])
            # Norm
            pix = np.where(np.abs(sky_spec.wavelength-skyline) < dwv)[0]
            f1 = np.sum(sky_spec.flux[pix])
            f2 = np.sum(arx_spec.flux[pix])
            norm = f1/f2
            # Plot
            ax.plot(sky_spec.wavelength[pix], sky_spec.flux[pix], 'k-', label='Obj',
                    drawstyle='steps-mid')
            pix2 = np.where(np.abs(arx_spec.wavelength-skyline) < dwv)[0]
            ax.plot(arx_spec.wavelength[pix2], arx_spec.flux[pix2]*norm, 'r-', label='Arx',
                    drawstyle='steps-mid')
            # Axes
            ax.xaxis.set_major_locator(plt.MultipleLocator(dwv.value))
            ax.set_xlabel('Wavelength')
            ax.set_ylabel('Counts')

        # Legend
        plt.legend(loc='upper left', scatterpoints=1, borderpad=0.3,
                   handletextpad=0.3, fontsize='small', numpoints=1)

        # Finish
        plt.savefig(outfile, dpi=400)
        plt.close()
        #plt.close()

    plt.rcdefaults()

